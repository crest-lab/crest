!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2021 - 2022 Philipp Pracht
!
! crest is free software: you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! crest is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with crest.  If not, see <https://www.gnu.org/licenses/>.
!================================================================================!

!> module turbom_sc
!> A module containing routines for
!> system calls to the code that uses Turbomole-like in & output conventions
!> TODO: for actual turbomole calculations the input files must be copied (cefine)

!=========================================================================================!
module turbom_sc
  use iso_fortran_env,only:wp => real64
  use strucrd
  use calc_type
  use iomod,only:makedir,directory_exist,remove,command,checkprog_silent
  use gradreader_module,only:rd_grad_tm
  implicit none
  !>--- private module variables and parameters
  private
  integer,parameter :: nf = 6
  character(len=*),parameter :: oldfiles(nf) = [&
          & 'energy       ','ceh.charges  ', &
          & 'output       ','.data        ', &
          & 'NOT_CONVERGED','gp3restart   ' ]
  character(len=*),parameter :: ridft = 'ridft' !> Turbomoles 'ridft'
  character(len=*),parameter :: xyzn = 'coord'  !> input coords must be in coord
  character(len=*),parameter :: ef = 'energy'   !> energy will be read from file energy
  character(len=*),parameter :: gf = 'gradient' !> gradient will be read from file gradient

  public :: turbom_engrad

!========================================================================================!
!========================================================================================!
contains  !>--- Module routines start here
!========================================================================================!
!========================================================================================!

  subroutine turbom_engrad(mol,calc,energy,grad,iostatus)
    implicit none
    type(coord) :: mol
    type(calculation_settings) :: calc

    real(wp),intent(inout) :: energy
    real(wp),intent(inout) :: grad(3,mol%nat)
    integer,intent(out) :: iostatus

    integer :: i,j,k,l,ich,och,io
    logical :: ex

    iostatus = 0

    !>--- setup system call information
    !$omp critical
    call turbom_setup(mol,calc)
    !$omp end critical

    !>--- do the systemcall
    call initsignal()
    call command(calc%systemcall,iostatus)
    if (iostatus /= 0) return

    !>--- read energy and gradient
    !$omp critical
    call rd_turbom_engrad(mol,calc,energy,grad,iostatus)
    !$omp end critical
    if (iostatus /= 0) return

    !>--- read WBOs?
    !!$omp critical
    !call rd_turbom_wbo(mol,calc,iostatus)
    !!$omp end critical
    !if (iostatus /= 0) return

    return
  end subroutine turbom_engrad

!========================================================================================!
  subroutine turbom_setup(mol,calc)
    implicit none
    type(coord) :: mol
    type(calculation_settings) :: calc

    integer :: l
    character(len=:),allocatable :: fname
    character(len=:),allocatable :: cpath
    character(len=10) :: num
    integer :: i,j,k,ich,och,io
    logical :: ex

    call initsignal()

    !>--- set default binary if not present
    if (.not.allocated(calc%binary)) then
      calc%binary = ridft
    end if

    !>--- check for the calculation space
    if (allocated(calc%calcspace)) then
      ex = directory_exist(calc%calcspace)
      if (.not.ex) then
        io = makedir(trim(calc%calcspace))
      end if
      cpath = calc%calcspace
    else
      cpath = ''
    end if
    !>--- cleanup old files
    do i = 1,nf
      call remove(trim(cpath)//sep//trim(oldfiles(i)))
    end do
    deallocate (cpath)

    !>--- construct path information and write coord file
    if (.not.allocated(calc%calcfile)) then
      if (allocated(calc%calcspace)) then
        l = len_trim(calc%calcspace)
        fname = trim(calc%calcspace)
        if (calc%calcspace(l:l) == sep) then
          fname = trim(fname)//xyzn
        else
          fname = trim(fname)//sep//xyzn
        end if
      else
        fname = xyzn
      end if
      calc%calcfile = fname
    else
      fname = calc%calcfile
    end if
    call mol%write(fname) !> should write a "coord" file, must be called for each SP
    deallocate (fname)

    !>--- write charge and uhf files
    if (calc%chrg /= 0) then
      call touch_chrg_tm(calc,calc%chrg)
    end if
    if (calc%uhf /= 0) then
      call touch_uhf_tm(calc,calc%uhf)
    endif

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<!
    !>--- if the systemcall was already set up, return
    if (allocated(calc%systemcall)) return
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<!

    !>--- construct path information for sys-call
    if (allocated(calc%calcspace)) then
      calc%systemcall = 'cd '//calc%calcspace//' &&'
      calc%systemcall = trim(calc%systemcall)//' '//trim(calc%binary)
    else
      calc%systemcall = trim(calc%binary)
    end if

    !>--- check if the binary exists
    call checkprog_silent(trim(calc%binary),.true.,io)
    if(io .ne. 0) error stop

    !>--- add other call information
    calc%systemcall = trim(calc%systemcall)//' '//xyzn
    !>--- user-set flags
    if (allocated(calc%other)) then
      calc%systemcall = trim(calc%systemcall)//' '//trim(calc%other)
    end if

    !>--- add printout information
    calc%systemcall = trim(calc%systemcall)//' '//'> output'
    calc%systemcall = trim(calc%systemcall)//dev0

    return
  end subroutine turbom_setup

!========================================================================================!
  subroutine rd_turbom_engrad(mol,calc,energy,grad,iostatus)
!*****************************************************
!* subroutine rd_turbom_engrad
!* read the energy and Cartesian gradient from file
!* Turbomole-style format is used for this
!*****************************************************
    implicit none
    type(coord) :: mol
    type(calculation_settings) :: calc
    real(wp),intent(inout) :: energy
    real(wp),intent(inout) :: grad(3,mol%nat)
    integer,intent(out) :: iostatus
    integer :: n,c
    real(wp) :: dum
    character(len=128) :: atmp

    integer :: i,j,k,ich,och,io
    logical :: ex

    call initsignal()

    iostatus = 0
    energy = 0.0_wp
    grad = 0.0_wp

!>--- read the energy file in TM format
    if (.not.allocated(calc%efile)) then
      if (allocated(calc%calcspace)) then
        calc%efile = trim(calc%calcspace)//sep//ef
      else
        calc%efile = ef
      end if
    end if

    inquire (file=calc%efile,exist=ex)
    if (.not.ex) then
      iostatus = 1
      return
    end if
    open (newunit=ich,file=calc%efile)
    call rd_energy_tm(ich,energy,iostatus)
    close (ich)

!>--- sometimes we don't want to calculate and read gradients for expensive calculations
!>--- where we are only interested in the total energy. Hence, reading it can be skipped
    if (calc%rdgrad) then

      if (.not.allocated(calc%gradfile)) then
        if (allocated(calc%calcspace)) then
          calc%gradfile = trim(calc%calcspace)//sep//gf
        else
          calc%gradfile = gf
        end if
      end if

      inquire (file=calc%gradfile,exist=ex)
      if (.not.ex) then
        iostatus = 1
        return
      end if

      open (newunit=ich,file=calc%gradfile)
      call rd_grad_tm(ich,mol%nat,energy,grad,iostatus)
      close (ich)
    end if

    return
  end subroutine rd_turbom_engrad

!========================================================================================!
!  subroutine rd_turbom_wbo(mol,calc,iostatus)
!    implicit none
!    type(coord) :: mol
!    type(calculation_settings) :: calc
!    integer,intent(out) :: iostatus
!
!    real(wp) :: dum
!    character(len=:),allocatable :: wbofile
!    character(len=128) :: atmp
!
!    integer :: i,j,k,l,ich,och,io
!    logical :: ex
!    call initsignal()
!
!    iostatus = 0
!
!    if (calc%rdwbo) then
!      if (allocated(calc%calcspace)) then
!        wbofile = trim(calc%calcspace)//sep//'wbo'
!      else
!        wbofile = 'wbo'
!      end if
!    else
!      return
!    end if
!
!    inquire (file=wbofile,exist=ex)
!    if (.not.ex) then
!      iostatus = 1
!      return
!    end if
!
!    if (allocated(calc%wbo)) deallocate (calc%wbo)
!    allocate (calc%wbo(mol%nat,mol%nat),source=0.0_wp)
!
!    open (newunit=ich,file=wbofile)
!    do
!      read (ich,'(a)',iostat=io) atmp
!      if (io < 0) exit
!      read (atmp,*) i,j,dum
!      calc%wbo(i,j) = dum
!      calc%wbo(j,i) = dum
!    end do
!    close (ich)
!
!  end subroutine rd_turbom_wbo

!========================================================================================!
  subroutine touch_chrg_tm(calc,chrg)
    implicit none
    type(calculation_settings) :: calc
    integer,intent(in) :: chrg
    character(len=:),allocatable :: fname
    integer :: ich
    if (allocated(calc%calcspace)) then
      fname = trim(calc%calcspace)//sep//'.CHRG'
    else
      fname = '.CHRG'
    end if
    open (newunit=ich,file=fname)
    write (ich,*) chrg
    close (ich)
    deallocate (fname)
  end subroutine touch_chrg_tm

  subroutine touch_uhf_tm(calc,uhf)
    implicit none
    type(calculation_settings) :: calc
    integer,intent(in) :: uhf
    character(len=:),allocatable :: fname
    integer :: ich
    if (allocated(calc%calcspace)) then
      fname = trim(calc%calcspace)//sep//'.UHF'
    else
      fname = '.UHF'
    end if
    open (newunit=ich,file=fname)
    write (ich,*) uhf
    close (ich)
    deallocate (fname)
  end subroutine touch_uhf_tm

  subroutine rd_energy_tm(ich,energy,iostat)
!*******************************************
!* read the 'energy' file
!* The most current energy should be the 
!* second entry in the second-to-last line
!*******************************************
     implicit none
     integer,intent(in) :: ich
     real(wp),intent(out) :: energy
     integer,intent(out) :: iostat
     integer :: i,j,nlines,io
     character(len=126) :: line,oldline
     real(wp) :: dum(2)
     energy = 0.0_wp
     iostat = 1
     nlines = 0
     do
       read (ich,'(a)',iostat=io) line
       if(io < 0 ) exit !> EOF exit
       if(index(line,'$end').ne.0)then
         read(oldline,*) j,energy,dum(1:2)
         iostat = 0
       endif
       oldline = line
     enddo 
  end subroutine rd_energy_tm
 
!========================================================================================!
!========================================================================================!
end module turbom_sc
