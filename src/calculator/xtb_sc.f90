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

!> module xtb_sc
!> A module containing routines for
!> system calls to the xtb code

!=========================================================================================!
module xtb_sc
  use iso_fortran_env,only:wp => real64
  use strucrd
  use calc_type
  use iomod,only:makedir,directory_exist,remove,command
  implicit none
  !>--- private module variables and parameters
  private
  integer,parameter :: nf = 3
  character(len=*),parameter :: xtbfiles(nf) = [&
          & 'charges    ','xtbinp.grad','xtbrestart ']
  character(len=3),parameter :: xtb = 'xtb'
  character(len=10),parameter :: xyzn = 'xtbinp.xyz'
  character(len=13),parameter :: gf = 'xtbinp.engrad'

  public :: xtb_engrad

!========================================================================================!
!========================================================================================!
contains  !>--- Module routines start here
!========================================================================================!
!========================================================================================!

  subroutine xtb_engrad(mol,calc,energy,grad,iostatus)
    use iso_fortran_env,only:wp => real64
    use strucrd
    use calc_type
    use iomod,only:makedir,directory_exist,remove

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
    call xtb_setup(mol,calc)
    !$omp end critical

    !>--- do the systemcall
    call initsignal()
    call command(calc%systemcall,iostatus)
    if (iostatus /= 0) return

    !>--- read energy and gradient
    !$omp critical
    call rd_xtb_engrad(mol,calc,energy,grad,iostatus)
    !$omp end critical
    if (iostatus /= 0) return

    !>--- read WBOs?
    !$omp critical
    call rd_xtb_wbo(mol,calc,iostatus)
    !$omp end critical
    if (iostatus /= 0) return

    return
  end subroutine xtb_engrad

!========================================================================================!
  subroutine xtb_setup(mol,calc)
    use iso_fortran_env,only:wp => real64
    use strucrd
    use calc_type
    use iomod,only:makedir,directory_exist,remove

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
      calc%binary = xtb
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
      !write(*,*) trim(cpath)//sep//trim(xtbfiles(i))
      call remove(trim(cpath)//sep//trim(xtbfiles(i)))
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
    call mol%write(fname)
    deallocate (fname)

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

    !>--- add other call information
    calc%systemcall = trim(calc%systemcall)//' '//xyzn
    !>--- chrg and uhf
    if (calc%chrg .ne. 0) then
      write (num,'(i0)') calc%chrg
      calc%systemcall = trim(calc%systemcall)//' '//'--chrg'
      calc%systemcall = trim(calc%systemcall)//' '//trim(num)
    end if
    if (calc%uhf .ne. 0) then
      write (num,'(i0)') calc%uhf
      calc%systemcall = trim(calc%systemcall)//' '//'--uhf'
      calc%systemcall = trim(calc%systemcall)//' '//trim(num)
    end if
    !>--- user-set flags
    if (allocated(calc%other)) then
      calc%systemcall = trim(calc%systemcall)//' '//trim(calc%other)
    end if
    !>--- don't miss the --grad flag!
    if (index(calc%systemcall,'-grad') .eq. 0) then
      calc%systemcall = trim(calc%systemcall)//' '//'--grad'
    end if

    !>--- add printout information
    calc%systemcall = trim(calc%systemcall)//' '//'> xtb.out'
    calc%systemcall = trim(calc%systemcall)//dev0

    !write (*,*) calc%systemcall
    return
  end subroutine xtb_setup

!========================================================================================!
! subroutine rd_xtb_engrad
! read xtb's energy and Cartesian gradient from file
! xtb's *.engrad format is used for this
  subroutine rd_xtb_engrad(mol,calc,energy,grad,iostatus)
    use iso_fortran_env,only:wp => real64
    use strucrd
    use calc_type
    use iomod,only:makedir,directory_exist,remove
    use gradreader_module,only:rd_grad_engrad

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

    c = 0
    open (newunit=ich,file=calc%gradfile)
    call rd_grad_engrad(ich,mol%nat,energy,grad,iostatus)
    close (ich)

    return
  end subroutine rd_xtb_engrad
!========================================================================================!
! subroutine rd_xtb_wbo
! helper routine ro read xtb WBOs
  subroutine rd_xtb_wbo(mol,calc,iostatus)
    implicit none
    type(coord) :: mol
    type(calculation_settings) :: calc
    integer,intent(out) :: iostatus

    real(wp) :: dum
    character(len=:),allocatable :: wbofile
    character(len=128) :: atmp

    integer :: i,j,k,l,ich,och,io
    logical :: ex
    call initsignal()

    iostatus = 0

    if (calc%rdwbo) then
      if (allocated(calc%calcspace)) then
        wbofile = trim(calc%calcspace)//sep//'wbo'
      else
        wbofile = 'wbo'
      end if
    else
      return
    end if

    inquire (file=wbofile,exist=ex)
    if (.not.ex) then
      iostatus = 1
      return
    end if

    if (allocated(calc%wbo)) deallocate (calc%wbo)
    allocate (calc%wbo(mol%nat,mol%nat),source=0.0_wp)

    open (newunit=ich,file=wbofile)
    do
      read (ich,'(a)',iostat=io) atmp
      if (io < 0) exit
      read (atmp,*) i,j,dum
      calc%wbo(i,j) = dum
      calc%wbo(j,i) = dum
    end do
    close (ich)

  end subroutine rd_xtb_wbo

!========================================================================================!
end module xtb_sc
