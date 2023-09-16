!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2021 - 2023 Philipp Pracht
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

!> module orca_sc
!> A module containing routines for ORCA subrpocess calls

!=========================================================================================!
module orca_sc
  use iso_fortran_env,only:wp => real64,stderr => error_unit
  use strucrd
  use calc_type
  use iomod,only:makedir,directory_exist,remove,command
  implicit none
  !>--- private module variables and parameters
  private
  integer,parameter :: nf = 2
  character(len=*),parameter :: orcafiles(nf) = [&
          & 'ORCA.in    ','ORCA.engrad']
  character(len=3),parameter :: orca = 'orca'
  character(len=8),parameter :: xyzn = 'ORCA.in'
  character(len=11),parameter :: gf = 'ORCA.engrad'

  public :: ORCA_engrad

!========================================================================================!
!========================================================================================!
contains  !>--- Module routines start here
!========================================================================================!
!========================================================================================!

  subroutine ORCA_engrad(mol,calc,energy,grad,iostatus)
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
    call ORCA_setup(mol,calc)
    !$omp end critical

    !>--- do the systemcall
    call initsignal()
    call command(calc%systemcall,iostatus)
    if (iostatus /= 0) return

    !>--- read energy and gradient
    !$omp critical
    call rd_ORCA_engrad(mol,calc,energy,grad,iostatus)
    !$omp end critical
    if (iostatus /= 0) return

    !>--- read WBOs?
    !$omp critical
    call rd_ORCA_wbo(mol,calc,iostatus)
    !$omp end critical
    if (iostatus /= 0) return

    return
  end subroutine ORCA_engrad

!========================================================================================!

  subroutine ORCA_setup(mol,calc)
!***********************************
!* Write ORCA inputs, and construct
!* subprocess call
!***********************************
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
      if(.not.allocated(calc%ORCA%cmd))then
       write(stdout,'(a)') "**WARNING** ORCA cmd line not specified, falling back to 'orca'"
       calc%binary = orca
      else
       calc%binary = calc%ORCA%cmd
      endif
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
      call remove(trim(cpath)//sep//trim(orcafiles(i)))
    end do
    deallocate (cpath)

    !>--- construct path information and write ORCA.in file
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
    if(calc%uhf < 1) calc%uhf = 1 !> ORCA uses multiplicity, not n_alpha - n_beta!
    call calc%ORCA%write(fname,mol,calc%chrg,calc%uhf)
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

    calc%systemcall = trim(calc%systemcall)//' ORCA.in' 
       
    if(calc%ORCA%mpi .and. index(calc%binary,'oversubscribe').eq.0)then
      calc%systemcall = trim(calc%systemcall)//' --oversubscribe'
    endif

    !>--- add printout information
    calc%systemcall = trim(calc%systemcall)//' '//'> ORCA.out'
    calc%systemcall = trim(calc%systemcall)//trim(dev0)

    !write (*,*) calc%systemcall
    return
  end subroutine ORCA_setup

!========================================================================================!

  subroutine rd_ORCA_engrad(mol,calc,energy,grad,iostatus)
!**************************************************
!* Read energy and gradient from ORCA output.
!* ORCA will always write to a .engrad file,
!* which with our naming conventions should always
!* be called ORCA.engrad
!**************************************************
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
  end subroutine rd_ORCA_engrad

!========================================================================================!

  subroutine rd_ORCA_wbo(mol,calc,iostatus)
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
    if(calc%rdwbo)then
    iostatus = 1
    write (stderr,'(a)') '**ERROR** ORCA bond order reader not yet implemented'
    endif
    return
  end subroutine rd_ORCA_wbo

!========================================================================================!
!========================================================================================!
end module orca_sc

