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

!====================================================!
! module generic_sc
! A module containing routines for
! system calls to a generic script
!====================================================!

module generic_sc

  use iso_fortran_env,only:wp => real64
  use strucrd
  use calc_type
  use iomod,only:makedir,directory_exist,remove
  implicit none

!=========================================================================================!
  !--- private module variables and parameters
  private
  integer :: i,j,k,l,ich,och,io
  logical :: ex

  integer,parameter :: nf = 2
  character(len=20),parameter :: rmfiles(nf) = [&
          & 'genericinp.grad','genericinp.engrad']
  character(len=6),parameter :: runscript = 'run.sh'
  character(len=14),parameter :: xyzn = 'genericinp.xyz'
  character(len=17),parameter :: gf = 'genericinp.engrad'

  public :: generic_engrad

contains
!========================================================================================!
  subroutine generic_engrad(mol,calc,energy,grad,iostatus)
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

    iostatus = 0
    
    !>--- setup system call information
    !$omp critical
    call generic_setup(mol,calc)
    !$omp end critical

    !>--- do the systemcall
    call initsignal()
    call execute_command_line(calc%systemcall,exitstat=iostatus)
    if (iostatus /= 0) return

    !>--- read energy and gradient
    !$omp critical
    call rd_generic_engrad(mol,calc,energy,grad,iostatus)
    !$omp end critical
    if (iostatus /= 0) return

    !>--- read WBOs?
    !!$omp critical
    !call rd_generic_wbo(mol,calc,iostatus)
    !!$omp end critical
    !if (iostatus /= 0) return

    return
  end subroutine generic_engrad

!========================================================================================!
  subroutine generic_setup(mol,calc)
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
    character(len=:),allocatable :: tmprelpath
    character(len=512) :: thispath
    character(len=10) :: num
    integer :: i,j,k,ich,och,io
    logical :: ex

    call initsignal()

    !>--- set default binary if not present
    if (.not. allocated(calc%binary)) then
      calc%binary = runscript
    end if

    !>--- important for generec calls:
    !     check if the requested binary/script is in the current working directory
    !     if so, convert to absolute path
    tmprelpath = '.'//sep//trim(calc%binary)
    inquire(file=tmprelpath,exist=ex)
    if(ex)then
      call getcwd(thispath)
      calc%binary = trim(thispath)//sep//trim(calc%binary)
    endif 
    deallocate(tmprelpath)

    !>--- check for the calculation space
    if (allocated(calc%calcspace)) then
      ex = directory_exist(calc%calcspace)
      if (.not. ex) then
        io = makedir(trim(calc%calcspace))
      end if
      cpath = calc%calcspace
    else
      cpath = ''
    end if
    !>--- cleanup old files
    do i = 1,nf
      !write(*,*) trim(cpath)//sep//trim(rmfiles(i)) 
      call remove(trim(cpath)//sep//trim(rmfiles(i)))
    end do
    deallocate (cpath)

    !>--- construct path information and write coord file
    if (.not. allocated(calc%calcfile)) then
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

    !>--- if the systemcall was already set up, return
    if (allocated(calc%systemcall)) return

    !>--- construct path information for sys-call
    if (allocated(calc%calcspace)) then
      calc%systemcall = 'cd '//calc%calcspace//' &&'
      calc%systemcall = trim(calc%systemcall)//' '//trim(calc%binary)
    else
      calc%systemcall = trim(calc%binary)
    end if

    !>--- add other call information
    calc%systemcall = trim(calc%systemcall)//' '//xyzn
    !>--- user-set flags
    if (allocated(calc%other)) then
      calc%systemcall = trim(calc%systemcall)//' '//trim(calc%other)
    end if

    !>--- add printout information
    calc%systemcall = trim(calc%systemcall)//' '//'> generic.out'
    calc%systemcall = trim(calc%systemcall)//dev0

    !write (*,*) calc%systemcall
    return
  end subroutine generic_setup

!========================================================================================!
! subroutine rd_generic_engrad
! read xtb's energy and Cartesian gradient from file
! xtb's *.engrad format is used for this
  subroutine rd_generic_engrad(mol,calc,energy,grad,iostatus)
    use iso_fortran_env,only:wp => real64
    use strucrd
    use calc_type
    use iomod,only:makedir,directory_exist,remove
    use gradreader_module
    implicit none
    type(coord) :: mol
    type(calculation_settings) :: calc
    real(wp),intent(inout) :: energy
    real(wp),intent(inout) :: grad(3,mol%nat)
    integer,intent(out) :: iostatus
    integer :: n,c
    real(wp) :: dum
    character(len=128) :: atmp
    character(len=512) :: tmpgradfile,tmpefile

    integer :: i,j,k,ich,och,io
    logical :: ex

    call initsignal()

    iostatus = 0
    tmpgradfile = ''
    tmpefile = ''
    if (.not. allocated(calc%gradfile)) then
      calc%gradfile = gf
    endif  
    if (allocated(calc%calcspace)) then
      tmpgradfile = trim(calc%calcspace)//sep//trim(calc%gradfile)
      if(allocated(calc%efile)) &
      & tmpefile  = trim(calc%calcspace)//sep//trim(calc%efile)
    else
      tmpgradfile = trim(calc%gradfile)
      if(allocated(calc%efile)) &
      & tmpefile  = trim(calc%efile)
    end if

    inquire (file=trim(tmpgradfile),exist=ex)
    if (.not. ex) then
      iostatus = 1
      return
    end if

    open (newunit=ich,file=trim(tmpgradfile))
    select case( calc%gradtype )
    case( gradtype%engrad )
      call rd_grad_engrad(ich,mol%nat,energy,grad,iostatus)
    case( gradtype%turbomole )
       
    case default
      if(.not.allocated(calc%efile))then
         iostatus = 1
         return
      endif
      call rd_efile(trim(tmpefile),energy,iostatus)
      if(allocated(calc%gradkey))then
        call rd_grad_generic(ich,mol%nat,grad, &
        & calc%gradkey,calc%gradfmt,iostatus)
      else
        call rd_grad_generic(ich,mol%nat,grad, &
        & '',calc%gradfmt,iostatus)
      endif 
    end select
    close (ich)

    return
  end subroutine rd_generic_engrad
!========================================================================================!
! subroutine rd_generic_wbo
! helper routine ro read generic WBOs
  subroutine rd_generic_wbo(mol,calc,iostatus)
    implicit none
    type(coord) :: mol
    type(calculation_settings) :: calc
    integer,intent(out) :: iostatus
    integer :: i,j
    real(wp) :: dum
    character(len=:),allocatable :: wbofile
    character(len=128) :: atmp

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
    if (.not. ex) then
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

  end subroutine rd_generic_wbo

!========================================================================================!
end module generic_sc
