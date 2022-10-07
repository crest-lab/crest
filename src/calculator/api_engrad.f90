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
! module api_engrad
! a collection of engrad calls for different apis
! this builds the communication between crests
! calculation settings and the api setups
!====================================================!

module api_engrad

  use iso_fortran_env,only:wp => real64,stdout=>output_unit
  use strucrd
  use calc_type
  use iomod,only:makedir,directory_exist,remove
  !> APIs
  use tblite_api

!=========================================================================================!
  implicit none
  !--- private module variables and parameters
  private
  integer :: i,j,k,l,ich,och,io
  logical :: ex


  public :: tblite_engrad

!=========================================================================================!
!=========================================================================================!
contains    !>--- Module routines start here
!=========================================================================================!
!=========================================================================================!
  subroutine tblite_engrad(mol,calc,energy,grad,iostatus)
    implicit none
    type(coord) :: mol
    type(calculation_settings) :: calc

    real(wp),intent(inout) :: energy
    real(wp),intent(inout) :: grad(3,mol%nat)
    integer,intent(out) :: iostatus

    character(len=:),allocatable :: cpath
    logical :: loadnew
    iostatus = 0
    
    !>--- setup system call information
    !$omp critical
    call tblite_init(calc,loadnew)
    !> tblite printout handling
    if(calc%ctx%unit .ne. stdout) close(calc%ctx%unit)
    if (allocated(calc%calcspace)) then
      ex = directory_exist(calc%calcspace)
      if (.not. ex) then
        io = makedir(trim(calc%calcspace))
      end if
      cpath = calc%calcspace//sep//'tblite.out'
    else
      cpath = 'tblite.out'
    end if
    open(newunit=calc%ctx%unit, file=cpath)
    deallocate (cpath)
    !> populate parameters and wavefunction
    if(loadnew)then
      call tblite_setup(mol,calc%chrg,calc%uhf,calc%tblitelvl,calc%etemp, &
      &    calc%ctx,calc%wfn,calc%tbcalc)
      call tblite_addsettings(calc%tbcalc,calc%maxscc)
    endif
    !$omp end critical

    !>--- do the engrad call
    call initsignal()
    call tblite_singlepoint(mol,calc%chrg,calc%uhf,calc%accuracy, &
    & calc%ctx,calc%wfn,calc%tbcalc,energy,grad,calc%tbres,iostatus)
    if(iostatus /= 0) return
     
    !>--- postprocessing, getting other data
    !> tbd

    return
  contains
    subroutine tblite_init(calc,loadnew)
      implicit none
      type(calculation_settings),intent(inout) :: calc
      logical,intent(out) :: loadnew
      loadnew = .false.
      if(.not.allocated(calc%wfn))then
      allocate(calc%wfn)
      loadnew = .true.
      endif
      if(.not.allocated(calc%tbcalc))then
      allocate(calc%tbcalc)
      loadnew=.true.
      endif
      if( .not.allocated(calc%ctx) )then
      allocate(calc%ctx)
      loadnew = .true.
      endif
      if( .not.allocated(calc%tbres) )then
      allocate(calc%tbres)
      loadnew=.true.
      endif
      if( calc%tbliteclean ) loadnew = .true.
    end subroutine tblite_init
  end subroutine tblite_engrad


!========================================================================================!



!========================================================================================!
end module api_engrad
