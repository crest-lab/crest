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

  use iso_fortran_env,only:wp => real64
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

    logical :: loadnew
    iostatus = 0
    
    !>--- setup system call information
    !$omp critical
    call tblite_init(calc,loadnew)
    if(loadnew)then
    call tblite_setup(mol,calc%chrg,calc%uhf,calc%tblitelvl,calc%etemp, &
    &    calc%ctx,calc%wfn,calc%tbcalc)
    endif
    !$omp end critical

    !>--- do the engrad call
    call initsignal()

  

    !>--- postprocessing, getting other data


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
      if( calc%tbliteclean ) loadnew = .true.
    end subroutine tblite_init
  end subroutine tblite_engrad


!========================================================================================!



!========================================================================================!
end module api_engrad
