!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2023 Philipp Pracht
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

module xhcff_api
  use iso_fortran_env,only:wp => real64,stdout => output_unit
  use strucrd
#ifdef WITH_XHCFF
  use xhcff_interface
#endif
  implicit none
  private

#ifndef WITH_XHCFF
  !> this is a placeholder if no xhcff module is used!
  type :: xhcff_calculator
    integer :: id = 0
  end type xhcff_calculator
#endif

  public :: xhcff_calculator  !> if compiled without(!!!) -DWITH_XHCFF=true this will export
                        !> the placeholder from above. Otherwise it will re-export
                        !> the type from xhcff_interface


  public :: xhcff_setup,xhcff_sp,xhcff_print

!========================================================================================!
!========================================================================================!
contains  !>--- Module routines start here
!========================================================================================!
!========================================================================================!

  subroutine xhcff_setup(mol, xhcff, pressure, gridpts, proberad, vdwset, pr, iunit, iostatus)
    implicit none
    type(coord),intent(in)  :: mol
    real(wp), intent(in) :: pressure !> pressure
    integer, intent(in) :: gridpts
    real(wp), intent(in) :: proberad
    integer, intent(in) :: vdwset
    type(xhcff_calculator),intent(inout) :: xhcff
    integer, intent(inout) :: iostatus
    logical, intent(in) :: pr
    integer,intent(in) :: iunit
#ifdef WITH_XHCFF
    !> initialize XHCFF
    call xhcff%init(mol%nat,mol%at,mol%xyz, &
     & pressure, gridpts, proberad, verbose=pr,iunit=iunit,vdwset=vdwset,iostat=iostatus)

#else /* WITH_XHCFF */
    write (stdout,*) 'Error: Compiled without XHCFF-lib support!'
    write (stdout,*) 'Use -DWITH_XHCFF=true in the setup to enable this function'
    error stop
#endif
  end subroutine xhcff_setup

!========================================================================================!

  subroutine xhcff_sp(mol,xhcff,energy,gradient,iostatus)
!********************************************************
!* The actual energy+gradient call to xhcff-lib.
!* Requires the xhcff_calculator object to be set up already.
!* Note that the original xhcff has no contribution to
!* the energy, only to the gradient
!********************************************************
    implicit none
    !> INPUT
    type(coord),intent(in)  :: mol
    type(xhcff_calculator),intent(inout) :: xhcff
    !> OUTPUT
    real(wp),intent(out) :: energy
    real(wp),intent(out) :: gradient(3,mol%nat)
    integer,intent(out) :: iostatus
    !> LOCAL
    logical :: fail
    energy = 0.0_wp
    gradient = 0.0_wp
    iostatus = 0
    fail = .false.
#ifdef WITH_XHCFF
!TODO update
     call xhcff%singlepoint(mol%nat,mol%at,mol%xyz,energy,gradient,iostatus)
#else
    write (stdout,*) 'Error: Compiled without XHCFF-lib support!'
    write (stdout,*) 'Use -DWITH_XHCFF=true in the setup to enable this function'
    error stop
#endif
  end subroutine xhcff_sp

!========================================================================================!

  subroutine xhcff_print(iunit,xhcff)
    implicit none
    integer,intent(in) :: iunit
    type(xhcff_calculator),intent(in) :: xhcff
#ifdef WITH_XHCFF
    call xhcff%info(iunit)
#endif
    return
  end subroutine xhcff_print

!========================================================================================!
!========================================================================================!
end module xhcff_api

