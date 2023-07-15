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
  type :: xhcff_data
    integer :: id = 0
  end type xhcff_data
#endif

  public :: xhcff_data  !> if compiled without(!!!) -DWITH_XHCFF=true this will export
                        !> the placeholder from above. Otherwise it will re-export
                        !> the type from xhcff_interface

  public :: xhcff_setup,xhcff_sp,xhcff_print

!========================================================================================!
!========================================================================================!
contains  !>--- Module routines start here
!========================================================================================!
!========================================================================================!

  subroutine xhcff_setup(mol,xhcff)
    implicit none
    type(coord),intent(in)  :: mol
    type(xhcff_data),intent(inout) :: xhcff
#ifdef WITH_XHCFF
!TODO pressure and grid size must probably be passed here
    !> initialize XHCFF
    call xhcff_initialize(mol%nat,mol%at,mol%xyz,xhcff)

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
!* Requires the xhcff_data object to be set up already.
!* Note that the original xhcff has no contribution to
!* the energy, only to the gradient
!********************************************************
    implicit none
    !> INPUT
    type(coord),intent(in)  :: mol
    type(xhcff_data),intent(inout) :: xhcff
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
      write(*,*) 'TODO: implement XHCFF'
!     call xhcff_singlepoint(mol%nat,mol%at,mol%xyz,xhcff, &
!     &    proberad,pressure,energy,gradient,.false.,iostatus)
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
    type(xhcff_data),intent(in) :: xhcff
#ifdef WITH_XHCFF
    call print_xhcff_results(xhcff,iunit)
#endif
    return
  end subroutine xhcff_print

!========================================================================================!
!========================================================================================!
end module xhcff_api

