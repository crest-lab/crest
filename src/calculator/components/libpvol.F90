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

module libpvol_api
  use iso_fortran_env,only:wp => real64,stdout => output_unit
  use strucrd
#ifdef WITH_LIBPVOL
  use libpvol_interface
#endif
  implicit none
  private

#ifndef WITH_LIBPVOL
  !> this is a placeholder if no libpvol module is used!
  type :: libpvol_calculator
    integer :: id = 0
  end type libpvol_calculator
#endif

  public :: libpvol_calculator  !> if compiled without(!!!) -DWITH_LIBPVOL=true this will export
  !> the placeholder from above. Otherwise it will re-export
  !> the type from libpvol_interface

  public :: libpvol_setup,libpvol_sp,libpvol_print

!========================================================================================!
!========================================================================================!
contains  !>--- Module routines start here
!========================================================================================!
!========================================================================================!

  subroutine libpvol_setup(mol,libpvol,pressure,model,gridpts, &
  &                        proberad,vdwset,radscal,pr,iunit,iostatus)
    implicit none
    type(coord),intent(in)  :: mol
    real(wp),intent(in)    :: pressure !> pressure
    integer,intent(in)      :: model
    integer,intent(in)     :: gridpts
    real(wp),intent(in)    :: proberad
    integer,intent(in)     :: vdwset
    real(wp),intent(in)     :: radscal
    type(libpvol_calculator),intent(inout) :: libpvol
    integer,intent(inout)  :: iostatus
    logical,intent(in)     :: pr
    integer,intent(in)      :: iunit
    integer :: selectedGrid
#ifdef WITH_LIBPVOL
    !> check grid size
    call match_lebedev(gridpts,selectedGrid)
    !> initialize LIBPVOL
    select case (model)
    case (1) !> PV
      call libpvol%init(mol%nat,mol%at,mol%xyz,pressure,'PV', &
       &  gridpts=selectedGrid,proberad=proberad,verbose=pr,iunit=iunit, &
       &  printlevel=2,scaling=radscal,vdwset=vdwset,iostat=iostatus)
    case (0) !> XHCFF
      call libpvol%init(mol%nat,mol%at,mol%xyz,pressure,'XHCFF', &
       &  gridpts=selectedGrid,proberad=proberad,verbose=pr,iunit=iunit, &
       &  printlevel=2,scaling=radscal,vdwset=vdwset,iostat=iostatus)
    case default
      error stop 'Unkown libpvol model type'
    end select

#else /* WITH_LIBPVOL */
    write (stdout,*) 'Error: Compiled without LIBPVOL-lib support!'
    write (stdout,*) 'Use -DWITH_LIBPVOL=true in the setup to enable this function'
    error stop
#endif
  end subroutine libpvol_setup

!========================================================================================!

  subroutine libpvol_sp(mol,libpvol,energy,gradient,iostatus)
!********************************************************
!* The actual energy+gradient call to libpvol-lib.
!* Requires the libpvol_calculator object to be set up already.
!* Note that the original libpvol has no contribution to
!* the energy, only to the gradient
!********************************************************
    implicit none
    !> INPUT
    type(coord),intent(in)  :: mol
    type(libpvol_calculator),intent(inout) :: libpvol
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
#ifdef WITH_LIBPVOL
!TODO update
    call libpvol%singlepoint(mol%nat,mol%at,mol%xyz,energy,gradient,iostat=iostatus)
#else
    write (stdout,*) 'Error: Compiled without LIBPVOL-lib support!'
    write (stdout,*) 'Use -DWITH_LIBPVOL=true in the setup to enable this function'
    error stop
#endif
  end subroutine libpvol_sp

!========================================================================================!

  subroutine libpvol_print(iunit,libpvol)
    implicit none
    integer,intent(in) :: iunit
    type(libpvol_calculator),intent(in) :: libpvol
#ifdef WITH_LIBPVOL
    call libpvol%info(iunit)
#endif
    return
  end subroutine libpvol_print

!========================================================================================!

  subroutine match_lebedev(ngrid,match)
!********************************************************************
!* Returns the Lebedev-Laikov grid size closest to a provided ngrid
!********************************************************************
    implicit none
    integer,intent(in)  :: ngrid
    integer,intent(out) :: match
    !> Available Lebedev-Laikov grids
    integer,parameter :: gridSize(32) = [ &
       &    6,14,26,38,50,74,86,110, &
       &  146,170,194,230,266,302,350,434, &
       &  590,770,974,1202,1454,1730,2030,2354, &
       & 2702,3074,3470,3890,4334,4802,5294,5810]
    integer :: i,dist,best_idx
    dist = abs(ngrid-gridSize(1))
    best_idx = 1
    do i = 2,size(gridSize)
      if (abs(ngrid-gridSize(i)) < dist) then
        dist = abs(ngrid-gridSize(i))
        best_idx = i
      end if
    end do
    match = gridSize(best_idx)
  end subroutine match_lebedev

!========================================================================================!
!========================================================================================!
end module libpvol_api

