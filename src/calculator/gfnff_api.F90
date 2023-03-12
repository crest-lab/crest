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

!====================================================!
! module gfnff_api
! An interface to GFN-FF standalone calculations
!====================================================!

module gfnff_api
  use iso_fortran_env,only:wp => real64,stdout => output_unit
  use strucrd
#ifdef WITH_GFNFF
  use gfnff_interface
#endif
  implicit none
  private

  !> link to subproject routines/datatypes
  public :: gfnff_data,print_gfnff_results
  !> routines from this file
  public :: gfnff_api_setup
  public :: gfnff_sp

#ifndef WITH_GFNFF
  !> these are placeholders if no gfnff module is used!
  type :: gfnff_data
    integer :: id = 0
  end type gfnff_data
#endif


!========================================================================================!
!========================================================================================!
contains  !> MODULE PROCEDURES START HERE
!========================================================================================!
!========================================================================================!

  subroutine gfnff_api_setup(mol,chrg,ff_dat,io,pr,iunit)
    implicit none
    type(coord),intent(in)      :: mol
    integer,intent(in)          :: chrg
    integer,intent(out)         :: io
    logical,intent(in),optional :: pr
    integer,intent(in),optional :: iunit
    type(gfnff_data),allocatable,intent(inout) :: ff_dat
#ifdef WITH_GFNFF
    io = 0

    !> initialize parametrization and topology of GFN-FF
    call gfnff_initialize(mol%nat,mol%at,mol%xyz,ff_dat,'no file',.false., &
    & ichrg=chrg,print=pr,iostat=io,iunit=iunit)

#else /* WITH_GFNFF */
    write (stdout,*) 'Error: Compiled without GFN-FF support!'
    write (stdout,*) 'Use -DWITH_GFNFF=true in the setup to enable this function'
    error stop
#endif
  end subroutine gfnff_api_setup

!========================================================================================!

  subroutine gfnff_sp(mol,ff_dat,energy,gradient,iostatus)
    implicit none
    !> INPUT
    type(coord),intent(in)  :: mol
    type(gfnff_data),allocatable,intent(inout) :: ff_dat
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
#ifdef WITH_GFNFF
    call gfnff_singlepoint(mol%nat,mol%at,mol%xyz,ff_dat, &
    & energy,gradient,iostat=iostatus)
#else
    write (stdout,*) 'Error: Compiled without GFN-FF support!'
    write (stdout,*) 'Use -DWITH_GFNFF=true in the setup to enable this function'
    error stop
#endif
  end subroutine gfnff_sp



!========================================================================================!
!========================================================================================!
end module gfnff_api

