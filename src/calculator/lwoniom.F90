!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2023 Patryk Wesołowski, Philipp Pracht
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

module lwoniom_module
!*************************************************************
!* Thin wrapper around the lwONIOM library.
!* Most of the ONIOM bookkeeping (job lists, fragment update,
!* gradient projection, energy/gradient/Hessian reconstruction)
!* is handled by lwONIOM itself via its job table and driver
!* layer; crest only provides the potentials through an
!* adapter extending the abstract lwoniom_calculator class
!* (see crest_oniom_calc in calculator.F90).
!*************************************************************
  use crest_parameters
  use strucrd
#ifdef WITH_LWONIOM
  use lwoniom_interface
#endif
  implicit none
  private

#ifndef WITH_LWONIOM
  !> placeholders if no lwONIOM module is used!
  type :: lwoniom_frag_placeholder
    integer,allocatable :: chrg
  end type lwoniom_frag_placeholder
  type :: lwoniom_input
    integer :: id = 0
  end type lwoniom_input
  type :: lwoniom_data
    integer :: id = 0
    integer :: calcids(2,2)
    integer :: nfrag = 0
    integer :: ncalcs = 0
    type(lwoniom_frag_placeholder),allocatable :: fragment(:)
  end type lwoniom_data
  type :: lwoniom_job
    integer :: id = 0
    integer :: fragid = 0
    integer :: level = 0
    integer :: theoryid = 0
    integer :: nat = 0
    integer :: nlink = 0
    integer,allocatable :: chrg
    integer,allocatable :: uhf
  end type lwoniom_job
  integer,parameter :: oniom_high = 1
  integer,parameter :: oniom_low = 2
  integer,parameter :: oniom_root = 3
#endif

  !> if compiled without(!!!) -DWITH_LWONIOM=true this will export
  !> the placeholders from above. Otherwise it will RE-export
  !> the types and driver layer from lwoniom_interface
  public :: lwoniom_input,lwoniom_data
  public :: lwoniom_job,lwoniom_get_jobs
  public :: oniom_high,oniom_low,oniom_root
#ifdef WITH_LWONIOM
  public :: lwoniom_calculator
  public :: lwoniom_engrad_driver,lwoniom_numhess_driver
  public :: lwoniom_initialize,lwoniom_get_jobgeo
#endif

  public :: ONIOM_read_toml
  public :: ONIOM_get_fraggrad
  public :: ONIOM_compile_error

!========================================================================================!
!========================================================================================!
contains  !> MODULE PROCEDURES START HERE
!========================================================================================!
!========================================================================================!

  subroutine ONIOM_read_toml(tomlfile,nat,at,xyz,ONIOM_data)
!********************************************
!* Read the [lwoniom] block from a toml file
!********************************************
    implicit none
    !> INPUT
    character(len=*),intent(in) :: tomlfile
    integer,intent(in) :: nat
    integer,intent(in) :: at(:)
    real(wp),intent(in) :: xyz(:,:)
    type(lwoniom_data),intent(out) :: ONIOM_data
    type(lwoniom_input),allocatable :: ONIOM_input
#ifdef WITH_LWONIOM
    allocate (ONIOM_input)
    call lwoniom_parse_inputfile(tomlfile,ONIOM_input,required=.false.,natoms=nat)
    ONIOM_input%at = at
    ONIOM_input%xyz = xyz*autoaa !> ONIOM_input needs to store coords in Angstroem rather than bohr
    call lwoniom_new_calculator(ONIOM_input,ONIOM_data) !> because this converts to Bohr
    deallocate (ONIOM_input)
    call ONIOM_data%dump_fragments()
#else
    call ONIOM_compile_error()
#endif
  end subroutine ONIOM_read_toml

!========================================================================================!

  subroutine ONIOM_get_fraggrad(ONIOM,F,gradient,highlow,energy)
!********************************************
!* get high or low level PROJECTED gradient
!* for fragment F, (optionally) with the
!* entire ONIOM energy.
!********************************************
    implicit none
    !> INPUT
    type(lwoniom_data),intent(inout) :: ONIOM
    integer,intent(in) :: F
    real(wp),intent(out) :: gradient(:,:)
    integer,intent(in) :: highlow
    real(wp),intent(out),optional :: energy
    integer :: natf,root_id
    gradient = 0.0_wp
#ifdef WITH_LWONIOM
    if (F > ONIOM%nfrag) error stop 'ONIOM fragment mismatch'
    select case (highlow)
    case (oniom_high)
      gradient = ONIOM%fragment(F)%gradient_high
    case default
      gradient = ONIOM%fragment(F)%gradient_low
    end select
    if (present(energy)) then
      root_id = ONIOM%root_id
      energy = ONIOM%fragment(root_id)%energy_qq
    end if
#else
    call ONIOM_compile_error()
#endif
  end subroutine ONIOM_get_fraggrad

!========================================================================================!

#ifndef WITH_LWONIOM
  subroutine lwoniom_get_jobs(dat,jobs)
!*******************************************************
!* placeholder for the lwONIOM job table setup routine
!*******************************************************
    implicit none
    type(lwoniom_data),intent(inout) :: dat
    type(lwoniom_job),allocatable,intent(out) :: jobs(:)
    call ONIOM_compile_error()
  end subroutine lwoniom_get_jobs
#endif

!========================================================================================!

  subroutine ONIOM_compile_error()
    write (stdout,*) 'Error: Compiled without lwONIOM support!'
    write (stdout,*) 'Use -DWITH_LWONIOM=true in the setup to enable this function'
    error stop
  end subroutine ONIOM_compile_error

!========================================================================================!
!========================================================================================!
end module lwoniom_module
