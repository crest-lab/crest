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
  use crest_parameters
  use strucrd
#ifdef WITH_LWONIOM
  use lwoniom_interface
#endif
  implicit none
  private

#ifndef WITH_LWONIOM
  !> this is a placeholder if no lwONIOM module is used!
  type :: lwoniom_input
    integer :: id = 0
  end type lwoniom_input
  type :: lwoniom_data
    integer :: id = 0
    integer :: calcids(2,2)
    integer :: nfrag = 0
  end type lwoniom_data
#endif

  !> if compiled without(!!!) -DWITH_LWONIOM=true this will export
  !> the placeholder from above. Otherwise it will RE-export
  !> the types from lwoniom_interface
  public :: lwoniom_input,lwoniom_data

  public :: ONIOM_read_toml
 
  public :: ONIOM_update_geo

  public :: ONIOM_associate_mol

  public :: ONIOM_get_fraggrad
  
  public :: ONIOM_engrad

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
    allocate(ONIOM_input)
    call lwoniom_parse_inputfile(tomlfile,ONIOM_input,required = .false., natoms=nat)
    ONIOM_input%at = at
    ONIOM_input%xyz = xyz
    call lwoniom_new_calculator( ONIOM_input, ONIOM_data )
    deallocate(ONIOM_input)
#else
    call ONIOM_compile_error()
#endif
  end subroutine ONIOM_read_toml


!========================================================================================!


  subroutine ONIOM_update_geo(ONIOM,mol,mollist)
!**********************************************************
!* Update all fragments of a previously set up ONIOM setup
!**********************************************************
    implicit none
    !> INPUT
    type(lwoniom_data),intent(inout) :: ONIOM
    type(coord),intent(in)  :: mol
    type(coord),intent(inout),optional :: mollist(ONIOM%nfrag)
    integer :: i
#ifdef WITH_LWONIOM
    call ONIOM%update( mol%xyz )  !> no point charge version

    !> optional, update a given list of coord-type molecules (the fragments)
    if(present(mollist))then
      do i=1,ONIOM%nfrag
       call ONIOM_get_mol(ONIOM,i,mollist(i))
      enddo
    endif
#else
    call ONIOM_compile_error()
#endif
  end subroutine ONIOM_update_geo


!========================================================================================!

  subroutine ONIOM_engrad(ONIOM,mol,energy,gradient)
!**********************************************
!* reconstruct ONIOM energy and gradient
!* energy and gradient is additive to whatever
!* is saved there already
!**********************************************
    implicit none
    !> INPUT
    type(lwoniom_data),intent(inout) :: ONIOM
    type(coord),intent(in)  :: mol

    real(wp),intent(inout) :: energy
    real(wp),intent(inout) :: gradient(3,mol%nat)
    integer :: i
#ifdef WITH_LWONIOM
    call lwoniom_singlepoint(mol%nat,ONIOM,energy,gradient)
 
#else
    call ONIOM_compile_error()
#endif
  end subroutine ONIOM_engrad




!========================================================================================!

  subroutine ONIOM_get_mol(ONIOM,F,mol)
!*****************************************
!* transfer data from fragment F into mol
!*****************************************
    implicit none
    !> INPUT
    type(lwoniom_data),intent(inout) :: ONIOM
    type(coord),intent(inout)  :: mol
    integer,intent(in) :: F 
    integer :: natf   
#ifdef WITH_LWONIOM
    if(F > ONIOM%nfrag ) error stop 'ONIOM fragment mismatch'

    natf =  ONIOM%fragment(F)%nat +  ONIOM%fragment(F)%nlink
    if(mol%nat /=  natf) call mol%deallocate()

    mol%nat = natf
    mol%at = reshape( [ONIOM%fragment(F)%at,ONIOM%fragment(F)%linkat], [natf])
    mol%xyz = reshape( [ONIOM%fragment(F)%xyz,ONIOM%fragment(F)%linkxyz], [3,natf]) 

#else
    call ONIOM_compile_error()
#endif
  end subroutine ONIOM_get_mol

!=========================================================================================!

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
    character(len=*),intent(in) :: highlow
    real(wp),intent(out),optional :: energy
    integer :: natf,root_id
#ifdef WITH_LWONIOM
    if(F > ONIOM%nfrag ) error stop 'ONIOM fragment mismatch'
    select case(highlow)
    case("high")
       gradient = ONIOM%fragment(F)%gradient_high
    case("low")
       gradient = ONIOM%fragment(F)%gradient_low
    end select
    if(present(energy))then
       root_id = ONIOM%root_id
       energy = ONIOM%fragment(root_id)%energy_qq
    endif
#else
    call ONIOM_compile_error()
#endif
  end subroutine ONIOM_get_fraggrad

!========================================================================================!

  subroutine ONIOM_associate_mol(mol,molptr)
    implicit none
    type(coord),target :: mol
    type(coord),pointer,intent(out) :: molptr
    molptr => mol
  end subroutine ONIOM_associate_mol

!========================================================================================!

  subroutine ONIOM_compile_error()
    write (stdout,*) 'Error: Compiled without lwONIOM support!'
    write (stdout,*) 'Use -DWITH_LWONIOM=true in the setup to enable this function'
    error stop
  end subroutine ONIOM_compile_error

!========================================================================================!
!========================================================================================!
end module lwoniom_module

