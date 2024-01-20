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
  type :: lwoniom_frag_placeholder
    integer,allocatable :: chrg 
  end  type lwoniom_frag_placeholder
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
    ONIOM_input%xyz = xyz*autoaa !> ONIOM_input needs to store coords in Angstroem rather than bohr
    call lwoniom_new_calculator( ONIOM_input, ONIOM_data ) !> because this converts to Bohr
    deallocate(ONIOM_input)
    call ONIOM_data%dump_fragments()
#else
    call ONIOM_compile_error()
#endif
  end subroutine ONIOM_read_toml


!========================================================================================!


  subroutine ONIOM_update_geo(ONIOM,mol,mollist,cmap)
!**********************************************************
!* Update all fragments of a previously set up ONIOM setup
!**********************************************************
    implicit none
    !> INPUT
    type(lwoniom_data),intent(inout) :: ONIOM
    type(coord),intent(in)  :: mol
    type(coord),intent(inout),optional :: mollist(ONIOM%ncalcs)
    integer,intent(in),optional :: cmap(ONIOM%ncalcs) !> calculation index mapping
    integer :: i,highlow,fragid,j
#ifdef WITH_LWONIOM
    call ONIOM%update( mol%xyz )  !> no point charge version

    !> optional, update a given list of coord-type molecules (the fragments)
    if(present(mollist).and.present(cmap))then
      do i=1,ONIOM%ncalcs
          j = cmap(i)
          call ONIOM_get_mapping(j, ONIOM%calcids, highlow, fragid)
          call ONIOM_get_mol(ONIOM,fragid,mollist(i),highlow)
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
    integer :: i,l
#ifdef WITH_LWONIOM
    call lwoniom_singlepoint(mol%nat,ONIOM,energy,gradient)
#else
    call ONIOM_compile_error()
#endif
  end subroutine ONIOM_engrad

!========================================================================================!

  subroutine ONIOM_get_mol(ONIOM,F,mol,highlow)
!*****************************************
!* transfer data from fragment F into mol
!*****************************************
    implicit none
    !> INPUT
    type(lwoniom_data),intent(inout) :: ONIOM
    type(coord),intent(inout)  :: mol
    integer,intent(in) :: F 
    integer,intent(in),optional :: highlow
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
    integer,intent(in) :: highlow
    real(wp),intent(out),optional :: energy
    integer :: natf,root_id
    gradient = 0.0_wp 
#ifdef WITH_LWONIOM
    if(F > ONIOM%nfrag ) error stop 'ONIOM fragment mismatch'
    select case(highlow)
    case(1)
       gradient = ONIOM%fragment(F)%gradient_high
    case default
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

  subroutine ONIOM_get_mapping(id, calcids, highlow, fragid)
     implicit none
     integer,intent(in) :: id !> calculation ID in CREST
     integer,intent(in) :: calcids(:,:) !> mappings form calculator setup
     integer,intent(out) :: highlow,fragid
     integer :: i,j,k,l
     highlow = 0
     fragid = 0
     iloop : do i=1,size(calcids,2)
       do j=1,2
         if( calcids(j,i) == id)then
           highlow = j
           fragid = i
           exit iloop
         endif
       enddo
     enddo iloop
   end subroutine ONIOM_get_mapping

!========================================================================================!

  subroutine ONIOM_compile_error()
    write (stdout,*) 'Error: Compiled without lwONIOM support!'
    write (stdout,*) 'Use -DWITH_LWONIOM=true in the setup to enable this function'
    error stop
  end subroutine ONIOM_compile_error

!========================================================================================!
!========================================================================================!
end module lwoniom_module

