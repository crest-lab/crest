!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2026 Philipp Pracht
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

module molecule_type_components
  use iso_c_binding
  use molecule_parameters
  implicit none
  private
! ══════════════════════════════════════════════════════════════════════════════

  !coord class. contains a single structure in the PDB format.
  !coordinates by definition are in Angstroem.
  type :: pdbdata
    !--- data
    integer :: nat = 0
    integer :: frag = 0
    !--- arrays
    integer,allocatable  :: athet(:) !ATOM (1) or HETATM (2)
    character(len=4),allocatable :: pdbat(:) !PDB atom specifier
    character(len=3),allocatable :: pdbas(:) !PDB amino acid specifier
    integer,allocatable :: pdbfrag(:) !PDB fragment specifier
    character(len=1),allocatable :: pdbgrp(:)  !PDB group specifier
    real(wp),allocatable :: pdbocc(:) !PDB occupancy
    real(wp),allocatable :: pdbtf(:)  !PDB temperature factor
  contains
    procedure :: deallocate => deallocate_pdb !clear memory space
    procedure :: allocate => allocate_pdb
  end type pdbdata
  public :: pdbdata

! ──────────────────────────────────────────────────────────────────────────────

  public :: signature,extxyz_signatures,extxyz_properties
  public :: parse_properties_tag,assemble_properties_tag
  public :: allocate_extxyz_properties_from_sigs

  ! Type representing a single property entry (e.g., pos:R:3)
  type :: signature
    character(len=32) :: name     ! Property name (e.g., "forces")
    character         :: p_type   ! 'R' for Real, 'S' for String, 'I' for Int
    integer           :: n_fields ! Number of columns (e.g., 3 for positions)
  end type signature

  ! Type representing the collection of all properties in the file
  type :: extxyz_signatures
    type(signature),allocatable :: signat(:)
    integer :: n_props = 0            ! Number of unique property keys
    integer :: total_fields = 0       ! Total sum of all n_fields (total columns)
  end type extxyz_signatures

  type :: extxyz_property
    integer :: natoms = 0
    type(signature) :: signat
    character(len=32),allocatable :: S(:,:)
    integer,allocatable           :: I(:,:)
    real(wp),allocatable          :: R(:,:)
  end type extxyz_property

  type :: extxyz_properties
    integer :: n_props = 0
    integer :: total_fields = 0
    type(extxyz_property),allocatable :: props(:)
  end type extxyz_properties

! ══════════════════════════════════════════════════════════════════════════════
contains  !> MODULE PROCEDURES START HERE
! ══════════════════════════════════════════════════════════════════════════════

!==================================================================!
! subroutine deallocate_pdb
! is used to clear memory for the pdbdata type
!==================================================================!
  subroutine deallocate_pdb(self)
    implicit none
    class(pdbdata) :: self
    self%nat = 0
    self%frag = 0
    if (allocated(self%athet)) deallocate (self%athet)
    if (allocated(self%pdbat)) deallocate (self%pdbat)
    if (allocated(self%pdbas)) deallocate (self%pdbas)
    if (allocated(self%pdbfrag)) deallocate (self%pdbfrag)
    if (allocated(self%pdbgrp)) deallocate (self%pdbgrp)
    if (allocated(self%pdbocc)) deallocate (self%pdbocc)
    if (allocated(self%pdbtf)) deallocate (self%pdbtf)
    return
  end subroutine deallocate_pdb

!==================================================================!
! subroutine allocate_pdb
! is used to clear memory for the pdbdata type
!==================================================================!
  subroutine allocate_pdb(self,nat)
    implicit none
    class(pdbdata) :: self
    integer :: nat
    call deallocate_pdb(self)
    self%nat = nat
    allocate (self%athet(nat))
    allocate (self%pdbat(nat))
    allocate (self%pdbas(nat))
    allocate (self%pdbfrag(nat))
    allocate (self%pdbgrp(nat))
    allocate (self%pdbocc(nat))
    allocate (self%pdbtf(nat))
    return
  end subroutine allocate_pdb

! ──────────────────────────────────────────────────────────────────────────────

  subroutine parse_properties_tag(prop_str,ext_sigs)
!*************************************************************************************
!*   Parses the "Properties" value string from an extXYZ comment line.               *
!*   Following the ASE (Atomic Simulation Environment) standard, it decomposes       *
!*   the string (format: "name:type:cols:name:type:cols...") into a structured       *
!*   array of 'signature' types.                                                     *
!*                                                                                   *
!* ARGUMENTS:                                                                        *
!*   prop_str  [IN]  : The raw string value of the Properties tag.                   *
!*                     Example: "species:S:1:pos:R:3:forces:R:3"                     *
!*   ext_sigs [OUT] : An instance of extxyz_signatures.                             *
!*                     - Allocates the 'props' array based on the number of triplets.*
!*                     - Calculates 'total_fields' for buffer allocation.            *
!*                                                                                   *
!* DATA MAPPING:                                                                     *
!*   - name    : String label of the property.                                       *
!*   - p_type  : 'R' (Real), 'S' (String), 'I' (Integer).                            *
!*   - n_fields: Integer representing the number of columns this property spans.     *
!*                                                                                   *
!* NOTES:                                                                            *
!*   - This routine assumes the input string is a valid series of triplets.          *
!*   - It handles both trailing colons and clean endings.                            *
!*************************************************************************************
    character(len=*),intent(in)        :: prop_str
    type(extxyz_signatures),intent(out) :: ext_sigs

    integer :: i,start_pos,end_pos,part_count,i_prop
    character(len=len_trim(prop_str)) :: buffer

    ! 1. Initial count of colons to determine array size
    ! Format is name:type:cols -> 2 colons per property, +1 at the end of parts
    part_count = 0
    do i = 1,len_trim(prop_str)
      if (prop_str(i:i) == ':') part_count = part_count+1
    end do

    ext_sigs%n_props = (part_count+1)/3
    allocate (ext_sigs%signat(ext_sigs%n_props))
    ext_sigs%total_fields = 0

    ! 2. Parse the triplets
    buffer = trim(prop_str)
    start_pos = 1

    do i_prop = 1,ext_sigs%n_props
      ! Extract Name
      end_pos = index(buffer(start_pos:),':')+start_pos-2
      ext_sigs%signat(i_prop)%name = buffer(start_pos:end_pos)
      start_pos = end_pos+2

      ! Extract Type (R/S/I)
      ext_sigs%signat(i_prop)%p_type = buffer(start_pos:start_pos)
      start_pos = start_pos+2 ! Skip char and following colon

      ! Extract Number of Fields
      end_pos = index(buffer(start_pos:),':')+start_pos-2
      if (end_pos < start_pos) end_pos = len_trim(buffer) ! Handle last element

      read (buffer(start_pos:end_pos),*) ext_sigs%signat(i_prop)%n_fields
      start_pos = end_pos+2

      ! Update global counter
      ext_sigs%total_fields = ext_sigs%total_fields+ext_sigs%signat(i_prop)%n_fields
    end do
  end subroutine parse_properties_tag

  subroutine assemble_properties_tag(ext_sigs,prop_str)
    implicit none
    type(extxyz_signatures),intent(in) :: ext_sigs
    character(len=*),intent(out)       :: prop_str

    integer :: i
    character(len=16) :: col_buffer ! Temporary buffer for integer conversion

    ! Initialize the string as empty
    prop_str = ""

    do i = 1,ext_sigs%n_props
      ! 1. Append the Name
      prop_str = trim(prop_str)//trim(ext_sigs%signat(i)%name)//":"

      ! 2. Append the Type (R/S/I)
      prop_str = trim(prop_str)//ext_sigs%signat(i)%p_type//":"

      ! 3. Append the Number of Columns
      write (col_buffer,'(I0)') ext_sigs%signat(i)%n_fields
      prop_str = trim(prop_str)//trim(col_buffer)

      ! 4. Add a colon separator UNLESS this is the last property
      if (i < ext_sigs%n_props) then
        prop_str = trim(prop_str)//":"
      end if
    end do
  end subroutine assemble_properties_tag

! ──────────────────────────────────────────────────────────────────────────────
  subroutine allocate_extxyz_properties_from_sigs(nat,ext_sigs,ext_props)
    implicit none
    integer,intent(in) :: nat
    type(extxyz_signatures),intent(in) :: ext_sigs
    type(extxyz_properties),intent(out) :: ext_props

    integer :: ii,n_props,n_fields

    n_props = ext_sigs%n_props
    ext_props%n_props = n_props
    ext_props%total_fields = ext_sigs%total_fields
    allocate (ext_props%props(n_props))
    do ii = 1,n_props
      associate (prop => ext_props%props(ii))
        prop%signat = ext_sigs%signat(ii)
        prop%natoms = nat
        n_fields = prop%signat%n_fields
        select case (prop%signat%p_type)
        case ('I')
          allocate (prop%I(n_fields,nat),source=0)
        case ('R')
          allocate (prop%R(n_fields,nat),source=0.0_wp)
        case ('S')
          allocate (prop%S(n_fields,nat),source=repeat(' ',32))
        end select
      end associate
    end do

  end subroutine allocate_extxyz_properties_from_sigs

! ══════════════════════════════════════════════════════════════════════════════
! ══════════════════════════════════════════════════════════════════════════════
end module molecule_type_components
