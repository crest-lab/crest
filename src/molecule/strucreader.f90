!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2020-2023 Philipp Pracht
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

!> Exports the "coord" type and I/O
module strucrd
  use molecule_type
  use molecule_type_components
  use molecule_type_ensemble
  use molecule_parameters, only: coordtype
  use molecule_io
  implicit none
  private

  !> RE-EXPORTS FROM THE ABOVE MODULES
! ══════════════════════════════════════════════════════════════════════════════
  public :: coord     !> coord type
  public :: ensemble  !> ensemble type (sparsely used, better use a list of coord objects)
  public :: mollist   !> list of coord objects

! ══════════════════════════════════════════════════════════════════════════════
  public :: i2e          !> function to convert atomic number to element symbol
  public :: asym         !> alterinative signature for i2e
  public :: e2i          !> function to convert element symbol into atomic number

  public :: checkcoordtype !> determine input coordinate file type (mostly via extension)
  public :: coordtype      !> Possible return types from checkcoordtype → e.g. coordtype%turbomole

  public :: rdnat        !> procedure to read number of atoms Nat
  public :: rdcoord      !> read an input file, determine format automatically
  public :: rdxmol       !> read a file in the Xmol (.xyz) format specifically
  public :: rdxmolselec  !> read only a certain structure in Xmol file

  !> NOTE, using coord%write() is safer than the ones below
  public :: wrc0     !> write file in turbomole format
  public :: wrcoord  !> write file by name, type via extension
  public :: wrxyz    !> write file in xyz format
  public :: wrsdf    !> write file in sdf format

  public :: xyz2coord
  public :: coord2xyz
  public :: rdensembleparam   !> read Nat and Nall for a XYZ trajectory
  public :: rdensemble        !> read a XYZ trajectory
  public :: wrensemble

  public :: coordline
  public :: grepenergy
  public :: get_atlist
  public :: sumform

! ══════════════════════════════════════════════════════════════════════════════
contains  !> MODULE PROCEDURES START HERE
! ══════════════════════════════════════════════════════════════════════════════

! ══════════════════════════════════════════════════════════════════════════════
end module strucrd
