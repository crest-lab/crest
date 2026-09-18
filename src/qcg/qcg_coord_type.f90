!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2025 Philipp Pracht
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
!
module qcg_coord_type
  use crest_parameters,only:wp
  use strucrd,only:coord
  implicit none
  public

  type,extends(coord) :: coord_qcg
    !> new components that are added to the coord type:
    integer   :: nmol          !> number of molecules
    real(wp)  :: cma(3)        !> center of mass
    real(wp)  :: aniso         !> anisotropy factor
    real(wp)  :: ell_abc(3)    !> ellipsoid axis
    real(wp)  :: atot          !> surface area
    real(wp)  :: vtot          !> volume
    real(wp)  :: rtot          !> radius
    real(wp)  :: mass          !> mass
    real(wp)  :: gt            !> gibbs free energy
    real(wp)  :: ht            !> enthalpy
    real(wp)  :: svib          !> vibrational entropy
    real(wp)  :: srot          !> rotational entropy
    real(wp)  :: stra          !> translational entropy
    real(wp)  :: eax(3)        !> molecular axis
    character(len=20) :: elevel = ''  !> level of theory at which %energy was evaluated
  contains
    procedure :: as_coord
    procedure :: from_coord
  end type coord_qcg

!==============================================================================!
contains  !> MODULE PROCEDURES START HERE
!==============================================================================!

  function as_coord(this) result(mol)
    class(coord_qcg),intent(in) :: this
    type(coord) :: mol

    mol%nat = this%nat
    if (allocated(this%at)) mol%at = this%at
    if (allocated(this%xyz)) mol%xyz = this%xyz

    mol%energy = this%energy
    if (allocated(this%comment)) mol%comment = this%comment
    mol%chrg = this%chrg
    mol%uhf = this%uhf
    mol%nbd = this%nbd
    if (allocated(this%bond)) mol%bond = this%bond
    if (allocated(this%lat)) mol%lat = this%lat
    if (allocated(this%qat)) mol%qat = this%qat
    mol%pdb = this%pdb

  end function as_coord

  subroutine from_coord(this,mol)
    class(coord_qcg),intent(inout) :: this
    type(coord),intent(in) :: mol

    this%nat = mol%nat
    if (allocated(mol%at)) this%at = mol%at
    if (allocated(mol%xyz)) this%xyz = mol%xyz

    this%energy = mol%energy
    if (allocated(mol%comment)) this%comment = mol%comment
    this%chrg = mol%chrg
    this%uhf = mol%uhf
    this%nbd = mol%nbd
    if (allocated(mol%bond)) this%bond = mol%bond
    if (allocated(mol%lat)) this%lat = mol%lat
    if (allocated(mol%qat)) this%qat = mol%qat
    this%pdb = mol%pdb
  end subroutine from_coord

end module qcg_coord_type

