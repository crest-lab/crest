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

!> Unit tests for the periodic (fixed-cell) CREGEN machinery:
!>   - coord%cellvol() lattice determinant
!>   - the alignment-free SVD fingerprint (pbc_fingerprint): invariance under
!>     atom permutation, global rotation/translation, and minimum-image wrapping
!>   - pbc_identical() duplicate vs. distinct decisions
!>   - minimum-image-aware coordination number (periodic topology)
!>   - the majority-vote periodic/molecular branch detection
!> None of these touch alignment/superposition, mirroring the periodic CREGEN
!> branch they protect.

module test_pbc_cregen
  use testdrive,only:new_unittest,unittest_type,error_type,check,test_failed
  use crest_parameters
  use strucrd
  use crest_testmol
  use pbc_fingerprint_module
  use crest_cn_module,only:calculate_CN
  use cregen_subroutines,only:cregen_majority_periodic
  implicit none
  private

  public :: collect_pbc_cregen

  real(wp),parameter :: cell = 50.0_wp  !> big cubic box (Bohr), no wrapping

!========================================================================================!
!========================================================================================!
contains  !> Unit tests for periodic CREGEN
!========================================================================================!
!========================================================================================!

  subroutine collect_pbc_cregen(testsuite)
    type(unittest_type),allocatable,intent(out) :: testsuite(:)
!&<
    testsuite = [ &
    new_unittest("cellvol cubic + triclinic       ",test_cellvol), &
    new_unittest("fingerprint permutation invar.  ",test_fp_permutation), &
    new_unittest("fingerprint rotation/translation",test_fp_rototrans), &
    new_unittest("fingerprint minimum-image wrap  ",test_fp_mic), &
    new_unittest("pbc_identical duplicate         ",test_identical_dup), &
    new_unittest("pbc_identical distinct          ",test_identical_distinct), &
    new_unittest("periodic CN across cell face    ",test_periodic_cn), &
    new_unittest("branch detection majority vote  ",test_branch_detect) &
    ]
!&>
  end subroutine collect_pbc_cregen

!========================================================================================!
!> cubic box of edge a → a³ ; upper-triangular triclinic → product of diagonal
  subroutine test_cellvol(error)
    type(error_type),allocatable,intent(out) :: error
    type(coord) :: mol
    real(wp),parameter :: a = 4.0_wp

    call get_testmol('cytosine',mol)
    if (allocated(mol%lat)) deallocate (mol%lat)
    allocate (mol%lat(3,3),source=0.0_wp)
    mol%lat(1,1) = a; mol%lat(2,2) = a; mol%lat(3,3) = a
    call check(error,mol%cellvol(),a**3,thr=1.0e-10_wp)
    if (allocated(error)) return

!>--- triclinic (columns = lattice vectors): det = 2*3*5 = 30
    mol%lat = 0.0_wp
    mol%lat(1,1) = 2.0_wp
    mol%lat(1,2) = 1.0_wp; mol%lat(2,2) = 3.0_wp
    mol%lat(2,3) = 1.0_wp; mol%lat(3,3) = 5.0_wp
    call check(error,mol%cellvol(),30.0_wp,thr=1.0e-10_wp)
  end subroutine test_cellvol

!========================================================================================!
!> permuting atoms permutes the columns of D → singular values unchanged
  subroutine test_fp_permutation(error)
    type(error_type),allocatable,intent(out) :: error
    type(coord) :: mol,molp
    real(wp),allocatable :: sig(:),sigp(:)
    integer :: i,nat,p

    call make_pbc_mol(mol)
    nat = mol%nat
    molp = mol
!>--- cyclic shift of the atom order (a non-trivial permutation)
    do i = 1,nat
      p = mod(i,nat)+1
      molp%at(p) = mol%at(i)
      molp%xyz(:,p) = mol%xyz(:,i)
    end do
    call pbc_fingerprint(mol,sig)
    call pbc_fingerprint(molp,sigp)
    call check(error,fp_distance(sig,sigp),0.0_wp,thr=1.0e-8_wp)
  end subroutine test_fp_permutation

!========================================================================================!
!> rigid rotation + translation leaves all distances unchanged → same fingerprint
  subroutine test_fp_rototrans(error)
    type(error_type),allocatable,intent(out) :: error
    type(coord) :: mol,molr
    real(wp),allocatable :: sig(:),sigr(:)
    real(wp) :: R(3,3),ang,shift(3)
    integer :: i

    call make_pbc_mol(mol)
    molr = mol
    ang = 0.7_wp
    R = 0.0_wp
    R(1,1) = cos(ang); R(1,2) = -sin(ang)
    R(2,1) = sin(ang); R(2,2) = cos(ang)
    R(3,3) = 1.0_wp
    shift = [3.1_wp,-2.4_wp,1.7_wp]
    do i = 1,mol%nat
      molr%xyz(:,i) = matmul(R,mol%xyz(:,i))+shift
    end do
    call pbc_fingerprint(mol,sig)
    call pbc_fingerprint(molr,sigr)
    call check(error,fp_distance(sig,sigr),0.0_wp,thr=1.0e-8_wp)
  end subroutine test_fp_rototrans

!========================================================================================!
!> shifting one atom by a full lattice vector is a no-op under the MIC
  subroutine test_fp_mic(error)
    type(error_type),allocatable,intent(out) :: error
    type(coord) :: mol,molw
    real(wp),allocatable :: sig(:),sigw(:)

    call make_pbc_mol(mol)
    molw = mol
    molw%xyz(:,1) = mol%xyz(:,1)+mol%lat(:,1)  !> wrap atom 1 by +a1
    call pbc_fingerprint(mol,sig)
    call pbc_fingerprint(molw,sigw)
    call check(error,fp_distance(sig,sigw),0.0_wp,thr=1.0e-9_wp)
  end subroutine test_fp_mic

!========================================================================================!
!> a permuted + rotated + wrapped copy at the same energy/volume is a duplicate
  subroutine test_identical_dup(error)
    type(error_type),allocatable,intent(out) :: error
    type(coord) :: mol,molc
    real(wp),allocatable :: sig(:),sigc(:)
    real(wp) :: R(3,3),ang
    integer :: i,nat,p
    real(wp),parameter :: ethr = 1.0e-4_wp

    call make_pbc_mol(mol)
    mol%energy = -42.0_wp
    nat = mol%nat
    molc = mol
    molc%energy = mol%energy+0.5e-4_wp   !> within ethr
!>--- rotate
    ang = -0.4_wp
    R = 0.0_wp
    R(1,1) = cos(ang); R(1,2) = -sin(ang)
    R(2,1) = sin(ang); R(2,2) = cos(ang)
    R(3,3) = 1.0_wp
    do i = 1,nat
      molc%xyz(:,i) = matmul(R,mol%xyz(:,i))
    end do
!>--- and permute (cyclic)
    block
      type(coord) :: tmp
      tmp = molc
      do i = 1,nat
        p = mod(i,nat)+1
        tmp%at(p) = molc%at(i)
        tmp%xyz(:,p) = molc%xyz(:,i)
      end do
      molc = tmp
    end block
    call pbc_fingerprint(mol,sig)
    call pbc_fingerprint(molc,sigc)
    call check(error,pbc_identical(mol,molc,sig,sigc,ethr))
  end subroutine test_identical_dup

!========================================================================================!
!> a genuinely different geometry (one atom displaced) is NOT a duplicate
  subroutine test_identical_distinct(error)
    type(error_type),allocatable,intent(out) :: error
    type(coord) :: mol,mold
    real(wp),allocatable :: sig(:),sigd(:)
    real(wp),parameter :: ethr = 1.0e-4_wp

    call make_pbc_mol(mol)
    mol%energy = -42.0_wp
    mold = mol            !> same energy and same cell volume ...
    mold%xyz(1,1) = mol%xyz(1,1)+1.0_wp  !> ... but a 1 Bohr displacement
    call pbc_fingerprint(mol,sig)
    call pbc_fingerprint(mold,sigd)
!>--- fingerprint must differ beyond the threshold
    call check(error,fp_distance(sig,sigd) > sigthr)
    if (allocated(error)) return
    call check(error,.not.pbc_identical(mol,mold,sig,sigd,ethr))
  end subroutine test_identical_distinct

!========================================================================================!
!> minimum-image CN bonds two atoms separated across a cell face; the plain
!> (non-periodic) CN does not.
  subroutine test_periodic_cn(error)
    type(error_type),allocatable,intent(out) :: error
    integer,parameter :: nat = 2
    integer :: at(nat)
    real(wp) :: xyz(3,nat),lat(3,3),cn(nat)
    real(wp),parameter :: box = 8.0_wp,bl = 2.85_wp

    at = [6,6]                       !> C-C
    xyz = 0.0_wp
    xyz(1,1) = 0.0_wp
    xyz(1,2) = box-bl                !> direct dist = box-bl (far); MIC = bl (bond)
    lat = 0.0_wp
    lat(1,1) = box; lat(2,2) = box; lat(3,3) = box

!>--- without lattice: atoms are far apart, no bond
    call calculate_CN(nat,at,xyz,cn)
    call check(error,cn(1) < 0.1_wp)
    if (allocated(error)) return

!>--- with lattice (MIC): atoms are bonded across the face
    call calculate_CN(nat,at,xyz,cn,lat=lat)
    call check(error,cn(1) > 0.5_wp)
  end subroutine test_periodic_cn

!========================================================================================!
!> majority of structures with a (non-zero) cell → periodic; otherwise molecular
  subroutine test_branch_detect(error)
    type(error_type),allocatable,intent(out) :: error
    type(coord),allocatable :: ens(:)
    integer :: i

    allocate (ens(3))
    do i = 1,3
      call get_testmol('cytosine',ens(i))
    end do
!>--- 2 of 3 carry a cell → periodic
    call add_cell(ens(1))
    call add_cell(ens(2))
    call check(error,cregen_majority_periodic(ens))
    if (allocated(error)) return

!>--- only 1 of 3 carries a cell → molecular
    if (allocated(ens(2)%lat)) deallocate (ens(2)%lat)
    call check(error,.not.cregen_majority_periodic(ens))
  end subroutine test_branch_detect

!========================================================================================!
!> helpers
!========================================================================================!

  subroutine make_pbc_mol(mol)
    type(coord),intent(out) :: mol
    call get_testmol('cytosine',mol)
    call add_cell(mol)
  end subroutine make_pbc_mol

  subroutine add_cell(mol)
    type(coord),intent(inout) :: mol
    if (allocated(mol%lat)) deallocate (mol%lat)
    allocate (mol%lat(3,3),source=0.0_wp)
    mol%lat(1,1) = cell
    mol%lat(2,2) = cell
    mol%lat(3,3) = cell
  end subroutine add_cell

!========================================================================================!
!========================================================================================!
end module test_pbc_cregen
