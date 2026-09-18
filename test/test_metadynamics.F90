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

!> Unit tests for the metadynamics (RMSD bias) potential evaluation.
!> These tests exercise calc_mtd() *directly* (no MD propagation): they set up
!> a static mtdpot with a few stored reference structures and check
!>   (1) the bias energy against a recorded baseline, and
!>   (2) the analytical Cartesian gradient against central finite differences
!>       of the same energy.
!> The finite-difference check is the regression guard that protects the bias
!> forces when the underlying RMSD routine is swapped out (e.g. flat ls_rmsd ->
!> coord-based irmsd wrapper).

module test_metadynamics
  use testdrive,only:new_unittest,unittest_type,error_type,check,test_failed
  use crest_parameters
  use strucrd
  use crest_testmol
  use metadynamics_module
  implicit none
  private

  public :: collect_metadynamics

  !> finite-difference step (Bohr) and gradient comparison tolerance
  real(wp),parameter :: fdstep = 1.0e-4_wp
  real(wp),parameter :: gthr = 1.0e-6_wp

!========================================================================================!
!========================================================================================!
contains  !> Unit tests for the metadynamics RMSD bias potential
!========================================================================================!
!========================================================================================!

!> Collect all exported unit tests
  subroutine collect_metadynamics(testsuite)
    type(unittest_type),allocatable,intent(out) :: testsuite(:)

!&<
    testsuite = [ &
    new_unittest("RMSD-MTD potential (all atoms)  ",test_mtd_rmsd_allatoms), &
    new_unittest("RMSD-MTD gradient (FD, all)     ",test_mtd_rmsd_grad_all), &
    new_unittest("RMSD-MTD gradient (FD, subset)  ",test_mtd_rmsd_grad_subset), &
    new_unittest("RMSD-MTD gradient (FD, PBC/MIC) ",test_mtd_rmsd_grad_pbc), &
    new_unittest("RMSD-MTD large box == gas phase ",test_mtd_rmsd_pbc_largebox_eq_gas), &
    new_unittest("RMSD-MTD minimum-image wrapping ",test_mtd_rmsd_pbc_wrap) &
    ]
!&>
  end subroutine collect_metadynamics

!========================================================================================!
!> Build a deterministic RMSD-bias potential around the given molecule.
!> The stored reference structures are fixed perturbations of mol so the bias
!> is reproducible and free of any RNG / MD state.
  subroutine setup_rmsd_pot(mol,pot,atinclude)
    type(coord),intent(in) :: mol
    type(mtdpot),intent(out) :: pot
    logical,intent(in),optional :: atinclude
    integer :: nref

    nref = 2
    pot%mtdtype = cv_rmsd
    pot%kpush = 0.05_wp
    pot%alpha = 0.60_wp
    pot%maxsave = nref
    pot%ncur = nref
    !> deterministic, "no damping" regime: ramp*cvdump large -> damp ~ 1
    pot%ramp = 0.030_wp
    pot%cvdump = 100000

    allocate (pot%cvxyz(3,mol%nat,nref),source=0.0_wp)
    !> reference 1: displace a couple of atoms
    pot%cvxyz(:,:,1) = mol%xyz
    pot%cvxyz(:,1,1) = pot%cvxyz(:,1,1)+[0.30_wp,-0.20_wp,0.10_wp]
    pot%cvxyz(:,min(3,mol%nat),1) = pot%cvxyz(:,min(3,mol%nat),1)+[0.00_wp,0.15_wp,-0.25_wp]
    !> reference 2: different displacements + a rigid translation (Kabsch must remove it)
    pot%cvxyz(:,:,2) = mol%xyz
    pot%cvxyz(:,2,2) = pot%cvxyz(:,2,2)+[-0.20_wp,0.10_wp,0.05_wp]
    pot%cvxyz(:,mol%nat,2) = pot%cvxyz(:,mol%nat,2)+[0.10_wp,0.10_wp,0.10_wp]
    pot%cvxyz(1,:,2) = pot%cvxyz(1,:,2)+0.05_wp   !> global x-shift

    if (present(atinclude)) then
      if (atinclude) then
        allocate (pot%atinclude(mol%nat),source=.false.)
        !> include roughly half the atoms (every second one)
        pot%atinclude(1:mol%nat:2) = .true.
      end if
    end if
  end subroutine setup_rmsd_pot

!========================================================================================!
!> Energy regression: the bias energy of the (perturbed) molecule against its
!> two stored references must match a recorded baseline value.
  subroutine test_mtd_rmsd_allatoms(error)
    type(error_type),allocatable,intent(out) :: error
    type(coord) :: mol
    type(mtdpot) :: pot
    real(wp),allocatable :: grd(:,:)
    real(wp) :: emtd
    !> baseline recorded from the reference implementation (quaternion RMSD)
    real(wp),parameter :: emtd_ref = 0.09934291441928428_wp

    call get_testmol('cytosine',mol)
    !> perturb the probe geometry so RMSD>0 to both references
    mol%xyz(:,1) = mol%xyz(:,1)+[0.12_wp,-0.07_wp,0.05_wp]

    call setup_rmsd_pot(mol,pot)
    allocate (grd(3,mol%nat),source=0.0_wp)
    call calc_mtd(mol,pot,emtd,grd)

    !> the bias must be positive (repulsive) and finite
    if (.not.(emtd > 0.0_wp .and. emtd < huge(1.0_wp))) then
      call test_failed(error,'RMSD-MTD bias energy is not a positive finite number')
      return
    end if
    !> regression against recorded baseline (see emtd_ref)
    call check(error,emtd,emtd_ref,thr=1.0e-8_wp)
  end subroutine test_mtd_rmsd_allatoms

!========================================================================================!
!> Gradient consistency (all atoms): analytical bias gradient vs central finite
!> differences of the bias energy. This guards the forces across RMSD-backend
!> changes independently of any recorded value.
  subroutine test_mtd_rmsd_grad_all(error)
    type(error_type),allocatable,intent(out) :: error
    type(coord) :: mol
    type(mtdpot) :: pot

    call get_testmol('cytosine',mol)
    mol%xyz(:,1) = mol%xyz(:,1)+[0.12_wp,-0.07_wp,0.05_wp]
    call setup_rmsd_pot(mol,pot)
    call check_mtd_gradient(error,mol,pot)
  end subroutine test_mtd_rmsd_grad_all

!========================================================================================!
!> Gradient consistency (atom subset): same check but with pot%atinclude set,
!> exercising the selected-atom branch of calc_rmsd_mtd.
  subroutine test_mtd_rmsd_grad_subset(error)
    type(error_type),allocatable,intent(out) :: error
    type(coord) :: mol
    type(mtdpot) :: pot

    call get_testmol('cytosine',mol)
    mol%xyz(:,1) = mol%xyz(:,1)+[0.12_wp,-0.07_wp,0.05_wp]
    call setup_rmsd_pot(mol,pot,atinclude=.true.)
    call check_mtd_gradient(error,mol,pot)
  end subroutine test_mtd_rmsd_grad_subset

!========================================================================================!
!> Gradient consistency under PBC: with a lattice present the probe is
!> minimum-image unwrapped onto each reference and the usual Kabsch fit
!> (translation + rotation removed) is applied. The analytical gradient must
!> reproduce central finite differences of the bias energy.
  subroutine test_mtd_rmsd_grad_pbc(error)
    type(error_type),allocatable,intent(out) :: error
    type(coord) :: mol
    type(mtdpot) :: pot

    call get_testmol('cytosine',mol)
    mol%xyz(:,1) = mol%xyz(:,1)+[0.12_wp,-0.07_wp,0.05_wp]
    call set_big_cell(mol)
    call setup_rmsd_pot(mol,pot)
    call check_mtd_gradient(error,mol,pot)
  end subroutine test_mtd_rmsd_grad_pbc

!========================================================================================!
!> Lab-frame equivalence: a flexible molecule in a box much larger than itself
!> must give the SAME bias as the gas-phase (non-periodic) calculation — i.e.
!> rigid-body rotation is still removed, so the bias drives conformational
!> change and not mere reorientation of the molecule within the box.
  subroutine test_mtd_rmsd_pbc_largebox_eq_gas(error)
    type(error_type),allocatable,intent(out) :: error
    type(coord) :: mol
    type(mtdpot) :: pot
    real(wp),allocatable :: grd(:,:)
    real(wp) :: e_gas,e_pbc

    call get_testmol('cytosine',mol)
    mol%xyz(:,1) = mol%xyz(:,1)+[0.12_wp,-0.07_wp,0.05_wp]
    call setup_rmsd_pot(mol,pot)
    allocate (grd(3,mol%nat),source=0.0_wp)

    !> gas-phase bias (no lattice)
    call calc_mtd(mol,pot,e_gas,grd)

    !> same molecule + references, now in a large periodic box
    call set_big_cell(mol)
    call calc_mtd(mol,pot,e_pbc,grd)

    call check(error,e_pbc,e_gas,thr=1.0e-10_wp)
  end subroutine test_mtd_rmsd_pbc_largebox_eq_gas

!========================================================================================!
!> Minimum-image invariance: shifting a reference structure rigidly by an exact
!> lattice vector must not change the bias (each atom's nearest image is
!> unchanged). Compares the bias for an unshifted reference (== mol) with one
!> shifted by one full cell vector.
  subroutine test_mtd_rmsd_pbc_wrap(error)
    type(error_type),allocatable,intent(out) :: error
    type(coord) :: mol
    type(mtdpot) :: pot0,pot1
    real(wp),allocatable :: grd(:,:)
    real(wp) :: e0,e1
    real(wp),parameter :: cell = 50.0_wp

    call get_testmol('cytosine',mol)
    call set_big_cell(mol)
    allocate (grd(3,mol%nat),source=0.0_wp)

    !> reference 0: identical to mol
    call setup_single_ref(mol,pot0,shiftvec=[0.0_wp,0.0_wp,0.0_wp])
    call calc_mtd(mol,pot0,e0,grd)

    !> reference 1: same, but rigidly translated by one lattice vector (+a1)
    call setup_single_ref(mol,pot1,shiftvec=[cell,0.0_wp,0.0_wp])
    call calc_mtd(mol,pot1,e1,grd)

    !> MIC must map the shifted reference back -> identical bias
    call check(error,e1,e0,thr=1.0e-9_wp)
  end subroutine test_mtd_rmsd_pbc_wrap

!========================================================================================!
!> attach a large cubic lattice (Bohr) to a molecule, marking it periodic
  subroutine set_big_cell(mol)
    type(coord),intent(inout) :: mol
    real(wp),parameter :: cell = 50.0_wp
    if (allocated(mol%lat)) deallocate (mol%lat)
    allocate (mol%lat(3,3),source=0.0_wp)
    mol%lat(1,1) = cell
    mol%lat(2,2) = cell
    mol%lat(3,3) = cell
  end subroutine set_big_cell

!========================================================================================!
!> single-reference RMSD pot: the lone reference is mol rigidly shifted by
!> shiftvec (used to probe minimum-image behavior)
  subroutine setup_single_ref(mol,pot,shiftvec)
    type(coord),intent(in) :: mol
    type(mtdpot),intent(out) :: pot
    real(wp),intent(in) :: shiftvec(3)
    integer :: j
    pot%mtdtype = cv_rmsd
    pot%kpush = 0.05_wp
    pot%alpha = 0.60_wp
    pot%maxsave = 1
    pot%ncur = 1
    pot%ramp = 0.030_wp
    pot%cvdump = 100000
    allocate (pot%cvxyz(3,mol%nat,1),source=0.0_wp)
    do j = 1,mol%nat
      pot%cvxyz(:,j,1) = mol%xyz(:,j)+shiftvec
    end do
  end subroutine setup_single_ref

!========================================================================================!
!> Shared helper: compare calc_mtd's analytical gradient to central finite
!> differences of calc_mtd's energy for every Cartesian degree of freedom.
  subroutine check_mtd_gradient(error,mol,pot)
    type(error_type),allocatable,intent(out) :: error
    type(coord),intent(in) :: mol
    type(mtdpot),intent(inout) :: pot
    type(coord) :: tmp
    real(wp),allocatable :: grd(:,:)
    real(wp) :: emtd,ep,em,gfd
    real(wp),allocatable :: gdum(:,:)
    integer :: i,j

    allocate (grd(3,mol%nat),gdum(3,mol%nat),source=0.0_wp)
    call calc_mtd(mol,pot,emtd,grd)

    tmp = mol
    do i = 1,mol%nat
      do j = 1,3
        tmp%xyz(j,i) = mol%xyz(j,i)+fdstep
        call calc_mtd(tmp,pot,ep,gdum)
        tmp%xyz(j,i) = mol%xyz(j,i)-fdstep
        call calc_mtd(tmp,pot,em,gdum)
        tmp%xyz(j,i) = mol%xyz(j,i)   !> restore
        gfd = (ep-em)/(2.0_wp*fdstep)
        call check(error,grd(j,i),gfd,thr=gthr)
        if (allocated(error)) return
      end do
    end do
  end subroutine check_mtd_gradient

!========================================================================================!
!========================================================================================!
end module test_metadynamics
