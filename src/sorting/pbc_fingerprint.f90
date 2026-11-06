!===============================================================================!
! This file is part of crest.
!
! Copyright (C) 2018-2026 Philipp Pracht
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
!===============================================================================!

module pbc_fingerprint_module
!********************************************************************************
!* Alignment-free duplicate detection for periodic structures.
!*
!* Implements the structure-comparison measure of
!*   P. Pracht, J. W. R. Morgan, D. J. Wales,
!*   J. Chem. Phys. 159, 064801 (2023), §II.C, Eq. (8).
!*
!* For a structure under full PBC, the symmetric N×N minimum-image distance
!* matrix D_ij = |r_i - r_j| is built. The singular values Σ of D (= |eigenvalues|
!* since D is symmetric) form a "fingerprint" that is invariant to atom
!* permutation (a simultaneous row/column permutation preserves the spectrum),
!* translation, rotation, and reflection — no superposition is required.
!*
!* NOTE: the fingerprint distance ΔΣ tracks structural similarity but is a much
!* tighter "is this literally the same structure" measure than RMSD; it is NOT a
!* drop-in conformer/rotamer metric. In the periodic CREGEN it is therefore used
!* as a *loose* pre-grouping filter (SIGTHR_GROUP), with the minimum-image RMSD
!* doing the precise duplicate decision. pbc_identical() (strict SIGTHR) remains
!* available as an exact-duplicate test.
!********************************************************************************
  use crest_parameters
  use strucrd,only:coord
  implicit none
  private

  public :: pbc_fingerprint
  public :: fp_distance
  public :: pbc_identical

!>--- thresholds (paper defaults)
!>    Σthr = 1e-5 (dimensionless fingerprint difference): strict, exact-duplicate
!>    test used by pbc_identical()
  real(wp),parameter,public :: sigthr = 1.0e-5_wp
!>    Σthr_group: *loose* fingerprint threshold used as a cheap skip-filter in the
!>    periodic CREGEN — a pair whose fingerprints differ by more than this cannot
!>    be a duplicate, so the MIC-RMSD is skipped. It must sit comfortably ABOVE
!>    the fingerprint distance of any genuine RMSD-duplicate (those are ~1e-2 or
!>    smaller), so it is deliberately generous; the MIC-RMSD (RTHR) makes the
!>    authoritative decision. Coarse, size-dependent, tunable per system.
  real(wp),parameter,public :: sigthr_group = 1.0_wp
!>    Ωthr = 0.5 Å³, converted to Bohr³
  real(wp),parameter,public :: volthr = 0.5_wp*(aatoau**3)

!========================================================================================!
!========================================================================================!
contains  !> MODULE PROCEDURES START HERE
!========================================================================================!
!========================================================================================!

  subroutine pbc_fingerprint(mol,sig)
!********************************************************************
!* Build the symmetric minimum-image distance matrix D_ij=|r_i-r_j|
!* and return its singular values (sig), the alignment-free
!* fingerprint Σ of the structure.
!*
!*   mol - input structure (coord); uses %lat for the MIC if allocated
!*   sig - output singular values, dimension (mol%nat)
!********************************************************************
    use strucrd,only:coord
    implicit none
    type(coord),intent(in) :: mol
    real(wp),allocatable,intent(out) :: sig(:)
    integer :: nat,i,j,info,lwork
    real(wp) :: latinv(3,3),sfrac(3),rij(3)
    real(wp),allocatable :: D(:,:),work(:)
    real(wp) :: udum(1,1),vtdum(1,1),wq(1)
    logical :: periodic

    nat = mol%nat
    allocate (sig(nat),source=0.0_wp)
    if (nat < 1) return

    periodic = allocated(mol%lat)
    if (periodic) latinv = inv3x3_fp(mol%lat)

! ── build the symmetric N×N MIC distance matrix D_ij = |r_i - r_j| ───────────
!    (EXPERIMENT: actual distance matrix, no per-column sorting. Symmetric, so
!     its singular values = |eigenvalues| are an isometry invariant directly.)
    allocate (D(nat,nat),source=0.0_wp)
    do i = 1,nat
      do j = 1,nat
        rij(:) = mol%xyz(:,j)-mol%xyz(:,i)
        if (periodic) then
          sfrac = matmul(latinv,rij)
          sfrac = sfrac-anint(sfrac)
          rij = matmul(mol%lat,sfrac)
        end if
        D(j,i) = sqrt(sum(rij**2))
      end do
    end do

! ── singular values only (JOBU=JOBVT='N') via LAPACK DGESVD ──────────────────
    lwork = -1
    call dgesvd('N','N',nat,nat,D,nat,sig,udum,1,vtdum,1,wq,lwork,info)
    lwork = max(1,nint(wq(1)))
    allocate (work(lwork))
    call dgesvd('N','N',nat,nat,D,nat,sig,udum,1,vtdum,1,work,lwork,info)
    deallocate (work,D)
  end subroutine pbc_fingerprint

!========================================================================================!

  function fp_distance(sigA,sigB) result(dsig)
!********************************************************************
!* Sum of squared differences between two fingerprints,
!* ΔΣ = Σ_i |Σ_{i,A} − Σ_{i,B}|²  (Eq. 8 comparison measure).
!* The two fingerprints must have the same length.
!********************************************************************
    implicit none
    real(wp),intent(in) :: sigA(:),sigB(:)
    real(wp) :: dsig
    integer :: n
    n = min(size(sigA),size(sigB))
    dsig = sum((sigA(1:n)-sigB(1:n))**2)
  end function fp_distance

!========================================================================================!

  function pbc_identical(molA,molB,sigA,sigB,ethr) result(same)
!********************************************************************
!* Decide whether two periodic structures are identical.
!* True iff ALL of the following hold:
!*   |ΔE|  < ethr    (energy, Hartree)
!*   |ΔΩ|  < volthr  (cell volume, Bohr³)
!*   ΔΣ    < sigthr  (fingerprint difference)
!********************************************************************
    use strucrd,only:coord
    implicit none
    type(coord),intent(in) :: molA,molB
    real(wp),intent(in) :: sigA(:),sigB(:)
    real(wp),intent(in) :: ethr
    logical :: same
    real(wp) :: de,dvol,dsig
    same = .false.
    de = abs(molA%energy-molB%energy)
    if (de >= ethr) return
    dvol = abs(molA%cellvol()-molB%cellvol())
    if (dvol >= volthr) return
    dsig = fp_distance(sigA,sigB)
    if (dsig >= sigthr) return
    same = .true.
  end function pbc_identical

!========================================================================================!
!> helpers
!========================================================================================!

  pure function inv3x3_fp(a) result(ainv)
!********************************************************************
!* Inverse of a 3x3 matrix (cofactor/adjugate) for the MIC.
!********************************************************************
    implicit none
    real(wp),intent(in) :: a(3,3)
    real(wp) :: ainv(3,3)
    real(wp) :: det,detinv
    det = a(1,1)*(a(2,2)*a(3,3)-a(2,3)*a(3,2)) &
      & -a(1,2)*(a(2,1)*a(3,3)-a(2,3)*a(3,1)) &
      & +a(1,3)*(a(2,1)*a(3,2)-a(2,2)*a(3,1))
    detinv = 1.0_wp/det
    ainv(1,1) = (a(2,2)*a(3,3)-a(2,3)*a(3,2))*detinv
    ainv(2,1) = (a(2,3)*a(3,1)-a(2,1)*a(3,3))*detinv
    ainv(3,1) = (a(2,1)*a(3,2)-a(2,2)*a(3,1))*detinv
    ainv(1,2) = (a(1,3)*a(3,2)-a(1,2)*a(3,3))*detinv
    ainv(2,2) = (a(1,1)*a(3,3)-a(1,3)*a(3,1))*detinv
    ainv(3,2) = (a(1,2)*a(3,1)-a(1,1)*a(3,2))*detinv
    ainv(1,3) = (a(1,2)*a(2,3)-a(1,3)*a(2,2))*detinv
    ainv(2,3) = (a(1,3)*a(2,1)-a(1,1)*a(2,3))*detinv
    ainv(3,3) = (a(1,1)*a(2,2)-a(1,2)*a(2,1))*detinv
  end function inv3x3_fp

!========================================================================================!
!========================================================================================!
end module pbc_fingerprint_module
