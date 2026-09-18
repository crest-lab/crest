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

module ir_spectrum
  use crest_parameters,only:wp,amutokg,metokg
  use atmasses,only:ams
  implicit none
  private
  public :: ir_intensities

contains

  subroutine ir_intensities(nat,at,nat3,evec,dipd,intens)
!***********************************************************************
!* Compute double-harmonic IR intensities from Cartesian dipole
!* gradients projected onto normal mode eigenvectors.
!* Intensities are returned in km/mol.
!*
!* Input:
!*   nat   - number of atoms
!*   at    - atomic numbers (nat)
!*   nat3  - 3*nat
!*   evec  - mass-weighted Hessian eigenvectors, shape (nat3,nat3),
!*           columns are normal modes (output of dsyevd via frequencies())
!*   dipd  - Cartesian dipole gradient ∂μ/∂x_i in a.u., shape (3,nat3)
!* Output:
!*   intens - IR intensities in km/mol, shape (nat3)
!***********************************************************************
    implicit none
    integer,intent(in)   :: nat,nat3
    integer,intent(in)   :: at(nat)
    real(wp),intent(in)  :: evec(nat3,nat3)
    real(wp),intent(in)  :: dipd(3,nat3)
    real(wp),intent(out) :: intens(nat3)
    !> conversion: |∂μ/∂Q|² in a.u. → km/mol
    real(wp),parameter   :: au_to_kmmol = 1.7770969e+6_wp
    real(wp),allocatable :: invmass(:)
    real(wp) :: trdip(3),sum2,amutoau
    integer  :: i,j,k,ii

! ── inverse square root mass vector, one entry per Cartesian coordinate ──
    amutoau = amutokg/metokg
    allocate(invmass(nat3),source=0.0_wp)
    do i = 1,nat
      do j = 1,3
        ii = (i-1)*3+j
        invmass(ii) = 1.0_wp/sqrt(ams(at(i))*amutoau)
      end do
    end do

! ── project dipole gradient onto each normal mode ────────────────────────
    intens = 0.0_wp
    do i = 1,nat3
      do k = 1,3
        sum2 = 0.0_wp
        do j = 1,nat3
          sum2 = sum2+dipd(k,j)*(evec(j,i)*invmass(j))
        end do
        trdip(k) = sum2
      end do
      intens(i) = au_to_kmmol*(trdip(1)**2+trdip(2)**2+trdip(3)**2)
    end do

    deallocate(invmass)
  end subroutine ir_intensities

end module ir_spectrum
