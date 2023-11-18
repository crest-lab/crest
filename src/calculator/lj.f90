!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2021 - 2022 Philipp Pracht
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

!> a small module containing the Lennard-Jones potential
!> (mainly for testing other implementations)

module lj

  use iso_fortran_env,only:wp => real64

  implicit none
  private

  public :: lj_engrad

!========================================================================================!
!========================================================================================!
contains  !> MODULE PROCEDURES START HERE
!========================================================================================!
!========================================================================================!

  subroutine lj_engrad(n,xyz,epsi,sigma,energy,grad)
!> Compute the total energy and corresponding forces for
!> a system of nPart interacting LJ particles.
!> The total energy is given by:
!>
!>  U(rⁿ) = ∑_i∑_j U( r_ij )    for i≠j
!>
!> Lennard-Jones potential U(r) =  4*ε*[ (σ/r)^12 - (σ/r)^6]
!> Here, we set σ=1 and ε=1 for reduced LJ units:
!>
!>  U(r) =  4*[ (1/r)^12 - (1/r)^6]
!>
!>  dU(r) = 4*(x/r)*[ 12*(1/r)^13 - 6*(1/r)^7]
!>        = 48 * x * (1/r)^8 *[(1/r)^6 - 0.5]
!>
    implicit none
    integer,intent(in) :: n
    real(wp),intent(in) :: xyz(3,n)
    real(wp),intent(in) :: sigma
    real(wp),intent(in) :: epsi

    real(wp),intent(inout) :: energy
    real(wp),intent(inout) :: grad(3,n)

    integer :: i,j,k,l,ich,och,io
    logical :: ex

    real(wp) :: rij,rij2,dr(3)
    real(wp) :: invR2,invR6,sig2
    real(wp) :: U,dU,dx,dy,dz

    energy = 0.0_wp
    grad = 0.0_wp

!>--- standard implementation
    do i = 1,n
      do j = 1,i-1
        rij = (xyz(1,i)-xyz(1,j))**2   &
        &   +(xyz(2,i)-xyz(2,j))**2   &
        &   +(xyz(3,i)-xyz(3,j))**2
        rij = sqrt(rij)
        !>--- energy
        U = 4.0_wp*epsi*((sigma/rij)**12-(sigma/rij)**6)
        energy = energy+U
        !>--- Cartesian derivative
        dU = 4.0_wp*epsi*(6.0_wp*(sigma**6)* &
        & (rij**-7)-12.0_wp*(sigma**12)*(rij**(-13)))
        dU = dU/rij
        dx = dU*(xyz(1,i)-xyz(1,j))
        dy = dU*(xyz(2,i)-xyz(2,j))
        dz = dU*(xyz(3,i)-xyz(3,j))
        grad(1,i) = grad(1,i)+dx
        grad(1,j) = grad(1,j)-dx
        grad(2,i) = grad(2,i)+dy
        grad(2,j) = grad(2,j)-dy
        grad(3,i) = grad(3,i)+dz
        grad(3,j) = grad(3,j)-dz
      end do
    end do

!>--- efficient implementation (needs debugging)
!    sig2 = sigma**2
!    !$omp parallel &
!    !$omp shared(xyz,energy,grad) &
!    !$omp private(i,j,dr,invR2,invR6,U,dU,dz,dx,dy) 
!    do i = 1,n
!      do j = 1,i-1
!        dr(:) = xyz(:,i)-xyz(:,j)
!        rij2 = dot_product(dr,dr)
!        invR2 = sig2/rij2
!        invR6 = invR2**3
!        !>--- energy
!        U = invR6*(invR6-1.0_wp)
!        !>--- Cartesian derivative
!        dU = (invR2**4)*(invR6-0.5d0)
!        dx = dU*dr(1)
!        dy = dU*dr(2)
!        dz = dU*dr(3)
!        !$omp critical
!        energy = energy + U
!        grad(1,i) = grad(1,i)+dx
!        grad(1,j) = grad(1,j)-dx
!        grad(2,i) = grad(2,i)+dy
!        grad(2,j) = grad(2,j)-dy
!        grad(3,i) = grad(3,i)+dz
!        grad(3,j) = grad(3,j)-dz
!        !$omp end critical
!      end do
!    end do
!    !$omp end parallel
!    energy = energy*(4.0_wp*epsi)
!    grad(:,:) = grad(:,:)*(48.0_wp*epsi)

    return
  end subroutine lj_engrad

!========================================================================================!
!========================================================================================!
end module lj
