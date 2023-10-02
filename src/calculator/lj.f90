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
    implicit none
    integer,intent(in) :: n
    real(wp),intent(in) :: xyz(3,n)
    real(wp),intent(in) :: sigma
    real(wp),intent(in) :: epsi

    real(wp),intent(inout) :: energy
    real(wp),intent(inout) :: grad(3,n)

    integer :: i,j,k,l,ich,och,io
    logical :: ex

    real(wp) :: rij
    real(wp) :: U,dU,dx,dy,dz

    energy = 0.0_wp
    grad = 0.0_wp

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

    return
  end subroutine lj_engrad

!========================================================================================!
!========================================================================================!
end module lj
