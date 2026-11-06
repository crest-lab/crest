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

module rmsdpot
  use strucrd
  use iso_fortran_env,only:wp => real64

  implicit none
  private

  type :: rmsdbias
    integer :: nbias = 0
    real(wp),allocatable :: alpha(:)
    real(wp),allocatable :: kpush(:)
    integer,allocatable  :: mult(:)
    type(coord),pointer :: ptr_structures(:)
  end type rmsdbias

  public :: rmsdbias

!========================================================================================!
!========================================================================================!
contains  !> MODULE PROCEDURES START HERE
!========================================================================================!
!========================================================================================!

  subroutine rmsd_push_engrad(mol,rbias,energy,grad)
!*************************************************************************
!* Compute a repulsive energy and corresponding forces for
!* the similarity match between the currnt mol and a list of references
!*************************************************************************
    implicit none
    type(coord),intent(in) :: mol
    type(rmsdbias) :: rbias

    real(wp),intent(inout) :: energy
    real(wp),intent(inout) :: grad(3,mol%nat)
    integer :: i,j,k,l

    real(wp) :: tmpe,ktot,rmssq 
    real(wp),allocatable :: tmpgrad(:,:)

    energy = 0.0_wp
    grad = 0.0_wp

    do i=1,rbias%nbias
      rmssq = 0.0_wp ** 2
      ktot = real(rbias%mult(i))*rbias%kpush(i)*real(mol%nat)
      tmpe = ktot*exp(-rbias%alpha(i)*rmssq )
       

    enddo

    return
  end subroutine rmsd_push_engrad

!========================================================================================!
!========================================================================================!
end module rmsdpot
