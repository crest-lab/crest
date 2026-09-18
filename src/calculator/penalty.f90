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

!> a small module for getting a penalty contribution, e.g. from the RMSD potential used in metadynamics

module penalty_module
  use crest_parameters
  use strucrd
  use irmsd_module
  implicit none
  private

  public :: penalty_params
  type :: penalty_params
    type(coord),pointer  :: biaslist(:)

    real(wp) :: alpha = 1.0_wp
    real(wp) :: kpush = 0.002_wp
    real(wp),allocatable :: ramp(:)
    real(wp),allocatable :: gradtmp(:,:)
    type(rmsd_core_cache) :: ccache

    character(len=:),allocatable :: biasfile
    type(coord),allocatable :: biastmp(:)
  end type penalty_params

  public :: rmsd_penalty_engrad

!========================================================================================!
!========================================================================================!
contains  !> MODULE PROCEDURES START HERE
!========================================================================================!
!========================================================================================!

  subroutine rmsd_penalty_engrad(mol,ppars,energy,grad,iostatus)
    type(coord),intent(in) :: mol
    type(penalty_params),intent(inout) :: ppars
    real(wp),intent(out) :: energy
    real(wp),intent(out) :: grad(:,:)
    integer,intent(out) :: iostatus
    integer :: nall,io,ii
    real(wp) :: etmp,rmsdval,dEdr,knat
    real(wp),parameter :: thr = sqrt(epsilon(thr))

    iostatus = 0
    energy = 0.0_wp
    grad(:,:) = 0.0_wp
    rmsdval = 0.0_wp
    nall = size(ppars%biaslist,1)
    knat = ppars%kpush*mol%nat

    do ii = 1,nall

      rmsdval = rmsd(mol,ppars%biaslist(ii),gradient=ppars%gradtmp,ccache=ppars%ccache)

      !> energy contribution
      call penalty_potential_gauss(knat,ppars%alpha,rmsdval,etmp,dEdr)
      energy = energy+etmp
      !> fallback: exactly matching structures will produce NaN gradients!
      if (rmsdval < thr) cycle
      !> gradient contribution
      grad(:,:) = grad(:,:)+dEdr*ppars%gradtmp(:,:)
    end do

  end subroutine rmsd_penalty_engrad

!========================================================================================!

  subroutine penalty_potential_gauss(k,a,r,etmp,dEdr)
    real(wp),intent(in) :: k,a,r
    real(wp),intent(out) :: etmp,dEdr
    etmp = k*exp(-a*r**2)
    dEdr = -2.0_wp*a*etmp*r
  end subroutine penalty_potential_gauss

!========================================================================================!
!========================================================================================!
end module penalty_module
