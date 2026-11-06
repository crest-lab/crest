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

!> a small module for getting a free energy as engrad call

module approxg_module
  use crest_parameters
  use modelhessian_core
  use thermochem_module
  use strucrd
  use optimize_maths,only:dhtosq
  implicit none
  private

  public :: approxg_params
  type :: approxg_params
     integer :: dim = 0
     logical :: pr = .false.
     real(wp) :: T = 298.15_wp
     real(wp),allocatable :: hess(:,:)
     real(wp),allocatable :: h(:)
     real(wp),allocatable :: freq(:)
     real(wp),allocatable :: xyz(:,:)
     real(wp) :: fscal = 1.0_wp
     real(wp) :: ithr = -50.0_wp
     real(wp) :: sthr = 50.0_wp
  end type approxg_params

  public :: modh_engrad

!========================================================================================!
!========================================================================================!
contains  !> MODULE PROCEDURES START HERE
!========================================================================================!
!========================================================================================!

  subroutine modh_engrad(mol,ag,dg,dggrad,iostatus)
    type(coord),intent(in) :: mol
    type(approxg_params),intent(inout) :: ag
    real(wp),intent(out) :: dg
    real(wp),intent(out) :: dggrad(:,:)
    integer,intent(out) :: iostatus
    integer :: n3,io
    integer,parameter :: nt = 1
    real(wp) :: temps(nt),et(nt),ht(nt),gt(nt),stot(nt)

    type(mhparam) :: mhset

    iostatus = 0
    dg = 0.0_wp
    dggrad(:,:) = 0.0_wp
    temps(1) = ag%T

    !> setup
    n3 = mol%nat*3
    call ddvopt(mol%xyz,mol%nat,ag%h,mol%at,mhset)

    call dhtosq(n3,ag%hess,ag%h)
    ag%h(:) = 0.0_wp

    call prj_mw_hess(mol%nat,mol%at,n3,mol%xyz, &
      &  ag%hess,ag%h)

    call frequencies(mol%nat,mol%at,mol%xyz,n3,ag%hess,ag%freq,io)
    iostatus = io
    if (iostatus .ne. 0) return

    ag%xyz(:,:) = mol%xyz(:,:)
    call calcthermo(mol%nat,mol%at,ag%xyz,ag%freq,ag%pr, &
                    ag%ithr,ag%fscal,ag%sthr,nt,temps,et,ht,gt,stot)

    dg = gt(1)
  end subroutine modh_engrad

!========================================================================================!
!========================================================================================!
end module approxg_module
