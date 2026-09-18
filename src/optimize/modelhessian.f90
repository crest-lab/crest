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
!
! Routines were adapted from the xtb code (github.com/grimme-lab/xtb)
! under the Open-source software LGPL-3.0 Licencse.
!================================================================================!
module modelhessian_module
  use iso_fortran_env,only:wp => real64,stdout => output_unit
  use crest_calculator,only:calcdata,constrhess
  use modelhessian_core
  implicit none

  public :: modhes

!==============================================================================!
contains  !> MODULE PROCEDURES START HERE
!==============================================================================!
!
  subroutine modhes(calc,modh,natoms,xyz,at,Hess,pr)
!**********************************************************
!* subroutine modhes
!* create a model Hessian for a given molecule
!*
!* Input:
!*     natoms - number of atoms
!*       xyz  - Cartesian coordinates
!*        at  - atom types as integers
!*      modh  - model Hessian settings (see above)
!*      calc  - calculation settings (for constraints)
!*        pr  - printout selection
!*
!* Output:
!*      Hess  - the (packed) model Hessian
!**********************************************************
    implicit none
    type(calcdata),intent(in) :: calc
    type(mhparam),intent(in) :: modh
    logical,intent(in) :: pr
    integer :: i
    integer :: nhess
    integer,intent(in) :: natoms
    real(wp),intent(in) :: xyz(3,natoms)
    real(wp),intent(out) :: hess((natoms*3)*((natoms*3)+1)/2)
    integer,intent(in) :: at(natoms)

!>  initialize
    nhess = 3*natoms
    Hess = 0.0_wp

    select case (modh%model)
    case (0)
      if (pr) write (stdout,'(a)') "Using Lindh-Hessian (1995)"
      call ddvopt(xyz,natoms,Hess,at,modh)
!> other model hessians currently not tested
    case (1)
      if (pr) write (stdout,'(a)') "Using Lindh-Hessian"
      call mh_lindh_d2(xyz,natoms,Hess,at,modh)
    case (2)
      if (pr) write (stdout,'(a)') "Using Lindh-Hessian (2007)"
      call mh_lindh(xyz,natoms,Hess,at,modh)
    case (3)
      if (pr) write (stdout,'(a)') "Using Swart-Hessian"
      call mh_swart(xyz,natoms,Hess,at,modh)
    end select

!> add user-set constraint contributions to modelhessian
    call constrhess(natoms,at,xyz,calc,Hess)

    return
  end subroutine modhes

!========================================================================================!
!########################################################################################!
!========================================================================================!
end module modelhessian_module
