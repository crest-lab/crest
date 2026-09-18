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

module coordinate_transform_module
!***********************************************************************
!* Module: coordinate_transform_module
!*
!* This module provides transformation routines for converting between
!* different coordinate representations. In particular, it includes routines
!* to transform 3D Cartesian coordinates (stored in a 'coord' type from the
!* strucrd module) into a 1D vector representation and vice versa. It also
!* provides routines to transform gradients between a 3D representation (as a
!* 2D array with dimensions 3 x nat) and a 1D vector. These routines are useful
!* in optimization contexts where a flattened variable representation is required.
!***********************************************************************
  use crest_parameters
  use strucrd
  implicit none
  private

  public :: compute_nvar,transform_mol,transform_grd

!========================================================================================!
!========================================================================================!
contains !> MODULE PROCEDURES START HERE
!========================================================================================!
!========================================================================================!

  function compute_nvar(mol) result(nvar)
    !***********************************************************************
    !* Function compute_nvar
    !* Computes the number of variables for a system with nat atoms.
    !* nvar is defined as 3 * nat.
    !***********************************************************************
    implicit none
    type(coord),intent(in) :: mol
    integer :: nvar

    nvar = 3*mol%nat
  end function compute_nvar

!========================================================================================!

  subroutine transform_mol(transformation_type,mol,nvar,vec)
    !***********************************************************************
    !* Subroutine transform_mol
    !* Wrapper routine for coordinate transformations on a molecule.
    !* Supported transformation types:
    !*   "cart2v" - Transforms mol%xyz (3D Cartesian coordinates) into a 1D vector.
    !*   "v2cart" - Transforms a 1D vector into mol%xyz.
    !*
    !* @param transformation_type  Character string specifying the transformation.
    !* @param mol                  Type(coord) variable containing Cartesian coordinates.
    !* @param vec                  1D real(wp) vector (input for "v2cart", output for "cart2v").
    !* @param nvar                 Integer, number of variables (nvar = 3*mol%nat).
    !***********************************************************************
    implicit none
    character(len=*),intent(in) :: transformation_type
    type(coord),intent(inout) :: mol
    real(wp),intent(inout) :: vec(nvar)
    integer,intent(in) :: nvar

    select case (trim(transformation_type))
    case ("cart2v")
      call cartesian_to_vector(mol,vec,nvar)
    case ("v2cart")
      call vector_to_cartesian(vec,nvar,mol)
    case default
      write(*,*)"Error: Transformation type not recognized in transform_mol."
      stop
    end select
  end subroutine transform_mol

!========================================================================================!

  subroutine transform_grd(transformation_type,mol,grd,nvar,vec)
    !***********************************************************************
    !* Subroutine transform_grd
    !* Wrapper routine for gradient transformations.
    !* Supported transformation types:
    !*   "grd2v" - Transforms a 3D gradient array grd(3, nat) into a 1D vector.
    !*   "v2grd" - Transforms a 1D gradient vector into a 3D array grd(3, nat).
    !*
    !* @param transformation_type  Character string specifying the transformation.
    !* @param grd                  3D gradient array (3 x nat); input for "grd2v"
    !*                             and output for "v2grd".
    !* @param vec                  1D real(wp) vector (output for "grd2v", input for "v2grd").
    !* @param nvar                 Integer, number of variables (nvar = 3*nat).
    !***********************************************************************
    implicit none
    character(len=*),intent(in) :: transformation_type
    type(coord),intent(inout) :: mol
    real(wp),intent(inout) :: grd(3,mol%nat)
    real(wp),intent(inout) :: vec(nvar)
    integer,intent(inout) :: nvar
    integer :: nat

    nat = mol%nat
    select case (trim(transformation_type))
    case ("cart2v")
      call gradient_to_vector(grd,nat,vec,nvar)
    case ("v2cart")
      call vector_to_gradient(vec,nvar,grd,nat)
    case default
      write(*,*)"Error: Transformation type not recognized in transform_grd."
      stop
    end select
  end subroutine transform_grd

!========================================================================================!

  subroutine cartesian_to_vector(mol,x,nvar)
    !***********************************************************************
    ! Subroutine cartesian_to_vector
    ! Transforms 3D Cartesian coordinates from mol%xyz into a 1D vector x.
    ! The number of variables is computed as nvar = 3 * mol%nat.
    !***********************************************************************
    implicit none
    type(coord),intent(in) :: mol
    real(wp),intent(out) :: x(nvar)
    integer,intent(in) :: nvar
    integer :: i,j,idx

    x = reshape(mol%xyz, [nvar])
    !idx = 0
    !do j = 1,mol%nat
    !  do i = 1,3
    !    idx = idx+1
    !    x(idx) = mol%xyz(i,j)
    !  end do
    !end do
  end subroutine cartesian_to_vector

  subroutine vector_to_cartesian(x,nvar,mol)
    !***********************************************************************
    ! Subroutine vector_to_cartesian
    ! Transforms a 1D vector x into 3D Cartesian coordinates stored in mol%xyz.
    ! It computes the number of atoms as nat = nvar / 3 and allocates mol%xyz.
    !***********************************************************************
    implicit none
    integer,intent(in) :: nvar
    real(wp),intent(in) :: x(nvar)
    type(coord),intent(inout) :: mol
    integer :: nat,i,j,idx

    if (mod(nvar,3) /= 0) then
      write(*,*)"Error: nvar must be a multiple of 3."
      stop
    end if

    mol%xyz = reshape(x,[3,mol%nat])
!    nat = mol%nat
!    idx = 0
!    do j = 1,nat
!      do i = 1,3
!        idx = idx+1
!        mol%xyz(i,j) = x(idx)
!      end do
!    end do
  end subroutine vector_to_cartesian

!========================================================================================!

  subroutine gradient_to_vector(grd,nat,g,nvar)
    !***********************************************************************
    !* Subroutine gradient_to_vector
    !* Transforms a 3D gradient array grd(3, nat) into a 1D vector g.
    !* The number of variables is computed as nvar = 3 * nat.
    !***********************************************************************
    implicit none
    integer,intent(in) :: nat,nvar
    real(wp),intent(in) :: grd(3,nat)
    real(wp),intent(out) :: g(nvar)
    integer :: i,j,idx

    g = reshape(grd, [nvar])
    !idx = 0
    !do j = 1,nat
    !  do i = 1,3
    !    idx = idx+1
    !    g(idx) = grd(i,j)
    !  end do
    !end do
  end subroutine gradient_to_vector

  subroutine vector_to_gradient(g,nvar,grd,nat)
    !***********************************************************************
    !* Subroutine vector_to_gradient
    !* Transforms a 1D gradient vector g into a 3D gradient array grd(3, nat).
    !* It computes the number of atoms as nat = nvar / 3 and allocates grd.
    !***********************************************************************
    implicit none
    integer,intent(in) :: nvar,nat
    real(wp),intent(in) :: g(nvar)
    real(wp),intent(out) :: grd(3,nat)
    integer :: i,j,idx

    if (mod(nvar,3) /= 0) then
      write(*,*)"Error: nvar must be a multiple of 3."
      stop
    end if

    grd = reshape(g, [3,nat])  
    !idx = 0
    !do j = 1,nat
    !  do i = 1,3
    !    idx = idx+1
    !    grd(i,j) = g(idx)
    !  end do
    !end do
  end subroutine vector_to_gradient

!========================================================================================!
!========================================================================================!
end module coordinate_transform_module

