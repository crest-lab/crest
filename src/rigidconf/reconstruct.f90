!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2023 Philipp Pracht, Christopher Zurek, Christoph Bannwarth
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

!========================================================================================!
!========================================================================================!
subroutine construct_new_zmat(nat,zmat,combi,ndieder,dvalues,dstep,ztod,zmat_new)
!************************************************************
!* Given a reference zmatrix and information about
!* dihedral angles, this routine constructs a new zmatrix
!*
!* Input:
!*    zmat,ndieder,dvalues,ztod,dstep,combi
!*
!* Output:
!*    zmat_new
!*
!************************************************************
  use crest_parameters
  use geo
  implicit none
  !> INPUT
  integer,intent(in)  :: nat
  real(wp),intent(in) :: zmat(3,nat)
  integer,intent(in)  :: ndieder
  integer(int8),intent(in) :: combi(ndieder)
  integer,intent(in)  :: dvalues(ndieder)
  real(wp),intent(in) :: dstep(ndieder)
  integer,intent(in)  :: ztod(nat)
  !> OUTPUT
  real(wp),intent(out) :: zmat_new(3,nat)
  !> LOCAL
  integer :: V,i,j,k,l,m
  real(wp) :: xref,xnew
!========================================================================================!

  zmat_new(:,:) = zmat(:,:)
  do i = 1,ndieder
    do j = 1,nat
      if (ztod(j) == i) then
        xref = zmat(3,j)+((float(combi(i))-1.0_wp)*dstep(i))
        xnew = angleshift(xref)
        zmat_new(3,j) = xnew
      end if
    end do
  end do

  return
end subroutine construct_new_zmat
!========================================================================================!
!========================================================================================!
subroutine reconstruct_zmat_to_mol(nat,at,zmat,na,nb,nc,mol)
  use crest_parameters
  use strucrd
  use INTERNALS_mod,only:GMETRY2
  implicit none
  !> INPUT
  integer,intent(in)  :: nat
  integer,intent(in)  :: at(nat)
  real(wp),intent(in) :: zmat(3,nat)
  integer,intent(in)  :: na(nat),nb(nat),nc(nat)
  !> OUTPUT
  type(coord),intent(out) :: mol

  call mol%deallocate()
  allocate (mol%at(nat))
  allocate (mol%xyz(3,nat))
  mol%nat = nat
  mol%at = at
  call GMETRY2(nat,zmat,mol%xyz,na,nb,nc)
end subroutine reconstruct_zmat_to_mol

