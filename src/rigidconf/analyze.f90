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

module rigidconf_analyze
  use crest_parameters
  use crest_data
  use strucrd
  use geo
  use INTERNALS_mod
  implicit none
  public

!========================================================================================!
!========================================================================================!
contains !> MODULE PROCEDURES START HERE
!========================================================================================!
!========================================================================================!

subroutine rigidconf_analyze_fallback(env,mol,zmat,na,nb,nc,wbo, &
&                                     ndieder,dvalues,dstep,ztod )
!************************************************************
!* Fallback routine for the dihedral setup.
!* Select all single bonds corresponding to a dihedral angle
!* in the zmatrix. 
!* This automatically excludes terminal atoms (e.g. hydrogen).
!* All the selected dihedrals get a dvalue of 3 and corresponding
!* dstep of 120°
!*
!* Input:
!*    env,mol,zmat,na,nb,nc,wbo,ndieder
!* 
!* Output: 
!*    dvalues, dstep, ztod
!*   
!************************************************************
  implicit none
  !> INPUT
  type(systemdata),intent(inout) :: env
  type(coord),intent(in) :: mol  !> by convention mol is in Bohrs
  real(wp),intent(in) :: zmat(3,mol%nat)
  integer,intent(in)  :: na(mol%nat),nb(mol%nat),nc(mol%nat)
  real(wp),intent(in) :: wbo(mol%nat,mol%nat)
  integer,intent(in)  :: ndieder
  !> OUTPUT 
  integer,intent(out)  :: dvalues( ndieder )
  real(wp),intent(out) :: dstep( ndieder )
  integer,intent(out)  :: ztod( mol%nat ) 
  !> LOCAL
  integer :: V,i,j,k,l,m
  integer,allocatable :: Amap(:,:)
!========================================================================================!
!>--- Init
  dvalues(:) = 3
  dstep(:) = 120.0_wp*degtorad
  ztod(:) = 0

!>--- Mapping
  V = mol%nat
  allocate(Amap(V,V), source=0) !> mapping matrix
  k = 0
  do i=1,V
    m = 0
    if( nc(i) > 0 ) then !> to skip the three incomplete zmate entries
      j=na(i)
      l=nb(i)
      if(Amap(j,l) > 0)then
        m = Amap(j,l)
      else
        if( nint(wbo(j,l)) == 1)then
         k = k + 1
         Amap(j,l) = k
         Amap(l,j) = k 
         m = k
        endif
      endif  
    endif
    ztod(i) = m
  enddo
 
  if(allocated(Amap)) deallocate(Amap)
  return
end subroutine rigidconf_analyze_fallback

!========================================================================================!
subroutine rigidconf_count_fallback(nat,na,nb,nc,wbo,ndieder,ztod)
!************************************************************
!* Count number of unique single-bond dihedral angles that
!* correspond to an entry in the zmat.
!************************************************************
  implicit none
  !> INPUT
  integer,intent(in)  :: nat
  integer,intent(in)  :: na(nat),nb(nat),nc(nat)
  real(wp),intent(in) :: wbo(nat,nat)
  !> OUTPUT
  integer,intent(out)  :: ndieder
  integer,intent(out),optional :: ztod(nat)
  !> LOCAL
  integer :: V,i,j,k,l,m
  integer,allocatable :: Amap(:,:)
!>--- init
  ndieder = 0
!>--- Mapping
  V = nat
  allocate(Amap(V,V), source=0) !> mapping matrix
  k = 0
  do i=1,V
    m = 0
    if( nc(i) > 0 ) then !> to skip the three incomplete zmate entries
      j=na(i)
      l=nb(i)
      if(Amap(j,l) > 0)then
        m = Amap(j,l)
      else
        if( nint(wbo(j,l)) == 1)then
         k = k + 1
         Amap(j,l) = k
         Amap(l,j) = k 
         m = k
        endif
      endif  
    endif
    if(present(ztod)) ztod(i) = m
  enddo
  ndieder = k
  if(allocated(Amap)) deallocate(Amap)
  return
end subroutine rigidconf_count_fallback


!========================================================================================!
subroutine prune_zmat_dihedrals(nat, xyz, zmat, na,nb,nc, ztod )
!********************************************************
!* Remove zmat entries that correspond
!* to the same bond and replace them with internal
!* dihedral angles. There you go, Christoph...
!********************************************************
   implicit none
   integer,intent(in)  :: nat
   real(wp),intent(in) :: xyz(3,nat)
   real(wp),intent(inout) :: zmat(3,nat)
   integer,intent(inout)  :: na(nat),nb(nat),nc(nat)
   integer,intent(inout)  :: ztod(nat)
   integer :: i,j,k,l
   integer :: maxgroup,nmembers,refi
  
   write(*,*) ztod
   maxgroup = maxval(ztod,1)
   do i=1,maxgroup 
     nmembers = count(ztod(:).eq.i)
     if(nmembers < 2) cycle
     do j=1,nat
       if(ztod(j)==i)then
         refi = j 
         exit
       endif 
     enddo
     do j=1,nat
        if(j==refi) cycle
        if(ztod(j) == i)then
           !nc(j) = nb(j)
           nc(j) = refi
           !call BANGLE2( xyz, j, na(j), nb(j), zmat(2,j) )
           call DIHED2( xyz, j, na(j), nb(j), nc(j), zmat(3,j) ) 
           ztod(j) = 0
        endif  
     enddo
   enddo
   write(*,*) ztod
end subroutine prune_zmat_dihedrals
!========================================================================================!
!========================================================================================!
end module rigidconf_analyze
