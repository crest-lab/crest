!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2020-2024 Philipp Pracht
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

module quicksort_interface
!********************************************************
!* module to load an interface to the quicksort routines
!* mandatory to handle optional input arguments
!********************************************************
  implicit none
  interface
    recursive subroutine quicksort(n,arr)
      implicit none
      integer :: n,arr(n)
    end subroutine quicksort

    recursive subroutine qsort(a,first,last,ind)
      use iso_fortran_env,only:wp => real64
      implicit none
      real(wp) :: a(*)
      integer :: ind(*)
      integer :: first,last
    end subroutine qsort

    recursive subroutine qqsort(a,first,last)
      use iso_fortran_env,only:wp => real64
      implicit none
      real(wp) :: a(*)
      integer :: first,last
    end subroutine qqsort

    recursive subroutine qqsorti(a,first,last)
      implicit none
      integer :: a(*)
      integer :: first,last
    end subroutine qqsorti

    recursive subroutine maskqsort(a,first,last,mask)
      use iso_fortran_env,only:wp => real64
      implicit none
      real(wp) :: a(*)
      integer :: first,last
      integer :: mask(*)
    end subroutine maskqsort

    recursive subroutine matqsort(adim,nall,a,adum,first,last,mask)
      use iso_fortran_env,only:wp => real64
      implicit none
      integer  :: adim,nall
      real(wp) :: a(adim,nall),adum(adim)
      integer :: first,last
      integer :: mask(nall)
    end subroutine matqsort

    recursive subroutine stringqsort(sdim,slen,strs,first,last,mask)
      implicit none
      integer,intent(in) :: sdim,slen
      character(len=slen) :: strs(sdim)
      integer :: first,last
      integer :: mask(sdim)
    end subroutine stringqsort

    subroutine maskinvert(nall,mask)
      implicit none
      integer :: nall
      integer :: mask(nall)
    end subroutine maskinvert

    recursive subroutine ensemble_qsort(nall,structures,first,last,mask)
      use crest_parameters
      use strucrd,only:coord
      implicit none
      integer,intent(in) :: nall
      type(coord),intent(inout) :: structures(nall)
      integer,intent(in) :: first,last
      integer,intent(inout),optional :: mask(nall)
    end subroutine ensemble_qsort

  end interface
end module quicksort_interface

!=============================================================!
! classical quicksort algorithm, sort LOW-to-HIGH
!=============================================================!
recursive subroutine quicksort(n,arr)
  implicit none
  integer :: n,arr(n),i,j,k,m
  integer :: pivot
  integer,allocatable :: R(:),L(:)
  integer :: rr,ll,rc,lc,pp

  if (n .le. 1) return

  pivot = arr(1)
  if (arr(2) .lt. arr(1)) pivot = arr(2)
  pp = 0
  do i = 1,n
    if (arr(i) .eq. pivot) pp = pp+1
  end do

  ll = 0
  do i = 1,n
    if (arr(i) .le. pivot) then
      ll = ll+1
    end if
  end do
  ll = ll-pp
  rr = n-ll-pp
  allocate (L(ll),R(rr))

  lc = 0
  rc = 0
  do j = 1,n
    if (arr(j) .lt. pivot) then
      lc = lc+1
      L(lc) = arr(j)
    else if (arr(j) .gt. pivot) then
      rc = rc+1
      R(rc) = arr(j)
    end if
  end do

  call quicksort(ll,L)
  call quicksort(rr,R)

  do i = 1,ll
    arr(i) = L(i)
  end do
  do k = 1,pp
    m = k+ll
    arr(m) = pivot
  end do
  do j = 1,rr
    m = j+ll+pp
    arr(m) = R(j)
  end do

  deallocate (R,L)
end subroutine quicksort

!=============================================================!
! other variant of quicksort algos
!=============================================================!
recursive subroutine qsort(a,first,last,ind)
  use iso_fortran_env,only:wp => real64
  implicit none
  real(wp) :: a(*)
  real(wp) :: x,t
  integer :: ind(*)
  integer :: first,last,i,j,ii

  x = a((first+last)/2)
  i = first
  j = last
  do
    do while (a(i) < x)
      i = i+1
    end do
    do while (x < a(j))
      j = j-1
    end do
    if (i >= j) exit
    t = a(i); a(i) = a(j); a(j) = t
    ii = ind(i); ind(i) = ind(j); ind(j) = ii
    i = i+1
    j = j-1
  end do
  if (first < i-1) call qsort(a,first,i-1,ind)
  if (j+1 < last) call qsort(a,j+1,last,ind)
end subroutine qsort

recursive subroutine qqsort(a,first,last)
  use iso_fortran_env,only:wp => real64
  implicit none
  real(wp) :: a(*)
  real(wp) :: x,t
  integer :: first,last,i,j

  x = a((first+last)/2)
  i = first
  j = last
  do
    do while (a(i) < x)
      i = i+1
    end do
    do while (x < a(j))
      j = j-1
    end do
    if (i >= j) exit
    t = a(i); a(i) = a(j); a(j) = t
    i = i+1
    j = j-1
  end do
  if (first < i-1) call qqsort(a,first,i-1)
  if (j+1 < last) call qqsort(a,j+1,last)
end subroutine qqsort

recursive subroutine qqsorti(a,first,last)
  implicit none
  integer :: a(*)
  integer :: x,t
  integer :: first,last,i,j

  x = a((first+last)/2)
  i = first
  j = last
  do
    do while (a(i) < x)
      i = i+1
    end do
    do while (x < a(j))
      j = j-1
    end do
    if (i >= j) exit
    t = a(i); a(i) = a(j); a(j) = t
    i = i+1
    j = j-1
  end do
  if (first < i-1) call qqsorti(a,first,i-1)
  if (j+1 < last) call qqsorti(a,j+1,last)
end subroutine qqsorti

recursive subroutine maskqsort(a,first,last,mask)
  use iso_fortran_env,only:wp => real64
  implicit none
  real(wp) :: a(*)
  real(wp) :: t
  integer :: x,first,last,i,j,ii
  integer :: mask(*)

  x = mask((first+last)/2)
  i = first
  j = last
  do
    do while (mask(i) < x)
      i = i+1
    end do
    do while (x < mask(j))
      j = j-1
    end do
    if (i >= j) exit
    t = a(i); a(i) = a(j); a(j) = t
    ii = mask(i); mask(i) = mask(j); mask(j) = ii
    i = i+1
    j = j-1
  end do
  if (first < i-1) call maskqsort(a,first,i-1,mask)
  if (j+1 < last) call maskqsort(a,j+1,last,mask)
end subroutine maskqsort

recursive subroutine matqsort(adim,nall,a,adum,first,last,mask)
  use iso_fortran_env,only:wp => real64
  implicit none
  integer  :: adim,nall
  real(wp) :: a(adim,nall),adum(adim)
  integer :: x,first,last,i,j,ii
  integer :: mask(nall)

  x = mask((first+last)/2)
  i = first
  j = last
  do
    do while (mask(i) < x)
      i = i+1
    end do
    do while (x < mask(j))
      j = j-1
    end do
    if (i >= j) exit
    adum(:) = a(:,i); a(:,i) = a(:,j); a(:,j) = adum(:)
    ii = mask(i); mask(i) = mask(j); mask(j) = ii
    i = i+1
    j = j-1
  end do
  if (first < i-1) call matqsort(adim,nall,a,adum,first,i-1,mask)
  if (j+1 < last) call matqsort(adim,nall,a,adum,j+1,last,mask)
end subroutine matqsort

recursive subroutine stringqsort(sdim,slen,strs,first,last,mask)
  implicit none
  integer,intent(in) :: sdim,slen
  character(len=slen) :: strs(sdim)
  character(len=slen) :: str
  integer :: x,first,last,i,j,ii
  integer :: mask(sdim)
  x = mask((first+last)/2)
  i = first
  j = last
  do
    do while (mask(i) < x)
      i = i+1
    end do
    do while (x < mask(j))
      j = j-1
    end do
    if (i >= j) exit
    str = strs(i); strs(i) = strs(j); strs(j) = str
    ii = mask(i); mask(i) = mask(j); mask(j) = ii
    i = i+1
    j = j-1
  end do
  if (first < i-1) call stringqsort(sdim,slen,strs,first,i-1,mask)
  if (j+1 < last) call stringqsort(sdim,slen,strs,j+1,last,mask)
end subroutine stringqsort

subroutine maskinvert(nall,mask)
  implicit none
  integer :: nall
  integer :: mask(nall)
  integer,allocatable :: imask(:)
  integer :: i
  allocate (imask(nall))
  do i = 1,nall
    imask(mask(i)) = i
  end do
  mask(:) = imask(:)
  deallocate (imask)
  return
end subroutine maskinvert

!========================================================================================!

recursive subroutine ensemble_qsort(nall,structures,first,last,mask)
  use crest_parameters
  use strucrd,only:coord
  implicit none
  integer,intent(in) :: nall
  type(coord),intent(inout) :: structures(nall)
  integer,intent(in) :: first,last
  integer,intent(inout),optional :: mask(nall)

  !> LOCAL
  type(coord),allocatable :: tmpmol
  integer :: i,j,mm,ii
  real(wp) :: ee

  if (present(mask)) then
!>--- sort according to a given mask (reference order)
    mm = mask((first+last)/2)
    i = first
    j = last
    do
      do while (mask(i) < mm)
        i = i+1
      end do
      do while (mm < mask(j))
        j = j-1
      end do
      if (i >= j) exit
      ii = mask(i); mask(i) = mask(j); mask(j) = ii
      allocate (tmpmol)
      tmpmol = structures(i); structures(i) = structures(j); structures(j) = tmpmol
      deallocate (tmpmol)
      i = i+1
      j = j-1
    end do
    if (first < i-1) call ensemble_qsort(nall,structures,first,i-1,mask)
    if (j+1 < last) call ensemble_qsort(nall,structures,j+1,last,mask)

  else
!>--- standard, sort according to energy of structures
    ee = structures((first+last)/2)%energy
    i = first
    j = last
    do
      do while (structures(i)%energy < ee)
        i = i+1
      end do
      do while (ee < structures(j)%energy)
        j = j-1
      end do
      if (i >= j) exit
      allocate (tmpmol)
      tmpmol = structures(i); structures(i) = structures(j); structures(j) = tmpmol
      deallocate (tmpmol)
      i = i+1
      j = j-1
    end do
    if (first < i-1) call ensemble_qsort(nall,structures,first,i-1)
    if (j+1 < last) call ensemble_qsort(nall,structures,j+1,last)
  end if
end subroutine ensemble_qsort

