!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2020 Philipp Pracht
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
! classical quicksort algorithm, sort HIGH-to-LOW
!=============================================================!
recursive subroutine revquicksort(n,arr)
  implicit none
  integer :: n,arr(n),i,j,k,m
  integer :: pivot
  integer,allocatable :: R(:),L(:)
  integer :: rr,ll,rc,lc,pp

  if (n .le. 1) return

  pivot = arr(1)
  pp = 0
  do i = 1,n
    if (arr(i) .eq. pivot) pp = pp+1
  end do

  ll = 0
  do i = 1,n
    if (arr(i) .ge. pivot) then
      ll = ll+1
    end if
  end do
  ll = ll-pp
  rr = n-ll-pp
  allocate (L(ll),R(rr))

  lc = 0
  rc = 0
  do j = 1,n
    if (arr(j) .gt. pivot) then
      lc = lc+1
      L(lc) = arr(j)
    else if (arr(j) .lt. pivot) then
      rc = rc+1
      R(rc) = arr(j)
    end if
  end do

  call revquicksort(ll,L)
  call revquicksort(rr,R)

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
end subroutine revquicksort

!=============================================================!
! other variant of quicksort algos
!=============================================================!
recursive subroutine qsort(a,first,last,ind)
  implicit none
  real*8 a(*),x,t
  integer ind(*)
  integer first,last
  integer i,j,ii

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
  implicit none
  real*8 a(*),x,t
  integer first,last
  integer i,j

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

recursive subroutine maskqsort(a,first,last,mask)
  implicit none
  real*8 a(*),t
  integer x
  integer mask(*)
  integer first,last
  integer i,j,ii

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
  implicit none
  integer :: adim,nall
  real*8 a(adim,nall),adum(adim)
  integer x
  integer mask(nall)
  integer first,last
  integer i,j,ii

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

recursive subroutine stringqsort(sdim,strs,first,last,mask)
  implicit none
  integer :: sdim
  character(len=*) :: strs(sdim)
  character(len=len(strs(1))) :: str
  integer x
  integer mask(sdim)
  integer first,last
  integer i,j,ii
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
  if (first < i-1) call stringqsort(sdim,strs,first,i-1,mask)
  if (j+1 < last) call stringqsort(sdim,strs,j+1,last,mask)
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
