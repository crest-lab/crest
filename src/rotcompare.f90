!===============================================================================!
! This file is part of crest.
!
! Copyright (C) 2018-2020 Philipp Pracht
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
!===============================================================================!

!===============================================================================!
! Routines for the comparison of two structures based
! on their rotational constants.
!===============================================================================!
!--- formerly in "eqrot.f"

logical function equalrot(i,j,nall,thr,rot)
  implicit none
  integer i,j,nall
  real*8 rot(3,nall),r,r1,thr
  equalrot = .false.
  r1 = rot(1,i)**2+rot(2,i)**2+rot(3,i)**2
  r = (rot(1,i)-rot(1,j))**2 &
 & +(rot(2,i)-rot(2,j))**2 &
 & +(rot(3,i)-rot(3,j))**2
  !write(*,*) sqrt(r)/sqrt(r1)
  if (sqrt(r)/sqrt(r1) .lt. thr) equalrot = .true.
end function equalrot

logical function equalrot2(i,j,nall,thr,rot)
  implicit none
  integer i,j,nall
  real*8 rot(3,0:nall),r,r1,r2,thr
  equalrot2 = .false.
  r1 = rot(1,i)**2+rot(2,i)**2+rot(3,i)**2
  r2 = rot(1,j)**2+rot(2,j)**2+rot(3,j)**2
  r = (rot(1,i)-rot(1,j))**2 &
 & +(rot(2,i)-rot(2,j))**2 &
 & +(rot(3,i)-rot(3,j))**2
  if (2.*sqrt(r)/(sqrt(r1)+sqrt(r2)) .lt. thr) equalrot2 = .true.
end function equalrot2

real*8 function rotdiff(i,j,nall,rot)
  implicit none
  integer i,j,nall
  real*8 rot(3,nall),r,r1,r2
  r1 = rot(1,i)**2+rot(2,i)**2+rot(3,i)**2
  r2 = rot(1,j)**2+rot(2,j)**2+rot(3,j)**2
  r = (rot(1,i)-rot(1,j))**2 &
 & +(rot(2,i)-rot(2,j))**2 &
 & +(rot(3,i)-rot(3,j))**2
  rotdiff = 2.*sqrt(r)/(sqrt(r1)+sqrt(r2))
end function rotdiff

!- use an absolute value as threshold
logical function equalrotall(i,j,nall,thr,rot)
  implicit none
  integer i,j,nall
  real*8 rot(3,nall),thr
  logical :: r1,r2,r3
  equalrotall = .false.
  r1 = abs(rot(1,i)-rot(1,j)) .le. thr
  r2 = abs(rot(2,i)-rot(2,j)) .le. thr
  r3 = abs(rot(3,i)-rot(3,j)) .le. thr
  equalrotall = r1.and.r2.and.r3
  return
end function equalrotall

!- use a relative value as threshold
logical function equalrotallrel(i,j,nall,thr,rot)
  implicit none
  integer i,j,nall
  real*8 rot(3,nall),thr
  logical :: r1,r2,r3
  equalrotallrel = .false.
  r1 = abs((rot(1,i)/rot(1,j))-1.0d0) .le. thr
  r2 = abs((rot(2,i)/rot(2,j))-1.0d0) .le. thr
  r3 = abs((rot(3,i)/rot(3,j))-1.0d0) .le. thr
  equalrotallrel = r1.and.r2.and.r3
  return
end function equalrotallrel

logical function equalrotmean(i,j,nall,thr,rot)
  implicit none
  integer i,j,nall
  real*8 rot(3,nall),r,thr
  equalrotmean = .false.
  r = abs(rot(1,i)-rot(1,j))
  r = r+abs(rot(2,i)-rot(2,j))
  r = r+abs(rot(3,i)-rot(3,j))
  r = r/3.0d0
  equalrotmean = r .le. thr
  return
end function equalrotmean

!===========================================================!
! anisotropy related functions

! calculate rot.const. anisotropy for a single structure
function rotaniso(i,nall,rot)
  use iso_fortran_env,wp => real64
  implicit none
  real(wp) :: rotaniso
  real(wp) :: rot(3,nall)
  integer :: i,nall
  real(wp) :: a,b,c,av
  a = rot(1,i); b = rot(2,i); c = rot(3,i)
  av = (a+b+c)/3.0_wp
  rotaniso = sqrt((a-av)**2+(b-av)**2+(c-av)**2)
  rotaniso = rotaniso/av
  rotaniso = rotaniso/(3.0_wp*sqrt(2.0_wp/3.0_wp))
  return
end function rotaniso

!the threshold used for the rotational constant comparison is
!is modified based on the anisotropy of the rot. constants
!the scaling function is an error function
function bthrerf(bthr,aniso,bthrmax,bthrshift) result(thr)
  use iso_fortran_env,wp => real64
  implicit none
  real(wp) :: bthr,bthrmax,bthrshift
  real(wp) :: aniso
  real(wp) :: thr
  real(wp) :: a,b,c,d
  thr = bthr
  c = ((bthrmax*100.0_wp)-(bthr*100.0_wp))/2.0_wp
  a = -erf(-2.5_wp)*c+(bthr*100.0_wp) ! the y-axis shift
  b = 4.0_wp/0.8_wp ! x-axis range from bthr to bthrmax
  d = bthrshift/0.15_wp
  thr = erf(aniso*b-d)*c+a
  thr = thr/100.0_wp
  return
end function bthrerf

! compare each rotational constant with a modified bthr threshold
! bthr is a relative value threshold
logical function equalrotaniso(i,j,nall,rot,bthr,bthrmax,bthrshift)
  use iso_fortran_env,wp => real64
  implicit none
  integer i,j,nall
  real(wp) :: rot(3,nall)
  real(wp) :: bthr
  real(wp) :: bthrmax,bthrshift
  real(wp) :: anisoi,anisoj,av
  real(wp) :: thr
  real(wp) :: rotaniso !this is a function
  real(wp) :: bthrerf  !this is a function
  logical :: r1,r2,r3
  equalrotaniso = .false.
  anisoi = rotaniso(i,nall,rot)
  anisoj = rotaniso(j,nall,rot)
  av = (anisoi+anisoj)/2.0d0
  !av=min(anisoi,anisoj)
  thr = bthrerf(bthr,av,bthrmax,bthrshift)
  r1 = abs((rot(1,i)/rot(1,j))-1.0d0) .le. thr
  r2 = abs((rot(2,i)/rot(2,j))-1.0d0) .le. thr
  r3 = abs((rot(3,i)/rot(3,j))-1.0d0) .le. thr
  equalrotaniso = r1.and.r2.and.r3
  return
end function equalrotaniso

