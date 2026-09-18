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

module rotaniso_mod
  use crest_parameters
  implicit none
  public

contains

  function rotaniso(rot) result(aniso)
!***********************************************************
!* calculate rot.const. anisotropy for a single structure
!***********************************************************
    implicit none
    real(wp) :: aniso
    real(wp),intent(in) :: rot(3)
    real(wp) :: a,b,c,av
    a = rot(1); b = rot(2); c = rot(3)
    !av = (a+b+c)/3.0_wp
    !aniso = sqrt((a-av)**2+(b-av)**2+(c-av)**2)
    !aniso = rotaniso/av
    !aniso = rotaniso/(3.0_wp*sqrt(2.0_wp/3.0_wp))

    !> the following is identical to the commented out part
    aniso = sqrt(a**2+b**2+c**2-a*b-a*c-b*c)/(a+b+c)
    return
  end function rotaniso

  function bthrerf(bthr,aniso,bthrmax,bthrshift) result(thr)
!***************************************************************
!* the threshold used for the rotational constant comparison is
!* is modified based on the anisotropy of the rot. constants
!* the scaling function is an error function
!***************************************************************
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

  logical function equalrotaniso(i,j,nall,rot,bthr,bthrmax,bthrshift)
!*********************************************************************
!* compare each rotational constant with a modified bthr threshold
!* bthr is a relative value (fraction) threshold
!*********************************************************************
    implicit none
    integer i,j,nall
    real(wp) :: rot(3,nall)
    real(wp) :: bthr
    real(wp) :: bthrmax,bthrshift
    real(wp) :: anisoi,anisoj,av
    real(wp) :: thr
    logical :: r1,r2,r3
    equalrotaniso = .false.
    anisoi = rotaniso(rot(:,i))
    anisoj = rotaniso(rot(:,j))
    av = (anisoi+anisoj)/2.0d0
    !av=min(anisoi,anisoj)
    thr = bthrerf(bthr,av,bthrmax,bthrshift)
    r1 = abs((rot(1,i)/rot(1,j))-1.0d0) .le. thr
    r2 = abs((rot(2,i)/rot(2,j))-1.0d0) .le. thr
    r3 = abs((rot(3,i)/rot(3,j))-1.0d0) .le. thr
    equalrotaniso = r1.and.r2.and.r3
    return
  end function equalrotaniso

end module rotaniso_mod

!===============================================================================!
! LEGACY Routines for the comparison of two structures based
! on their rotational constants.
!===============================================================================!

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

