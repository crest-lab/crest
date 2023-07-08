!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2023 Philipp Pracht
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

module bessel_functions
  use iso_fortran_env,only:wp => real64
  implicit none

  private

  real(wp),parameter :: pi = acos(0.0_wp)*2.0_wp

  public :: bessel_In,bessel_I0,bessel_I1

!========================================================================================!
!========================================================================================!
contains  !> MODULE PROCEDURES START HERE
!========================================================================================!
!========================================================================================!

!************************************************************
!* Modified Bessel function of first kind routines
!* adapted from 
!*
!*      Fortran Numerical Recepies, Chapter 6.6
!*
!************************************************************
  function bessel_I0(x) result(bessel)
!************************************************************
!* Calculates the modified Bessel function of the first
!* kind Iₙ(x), where n=0 and x is a real number
!************************************************************
    implicit none
    real(wp),intent(in) :: x
    real(wp) :: bessel
    real(wp) :: absx,tmp,y
    real(wp) :: p1,p2,p3,p4,p5,p6,p7
    real(wp) :: q1,q2,q3,q4,q5,q6,q7,q8,q9
    SAVE p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7,q8,q9
    DATA p1,p2,p3,p4,p5,p6,p7/1.0d0,3.5156229d0,3.0899424d0,1.2067492d0, &
    & 0.2659732d0,0.360768d-1,0.45813d-2/
    DATA q1,q2,q3,q4,q5,q6,q7,q8,q9/0.39894228d0,0.1328592d-1, &
    & 0.225319d-2,-0.157565d-2,0.916281d-2,-0.2057706d-1,      &
    & 0.2635537d-1,-0.1647633d-1,0.392377d-2/
    absx = abs(x)
    if (absx < 3.75_wp) then
      y = x/3.75_wp
      y = y*y
      bessel = p1+y*(p2+y*(p3+y*(p4+y*(p5+y*(p6+y*p7)))))
    else
      y = 3.75_wp/absx
      tmp = exp(absx)/sqrt(absx)
      bessel = tmp*(q1+y*(q2+y*(q3+y*(q4 &
      &        +y*(q5+y*(q6+y*(q7+y*(q8+y*q9))))))))
    end if
  end function bessel_I0
!========================================================================================!
  function bessel_I1(x) result(bessel)
!************************************************************
!* Calculates the modified Bessel function of the first
!* kind Iₙ(x), where n=0 and x is a real number
!************************************************************
    implicit none
    real(wp),intent(in) :: x
    real(wp) :: bessel
    real(wp) :: absx,y,tmp,tmp2
    real(wp) :: p1,p2,p3,p4,p5,p6,p7
    real(wp) :: q1,q2,q3,q4,q5,q6,q7,q8,q9
    save p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7,q8,q9
    DATA p1,p2,p3,p4,p5,p6,p7/0.5d0,0.87890594d0,0.51498869d0, &
    & 0.15084934d0,0.2658733d-1,0.301532d-2,0.32411d-3/
    DATA q1,q2,q3,q4,q5,q6,q7,q8,q9/0.39894228d0,-0.3988024d-1, &
    & -0.362018d-2,0.163801d-2,-0.1031555d-1,0.2282967d-1,      &
    & -0.2895312d-1,0.1787654d-1,-0.420059d-2/
    absx = abs(x)
    if (absx < 3.75_wp) then
      y = x/3.75_wp
      y = y*y
      bessel = x*(p1+y*(p2+y*(p3+y*(p4+y*(p5+y*(p6+y*p7))))))
    else
      y = 3.75_wp/absx
      tmp = (exp(absx)/sqrt(absx))
      bessel = tmp*(q1+y*(q2+y*(q3+y*(q4+ &
      &       y*(q5+y*(q6+y*(q7+y*(q8+y*q9))))))))
      if (x < 0.0_wp) bessel = -bessel
    end if
    return
  end function bessel_I1
!========================================================================================!
  function bessel_In(n,x) result(bessel)
!************************************************************
!* Calculates the modified Bessel function of the first
!* kind Iₙ(x), where n is an integer and x is a real number
!************************************************************
    implicit none
    integer,intent(in)  :: n
    real(wp),intent(in) :: x
    real(wp) :: bessel
    real(wp) :: absx,tox,bip,bi,bim
    integer :: m,j
    real(wp),parameter :: zerotol = 1.0d-12
    real(wp),parameter :: bigno = 1.0e10
    real(wp),parameter :: bigni = 1.0e-10
    integer,parameter  :: iacc = 40

    absx = abs(x)
    if (n == 0) then
      bessel = bessel_I0(x)
    elseif (n == 1) then
      bessel = bessel_I1(x)
    elseif (absx < zerotol) then
      bessel = 0.0_wp
    else
      tox = 2.0_wp/absx
      bip = 0.0_wp
      bi = 1.0_wp
      bessel = 0.0_wp
      !> Downward recurrence from even m
      !> make IACC larger to increase accuracy
      m = 2*((n+int(sqrt(float(IACC*n)))))
      do j = m,1,-1
        bim = bip+float(j)*tox*bi
        bip = bi
        bi = bim
        !> Renormalize to prevent overflows
        if (abs(bi) > bigno) then
          bessel = bessel*bigni
          bi = bi*bigni
          bip = bip*bigni
        end if
        if (j .eq. n) bessel = bip
      end do
      !> Normalize with bessel_I0
      bessel = bessel*bessel_I0(x)/bi
      if (x < 0.0_wp.and.mod(n,2) .eq. 1) bessel = -bessel
    end if
  end function bessel_In

!========================================================================================!
!========================================================================================!
end module bessel_functions
