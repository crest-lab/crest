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

module probabilities_module
  use iso_fortran_env,only:wp => real64
  use bessel_functions
  implicit none

  private

  real(wp),parameter :: pi = acos(0.0_wp)*2.0_wp
  real(wp),parameter,public :: radtodeg = 180.0_wp / pi
  real(wp),parameter,public :: degtorad = 1.0_wp / radtodeg

  public :: normal_dist

  public :: vonMises
  interface vonMises
    module procedure vonMises_fraction
    module procedure vonMises_avg
  end interface vonMises
  public :: vonMises_plot

!========================================================================================!
!========================================================================================!
contains  !> MODULE PROCEDURES START HERE
!========================================================================================!
!========================================================================================!

  subroutine normal_dist(x,mu,sigma,pdf,d_pdf,d2_pdf)
!*******************************************************************
!* This subroutine computes the Gaussian probability density 
!*
!*                1
!* ρ(x|μ,σ) = ──────────  exp(-1/2 * ((x -  μ)/σ)²)
!*            σ sqrt(2π)
!*
!* at x, where σ is the standard deviation, σ² is the variance, and
!* μ is the mean.
!*******************************************************************
    real(wp),intent(in) :: x        !> Input value
    real(wp),intent(in) :: mu       !> Mean
    real(wp),intent(in) :: sigma    !> Standard deviation
    real(wp),intent(out) :: pdf     !> Probability density function value
    real(wp),intent(out),optional :: d_pdf   !> First derivative of the PDF
    real(wp),intent(out),optional :: d2_pdf  !> Second derivative of the PDF

    real(wp) :: expo,coefficient

    expo = -((x-mu)**2)/(2.0_wp*sigma**2)
    coefficient = 1.0_wp/(sigma*sqrt(2.0_wp*pi))

    pdf = coefficient*exp(expo)
    if (present(d_pdf)) then
      d_pdf = -((x-mu)/sigma**2)*pdf
    end if
    if (present(d2_pdf)) then
      d2_pdf = ((x-mu)**2/sigma**4-1.0_wp/sigma**2)*pdf
    end if
  end subroutine normal_dist

!=======================================================================================!

  subroutine vonMises_fraction(theta,kappa,mu,N,p,dpdt,dpdt2)
!********************************************************************
!* The subroutine computes a smooth probability density
!* ρ(θ|μ,κ) as the average of N von-Mises kernels, which are
!* continuous probability distributions around a specific angle.
!* It takes into account both the central location parameter,
!* i.e. the N observations μ of the angle θ, and the concentration
!* parameter κ.
!* The latter is obtained from the "von Mises-scale plug-in rule".
!* ρ(θ|μ,κ) is given by
!*
!*                1
!* ρ(θ|μ,κ) = ─────────  exp(κ' cos(θ -  μᵢ))
!*            2πNI₀(κ')
!*
!*
!* where I₀ is the modiﬁed Bessel function of order 0 and κ' is
!* computed from κ via the "scale plug-in rule" (see below).
!*
!* This routine is the version computes the N-th fraction of the
!* average and adds up the contribution to p.
!*******************************************************************
    implicit none
    !> INPUT
    real(wp),intent(in) :: theta
    real(wp),intent(in) :: kappa
    real(wp),intent(in) :: mu
    integer,intent(in)  :: N
    !> OUTPUT
    real(wp),intent(inout) :: p
    real(wp),intent(inout),optional :: dpdt  !> first derivative wrt θ
    real(wp),intent(inout),optional :: dpdt2 !> second derivcative wrt θ
    !> LOCAL
    integer :: i,j
    real(wp) :: Nf,skappa
    real(wp) :: tmp,tmp2,tmp3
    real(wp) :: tmpdpdt,tmpdpdt2

    Nf = float(N)

    tmpdpdt = 0.0_wp
    tmpdpdt2 = 0.0_wp

    skappa = vonMises_scale(kappa,Nf)
    tmp = 2.0_wp*pi*Nf*bessel_In(0,skappa)
    tmp = 1.0_wp/tmp

    tmp2 = exp(skappa*cos(theta-mu))
    p = p+tmp*tmp2
    !> first derivative
    tmp3 = tmp*skappa*sin(mu-theta)*exp(skappa*cos(mu-theta))
    tmpdpdt = tmpdpdt+tmp3
    !> second derivative
    tmp3 = tmp*skappa*exp(skappa*cos(mu-theta))*(-cos(mu-theta) &
    &      +skappa*(sin(mu-theta)**2))
    tmpdpdt2 = tmpdpdt2+tmp3

    if (present(dpdt)) dpdt = dpdt+tmpdpdt
    if (present(dpdt2)) dpdt2 = dpdt2+tmpdpdt2

  end subroutine vonMises_fraction

!========================================================================================!
  subroutine vonMises_avg(theta,kappa,mu,p,dpdt,dpdt2)
!********************************************************************
!* The subroutine computes a smooth probability density
!* ρ(θ|μ,κ) as the average of N von-Mises kernels, which are
!* continuous probability distributions around a specific angle.
!* It takes into account both the central location parameter,
!* i.e. the N observations μ of the angle θ, and the concentration
!* parameter κ.
!* The latter is obtained from the "von Mises-scale plug-in rule".
!* ρ(θ|μ,κ) is given by
!*
!*                1       N
!* ρ(θ|μ,κ) = ─────────  ∑ exp(κ' cos(θ -  μᵢ))
!*            2πNI₀(κ')   i
!*
!* where I₀ is the modiﬁed Bessel function of order 0 and κ' is
!* computed from κ via the "scale plug-in rule" (see below).
!*******************************************************************
    implicit none
    !> INPUT
    real(wp),intent(in) :: theta
    real(wp),intent(in) :: kappa
    real(wp),intent(in) :: mu(:)
    !> OUTPUT
    real(wp),intent(out) :: p
    real(wp),intent(out),optional :: dpdt  !> first derivative wrt θ
    real(wp),intent(out),optional :: dpdt2 !> second derivcative wrt θ
    !> LOCAL
    integer :: N
    integer :: i,j
    real(wp) :: Nf,skappa
    real(wp) :: tmp,tmp2,tmp3
    real(wp) :: tmpdpdt,tmpdpdt2
    real(wp) :: t1,t2,t3

    N = size(mu,1)
    Nf = float(N)

    p = 0.0_wp
    if (present(dpdt)) dpdt = 0.0_wp
    if (present(dpdt2)) dpdt2 = 0.0_wp
    tmpdpdt = 0.0_wp
    tmpdpdt2 = 0.0_wp

    skappa = vonMises_scale(kappa,Nf)
    tmp = 2.0_wp*pi*Nf*bessel_I0(skappa)
    tmp = 1.0_wp/tmp

    do i = 1,N
      tmp2 = exp(skappa*cos(theta-mu(i)))
      p = p+tmp*tmp2
      !> first derivative
      t1 = cos(mu(i)-theta)
      t2 = sin(mu(i)-theta)
      tmp3 = tmp*skappa*t2*exp(skappa*t1)
      tmpdpdt = tmpdpdt+tmp3
      !> second derivative
      tmp3 = tmp*skappa*exp(skappa*t1)*(-t1+skappa*(t2*t2))
      tmpdpdt2 = tmpdpdt2+tmp3
    end do

    if (present(dpdt)) dpdt = tmpdpdt
    if (present(dpdt2)) dpdt2 = tmpdpdt2

  end subroutine vonMises_avg

!========================================================================================!
  function vonMises_scale(kappa,Nf) result(skappa)
!*****************************************************************
!* Computes the modified concentration parameter κ' from κ
!* and the number of data points N (here expected as a float Nf)
!* via the "von Mises scale plug-in rule".
!* See https://doi.org/10.1016/j.csda.2007.11.003
!*
!*  κ' = [3*N*κ²*I₂(2κ)*{4π^(1/2)*I₀(κ)²}⁻¹]^(2/5)
!*
!*****************************************************************
    implicit none
    !> INPUT
    real(wp),intent(in) :: kappa
    real(wp),intent(in) :: Nf
    !> OUTPUT
    real(wp) :: skappa
    !> LOCAL
    real(wp) :: I0,I2
    real(wp) :: t1

    !> modified Bessel
    I0 = bessel_In(0,kappa)
    I2 = bessel_In(2,2.0_wp*kappa)

    !> calculate κ'
    skappa = 3.0_wp*Nf*(kappa**2)*I2
    t1 = 4.0_wp*sqrt(pi)*(I0**2)
    skappa = skappa*(1.0_wp/t1)
    skappa = skappa**(2.0_wp/5.0_wp)
  end function vonMises_scale


!========================================================================================!

  subroutine vonMises_plot(kappa,mu,n)
!*******************************************************
!* Plot an example von Mises distribution with n points
!*******************************************************  
  implicit none
  integer,intent(in) :: n
  real(wp),intent(in) :: mu(:)
  real(wp),intent(in) ::  kappa
  real(wp) :: theta,p,dpdt,dpdt2,numdpdt,pp,pm,tmp
  real(wp) :: d
  integer :: i,j,k,l,NM
  integer :: ich

  d = (2.0_wp*pi)/float(n)
  theta = 0.0_wp
  open (newunit=ich,file='vonmises.txt')
  do i = 1,n
    call vonMises(theta,kappa,mu,p,dpdt,dpdt2)
    write (ich,'(2f16.8)') theta,p
    theta = theta+d
  end do
  close (ich)

end subroutine vonMises_plot

  

!========================================================================================!
!========================================================================================!
end module probabilities_module
