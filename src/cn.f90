! This file is part of crest.
!
! Copyright (C) 2024 Philipp Pracht
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

module crest_cn_module
  use crest_parameters
  use miscdata,only:RCOV,PAULING_EN
  implicit none
!=========================================================================================!
!>--- private module variables and parameters
  private
  real(wp),parameter :: k1 = 16.0_wp
  real(wp),parameter :: sqrtpi = sqrt(pi)

  !> CN type enumerator
  type,private:: enum_cntype
    integer :: exp = 1
    integer :: erf = 2
    integer :: erf_en = 3
    integer :: gfn = 4
    integer :: d4 = 5
  end type enum_cntype
  type(enum_cntype),parameter,private :: cn_ver = enum_cntype()

  public :: calculate_CN
  public :: calc_ncoord
  interface calc_ncoord
    module procedure cn_ncoord
    module procedure cn_xcoord
    module procedure cn_xcoord2
  end interface calc_ncoord
  public :: cn_ycoord
  public :: cn_ycoord2

!========================================================================================!
!========================================================================================!
contains  !> MODULE PROCEDURES START HERE
!========================================================================================!
!========================================================================================!

  subroutine calculate_CN(nat,at,xyz,cn,cnthr,cntype,dcndr,bond)
!*********************************************************
!* Universal CN calculator with several optional settings
!* for customisation.
!* The routine requires at least the following arguments:
!*   nat - number of atoms
!*   at  - atomic numbers
!*   xyz - coordinates in Bohr
!* Optional arguments:
!*   cnthr  - pair distance cutoff in Bohr
!*   cntype - string to select CN type (exp,erf,erf_en)
!*   dcndr  - optional output cn derivatives
!*   bond   - optional output "bond" connectivity matrix
!*********************************************************
    implicit none
    !> INPUT
    integer,intent(in) :: nat
    integer,intent(in) :: at(nat)
    real(wp),intent(in) :: xyz(3,nat)
    !> OUTPUT
    real(wp),intent(out) :: cn(nat)
    !> OPTIONAL
    real(wp),intent(in),optional :: cnthr
    character(len=*),intent(in),optional :: cntype
    real(wp),intent(out),optional :: dcndr(:,:,:)
    real(wp),intent(out),allocatable,optional :: bond(:,:)
    !> LOCAL
    real(wp) :: cn_thr,cn_direct
    integer :: cn_type,ati,atj
    real(wp) :: rcovi,rcovj,rij(3),r2,r,rco
    real(wp) :: damp,ddamp(3),den
    integer :: i,j,k,l,iat,jat
    logical :: deriv,getbond

!>--- check options and defaults
    if (present(cnthr)) then
      cn_thr = cnthr**2  !> note the value will be squared!
    else
      cn_thr = 625.0_wp  !> 25.0^2 as in tblite
    end if

    cn_direct = 1.0_wp
    if (present(cntype)) then
      select case (cntype)
      case ('exp')
        cn_type = cn_ver%exp
      case ('erf')
        cn_type = cn_ver%erf
      case ('erf_en') !> DO NOT USE, SPECIAL FOR CEH
        cn_type = cn_ver%erf_en
        cn_direct = -1.0_wp
      case ('gfn')
        cn_type = cn_ver%gfn
      case ('d4','cov')
        cn_type = cn_ver%d4
      case default
        cn_type = cn_ver%exp
      end select
    else
      cn_type = cn_ver%exp
    end if

    if (present(dcndr)) then
      deriv = .true.
      dcndr(:,:,:) = 0.0_wp
    else
      deriv = .false.
    end if

    if(present(bond))then
       getbond=.true.
       allocate(bond(nat,nat),source=0.0_wp)
    else
       getbond=.false.
    endif


    cn(:) = 0.0_wp
!>--- actual calculation

    do i = 1,nat
      ati = at(i)
      rcovi = RCOV(ati)

      do j = 1,i-1
        rij(:) = xyz(:,i)-xyz(:,j)
        r2 = sum(rij**2)
!>--- cycle cutoff
        if (r2 > cn_thr) cycle
        r = sqrt(r2)
        atj = at(j)
        rcovj = RCOV(atj)
        rco = rcovi+rcovj

!>--- select the correct CN version
        select case (cn_type)
        case (cn_ver%exp)
          damp = cn_damp_exp(rco,r)
          if (deriv) ddamp(:) = dcn_damp_exp(rco,r)*(rij(:)/r)
        case (cn_ver%erf)
          damp = cn_damp_erf(rco,r)
          if (deriv) ddamp(:) = dcn_damp_erf(rco,r)*(rij(:)/r)
        case (cn_ver%erf_en) !> DO NOT USE, SPECIAL FOR CEH
          den = PAULING_EN(atj)-PAULING_EN(ati)
          damp = cn_damp_erf_en(rco,r,den)
          if (deriv) ddamp(:) = dcn_damp_erf_en(rco,r,den)*(rij(:)/r)
        case (cn_ver%d4)
          den = PAULING_EN(atj)-PAULING_EN(ati)
          damp = cn_damp_d4(rco,r,den)
          if (deriv) ddamp(:) = dcn_damp_d4(rco,r,den)*(rij(:)/r)
        case (cn_ver%gfn)
          damp = cn_damp_gfn(rco,r)
          if (deriv) ddamp(:) = dcn_damp_gfn(rco,r)*(rij(:)/r)
        end select
        
        if(getbond)then
           bond(j,i) = damp
           bond(i,j) = damp 
        endif

        cn(i) = cn(i)+damp
        if (i /= j) then
          cn(j) = cn(j)+damp*cn_direct
        end if
        if (deriv) then
          dcndr(:,i,i) = dcndr(:,i,i)+ddamp(:)
          dcndr(:,i,j) = dcndr(:,i,j)+ddamp(:)*cn_direct
          dcndr(:,j,i) = dcndr(:,j,i)-ddamp(:)
          dcndr(:,j,j) = dcndr(:,j,j)-ddamp(:)*cn_direct
        end if

      end do
    end do

  end subroutine calculate_CN

!=========================================================================================!
!> old ncoord-based routines for compatibility
!=========================================================================================!

  subroutine cn_ncoord(nat,rcovin,iz,xyz,cn,cn_thr)
!****************************************************************
!* Basic version taking covalent radii and cn threshold as input
!*****************************************************************
    implicit none
    !> INPUT
    integer,intent(in) :: nat
    real(wp),intent(in) :: rcovin(*)
    integer,intent(in) :: iz(nat)
    real(wp),intent(in) :: xyz(3,nat)
    real(wp),intent(in) :: cn_thr
    !> OUTPUT
    real(wp),intent(out) :: cn(nat)
    integer :: iat,i
    real(wp) :: r,damp,xn,rr,rco,r2,rcovi,rcovj

    do i = 1,nat
      xn = 0.0d0
      rcovi = rcovin(iz(i))
      do iat = 1,nat
        if (iat .ne. i) then
          r2 = cn_help_rdist(nat,xyz,iat,i)
          if (r2 .gt. cn_thr) cycle
          r = sqrt(r2)
          rcovj = rcovin(iz(iat))
!> covalent distance in Bohr
          rco = rcovi+rcovj
          rr = rco/r
!> counting function exponential has a better long-range behavior than MHGs inverse damping
          damp = 1.d0/(1.d0+exp(-k1*(rr-1.0d0)))
          xn = xn+damp
        end if
      end do
      cn(i) = xn
    end do
  end subroutine cn_ncoord

!=========================================================================================!

  subroutine cn_ycoord(nat,rcovin,iz,xyz,cn,cn_thr)
!**********************************************************
!* modified verison that multiplies the damping factor
!* by the atomic number. Created for the CONFCROSS routine
!*
!* IMPORTANT: THIS IS NOT AN ACTUAL CN!
!**********************************************************
    implicit none
    !> INPUT
    integer,intent(in) :: nat
    real(wp),intent(in) :: rcovin(*)
    integer,intent(in) :: iz(nat)
    real(wp),intent(in) :: xyz(3,nat)
    real(wp),intent(in) :: cn_thr
    !> OUTPUT
    real(wp),intent(out) :: cn(nat)
    !> LOCAL
    integer  :: iat,i
    real(wp) :: r,damp,xn,rr,rco,r2,rcovi,rcovj

    do i = 1,nat
      xn = 0.0d0
      rcovi = rcovin(iz(i))
      do iat = 1,nat
        if (iat .ne. i) then
          r2 = cn_help_rdist(nat,xyz,iat,i)
          if (r2 .gt. cn_thr) cycle
          r = sqrt(r2)
          rcovj = rcovin(iz(iat))
!> covalent distance in Bohr
          rco = rcovi+rcovj
          rr = rco/r
!> counting function exponential has a better long-range behavior than MHGs inverse damping
          damp = 1.d0/(1.d0+exp(-k1*(rr-1.0d0)))
          xn = xn+damp*real(iz(iat))
        end if
      end do
      cn(i) = xn
    end do
  end subroutine cn_ycoord

  subroutine cn_ycoord2(nat,rcovin,iz,xyz,cn,cn_thr,cthr,clash)
!**********************************************************
!* variation to the ycoord version that checks for clashes
!* if a clash is detected, the routine exits early!
!* Note the "CN" is NOT an output, but an INPUT argument,
!* to be used in tandem with ycoord
!**********************************************************
    implicit none
    !> INPUT
    integer,intent(in) :: nat
    real(wp),intent(in) :: rcovin(*)
    integer,intent(in) :: iz(nat)
    real(wp),intent(in) :: xyz(3,nat)
    real(wp),intent(in) :: cn_thr
    real(wp),intent(in) :: cthr
    real(wp),intent(in) :: cn(nat)
    !> OUTPUT
    logical,intent(out)  :: clash
    !> LOCAL
    integer  :: iat,i
    real(wp) :: r,damp,xn,rr,rco,r2,rcovi,rcovj

    clash = .false.
    do i = 1,nat
      xn = 0.0d0
      rcovi = rcovin(iz(i))
      do iat = 1,nat
        if (iat .ne. i) then
          r2 = cn_help_rdist(nat,xyz,iat,i)
          if (r2 .gt. cn_thr) cycle
          r = sqrt(r2)
          rcovj = rcovin(iz(iat))
!> covalent distance in Bohr
          rco = rcovi+rcovj
          rr = rco/r
!> counting function exponential has a better long-range behavior than MHGs inverse damping
          damp = 1.d0/(1.d0+exp(-k1*(rr-1.0d0)))
          xn = xn+damp*iz(iat)
        end if
      end do
      if (abs(cn(i)-xn) .gt. cthr) then                  !   clash check
        clash = .true.
        return
      end if
    end do
  end subroutine cn_ycoord2

!=========================================================================================!

  subroutine cn_xcoord(nat,iz,xyz,cn,bond)
!*********************************************
!* simple version of CN calculation routine
!* using a predefined CN cutoff of 400.0
!* Also outputs a "bond" matrix
!********************************************
    implicit none
    !> INPUT
    integer,intent(in) :: nat
    integer,intent(in) :: iz(nat)
    real(wp),intent(in) :: xyz(3,nat)
    !> OUTPUT
    real(wp),intent(out) :: cn(nat)
    real(wp),intent(out) :: bond(nat,nat)
    !> LOCAL
    real(wp) :: cn_thr,dx,dy,dz,r,damp
    real(wp) :: xn,rr,rco,r2,rcovi,rcovj
    integer  :: iat,i
    !> take a fixed cn_thr
    !> and the global
    cn_thr = 400.0d0
    call cn_xcoord2(nat,iz,xyz,RCOV,cn,cn_thr,bond)
  end subroutine cn_xcoord

  subroutine cn_xcoord2(nat,iz,xyz,rcovin,cn,cn_thr,bond)
!***********************************************
!* A clean version using input arguments for
!* rcov, cn_thr and outputting a "bond" matrix
!***********************************************
    implicit none
    !> INPUT
    integer,intent(in) :: nat
    integer,intent(in) :: iz(nat)
    real(wp),intent(in) :: xyz(3,nat)
    real(wp),intent(in)  :: cn_thr
    real(wp),intent(in)  :: rcovin(*)
    !> OUTPUT
    real(wp),intent(out) :: cn(nat)
    real(wp),intent(out) :: bond(nat,nat)
    !> LOCAL
    integer :: iat,i
    real(wp) :: r,damp,xn,rr,rco,r2,rcovi,rcovj
    bond = 0.0d0
    cn = 0.0d0
    do i = 1,nat
      xn = 0.0d0
      rcovi = rcovin(iz(i))
      do iat = 1,nat
        if (iat .ne. i) then
          r2 = cn_help_rdist(nat,xyz,iat,i)
          if (r2 .gt. cn_thr) cycle
          r = sqrt(r2)
          rcovj = rcovin(iz(iat))
!> covalent distance in Bohr
          rco = (rcovi+rcovj)
          rr = rco/r
!> counting function exponential has a better long-range behavior than MHGs inverse damping
          damp = 1.d0/(1.d0+exp(-k1*(rr-1.0d0)))
          bond(iat,i) = damp
          xn = xn+damp
        end if
      end do
      cn(i) = xn
    end do
    return
  end subroutine cn_xcoord2

!========================================================================================!
!> Helper functions and damping versions
!========================================================================================!

  function cn_help_rdist(nat,xyz,iat,jat) result(r2)
    implicit none
    integer,intent(in) :: nat,iat,jat
    real(wp),intent(in) :: xyz(3,nat)
    real(wp) :: dx,dy,dz,r2
    dx = xyz(1,iat)-xyz(1,jat)
    dy = xyz(2,iat)-xyz(2,jat)
    dz = xyz(3,iat)-xyz(3,jat)
    r2 = dx*dx+dy*dy+dz*dz
  end function cn_help_rdist

  function cn_damp_exp(rco,r) result(damp)
!*****************************************
!* classic exponential CN damping factor
!*****************************************
    implicit none
    real(wp) :: damp
    real(wp),intent(in) :: rco,r
    real(wp),parameter :: kcn = 16.0_wp
    damp = 1.0_wp/(1.0_wp+exp(-kcn*(rco/r-1.0_wp)))
  end function cn_damp_exp
  function dcn_damp_exp(rco,r) result(ddamp)
!**********************************************
!* derivative of exponential CN damping factor
!***********************************************
    implicit none
    real(wp) :: ddamp
    real(wp),intent(in) :: rco,r
    real(wp) :: tmp
    real(wp),parameter :: kcn = 16.0_wp
    tmp = exp(-kcn*(rco/r-1.0_wp))
    ddamp = (-kcn*rco*tmp)/(r**2*((tmp+1.0_wp)**2))
  end function dcn_damp_exp

  function cn_damp_erf(rco,r) result(damp)
!*************************
!* erf-CN damping factor
!*************************
    implicit none
    real(wp) :: damp
    real(wp),intent(in) :: rco,r
    real(wp),parameter :: kcn = 2.60_wp
    damp = 0.5_wp*(1.0_wp+erf(-kcn*(r-rco)/rco))
  end function cn_damp_erf
  function dcn_damp_erf(rco,r) result(ddamp)
!**************************************
!* derivative of erf-CN damping factor
!**************************************
    implicit none
    real(wp) :: ddamp
    real(wp),intent(in) :: rco,r
    real(wp) :: tmp
    real(wp),parameter :: kcn = 2.60_wp
    tmp = exp(-(kcn*(r-rco)/rco)**2)
    ddamp = (-kcn*tmp)/(rco*sqrtpi)
  end function dcn_damp_erf

  function cn_damp_gfn(rco,r) result(damp)
!***************************************
!* double-exponential CN damping factor
!***************************************
    implicit none
    real(wp) :: damp
    real(wp),intent(in) :: rco,r
    real(wp) :: dampa,dampb,rcob
    !> Steepness of first counting function
    real(wp),parameter :: ka = 10.0_wp
    !> Steepness of second counting function
    real(wp),parameter :: kb = 20.0_wp
    !> Offset of the second counting function
    real(wp),parameter :: r_shift = 2.0_wp
    dampa = 1.0_wp/(1.0_wp+exp(-ka*(rco/r-1.0_wp)))
    rcob = rco+r_shift
    dampb = 1.0_wp/(1.0_wp+exp(-kb*(rcob/r-1.0_wp)))
    damp = dampa*dampb
  end function cn_damp_gfn
  function dcn_damp_gfn(rco,r) result(ddamp)
!*****************************************************
!* derivative of double-exponential CN damping factor
!*****************************************************
    implicit none
    real(wp) :: ddamp
    real(wp),intent(in) :: rco,r
    real(wp) :: tmp
    real(wp) :: dampa,dampb,rcob
    real(wp) :: ddampa,ddampb
    !> Steepness of first counting function
    real(wp),parameter :: ka = 10.0_wp
    !> Steepness of second counting function
    real(wp),parameter :: kb = 20.0_wp
    !> Offset of the second counting function
    real(wp),parameter :: r_shift = 2.0_wp
    dampa = 1.0_wp/(1.0_wp+exp(-ka*(rco/r-1.0_wp)))
    rcob = rco+r_shift
    dampb = 1.0_wp/(1.0_wp+exp(-kb*(rcob/r-1.0_wp)))
    tmp = exp(-ka*(rco/r-1._wp))
    ddampa = (-ka*rco*tmp)/(r**2*((tmp+1.0_wp)**2))
    tmp = exp(-kb*(rcob/r-1._wp))
    ddampb = (-kb*rcob*tmp)/(r**2*((tmp+1.0_wp)**2))
    ddamp = ddampa*dampb+dampa*ddampb
  end function dcn_damp_gfn

  function cn_damp_erf_en(rco,r,den) result(damp)
!*********************************************************
!* EN-dependent erf-CN damping factor (developed for CEH)
!*********************************************************
    implicit none
    real(wp) :: damp
    real(wp),intent(in) :: rco,r,den
    damp = den*cn_damp_erf(rco,r)
  end function cn_damp_erf_en
  function dcn_damp_erf_en(rco,r,den) result(ddamp)
!***********************************************************************
!* derivative of EN-dependent erf-CN damping factor (developed for CEH)
!***********************************************************************
    implicit none
    real(wp) :: ddamp
    real(wp),intent(in) :: rco,r,den
    ddamp = den*dcn_damp_erf(rco,r)
  end function dcn_damp_erf_en

  function cn_damp_d4(rco,r,den) result(damp)
!*********************************************
!* EN-dependent erf-CN damping factor from D4
!*********************************************
    implicit none
    real(wp) :: damp
    real(wp),intent(in) :: rco,r,den
    real(wp) :: tmp,tmperf
    real(wp),parameter :: k4 = 4.10451_wp
    real(wp),parameter :: k5 = 19.08857_wp
    real(wp),parameter :: k6 = 2*11.28174_wp**2
    real(wp),parameter :: kcn = 7.50_wp
    tmp = k4*exp(-((abs(den)+k5)**2)/k6)
    tmperf = 0.5_wp*(1.0_wp+erf(-kcn*(r-rco)/rco))
    damp = tmp*tmperf
!    damp = tmp*cn_damp_erf(rco,r)
  end function cn_damp_d4
  function dcn_damp_d4(rco,r,den) result(ddamp)
!***********************************************************
!* derivative of EN-dependent erf-CN damping factor from D4
!***********************************************************
    implicit none
    real(wp) :: ddamp
    real(wp),intent(in) :: rco,r,den
    real(wp) :: tmp,dtmperf,tmp2
    real(wp),parameter :: k4 = 4.10451_wp
    real(wp),parameter :: k5 = 19.08857_wp
    real(wp),parameter :: k6 = 2*11.28174_wp**2
    real(wp),parameter :: kcn = 7.50_wp
    tmp = k4*exp(-(abs(den)+k5)**2/k6)
    tmp2 = exp(-(kcn*(r-rco)/rco)**2)
    dtmperf = (-kcn*tmp2)/(rco*sqrtpi)
    ddamp = tmp*dtmperf
  end function dcn_damp_d4

!========================================================================================!
!========================================================================================!
end module crest_cn_module
