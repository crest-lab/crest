!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2022 Philipp Pracht
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

module nonadiabatic_module

  use iso_fortran_env,only:wp => real64
  use constraints
  implicit none

!=========================================================================================!
  !--- private module variables and parameters
  private
  !--- some constants and name mappings
  real(wp),parameter :: bohr = 0.52917726_wp
  real(wp),parameter :: autokcal = 627.509541_wp

  public :: engrad_mean
  public :: calc_nonadiabatic_constraint

contains

!========================================================================================!
!> subrotuine engrad_mean
!> given two energies and the respective gradients,
!> the routine returns the arithmetic mean between
!> those two.
!>---------------------------------------------------
  subroutine engrad_mean(nat,e1,e2,grd1,grd2,em,grdm)
    implicit none
    integer,intent(in) :: nat
    real(wp),intent(in) :: e1,e2
    real(wp),intent(in) :: grd1(3,nat)
    real(wp),intent(in) :: grd2(3,nat)
    real(wp),intent(out) :: em
    real(wp),intent(out) :: grdm(3,nat)

    em = 0.5_wp * (e1 + e2)

    grdm = 0.5_wp * (grd1 + grd2)

    return
  end subroutine engrad_mean



!========================================================================================!
!> subroutine calc_nonadiabatic_constraint
!> inteface to constraints of the non-adiabatic surfaces
  subroutine calc_nonadiabatic_constraint(n,constr,d,energies,grds,efix,gfix)
    implicit none
    integer,intent(in) :: n
    type(constraint) :: constr
    integer,intent(in) :: d
    real(wp),intent(in) :: energies(d)
    real(wp),intent(in) :: grds(3,n,d)
    real(wp),intent(out) :: efix
    real(wp),intent(out) :: gfix(3,n)

    real(wp) :: sigm,alph,c

    efix = 0.0_wp
    gfix = 0.0_wp

    sigm = constr%fc(1)
    alph = constr%fc(2)
    c    = constr%fc(3)

    select case (constr%type)
    case (na_gapdiff)
      call gapdiff_constraint(n,sigm,alph,energies(1),energies(2), &
      &    grds(:,:,1),grds(:,:,2),efix,gfix)
    case (na_gapdiff2)
      call gapdiff_constraint2(n,sigm,alph,c,energies(1),energies(2), &
      &    grds(:,:,1),grds(:,:,2),efix,gfix)
    case default
      return
    end select
 
    !write(*,*) energies(2)-energies(1), efix, gfix(1,1)

    return
  end subroutine calc_nonadiabatic_constraint

!========================================================================================!
!> subroutine gapdiff_constraint
!> construct an constraint that minimizes minimizes the
!> gap between two given energy surfaces
!> See 
!>     B. Levine, J. Coe, T. Martinez
!>     J.Phys.Chem.B, 2008, 112, 405-413
!>-------------------------------------------------------------------
  subroutine gapdiff_constraint(nat,sigm,alph,e1,e2,grd1,grd2,efix,gfix)
    implicit none
    integer,intent(in) :: nat
    real(wp),intent(in) :: e1,e2,alph,sigm
    real(wp),intent(in) :: grd1(3,nat),grd2(3,nat)
    real(wp),intent(out) :: efix
    real(wp),intent(out) :: gfix(3,nat)
    integer :: i,j
    real(wp) :: gap,efixgrd

    efix = 0.0_wp
    gfix = 0.0_wp

    gap = e2 - e1
    
    efix = sigm * ( gap**2 / (gap + alph))
    efixgrd = sigm * ( (gap**2 + 2.0_wp*alph*gap)  / (gap + alph)**2)    

    do i=1,nat
      gfix(:,i) = efixgrd * ( grd2(:,i) - grd1(:,i))
    enddo

    return
  end subroutine gapdiff_constraint


!========================================================================================!
!> subroutine gapdiff_constraint2
!> construct an constraint that minimizes minimizes the
!> gap between two given energy surfaces
!> new version with Gaussian potential as prefactor 
!>-------------------------------------------------------------------
  subroutine gapdiff_constraint2(nat,sigm,alph,c,e1,e2,grd1,grd2,efix,gfix)
    implicit none
    integer,intent(in) :: nat
    real(wp),intent(in) :: e1,e2,alph,sigm,c
    real(wp),intent(in) :: grd1(3,nat),grd2(3,nat)
    real(wp),intent(out) :: efix
    real(wp),intent(out) :: gfix(3,nat)
    integer :: i,j
    real(wp) :: b
    real(wp) :: gap,efixgrd,gss,gssgrd,absfct
    real(wp) :: vfix,vgrd
    real(wp),parameter :: autoev = 27.2114_wp

    efix = 0.0_wp
    gfix = 0.0_wp

    gap = e2 - e1

    !> Gaussian function (note: abs(gap) = sqrt(gap)²)
    b = autoev
    gss = sigm * (exp(-b*abs(gap)) + c)
    gssgrd = -b*sigm*exp(-b*abs(gap))

    !> Gap-based potential
    vfix = ( gap**2 / (abs(gap) + alph))
    vgrd = ( (gap**2 + 2.0_wp*alph*abs(gap))  / (abs(gap) + alph)**2)

    !> factor to assign the correct prefactor to gradient
    absfct = gap/abs(gap)

    !> construct bias energy and gradient
    efix = gss * vfix
    efixgrd = gss*vgrd + gssgrd*vfix

    do i=1,nat
      gfix(:,i) = efixgrd * absfct * ( grd2(:,i) - grd1(:,i))
    enddo

    return
  end subroutine gapdiff_constraint2

end module nonadiabatic_module
