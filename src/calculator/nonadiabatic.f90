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
  interface engrad_mean
    module procedure :: engrad_mean_2
    module procedure :: engrad_mean_multi
    module procedure :: engrad_mean_weight
  end interface engrad_mean

  public :: calc_nonadiabatic_constraint

!========================================================================================!
!========================================================================================!
contains  !>--- Module routines start here
!========================================================================================!
!========================================================================================!

!========================================================================================!
!> subrotuine engrad_mean
!> given energies and the respective gradients,
!> the routine returns the arithmetic mean between
!> those.
!>---------------------------------------------------
  subroutine engrad_mean_2(nat,e1,e2,grd1,grd2,em,grdm)
    implicit none
    integer,intent(in) :: nat
    real(wp),intent(in) :: e1,e2
    real(wp),intent(in) :: grd1(3,nat)
    real(wp),intent(in) :: grd2(3,nat)
    real(wp),intent(out) :: em
    real(wp),intent(out) :: grdm(3,nat)

    em = 0.5_wp*(e1+e2)

    grdm = 0.5_wp*(grd1+grd2)

    return
  end subroutine engrad_mean_2
!===============!
  subroutine engrad_mean_multi(nat,nlev,e,grd,em,grdm)
    implicit none
    integer,intent(in) :: nat
    integer,intent(in) :: nlev
    real(wp),intent(in) :: e(nlev)
    real(wp),intent(in) :: grd(3,nat,nlev)
    real(wp),intent(out) :: em
    real(wp),intent(out) :: grdm(3,nat)
    integer :: l,i
    real(wp) :: dum

    em = 0.0_wp
    grdm = 0.0_wp

    l = nlev
    dum = 1.0_wp/float(l)

    do i = 1,l
      em = em+dum*e(i)
      grdm = grdm+dum*grd(:,:,i)
    end do
    return
  end subroutine engrad_mean_multi
!===============!
  subroutine engrad_mean_weight(nat,nlev,weights,e,grd,em,grdm)
    implicit none
    integer,intent(in) :: nat
    integer,intent(in) :: nlev
    real(wp),intent(in) :: weights(nlev)
    real(wp),intent(in) :: e(nlev)
    real(wp),intent(in) :: grd(3,nat,nlev)
    real(wp),intent(out) :: em
    real(wp),intent(out) :: grdm(3,nat)
    integer :: l,i
    real(wp) :: dum

    em = 0.0_wp
    grdm = 0.0_wp
    l = nlev
    do i = 1,l
      dum = weights(i)
      em = em+dum*e(i)
      grdm = grdm+dum*grd(:,:,i)
    end do
    return
  end subroutine engrad_mean_weight

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
    c = constr%fc(3)

    select case (constr%type)
    case (na_gapdiff)
      call gapdiff_constraint(n,sigm,alph,energies(1),energies(2), &
      &    grds(:,:,1),grds(:,:,2),efix,gfix)
    case (na_gapdiff2)
      !call gapdiff_constraint2(n,sigm,alph,c,energies(1),energies(2), &
      !&    grds(:,:,1),grds(:,:,2),efix,gfix)
      call gapdiff_constraint2_pairwise(n,sigm,alph,c,d,energies, &
      &    grds,efix,gfix)
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

    gap = e2-e1

    efix = sigm*(gap**2/(gap+alph))
    efixgrd = sigm*((gap**2+2.0_wp*alph*gap)/(gap+alph)**2)

    do i = 1,nat
      gfix(:,i) = efixgrd*(grd2(:,i)-grd1(:,i))
    end do

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

    gap = e2-e1

    !> Gaussian function (note: abs(gap) = sqrt(gap)²)
    b = autoev
    gss = sigm*(exp(-b*abs(gap))+c)
    gssgrd = -b*sigm*exp(-b*abs(gap))

    !> Gap-based potential
    vfix = (gap**2/(abs(gap)+alph))
    vgrd = ((gap**2+2.0_wp*alph*abs(gap))/(abs(gap)+alph)**2)

    !> factor to assign the correct prefactor to gradient
    absfct = gap/abs(gap)

    !> construct bias energy and gradient
    efix = gss*vfix
    efixgrd = gss*vgrd+gssgrd*vfix

    do i = 1,nat
      gfix(:,i) = efixgrd*absfct*(grd2(:,i)-grd1(:,i))
    end do

    return
  end subroutine gapdiff_constraint2

!========================================================================================!
!> subroutine gapdiff_constraint2
!> construct an constraint that minimizes minimizes the
!> gap between two given energy surfaces
!> new version with Gaussian potential as prefactor
!> Pairwise version
!>-------------------------------------------------------------------
  subroutine gapdiff_constraint2_pairwise(nat,sigm,alph,c,d,energies,grads,efix,gfix)
    implicit none
    integer,intent(in) :: nat,d
    real(wp),intent(in) :: alph,sigm,c
    real(wp),intent(in) :: energies(d)
    real(wp),intent(in) :: grads(3,nat,d)
    real(wp),intent(out) :: efix
    real(wp),intent(out) :: gfix(3,nat)
    integer :: i,j,ei,ej
    real(wp) :: b
    real(wp) :: gap,efixgrd,gss,gssgrd,absfct
    real(wp) :: vfix,vgrd
    real(wp),parameter :: autoev = 27.2114_wp

    efix = 0.0_wp
    gfix = 0.0_wp
    b = autoev

    do ei = 1,d
      do ej = 1,ei-1
        gap = energies(ei) - energies(ej)
        !> Gaussian function (note: abs(gap) = sqrt(gap)²)
        gss = sigm*(exp(-b*abs(gap))+c)
        gssgrd = -b*sigm*exp(-b*abs(gap))

        !> Gap-based potential
        vfix = (gap**2/(abs(gap)+alph))
        vgrd = ((gap**2+2.0_wp*alph*abs(gap))/(abs(gap)+alph)**2)

        !> factor to assign the correct prefactor to gradient
        absfct = gap/abs(gap)

        !> construct bias energy and gradient
        efix = efix + gss*vfix
        efixgrd = gss*vgrd+gssgrd*vfix
        do i = 1,nat
          gfix(:,i) = gfix(:,i) + efixgrd*absfct*(grads(:,i,ei)-grads(:,i,ej))
        end do

      end do
    end do

    return
  end subroutine gapdiff_constraint2_pairwise

end module nonadiabatic_module
