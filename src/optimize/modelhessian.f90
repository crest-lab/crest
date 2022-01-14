!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2021 - 2022 Philipp Pracht
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
module modelhessian_module
  use iso_fortran_env,only:wp => real64
  use calc_type
  use calc_module,only:constrhess
  implicit none

!> a modelhessian type to save settings
  type :: mhparam
    integer :: model = 0       !> model hessian selection
    real(wp) :: s6 = 20.0_wp   !> dispersion scaling
    real(wp) :: rcut = 70.0_wp !> cutoff parameter
    !> force constants
    real(wp) :: kr = 0.4000_wp
    real(wp) :: kf = 0.1300_wp
    real(wp) :: kt = 0.0075_wp
    real(wp) :: ko = 0.0000_wp
    real(wp) :: kd = 0.0000_wp
    real(wp) :: kq = 0.0000_wp
  end type mhparam

!> Parameters & constants
  real(wp),parameter :: bohr = 0.52917726_wp
  real(wp),parameter :: pi = 3.141592653589793_wp
  real(wp),parameter :: Zero = 0.0_wp
  real(wp),parameter :: One = 1.0_wp
  real(wp),parameter :: Two = 2.0_wp
  real(wp),parameter :: Three = 3.0_wp
  real(wp),parameter :: Four = 4.0_wp
  real(wp),parameter :: Five = 5.0_wp
  real(wp),parameter :: Six = 6.0_wp
  real(wp),parameter :: Seven = 7.0_wp
  real(wp),parameter :: Eight = 8.0_wp
  real(wp),parameter :: RNine = 9.0_wp
  real(wp),parameter :: Ten = 10.0_wp
  real(wp),parameter :: Half = 0.5_wp
  real(wp),parameter :: SqrtP2 = 0.8862269254527579_wp
  real(wp),parameter :: TwoP34 = 0.2519794355383808_wp
  real(wp),parameter :: TwoP54 = 5.914967172795612_wp
  real(wp),parameter :: One2C2 = 0.2662567690426443D-04

  !>  van-der-Waals radii used in the D2 model (NOTE: here not in a.u.)
  real(wp),parameter :: vander(86) = (/ &
  & 0.91_wp,0.92_wp, & ! H, He
  & 0.75_wp,1.28_wp,1.35_wp,1.32_wp,1.27_wp,1.22_wp,1.17_wp,1.13_wp, & ! Li-Ne
  & 1.04_wp,1.24_wp,1.49_wp,1.56_wp,1.55_wp,1.53_wp,1.49_wp,1.45_wp, & ! Na-Ar
  & 1.35_wp,1.34_wp, & ! K, Ca
  & 1.42_wp,1.42_wp,1.42_wp,1.42_wp,1.42_wp, & ! Sc-Zn
  & 1.42_wp,1.42_wp,1.42_wp,1.42_wp,1.42_wp, &
  & 1.50_wp,1.57_wp,1.60_wp,1.61_wp,1.59_wp,1.57_wp, & ! Ga-Kr
  & 1.48_wp,1.46_wp, & ! Rb, Sr
  & 1.49_wp,1.49_wp,1.49_wp,1.49_wp,1.49_wp, & ! Y-Cd
  & 1.49_wp,1.49_wp,1.49_wp,1.49_wp,1.49_wp, &
  & 1.52_wp,1.64_wp,1.71_wp,1.72_wp,1.72_wp,1.71_wp, & ! In-Xe
  & 2.00_wp,2.00_wp, &
  & 2.00_wp,2.00_wp,2.00_wp,2.00_wp,2.00_wp,2.00_wp,2.00_wp, & ! La-Yb
  & 2.00_wp,2.00_wp,2.00_wp,2.00_wp,2.00_wp,2.00_wp,2.00_wp, &
  & 2.00_wp,2.00_wp,2.00_wp,2.00_wp,2.00_wp, & ! Lu-Hg
  & 2.00_wp,2.00_wp,2.00_wp,2.00_wp,2.00_wp, &
  & 2.00_wp,2.00_wp,2.00_wp,2.00_wp,2.00_wp,2.00_wp/) ! Tl-Rn
  !>  C6 coefficients used in the D2 model
  real(wp),parameter :: c6(86) = (/ &
  & 0.14_wp,0.08_wp, & ! H,He
  & 1.61_wp,1.61_wp,3.13_wp,1.75_wp,1.23_wp,0.70_wp,0.75_wp,0.63_wp, &
  & 5.71_wp,5.71_wp,10.79_wp,9.23_wp,7.84_wp,5.57_wp,5.07_wp,4.61_wp, &
  & 10.80_wp,10.80_wp, & ! K,Ca
  & 10.80_wp,10.80_wp,10.80_wp,10.80_wp,10.80_wp, & ! Sc-Zn
  & 10.80_wp,10.80_wp,10.80_wp,10.80_wp,10.80_wp, &
  & 16.99_wp,17.10_wp,16.37_wp,12.64_wp,12.47_wp,12.01_wp, & ! Ga-Kr
  & 24.67_wp,24.67_wp, & ! Rb,Sr
  & 24.67_wp,24.67_wp,24.67_wp,24.67_wp,24.67_wp, & ! Y-Cd
  & 24.67_wp,24.67_wp,24.67_wp,24.67_wp,24.67_wp, &
  & 37.32_wp,38.71_wp,38.44_wp,31.74_wp,31.50_wp,29.99_wp, & ! In-Xe
  & 50.00_wp,50.00_wp, & ! Cs,Ba
  & 50.00_wp,50.00_wp,50.00_wp,50.00_wp,50.00_wp,50.00_wp,50.00_wp, & ! La-Yb
  & 50.00_wp,50.00_wp,50.00_wp,50.00_wp,50.00_wp,50.00_wp,50.00_wp, &
  & 50.00_wp,50.00_wp,50.00_wp,50.00_wp,50.00_wp, & ! Lu-Hg
  & 50.00_wp,50.00_wp,50.00_wp,50.00_wp,50.00_wp, &
  & 50.00_wp,50.00_wp,50.00_wp,50.00_wp,50.00_wp,50.00_wp/) ! Tl-Rn

  public :: modhes

contains
!========================================================================================!
!> subroutine modhes
!> create a model Hessian for a given molecule
!>
!> Input:
!>     natoms - number of atoms
!>       xyz  - Cartesian coordinates
!>        at  - atom types as integers
!>      modh  - model Hessian settings (see above)
!>      calc  - calculation settings (for constraints)
!>        pr  - printout selection
!>
!> Output:
!>      Hess  - the (packed) model Hessian
!>--------------------------------------------------------
  subroutine modhes(calc,modh,natoms,xyz,at,Hess,pr)
    implicit none
    type(calcdata),intent(in) :: calc
    type(mhparam),intent(in) :: modh
    logical,intent(in) :: pr
    integer :: i
    integer :: nhess
    integer,intent(in) :: natoms
    real(wp),intent(in) :: xyz(3,natoms)
    real(wp),intent(out) :: hess((natoms * 3) * ((natoms * 3) + 1) / 2)
    integer,intent(in) :: at(natoms)

!>  initialize
    nhess = 3 * natoms
    Hess = 0.0_wp

    select case (modh%model)
    case default  !case (p_modh_old)
      if (pr) write (*,'(a)') "Using Lindh-Hessian (1995)"
      call ddvopt(xyz,natoms,Hess,at,modh)
!> other model hessians currently not implemented
      !case (p_modh_lindh_d2)
      !  if (pr) write (*,'(a)') "Using Lindh-Hessian"
      !  call mh_lindh_d2(xyz,natoms,Hess,at,modh)
      !case (p_modh_lindh)
      !  if (pr) write (*,'(a)') "Using Lindh-Hessian (2007)"
      !  call mh_lindh(xyz,natoms,Hess,at,modh)
      !case (p_modh_swart)
      !  if (pr) write (*,'(a)') "Using Swart-Hessian"
      !  call mh_swart(xyz,natoms,Hess,at,modh)
    end select

!> add user-set constraint contributions to modelhessian
    call constrhess(natoms,at,xyz,calc,Hess)

    return
  end subroutine modhes

!========================================================================================!
!> subroutine ddvopt
!> generates a Lindh Model Hessian
!> Chem. Phys. Let. 241(1995) 423-428
!>
!> Input:
!>     Cart  - cartesian coordinates
!>    nAtoms - number of atoms
!>     iANr  - atom types as integers
!>    mhset  - model Hessian parameters
!>
!> Output:
!>     Hess  - the (packed) model Hessian
!>------------------------------------------------------
  subroutine ddvopt(Cart,nAtoms,Hess,iANr,mhset)
    Implicit Integer(i - n)
    Implicit Real(wp) (a - h,o - z)
    type(mhparam) :: mhset

    real(wp) :: s6
    real(wp) :: rcut

    Real(wp) :: Cart(3,nAtoms),rij(3),rjk(3),rkl(3), &
   &       Hess((3 * nAtoms) * (3 * nAtoms + 1) / 2),si(3),sj(3),sk(3), &
   &       sl(3),sm(3),x(2),y(2),z(2), &
   &       xyz(3,4),C(3,4),Dum(3,4,3,4)
    Integer iANr(nAtoms)

! include  "common/ddvdt.inc" (molpro 2002.6)
    Real(wp) :: rAV(3,3),aAV(3,3), &
   &       B_Str(6),A_Bend(2),A_Trsn(2),A_StrH(2), &
   &       rkr,rkf,A_Str,RF_Const, &
   &       wthr

    Data rAv/1.3500d+00,2.1000d+00,2.5300d+00, &
   &         2.1000d+00,2.8700d+00,3.4000d+00, &
   &         2.5300d+00,3.4000d+00,3.4000d+00/
    Data aAv/1.0000d+00,0.3949d+00,0.3949d+00, &
   &         0.3949d+00,0.2800d+00,0.2800d+00, &
   &         0.3949d+00,0.2800d+00,0.2800d+00/
!org  Data rkr,rkf,rkt/0.4500D+00,0.1500D+00,0.5000D-02/
    Data rkr,rkf,rkt/0.4000D+00,0.1300D+00,0.7500D-02/
    Data A_Str/1.734d0/
    Data B_Str/-.244d0,0.352d0,1.085d0,0.660d0,1.522d0,2.068d0/
    Data A_Bend/0.160d0,0.250d0/
    Data A_Trsn/0.0023d0,0.07d0/
    Data A_StrH/0.3601d0,1.944d0/
    Data RF_Const/1.0D-2/
    Data wthr/0.2/

!cc VDWx-Parameters (Grimme) used for vdw-correction of model hessian
    real(wp) :: alphavdw,damp,c6k,c6l,c66,vdw(3,3),dr(3)
    integer :: kxyz,lxyz
!cc End: VDWx ccccccccccccccccc

    !> BLAS
    external :: dcopy

    s6 = mhset%s6
    rcut = mhset%rcut

!
!------- Statement functions
!
!      ixyz(i,iAtom) = (iAtom-1)*3 + i
!      Jnd(i,j) = i*(i-1)/2 +j
!      Ind(i,iAtom,j,jAtom)=Jnd(Max(ixyz(i,iAtom),ixyz(j,jAtom)), &
!     &                         Min(ixyz(i,iAtom),ixyz(j,jAtom)))
!end

    Fact = One
!hjw threshold reduced
    rZero = 1.0d-10
    n3 = 3 * nAtoms
    Hess = 0.0d0

!
!     Hessian for tension
!
    Do kAtom = 1,nAtoms
      kr = iTabRow(iANr(kAtom))
!        If (kr.eq.0) Go To 5

      Do lAtom = 1,kAtom - 1
        lr = iTabRow(iANr(lAtom))
!           If (lr.eq.0) Go To 10
        xkl = Cart(1,kAtom) - Cart(1,lAtom)
        ykl = Cart(2,kAtom) - Cart(2,lAtom)
        zkl = Cart(3,kAtom) - Cart(3,lAtom)
        rkl2 = xkl**2 + ykl**2 + zkl**2
        r0 = rAv(kr,lr)
        alpha = aAv(kr,lr)

!cccccc VDWx ccccccccccccccccccccccccccccccccc
        c6k = c6(iANr(katom))
        c6l = c6(iANr(latom))
        c66 = sqrt(c6k * c6l)
        Rv = (vander(iANr(katom)) + vander(iANr(latom))) / bohr

        call getvdwxx(xkl,ykl,zkl,c66,s6,Rv,vdw(1,1))
        call getvdwxy(xkl,ykl,zkl,c66,s6,Rv,vdw(1,2))
        call getvdwxy(xkl,zkl,ykl,c66,s6,Rv,vdw(1,3))
        call getvdwxx(ykl,xkl,zkl,c66,s6,Rv,vdw(2,2))
        call getvdwxy(ykl,zkl,xkl,c66,s6,Rv,vdw(2,3))
        call getvdwxx(zkl,xkl,ykl,c66,s6,Rv,vdw(3,3))
!cccccc Ende VDWx ccccccccccccccccccccccccccccccc

        gamma = rkr * Exp(alpha * r0**2)
! not better: *sqrt(abs(wb(kAtom,lAtom)))
        gmm = gamma * Exp(-alpha * rkl2)
        Hxx = gmm * xkl * xkl / rkl2 - vdw(1,1)
        Hxy = gmm * xkl * ykl / rkl2 - vdw(1,2)
        Hxz = gmm * xkl * zkl / rkl2 - vdw(1,3)
        Hyy = gmm * ykl * ykl / rkl2 - vdw(2,2)
        Hyz = gmm * ykl * zkl / rkl2 - vdw(2,3)
        Hzz = gmm * zkl * zkl / rkl2 - vdw(3,3)

!
        Hess(Ind(1,kAtom,1,kAtom)) = Hess(Ind(1,kAtom,1,kAtom)) + Hxx
        Hess(Ind(2,kAtom,1,kAtom)) = Hess(Ind(2,kAtom,1,kAtom)) + Hxy
        Hess(Ind(2,kAtom,2,kAtom)) = Hess(Ind(2,kAtom,2,kAtom)) + Hyy
        Hess(Ind(3,kAtom,1,kAtom)) = Hess(Ind(3,kAtom,1,kAtom)) + Hxz
        Hess(Ind(3,kAtom,2,kAtom)) = Hess(Ind(3,kAtom,2,kAtom)) + Hyz
        Hess(Ind(3,kAtom,3,kAtom)) = Hess(Ind(3,kAtom,3,kAtom)) + Hzz
!
        Hess(Ind(1,kAtom,1,lAtom)) = Hess(Ind(1,kAtom,1,lAtom)) - Hxx
        Hess(Ind(1,kAtom,2,lAtom)) = Hess(Ind(1,kAtom,2,lAtom)) - Hxy
        Hess(Ind(1,kAtom,3,lAtom)) = Hess(Ind(1,kAtom,3,lAtom)) - Hxz
        Hess(Ind(2,kAtom,1,lAtom)) = Hess(Ind(2,kAtom,1,lAtom)) - Hxy
        Hess(Ind(2,kAtom,2,lAtom)) = Hess(Ind(2,kAtom,2,lAtom)) - Hyy
        Hess(Ind(2,kAtom,3,lAtom)) = Hess(Ind(2,kAtom,3,lAtom)) - Hyz
        Hess(Ind(3,kAtom,1,lAtom)) = Hess(Ind(3,kAtom,1,lAtom)) - Hxz
        Hess(Ind(3,kAtom,2,lAtom)) = Hess(Ind(3,kAtom,2,lAtom)) - Hyz
        Hess(Ind(3,kAtom,3,lAtom)) = Hess(Ind(3,kAtom,3,lAtom)) - Hzz
!
        Hess(Ind(1,lAtom,1,lAtom)) = Hess(Ind(1,lAtom,1,lAtom)) + Hxx
        Hess(Ind(2,lAtom,1,lAtom)) = Hess(Ind(2,lAtom,1,lAtom)) + Hxy
        Hess(Ind(2,lAtom,2,lAtom)) = Hess(Ind(2,lAtom,2,lAtom)) + Hyy
        Hess(Ind(3,lAtom,1,lAtom)) = Hess(Ind(3,lAtom,1,lAtom)) + Hxz
        Hess(Ind(3,lAtom,2,lAtom)) = Hess(Ind(3,lAtom,2,lAtom)) + Hyz
        Hess(Ind(3,lAtom,3,lAtom)) = Hess(Ind(3,lAtom,3,lAtom)) + Hzz
!
10      Continue
      End Do

5     Continue
    End Do

!
!     Hessian for bending
!
    Do mAtom = 1,nAtoms
      mr = iTabRow(iANr(mAtom))
!        If (mr.eq.0) Go To 20
      Do iAtom = 1,nAtoms
        If (iAtom .eq. mAtom) Go To 30
        ir = iTabRow(iANr(iAtom))
!          If (ir.eq.0) Go To 30
        if (rcutoff(cart,iatom,matom,rcut)) cycle
!          if(wb(iatom,matom).lt.wthr) cycle
        Do jAtom = 1,iAtom - 1
          If (jAtom .eq. mAtom) Go To 40
          jr = iTabRow(iANr(jAtom))
!           If (jr.eq.0) Go To 40
          if (rcutoff(cart,jatom,iatom,rcut)) cycle
          if (rcutoff(cart,jatom,matom,rcut)) cycle
!           if(wb(jatom,iatom).lt.wthr) cycle
!           if(wb(jatom,matom).lt.wthr) cycle

          xmi = (Cart(1,iAtom) - Cart(1,mAtom))
          ymi = (Cart(2,iAtom) - Cart(2,mAtom))
          zmi = (Cart(3,iAtom) - Cart(3,mAtom))
          rmi2 = xmi**2 + ymi**2 + zmi**2
          rmi = sqrt(rmi2)
          r0mi = rAv(mr,ir)
          ami = aAv(mr,ir)
!
          xmj = (Cart(1,jAtom) - Cart(1,mAtom))
          ymj = (Cart(2,jAtom) - Cart(2,mAtom))
          zmj = (Cart(3,jAtom) - Cart(3,mAtom))
          rmj2 = xmj**2 + ymj**2 + zmj**2
          rmj = sqrt(rmj2)
          r0mj = rAv(mr,jr)
          amj = aAv(mr,jr)
!
!---------- Test if zero angle
!
          Test = xmi * xmj + ymi * ymj + zmi * zmj
          Test = Test / (rmi * rmj)
          If (Test .eq. One) Go To 40
!
          xij = (Cart(1,jAtom) - Cart(1,iAtom))
          yij = (Cart(2,jAtom) - Cart(2,iAtom))
          zij = (Cart(3,jAtom) - Cart(3,iAtom))
          rij2 = xij**2 + yij**2 + zij**2
          rrij = sqrt(rij2)
!
          alpha = rkf * exp((ami * r0mi**2 + amj * r0mj**2))
!
          r = sqrt(rmj2 + rmi2)
          gij = alpha * exp(-(ami * rmi2 + amj * rmj2))
!           Write (*,*) ' gij=',gij
          rL2 = (ymi * zmj - zmi * ymj)**2 + (zmi * xmj - xmi * zmj)**2 + &
   &         (xmi * ymj - ymi * xmj)**2
!hjw modified
          if (rL2 .lt. 1.d-14) then
            rL = 0
          else
            rL = sqrt(rL2)
          end if
!
          if ((rmj .gt. rZero) .and. (rmi .gt. rZero) .and. &
   &                                (rrij .gt. rZero)) Then
            SinPhi = rL / (rmj * rmi)
            rmidotrmj = xmi * xmj + ymi * ymj + zmi * zmj
            CosPhi = rmidotrmj / (rmj * rmi)
!
!-------------None linear case
!
            If (SinPhi .gt. rZero) Then
!               Write (*,*) ' None linear case'
              si(1) = (xmi / rmi * cosphi - xmj / rmj) / (rmi * sinphi)
              si(2) = (ymi / rmi * cosphi - ymj / rmj) / (rmi * sinphi)
              si(3) = (zmi / rmi * cosphi - zmj / rmj) / (rmi * sinphi)
              sj(1) = (cosphi * xmj / rmj - xmi / rmi) / (rmj * sinphi)
              sj(2) = (cosphi * ymj / rmj - ymi / rmi) / (rmj * sinphi)
              sj(3) = (cosphi * zmj / rmj - zmi / rmi) / (rmj * sinphi)
              sm(1) = -si(1) - sj(1)
              sm(2) = -si(2) - sj(2)
              sm(3) = -si(3) - sj(3)
              Do icoor = 1,3
                Do jCoor = 1,3
                  If (mAtom .gt. iAtom) Then
                    Hess(Ind(icoor,mAtom,jcoor,iAtom)) = &
  &                        Hess(Ind(icoor,mAtom,jcoor,iAtom)) &
  &                        + gij * sm(icoor) * si(jcoor)
                  else
                    Hess(Ind(icoor,iAtom,jcoor,mAtom)) = &
   &                        Hess(Ind(icoor,iAtom,jcoor,mAtom)) &
   &                        + gij * si(icoor) * sm(jcoor)
                  End If
                  If (mAtom .gt. jAtom) Then
                    Hess(Ind(icoor,mAtom,jcoor,jAtom)) = &
 &                        Hess(Ind(icoor,mAtom,jcoor,jAtom)) &
 &                        + gij * sm(icoor) * sj(jcoor)
                  else
                    Hess(Ind(icoor,jAtom,jcoor,mAtom)) = &
   &                        Hess(Ind(icoor,jAtom,jcoor,mAtom)) &
   &                        + gij * sj(icoor) * sm(jcoor)
                  End If
                  If (iAtom .gt. jAtom) Then
                    Hess(Ind(icoor,iAtom,jcoor,jAtom)) = &
 &                        Hess(Ind(icoor,iAtom,jcoor,jAtom)) &
 &                        + gij * si(icoor) * sj(jcoor)
                  else
                    Hess(Ind(icoor,jAtom,jcoor,iAtom)) = &
 &                        Hess(Ind(icoor,jAtom,jcoor,iAtom)) &
 &                        + gij * sj(icoor) * si(jcoor)
                  End If
                End Do
              End Do
              Do icoor = 1,3
                Do jCoor = 1,icoor
                  Hess(Ind(icoor,iAtom,jcoor,iAtom)) = &
   &                        Hess(Ind(icoor,iAtom,jcoor,iAtom)) &
   &                        + gij * si(icoor) * si(jcoor)
                  Hess(Ind(icoor,mAtom,jcoor,mAtom)) = &
   &                        Hess(Ind(icoor,mAtom,jcoor,mAtom)) &
   &                        + gij * sm(icoor) * sm(jcoor)
                  Hess(Ind(icoor,jAtom,jcoor,jAtom)) = &
   &                        Hess(Ind(icoor,jAtom,jcoor,jAtom)) &
   &                        + gij * sj(icoor) * sj(jcoor)

!
                End Do
              End Do
            Else
!
!----------------Linear case
!
              if ((abs(ymi) .gt. rZero) .or. &
&                 (abs(xmi) .gt. rZero)) Then
                x(1) = -ymi
                y(1) = xmi
                z(1) = Zero
                x(2) = -xmi * zmi
                y(2) = -ymi * zmi
                z(2) = xmi * xmi + ymi * ymi
              Else
                x(1) = One
                y(1) = Zero
                z(1) = Zero
                x(2) = Zero
                y(2) = One
                z(2) = Zero
              End If
              Do i = 1,2
                r1 = sqrt(x(i)**2 + y(i)**2 + z(i)**2)
                cosThetax = x(i) / r1
                cosThetay = y(i) / r1
                cosThetaz = z(i) / r1
                si(1) = -cosThetax / rmi
                si(2) = -cosThetay / rmi
                si(3) = -cosThetaz / rmi
                sj(1) = -cosThetax / rmj
                sj(2) = -cosThetay / rmj
                sj(3) = -cosThetaz / rmj
                sm(1) = -(si(1) + sj(1))
                sm(2) = -(si(2) + sj(2))
                sm(3) = -(si(3) + sj(3))
!
                Do icoor = 1,3
                  Do jCoor = 1,3
                    If (mAtom .gt. iAtom) Then
                      Hess(Ind(icoor,mAtom,jcoor,iAtom)) = &
 &                        Hess(Ind(icoor,mAtom,jcoor,iAtom)) &
 &                         + gij * sm(icoor) * si(jcoor)
                    else
                      Hess(Ind(icoor,iAtom,jcoor,mAtom)) = &
&                        Hess(Ind(icoor,iAtom,jcoor,mAtom)) &
&                         + gij * si(icoor) * sm(jcoor)
                    End If
                    If (mAtom .gt. jAtom) Then
                      Hess(Ind(icoor,mAtom,jcoor,jAtom)) = &
 &                        Hess(Ind(icoor,mAtom,jcoor,jAtom)) &
 &                         + gij * sm(icoor) * sj(jcoor)
                    else
                      Hess(Ind(icoor,jAtom,jcoor,mAtom)) = &
 &                        Hess(Ind(icoor,jAtom,jcoor,mAtom)) &
 &                         + gij * sj(icoor) * sm(jcoor)
                    End If
                    If (iAtom .gt. jAtom) Then
                      Hess(Ind(icoor,iAtom,jcoor,jAtom)) = &
&                        Hess(Ind(icoor,iAtom,jcoor,jAtom)) &
&                         + gij * si(icoor) * sj(jcoor)
                    else
                      Hess(Ind(icoor,jAtom,jcoor,iAtom)) = &
&                        Hess(Ind(icoor,jAtom,jcoor,iAtom)) &
&                         + gij * sj(icoor) * si(jcoor)
                    End If
                  End Do
                End Do
                Do icoor = 1,3
                  Do jCoor = 1,icoor
                    Hess(Ind(icoor,iAtom,jcoor,iAtom)) = &
&                        Hess(Ind(icoor,iAtom,jcoor,iAtom)) &
&                         + gij * si(icoor) * si(jcoor)
                    Hess(Ind(icoor,mAtom,jcoor,mAtom)) = &
&                        Hess(Ind(icoor,mAtom,jcoor,mAtom)) &
&                         + gij * sm(icoor) * sm(jcoor)
                    Hess(Ind(icoor,jAtom,jcoor,jAtom)) = &
&                         Hess(Ind(icoor,jAtom,jcoor,jAtom)) &
&                         + gij * sj(icoor) * sj(jcoor)
                  End Do
                End Do
              End Do
            End If
          End If
!
40        Continue
        End Do
30      Continue
      End Do
20    Continue
    End Do
!
!     Hessian for torsion
!
    Do jAtom = 1,nAtoms
      jr = iTabRow(iANr(jAtom))
!       If (jr.eq.0) Go To 444
!
      Call DCopy(3,Cart(1,jAtom),1,xyz(1,2),1)
!
      Do kAtom = 1,nAtoms
        If (kAtom .eq. jAtom) Go To 111
        kr = iTabRow(iANr(kAtom))
!          If (kr.eq.0) Go To 111

        if (rcutoff(cart,katom,jatom,rcut)) cycle
!          if(wb(katom,jatom).lt.wthr) cycle
!
        Call DCopy(3,Cart(1,kAtom),1,xyz(1,3),1)
!
        Do iAtom = 1,nAtoms
          ij_ = nAtoms * (jAtom - 1) + iAtom
          If (iAtom .eq. jAtom) Go To 333
          If (iAtom .eq. kAtom) Go To 333
          ir = iTabRow(iANr(iAtom))
!             If (ir.eq.0) Go To 333
!
          if (rcutoff(cart,iatom,katom,rcut)) cycle
          if (rcutoff(cart,iatom,jatom,rcut)) cycle
!             if(wb(iatom,katom).lt.wthr) cycle
!             if(wb(iatom,jatom).lt.wthr) cycle

          Call DCopy(3,Cart(1,iAtom),1,xyz(1,1),1)
!
          Do lAtom = 1,nAtoms
            lk_ = nAtoms * (kAtom - 1) + lAtom
            If (ij_ .le. lk_) Go To 222
            If (lAtom .eq. iAtom) Go To 222
            If (lAtom .eq. jAtom) Go To 222
            If (lAtom .eq. kAtom) Go To 222
            lr = iTabRow(iANr(lAtom))
!                If (lr.eq.0) Go To 222
!
            if (rcutoff(cart,latom,iatom,rcut)) cycle
            if (rcutoff(cart,latom,katom,rcut)) cycle
            if (rcutoff(cart,latom,jatom,rcut)) cycle
!                if(wb(latom,iatom).lt.wthr) cycle
!                if(wb(latom,katom).lt.wthr) cycle
!                if(wb(latom,jatom).lt.wthr) cycle

            Call DCopy(3,Cart(1,lAtom),1,xyz(1,4),1)
!
            rij(1) = Cart(1,iAtom) - Cart(1,jAtom)
            rij(2) = Cart(2,iAtom) - Cart(2,jAtom)
            rij(3) = Cart(3,iAtom) - Cart(3,jAtom)
            rij0 = rAv(ir,jr)**2
            aij = aAv(ir,jr)
!
            rjk(1) = Cart(1,jAtom) - Cart(1,kAtom)
            rjk(2) = Cart(2,jAtom) - Cart(2,kAtom)
            rjk(3) = Cart(3,jAtom) - Cart(3,kAtom)
            rjk0 = rAv(jr,kr)**2
            ajk = aAv(jr,kr)
!
            rkl(1) = Cart(1,kAtom) - Cart(1,lAtom)
            rkl(2) = Cart(2,kAtom) - Cart(2,lAtom)
            rkl(3) = Cart(3,kAtom) - Cart(3,lAtom)
            rkl0 = rAv(kr,lr)**2
            akl = aAv(kr,lr)
!
            rij2 = rij(1)**2 + rij(2)**2 + rij(3)**2
            rjk2 = rjk(1)**2 + rjk(2)**2 + rjk(3)**2
            rkl2 = rkl(1)**2 + rkl(2)**2 + rkl(3)**2
!              Allow only angles in the range of 35-145
            A35 = (35.0D0 / 180.D0) * Pi
            CosFi_Max = Cos(A35)
            CosFi2 = (rij(1) * rjk(1) + rij(2) * rjk(2) + rij(3) * rjk(3)) &
  &               / Sqrt(rij2 * rjk2)
            If (Abs(CosFi2) .gt. CosFi_Max) Go To 222
            CosFi3 = (rkl(1) * rjk(1) + rkl(2) * rjk(2) + rkl(3) * rjk(3)) &
  &               / Sqrt(rkl2 * rjk2)
            If (Abs(CosFi3) .gt. CosFi_Max) Go To 222

            beta = rkt * &
  &                       exp((aij * rij0 + ajk * rjk0 + akl * rkl0))
            tij = beta * exp(-(aij * rij2 + ajk * rjk2 + akl * rkl2))

            Call Trsn(xyz,4,Tau,C,.False.,.False.,'        ', &
  &                  Dum,.False.)
            Call DCopy(3,C(1,1),1,si,1)
            Call DCopy(3,C(1,2),1,sj,1)
            Call DCopy(3,C(1,3),1,sk,1)
            Call DCopy(3,C(1,4),1,sl,1)
!
!-------------Off diagonal block
!
            Do icoor = 1,3
              Do jCoor = 1,3
                Hess(Ind(icoor,iAtom,jcoor,jAtom)) = &
    &           Hess(Ind(icoor,iAtom,jcoor,jAtom)) &
    &            + tij * si(icoor) * sj(jcoor)
                Hess(Ind(icoor,iAtom,jcoor,kAtom)) = &
    &           Hess(Ind(icoor,iAtom,jcoor,kAtom)) &
    &            + tij * si(icoor) * sk(jcoor)
                Hess(Ind(icoor,iAtom,jcoor,lAtom)) = &
    &           Hess(Ind(icoor,iAtom,jcoor,lAtom)) &
    &            + tij * si(icoor) * sl(jcoor)
                Hess(Ind(icoor,jAtom,jcoor,kAtom)) = &
    &           Hess(Ind(icoor,jAtom,jcoor,kAtom)) &
    &            + tij * sj(icoor) * sk(jcoor)
                Hess(Ind(icoor,jAtom,jcoor,lAtom)) = &
    &           Hess(Ind(icoor,jAtom,jcoor,lAtom)) &
    &            + tij * sj(icoor) * sl(jcoor)
                Hess(Ind(icoor,kAtom,jcoor,lAtom)) = &
    &           Hess(Ind(icoor,kAtom,jcoor,lAtom)) &
    &            + tij * sk(icoor) * sl(jcoor)

              End Do
            End Do
!
!-------------Diagonal block
!
            Do icoor = 1,3
              Do jCoor = 1,icoor
                Hess(Ind(icoor,iAtom,jcoor,iAtom)) = &
    &           Hess(Ind(icoor,iAtom,jcoor,iAtom)) &
    &            + tij * si(icoor) * si(jcoor)
                Hess(Ind(icoor,jAtom,jcoor,jAtom)) = &
    &           Hess(Ind(icoor,jAtom,jcoor,jAtom)) &
    &            + tij * sj(icoor) * sj(jcoor)
                Hess(Ind(icoor,kAtom,jcoor,kAtom)) = &
    &           Hess(Ind(icoor,kAtom,jcoor,kAtom)) &
    &            + tij * sk(icoor) * sk(jcoor)
                Hess(Ind(icoor,lAtom,jcoor,lAtom)) = &
    &           Hess(Ind(icoor,lAtom,jcoor,lAtom)) &
    &            + tij * sl(icoor) * sl(jcoor)

!
              End Do
            End Do
222         Continue
          End Do        ! lAtom
333       Continue
        End Do          ! iAtom
111     Continue
      End Do             ! kAtom
444   Continue
    End Do               ! jAtom
    Return

  contains
    function ixyz(i,iatom)
      integer :: ixyz
      integer,intent(in) :: i,iatom
      ixyz = (iatom - 1) * 3 + i
    end function ixyz
    function jnd(i,j)
      integer :: jnd
      integer,intent(in) :: i,j
      jnd = i * (i - 1) / 2 + j
    end function jnd
    function ind(i,iatom,j,jatom)
      integer :: ind
      integer,intent(in) :: i,iatom,j,jatom
      ind = jnd(max(ixyz(i,iatom),ixyz(j,jatom)),min(ixyz(i,iatom),ixyz(j,jatom)))
    end function ind
  end subroutine ddvopt

!========================================================================================!
  logical function rcutoff(cart,katom,latom,rcut)
    implicit none
    Real(wp) :: Cart(3,*),xkl,ykl,zkl,rkl2
    real(wp) :: rcut
    integer :: katom,latom
    rcutoff = .false.
    xkl = Cart(1,kAtom) - Cart(1,lAtom)
    ykl = Cart(2,kAtom) - Cart(2,lAtom)
    zkl = Cart(3,kAtom) - Cart(3,lAtom)
    rkl2 = xkl**2 + ykl**2 + zkl**2
    if (rkl2 .gt. rcut) rcutoff = .true.
  end function

  function itabrow(i)
    integer :: itabrow
    integer,intent(in) :: i

    itabrow = 0
    if (i .gt. 0 .and. i .le. 2) then
      itabrow = 1
    else if (i .gt. 2 .and. i .le. 10) then
      itabrow = 2
    else if (i .gt. 10 .and. i .le. 18) then
      itabrow = 3
    else if (i .gt. 18 .and. i .le. 36) then
      itabrow = 3
    else if (i .gt. 36 .and. i .le. 54) then
      itabrow = 3
    else if (i .gt. 54 .and. i .le. 86) then
      itabrow = 3
    else if (i .gt. 86) then
      itabrow = 3
    end if

    return
  end function itabrow

  subroutine getvdwxy(rx,ry,rz,c66,s6,r0,vdw)
    !cc Ableitung nach rx und ry
    Implicit Real * 8(a - h,o - z)
    integer k,l
    !    write(*,*) 's6:', s6
    avdw = 20.0
    t1 = s6 * C66
    t2 = rx**2
    t3 = ry**2
    t4 = rz**2
    t5 = t2 + t3 + t4
    t6 = t5**2
    t7 = t6**2
    t11 = sqrt(t5)
    t12 = 0.1D1 / r0
    t16 = exp(-avdw * (t11 * t12 - 0.1D1))
    t17 = 0.1D1 + t16
    t25 = t17**2
    t26 = 0.1D1 / t25
    t35 = 0.1D1 / t7
    t40 = avdw**2
    t41 = r0**2
    t43 = t40 / t41
    t44 = t16**2
    t56 = -0.48D2 * t1 / t7 / t5 / t17 * rx * ry + 0.13D2 * t1 / t11 /&
       & t7 * t26 * rx * avdw * t12 * ry * t16 - 0.2D1 * t1 * t35 / t25 /&
       &t17 * t43 * rx * t44 * ry + t1 * t35 * t26 * t43 * rx * ry * t16
    vdw = t56
    return
  end subroutine getvdwxy

  subroutine getvdwxx(rx,ry,rz,c66,s6,r0,vdw)
    !cc Ableitung nach rx und rx
    Implicit Real * 8(a - h,o - z)
    avdw = 20.0
    !      write(*,*) 's6:', s6
    t1 = s6 * C66
    t2 = rx**2
    t3 = ry**2
    t4 = rz**2
    t5 = t2 + t3 + t4
    t6 = t5**2
    t7 = t6**2
    t10 = sqrt(t5)
    t11 = 0.1D1 / r0
    t15 = exp(-avdw * (t10 * t11 - 0.1D1))
    t16 = 0.1D1 + t15
    t17 = 0.1D1 / t16
    t24 = t16**2
    t25 = 0.1D1 / t24
    t29 = t11 * t15
    t33 = 0.1D1 / t7
    t41 = avdw**2
    t42 = r0**2
    t44 = t41 / t42
    t45 = t15**2
    t62 = -0.48D2 * t1 / t7 / t5 * t17 * t2 + 0.13D2 * t1 / t10 / t7 *&
       & t25 * t2 * avdw * t29 + 0.6D1 * t1 * t33 * t17 - 0.2D1 * t1 * t33&
       & / t24 / t16 * t44 * t2 * t45 - t1 / t10 / t6 / t5 * t25 * avdw *&
       &t29 + t1 * t33 * t25 * t44 * t2 * t15
    vdw = t62
  end subroutine getvdwxx

  pure subroutine trsn2(xyz,tau,bt)
    implicit none
    real(wp),intent(out) :: bt(3,4)
    real(wp),intent(out) :: tau
    real(wp),intent(in)  :: xyz(3,4)
    real(wp) :: rij(3),rij1,brij(3,2)
    real(wp) :: rjk(3),rjk1,brjk(3,2)
    real(wp) :: rkl(3),rkl1,brkl(3,2)
    real(wp) :: bf2(3,3),fi2,sinfi2,cosfi2
    real(wp) :: bf3(3,3),fi3,sinfi3,cosfi3
    real(wp) :: costau,sintau
    integer  :: ix,iy,iz
    call strtch2(xyz(1,1),rij1,brij)
    call strtch2(xyz(1,2),rjk1,brjk)
    call strtch2(xyz(1,3),rkl1,brkl)
    call bend2(xyz(1,1),fi2,bf2)
    sinfi2 = sin(fi2)
    cosfi2 = cos(fi2)
    call bend2(xyz(1,2),fi3,bf3)
    sinfi3 = sin(fi3)
    cosfi3 = cos(fi3)
    costau = ((brij(2,1) * brjk(3,2) - brij(3,1) * brjk(2,2)) * &
              (brjk(2,1) * brkl(3,2) - brjk(3,1) * brkl(2,2)) + &
              (brij(3,1) * brjk(1,2) - brij(1,1) * brjk(3,2)) * &
              (brjk(3,1) * brkl(1,2) - brjk(1,1) * brkl(3,2)) + &
              (brij(1,1) * brjk(2,2) - brij(2,1) * brjk(1,2)) * &
              (brjk(1,1) * brkl(2,2) - brjk(2,1) * brkl(1,2))) &
             / (sinfi2 * sinfi3)
    sintau = (brij(1,2) * (brjk(2,1) * brkl(3,2) - brjk(3,1) * brkl(2,2)) &
              + brij(2,2) * (brjk(3,1) * brkl(1,2) - brjk(1,1) * brkl(3,2)) &
              + brij(3,2) * (brjk(1,1) * brkl(2,2) - brjk(2,1) * brkl(1,2))) &
             / (sinfi2 * sinfi3)
    tau = atan2(sintau,costau)
    if (abs(tau) .eq. pi) tau = pi
    do ix = 1,3
      iy = ix + 1
      if (iy .gt. 3) iy = iy - 3
      iz = iy + 1
      if (iz .gt. 3) iz = iz - 3
      bt(ix,1) = (brij(iy,2) * brjk(iz,2) - brij(iz,2) * brjk(iy,2)) &
         &           / (rij1 * sinfi2**2)
      bt(ix,4) = (brkl(iy,1) * brjk(iz,1) - brkl(iz,1) * brjk(iy,1)) &
         &           / (rkl1 * sinfi3**2)
      bt(ix,2) = -((rjk1 - rij1 * cosfi2) * bt(ix,1) &
         &             + rkl1 * cosfi3 * bt(ix,4)) / rjk1
      bt(ix,3) = -(bt(ix,1) + bt(ix,2) + bt(ix,4))
    end do
  contains
    pure subroutine strtch2(xyz,avst,b)
      implicit none
      real(wp),intent(out) :: b(3,2)
      real(wp),intent(in)  :: xyz(3,2)
      real(wp) :: r(3)
      real(wp) :: rr
      real(wp),intent(out) :: avst
      r = xyz(:,2) - xyz(:,1)
      rr = norm2(r)
      avst = rr
      b(:,1) = -r / rr
      b(:,2) = -b(:,1)
    end subroutine strtch2
    pure subroutine bend2(xyz,fir,bf)
      implicit none
      real(wp),intent(out) :: bf(3,3)
      real(wp),intent(in)  :: xyz(3,3)
      real(wp) :: brij(3,2)
      real(wp) :: brjk(3,2)
      real(wp) :: co,crap
      real(wp),intent(out) :: fir
      real(wp) :: si
      real(wp) :: rij1,rjk1
      integer  :: i
      call strtch2(xyz(1,1),rij1,brij)
      call strtch2(xyz(1,2),rjk1,brjk)
      co = 0.0_wp
      crap = 0.0_wp
      do i = 1,3
        co = co + brij(i,1) * brjk(i,2)
        crap = crap + (brjk(i,2) + brij(i,1))**2
      end do
      if (sqrt(crap) .lt. 1.0d-6) then
        fir = pi - asin(sqrt(crap))
        si = sqrt(crap)
      else
        fir = acos(co)
        si = sqrt(1.0_wp - co**2)
      end if
      if (abs(fir - pi) .lt. 1.0d-13) then
        fir = pi
        return
      end if
      do i = 1,3
        bf(i,1) = (co * brij(i,1) - brjk(i,2)) / (si * rij1)
        bf(i,3) = (co * brjk(i,2) - brij(i,1)) / (si * rjk1)
        bf(i,2) = -(bf(i,1) + bf(i,3))
      end do
    end subroutine bend2
  end subroutine trsn2

  Subroutine Trsn(xyz,nCent,Tau,Bt,lWrite,lWarn,Label,dBt,ldB)
    !***********************************************************************
    !                                                                      *
    ! Reference: Molecular Vibrations, E. Bright Wilson, Jr, J. C. Decicius*
    !             nd Paul C. Cross, Sec. 4-1, Eq. 20-24                    *
    !                                                                      *
    ! R.Lindh May-June '96                                                 *
    !***********************************************************************
    Implicit Real(wp) (a - h,o - z)

    integer :: nCent,mCent,i,j,ix,iy,iz,jx,jy,jz
    Real(wp) Bt(3,nCent),xyz(3,nCent),Rij(3),Eij(3),Rjk(3),Ejk(3),&
       &       Rkl(3),Ekl(3),Rijk(3),Eijk(3),dBt(3,nCent,3,nCent),&
       &       BRij(3,2),dBRij(3,2,3,2),BRjk(3,2),dBRjk(3,2,3,2),&
       &       BRkl(3,2),dBRkl(3,2,3,2),Bf2(3,3),dum(3,4,3,4),&
       &       Bf3(3,3)
    Logical :: lWrite,lWarn,ldB
    Character(len=8) :: Label
    !
    !     Call qEnter('Trsn')
    mCent = 2
    Call Strtch(xyz(1,1),mCent,Rij1,BRij,.False.,Label,dBRij,ldB)
    Call Strtch(xyz(1,2),mCent,Rjk1,BRjk,.False.,Label,dBRjk,ldB)
    Call Strtch(xyz(1,3),mCent,Rkl1,BRkl,.False.,Label,dBRkl,ldB)
    mCent = 3
    Call Bend(xyz(1,1),mCent,Fi2,Bf2,.False.,.False.,Label,Dum,&
       &          .False.)
    SinFi2 = Sin(Fi2)
    CosFi2 = Cos(Fi2)
    Call Bend(xyz(1,2),mCent,Fi3,Bf3,.False.,.False.,Label,Dum,&
       &          .False.)
    SinFi3 = Sin(Fi3)
    CosFi3 = Cos(Fi3)
    !
    !     Get the angle between the two planes, i.e. the
    !     angle between the normal vectors.
    !
    !     r123 * r234 = CosTau
    !
    CosTau = ((BRij(2,1) * BRjk(3,2) - BRij(3,1) * BRjk(2,2)) *&
       &           (BRjk(2,1) * BRkl(3,2) - BRjk(3,1) * BRkl(2,2)) +&
       &           (BRij(3,1) * BRjk(1,2) - BRij(1,1) * BRjk(3,2)) *&
       &           (BRjk(3,1) * BRkl(1,2) - BRjk(1,1) * BRkl(3,2)) +&
       &           (BRij(1,1) * BRjk(2,2) - BRij(2,1) * BRjk(1,2)) *&
       &           (BRjk(1,1) * BRkl(2,2) - BRjk(2,1) * BRkl(1,2)))&
       &         / (SinFi2 * SinFi3)
    !
    !     For the vector product of the two vectors. This
    !     will give a vector parallell to e23. The direction
    !     relative to e23 defines the sign.
    !
    !     e123 X e234 = SinTau * e23
    !
    SinTau = (BRij(1,2) * (BRjk(2,1) * BRkl(3,2) - BRjk(3,1) * BRkl(2,2))&
       &         + BRij(2,2) * (BRjk(3,1) * BRkl(1,2) - BRjk(1,1) * BRkl(3,2))&
       &         + BRij(3,2) * (BRjk(1,1) * BRkl(2,2) - BRjk(2,1) * BRkl(1,2)))&
       &         / (SinFi2 * SinFi3)
    !
    !     (-Pi < Tau <= Pi)
    !
    Tau = ATan2(SinTau,CosTau)
    If (Abs(Tau) .eq. Pi) Tau = Pi
    !
    dTau = 180.0D+00 * Tau / Pi
    dFi2 = 180.0D+00 * Fi2 / Pi
    dFi3 = 180.0D+00 * Fi3 / Pi
    If (lWarn) Then
      If (dTau .gt. 177.5 .or. dTau .lt. -177.5) Then
        Write (*,*) ' Warning: dihedral angle close to'&
           &         //' end of range'
      End If
      If (dFi2 .gt. 177.5 .or. dFi2 .lt. 2.5) Then
        Write (*,*) ' Warning: bond angle close to'&
           &         //' end of range'
      End If
      If (dFi3 .gt. 177.5 .or. dFi3 .lt. 2.5) Then
        Write (*,*) ' Warning: bond angle close to'&
           &         //' end of range'
      End If
    End If
    If (LWRITE) Write (*,1) Label,dTau,Tau
1   FORMAT(1X,A,' : Dihedral Angle=',F10.4,&
       & '/degree,',F10.4,'/rad')
    !
    !---- Compute the WDC matrix.
    !
    Do ix = 1,3
      iy = ix + 1
      If (iy .gt. 3) iy = iy - 3
      iz = iy + 1
      If (iz .gt. 3) iz = iz - 3
      Bt(ix,1) = (BRij(iy,2) * BRjk(iz,2) - BRij(iz,2) * BRjk(iy,2))&
         &           / (Rij1 * SinFi2**2)
      Bt(ix,4) = (BRkl(iy,1) * BRjk(iz,1) - BRkl(iz,1) * BRjk(iy,1))&
         &           / (Rkl1 * SinFi3**2)
      Bt(ix,2) = -((Rjk1 - Rij1 * CosFi2) * Bt(ix,1)&
         &             + Rkl1 * CosFi3 * Bt(ix,4)) / Rjk1
      Bt(ix,3) = -(Bt(ix,1) + Bt(ix,2) + Bt(ix,4))
    End Do
    !
    If (ldB) Then
      !
      !------- Compute the derivative of the WDC matrix.
      !
      Do ix = 1,3
        iy = ix + 1
        If (iy .gt. 3) iy = iy - 3
        iz = iy + 1
        If (iz .gt. 3) iz = iz - 3
        Do jx = 1,ix
          jy = jx + 1
          If (jy .gt. 3) jy = jy - 3
          jz = jy + 1
          If (jz .gt. 3) jz = jz - 3
          !
          dBt(ix,1,jx,1) = (dBRij(ix,1,jy,2) * BRjk(jz,2)&
             &                       - dBRij(ix,1,jz,2) * BRjk(jy,2)&
             &                       - Bt(jx,1) * (BRij(ix,1) * SinFi2**2&
             &                       + Rij1 * Two * SinFi2 * CosFi2 * Bf2(ix,1)))&
             &                       / (Rij1 * SinFi2**2)
          dBt(ix,1,jx,2) = -((-BRij(ix,1) * CosFi2&
             &                         + Rij1 * SinFi2 * Bf2(ix,1)) * Bt(jx,1)&
             &                         + (Rjk1 - Rij1 * CosFi2) * dBt(ix,1,jx,1))&
             &                       / Rjk1
          dBt(jx,2,ix,1) = dBt(ix,1,jx,2)
          dBt(ix,1,jx,4) = Zero
          dBt(jx,4,ix,1) = dBt(ix,1,jx,4)
          dBt(ix,1,jx,3) = -(dBt(ix,1,jx,1) + dBt(ix,1,jx,2))
          dBt(jx,3,ix,1) = dBt(ix,1,jx,3)
          dBt(ix,4,jx,4) = (dBRkl(ix,2,jy,1) * BRjk(jz,1)&
             &                       - dBRkl(ix,2,jz,1) * BRjk(jy,1)&
             &                       - Bt(jx,4) * (BRkl(ix,2) * SinFi3**2&
             &                       + Rkl1 * Two * SinFi3 * CosFi3 * Bf3(ix,3)))&
             &                       / (Rkl1 * SinFi3**2)
          dBt(ix,4,jx,3) = -((-BRkl(ix,2) * CosFi3&
             &                         + Rkl1 * SinFi3 * Bf3(ix,3)) * Bt(jx,4)&
             &                         + (Rjk1 - Rkl1 * CosFi3) * dBt(ix,4,jx,4))&
             &                       / Rjk1
          dBt(jx,3,ix,4) = dBt(ix,4,jx,3)
          dBt(ix,4,jx,2) = -(dBt(ix,4,jx,4) + dBt(ix,4,jx,3))
          dBt(jx,2,ix,4) = dBt(ix,4,jx,2)
          If (ix .ne. jx) Then
            dBt(jx,1,ix,1) = dBt(ix,1,jx,1)
            dBt(ix,4,jx,1) = Zero
            dBt(jx,4,ix,4) = dBt(ix,4,jx,4)
            dBt(jx,1,ix,4) = dBt(ix,4,jx,1)
            dBt(jx,1,ix,2) = -((-BRij(jx,1) * CosFi2&
               &                            + Rij1 * SinFi2 * Bf2(jx,1)) * Bt(ix,1)&
               &                            + (Rjk1 - Rij1 * CosFi2) * dBt(jx,1,ix,1))&
               &                          / Rjk1
            dBt(ix,2,jx,1) = dBt(jx,1,ix,2)
            dBt(ix,3,jx,1) = -(dBt(ix,1,jx,1) + dBt(ix,2,jx,1)&
               &                          + dBt(ix,4,jx,1))
            dBt(jx,1,ix,3) = dBt(ix,3,jx,1)
            dBt(jx,4,ix,3) = -((-BRkl(jx,2) * CosFi3&
               &                            + Rkl1 * SinFi3 * Bf3(jx,3)) * Bt(ix,4)&
               &                            + (Rjk1 - Rkl1 * CosFi3) * dBt(jx,4,ix,4))&
               &                          / Rjk1
            dBt(ix,3,jx,4) = dBt(jx,4,ix,3)
            dBt(ix,2,jx,4) = -(dBt(ix,4,jx,4) + dBt(ix,3,jx,4))
            dBt(jx,4,ix,2) = dBt(ix,2,jx,4)
          End If
          dBt(ix,2,jx,3) = -((BRjk(ix,1)&
             &                           + Rkl1 * SinFi3 * Bf3(ix,1)) * Bt(jx,4)&
             &                         + (Rjk1 - Rkl1 * CosFi3) * dBt(ix,2,jx,4)&
             &                         + (BRij(ix,2) * CosFi2&
             &                           - Rij1 * SinFi2 * Bf2(ix,2)) * Bt(jx,1)&
             &                         + Rij1 * CosFi2 * dBt(ix,2,jx,1)&
             &                         + Bt(jx,3) * BRjk(ix,1)) / Rjk1
          dBt(jx,3,ix,2) = dBt(ix,2,jx,3)
          dBt(ix,2,jx,2) = -(dBt(ix,2,jx,1) + dBt(ix,2,jx,4)&
             &                         + dBt(ix,2,jx,3))
          dBt(ix,3,jx,3) = -(dBt(ix,2,jx,3) + dBt(ix,1,jx,3)&
             &                         + dBt(ix,4,jx,3))
          If (ix .ne. jx) Then
            dBt(ix,3,jx,2) = -(dBt(ix,2,jx,2) + dBt(ix,1,jx,2)&
               &                            + dBt(ix,4,jx,2))
            dBt(jx,2,ix,3) = dBt(ix,3,jx,2)
            dBt(jx,2,ix,2) = dBt(ix,2,jx,2)
            dBt(jx,3,ix,3) = dBt(ix,3,jx,3)
          End If
          !
        End Do
      End Do
      !
    End If
    !     Call qExit('Trsn')
    Return
  contains
    Subroutine Strtch(xyz,nCent,Avst,B,lWrite,Label,dB,ldB)
      Implicit Real(wp) (a - h,o - z)
      !      include "common/real.inc"
      !comdeck real.inc $Revision: 2002.3 $
      Real(wp) :: Zero,One,Two,Three,Four,Five,Six,Seven,&
         &       Eight,RNine,Ten,Half,Pi,SqrtP2,TwoP34,&
         &       TwoP54,One2C2
      Parameter(Zero=0.0D0,One=1.0D0,Two=2.0D0,Three=3.0D0,&
         &          Four=4.0D0,Five=5.0D0,Six=6.0D0,Seven=7.0D0,&
         &          Eight=8.0D0,rNine=9.0D0,Ten=1.0D1,Half=0.5D0,&
         &          Pi=3.141592653589793D0,&
         &          SqrtP2=0.8862269254527579D0,&
         &          TwoP34=0.2519794355383808D0,&
         &          TwoP54=5.914967172795612D0,&
         &          One2C2=0.2662567690426443D-04)

      integer :: nCent
      Real(wp) :: B(3,nCent),xyz(3,nCent),dB(3,nCent,3,nCent),R(3)
      Logical :: lWrite,ldB
      Character(len=8) :: Label
      !      include "common/angstr.inc"
      !comdeck angstr.inc $Revision: 2002.3 $
      !
      !     Conversion factor angstrom to bohr from the IUPAC
      !     publication
      !     .529177249(24) angstrom / bohr
      !     "Quantities, Units and Symbols in Physical Chemistry"
      !     I. Mills, T. Cvitas, K. Homann, N. Kallay and
      !     K. Kuchitsu, Blackwell Scientific Publications,
      !     Oxford, 1988.
      !
      Data Angstr/0.529177249D+00/
      !
      R(1) = xyz(1,2) - xyz(1,1)
      R(2) = xyz(2,2) - xyz(2,1)
      R(3) = xyz(3,2) - xyz(3,1)
      R2 = R(1)**2 + R(2)**2 + R(3)**2
      RR = Sqrt(R2)
      Avst = RR
      !
      aRR = RR * Angstr
      If (lWrite) Write (*,'(1X,A,A,2(F10.6,A))') Label,&
         &      ' : Bond Length=',aRR,' / Angstrom',RR,' / bohr'
      !
      !---- Compute the WDC B-matrix.
      !
      B(1,1) = -R(1) / RR
      B(2,1) = -R(2) / RR
      B(3,1) = -R(3) / RR
      !.... Utilize translational invariance.
      B(1,2) = -B(1,1)
      B(2,2) = -B(2,1)
      B(3,2) = -B(3,1)
      !
      !---- Compute the cartesian derivative of the B-matrix.
      !
      If (ldB) Then
        !
        Do i = 1,3
          Do j = 1,i
            If (i .eq. j) Then
              dB(i,1,j,1) = (One - B(j,1) * B(i,1)) / RR
            Else
              dB(i,1,j,1) = (-B(j,1) * B(i,1)) / RR
            End If
            dB(j,1,i,1) = dB(i,1,j,1)
            !
            dB(i,2,j,1) = -dB(i,1,j,1)
            dB(j,1,i,2) = dB(i,2,j,1)
            !
            dB(i,1,j,2) = -dB(i,1,j,1)
            dB(j,2,i,1) = dB(i,1,j,2)
            !
            dB(i,2,j,2) = -dB(i,2,j,1)
            dB(j,2,i,2) = dB(i,2,j,2)
          End Do
        End Do
        !
      End If
      !     Call qExit('Strtch')
      !     Call GetMem('Exit Strtch','Chec','Real',ipMass,2*msAtom)
      Return
    End subroutine strtch
    Subroutine Bend(xyz,nCent,Fir,Bf,lWrite,lWarn,Label,dBf,ldB)
      Implicit Real(wp) (a - h,o - z)

      integer :: nCent
      !Real(wp) ::   Bf(3,nCent),xyz(3,nCent),dBf(3,nCent,3,nCent),&
      Real(wp) ::   Bf(3,3),xyz(3,nCent),dBf(3,nCent,3,nCent),&
         &        BRij(3,2),dBRij(3,2,3,2),&
         &        BRjk(3,2),dBRjk(3,2,3,2)
      Logical lWrite,ldB,lWarn
      Character(len=8) :: Label
      !
      !     Call QEnter('Bend')
      !
      mCent = 2
      Call Strtch(xyz(1,1),mCent,Rij1,BRij,.False.,Label,dBRij,ldB)
      Call Strtch(xyz(1,2),mCent,Rjk1,BRjk,.False.,Label,dBRjk,ldB)
      Co = Zero
      Crap = Zero
      Do i = 1,3
        Co = Co + BRij(i,1) * BRjk(i,2)
        Crap = Crap + (BRjk(i,2) + BRij(i,1))**2
      End Do
      !
      !.... Special care for cases close to linearity
      !
      If (Sqrt(Crap) .lt. 1.0D-6) Then
        Fir = Pi - ArSin(Sqrt(Crap))
        Si = Sqrt(Crap)
      Else
        Fir = ArCos(Co)
        Si = Sqrt(One - Co**2)
      End If
      !
      If (Abs(Fir - Pi) .lt. 1.0d-13) Then
        Fir = Pi
        Return
      End If
      dFir = 180.0D0 * Fir / Pi
      If ((Abs(dFir) .gt. 177.5 .or. Abs(dFir) .lt. 2.5) .and. lWarn)&
         &   Write (*,*) ' Valence angle close to end in '//&
         &               'range of definition'
      If (lWrite) Write (*,'(1X,A,A,F10.4,A,F10.6,A)') Label,&
         &            ' : Angle=',dFir,'/degree, ',Fir,'/rad'
      !
      !---- Compute the WDC B-matrix
      !
      !     Bf=-11.1111
      Do i = 1,3
        Bf(i,1) = (Co * BRij(i,1) - BRjk(i,2)) / (Si * Rij1)
        Bf(i,3) = (Co * BRjk(i,2) - BRij(i,1)) / (Si * Rjk1)
        !....... Utilize translational invariance.
        Bf(i,2) = -(Bf(i,1) + Bf(i,3))
      End Do
      !     Call RecPrt('Bf',' ',Bf,9,1)
      !
      !---- Compute the cartesian derivative of the B-Matrix.
      !
      If (ldB) Then
        !
        !        dBf=-11.11111
        Do i = 1,3
          Do j = 1,i
            dBf(i,1,j,1) = (-Si * Bf(i,1) * BRij(j,1)&
               &                        + Co * dBRij(i,1,j,1)&
               &                        - Bf(j,1) * (Co * Bf(i,1) * Rij1&
               &                        + Si * BRij(i,1))) / (Si * Rij1)
            dBf(i,1,j,3) = (-Si * Bf(i,1) * BRjk(j,2)&
               &                       + dBRij(i,1,j,2)&
               &                       - Bf(j,3) * Co * Bf(i,1) * Rjk1)&
               &                       / (Si * Rjk1)
            !              Write (*,*) '13',dBf(i,1,j,3), i, j
            dBf(i,3,j,1) = (-Si * Bf(i,3) * BRij(j,1)&
               &                       + dBRjk(i,2,j,1)&
               &                       - Bf(j,1) * Co * Bf(i,3) * Rij1)&
               &                       / (Si * Rij1)
            dBf(i,3,j,3) = (-Si * Bf(i,3) * BRjk(j,2)&
               &                        + Co * dBRjk(i,2,j,2)&
               &                        - Bf(j,3) * (Co * Bf(i,3) * Rjk1&
               &                        + Si * BRjk(i,2))) / (Si * Rjk1)
            !
            dBf(j,1,i,1) = dBf(i,1,j,1)
            dBf(j,3,i,1) = dBf(i,1,j,3)
            dBf(j,1,i,3) = dBf(i,3,j,1)
            dBf(j,3,i,3) = dBf(i,3,j,3)
            !
            dBf(i,1,j,2) = -(dBf(i,1,j,1) + dBf(i,1,j,3))
            dBf(j,2,i,1) = dBf(i,1,j,2)
            dBf(j,1,i,2) = -(dBf(j,1,i,1) + dBf(j,1,i,3))
            dBf(i,2,j,1) = dBf(j,1,i,2)
            dBf(i,3,j,2) = -(dBf(i,3,j,1) + dBf(i,3,j,3))
            dBf(j,2,i,3) = dBf(i,3,j,2)
            dBf(j,3,i,2) = -(dBf(j,3,i,1) + dBf(j,3,i,3))
            dBf(i,2,j,3) = dBf(j,3,i,2)
            !
            dBf(i,2,j,2) = -(dBf(i,2,j,1) + dBf(i,2,j,3))
            dBf(j,2,i,2) = dBf(i,2,j,2)
            !
          End Do
        End Do
        !        Call RecPrt('dBf','(9F9.1)',dBf,9,9)
        !
      End If
      !
      !     Call QExit('Bend')
      Return
    End subroutine bend
    Function arSin(Arg)
      Implicit Real * 8(a - h,o - z)
      Real * 8 ArSin

      A = Arg
      IF (ABS(A) .GT. One) Then
        PRINT 3,A
3       FORMAT(1X,'Warning argument of aSin= ',1F21.18)
        A = Sign(One,A)
      End If
      !
      ArSin = ASin(A)
      Return
    End function arSin
    Function arCos(Arg)
      Implicit Real(wp) (a - h,o - z)
      Real(wp) :: ArCos
      A = Arg
      IF (ABS(A) .GT. One) Then
        A = Sign(One,A)
      End If
      ArCos = ACos(A)
      Return
    End function arCos
  End subroutine trsn

!========================================================================================!
end module modelhessian_module
