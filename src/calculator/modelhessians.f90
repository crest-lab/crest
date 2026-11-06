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
module modelhessian_core
  use iso_fortran_env,only:wp => real64,stdout => output_unit
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
  real(wp),parameter :: aatoau = 1.0/bohr
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

!&<
  integer, private, parameter :: max_elem = 118
  !> covalent radii (taken from Pyykko and Atsumi, Chem. Eur. J. 15, 2009,
  !  188-197), values for metals decreased by 10 %
  real(wp),parameter :: covrad_2009(max_elem) = aatoau * [ &
  & 0.32_wp,0.46_wp, & ! H,He
  & 1.20_wp,0.94_wp,0.77_wp,0.75_wp,0.71_wp,0.63_wp,0.64_wp,0.67_wp, & ! Li-Ne
  & 1.40_wp,1.25_wp,1.13_wp,1.04_wp,1.10_wp,1.02_wp,0.99_wp,0.96_wp, & ! Na-Ar
  & 1.76_wp,1.54_wp, & ! K,Ca
  &                 1.33_wp,1.22_wp,1.21_wp,1.10_wp,1.07_wp, & ! Sc-
  &                 1.04_wp,1.00_wp,0.99_wp,1.01_wp,1.09_wp, & ! -Zn
  &                 1.12_wp,1.09_wp,1.15_wp,1.10_wp,1.14_wp,1.17_wp, & ! Ga-Kr
  & 1.89_wp,1.67_wp, & ! Rb,Sr
  &                 1.47_wp,1.39_wp,1.32_wp,1.24_wp,1.15_wp, & ! Y-
  &                 1.13_wp,1.13_wp,1.08_wp,1.15_wp,1.23_wp, & ! -Cd
  &                 1.28_wp,1.26_wp,1.26_wp,1.23_wp,1.32_wp,1.31_wp, & ! In-Xe
  & 2.09_wp,1.76_wp, & ! Cs,Ba
  &         1.62_wp,1.47_wp,1.58_wp,1.57_wp,1.56_wp,1.55_wp,1.51_wp, & ! La-Eu
  &         1.52_wp,1.51_wp,1.50_wp,1.49_wp,1.49_wp,1.48_wp,1.53_wp, & ! Gd-Yb
  &                 1.46_wp,1.37_wp,1.31_wp,1.23_wp,1.18_wp, & ! Lu-
  &                 1.16_wp,1.11_wp,1.12_wp,1.13_wp,1.32_wp, & ! -Hg
  &                 1.30_wp,1.30_wp,1.36_wp,1.31_wp,1.38_wp,1.42_wp, & ! Tl-Rn
  & 2.01_wp,1.81_wp, & ! Fr,Ra
  &      1.67_wp,1.58_wp,1.52_wp,1.53_wp,1.54_wp,1.55_wp,1.49_wp, & ! Ac-Am
  &      1.49_wp,1.51_wp,1.51_wp,1.48_wp,1.50_wp,1.56_wp,1.58_wp, & ! Cm-No
  &                 1.45_wp,1.41_wp,1.34_wp,1.29_wp,1.27_wp, & ! Lr-
  &                 1.21_wp,1.16_wp,1.15_wp,1.09_wp,1.22_wp, & ! -Cn
  &                 1.36_wp,1.43_wp,1.46_wp,1.58_wp,1.48_wp,1.57_wp ] ! Nh-Og
!&>

  public :: ddvopt,mh_lindh_d2,mh_lindh,mh_swart

!==============================================================================!
contains  !> MODULE PROCEDURES START HERE
!==============================================================================!

  subroutine ddvopt(Cart,nAtoms,Hess,iANr,mhset)
!***********************************************************
!* subroutine ddvopt
!* generates a Lindh Model Hessian
!* Chem. Phys. Let. 241(1995) 423-428
!*
!* Input:
!*     Cart  - cartesian coordinates
!*    nAtoms - number of atoms
!*     iANr  - atom types as integers
!*    mhset  - model Hessian parameters
!*
!* Output:
!*     Hess  - the (packed) model Hessian
!**********************************************************
    Implicit Integer(i-n)
    Implicit Real(wp) (a-h,o-z)
    type(mhparam) :: mhset

    real(wp) :: s6
    real(wp) :: rcut

    Real(wp) :: Cart(3,nAtoms),rij(3),rjk(3),rkl(3), &
   &       Hess((3*nAtoms)*(3*nAtoms+1)/2),si(3),sj(3),sk(3), &
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
    n3 = 3*nAtoms
    Hess = 0.0d0

!
!     Hessian for tension
!
    Do kAtom = 1,nAtoms
      kr = iTabRow(iANr(kAtom))
!        If (kr.eq.0) Go To 5

      Do lAtom = 1,kAtom-1
        lr = iTabRow(iANr(lAtom))
!           If (lr.eq.0) Go To 10
        xkl = Cart(1,kAtom)-Cart(1,lAtom)
        ykl = Cart(2,kAtom)-Cart(2,lAtom)
        zkl = Cart(3,kAtom)-Cart(3,lAtom)
        rkl2 = xkl**2+ykl**2+zkl**2
        r0 = rAv(kr,lr)
        alpha = aAv(kr,lr)

!cccccc VDWx ccccccccccccccccccccccccccccccccc
        c6k = c6(iANr(katom))
        c6l = c6(iANr(latom))
        c66 = sqrt(c6k*c6l)
        Rv = (vander(iANr(katom))+vander(iANr(latom)))/bohr

        call getvdwxx(xkl,ykl,zkl,c66,s6,Rv,vdw(1,1))
        call getvdwxy(xkl,ykl,zkl,c66,s6,Rv,vdw(1,2))
        call getvdwxy(xkl,zkl,ykl,c66,s6,Rv,vdw(1,3))
        call getvdwxx(ykl,xkl,zkl,c66,s6,Rv,vdw(2,2))
        call getvdwxy(ykl,zkl,xkl,c66,s6,Rv,vdw(2,3))
        call getvdwxx(zkl,xkl,ykl,c66,s6,Rv,vdw(3,3))
!cccccc Ende VDWx ccccccccccccccccccccccccccccccc

        gamma = rkr*Exp(alpha*r0**2)
! not better: *sqrt(abs(wb(kAtom,lAtom)))
        gmm = gamma*Exp(-alpha*rkl2)
        Hxx = gmm*xkl*xkl/rkl2-vdw(1,1)
        Hxy = gmm*xkl*ykl/rkl2-vdw(1,2)
        Hxz = gmm*xkl*zkl/rkl2-vdw(1,3)
        Hyy = gmm*ykl*ykl/rkl2-vdw(2,2)
        Hyz = gmm*ykl*zkl/rkl2-vdw(2,3)
        Hzz = gmm*zkl*zkl/rkl2-vdw(3,3)

!
        Hess(Ind(1,kAtom,1,kAtom)) = Hess(Ind(1,kAtom,1,kAtom))+Hxx
        Hess(Ind(2,kAtom,1,kAtom)) = Hess(Ind(2,kAtom,1,kAtom))+Hxy
        Hess(Ind(2,kAtom,2,kAtom)) = Hess(Ind(2,kAtom,2,kAtom))+Hyy
        Hess(Ind(3,kAtom,1,kAtom)) = Hess(Ind(3,kAtom,1,kAtom))+Hxz
        Hess(Ind(3,kAtom,2,kAtom)) = Hess(Ind(3,kAtom,2,kAtom))+Hyz
        Hess(Ind(3,kAtom,3,kAtom)) = Hess(Ind(3,kAtom,3,kAtom))+Hzz
!
        Hess(Ind(1,kAtom,1,lAtom)) = Hess(Ind(1,kAtom,1,lAtom))-Hxx
        Hess(Ind(1,kAtom,2,lAtom)) = Hess(Ind(1,kAtom,2,lAtom))-Hxy
        Hess(Ind(1,kAtom,3,lAtom)) = Hess(Ind(1,kAtom,3,lAtom))-Hxz
        Hess(Ind(2,kAtom,1,lAtom)) = Hess(Ind(2,kAtom,1,lAtom))-Hxy
        Hess(Ind(2,kAtom,2,lAtom)) = Hess(Ind(2,kAtom,2,lAtom))-Hyy
        Hess(Ind(2,kAtom,3,lAtom)) = Hess(Ind(2,kAtom,3,lAtom))-Hyz
        Hess(Ind(3,kAtom,1,lAtom)) = Hess(Ind(3,kAtom,1,lAtom))-Hxz
        Hess(Ind(3,kAtom,2,lAtom)) = Hess(Ind(3,kAtom,2,lAtom))-Hyz
        Hess(Ind(3,kAtom,3,lAtom)) = Hess(Ind(3,kAtom,3,lAtom))-Hzz
!
        Hess(Ind(1,lAtom,1,lAtom)) = Hess(Ind(1,lAtom,1,lAtom))+Hxx
        Hess(Ind(2,lAtom,1,lAtom)) = Hess(Ind(2,lAtom,1,lAtom))+Hxy
        Hess(Ind(2,lAtom,2,lAtom)) = Hess(Ind(2,lAtom,2,lAtom))+Hyy
        Hess(Ind(3,lAtom,1,lAtom)) = Hess(Ind(3,lAtom,1,lAtom))+Hxz
        Hess(Ind(3,lAtom,2,lAtom)) = Hess(Ind(3,lAtom,2,lAtom))+Hyz
        Hess(Ind(3,lAtom,3,lAtom)) = Hess(Ind(3,lAtom,3,lAtom))+Hzz
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
        Do jAtom = 1,iAtom-1
          If (jAtom .eq. mAtom) Go To 40
          jr = iTabRow(iANr(jAtom))
!           If (jr.eq.0) Go To 40
          if (rcutoff(cart,jatom,iatom,rcut)) cycle
          if (rcutoff(cart,jatom,matom,rcut)) cycle
!           if(wb(jatom,iatom).lt.wthr) cycle
!           if(wb(jatom,matom).lt.wthr) cycle

          xmi = (Cart(1,iAtom)-Cart(1,mAtom))
          ymi = (Cart(2,iAtom)-Cart(2,mAtom))
          zmi = (Cart(3,iAtom)-Cart(3,mAtom))
          rmi2 = xmi**2+ymi**2+zmi**2
          rmi = sqrt(rmi2)
          r0mi = rAv(mr,ir)
          ami = aAv(mr,ir)
!
          xmj = (Cart(1,jAtom)-Cart(1,mAtom))
          ymj = (Cart(2,jAtom)-Cart(2,mAtom))
          zmj = (Cart(3,jAtom)-Cart(3,mAtom))
          rmj2 = xmj**2+ymj**2+zmj**2
          rmj = sqrt(rmj2)
          r0mj = rAv(mr,jr)
          amj = aAv(mr,jr)
!
!---------- Test if zero angle
!
          Test = xmi*xmj+ymi*ymj+zmi*zmj
          Test = Test/(rmi*rmj)
          If (Test .eq. One) Go To 40
!
          xij = (Cart(1,jAtom)-Cart(1,iAtom))
          yij = (Cart(2,jAtom)-Cart(2,iAtom))
          zij = (Cart(3,jAtom)-Cart(3,iAtom))
          rij2 = xij**2+yij**2+zij**2
          rrij = sqrt(rij2)
!
          alpha = rkf*exp((ami*r0mi**2+amj*r0mj**2))
!
          r = sqrt(rmj2+rmi2)
          gij = alpha*exp(-(ami*rmi2+amj*rmj2))
!           Write (*,*) ' gij=',gij
          rL2 = (ymi*zmj-zmi*ymj)**2+(zmi*xmj-xmi*zmj)**2+ &
   &         (xmi*ymj-ymi*xmj)**2
!hjw modified
          if (rL2 .lt. 1.d-14) then
            rL = 0
          else
            rL = sqrt(rL2)
          end if
!
          if ((rmj .gt. rZero).and.(rmi .gt. rZero).and. &
   &                                (rrij .gt. rZero)) Then
            SinPhi = rL/(rmj*rmi)
            rmidotrmj = xmi*xmj+ymi*ymj+zmi*zmj
            CosPhi = rmidotrmj/(rmj*rmi)
!
!-------------None linear case
!
            If (SinPhi .gt. rZero) Then
!               Write (*,*) ' None linear case'
              si(1) = (xmi/rmi*cosphi-xmj/rmj)/(rmi*sinphi)
              si(2) = (ymi/rmi*cosphi-ymj/rmj)/(rmi*sinphi)
              si(3) = (zmi/rmi*cosphi-zmj/rmj)/(rmi*sinphi)
              sj(1) = (cosphi*xmj/rmj-xmi/rmi)/(rmj*sinphi)
              sj(2) = (cosphi*ymj/rmj-ymi/rmi)/(rmj*sinphi)
              sj(3) = (cosphi*zmj/rmj-zmi/rmi)/(rmj*sinphi)
              sm(1) = -si(1)-sj(1)
              sm(2) = -si(2)-sj(2)
              sm(3) = -si(3)-sj(3)
              Do icoor = 1,3
                Do jCoor = 1,3
                  If (mAtom .gt. iAtom) Then
                    Hess(Ind(icoor,mAtom,jcoor,iAtom)) = &
  &                        Hess(Ind(icoor,mAtom,jcoor,iAtom)) &
  &                        +gij*sm(icoor)*si(jcoor)
                  else
                    Hess(Ind(icoor,iAtom,jcoor,mAtom)) = &
   &                        Hess(Ind(icoor,iAtom,jcoor,mAtom)) &
   &                        +gij*si(icoor)*sm(jcoor)
                  End If
                  If (mAtom .gt. jAtom) Then
                    Hess(Ind(icoor,mAtom,jcoor,jAtom)) = &
 &                        Hess(Ind(icoor,mAtom,jcoor,jAtom)) &
 &                        +gij*sm(icoor)*sj(jcoor)
                  else
                    Hess(Ind(icoor,jAtom,jcoor,mAtom)) = &
   &                        Hess(Ind(icoor,jAtom,jcoor,mAtom)) &
   &                        +gij*sj(icoor)*sm(jcoor)
                  End If
                  If (iAtom .gt. jAtom) Then
                    Hess(Ind(icoor,iAtom,jcoor,jAtom)) = &
 &                        Hess(Ind(icoor,iAtom,jcoor,jAtom)) &
 &                        +gij*si(icoor)*sj(jcoor)
                  else
                    Hess(Ind(icoor,jAtom,jcoor,iAtom)) = &
 &                        Hess(Ind(icoor,jAtom,jcoor,iAtom)) &
 &                        +gij*sj(icoor)*si(jcoor)
                  End If
                End Do
              End Do
              Do icoor = 1,3
                Do jCoor = 1,icoor
                  Hess(Ind(icoor,iAtom,jcoor,iAtom)) = &
   &                        Hess(Ind(icoor,iAtom,jcoor,iAtom)) &
   &                        +gij*si(icoor)*si(jcoor)
                  Hess(Ind(icoor,mAtom,jcoor,mAtom)) = &
   &                        Hess(Ind(icoor,mAtom,jcoor,mAtom)) &
   &                        +gij*sm(icoor)*sm(jcoor)
                  Hess(Ind(icoor,jAtom,jcoor,jAtom)) = &
   &                        Hess(Ind(icoor,jAtom,jcoor,jAtom)) &
   &                        +gij*sj(icoor)*sj(jcoor)

!
                End Do
              End Do
            Else
!
!----------------Linear case
!
              if ((abs(ymi) .gt. rZero).or. &
&                 (abs(xmi) .gt. rZero)) Then
                x(1) = -ymi
                y(1) = xmi
                z(1) = Zero
                x(2) = -xmi*zmi
                y(2) = -ymi*zmi
                z(2) = xmi*xmi+ymi*ymi
              Else
                x(1) = One
                y(1) = Zero
                z(1) = Zero
                x(2) = Zero
                y(2) = One
                z(2) = Zero
              End If
              Do i = 1,2
                r1 = sqrt(x(i)**2+y(i)**2+z(i)**2)
                cosThetax = x(i)/r1
                cosThetay = y(i)/r1
                cosThetaz = z(i)/r1
                si(1) = -cosThetax/rmi
                si(2) = -cosThetay/rmi
                si(3) = -cosThetaz/rmi
                sj(1) = -cosThetax/rmj
                sj(2) = -cosThetay/rmj
                sj(3) = -cosThetaz/rmj
                sm(1) = -(si(1)+sj(1))
                sm(2) = -(si(2)+sj(2))
                sm(3) = -(si(3)+sj(3))
!
                Do icoor = 1,3
                  Do jCoor = 1,3
                    If (mAtom .gt. iAtom) Then
                      Hess(Ind(icoor,mAtom,jcoor,iAtom)) = &
 &                        Hess(Ind(icoor,mAtom,jcoor,iAtom)) &
 &                         +gij*sm(icoor)*si(jcoor)
                    else
                      Hess(Ind(icoor,iAtom,jcoor,mAtom)) = &
&                        Hess(Ind(icoor,iAtom,jcoor,mAtom)) &
&                         +gij*si(icoor)*sm(jcoor)
                    End If
                    If (mAtom .gt. jAtom) Then
                      Hess(Ind(icoor,mAtom,jcoor,jAtom)) = &
 &                        Hess(Ind(icoor,mAtom,jcoor,jAtom)) &
 &                         +gij*sm(icoor)*sj(jcoor)
                    else
                      Hess(Ind(icoor,jAtom,jcoor,mAtom)) = &
 &                        Hess(Ind(icoor,jAtom,jcoor,mAtom)) &
 &                         +gij*sj(icoor)*sm(jcoor)
                    End If
                    If (iAtom .gt. jAtom) Then
                      Hess(Ind(icoor,iAtom,jcoor,jAtom)) = &
&                        Hess(Ind(icoor,iAtom,jcoor,jAtom)) &
&                         +gij*si(icoor)*sj(jcoor)
                    else
                      Hess(Ind(icoor,jAtom,jcoor,iAtom)) = &
&                        Hess(Ind(icoor,jAtom,jcoor,iAtom)) &
&                         +gij*sj(icoor)*si(jcoor)
                    End If
                  End Do
                End Do
                Do icoor = 1,3
                  Do jCoor = 1,icoor
                    Hess(Ind(icoor,iAtom,jcoor,iAtom)) = &
&                        Hess(Ind(icoor,iAtom,jcoor,iAtom)) &
&                         +gij*si(icoor)*si(jcoor)
                    Hess(Ind(icoor,mAtom,jcoor,mAtom)) = &
&                        Hess(Ind(icoor,mAtom,jcoor,mAtom)) &
&                         +gij*sm(icoor)*sm(jcoor)
                    Hess(Ind(icoor,jAtom,jcoor,jAtom)) = &
&                         Hess(Ind(icoor,jAtom,jcoor,jAtom)) &
&                         +gij*sj(icoor)*sj(jcoor)
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
          ij_ = nAtoms*(jAtom-1)+iAtom
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
            lk_ = nAtoms*(kAtom-1)+lAtom
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
            rij(1) = Cart(1,iAtom)-Cart(1,jAtom)
            rij(2) = Cart(2,iAtom)-Cart(2,jAtom)
            rij(3) = Cart(3,iAtom)-Cart(3,jAtom)
            rij0 = rAv(ir,jr)**2
            aij = aAv(ir,jr)
!
            rjk(1) = Cart(1,jAtom)-Cart(1,kAtom)
            rjk(2) = Cart(2,jAtom)-Cart(2,kAtom)
            rjk(3) = Cart(3,jAtom)-Cart(3,kAtom)
            rjk0 = rAv(jr,kr)**2
            ajk = aAv(jr,kr)
!
            rkl(1) = Cart(1,kAtom)-Cart(1,lAtom)
            rkl(2) = Cart(2,kAtom)-Cart(2,lAtom)
            rkl(3) = Cart(3,kAtom)-Cart(3,lAtom)
            rkl0 = rAv(kr,lr)**2
            akl = aAv(kr,lr)
!
            rij2 = rij(1)**2+rij(2)**2+rij(3)**2
            rjk2 = rjk(1)**2+rjk(2)**2+rjk(3)**2
            rkl2 = rkl(1)**2+rkl(2)**2+rkl(3)**2
!              Allow only angles in the range of 35-145
            A35 = (35.0D0/180.D0)*Pi
            CosFi_Max = Cos(A35)
            CosFi2 = (rij(1)*rjk(1)+rij(2)*rjk(2)+rij(3)*rjk(3)) &
  &               /Sqrt(rij2*rjk2)
            If (Abs(CosFi2) .gt. CosFi_Max) Go To 222
            CosFi3 = (rkl(1)*rjk(1)+rkl(2)*rjk(2)+rkl(3)*rjk(3)) &
  &               /Sqrt(rkl2*rjk2)
            If (Abs(CosFi3) .gt. CosFi_Max) Go To 222

            beta = rkt* &
  &                       exp((aij*rij0+ajk*rjk0+akl*rkl0))
            tij = beta*exp(-(aij*rij2+ajk*rjk2+akl*rkl2))

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
    &            +tij*si(icoor)*sj(jcoor)
                Hess(Ind(icoor,iAtom,jcoor,kAtom)) = &
    &           Hess(Ind(icoor,iAtom,jcoor,kAtom)) &
    &            +tij*si(icoor)*sk(jcoor)
                Hess(Ind(icoor,iAtom,jcoor,lAtom)) = &
    &           Hess(Ind(icoor,iAtom,jcoor,lAtom)) &
    &            +tij*si(icoor)*sl(jcoor)
                Hess(Ind(icoor,jAtom,jcoor,kAtom)) = &
    &           Hess(Ind(icoor,jAtom,jcoor,kAtom)) &
    &            +tij*sj(icoor)*sk(jcoor)
                Hess(Ind(icoor,jAtom,jcoor,lAtom)) = &
    &           Hess(Ind(icoor,jAtom,jcoor,lAtom)) &
    &            +tij*sj(icoor)*sl(jcoor)
                Hess(Ind(icoor,kAtom,jcoor,lAtom)) = &
    &           Hess(Ind(icoor,kAtom,jcoor,lAtom)) &
    &            +tij*sk(icoor)*sl(jcoor)

              End Do
            End Do
!
!-------------Diagonal block
!
            Do icoor = 1,3
              Do jCoor = 1,icoor
                Hess(Ind(icoor,iAtom,jcoor,iAtom)) = &
    &           Hess(Ind(icoor,iAtom,jcoor,iAtom)) &
    &            +tij*si(icoor)*si(jcoor)
                Hess(Ind(icoor,jAtom,jcoor,jAtom)) = &
    &           Hess(Ind(icoor,jAtom,jcoor,jAtom)) &
    &            +tij*sj(icoor)*sj(jcoor)
                Hess(Ind(icoor,kAtom,jcoor,kAtom)) = &
    &           Hess(Ind(icoor,kAtom,jcoor,kAtom)) &
    &            +tij*sk(icoor)*sk(jcoor)
                Hess(Ind(icoor,lAtom,jcoor,lAtom)) = &
    &           Hess(Ind(icoor,lAtom,jcoor,lAtom)) &
    &            +tij*sl(icoor)*sl(jcoor)

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
      ixyz = (iatom-1)*3+i
    end function ixyz
    function jnd(i,j)
      integer :: jnd
      integer,intent(in) :: i,j
      jnd = i*(i-1)/2+j
    end function jnd
    function ind(i,iatom,j,jatom)
      integer :: ind
      integer,intent(in) :: i,iatom,j,jatom
      ind = jnd(max(ixyz(i,iatom),ixyz(j,jatom)),min(ixyz(i,iatom),ixyz(j,jatom)))
    end function ind
  end subroutine ddvopt

!========================================================================================!
!########################################################################################!
!========================================================================================!

  subroutine mh_swart(xyz,n,hess,at,modh)
!****************************************************************************
!* Swart's Model Hessian augmented with D2
!* ------------------------------------------------------------------------
!* Implemented after:
!* M. Swart, F. M. Bickelhaupt, Int. J. Quantum Chem., 2006, 106, 2536–2544.
!* DOI:10.1002/qua.21049
!*
!* gij = exp[-(Rij/Cij-1)]
!* kij   = rkr·gij
!* kijk  = rkf·gij·gjk
!* kijkl = rkt·gij·gjk·gkl
!*
!* The proposed force constants by Swart are:
!* rkr = 0.35, rkf = 0.15, rkt = 0.005
!*
!* This Hessian is additionally augmented with D2, please note that D2
!* is not implemented in atomic units and requires some magical conversion
!* factor somewhere hidden in the implementation below.
!****************************************************************************
    implicit none

    integer,intent(in)  :: n
    real(wp),intent(in)  :: xyz(3,n)
    real(wp),intent(out) :: hess((3*n)*(3*n+1)/2)
    integer,intent(in)  :: at(n)
    type(mhparam),intent(in) :: modh

    integer  :: n3
    real(wp),parameter :: rzero = 1.0e-10_wp
    logical,allocatable :: lcutoff(:,:)
    real(wp) :: kd

    allocate (lcutoff(n,n),source=.false.)

    n3 = 3*n
    hess = 0.0d0

!  the dispersion force constant is used relative to the stretch force constant
    kd = modh%kd/modh%kr

    associate (rad => covrad_2009)

      call mh_swart_stretch(n,at,xyz,hess,modh%kr,kd,modh%s6,rad,rad,lcutoff,modh%rcut)
      if (modh%kf .ne. 0.0_wp) &
        call mh_swart_bend(n,at,xyz,hess,modh%kf,kd,rad,rad,lcutoff)
      if (modh%kt .ne. 0.0_wp) &
        call mh_swart_torsion(n,at,xyz,hess,modh%kt,kd,rad,rad,lcutoff)
      if (modh%ko .ne. 0.0_wp) &
        call mh_swart_outofp(n,at,xyz,hess,modh%ko,kd,rad,rad,lcutoff)
      if (modh%kq .ne. 0.0_wp) then
!      call new_charge_model_2019(chrgeq,n,at)
        call mh_eeq(n,at,xyz,0.0_wp,modh%kq,hess)
      end if

    end associate

  end subroutine mh_swart

  pure subroutine mh_swart_stretch(n,at,xyz,hess,kr,kd,s6,rcov,rvdw,lcutoff,rcut)
    implicit none

    integer,intent(in)    :: n
    integer,intent(in)    :: at(n)
    real(wp),intent(in)    :: xyz(3,n)
    real(wp),intent(inout) :: hess((3*n)*(3*n+1)/2)
    real(wp),intent(in)    :: kr
    real(wp),intent(in)    :: kd
    real(wp),intent(in)    :: s6
    real(wp),intent(in)    :: rcov(:)
    real(wp),intent(in)    :: rvdw(:)
    logical,intent(out)   :: lcutoff(n,n)
    real(wp),intent(in)    :: rcut

    integer  :: i,j
    real(wp) :: xij,yij,zij,rij2,r0,d0
    real(wp) :: gmm
    real(wp) :: c6i,c6j,c6ij,rv
    real(wp) :: hxx,hxy,hxz,hyy,hyz,hzz
    real(wp) :: vdw(3,3)

!! ------------------------------------------------------------------------
!  Hessian for stretch
!! ------------------------------------------------------------------------
    stretch_iAt: do i = 1,n

      stretch_jAt: do j = 1,i-1

        ! save for later
        lcutoff(i,j) = rcutoff(xyz,i,j,rcut)
        lcutoff(j,i) = lcutoff(i,j)

        xij = xyz(1,i)-xyz(1,j)
        yij = xyz(2,i)-xyz(2,j)
        zij = xyz(3,i)-xyz(3,j)
        rij2 = xij**2+yij**2+zij**2
        r0 = rcov(at(i))+rcov(at(j))
        d0 = rvdw(at(i))+rvdw(at(j))

        !cccccc vdwx ccccccccccccccccccccccccccccccccc
        c6i = c6(at(i))
        c6j = c6(at(j))
        c6ij = sqrt(c6i*c6j)
        rv = (vander(at(i))+vander(at(j)))*aatoau

        call getvdwxx(xij,yij,zij,c6ij,s6,rv,vdw(1,1))
        call getvdwxy(xij,yij,zij,c6ij,s6,rv,vdw(1,2))
        call getvdwxy(xij,zij,yij,c6ij,s6,rv,vdw(1,3))
        call getvdwxx(yij,xij,zij,c6ij,s6,rv,vdw(2,2))
        call getvdwxy(yij,zij,xij,c6ij,s6,rv,vdw(2,3))
        call getvdwxx(zij,xij,yij,c6ij,s6,rv,vdw(3,3))
        !cccccc ende vdwx ccccccccccccccccccccccccccccccc

        gmm = kr*fk_swart(1.0_wp,r0,rij2) &
              +kr*kd*fk_vdw(5.0_wp,d0,rij2)

        !gmm = max(gmm,min_fk)

        hxx = gmm*xij*xij/rij2-vdw(1,1)
        hxy = gmm*xij*yij/rij2-vdw(1,2)
        hxz = gmm*xij*zij/rij2-vdw(1,3)
        hyy = gmm*yij*yij/rij2-vdw(2,2)
        hyz = gmm*yij*zij/rij2-vdw(2,3)
        hzz = gmm*zij*zij/rij2-vdw(3,3)

        ! save diagonal elements for atom i
        hess(ind(1,i,1,i)) = hess(ind(1,i,1,i))+hxx
        hess(ind(2,i,1,i)) = hess(ind(2,i,1,i))+hxy
        hess(ind(2,i,2,i)) = hess(ind(2,i,2,i))+hyy
        hess(ind(3,i,1,i)) = hess(ind(3,i,1,i))+hxz
        hess(ind(3,i,2,i)) = hess(ind(3,i,2,i))+hyz
        hess(ind(3,i,3,i)) = hess(ind(3,i,3,i))+hzz
        ! save elements between atom i and atom j
        hess(ind(1,i,1,j)) = hess(ind(1,i,1,j))-hxx
        hess(ind(1,i,2,j)) = hess(ind(1,i,2,j))-hxy
        hess(ind(1,i,3,j)) = hess(ind(1,i,3,j))-hxz
        hess(ind(2,i,1,j)) = hess(ind(2,i,1,j))-hxy
        hess(ind(2,i,2,j)) = hess(ind(2,i,2,j))-hyy
        hess(ind(2,i,3,j)) = hess(ind(2,i,3,j))-hyz
        hess(ind(3,i,1,j)) = hess(ind(3,i,1,j))-hxz
        hess(ind(3,i,2,j)) = hess(ind(3,i,2,j))-hyz
        hess(ind(3,i,3,j)) = hess(ind(3,i,3,j))-hzz
        ! save diagonal elements for atom j
        hess(ind(1,j,1,j)) = hess(ind(1,j,1,j))+hxx
        hess(ind(2,j,1,j)) = hess(ind(2,j,1,j))+hxy
        hess(ind(2,j,2,j)) = hess(ind(2,j,2,j))+hyy
        hess(ind(3,j,1,j)) = hess(ind(3,j,1,j))+hxz
        hess(ind(3,j,2,j)) = hess(ind(3,j,2,j))+hyz
        hess(ind(3,j,3,j)) = hess(ind(3,j,3,j))+hzz

      end do stretch_jAt
    end do stretch_iAt

  end subroutine mh_swart_stretch

  pure subroutine mh_swart_bend(n,at,xyz,hess,kf,kd,rcov,rvdw,lcutoff)
    implicit none

    integer,intent(in)    :: n
    integer,intent(in)    :: at(n)
    real(wp),intent(in)    :: xyz(3,n)
    real(wp),intent(inout) :: hess((3*n)*(3*n+1)/2)
    real(wp),intent(in)    :: kf
    real(wp),intent(in)    :: kd
    real(wp),intent(in)    :: rcov(:)
    real(wp),intent(in)    :: rvdw(:)
    logical,intent(in)    :: lcutoff(n,n)

    integer  :: i,j,m,ic,jc,ii
    real(wp),parameter :: rzero = 1.0e-10_wp
    real(wp) :: xij,yij,zij,rij2,rrij,r1
    real(wp) :: xmi,ymi,zmi,rmi2,rmi,r0mi,d0mj,gmi
    real(wp) :: xmj,ymj,zmj,rmj2,rmj,r0mj,d0mi,gmj
    real(wp) :: test,gij,rl2,rl,rmidotrmj
    real(wp) :: sinphi,cosphi,costhetax,costhetay,costhetaz
    real(wp) :: alpha
    real(wp) :: si(3),sj(3),sm(3),x(2),y(2),z(2)

!! ------------------------------------------------------------------------
!  Hessian for bending
!! ------------------------------------------------------------------------
    bend_mAt: do m = 1,n
      bend_iAt: do i = 1,n
        if (i .eq. m) cycle bend_iAt
        if (lcutoff(i,m)) cycle bend_iAt

        xmi = (xyz(1,i)-xyz(1,m))
        ymi = (xyz(2,i)-xyz(2,m))
        zmi = (xyz(3,i)-xyz(3,m))
        rmi2 = xmi**2+ymi**2+zmi**2
        rmi = sqrt(rmi2)
        r0mi = rcov(at(m))+rcov(at(i))
        d0mi = rvdw(at(m))+rvdw(at(i))

        bend_jAt: do j = 1,i-1
          if (j .eq. m) cycle bend_jAt
          if (lcutoff(j,i)) cycle bend_jAt
          if (lcutoff(j,m)) cycle bend_jAt

          xmj = (xyz(1,j)-xyz(1,m))
          ymj = (xyz(2,j)-xyz(2,m))
          zmj = (xyz(3,j)-xyz(3,m))
          rmj2 = xmj**2+ymj**2+zmj**2
          rmj = sqrt(rmj2)
          r0mj = rcov(at(m))+rcov(at(j))
          d0mj = rvdw(at(m))+rvdw(at(j))

          ! test if zero angle
          test = xmi*xmj+ymi*ymj+zmi*zmj
          test = test/(rmi*rmj)
          if (abs(test-1.0_wp) .lt. 1.0e-12_wp) cycle bend_jAt

          xij = (xyz(1,j)-xyz(1,i))
          yij = (xyz(2,j)-xyz(2,i))
          zij = (xyz(3,j)-xyz(3,i))
          rij2 = xij**2+yij**2+zij**2
          rrij = sqrt(rij2)

          gmi = fk_swart(1.0_wp,r0mi,rmi2) &
                +0.5_wp*kd*fk_vdw(5.0_wp,d0mi,rmi2)
          gmj = fk_swart(1.0_wp,r0mj,rmj2) &
                +0.5_wp*kd*fk_vdw(5.0_wp,d0mj,rmj2)

          gij = kf*gmi*gmj

          rl2 = (ymi*zmj-zmi*ymj)**2+(zmi*xmj-xmi*zmj)**2+(xmi*ymj-ymi*xmj)**2

          if (rl2 .lt. 1.e-14_wp) then
            rl = 0.0_wp
          else
            rl = sqrt(rl2)
          end if

          !gij = max(gij,min_fk)

          if ((rmj .gt. rzero).and.(rmi .gt. rzero).and.(rrij .gt. rzero)) then
            sinphi = rl/(rmj*rmi)
            rmidotrmj = xmi*xmj+ymi*ymj+zmi*zmj
            cosphi = rmidotrmj/(rmj*rmi)
            ! none linear case
            if (sinphi .gt. rzero) then
              si(1) = (xmi/rmi*cosphi-xmj/rmj)/(rmi*sinphi)
              si(2) = (ymi/rmi*cosphi-ymj/rmj)/(rmi*sinphi)
              si(3) = (zmi/rmi*cosphi-zmj/rmj)/(rmi*sinphi)
              sj(1) = (cosphi*xmj/rmj-xmi/rmi)/(rmj*sinphi)
              sj(2) = (cosphi*ymj/rmj-ymi/rmi)/(rmj*sinphi)
              sj(3) = (cosphi*zmj/rmj-zmi/rmi)/(rmj*sinphi)
              sm(1) = -si(1)-sj(1)
              sm(2) = -si(2)-sj(2)
              sm(3) = -si(3)-sj(3)
              do ic = 1,3
                do jc = 1,3
                  if (m .gt. i) then
                    hess(ind(ic,m,jc,i)) = hess(ind(ic,m,jc,i)) &
                                           +gij*sm(ic)*si(jc)
                  else
                    hess(ind(ic,i,jc,m)) = hess(ind(ic,i,jc,m)) &
                                           +gij*si(ic)*sm(jc)
                  end if
                  if (m .gt. j) then
                    hess(ind(ic,m,jc,j)) = hess(ind(ic,m,jc,j)) &
                                           +gij*sm(ic)*sj(jc)
                  else
                    hess(ind(ic,j,jc,m)) = hess(ind(ic,j,jc,m)) &
                                           +gij*sj(ic)*sm(jc)
                  end if
                  if (i .gt. j) then
                    hess(ind(ic,i,jc,j)) = hess(ind(ic,i,jc,j)) &
                                           +gij*si(ic)*sj(jc)
                  else
                    hess(ind(ic,j,jc,i)) = hess(ind(ic,j,jc,i)) &
                                           +gij*sj(ic)*si(jc)
                  end if
                end do
              end do
              do ic = 1,3
                do jc = 1,ic
                  hess(ind(ic,i,jc,i)) = hess(ind(ic,i,jc,i))+gij*si(ic)*si(jc)
                  hess(ind(ic,m,jc,m)) = hess(ind(ic,m,jc,m))+gij*sm(ic)*sm(jc)
                  hess(ind(ic,j,jc,j)) = hess(ind(ic,j,jc,j))+gij*sj(ic)*sj(jc)
                end do
              end do
            else
              ! linear case
              if ((abs(ymi) .gt. rzero).or.(abs(xmi) .gt. rzero)) then
                x(1) = -ymi
                y(1) = xmi
                z(1) = 0.0_wp
                x(2) = -xmi*zmi
                y(2) = -ymi*zmi
                z(2) = xmi*xmi+ymi*ymi
              else
                x(1) = 1.0_wp
                y(1) = 0.0_wp
                z(1) = 0.0_wp
                x(2) = 0.0_wp
                y(2) = 1.0_wp
                z(2) = 0.0_wp
              end if
              do ii = 1,2
                r1 = sqrt(x(ii)**2+y(ii)**2+z(ii)**2)
                costhetax = x(ii)/r1
                costhetay = y(ii)/r1
                costhetaz = z(ii)/r1
                si(1) = -costhetax/rmi
                si(2) = -costhetay/rmi
                si(3) = -costhetaz/rmi
                sj(1) = -costhetax/rmj
                sj(2) = -costhetay/rmj
                sj(3) = -costhetaz/rmj
                sm(1) = -(si(1)+sj(1))
                sm(2) = -(si(2)+sj(2))
                sm(3) = -(si(3)+sj(3))
                !
                do ic = 1,3
                  do jc = 1,3
                    if (m .gt. i) then
                      hess(ind(ic,m,jc,i)) = hess(ind(ic,m,jc,i)) &
                                             +gij*sm(ic)*si(jc)
                    else
                      hess(ind(ic,i,jc,m)) = hess(ind(ic,i,jc,m)) &
                                             +gij*si(ic)*sm(jc)
                    end if
                    if (m .gt. j) then
                      hess(ind(ic,m,jc,j)) = hess(ind(ic,m,jc,j)) &
                                             +gij*sm(ic)*sj(jc)
                    else
                      hess(ind(ic,j,jc,m)) = hess(ind(ic,j,jc,m)) &
                                             +gij*sj(ic)*sm(jc)
                    end if
                    if (i .gt. j) then
                      hess(ind(ic,i,jc,j)) = hess(ind(ic,i,jc,j)) &
                                             +gij*si(ic)*sj(jc)
                    else
                      hess(ind(ic,j,jc,i)) = hess(ind(ic,j,jc,i)) &
                                             +gij*sj(ic)*si(jc)
                    end if
                  end do
                end do
                do ic = 1,3
                  do jc = 1,ic
                    hess(ind(ic,i,jc,i)) = hess(ind(ic,i,jc,i)) &
                                           +gij*si(ic)*si(jc)
                    hess(ind(ic,m,jc,m)) = hess(ind(ic,m,jc,m)) &
                                           +gij*sm(ic)*sm(jc)
                    hess(ind(ic,j,jc,j)) = hess(ind(ic,j,jc,j)) &
                                           +gij*sj(ic)*sj(jc)
                  end do
                end do
              end do

            end if
          end if

        end do bend_jAt
      end do bend_iAt
    end do bend_mAt

  end subroutine mh_swart_bend

  pure subroutine mh_swart_torsion(n,at,xyz,hess,kt,kd,rcov,rvdw,lcutoff)
    implicit none

    integer,intent(in)    :: n
    integer,intent(in)    :: at(n)
    real(wp),intent(in)    :: xyz(3,n)
    real(wp),intent(inout) :: hess((3*n)*(3*n+1)/2)
    real(wp),intent(in)    :: kt
    real(wp),intent(in)    :: kd
    real(wp),intent(in)    :: rcov(:)
    real(wp),intent(in)    :: rvdw(:)
    logical,intent(in)    :: lcutoff(n,n)

    integer  :: i,j,k,l,ic,jc,ij,kl
!  allow only angles in the range of 35-145
    real(wp),parameter :: a35 = (35.0d0/180.d0)*pi
    real(wp),parameter :: cosfi_max = cos(a35)
    real(wp) :: txyz(3,4),c(3,4)
    real(wp) :: rij(3),rij0,aij,rij2,d0ij,gij
    real(wp) :: rjk(3),rjk0,ajk,rjk2,d0jk,gjk
    real(wp) :: rkl(3),rkl0,akl,rkl2,d0kl,gkl
    real(wp) :: cosfi2,cosfi3,cosfi4
    real(wp) :: beta,tij,tau
    real(wp) :: si(3),sj(3),sk(3),sl(3)

!! ------------------------------------------------------------------------
!  Hessian for torsion
!! ------------------------------------------------------------------------
    torsion_jAt: do j = 1,n
      txyz(:,2) = xyz(:,j)
      torsion_kAt: do k = 1,n
        if (k .eq. j) cycle torsion_kAt
        if (lcutoff(k,j)) cycle torsion_kAt
        txyz(:,3) = xyz(:,k)
        torsion_iAt: do i = 1,n
          ij = n*(j-1)+i
          if (i .eq. j) cycle torsion_iAt
          if (i .eq. k) cycle torsion_iAt
          if (lcutoff(i,k)) cycle torsion_iAt
          if (lcutoff(i,j)) cycle torsion_iAt

          txyz(:,1) = xyz(:,i)
          torsion_lAt: do l = 1,n
            kl = n*(l-1)+k
            if (ij .le. kl) cycle torsion_lAt
            if (l .eq. i) cycle torsion_lAt
            if (l .eq. j) cycle torsion_lAt
            if (l .eq. k) cycle torsion_lAt
!
            if (lcutoff(l,i)) cycle torsion_lAt
            if (lcutoff(l,k)) cycle torsion_lAt
            if (lcutoff(l,j)) cycle torsion_lAt

            txyz(:,4) = xyz(:,l)

            rij = xyz(:,i)-xyz(:,j)
            d0ij = rvdw(at(i))+rvdw(at(j))
            rij0 = rcov(at(i))+rcov(at(j))

            rjk = xyz(:,j)-xyz(:,k)
            d0jk = rvdw(at(j))+rvdw(at(k))
            rjk0 = rcov(at(j))+rcov(at(k))

            rkl = xyz(:,k)-xyz(:,l)
            d0kl = rvdw(at(k))+rvdw(at(l))
            rkl0 = rcov(at(k))+rcov(at(l))

            rij2 = sum(rij**2)
            rjk2 = sum(rjk**2)
            rkl2 = sum(rjk**2)

            cosfi2 = dot_product(rij,rjk)/sqrt(rij2*rjk2)
            if (abs(cosfi2) .gt. cosfi_max) cycle
            cosfi3 = dot_product(rkl,rjk)/sqrt(rkl2*rjk2)
            if (abs(cosfi3) .gt. cosfi_max) cycle

            gij = fk_swart(1.0_wp,rij0,rij2) &
                  +0.5_wp*kd*fk_vdw(5.0_wp,d0ij,rij2)
            gjk = fk_swart(1.0_wp,rjk0,rjk2) &
                  +0.5_wp*kd*fk_vdw(5.0_wp,d0jk,rjk2)
            gkl = fk_swart(1.0_wp,rkl0,rkl2) &
                  +0.5_wp*kd*fk_vdw(5.0_wp,d0kl,rkl2)

            tij = kt*gij*gjk*gkl

            !tij = max(tij,10*min_fk)

            call trsn2(txyz,tau,c)
            si = c(:,1)
            sj = c(:,2)
            sk = c(:,3)
            sl = c(:,4)

            ! off diagonal block
            do ic = 1,3
              do jc = 1,3
                hess(ind(ic,i,jc,j)) = hess(ind(ic,i,jc,j))+tij*si(ic)*sj(jc)
                hess(ind(ic,i,jc,k)) = hess(ind(ic,i,jc,k))+tij*si(ic)*sk(jc)
                hess(ind(ic,i,jc,l)) = hess(ind(ic,i,jc,l))+tij*si(ic)*sl(jc)
                hess(ind(ic,j,jc,k)) = hess(ind(ic,j,jc,k))+tij*sj(ic)*sk(jc)
                hess(ind(ic,j,jc,l)) = hess(ind(ic,j,jc,l))+tij*sj(ic)*sl(jc)
                hess(ind(ic,k,jc,l)) = hess(ind(ic,k,jc,l))+tij*sk(ic)*sl(jc)
              end do
            end do

            ! diagonal block
            do ic = 1,3
              do jc = 1,ic
                hess(ind(ic,i,jc,i)) = hess(ind(ic,i,jc,i))+tij*si(ic)*si(jc)
                hess(ind(ic,j,jc,j)) = hess(ind(ic,j,jc,j))+tij*sj(ic)*sj(jc)
                hess(ind(ic,k,jc,k)) = hess(ind(ic,k,jc,k))+tij*sk(ic)*sk(jc)
                hess(ind(ic,l,jc,l)) = hess(ind(ic,l,jc,l))+tij*sl(ic)*sl(jc)
              end do
            end do

          end do torsion_lAt
        end do torsion_iAt
      end do torsion_kAt
    end do torsion_jAt

  end subroutine mh_swart_torsion

  pure subroutine mh_swart_outofp(n,at,xyz,hess,ko,kd,rcov,rvdw,lcutoff)
    implicit none

    integer,intent(in)    :: n
    integer,intent(in)    :: at(n)
    real(wp),intent(in)    :: xyz(3,n)
    real(wp),intent(inout) :: hess((3*n)*(3*n+1)/2)
    real(wp),intent(in)    :: ko
    real(wp),intent(in)    :: kd
    real(wp),intent(in)    :: rcov(:)
    real(wp),intent(in)    :: rvdw(:)
    logical,intent(in)    :: lcutoff(n,n)

    integer  :: i,ir,j,jr,k,kr,l,lr,ic,jc
    real(wp) :: txyz(3,4),c(3,4)
    real(wp) :: rij(3),rij0,d0ij,rij2,gij
    real(wp) :: rik(3),rik0,d0ik,rik2,gik
    real(wp) :: ril(3),ril0,d0il,ril2,gil
    real(wp) :: cosfi2,cosfi3,cosfi4
    real(wp) :: beta,tij,tau
    real(wp) :: si(3),sj(3),sk(3),sl(3)

!! ------------------------------------------------------------------------
!  Hessian for out-of-plane
!! ------------------------------------------------------------------------
    outofplane_iAt: do i = 1,n
      txyz(:,4) = xyz(:,i)
      outofplane_jAt: do j = 1,n
        if (j .eq. i) cycle outofplane_jAt
        if (lcutoff(j,i)) cycle outofplane_jAt
        txyz(:,1) = xyz(:,j)
        outofplane_kAt: do k = 1,n
          if (i .eq. k) cycle outofplane_kAt
          if (j .eq. k) cycle outofplane_kat
          if (lcutoff(k,i)) cycle outofplane_kAt
          if (lcutoff(k,j)) cycle outofplane_kAt
          txyz(:,2) = xyz(:,k)
          outofplane_lAt: do l = 1,n
            txyz(:,3) = xyz(:,l)
            if (l .eq. i) cycle outofplane_lAt
            if (l .eq. j) cycle outofplane_lAt
            if (l .eq. k) cycle outofplane_lAt
            if (lcutoff(l,i)) cycle outofplane_lAt
            if (lcutoff(l,k)) cycle outofplane_lAt
            if (lcutoff(l,j)) cycle outofplane_lAt

            rij = xyz(:,i)-xyz(:,j)
            rij0 = rcov(at(i))+rcov(at(j))
            d0ij = rvdw(at(i))+rvdw(at(j))

            rik = xyz(:,i)-xyz(:,k)
            rik0 = rcov(at(i))+rcov(at(k))
            d0ik = rvdw(at(i))+rvdw(at(k))

            ril = xyz(:,i)-xyz(:,l)
            ril0 = rcov(at(i))+rcov(at(l))
            d0il = rvdw(at(i))+rvdw(at(l))

            rij2 = sum(rij**2)
            rik2 = sum(rik**2)
            ril2 = sum(ril**2)

            cosfi2 = dot_product(rij,rik)/sqrt(rij2*rik2)
            if (abs(abs(cosfi2)-1.0_wp) .lt. 1.0e-1_wp) cycle
            cosfi3 = dot_product(rij,ril)/sqrt(rij2*ril2)
            if (abs(abs(cosfi3)-1.0_wp) .lt. 1.0e-1_wp) cycle
            cosfi4 = dot_product(rik,ril)/sqrt(rik2*ril2)
            if (abs(abs(cosfi4)-1.0_wp) .lt. 1.0e-1_wp) cycle

            gij = fk_swart(1.0_wp,rij0,rij2) &
                  +0.5_wp*kd*fk_vdw(5.0_wp,d0ij,rij2)
            gik = fk_swart(1.0_wp,rik0,rik2) &
                  +0.5_wp*kd*fk_vdw(5.0_wp,d0ik,rik2)
            gil = fk_swart(1.0_wp,ril0,ril2) &
                  +0.5_wp*kd*fk_vdw(5.0_wp,d0il,ril2)

            tij = ko*gij*gik*gil

            !tij = max(tij,10*min_fk)

            call outofp2(xyz,tau,c)
            If (abs(tau) .gt. 45.0d0*(pi/180.d0)) cycle

            si = c(:,4)
            sj = c(:,1)
            sk = c(:,2)
            sl = c(:,3)

            ! off diagonal block
            do ic = 1,3
              do jc = 1,3
                hess(ind(ic,i,jc,j)) = hess(ind(ic,i,jc,j))+tij*si(ic)*sj(jc)
                hess(ind(ic,i,jc,k)) = hess(ind(ic,i,jc,k))+tij*si(ic)*sk(jc)
                hess(ind(ic,i,jc,l)) = hess(ind(ic,i,jc,l))+tij*si(ic)*sl(jc)
                hess(ind(ic,j,jc,k)) = hess(ind(ic,j,jc,k))+tij*sj(ic)*sk(jc)
                hess(ind(ic,j,jc,l)) = hess(ind(ic,j,jc,l))+tij*sj(ic)*sl(jc)
                hess(ind(ic,k,jc,l)) = hess(ind(ic,k,jc,l))+tij*sk(ic)*sl(jc)
              end do
            end do

            ! diagonal block
            do ic = 1,3
              do jc = 1,ic
                hess(ind(ic,i,jc,i)) = hess(ind(ic,i,jc,i))+tij*si(ic)*si(jc)
                hess(ind(ic,j,jc,j)) = hess(ind(ic,j,jc,j))+tij*sj(ic)*sj(jc)
                hess(ind(ic,k,jc,k)) = hess(ind(ic,k,jc,k))+tij*sk(ic)*sk(jc)
                hess(ind(ic,l,jc,l)) = hess(ind(ic,l,jc,l))+tij*sl(ic)*sl(jc)
              end do
            end do

          end do outofplane_lAt
        end do outofplane_kAt
      end do outofplane_jAt
    end do outofplane_iAt

  end subroutine mh_swart_outofp

!========================================================================================!
!########################################################################################!
!========================================================================================!

  subroutine mh_lindh(xyz,n,hess,at,modh)
!**************************************************************************
!*  Lindh's Model Hessian updated around 2007
!* ------------------------------------------------------------------------
!*  R. Lindh, personal communication.
!*
!*  gij = exp[αij(R²ref - R²ij)]
!*  dij = exp[-4·(Rvdw - Rij)²]
!*  kij   = rkr·gij + rkd·dij
!*  kijk  = rkf·(gij+½·rkd/rkr·dij)·(gjk+½·rkd/rkr·djk)
!*  kijkl = rkt·(gij+½·rkd/rkr·dij)·(gjk+½·rkd/rkr·djk)·(gkl+½·rkd/rkr·dkl)
!*
!*  parameters tweaked by R. Lindh in 2007:
!*  rkr = 0.45, rkf = 0.10, rkt = 0.0025, rko = 0.16, rkd = 0.05
!*
!*  the reference distances are divided by rows in the PSE:
!*  rAv:        1        2        3       aAv:        1        2        3
!*    1    1.3500   2.1000   2.5300         1    1.0000   0.3949   0.3949
!*    2    2.1000   2.8700   3.8000         2    0.3949   0.2800   0.1200
!*    3    2.5300   3.8000   4.5000         3    0.3949   0.1200   0.0600
!*
!*  dAv:        1        2        3
!*    1    0.0000   3.6000   3.6000
!*    2    3.6000   5.3000   5.3000
!*    3    3.6000   5.3000   5.3000
!*
!**************************************************************************
    implicit none

    integer,intent(in)  :: n
    real(wp),intent(in)  :: xyz(3,n)
    real(wp),intent(out) :: hess((3*n)*(3*n+1)/2)
    integer,intent(in)  :: at(n)
    type(mhparam),intent(in) :: modh

    real(wp),parameter :: rAv(3,3) = reshape( &
                          (/1.3500_wp,2.1000_wp,2.5300_wp, &
                            2.1000_wp,2.8700_wp,3.8000_wp, &
                            2.5300_wp,3.8000_wp,4.5000_wp/),shape(rAv))
    real(wp),parameter :: aAv(3,3) = reshape( &
                          (/1.0000_wp,0.3949_wp,0.3949_wp, &
                            0.3949_wp,0.2800_wp,0.1200_wp, &
                            0.3949_wp,0.1200_wp,0.0600_wp/),shape(aAv))
    real(wp),parameter :: dAv(3,3) = reshape( &
                          (/0.0000_wp,3.6000_wp,3.6000_wp, &
                            3.6000_wp,5.3000_wp,5.3000_wp, &
                            3.6000_wp,5.3000_wp,5.3000_wp/),shape(aAv))

    integer  :: n3
    real(wp) :: kd
    logical,allocatable :: lcutoff(:,:)
    !type(chrg_parameter) :: chrgeq

    allocate (lcutoff(n,n),source=.false.)

    n3 = 3*n
    hess = 0.0d0

!  the dispersion force constant is used relative to the stretch force constant
    kd = modh%kd/modh%kr

    call mh_lindh_stretch(n,at,xyz,hess,modh%kr,kd,modh%s6,aav,rav,dav,lcutoff,modh%rcut)
    if (modh%kf .ne. 0.0_wp) &
      call mh_lindh_bend(n,at,xyz,hess,modh%kf,kd,aav,rav,dav,lcutoff)
    if (modh%kt .ne. 0.0_wp) &
      call mh_lindh_torsion(n,at,xyz,hess,modh%kt,kd,aav,rav,dav,lcutoff)
    if (modh%ko .ne. 0.0_wp) &
      call mh_lindh_outofp(n,at,xyz,hess,modh%ko,0.0_wp,aav,rav,dav,lcutoff)
    if (modh%kq .ne. 0.0_wp) then
      !call new_charge_model_2019(chrgeq,n,at)
      call mh_eeq(n,at,xyz,0.0_wp,modh%kq,hess)
    end if

  end subroutine mh_lindh

  subroutine mh_lindh_d2(xyz,n,hess,at,modh)
!**************************************************************************
!*  Lindh's Model Hessian augmented with D2
!* ------------------------------------------------------------------------
!*  Implemented after:
!*  Lindh, R., Bernhardsson, A., Karlström, G., & Malmqvist, P.-Å. (1995).
!*  On the use of a Hessian model function in molecular geometry optimizations.
!*  Chem. Phys. Lett., 241(4), 423–428. doi:10.1016/0009-2614(95)00646-l
!*
!*  gij = exp[αij(R²ref - R²ij)]
!*  kij   = rkr·gij
!*  kijk  = rkf·gij·gjk
!*  kijkl = rkt·gij·gjk·gkl
!*
!*  Originally Lindh proposed (we tweaked those a little bit):
!*  rkr = 0.45, rkf = 0.15, rkt = 0.005
!*
!*  the reference distances are divided by rows in the PSE:
!*  rAv:        1        2        3       aAv:        1        2        3
!*    1    1.3500   2.1000   2.5300         1    1.0000   0.3949   0.3949
!*    2    2.1000   2.8700   3.4000         2    0.3949   0.2800   0.2800
!*    3    2.5300   3.4000   3.4000         3    0.3949   0.2800   0.2800
!*
!*  This Hessian is additionally augmented with D2, please note that D2
!*  is not implemented in atomic units and requires some magical conversion
!*  factor somewhere hidden in the implementation below.
!*************************************************************************
    implicit none

    integer,intent(in)  :: n
    real(wp),intent(in)  :: xyz(3,n)
    real(wp),intent(out) :: hess((3*n)*(3*n+1)/2)
    integer,intent(in)  :: at(n)
    type(mhparam),intent(in) :: modh

    real(wp),parameter :: rAv(3,3) = reshape( &
                          (/1.3500_wp,2.1000_wp,2.5300_wp, &
                            2.1000_wp,2.8700_wp,3.4000_wp, &
                            2.5300_wp,3.4000_wp,3.4000_wp/),shape(rAv))
    real(wp),parameter :: aAv(3,3) = reshape( &
                          (/1.0000_wp,0.3949_wp,0.3949_wp, &
                            0.3949_wp,0.2800_wp,0.2800_wp, &
                            0.3949_wp,0.2800_wp,0.2800_wp/),shape(aAv))
    real(wp),parameter :: dAv(3,3) = reshape( &
                          (/0.0000_wp,0.0000_wp,0.0000_wp, &
                            0.0000_wp,0.0000_wp,0.0000_wp, &
                            0.0000_wp,0.0000_wp,0.0000_wp/),shape(aAv))

    integer  :: n3
    real(wp) :: kd
    logical,allocatable :: lcutoff(:,:)

    allocate (lcutoff(n,n),source=.false.)

    n3 = 3*n
    hess = 0.0d0

!  the dispersion force constant is used relative to the stretch force constant
    kd = modh%kd/modh%kr

    call mh_lindh_stretch(n,at,xyz,hess,modh%kr,kd,modh%s6,aav,rav,dav,lcutoff,modh%rcut)
    if (modh%kf .ne. 0.0_wp) &
      call mh_lindh_bend(n,at,xyz,hess,modh%kf,kd,aav,rav,dav,lcutoff)
    if (modh%kt .ne. 0.0_wp) &
      call mh_lindh_torsion(n,at,xyz,hess,modh%kt,kd,aav,rav,dav,lcutoff)
    if (modh%ko .ne. 0.0_wp) &
      call mh_lindh_outofp(n,at,xyz,hess,modh%ko,kd,aav,rav,dav,lcutoff)
    if (modh%kq .ne. 0.0_wp) then
      call mh_eeq(n,at,xyz,0.0_wp,modh%kq,hess)
    end if
  end subroutine mh_lindh_d2

  pure subroutine mh_lindh_stretch(n,at,xyz,hess,kr,kd,s6,aav,rav,dav,lcutoff,rcut)
    implicit none

    integer,intent(in)    :: n
    integer,intent(in)    :: at(n)
    real(wp),intent(in)    :: xyz(3,n)
    real(wp),intent(inout) :: hess((3*n)*(3*n+1)/2)
    real(wp),intent(in)    :: kr
    real(wp),intent(in)    :: kd
    real(wp),intent(in)    :: s6
    real(wp),intent(in)    :: aav(3,3)
    real(wp),intent(in)    :: rav(3,3)
    real(wp),intent(in)    :: dav(3,3)
    logical,intent(out)   :: lcutoff(n,n)
    real(wp),intent(in)    :: rcut

    integer  :: i,ir,j,jr
    real(wp) :: xij,yij,zij,rij2,r0,d0
    real(wp) :: alpha,gmm
    real(wp) :: c6i,c6j,c6ij,rv
    real(wp) :: hxx,hxy,hxz,hyy,hyz,hzz
    real(wp) :: vdw(3,3)

!! ------------------------------------------------------------------------
!  Hessian for stretch
!! ------------------------------------------------------------------------
    stretch_iAt: do i = 1,n
      ir = itabrow(at(i))

      stretch_jAt: do j = 1,i-1
        jr = itabrow(at(j))

        ! save for later
        lcutoff(i,j) = rcutoff(xyz,i,j,rcut)
        lcutoff(j,i) = lcutoff(i,j)

        xij = xyz(1,i)-xyz(1,j)
        yij = xyz(2,i)-xyz(2,j)
        zij = xyz(3,i)-xyz(3,j)
        rij2 = xij**2+yij**2+zij**2
        r0 = rav(ir,jr)
        d0 = dav(ir,jr)
        alpha = aav(ir,jr)

        !cccccc vdwx ccccccccccccccccccccccccccccccccc
        c6i = c6(at(i))
        c6j = c6(at(j))
        c6ij = sqrt(c6i*c6j)
        rv = (vander(at(i))+vander(at(j)))*aatoau

        call getvdwxx(xij,yij,zij,c6ij,s6,rv,vdw(1,1))
        call getvdwxy(xij,yij,zij,c6ij,s6,rv,vdw(1,2))
        call getvdwxy(xij,zij,yij,c6ij,s6,rv,vdw(1,3))
        call getvdwxx(yij,xij,zij,c6ij,s6,rv,vdw(2,2))
        call getvdwxy(yij,zij,xij,c6ij,s6,rv,vdw(2,3))
        call getvdwxx(zij,xij,yij,c6ij,s6,rv,vdw(3,3))
        !cccccc ende vdwx ccccccccccccccccccccccccccccccc

        gmm = kr*fk_lindh(alpha,r0,rij2) &
              +kr*kd*fk_vdw(4.0_wp,d0,rij2)

        !gmm = max(gmm,min_fk)

        hxx = gmm*xij*xij/rij2-vdw(1,1)
        hxy = gmm*xij*yij/rij2-vdw(1,2)
        hxz = gmm*xij*zij/rij2-vdw(1,3)
        hyy = gmm*yij*yij/rij2-vdw(2,2)
        hyz = gmm*yij*zij/rij2-vdw(2,3)
        hzz = gmm*zij*zij/rij2-vdw(3,3)

        ! save diagonal elements for atom i
        hess(ind(1,i,1,i)) = hess(ind(1,i,1,i))+hxx
        hess(ind(2,i,1,i)) = hess(ind(2,i,1,i))+hxy
        hess(ind(2,i,2,i)) = hess(ind(2,i,2,i))+hyy
        hess(ind(3,i,1,i)) = hess(ind(3,i,1,i))+hxz
        hess(ind(3,i,2,i)) = hess(ind(3,i,2,i))+hyz
        hess(ind(3,i,3,i)) = hess(ind(3,i,3,i))+hzz
        ! save elements between atom i and atom j
        hess(ind(1,i,1,j)) = hess(ind(1,i,1,j))-hxx
        hess(ind(1,i,2,j)) = hess(ind(1,i,2,j))-hxy
        hess(ind(1,i,3,j)) = hess(ind(1,i,3,j))-hxz
        hess(ind(2,i,1,j)) = hess(ind(2,i,1,j))-hxy
        hess(ind(2,i,2,j)) = hess(ind(2,i,2,j))-hyy
        hess(ind(2,i,3,j)) = hess(ind(2,i,3,j))-hyz
        hess(ind(3,i,1,j)) = hess(ind(3,i,1,j))-hxz
        hess(ind(3,i,2,j)) = hess(ind(3,i,2,j))-hyz
        hess(ind(3,i,3,j)) = hess(ind(3,i,3,j))-hzz
        ! save diagonal elements for atom j
        hess(ind(1,j,1,j)) = hess(ind(1,j,1,j))+hxx
        hess(ind(2,j,1,j)) = hess(ind(2,j,1,j))+hxy
        hess(ind(2,j,2,j)) = hess(ind(2,j,2,j))+hyy
        hess(ind(3,j,1,j)) = hess(ind(3,j,1,j))+hxz
        hess(ind(3,j,2,j)) = hess(ind(3,j,2,j))+hyz
        hess(ind(3,j,3,j)) = hess(ind(3,j,3,j))+hzz

      end do stretch_jAt
    end do stretch_iAt

  end subroutine mh_lindh_stretch

  pure subroutine mh_lindh_bend(n,at,xyz,hess,kf,kd,aav,rav,dav,lcutoff)
    implicit none

    integer,intent(in)    :: n
    integer,intent(in)    :: at(n)
    real(wp),intent(in)    :: xyz(3,n)
    real(wp),intent(inout) :: hess((3*n)*(3*n+1)/2)
    real(wp),intent(in)    :: kf
    real(wp),intent(in)    :: kd
    real(wp),intent(in)    :: aav(3,3)
    real(wp),intent(in)    :: rav(3,3)
    real(wp),intent(in)    :: dav(3,3)
    logical,intent(in)    :: lcutoff(n,n)

    integer  :: i,ir,j,jr,m,mr,ic,jc,ii
    real(wp),parameter :: rzero = 1.0e-10_wp
    real(wp) :: xij,yij,zij,rij2,rrij,r1
    real(wp) :: xmi,ymi,zmi,rmi2,rmi,r0mi,ami,d0mj,gmi
    real(wp) :: xmj,ymj,zmj,rmj2,rmj,r0mj,amj,d0mi,gmj
    real(wp) :: test,gij,rl2,rl,rmidotrmj
    real(wp) :: sinphi,cosphi,costhetax,costhetay,costhetaz
    real(wp) :: alpha
    real(wp) :: si(3),sj(3),sm(3),x(2),y(2),z(2)

!! ------------------------------------------------------------------------
!  Hessian for bending
!! ------------------------------------------------------------------------
    bend_mAt: do m = 1,n
      mr = itabrow(at(m))
      bend_iAt: do i = 1,n
        if (i .eq. m) cycle bend_iAt
        ir = itabrow(at(i))
        if (lcutoff(i,m)) cycle bend_iAt

        xmi = (xyz(1,i)-xyz(1,m))
        ymi = (xyz(2,i)-xyz(2,m))
        zmi = (xyz(3,i)-xyz(3,m))
        rmi2 = xmi**2+ymi**2+zmi**2
        rmi = sqrt(rmi2)
        r0mi = rav(mr,ir)
        d0mi = dav(mr,ir)
        ami = aav(mr,ir)

        bend_jAt: do j = 1,i-1
          if (j .eq. m) cycle bend_jAt
          jr = itabrow(at(j))
          if (lcutoff(j,i)) cycle bend_jAt
          if (lcutoff(j,m)) cycle bend_jAt

          xmj = (xyz(1,j)-xyz(1,m))
          ymj = (xyz(2,j)-xyz(2,m))
          zmj = (xyz(3,j)-xyz(3,m))
          rmj2 = xmj**2+ymj**2+zmj**2
          rmj = sqrt(rmj2)
          r0mj = rav(mr,jr)
          d0mj = dav(mr,jr)
          amj = aav(mr,jr)

          ! test if zero angle
          test = xmi*xmj+ymi*ymj+zmi*zmj
          test = test/(rmi*rmj)
          if (abs(test-1.0_wp) .lt. 1.0e-12_wp) cycle bend_jAt

          xij = (xyz(1,j)-xyz(1,i))
          yij = (xyz(2,j)-xyz(2,i))
          zij = (xyz(3,j)-xyz(3,i))
          rij2 = xij**2+yij**2+zij**2
          rrij = sqrt(rij2)

          gmi = fk_lindh(ami,r0mi,rmi2) &
                +0.5_wp*kd*fk_vdw(4.0_wp,d0mi,rmi2)
          gmj = fk_lindh(amj,r0mj,rmj2) &
                +0.5_wp*kd*fk_vdw(4.0_wp,d0mj,rmj2)

          gij = kf*gmi*gmj

          rl2 = (ymi*zmj-zmi*ymj)**2+(zmi*xmj-xmi*zmj)**2+(xmi*ymj-ymi*xmj)**2

          if (rl2 .lt. 1.e-14_wp) then
            rl = 0.0_wp
          else
            rl = sqrt(rl2)
          end if

          !gij = max(gij,min_fk)

          if ((rmj .gt. rzero).and.(rmi .gt. rzero).and.(rrij .gt. rzero)) then
            sinphi = rl/(rmj*rmi)
            rmidotrmj = xmi*xmj+ymi*ymj+zmi*zmj
            cosphi = rmidotrmj/(rmj*rmi)
            ! none linear case
            if (sinphi .gt. rzero) then
              si(1) = (xmi/rmi*cosphi-xmj/rmj)/(rmi*sinphi)
              si(2) = (ymi/rmi*cosphi-ymj/rmj)/(rmi*sinphi)
              si(3) = (zmi/rmi*cosphi-zmj/rmj)/(rmi*sinphi)
              sj(1) = (cosphi*xmj/rmj-xmi/rmi)/(rmj*sinphi)
              sj(2) = (cosphi*ymj/rmj-ymi/rmi)/(rmj*sinphi)
              sj(3) = (cosphi*zmj/rmj-zmi/rmi)/(rmj*sinphi)
              sm(1) = -si(1)-sj(1)
              sm(2) = -si(2)-sj(2)
              sm(3) = -si(3)-sj(3)
              do ic = 1,3
                do jc = 1,3
                  if (m .gt. i) then
                    hess(ind(ic,m,jc,i)) = hess(ind(ic,m,jc,i)) &
                                           +gij*sm(ic)*si(jc)
                  else
                    hess(ind(ic,i,jc,m)) = hess(ind(ic,i,jc,m)) &
                                           +gij*si(ic)*sm(jc)
                  end if
                  if (m .gt. j) then
                    hess(ind(ic,m,jc,j)) = hess(ind(ic,m,jc,j)) &
                                           +gij*sm(ic)*sj(jc)
                  else
                    hess(ind(ic,j,jc,m)) = hess(ind(ic,j,jc,m)) &
                                           +gij*sj(ic)*sm(jc)
                  end if
                  if (i .gt. j) then
                    hess(ind(ic,i,jc,j)) = hess(ind(ic,i,jc,j)) &
                                           +gij*si(ic)*sj(jc)
                  else
                    hess(ind(ic,j,jc,i)) = hess(ind(ic,j,jc,i)) &
                                           +gij*sj(ic)*si(jc)
                  end if
                end do
              end do
              do ic = 1,3
                do jc = 1,ic
                  hess(ind(ic,i,jc,i)) = hess(ind(ic,i,jc,i))+gij*si(ic)*si(jc)
                  hess(ind(ic,m,jc,m)) = hess(ind(ic,m,jc,m))+gij*sm(ic)*sm(jc)
                  hess(ind(ic,j,jc,j)) = hess(ind(ic,j,jc,j))+gij*sj(ic)*sj(jc)
                end do
              end do
            else
              ! linear case
              if ((abs(ymi) .gt. rzero).or.(abs(xmi) .gt. rzero)) then
                x(1) = -ymi
                y(1) = xmi
                z(1) = 0.0_wp
                x(2) = -xmi*zmi
                y(2) = -ymi*zmi
                z(2) = xmi*xmi+ymi*ymi
              else
                x(1) = 1.0_wp
                y(1) = 0.0_wp
                z(1) = 0.0_wp
                x(2) = 0.0_wp
                y(2) = 1.0_wp
                z(2) = 0.0_wp
              end if
              do ii = 1,2
                r1 = sqrt(x(ii)**2+y(ii)**2+z(ii)**2)
                costhetax = x(ii)/r1
                costhetay = y(ii)/r1
                costhetaz = z(ii)/r1
                si(1) = -costhetax/rmi
                si(2) = -costhetay/rmi
                si(3) = -costhetaz/rmi
                sj(1) = -costhetax/rmj
                sj(2) = -costhetay/rmj
                sj(3) = -costhetaz/rmj
                sm(1) = -(si(1)+sj(1))
                sm(2) = -(si(2)+sj(2))
                sm(3) = -(si(3)+sj(3))
                !
                do ic = 1,3
                  do jc = 1,3
                    if (m .gt. i) then
                      hess(ind(ic,m,jc,i)) = hess(ind(ic,m,jc,i)) &
                                             +gij*sm(ic)*si(jc)
                    else
                      hess(ind(ic,i,jc,m)) = hess(ind(ic,i,jc,m)) &
                                             +gij*si(ic)*sm(jc)
                    end if
                    if (m .gt. j) then
                      hess(ind(ic,m,jc,j)) = hess(ind(ic,m,jc,j)) &
                                             +gij*sm(ic)*sj(jc)
                    else
                      hess(ind(ic,j,jc,m)) = hess(ind(ic,j,jc,m)) &
                                             +gij*sj(ic)*sm(jc)
                    end if
                    if (i .gt. j) then
                      hess(ind(ic,i,jc,j)) = hess(ind(ic,i,jc,j)) &
                                             +gij*si(ic)*sj(jc)
                    else
                      hess(ind(ic,j,jc,i)) = hess(ind(ic,j,jc,i)) &
                                             +gij*sj(ic)*si(jc)
                    end if
                  end do
                end do
                do ic = 1,3
                  do jc = 1,ic
                    hess(ind(ic,i,jc,i)) = hess(ind(ic,i,jc,i)) &
                                           +gij*si(ic)*si(jc)
                    hess(ind(ic,m,jc,m)) = hess(ind(ic,m,jc,m)) &
                                           +gij*sm(ic)*sm(jc)
                    hess(ind(ic,j,jc,j)) = hess(ind(ic,j,jc,j)) &
                                           +gij*sj(ic)*sj(jc)
                  end do
                end do
              end do

            end if
          end if

        end do bend_jAt
      end do bend_iAt
    end do bend_mAt

  end subroutine mh_lindh_bend

  subroutine mh_lindh_torsion(n,at,xyz,hess,kt,kd,aav,rav,dav,lcutoff)
    implicit none

    integer,intent(in)    :: n
    integer,intent(in)    :: at(n)
    real(wp),intent(in)    :: xyz(3,n)
    real(wp),intent(inout) :: hess((3*n)*(3*n+1)/2)
    real(wp),intent(in)    :: kt
    real(wp),intent(in)    :: kd
    real(wp),intent(in)    :: aav(3,3)
    real(wp),intent(in)    :: rav(3,3)
    real(wp),intent(in)    :: dav(3,3)
    logical,intent(in)    :: lcutoff(n,n)

    integer  :: i,ir,j,jr,k,kr,l,lr,ic,jc,ij,kl
!  allow only angles in the range of 35-145
    real(wp),parameter :: a35 = (35.0d0/180.d0)*pi
    real(wp),parameter :: cosfi_max = cos(a35)
    real(wp) :: txyz(3,4),c(3,4)
    real(wp) :: rij(3),rij0,aij,rij2,d0ij,gij
    real(wp) :: rjk(3),rjk0,ajk,rjk2,d0jk,gjk
    real(wp) :: rkl(3),rkl0,akl,rkl2,d0kl,gkl
    real(wp) :: cosfi2,cosfi3,cosfi4
    real(wp) :: beta,tij,tau,dum(3,4,3,4)
    real(wp) :: si(3),sj(3),sk(3),sl(3)

!! ------------------------------------------------------------------------
!  Hessian for torsion
!! ------------------------------------------------------------------------
    torsion_jAt: do j = 1,n
      jr = itabrow(at(j))
      txyz(:,2) = xyz(:,j)
      torsion_kAt: do k = 1,n
        if (k .eq. j) cycle torsion_kAt
        kr = itabrow(at(k))
        if (lcutoff(k,j)) cycle torsion_kAt
        txyz(:,3) = xyz(:,k)
        torsion_iAt: do i = 1,n
          ij = n*(j-1)+i
          if (i .eq. j) cycle torsion_iAt
          if (i .eq. k) cycle torsion_iAt
          ir = itabrow(at(i))
          if (lcutoff(i,k)) cycle torsion_iAt
          if (lcutoff(i,j)) cycle torsion_iAt

          txyz(:,1) = xyz(:,i)
          torsion_lAt: do l = 1,n
            kl = n*(k-1)+l
            if (ij .le. kl) cycle torsion_lAt
            if (l .eq. i) cycle torsion_lAt
            if (l .eq. j) cycle torsion_lAt
            if (l .eq. k) cycle torsion_lAt
            lr = itabrow(at(l))
!
            if (lcutoff(l,i)) cycle torsion_lAt
            if (lcutoff(l,k)) cycle torsion_lAt
            if (lcutoff(l,j)) cycle torsion_lAt

            txyz(:,4) = xyz(:,l)

            rij = xyz(:,i)-xyz(:,j)
            d0ij = dav(ir,jr)
            rij0 = rav(ir,jr)
            aij = aav(ir,jr)

            rjk = xyz(:,j)-xyz(:,k)
            d0jk = dav(jr,kr)
            rjk0 = rav(jr,kr)
            ajk = aav(jr,kr)

            rkl = xyz(:,k)-xyz(:,l)
            d0kl = dav(kr,lr)
            rkl0 = rav(kr,lr)
            akl = aav(kr,lr)

            rij2 = sum(rij**2)
            rjk2 = sum(rjk**2)
            rkl2 = sum(rjk**2)

            cosfi2 = dot_product(rij,rjk)/sqrt(rij2*rjk2)
            if (abs(cosfi2) .gt. cosfi_max) cycle
            cosfi3 = dot_product(rkl,rjk)/sqrt(rkl2*rjk2)
            if (abs(cosfi3) .gt. cosfi_max) cycle

            gij = fk_lindh(aij,rij0,rij2) &
                  +0.5_wp*kd*fk_vdw(4.0_wp,d0ij,rij2)
            gjk = fk_lindh(ajk,rjk0,rjk2) &
                  +0.5_wp*kd*fk_vdw(4.0_wp,d0jk,rjk2)
            gkl = fk_lindh(akl,rkl0,rkl2) &
                  +0.5_wp*kd*fk_vdw(4.0_wp,d0kl,rkl2)

            tij = kt*gij*gjk*gkl

            !tij = max(tij,10*min_fk)

            !call trsn2(txyz,tau,c)
            Call Trsn(txyz,4,Tau,C,.False.,.False.,'        ', &
   &                  Dum,.False.)
            si = c(:,1)
            sj = c(:,2)
            sk = c(:,3)
            sl = c(:,4)

            ! off diagonal block
            do ic = 1,3
              do jc = 1,3
                hess(ind(ic,i,jc,j)) = hess(ind(ic,i,jc,j))+tij*si(ic)*sj(jc)
                hess(ind(ic,i,jc,k)) = hess(ind(ic,i,jc,k))+tij*si(ic)*sk(jc)
                hess(ind(ic,i,jc,l)) = hess(ind(ic,i,jc,l))+tij*si(ic)*sl(jc)
                hess(ind(ic,j,jc,k)) = hess(ind(ic,j,jc,k))+tij*sj(ic)*sk(jc)
                hess(ind(ic,j,jc,l)) = hess(ind(ic,j,jc,l))+tij*sj(ic)*sl(jc)
                hess(ind(ic,k,jc,l)) = hess(ind(ic,k,jc,l))+tij*sk(ic)*sl(jc)
              end do
            end do

            ! diagonal block
            do ic = 1,3
              do jc = 1,ic
                hess(ind(ic,i,jc,i)) = hess(ind(ic,i,jc,i))+tij*si(ic)*si(jc)
                hess(ind(ic,j,jc,j)) = hess(ind(ic,j,jc,j))+tij*sj(ic)*sj(jc)
                hess(ind(ic,k,jc,k)) = hess(ind(ic,k,jc,k))+tij*sk(ic)*sk(jc)
                hess(ind(ic,l,jc,l)) = hess(ind(ic,l,jc,l))+tij*sl(ic)*sl(jc)
              end do
            end do

          end do torsion_lAt
        end do torsion_iAt
      end do torsion_kAt
    end do torsion_jAt

  end subroutine mh_lindh_torsion

  pure subroutine mh_lindh_outofp(n,at,xyz,hess,ko,kd,aav,rav,dav,lcutoff)
    implicit none

    integer,intent(in)    :: n
    integer,intent(in)    :: at(n)
    real(wp),intent(in)    :: xyz(3,n)
    real(wp),intent(inout) :: hess((3*n)*(3*n+1)/2)
    real(wp),intent(in)    :: ko
    real(wp),intent(in)    :: kd
    real(wp),intent(in)    :: aav(3,3)
    real(wp),intent(in)    :: rav(3,3)
    real(wp),intent(in)    :: dav(3,3)
    logical,intent(in)    :: lcutoff(n,n)

    integer  :: i,ir,j,jr,k,kr,l,lr,ic,jc
    real(wp) :: txyz(3,4),c(3,4)
    real(wp) :: rij(3),rij0,aij,rij2,gij,d0ij
    real(wp) :: rik(3),rik0,aik,rik2,gik,d0ik
    real(wp) :: ril(3),ril0,ail,ril2,gil,d0il
    real(wp) :: cosfi2,cosfi3,cosfi4
    real(wp) :: beta,tij,tau
    real(wp) :: si(3),sj(3),sk(3),sl(3)

!! ------------------------------------------------------------------------
!  Hessian for out-of-plane
!! ------------------------------------------------------------------------
    outofplane_iAt: do i = 1,n
      ir = itabrow(at(i))
      txyz(:,4) = xyz(:,i)
      outofplane_jAt: do j = 1,n
        if (j .eq. i) cycle outofplane_jAt
        if (lcutoff(j,i)) cycle outofplane_jAt
        jr = itabrow(at(j))
        txyz(:,1) = xyz(:,j)
        outofplane_kAt: do k = 1,n
          if (i .eq. k) cycle outofplane_kAt
          if (j .eq. k) cycle outofplane_kat
          if (lcutoff(k,i)) cycle outofplane_kAt
          if (lcutoff(k,j)) cycle outofplane_kAt
          kr = itabrow(at(k))
          txyz(:,2) = xyz(:,k)
          outofplane_lAt: do l = 1,n
            lr = itabrow(at(l))
            txyz(:,3) = xyz(:,l)
            if (l .eq. i) cycle outofplane_lAt
            if (l .eq. j) cycle outofplane_lAt
            if (l .eq. k) cycle outofplane_lAt
            if (lcutoff(l,i)) cycle outofplane_lAt
            if (lcutoff(l,k)) cycle outofplane_lAt
            if (lcutoff(l,j)) cycle outofplane_lAt

            rij = xyz(:,i)-xyz(:,j)
            d0ij = dav(ir,jr)
            rij0 = rav(ir,jr)
            aij = aav(ir,jr)

            rik = xyz(:,i)-xyz(:,k)
            d0ik = dav(ir,kr)
            rik0 = rav(ir,kr)
            aik = aav(ir,kr)

            ril = xyz(:,i)-xyz(:,l)
            d0il = dav(ir,lr)
            ril0 = rav(ir,lr)
            ail = aav(ir,lr)

            rij2 = sum(rij**2)
            rik2 = sum(rik**2)
            ril2 = sum(ril**2)

            cosfi2 = dot_product(rij,rik)/sqrt(rij2*rik2)
            if (abs(abs(cosfi2)-1.0_wp) .lt. 1.0e-1_wp) cycle
            cosfi3 = dot_product(rij,ril)/sqrt(rij2*ril2)
            if (abs(abs(cosfi3)-1.0_wp) .lt. 1.0e-1_wp) cycle
            cosfi4 = dot_product(rik,ril)/sqrt(rik2*ril2)
            if (abs(abs(cosfi4)-1.0_wp) .lt. 1.0e-1_wp) cycle

            gij = fk_lindh(aij,rij0,rij2) &
                  +0.5_wp*kd*fk_vdw(4.0_wp,d0ij,rij2)
            gik = fk_lindh(aik,rik0,rik2) &
                  +0.5_wp*kd*fk_vdw(4.0_wp,d0ik,rik2)
            gil = fk_lindh(ail,ril0,ril2) &
                  +0.5_wp*kd*fk_vdw(4.0_wp,d0il,ril2)

            tij = ko*gij*gik*gil

            !tij = max(tij,10*min_fk)

            call outofp2(xyz,tau,c)
            If (abs(tau) .gt. 45.0d0*(pi/180.d0)) cycle

            si = c(:,4)
            sj = c(:,1)
            sk = c(:,2)
            sl = c(:,3)

            ! off diagonal block
            do ic = 1,3
              do jc = 1,3
                hess(ind(ic,i,jc,j)) = hess(ind(ic,i,jc,j))+tij*si(ic)*sj(jc)
                hess(ind(ic,i,jc,k)) = hess(ind(ic,i,jc,k))+tij*si(ic)*sk(jc)
                hess(ind(ic,i,jc,l)) = hess(ind(ic,i,jc,l))+tij*si(ic)*sl(jc)
                hess(ind(ic,j,jc,k)) = hess(ind(ic,j,jc,k))+tij*sj(ic)*sk(jc)
                hess(ind(ic,j,jc,l)) = hess(ind(ic,j,jc,l))+tij*sj(ic)*sl(jc)
                hess(ind(ic,k,jc,l)) = hess(ind(ic,k,jc,l))+tij*sk(ic)*sl(jc)
              end do
            end do

            ! diagonal block
            do ic = 1,3
              do jc = 1,ic
                hess(ind(ic,i,jc,i)) = hess(ind(ic,i,jc,i))+tij*si(ic)*si(jc)
                hess(ind(ic,j,jc,j)) = hess(ind(ic,j,jc,j))+tij*sj(ic)*sj(jc)
                hess(ind(ic,k,jc,k)) = hess(ind(ic,k,jc,k))+tij*sk(ic)*sk(jc)
                hess(ind(ic,l,jc,l)) = hess(ind(ic,l,jc,l))+tij*sl(ic)*sl(jc)
              end do
            end do

          end do outofplane_lAt
        end do outofplane_kAt
      end do outofplane_jAt
    end do outofplane_iAt

  end subroutine mh_lindh_outofp

!========================================================================================!
!########################################################################################!
!========================================================================================!

  pure function rcutoff(xyz,katom,latom,rcut)
    implicit none
    logical  :: rcutoff
    real(wp),intent(in) :: xyz(3,*)
    real(wp),intent(in) :: rcut
    real(wp) :: rkl(3),rkl2
    integer,intent(in) :: katom,latom
    rcutoff = .false.
    rkl = xyz(:,kAtom)-xyz(:,lAtom)
    rkl2 = sum(rkl**2)
    if (rkl2 .gt. rcut) rcutoff = .true.
  end function rcutoff

  pure elemental function itabrow(i)
    integer :: itabrow
    integer,intent(in) :: i

    itabrow = 0
    if (i .gt. 0.and.i .le. 2) then
      itabrow = 1
    else if (i .gt. 2.and.i .le. 10) then
      itabrow = 2
    else if (i .gt. 10.and.i .le. 18) then
      itabrow = 3
    else if (i .gt. 18.and.i .le. 36) then
      itabrow = 3
    else if (i .gt. 36.and.i .le. 54) then
      itabrow = 3
    else if (i .gt. 54.and.i .le. 86) then
      itabrow = 3
    else if (i .gt. 86) then
      itabrow = 3
    end if

    return
  end function itabrow

  pure subroutine getvdwxy(rx,ry,rz,c66,s6,r0,vdw)
    !cc Ableitung nach rx und ry
    implicit none
    real(wp),intent(in)  :: rx,ry,rz,c66,s6,r0
    real(wp),intent(out) :: vdw
    real(wp) :: t1,t2,t3,t4,t5,t6,t7,t11,t12,t16,t17,t25,t26,t35
    real(wp) :: t40,t41,t43,t44,t56,avdw

    !    write(*,*) 's6:', s6
    avdw = 20.0
    t1 = s6*C66
    t2 = rx**2
    t3 = ry**2
    t4 = rz**2
    t5 = t2+t3+t4
    t6 = t5**2
    t7 = t6**2
    t11 = sqrt(t5)
    t12 = 0.1D1/r0
    t16 = exp(-avdw*(t11*t12-0.1D1))
    t17 = 0.1D1+t16
    t25 = t17**2
    t26 = 0.1D1/t25
    t35 = 0.1D1/t7
    t40 = avdw**2
    t41 = r0**2
    t43 = t40/t41
    t44 = t16**2
    t56 = -0.48D2*t1/t7/t5/t17*rx*ry+0.13D2*t1/t11/&
       & t7*t26*rx*avdw*t12*ry*t16-0.2D1*t1*t35/t25/&
       &t17*t43*rx*t44*ry+t1*t35*t26*t43*rx*ry*t16
    vdw = t56
    return
  end subroutine getvdwxy

  pure subroutine getvdwxx(rx,ry,rz,c66,s6,r0,vdw)
    !cc Ableitung nach rx und rx
    Implicit none
    real(wp),intent(in)  :: rx,ry,rz,c66,s6,r0
    real(wp),intent(out) :: vdw
    real(wp) :: t1,t2,t3,t4,t5,t6,t7,t10,t11,t15,t16,t17,t24,t25,t29
    real(wp) :: t33,t41,t42,t44,t45,t62,avdw
    avdw = 20.0
    !      write(*,*) 's6:', s6
    t1 = s6*C66
    t2 = rx**2
    t3 = ry**2
    t4 = rz**2
    t5 = t2+t3+t4
    t6 = t5**2
    t7 = t6**2
    t10 = sqrt(t5)
    t11 = 0.1D1/r0
    t15 = exp(-avdw*(t10*t11-0.1D1))
    t16 = 0.1D1+t15
    t17 = 0.1D1/t16
    t24 = t16**2
    t25 = 0.1D1/t24
    t29 = t11*t15
    t33 = 0.1D1/t7
    t41 = avdw**2
    t42 = r0**2
    t44 = t41/t42
    t45 = t15**2
    t62 = -0.48D2*t1/t7/t5*t17*t2+0.13D2*t1/t10/t7*&
       & t25*t2*avdw*t29+0.6D1*t1*t33*t17-0.2D1*t1*t33&
       & /t24/t16*t44*t2*t45-t1/t10/t6/t5*t25*avdw*&
       &t29+t1*t33*t25*t44*t2*t15
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
    costau = ((brij(2,1)*brjk(3,2)-brij(3,1)*brjk(2,2))* &
              (brjk(2,1)*brkl(3,2)-brjk(3,1)*brkl(2,2))+ &
              (brij(3,1)*brjk(1,2)-brij(1,1)*brjk(3,2))* &
              (brjk(3,1)*brkl(1,2)-brjk(1,1)*brkl(3,2))+ &
              (brij(1,1)*brjk(2,2)-brij(2,1)*brjk(1,2))* &
              (brjk(1,1)*brkl(2,2)-brjk(2,1)*brkl(1,2))) &
             /(sinfi2*sinfi3)
    sintau = (brij(1,2)*(brjk(2,1)*brkl(3,2)-brjk(3,1)*brkl(2,2)) &
              +brij(2,2)*(brjk(3,1)*brkl(1,2)-brjk(1,1)*brkl(3,2)) &
              +brij(3,2)*(brjk(1,1)*brkl(2,2)-brjk(2,1)*brkl(1,2))) &
             /(sinfi2*sinfi3)
    tau = atan2(sintau,costau)
    if (abs(tau) .eq. pi) tau = pi
    do ix = 1,3
      iy = ix+1
      if (iy .gt. 3) iy = iy-3
      iz = iy+1
      if (iz .gt. 3) iz = iz-3
      bt(ix,1) = (brij(iy,2)*brjk(iz,2)-brij(iz,2)*brjk(iy,2)) &
         &           /(rij1*sinfi2**2)
      bt(ix,4) = (brkl(iy,1)*brjk(iz,1)-brkl(iz,1)*brjk(iy,1)) &
         &           /(rkl1*sinfi3**2)
      bt(ix,2) = -((rjk1-rij1*cosfi2)*bt(ix,1) &
         &             +rkl1*cosfi3*bt(ix,4))/rjk1
      bt(ix,3) = -(bt(ix,1)+bt(ix,2)+bt(ix,4))
    end do
  end subroutine trsn2
  pure subroutine strtch2(xyz,avst,b)
    implicit none
    real(wp),intent(out) :: b(3,2)
    real(wp),intent(in)  :: xyz(3,2)
    real(wp) :: r(3)
    real(wp) :: rr
    real(wp),intent(out) :: avst
    r = xyz(:,2)-xyz(:,1)
    rr = norm2(r)
    avst = rr
    b(:,1) = -r/rr
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
      co = co+brij(i,1)*brjk(i,2)
      crap = crap+(brjk(i,2)+brij(i,1))**2
    end do
    if (sqrt(crap) .lt. 1.0d-6) then
      fir = pi-asin(sqrt(crap))
      si = sqrt(crap)
    else
      fir = acos(co)
      si = sqrt(1.0_wp-co**2)
    end if
    if (abs(fir-pi) .lt. 1.0d-13) then
      fir = pi
      return
    end if
    do i = 1,3
      bf(i,1) = (co*brij(i,1)-brjk(i,2))/(si*rij1)
      bf(i,3) = (co*brjk(i,2)-brij(i,1))/(si*rjk1)
      bf(i,2) = -(bf(i,1)+bf(i,3))
    end do
  end subroutine bend2

  pure subroutine outofp2(xyz,teta,bt)
    implicit none
    real(wp),intent(out) :: teta
    real(wp),intent(out) :: bt(3,4)
    real(wp),intent(in)  :: xyz(3,4)
    real(wp) :: r1(3),r2(3),r3(3)
    real(wp) :: q41,q42,q43,e41(3),e42(3),e43(3)
    real(wp) :: cosfi1,fi1,dfi1,cosfi2,fi2,dfi2,cosfi3,fi3,dfi3
    real(wp) :: c14(3,3),br14(3,3)
    real(wp) :: r42(3),r43(3)
    integer  :: ix,iy,iz
!  4 -> 1 (bond)
    r1 = xyz(:,1)-xyz(:,4)
    q41 = norm2(r1)
    e41 = r1/q41
!  4 -> 2 (bond in plane)
    r2 = xyz(:,2)-xyz(:,4)
    q42 = norm2(r2)
    e42 = r2/q42
!  4 -> 3 (bond in plane)
    r3 = xyz(:,3)-xyz(:,4)
    q43 = norm2(r3)
    e43 = r3/q43
!
!  get the angle between e43 and e42
!
    cosfi1 = dot_product(e43,e42)

    fi1 = acos(cosfi1)
    dfi1 = 180.d0*fi1/pi
!
!  dirty exit! this happens when an earlier structure is ill defined.
!
    if (abs(fi1-pi) .lt. 1.0d-13) then
      teta = 0.0_wp
      bt = 0.0_wp
      return
    end if
!
!  get the angle between e41 and e43
!
    cosfi2 = dot_product(e41,e43)

    fi2 = acos(cosfi2)
    dfi2 = 180.d0*fi2/pi
!
!  get the angle between e41 and e42
!
    cosfi3 = dot_product(e41,e42)

    fi3 = acos(cosfi3)
    dfi3 = 180.d0*fi3/pi
!
!  the first two centers are trivially
!
    c14(:,1) = xyz(:,1)
    c14(:,2) = xyz(:,4)
!
!  the 3rd is
!
    r42 = xyz(:,2)-xyz(:,4)
    r43 = xyz(:,3)-xyz(:,4)
    c14(1,3) = r42(2)*r43(3)-r42(3)*r43(2)
    c14(2,3) = r42(3)*r43(1)-r42(1)*r43(3)
    c14(3,3) = r42(1)*r43(2)-r42(2)*r43(1)
!
!  exit if 2-3-4 are collinear
!  (equivalent to the above check, but this is more concrete)
!
    if ((c14(1,3)**2+c14(2,3)**2+c14(3,3)**2) .lt. 1.0d-10) then
      teta = 0.0d0
      bt = 0.0_wp
      return
    end if
    c14(1,3) = c14(1,3)+xyz(1,4)
    c14(2,3) = c14(2,3)+xyz(2,4)
    c14(3,3) = c14(3,3)+xyz(3,4)

    call bend2(c14,teta,br14)

    teta = teta-0.5_wp*pi
!
!--compute the wdc matrix
!
    do ix = 1,3
      iy = mod(ix+1,4)+(ix+1)/4
      iz = mod(iy+1,4)+(iy+1)/4

      bt(ix,1) = -br14(ix,1)
      bt(ix,2) = r43(iz)*br14(iy,3)-r43(iy)*br14(iz,3)
      bt(ix,3) = -r42(iz)*br14(iy,3)+r42(iy)*br14(iz,3)

      bt(ix,4) = -(bt(ix,1)+bt(ix,2)+bt(ix,3))

    end do

    bt = -bt
  end subroutine outofp2

  Subroutine Trsn(xyz,nCent,Tau,Bt,lWrite,lWarn,Label,dBt,ldB)
!************************************************************************
!*                                                                      *
!* Reference: Molecular Vibrations, E. Bright Wilson, Jr, J. C. Decicius*
!*             nd Paul C. Cross, Sec. 4-1, Eq. 20-24                    *
!*                                                                      *
!* R.Lindh May-June '96                                                 *
!************************************************************************
    Implicit Real(wp) (a-h,o-z)

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
    CosTau = ((BRij(2,1)*BRjk(3,2)-BRij(3,1)*BRjk(2,2))*&
       &           (BRjk(2,1)*BRkl(3,2)-BRjk(3,1)*BRkl(2,2))+&
       &           (BRij(3,1)*BRjk(1,2)-BRij(1,1)*BRjk(3,2))*&
       &           (BRjk(3,1)*BRkl(1,2)-BRjk(1,1)*BRkl(3,2))+&
       &           (BRij(1,1)*BRjk(2,2)-BRij(2,1)*BRjk(1,2))*&
       &           (BRjk(1,1)*BRkl(2,2)-BRjk(2,1)*BRkl(1,2)))&
       &         /(SinFi2*SinFi3)
    !
    !     For the vector product of the two vectors. This
    !     will give a vector parallell to e23. The direction
    !     relative to e23 defines the sign.
    !
    !     e123 X e234 = SinTau * e23
    !
    SinTau = (BRij(1,2)*(BRjk(2,1)*BRkl(3,2)-BRjk(3,1)*BRkl(2,2))&
       &         +BRij(2,2)*(BRjk(3,1)*BRkl(1,2)-BRjk(1,1)*BRkl(3,2))&
       &         +BRij(3,2)*(BRjk(1,1)*BRkl(2,2)-BRjk(2,1)*BRkl(1,2)))&
       &         /(SinFi2*SinFi3)
    !
    !     (-Pi < Tau <= Pi)
    !
    Tau = ATan2(SinTau,CosTau)
    If (Abs(Tau) .eq. Pi) Tau = Pi
    !
    dTau = 180.0D+00*Tau/Pi
    dFi2 = 180.0D+00*Fi2/Pi
    dFi3 = 180.0D+00*Fi3/Pi
    If (lWarn) Then
      If (dTau .gt. 177.5.or.dTau .lt. -177.5) Then
        Write (*,*) ' Warning: dihedral angle close to'&
           &         //' end of range'
      End If
      If (dFi2 .gt. 177.5.or.dFi2 .lt. 2.5) Then
        Write (*,*) ' Warning: bond angle close to'&
           &         //' end of range'
      End If
      If (dFi3 .gt. 177.5.or.dFi3 .lt. 2.5) Then
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
      iy = ix+1
      If (iy .gt. 3) iy = iy-3
      iz = iy+1
      If (iz .gt. 3) iz = iz-3
      Bt(ix,1) = (BRij(iy,2)*BRjk(iz,2)-BRij(iz,2)*BRjk(iy,2))&
         &           /(Rij1*SinFi2**2)
      Bt(ix,4) = (BRkl(iy,1)*BRjk(iz,1)-BRkl(iz,1)*BRjk(iy,1))&
         &           /(Rkl1*SinFi3**2)
      Bt(ix,2) = -((Rjk1-Rij1*CosFi2)*Bt(ix,1)&
         &             +Rkl1*CosFi3*Bt(ix,4))/Rjk1
      Bt(ix,3) = -(Bt(ix,1)+Bt(ix,2)+Bt(ix,4))
    End Do
    !
    If (ldB) Then
      !
      !------- Compute the derivative of the WDC matrix.
      !
      Do ix = 1,3
        iy = ix+1
        If (iy .gt. 3) iy = iy-3
        iz = iy+1
        If (iz .gt. 3) iz = iz-3
        Do jx = 1,ix
          jy = jx+1
          If (jy .gt. 3) jy = jy-3
          jz = jy+1
          If (jz .gt. 3) jz = jz-3
          !
          dBt(ix,1,jx,1) = (dBRij(ix,1,jy,2)*BRjk(jz,2)&
             &                       -dBRij(ix,1,jz,2)*BRjk(jy,2)&
             &                       -Bt(jx,1)*(BRij(ix,1)*SinFi2**2&
             &                       +Rij1*Two*SinFi2*CosFi2*Bf2(ix,1)))&
             &                       /(Rij1*SinFi2**2)
          dBt(ix,1,jx,2) = -((-BRij(ix,1)*CosFi2&
             &                         +Rij1*SinFi2*Bf2(ix,1))*Bt(jx,1)&
             &                         +(Rjk1-Rij1*CosFi2)*dBt(ix,1,jx,1))&
             &                       /Rjk1
          dBt(jx,2,ix,1) = dBt(ix,1,jx,2)
          dBt(ix,1,jx,4) = Zero
          dBt(jx,4,ix,1) = dBt(ix,1,jx,4)
          dBt(ix,1,jx,3) = -(dBt(ix,1,jx,1)+dBt(ix,1,jx,2))
          dBt(jx,3,ix,1) = dBt(ix,1,jx,3)
          dBt(ix,4,jx,4) = (dBRkl(ix,2,jy,1)*BRjk(jz,1)&
             &                       -dBRkl(ix,2,jz,1)*BRjk(jy,1)&
             &                       -Bt(jx,4)*(BRkl(ix,2)*SinFi3**2&
             &                       +Rkl1*Two*SinFi3*CosFi3*Bf3(ix,3)))&
             &                       /(Rkl1*SinFi3**2)
          dBt(ix,4,jx,3) = -((-BRkl(ix,2)*CosFi3&
             &                         +Rkl1*SinFi3*Bf3(ix,3))*Bt(jx,4)&
             &                         +(Rjk1-Rkl1*CosFi3)*dBt(ix,4,jx,4))&
             &                       /Rjk1
          dBt(jx,3,ix,4) = dBt(ix,4,jx,3)
          dBt(ix,4,jx,2) = -(dBt(ix,4,jx,4)+dBt(ix,4,jx,3))
          dBt(jx,2,ix,4) = dBt(ix,4,jx,2)
          If (ix .ne. jx) Then
            dBt(jx,1,ix,1) = dBt(ix,1,jx,1)
            dBt(ix,4,jx,1) = Zero
            dBt(jx,4,ix,4) = dBt(ix,4,jx,4)
            dBt(jx,1,ix,4) = dBt(ix,4,jx,1)
            dBt(jx,1,ix,2) = -((-BRij(jx,1)*CosFi2&
               &                            +Rij1*SinFi2*Bf2(jx,1))*Bt(ix,1)&
               &                            +(Rjk1-Rij1*CosFi2)*dBt(jx,1,ix,1))&
               &                          /Rjk1
            dBt(ix,2,jx,1) = dBt(jx,1,ix,2)
            dBt(ix,3,jx,1) = -(dBt(ix,1,jx,1)+dBt(ix,2,jx,1)&
               &                          +dBt(ix,4,jx,1))
            dBt(jx,1,ix,3) = dBt(ix,3,jx,1)
            dBt(jx,4,ix,3) = -((-BRkl(jx,2)*CosFi3&
               &                            +Rkl1*SinFi3*Bf3(jx,3))*Bt(ix,4)&
               &                            +(Rjk1-Rkl1*CosFi3)*dBt(jx,4,ix,4))&
               &                          /Rjk1
            dBt(ix,3,jx,4) = dBt(jx,4,ix,3)
            dBt(ix,2,jx,4) = -(dBt(ix,4,jx,4)+dBt(ix,3,jx,4))
            dBt(jx,4,ix,2) = dBt(ix,2,jx,4)
          End If
          dBt(ix,2,jx,3) = -((BRjk(ix,1)&
             &                           +Rkl1*SinFi3*Bf3(ix,1))*Bt(jx,4)&
             &                         +(Rjk1-Rkl1*CosFi3)*dBt(ix,2,jx,4)&
             &                         +(BRij(ix,2)*CosFi2&
             &                           -Rij1*SinFi2*Bf2(ix,2))*Bt(jx,1)&
             &                         +Rij1*CosFi2*dBt(ix,2,jx,1)&
             &                         +Bt(jx,3)*BRjk(ix,1))/Rjk1
          dBt(jx,3,ix,2) = dBt(ix,2,jx,3)
          dBt(ix,2,jx,2) = -(dBt(ix,2,jx,1)+dBt(ix,2,jx,4)&
             &                         +dBt(ix,2,jx,3))
          dBt(ix,3,jx,3) = -(dBt(ix,2,jx,3)+dBt(ix,1,jx,3)&
             &                         +dBt(ix,4,jx,3))
          If (ix .ne. jx) Then
            dBt(ix,3,jx,2) = -(dBt(ix,2,jx,2)+dBt(ix,1,jx,2)&
               &                            +dBt(ix,4,jx,2))
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
      Implicit Real(wp) (a-h,o-z)
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
      R(1) = xyz(1,2)-xyz(1,1)
      R(2) = xyz(2,2)-xyz(2,1)
      R(3) = xyz(3,2)-xyz(3,1)
      R2 = R(1)**2+R(2)**2+R(3)**2
      RR = Sqrt(R2)
      Avst = RR
      !
      aRR = RR*Angstr
      If (lWrite) Write (*,'(1X,A,A,2(F10.6,A))') Label,&
         &      ' : Bond Length=',aRR,' / Angstrom',RR,' / bohr'
      !
      !---- Compute the WDC B-matrix.
      !
      B(1,1) = -R(1)/RR
      B(2,1) = -R(2)/RR
      B(3,1) = -R(3)/RR
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
              dB(i,1,j,1) = (One-B(j,1)*B(i,1))/RR
            Else
              dB(i,1,j,1) = (-B(j,1)*B(i,1))/RR
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
      Implicit Real(wp) (a-h,o-z)

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
        Co = Co+BRij(i,1)*BRjk(i,2)
        Crap = Crap+(BRjk(i,2)+BRij(i,1))**2
      End Do
      !
      !.... Special care for cases close to linearity
      !
      If (Sqrt(Crap) .lt. 1.0D-6) Then
        Fir = Pi-ArSin(Sqrt(Crap))
        Si = Sqrt(Crap)
      Else
        Fir = ArCos(Co)
        Si = Sqrt(One-Co**2)
      End If
      !
      If (Abs(Fir-Pi) .lt. 1.0d-13) Then
        Fir = Pi
        Return
      End If
      dFir = 180.0D0*Fir/Pi
      If ((Abs(dFir) .gt. 177.5.or.Abs(dFir) .lt. 2.5).and.lWarn)&
         &   Write (*,*) ' Valence angle close to end in '//&
         &               'range of definition'
      If (lWrite) Write (*,'(1X,A,A,F10.4,A,F10.6,A)') Label,&
         &            ' : Angle=',dFir,'/degree, ',Fir,'/rad'
      !
      !---- Compute the WDC B-matrix
      !
      !     Bf=-11.1111
      Do i = 1,3
        Bf(i,1) = (Co*BRij(i,1)-BRjk(i,2))/(Si*Rij1)
        Bf(i,3) = (Co*BRjk(i,2)-BRij(i,1))/(Si*Rjk1)
        !....... Utilize translational invariance.
        Bf(i,2) = -(Bf(i,1)+Bf(i,3))
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
            dBf(i,1,j,1) = (-Si*Bf(i,1)*BRij(j,1)&
               &                        +Co*dBRij(i,1,j,1)&
               &                        -Bf(j,1)*(Co*Bf(i,1)*Rij1&
               &                        +Si*BRij(i,1)))/(Si*Rij1)
            dBf(i,1,j,3) = (-Si*Bf(i,1)*BRjk(j,2)&
               &                       +dBRij(i,1,j,2)&
               &                       -Bf(j,3)*Co*Bf(i,1)*Rjk1)&
               &                       /(Si*Rjk1)
            !              Write (*,*) '13',dBf(i,1,j,3), i, j
            dBf(i,3,j,1) = (-Si*Bf(i,3)*BRij(j,1)&
               &                       +dBRjk(i,2,j,1)&
               &                       -Bf(j,1)*Co*Bf(i,3)*Rij1)&
               &                       /(Si*Rij1)
            dBf(i,3,j,3) = (-Si*Bf(i,3)*BRjk(j,2)&
               &                        +Co*dBRjk(i,2,j,2)&
               &                        -Bf(j,3)*(Co*Bf(i,3)*Rjk1&
               &                        +Si*BRjk(i,2)))/(Si*Rjk1)
            !
            dBf(j,1,i,1) = dBf(i,1,j,1)
            dBf(j,3,i,1) = dBf(i,1,j,3)
            dBf(j,1,i,3) = dBf(i,3,j,1)
            dBf(j,3,i,3) = dBf(i,3,j,3)
            !
            dBf(i,1,j,2) = -(dBf(i,1,j,1)+dBf(i,1,j,3))
            dBf(j,2,i,1) = dBf(i,1,j,2)
            dBf(j,1,i,2) = -(dBf(j,1,i,1)+dBf(j,1,i,3))
            dBf(i,2,j,1) = dBf(j,1,i,2)
            dBf(i,3,j,2) = -(dBf(i,3,j,1)+dBf(i,3,j,3))
            dBf(j,2,i,3) = dBf(i,3,j,2)
            dBf(j,3,i,2) = -(dBf(j,3,i,1)+dBf(j,3,i,3))
            dBf(i,2,j,3) = dBf(j,3,i,2)
            !
            dBf(i,2,j,2) = -(dBf(i,2,j,1)+dBf(i,2,j,3))
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
      Implicit Real*8(a-h,o-z)
      Real*8 ArSin

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
      Implicit Real(wp) (a-h,o-z)
      Real(wp) :: ArCos
      A = Arg
      IF (ABS(A) .GT. One) Then
        A = Sign(One,A)
      End If
      ArCos = ACos(A)
      Return
    End function arCos
  End subroutine trsn

  pure elemental function ixyz(i,iatom)
    integer :: ixyz
    integer,intent(in) :: i,iatom
    ixyz = (iatom-1)*3+i
  end function ixyz
  pure elemental function jnd(i,j)
    integer :: jnd
    integer,intent(in) :: i,j
    jnd = i*(i-1)/2+j
  end function jnd
  pure elemental function ind(i,iatom,j,jatom)
    integer :: ind
    integer,intent(in) :: i,iatom,j,jatom
    ind = jnd(max(ixyz(i,iatom),ixyz(j,jatom)),min(ixyz(i,iatom),ixyz(j,jatom)))
  end function ind

  pure elemental function fk_lindh(alpha,r0,r2) result(gmm)
    implicit none
    real(wp),intent(in) :: alpha,r0,r2
    real(wp) :: gmm
    gmm = exp(alpha*(r0**2-r2))
  end function fk_lindh

  pure elemental function fk_swart(alpha,r0,r2) result(gmm)
    implicit none
    real(wp),intent(in) :: alpha,r0,r2
    real(wp) :: gmm
    gmm = exp(-alpha*(sqrt(r2)/r0-1.0_wp))
  end function fk_swart

  pure elemental function fk_vdw(alpha,r0,r2) result(gmm)
    implicit none
    real(wp),intent(in) :: alpha,r0,r2
    real(wp) :: gmm
    gmm = exp(-alpha*(r0-sqrt(r2))**2)
  end function fk_vdw

!========================================================================================!
!########################################################################################!
!========================================================================================!

  subroutine mh_eeq(n,at,xyz,chrg,kq,hess)
    implicit none

!! ------------------------------------------------------------------------
!  Input
!! ------------------------------------------------------------------------
    integer,intent(in)     :: n                ! number of atoms
    integer,intent(in)     :: at(n)            ! ordinal numbers
    real(wp),intent(in)    :: xyz(3,n)         ! geometry
    real(wp),intent(in)    :: chrg             ! total charge
    real(wp),intent(in)    :: kq               ! scaling parameter
!    type(chrg_parameter),intent(in) :: chrgeq  ! charge model
!! ------------------------------------------------------------------------
!  Output
!! ------------------------------------------------------------------------
    real(wp),intent(out)   :: hess((3*n)*(3*n+1)/2)
    real(wp),allocatable   :: hessian(:,:,:,:) ! molecular hessian of IES

!  π itself
    real(wp),parameter :: pi = 3.1415926535897932384626433832795029_wp
!  √π
    real(wp),parameter :: sqrtpi = sqrt(pi)
!  √(2/π)
    real(wp),parameter :: sqrt2pi = sqrt(2.0_wp/pi)
!
!! ------------------------------------------------------------------------
!  charge model
!! ------------------------------------------------------------------------
    integer  :: m ! dimension of the Lagrangian
    real(wp),allocatable :: Amat(:,:)
    real(wp),allocatable :: Xvec(:)
    real(wp),allocatable :: Ainv(:,:)
    real(wp),allocatable :: dAmat(:,:,:)
    real(wp),allocatable :: dqdr(:,:,:)

!! ------------------------------------------------------------------------
!  local variables
!! ------------------------------------------------------------------------
    integer  :: i,j,k,l
    real(wp) :: r,rij(3),r2
    real(wp) :: gamij,gamij2
    real(wp) :: arg,arg2,tmp,dtmp
    real(wp) :: lambda
    real(wp) :: es,expterm,erfterm
    real(wp) :: htmp,rxr(3,3)
    real(wp) :: rcovij,rr

!! ------------------------------------------------------------------------
!  scratch variables
!! ------------------------------------------------------------------------
    real(wp),allocatable :: alpha(:)
    real(wp),allocatable :: xtmp(:)
    real(wp),allocatable :: atmp(:,:)

!! ------------------------------------------------------------------------
!  Lapack work variables
!! ------------------------------------------------------------------------
    integer,allocatable :: ipiv(:)
    real(wp),allocatable :: temp(:)
    real(wp),allocatable :: work(:)
    integer  :: lwork
    integer  :: info
    real(wp) :: test(1)

!! ------------------------------------------------------------------------
!  EEQ parameters
!  PARAMETRISATION BY S. SPICHER (Fri, 14 Dec 2018 16:13:08 +0100)
!! ------------------------------------------------------------------------
    integer,parameter :: max_elem = 86
!&<
    real(wp),parameter :: enparam(max_elem) = (/ &
     1.23695041_wp, 1.26590957_wp, 0.54341808_wp, 0.99666991_wp, 1.26691604_wp, &
     1.40028282_wp, 1.55819364_wp, 1.56866440_wp, 1.57540015_wp, 1.15056627_wp, &
     0.55936220_wp, 0.72373742_wp, 1.12910844_wp, 1.12306840_wp, 1.52672442_wp, &
     1.40768172_wp, 1.48154584_wp, 1.31062963_wp, 0.40374140_wp, 0.75442607_wp, &
     0.76482096_wp, 0.98457281_wp, 0.96702598_wp, 1.05266584_wp, 0.93274875_wp, &
     1.04025281_wp, 0.92738624_wp, 1.07419210_wp, 1.07900668_wp, 1.04712861_wp, &
     1.15018618_wp, 1.15388455_wp, 1.36313743_wp, 1.36485106_wp, 1.39801837_wp, &
     1.18695346_wp, 0.36273870_wp, 0.58797255_wp, 0.71961946_wp, 0.96158233_wp, &
     0.89585296_wp, 0.81360499_wp, 1.00794665_wp, 0.92613682_wp, 1.09152285_wp, &
     1.14907070_wp, 1.13508911_wp, 1.08853785_wp, 1.11005982_wp, 1.12452195_wp, &
     1.21642129_wp, 1.36507125_wp, 1.40340000_wp, 1.16653482_wp, 0.34125098_wp, &
     0.58884173_wp, 0.68441115_wp, 0.56999999_wp, 0.56999999_wp, 0.56999999_wp, &
     0.56999999_wp, 0.56999999_wp, 0.56999999_wp, 0.56999999_wp, 0.56999999_wp, &
     0.56999999_wp, 0.56999999_wp, 0.56999999_wp, 0.56999999_wp, 0.56999999_wp, &
     0.56999999_wp, 0.87936784_wp, 1.02761808_wp, 0.93297476_wp, 1.10172128_wp, &
     0.97350071_wp, 1.16695666_wp, 1.23997927_wp, 1.18464453_wp, 1.14191734_wp, &
     1.12334192_wp, 1.01485321_wp, 1.12950808_wp, 1.30804834_wp, 1.33689961_wp, &
     1.27465977_wp /)
    real(wp),parameter :: gamparam(max_elem) = (/ &
    -0.35015861_wp, 1.04121227_wp, 0.09281243_wp, 0.09412380_wp, 0.26629137_wp, &
     0.19408787_wp, 0.05317918_wp, 0.03151644_wp, 0.32275132_wp, 1.30996037_wp, &
     0.24206510_wp, 0.04147733_wp, 0.11634126_wp, 0.13155266_wp, 0.15350650_wp, &
     0.15250997_wp, 0.17523529_wp, 0.28774450_wp, 0.42937314_wp, 0.01896455_wp, &
     0.07179178_wp,-0.01121381_wp,-0.03093370_wp, 0.02716319_wp,-0.01843812_wp, &
    -0.15270393_wp,-0.09192645_wp,-0.13418723_wp,-0.09861139_wp, 0.18338109_wp, &
     0.08299615_wp, 0.11370033_wp, 0.19005278_wp, 0.10980677_wp, 0.12327841_wp, &
     0.25345554_wp, 0.58615231_wp, 0.16093861_wp, 0.04548530_wp,-0.02478645_wp, &
     0.01909943_wp, 0.01402541_wp,-0.03595279_wp, 0.01137752_wp,-0.03697213_wp, &
     0.08009416_wp, 0.02274892_wp, 0.12801822_wp,-0.02078702_wp, 0.05284319_wp, &
     0.07581190_wp, 0.09663758_wp, 0.09547417_wp, 0.07803344_wp, 0.64913257_wp, &
     0.15348654_wp, 0.05054344_wp, 0.11000000_wp, 0.11000000_wp, 0.11000000_wp, &
     0.11000000_wp, 0.11000000_wp, 0.11000000_wp, 0.11000000_wp, 0.11000000_wp, &
     0.11000000_wp, 0.11000000_wp, 0.11000000_wp, 0.11000000_wp, 0.11000000_wp, &
     0.11000000_wp,-0.02786741_wp, 0.01057858_wp,-0.03892226_wp,-0.04574364_wp, &
    -0.03874080_wp,-0.03782372_wp,-0.07046855_wp, 0.09546597_wp, 0.21953269_wp, &
     0.02522348_wp, 0.15263050_wp, 0.08042611_wp, 0.01878626_wp, 0.08715453_wp, &
     0.10500484_wp /)
    real(wp),parameter :: kappa(max_elem) = (/ &
     0.04916110_wp, 0.10937243_wp,-0.12349591_wp,-0.02665108_wp,-0.02631658_wp, &
     0.06005196_wp, 0.09279548_wp, 0.11689703_wp, 0.15704746_wp, 0.07987901_wp, &
     -0.10002962_wp,-0.07712863_wp,-0.02170561_wp,-0.04964052_wp, 0.14250599_wp, &
     0.07126660_wp, 0.13682750_wp, 0.14877121_wp,-0.10219289_wp,-0.08979338_wp, &
    -0.08273597_wp,-0.01754829_wp,-0.02765460_wp,-0.02558926_wp,-0.08010286_wp, &
    -0.04163215_wp,-0.09369631_wp,-0.03774117_wp,-0.05759708_wp, 0.02431998_wp, &
    -0.01056270_wp,-0.02692862_wp, 0.07657769_wp, 0.06561608_wp, 0.08006749_wp, &
     0.14139200_wp,-0.05351029_wp,-0.06701705_wp,-0.07377246_wp,-0.02927768_wp, &
    -0.03867291_wp,-0.06929825_wp,-0.04485293_wp,-0.04800824_wp,-0.01484022_wp, &
     0.07917502_wp, 0.06619243_wp, 0.02434095_wp,-0.01505548_wp,-0.03030768_wp, &
     0.01418235_wp, 0.08953411_wp, 0.08967527_wp, 0.07277771_wp,-0.02129476_wp, &
    -0.06188828_wp,-0.06568203_wp,-0.11000000_wp,-0.11000000_wp,-0.11000000_wp, &
    -0.11000000_wp,-0.11000000_wp,-0.11000000_wp,-0.11000000_wp,-0.11000000_wp, &
    -0.11000000_wp,-0.11000000_wp,-0.11000000_wp,-0.11000000_wp,-0.11000000_wp, &
    -0.11000000_wp,-0.03585873_wp,-0.03132400_wp,-0.05902379_wp,-0.02827592_wp, &
    -0.07606260_wp,-0.02123839_wp, 0.03814822_wp, 0.02146834_wp, 0.01580538_wp, &
    -0.00894298_wp,-0.05864876_wp,-0.01817842_wp, 0.07721851_wp, 0.07936083_wp, &
     0.05849285_wp /)
    real(wp),parameter :: alphaparam(max_elem) = (/ &
     0.55159092_wp, 0.66205886_wp, 0.90529132_wp, 1.51710827_wp, 2.86070364_wp, &
     1.88862966_wp, 1.32250290_wp, 1.23166285_wp, 1.77503721_wp, 1.11955204_wp, &
     1.28263182_wp, 1.22344336_wp, 1.70936266_wp, 1.54075036_wp, 1.38200579_wp, &
     2.18849322_wp, 1.36779065_wp, 1.27039703_wp, 1.64466502_wp, 1.58859404_wp, &
     1.65357953_wp, 1.50021521_wp, 1.30104175_wp, 1.46301827_wp, 1.32928147_wp, &
     1.02766713_wp, 1.02291377_wp, 0.94343886_wp, 1.14881311_wp, 1.47080755_wp, &
     1.76901636_wp, 1.98724061_wp, 2.41244711_wp, 2.26739524_wp, 2.95378999_wp, &
     1.20807752_wp, 1.65941046_wp, 1.62733880_wp, 1.61344972_wp, 1.63220728_wp, &
     1.60899928_wp, 1.43501286_wp, 1.54559205_wp, 1.32663678_wp, 1.37644152_wp, &
     1.36051851_wp, 1.23395526_wp, 1.65734544_wp, 1.53895240_wp, 1.97542736_wp, &
     1.97636542_wp, 2.05432381_wp, 3.80138135_wp, 1.43893803_wp, 1.75505957_wp, &
     1.59815118_wp, 1.76401732_wp, 1.63999999_wp, 1.63999999_wp, 1.63999999_wp, &
     1.63999999_wp, 1.63999999_wp, 1.63999999_wp, 1.63999999_wp, 1.63999999_wp, &
     1.63999999_wp, 1.63999999_wp, 1.63999999_wp, 1.63999999_wp, 1.63999999_wp, &
     1.63999999_wp, 1.47055223_wp, 1.81127084_wp, 1.40189963_wp, 1.54015481_wp, &
     1.33721475_wp, 1.57165422_wp, 1.04815857_wp, 1.78342098_wp, 2.79106396_wp, &
     1.78160840_wp, 2.47588882_wp, 2.37670734_wp, 1.76613217_wp, 2.66172302_wp, &
     2.82773085_wp /)
!&>

!! ------------------------------------------------------------------------
!  initizialization
!! ------------------------------------------------------------------------
    m = n+1
    allocate (ipiv(m),source=0)
    allocate (Amat(m,m),Xvec(m),alpha(n),dqdr(3,n,m),source=0.0_wp)

!! ------------------------------------------------------------------------
!  set up the A matrix and X vector
!! ------------------------------------------------------------------------
!  αi -> alpha(i), ENi -> xi(i), κi -> kappa(i), Jii -> gam(i)
!  γij = 1/√(αi+αj)
!  Xi  = -ENi + κi·√CNi
!  Aii = Jii + 2/√π·γii
!  Aij = erf(γij·Rij)/Rij = 2/√π·F0(γ²ij·R²ij)
!! ------------------------------------------------------------------------
!  prepare some arrays
!$omp parallel default(none) &
!!$omp shared(n,at,chrgeq) &
!$omp shared(n,at) &
!$omp private(i) &
!$omp shared(Xvec,alpha)
!$omp do schedule(dynamic)
    do i = 1,n
!      Xvec(i) = -chrgeq%en(i)
!      alpha(i) = chrgeq%alpha(i)**2
      Xvec(i) = -enparam(at(i))
      alpha(i) = alphaparam(at(i))**2
    end do
!$omp enddo
!$omp endparallel

!$omp parallel default(none) &
!!$omp shared(n,at,xyz,chrgeq,alpha) &
!$omp shared(n,at,xyz,alpha) &
!$omp private(i,j,r,gamij) &
!$omp shared(Amat)
!$omp do schedule(dynamic)
    ! prepare A matrix
    do i = 1,n
      ! EN of atom i
      do j = 1,i-1
        r = sqrt(sum((xyz(:,j)-xyz(:,i))**2))
        gamij = 1.0_wp/sqrt(alpha(i)+alpha(j))
        Amat(j,i) = erf(gamij*r)/r
        Amat(i,j) = Amat(j,i)
      end do
!      Amat(i,i) = chrgeq%gam(i)+sqrt2pi/sqrt(alpha(i))
      Amat(i,i) = gamparam(at(i))+sqrt2pi/sqrt(alpha(i))
    end do
!$omp enddo
!$omp endparallel

!! ------------------------------------------------------------------------
!  solve the linear equations to obtain partial charges
!! ------------------------------------------------------------------------
    Amat(m,1:m) = 1.0_wp
    Amat(1:m,m) = 1.0_wp
    Amat(m,m) = 0.0_wp
    Xvec(m) = chrg
    ! generate temporary copy
    allocate (Atmp(m,m),source=Amat)
    allocate (Xtmp(m),source=Xvec)

    ! assume work space query, set best value to test after first dsysv call
    call dsysv('u',m,1,Atmp,m,ipiv,Xtmp,m,test,-1,info)
    lwork = int(test(1))
    allocate (work(lwork),source=0.0_wp)

    call dsysv('u',m,1,Atmp,m,ipiv,Xtmp,m,work,lwork,info)
    if (info > 0) error stop '** ERROR ** (goedecker_solve) DSYSV failed'

    if (abs(sum(Xtmp(:n))-chrg) > 1.e-6_wp) &
      error stop '** ERROR ** (goedecker_solve) charge constrain error'
    !print'(3f20.14)',Xtmp

!! ------------------------------------------------------------------------
!  calculate isotropic electrostatic (IES) energy
!! ------------------------------------------------------------------------
!  E = ∑i (ENi - κi·√CNi)·qi + ∑i (Jii + 2/√π·γii)·q²i
!      + ½ ∑i ∑j,j≠i qi·qj·2/√π·F0(γ²ij·R²ij)
!    = q·(½A·q - X)
!! ------------------------------------------------------------------------
!   work(:m) = Xvec
!   call dsymv('u',m,0.5_wp,Amat,m,Xtmp,1,-1.0_wp,work,1)
!   es = dot_product(Xtmp,work(:m))
!   energy = es + energy

!! ------------------------------------------------------------------------
!  calculate molecular gradient of the IES energy
!! ------------------------------------------------------------------------
!  dE/dRj -> g(:,j), ∂Xi/∂Rj -> -dcn(:,i,j), ½∂Aij/∂Rj -> dAmat(:,j,i)
!  dE/dR = (½∂A/∂R·q - ∂X/∂R)·q
!  ∂Aij/∂Rj = ∂Aij/∂Ri
!! ------------------------------------------------------------------------
    allocate (dAmat(3,n,m),source=0.0_wp)
!$omp parallel default(none) &
!$omp shared(n,xyz,alpha,Amat,Xtmp) &
!$omp private(i,j,rij,r2,gamij,arg,dtmp) &
!$omp reduction(+:dAmat)
!$omp do schedule(dynamic)
    do i = 1,n
      do j = 1,i-1
        rij = xyz(:,i)-xyz(:,j)
        r2 = sum(rij**2)
        gamij = 1.0_wp/sqrt(alpha(i)+alpha(j))
        arg = gamij**2*r2
        dtmp = 2.0_wp*gamij*exp(-arg)/(sqrtpi*r2)-Amat(j,i)/r2
        dAmat(:,i,i) = +dtmp*rij*Xtmp(j)+dAmat(:,i,i)
        dAmat(:,j,j) = -dtmp*rij*Xtmp(i)+dAmat(:,j,j)
        dAmat(:,i,j) = +dtmp*rij*Xtmp(i)
        dAmat(:,j,i) = -dtmp*rij*Xtmp(j)
      end do
    end do
!$omp enddo
!$omp endparallel

!! ------------------------------------------------------------------------
!  invert the A matrix using a Bunch-Kaufman factorization
!  A⁻¹ = (L·D·L^T)⁻¹ = L^T·D⁻¹·L
!! ------------------------------------------------------------------------
    allocate (Ainv(m,m),source=Amat)

    ! assume work space query, set best value to test after first dsytrf call
    call dsytrf('L',m,Ainv,m,ipiv,test,-1,info)
    if (int(test(1)) > lwork) then
      deallocate (work)
      lwork = int(test(1))
      allocate (work(lwork),source=0.0_wp)
    end if

    ! Bunch-Kaufman factorization A = L*D*L**T
    call dsytrf('L',m,Ainv,m,ipiv,work,lwork,info)
    if (info > 0) then
      error stop '** ERROR ** (goedecker_inversion) DSYTRF failed'

    end if

    ! A⁻¹ from factorized L matrix, save lower part of A⁻¹ in Ainv matrix
    ! Ainv matrix is overwritten with lower triangular part of A⁻¹
    call dsytri('L',m,Ainv,m,ipiv,work,info)
    if (info > 0) then
      error stop '** ERROR ** (goedecker_inversion) DSYTRI failed'
    end if

    ! symmetrizes A⁻¹ matrix from lower triangular part of inverse matrix
    do i = 1,m
      do j = i+1,m
        Ainv(i,j) = Ainv(j,i)
      end do
    end do

!! ------------------------------------------------------------------------
!  calculate gradient of the partial charge w.r.t. the nuclear coordinates
!! ------------------------------------------------------------------------
    !call dsymm('r','l',3*n,m,-1.0_wp,Ainv,m,dAmat,3*n,1.0_wp,dqdr,3*n)
    call dgemm('n','n',3*n,m,m,-1.0_wp,dAmat,3*n,Ainv,m,1.0_wp,dqdr,3*n)
    !print'(/,"analytical gradient")'
    !print'(3f20.14)',dqdr(:,:,:n)

!! ------------------------------------------------------------------------
!  molecular Hessian calculation
!! ------------------------------------------------------------------------
    do i = 1,n
      do j = 1,i-1
        rij = xyz(:,j)-xyz(:,i)
        r2 = sum(rij**2)
        r = sqrt(r2)
        gamij = 1.0_wp/sqrt(alpha(i)+alpha(j))
        gamij2 = gamij**2
        arg2 = gamij2*r2
        arg = sqrt(arg2)
        erfterm = Xtmp(i)*Xtmp(j)*erf(arg)/r
        expterm = Xtmp(i)*Xtmp(j)*2*gamij*exp(-arg2)/sqrtpi
        ! ∂²(qAq)/(∂Ri∂Rj):
        ! ∂²(qAq)/(∂Xi∂Xi) = (1-3X²ij/R²ij-2γ²ijX²ij) 2γij/√π exp[-γ²ij·R²ij]/R²ij
        !                  - (R²ij-3X²ij) erf[γij·Rij]/R⁵ij
        ! ∂²(qAq)/(∂Xi∂Xj) = (R²ij-3X²ij) erf[γij·Rij]/R⁵ij
        !                  - (1-3X²ij/R²ij-2γ²ijX²ij) 2γij/√π exp[-γ²ij·R²ij]/R²ij
        ! ∂²(qAq)/(∂Xi∂Yi) = 3X²ij erf[γij·Rij]/R⁵ij
        !                  - (3X²ij/R²ij+2γ²ijX²ij) 2γij/√π exp[-γ²ij·R²ij]/R²ij
        ! ∂²(qAq)/(∂Xi∂Yj) = (3X²ij/R²ij+2γ²ijX²ij) 2γij/√π exp[-γ²ij·R²ij]/R²ij
        !                  - 3X²ij erf[γij·Rij]/R⁵ij
        rxr(1,1) = erfterm*(3*rij(1)**2/r2**2-1.0_wp/r2) &
                   -expterm*(3*rij(1)**2/r2**2+2*gamij2*rij(1)**2/r2-1/r2)
        rxr(2,2) = erfterm*(3*rij(2)**2/r2**2-1.0_wp/r2) &
                   -expterm*(3*rij(2)**2/r2**2+2*gamij2*rij(2)**2/r2-1/r2)
        rxr(3,3) = erfterm*(3*rij(3)**2/r2**2-1.0_wp/r2) &
                   -expterm*(3*rij(3)**2/r2**2+2*gamij2*rij(3)**2/r2-1/r2)
        rxr(2,1) = erfterm*3*rij(2)*rij(1)/r2**2 &
                   -expterm*(3*rij(2)*rij(1)/r2**2+2*gamij2*rij(2)*rij(1)/r2)
        rxr(3,1) = erfterm*3*rij(3)*rij(1)/r2**2 &
                   -expterm*(3*rij(3)*rij(1)/r2**2+2*gamij2*rij(3)*rij(1)/r2)
        rxr(3,2) = erfterm*3*rij(3)*rij(2)/r2**2 &
                   -expterm*(3*rij(3)*rij(2)/r2**2+2*gamij2*rij(3)*rij(2)/r2)

        do k = 1,m
          rxr(1,1) = rxr(1,1)+0.5_wp*dqdr(1,i,k)*dAmat(1,j,k) &
                     +0.5_wp*dqdr(1,j,k)*dAmat(1,i,k)
          rxr(2,1) = rxr(2,1)+0.5_wp*dqdr(2,i,k)*dAmat(1,j,k) &
                     +0.5_wp*dqdr(2,j,k)*dAmat(1,i,k)
          rxr(3,1) = rxr(3,1)+0.5_wp*dqdr(3,i,k)*dAmat(1,j,k) &
                     +0.5_wp*dqdr(3,j,k)*dAmat(1,i,k)
          rxr(2,2) = rxr(2,2)+0.5_wp*dqdr(2,i,k)*dAmat(2,j,k) &
                     +0.5_wp*dqdr(2,j,k)*dAmat(2,i,k)
          rxr(3,2) = rxr(3,2)+0.5_wp*dqdr(3,i,k)*dAmat(2,j,k) &
                     +0.5_wp*dqdr(3,j,k)*dAmat(2,i,k)
          rxr(3,3) = rxr(3,3)+0.5_wp*dqdr(3,i,k)*dAmat(3,j,k) &
                     +0.5_wp*dqdr(3,j,k)*dAmat(3,i,k)
        end do
        ! symmetrize
        rxr(1,2) = rxr(2,1)
        rxr(1,3) = rxr(3,1)
        rxr(2,3) = rxr(3,2)

        ! save diagonal elements for atom i
        hess(ind(1,i,1,i)) = hess(ind(1,i,1,i))+kq*rxr(1,1)
        hess(ind(2,i,1,i)) = hess(ind(2,i,1,i))+kq*rxr(2,1)
        hess(ind(2,i,2,i)) = hess(ind(2,i,2,i))+kq*rxr(2,2)
        hess(ind(3,i,1,i)) = hess(ind(3,i,1,i))+kq*rxr(3,1)
        hess(ind(3,i,2,i)) = hess(ind(3,i,2,i))+kq*rxr(3,2)
        hess(ind(3,i,3,i)) = hess(ind(3,i,3,i))+kq*rxr(3,3)
        ! save elements between atom i and atom j
        hess(ind(1,i,1,j)) = hess(ind(1,i,1,j))-kq*rxr(1,1)
        hess(ind(1,i,2,j)) = hess(ind(1,i,2,j))-kq*rxr(2,1)
        hess(ind(1,i,3,j)) = hess(ind(1,i,3,j))-kq*rxr(3,1)
        hess(ind(2,i,1,j)) = hess(ind(2,i,1,j))-kq*rxr(2,1)
        hess(ind(2,i,2,j)) = hess(ind(2,i,2,j))-kq*rxr(2,2)
        hess(ind(2,i,3,j)) = hess(ind(2,i,3,j))-kq*rxr(3,2)
        hess(ind(3,i,1,j)) = hess(ind(3,i,1,j))-kq*rxr(3,1)
        hess(ind(3,i,2,j)) = hess(ind(3,i,2,j))-kq*rxr(3,2)
        hess(ind(3,i,3,j)) = hess(ind(3,i,3,j))-kq*rxr(3,3)
        ! save diagonal elements for atom j
        hess(ind(1,j,1,j)) = hess(ind(1,j,1,j))+kq*rxr(1,1)
        hess(ind(2,j,1,j)) = hess(ind(2,j,1,j))+kq*rxr(2,1)
        hess(ind(2,j,2,j)) = hess(ind(2,j,2,j))+kq*rxr(2,2)
        hess(ind(3,j,1,j)) = hess(ind(3,j,1,j))+kq*rxr(3,1)
        hess(ind(3,j,2,j)) = hess(ind(3,j,2,j))+kq*rxr(3,2)
        hess(ind(3,j,3,j)) = hess(ind(3,j,3,j))+kq*rxr(3,3)
      end do
    end do

    ! ∂²(qA)/(∂Ri∂q)·∂q/∂Rj
    ! hessian = hessian + reshape(matmul(reshape(dqdr,(/3*n,m/)),&
    !    transpose(reshape(dAmat,(/3*n,m/)))),(/3,n,3,n/))
    !call dgemm('n','t',3*n,m,3*n,+1.0_wp,dqdr,3*n,dAmat,3*n,1.0_wp,hessian,3*n)
    !call dgemm('n','t',3*n,m,3*n,+1.0_wp,dAmat,3*n,dqdr,3*n,1.0_wp,hessian,3*n)

  end subroutine mh_eeq

!========================================================================================!
!########################################################################################!
!========================================================================================!
end module modelhessian_core
