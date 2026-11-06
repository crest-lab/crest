!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2022-2026 Christoph Plett, Philipp Pracht
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
! Ported from the xtb IFF submodule (github.com/grimme-lab/xtb,
! src/iff/iff_ini.f90), LGPL-3.0-or-later.
!================================================================================!

!> Setup of the IFF pair data (D4 C6, Drude charges, Gaussian densities)
module aiss_iff_ini
  use crest_parameters,only:wp,stdout
  use aiss_param
  use aiss_iff_lmo,only:setlmo
#ifdef WITH_TBLITE
  use mctc_env,only:error_type
  use mctc_io,only:structure_type,new
  use dftd4_model,only:d4_model,new_d4_model,d4_qmod
  use dftd4_ncoord,only:get_coordination_number
#endif
  implicit none
  private

  public :: init_iff

  !> D3 C6AA(CNmax) coefficients
  real(wp),parameter :: c6d3(94) = [ &
  &  3.03_wp,1.56_wp,85.32_wp,55.14_wp,28.03_wp,18.21_wp,15.58_wp,10.37_wp, &
  &  7.13_wp,6.29_wp,186.11_wp,175.56_wp,153.59_wp,149.77_wp,151.69_wp,125.81_wp, &
  & 90.40_wp,64.65_wp,338.02_wp,343.33_wp,313.98_wp,321.59_wp,269.33_wp,153.34_wp, &
  &113.85_wp,109.50_wp,118.29_wp,141.88_wp,175.00_wp,153.45_wp,172.32_wp,180.29_wp, &
  &224.29_wp,210.52_wp,169.24_wp,130.40_wp,476.26_wp,491.74_wp,451.08_wp,458.69_wp, &
  &428.90_wp,337.11_wp,265.87_wp,239.01_wp,244.66_wp,265.87_wp,268.61_wp,263.45_wp, &
  &310.83_wp,336.10_wp,417.38_wp,418.58_wp,358.50_wp,290.22_wp,649.09_wp,716.55_wp, &
  &766.02_wp,688.04_wp,647.02_wp,620.52_wp,581.74_wp,587.79_wp,549.03_wp,595.36_wp, &
  &524.92_wp,544.44_wp,472.80_wp,446.23_wp,440.52_wp,487.08_wp,411.57_wp,420.98_wp, &
  &458.77_wp,398.10_wp,332.14_wp,297.83_wp,305.53_wp,336.87_wp,317.16_wp,298.86_wp, &
  &389.61_wp,455.29_wp,533.88_wp,558.97_wp,495.45_wp,412.83_wp,822.05_wp,867.56_wp, &
  &906.04_wp,1224.59_wp,980.38_wp,753.86_wp,698.16_wp,1051.25_wp]

!========================================================================================!
!========================================================================================!
contains  !> MODULE PROCEDURES START HERE
!========================================================================================!
!========================================================================================!

  subroutine init_iff(d,nfrag1,icoord,icoord0,iunit)
    !*****************************************************************
    !* Precompute all IFF pair quantities from the fragment data set
    !* by the LMO precomputation (at, xyz, q, LMOs, CT populations):
    !* D4 C6/alpha(0), BJ radii, neighbor lists, LMO point charges,
    !* Drude charges and Gaussian atomic densities. Both fragments
    !* are shifted to their CMA (A: CMA of the first nfrag1 atoms if
    !* nfrag1 > 0); icoord/icoord0 receive the relative placement.
    !*****************************************************************
    type(iff_data),intent(inout) :: d
    integer,intent(in)     :: nfrag1
    real(wp),intent(inout) :: icoord(6),icoord0(6)
    integer,intent(in)     :: iunit

    integer  :: i,j,k,n,n1,n2,np
    real(wp) :: distance,max_distanceA,max_distanceB
    real(wp) :: r(3),rr,rrr,r42,c6d4,c6aad4,tmp,thr,a1,a2,s8,ai,aj,p1,fc
    real(wp) :: f,fall,fhyd,r0,aa,e1(4),c1(4)
    real(wp),allocatable :: cn(:),alp0(:),dum1(:),dum2(:)

    n = d%n
    n1 = d%n1
    n2 = d%n2
    allocate (cn(n),alp0(n),dum1(n1),dum2(n2),source=0.0_wp)
    d%c6ab = 0.0_wp

!>-- combined system for D4
    do i = 1,n1
      d%xyz(1:3,i) = d%xyz1(1:3,i)
      d%at(i) = d%at1(i)
      d%q(i) = d%q1(i)
    end do
    do i = 1,n2
      d%xyz(1:3,i+n1) = d%xyz2(1:3,i)
      d%at(i+n1) = d%at2(i)
      d%q(i+n1) = d%q2(i)
    end do
    d%chrg(1) = sum(d%q1)+1.d-8
    d%chrg(2) = sum(d%q2)+1.d-8

!>-- shift B far away so that the D4 coefficients are those of separated fragments
    max_distanceA = 0.0_wp
    do i = 1,n1
      do j = 1,i-1
        distance = sqrt(sum((d%xyz1(:,j)-d%xyz1(:,i))**2))
        max_distanceA = max(max_distanceA,distance)
      end do
    end do
    max_distanceB = 0.0_wp
    do i = 1,n2
      do j = 1,i-1
        distance = sqrt(sum((d%xyz2(:,j)-d%xyz2(:,i))**2))
        max_distanceB = max(max_distanceB,distance)
      end do
    end do
    d%shift_geo = max_distanceA+max_distanceB+1.d+6
    do i = 1,n2
      d%xyz(1:3,i+n1) = d%xyz2(1:3,i)+d%shift_geo
    end do

!>-- free-atom radii and BJ damping terms
    a1 = par_d3_a1
    a2 = par_d3_a2
    s8 = par_d3_s8
    do i = 1,86
      call free_atom_alpha(i,tmp)
      r0_atom(i) = tmp**(1./3.)
      do j = 1,i
        r42 = r2r4(i)*r2r4(j)
        rrab(j,i) = r42*3.0d0*s8
        rrab(i,j) = r42*3.0d0*s8
        r0ab6(j,i) = (a1*dsqrt(3.0d0*r42)+a2)**6
        r0ab8(j,i) = (a1*dsqrt(3.0d0*r42)+a2)**8
        r0ab6(i,j) = r0ab6(j,i)
        r0ab8(i,j) = r0ab8(j,i)
      end do
    end do

    call c6_alp0(n,d%at,d%xyz,d%q,cn,alp0,d%c6ab)
    c6aad4 = sum(d%c6ab)

    c6d4 = 0
    do i = 1,n1
      d%cprob(i) = sqrt(d%c6ab(i,i)*c6d3(probe_atom_type))
      ai = alp0(i)**(1./3)
      do j = 1,n2
        aj = alp0(j+n1)**(1./3)
        !> converts to the r0ab scale so that element factors are about one
        d%alpab(j,i) = 2.6*(ai+aj)*0.5
        c6d4 = c6d4+d%c6ab(j+n1,i)
      end do
    end do
    if (d%debug) then
      write (stdout,'(''molecular C6AA D4 /au    :'',2F12.1)') c6aad4
      write (stdout,'(''molecular C6AB D4 /au    :'',2F12.1)') c6d4
    end if

!>-- ATM neighbor lists (intra cut-off thr) and Gamma matrices
    thr = 25.0
    do i = 1,n1
      d%alp1(i) = alp0(i)
      d%cn1(i) = cn(i)
      d%z1(i) = val_e(d%at1(i))
      k = 0
      do j = 1,n1
        rr = sum((d%xyz1(:,j)-d%xyz1(:,i))**2)
        r0 = sqrt(rr)
        d%gab1(j,i) = 1.0d0/(r0+1./sqrt(gam(d%at1(i))*gam(d%at1(j))))
        if (j .eq. i) cycle
        ai = alp0(i)**(1./3)
        aj = alp0(j)**(1./3)
        rrr = 1.13*2.6*(ai+aj)*0.5/r0
        tmp = 1.0d0/(1.d0+6.d0*rrr**14)
        if (tmp .gt. 0.90.and.r0 .lt. thr) then
          k = k+1
          d%neigh(k,i) = j
        end if
      end do
      d%neigh(0,i) = k
    end do
    do i = 1,n2
      d%alp2(i) = alp0(i+n1)
      d%cn2(i) = cn(i+n1)
      d%z2(i) = val_e(d%at2(i))
      k = 0
      do j = 1,n2
        rr = sum((d%xyz2(:,j)-d%xyz2(:,i))**2)
        r0 = sqrt(rr)
        d%gab2(j,i) = 1.0d0/(r0+1./sqrt(gam(d%at2(i))*gam(d%at2(j))))
        if (j .eq. i) cycle
        ai = alp0(i+n1)**(1./3)
        aj = alp0(j+n1)**(1./3)
        rrr = 1.13*2.6*(ai+aj)*0.5/r0
        tmp = 1.0d0/(1.d0+6.d0*rrr**14)
        if (tmp .gt. 0.90.and.r0 .lt. thr) then
          k = k+1
          d%neigh(k,i+n1) = j
        end if
      end do
      d%neigh(0,i+n1) = k
    end do

    call setlmo(1,n1,d%nlmo1,d%at1,d%xyz1,d%alp1,d%q1,d%cn1,d%rlmo1,d%lmo1,d%dipol(1),d%debug)
    call setlmo(2,n2,d%nlmo2,d%at2,d%xyz2,d%alp2,d%q2,d%cn2,d%rlmo2,d%lmo2,d%dipol(2),d%debug)

!>-- shift to the CMA; xyz1 stays fixed from here on
    icoord = 0
    if (nfrag1 .eq. 0) then
      call cmadock(n1,n1,d%at1,d%xyz1,r)
    else
      call cmadock(nfrag1,n1,d%at1,d%xyz1,r)
    end if
    icoord(1:3) = -r(1:3)
    do i = 1,3
      d%xyz1(i,1:n1) = d%xyz1(i,1:n1)-r(i)
      d%rlmo1(i,1:d%nlmo1) = d%rlmo1(i,1:d%nlmo1)-r(i)
    end do
    call cmadock(n2,n2,d%at2,d%xyz2,r)
    icoord(1:3) = icoord(1:3)+r(1:3)
    do i = 1,3
      d%xyz2(i,1:n2) = d%xyz2(i,1:n2)-r(i)
      d%rlmo2(i,1:d%nlmo2) = d%rlmo2(i,1:d%nlmo2)-r(i)
    end do
    icoord0 = icoord

!>-- Drude charges and 4-primitive Gaussian atomic densities
    np = 4
    fall = 0.50d0
    fhyd = 3.00d0
    p1 = 4.00
    do i = 1,n1
      f = fall
      if (d%at1(i) .le. 2) f = fhyd
      fc = par_drude_fc*gam(d%at1(i))
      d%qdr1(i) = -f*(fc*d%alp1(i))**0.50d0
      d%alp1(i) = d%alp1(i)**(1./3.)
      aa = p1*r0scal(d%at1(i))/d%alp1(i)
      if (d%at1(i) .eq. 1) then
        aa = aa*(1.+d%cn1(i))**(1./2.)
      else
        aa = aa*(1.+d%cn1(i))**(-1./3.)
      end if
      dum1(i) = aa
      call setsto4(np,1,1,aa,e1,c1)
      d%den1(1,1:np,i) = e1(1:np)
      d%den1(2,1:np,i) = c1(1:np)
    end do
    do i = 1,n2
      f = fall
      if (d%at2(i) .le. 2) f = fhyd
      fc = par_drude_fc*gam(d%at2(i))
      d%qdr2(i) = -f*(fc*d%alp2(i))**0.50d0
      d%alp2(i) = d%alp2(i)**(1./3.)
      aa = p1*r0scal(d%at2(i))/d%alp2(i)
      if (d%at2(i) .eq. 1) then
        aa = aa*(1.+d%cn2(i))**(1./2.)
      else
        aa = aa*(1.+d%cn2(i))**(-1./3.)
      end if
      dum2(i) = aa
      call setsto4(np,1,1,aa,e1,c1)
      d%den2(1,1:np,i) = e1(1:np)
      d%den2(2,1:np,i) = c1(1:np)
    end do

    if (d%debug) then
      write (stdout,*) '  #   Z           coordinates (CMA)   CN(D4)    qM/CM5','   alp   dexpo'
      do i = 1,n1
        write (stdout,'(2i4,3f10.5,3f6.2,3f7.3)') &
        & i,d%at1(i),d%xyz1(1:3,i),d%cn1(i),d%q1(i),d%q(i),d%alp1(i),dum1(i)
      end do
      write (stdout,*) 'LMOcent on atom     coordinates    charge'
      do i = 1,d%nlmo1
        write (stdout,'(2i7,3f10.5,5x,f7.3)') i,d%lmo1(i),d%rlmo1(1:4,i)
      end do
      write (stdout,*) '  #   Z           coordinates (CMA)   CN(D4)     qM/CM5','   alp   dexpo'
      do i = 1,n2
        write (stdout,'(2i4,3f10.5,3f6.2,3f7.3)') &
        & i,d%at2(i),d%xyz2(1:3,i),d%cn2(i),d%q2(i),d%q(i+n1),d%alp2(i),dum2(i)
      end do
      write (stdout,*) 'LMOcent on atom     coordinates    charge'
      do i = 1,d%nlmo2
        write (stdout,'(2i7,3f10.5,5x,f7.3)') i,d%lmo2(i),d%rlmo2(1:4,i)
      end do
    end if

    if (iunit /= 0) then
      write (iunit,'(''  LUMO energy 1 (read)           :'',F8.3)') d%elumo(1)*27.212
      write (iunit,'(''  HOMO energy 1 (read)           :'',F8.3)') d%ehomo(1)*27.212
      write (iunit,'(''  LUMO energy 2 (read)           :'',F8.3)') d%elumo(2)*27.212
      write (iunit,'(''  HOMO energy 2 (read)           :'',F8.3)') d%ehomo(2)*27.212
    end if
  end subroutine init_iff

!========================================================================================!

  subroutine setsto4(nprim,n,l,zeta,expo,cont)
    !*****************************************************************
    !* STO-4G expansion (normalized primitives) of an n,l Slater
    !* function with exponent zeta
    !*****************************************************************
    integer,intent(in)    :: n,l
    integer,intent(inout) :: nprim
    real(wp),intent(in)   :: zeta
    real(wp),intent(out)  :: expo(nprim),cont(nprim)
    real(wp) :: allc(4,5,4),allz(4,5,4),dex(-1:96),xnorm
    integer  :: iam,i,j,dex2

    nprim = 4
    do i = -1,10
      dex2 = 1
      if (i .ge. 2) then
        do j = 1,i,2
          dex2 = dex2*j
        end do
      end if
      dex(i) = dex2
    end do
    allz = 0.0_wp
    allc = 0.0_wp
    !> 1s
    allz(:,1,1) = [5.216844534_wp,9.546182760D-1,2.652034102D-1,8.801862774D-2]
    allc(:,1,1) = [5.675242080D-2,2.601413550D-1,5.328461143D-1,2.916254405D-1]
    !> 2s,2p
    allz(:,2,1) = [1.161525551D1,2.000243111_wp,1.607280687D-1,6.125744532D-2]
    allc(:,2,1) = [-1.198411747D-2,-5.472052539D-2,5.805587176D-1,4.770079976D-1]
    allz(:,2,2) = [1.798260992_wp,4.662622228D-1,1.643718620D-1,6.543927065D-2]
    allc(:,2,2) = [5.713170255D-2,2.857455515D-1,5.517873105D-1,2.632314924D-1]
    !> 3s,3p,3d
    allz(:,3,1) = [1.513265591_wp,4.262497508D-1,7.643320863D-2,3.760545063D-2]
    allc(:,3,1) = [-3.295496352D-2,-1.724516959D-1,7.518511194D-1,3.589627317D-1]
    allz(:,3,2) = [1.853180239_wp,1.915075719D-1,8.655487938D-2,4.184253862D-2]
    allc(:,3,2) = [-1.434249391D-2,2.755177580D-1,5.846750879D-1,2.144986514D-1]
    allz(:,3,3) = [9.185846715D-1,2.920461109D-1,1.187568890D-1,5.286755896D-2]
    allc(:,3,3) = [5.799057705D-2,3.045581349D-1,5.601358038D-1,2.432423313D-1]
    !> 4s,4p,4d,4f
    allz(:,4,1) = [3.242212833D-1,1.663217177D-1,5.081097451D-2,2.829066600D-2]
    allc(:,4,1) = [-1.120682822D-1,-2.845426863D-1,8.909873788D-1,3.517811205D-1]
    allz(:,4,2) = [1.492607880_wp,4.327619272D-1,7.553156064D-2,3.706272183D-2]
    allc(:,4,2) = [-6.035216774D-3,-6.013310874D-2,6.451518200D-1,4.117923820D-1]
    allz(:,4,3) = [1.995825422_wp,1.823461280D-1,8.197240896D-2,4.000634951D-2]
    allc(:,4,3) = [-2.816702620D-3,2.177095871D-1,6.058047348D-1,2.717811257D-1]
    !> 4f carries an extra factor, the normalization below is wrong for f
    allz(:,4,4) = [0.05691670217_wp,0.2074585819_wp,0.09298346885_wp,0.04473508853_wp]
    allc(:,4,4) = [0.05902730589_wp,0.3191828952_wp,0.5639423893_wp,0.2284796537_wp] &
    &             /1.063832358490576_wp**0.5
    !> 5s,5p,5d
    allz(:,5,1) = [8.602284252D-1,1.189050200D-1,3.446076176D-2,1.974798796D-2]
    allc(:,5,1) = [1.103657561D-2,-5.606519023D-1,1.179429987_wp,1.734974376D-1]
    allz(:,5,2) = [3.962838833D-1,1.838858552D-1,4.943555157D-2,2.750222273D-2]
    allc(:,5,2) = [-1.801459207D-2,-1.360777372D-1,7.533973719D-1,3.409304859D-1]
    allz(:,5,3) = [4.230617826D-1,8.293863702D-2,4.590326388D-2,2.628744797D-2]
    allc(:,5,3) = [-2.421626009D-2,3.937644956D-1,5.489520286D-1,1.190436963D-1]

    iam = l-1
    do j = 1,nprim
      cont(j) = allc(j,n,l)
      expo(j) = allz(j,n,l)*zeta**2
      xnorm = (2.d0*expo(j)/pi)**0.75d0*(4.d0*expo(j))**(iam/2.d0)/sqrt(dex(2*iam-1))
      cont(j) = cont(j)*xnorm
    end do
  end subroutine setsto4

!========================================================================================!

  subroutine c6_alp0(nat,at,xyz,q,cn,azero,c6ab)
    !*****************************************************************
    !* D4 (GFN2 reference charges) coordination numbers cn, pair
    !* C6 coefficients c6ab and static polarizabilities azero for
    !* the given atomic charges q (dftd4 library replaces the
    !* xtb-internal D4 model)
    !*****************************************************************
    integer,intent(in)   :: nat,at(nat)
    real(wp),intent(in)  :: xyz(3,nat),q(nat)
    real(wp),intent(out) :: cn(nat),azero(nat),c6ab(nat,nat)
#ifdef WITH_TBLITE
    type(structure_type) :: mol
    type(d4_model) :: d4
    type(error_type),allocatable :: error
    real(wp),allocatable :: gwvec(:,:,:)
    real(wp) :: trans(3,1)

    call new(mol,at,xyz)
    call new_d4_model(error,d4,mol,ga=3.0_wp,gc=2.0_wp,wf=6.0_wp,qmod=d4_qmod%gfn2)
    if (allocated(error)) then
      write (stdout,'(a,1x,a)') '[Error] D4 setup for IFF:',error%message
      error stop
    end if
    !> xtb used a squared distance threshold of 500 Bohr^2
    trans = 0.0_wp
    call get_coordination_number(mol,trans,sqrt(500.0_wp),d4%rcov,d4%en,cn)
    allocate (gwvec(maxval(d4%ref),nat,d4%ncoup),source=0.0_wp)
    call d4%weight_references(mol,cn,q,gwvec)
    call d4%get_atomic_c6(mol,gwvec,c6=c6ab)
    call d4%get_polarizabilities(mol,gwvec,alpha=azero)
#else
    cn = 0.0_wp
    azero = 0.0_wp
    c6ab = 0.0_wp
    write (stdout,'(a)') 'Error: the aISS docking requires CREST compiled with tblite (dftd4)'
    error stop
#endif
  end subroutine c6_alp0

  subroutine free_atom_alpha(iat,alpha)
    !> static polarizability of the free neutral atom (CN = 0, q = 0)
    integer,intent(in)   :: iat
    real(wp),intent(out) :: alpha
#ifdef WITH_TBLITE
    type(structure_type) :: mol
    type(d4_model) :: d4
    type(error_type),allocatable :: error
    real(wp),allocatable :: gwvec(:,:,:)
    real(wp) :: xyz(3,1),a(1)

    xyz = 0.0_wp
    call new(mol,[iat],xyz)
    call new_d4_model(error,d4,mol,ga=3.0_wp,gc=2.0_wp,wf=6.0_wp,qmod=d4_qmod%gfn2)
    if (allocated(error)) then
      write (stdout,'(a,1x,a)') '[Error] D4 setup for IFF:',error%message
      error stop
    end if
    allocate (gwvec(maxval(d4%ref),1,d4%ncoup),source=0.0_wp)
    call d4%weight_references(mol,[0.0_wp],[0.0_wp],gwvec)
    call d4%get_polarizabilities(mol,gwvec,alpha=a)
    alpha = a(1)
#else
    alpha = 0.0_wp
    write (stdout,'(a)') 'Error: the aISS docking requires CREST compiled with tblite (dftd4)'
    error stop
#endif
  end subroutine free_atom_alpha

!========================================================================================!
end module aiss_iff_ini
