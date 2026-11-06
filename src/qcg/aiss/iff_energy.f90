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
! src/iff/iff_energy.f90), LGPL-3.0-or-later.
!================================================================================!

!> Intermolecular force-field (xTB-IFF) energy of a docked A-B pair
module aiss_iff_energy
  use crest_parameters,only:wp,stdout
  use aiss_param
  implicit none
  private

  public :: iff_e,intermole_probe,alignmol

!========================================================================================!
!========================================================================================!
contains  !> MODULE PROCEDURES START HERE
!========================================================================================!
!========================================================================================!

  subroutine iff_e(d,pr,mode,e,icoord)
    !*****************************************************************
    !* IFF interaction energy e (Eh) of fragment A and fragment B,
    !* the latter rotated/shifted by icoord (translation 1:3, Euler
    !* angles 4:6). Sums Pauli, D3-type dispersion (+ATM for pr,
    !* d%fulle or mode == 0), density-based and LMO electrostatics,
    !* Drude induction, charge transfer, wall and directed-docking
    !* potentials. Thread-safe; called inside OpenMP loops.
    !*****************************************************************
    type(iff_data),intent(in) :: d
    logical,intent(in)   :: pr
    integer,intent(in)   :: mode
    real(wp),intent(out) :: e
    real(wp),intent(in),optional :: icoord(6)

    real(wp) :: ed,es,ep,esl,ei,eabc,ect,ep2,esph,gsolv
    real(wp) :: r2,r6,r8,r42,qq,rmin,dc6,p72,nn12
    real(wp) :: d2(3),t1,t2,t3,c9,ang,damp,tmp2,zz1,zz2,nn1,nn2,ees
    real(wp) :: rotm(3,3),r,rrr,dip(3),rabc(3)
    real(wp),allocatable :: a2(:,:),al2(:,:),q1(:),q2(:)
    real(wp),allocatable :: rab(:,:),r2ab(:,:),r0tmp(:,:),sab(:,:)
    real(wp),allocatable :: rab_solu_solv(:,:),xyz(:,:)
    integer,allocatable  :: pair(:,:)
    integer :: npair,n,n1,n2,i,j,k,kk,iat,jat,metal(86),counter,iw
    logical :: atm,qcgpot

    n = d%n
    n1 = d%n1
    n2 = d%n2
    metal = 1
    metal(1:2) = 0
    metal(6:10) = 0
    metal(14:18) = 0
    metal(32:36) = 0
    metal(50:54) = 0
    metal(82:86) = 0

    allocate (a2(3,n2),al2(4,n2*10),q1(n1),q2(n2))
    allocate (rab(n2,n1),r2ab(n2,n1),r0tmp(n2,n1),sab(n2,n1))
    allocate (pair(2,n1*n2),source=0)

    if (present(icoord)) then
      call rotmat(icoord,rotm)
      call move(n2,d%nlmo2,d%xyz2,d%rlmo2,a2,al2,rotm,icoord)
    else
      a2 = d%xyz2
      al2 = d%rlmo2
    end if

    atm = (pr.or.d%fulle.or.mode .eq. 0)
    qcgpot = (d%ndirected > 0.and.d%directed_type == p_atom_qcg)

    ed = 0
    ep = 0
    es = 0
    esl = 0
    ei = 0
    eabc = 0
    esph = 0
    gsolv = 0
    npair = 0
    counter = 0

    if (qcgpot) allocate (rab_solu_solv(d%ndirected,n2),source=0.0_wp)

!>-- pair distances and induction pair list
    do i = 1,n1
      if (qcgpot) then
        if (any(i == d%directed_atoms)) counter = counter+1
      end if
      do j = 1,n2
        r2 = (d%xyz1(1,i)-a2(1,j))**2+(d%xyz1(2,i)-a2(2,j))**2+(d%xyz1(3,i)-a2(3,j))**2
        r2ab(j,i) = r2+1.d-12
        rab(j,i) = sqrt(r2)+1.d-12
        if (qcgpot) then
          !> solute-solvent distances in Angstrom
          if (any(i == d%directed_atoms)) rab_solu_solv(counter,j) = sqrt(r2)*autoang
        end if
        if (r2 .gt. 12000.0d0.or.r2 .lt. 5.0d0) cycle
        if (metal(d%at1(i)) .eq. 1.or.metal(d%at2(j)) .eq. 1) cycle
        npair = npair+1
        pair(1,npair) = i
        pair(2,npair) = j
      end do
    end do

    call ovlp(n1,n2,r2ab,d%den1,d%den2,sab)

    q1 = d%q1
    q2 = d%q2
    call chrgtransfer(d,n1,n2,d%at1,d%at2,q1,q2,rab,d%debug.and.pr,ect)

!>-- atomic electrostatics: point charges beyond 8 Bohr, else via densities
    p72 = par_es_damp
    do i = 1,n1
      do j = 1,n2
        r = rab(j,i)
        r0tmp(j,i) = p72*d%alpab(j,i)
        if (r .gt. 8) then
          es = es+q1(i)*q2(j)/r
        else
          zz1 = d%z1(i)
          zz2 = d%z2(j)
          nn1 = zz1-q1(i)
          nn2 = zz2-q2(j)
          call truees(n1,n2,i,j,r,r2ab(j,i),d%den1,d%den2,zz1,zz2,nn1,nn2,ees)
          es = es+ees
        end if
      end do
    end do

!>-- D3 dispersion and density-overlap Pauli repulsion
    ep2 = 0
    do i = 1,n1
      iat = d%at1(i)
      nn1 = d%z1(i)-q1(i)
      do j = 1,n2
        r6 = r2ab(j,i)**3
        r8 = r2ab(j,i)*r6
        jat = d%at2(j)
        r42 = rrab(jat,iat)
        dc6 = 1./(r6+r0ab6(jat,iat))+r42/(r8+r0ab8(jat,iat))
        ed = ed+dc6*d%c6ab(j+n1,i)
        if (d%ndirected > 0.and.d%directed_type == p_atom_att) then
          ed = ed-d%directed_val(i)*(exp(-0.01*((rab(j,i)-2)**2)))
        end if
        if (rab(j,i) .gt. 25) cycle
        nn12 = nn1*(d%z2(j)-q2(j))
        ep = ep+nn12*sab(j,i)/rab(j,i)
        if (iat .eq. 1.or.jat .eq. 1) &
        & ep2 = ep2+sqrt(nn12)*exp(-par_xh1*rab(j,i)*d%alpab(j,i))
        if (d%ndirected > 0.and.d%directed_type == p_atom_pot) then
          ep = ep+d%directed_val(i)*(1/rab(j,i)**pot_expo)*pot_scal
        end if
      end do
    end do

!>-- ATM three-body dispersion (no damping on i..k and j..k)
    if (atm) then
      do i = 1,n1
        do j = 1,n2
          r2 = r2ab(j,i)
          if (r2 .gt. 500.0d0) cycle
          d2(1) = r2
          r = rab(j,i)
          rrr = r/d%alpab(j,i)
          damp = 1.0d0/(1.d0+6.d0*rrr**(-8))
          do kk = 1,d%neigh(0,i)
            k = d%neigh(kk,i)
            c9 = d%c6ab(j+n1,i)*d%c6ab(i,k)*d%c6ab(j+n1,k)
            d2(2) = r2ab(j,k)
            d2(3) = sum((d%xyz1(:,i)-d%xyz1(:,k))**2)
            t1 = (d2(1)+d2(2)-d2(3))
            t2 = (d2(1)+d2(3)-d2(2))
            t3 = (d2(3)+d2(2)-d2(1))
            tmp2 = d2(1)*d2(2)*d2(3)
            ang = (0.375d0*t1*t2*t3/tmp2+1.0d0)/tmp2**1.5d0
            eabc = eabc+damp*ang*sqrt(c9)
          end do
          do kk = 1,d%neigh(0,j+n1)
            k = d%neigh(kk,j+n1)
            c9 = d%c6ab(j+n1,i)*d%c6ab(j+n1,k+n1)*d%c6ab(i,k+n1)
            d2(2) = sum((a2(:,j)-a2(:,k))**2)
            d2(3) = r2ab(k,i)
            t1 = (d2(1)+d2(2)-d2(3))
            t2 = (d2(1)+d2(3)-d2(2))
            t3 = (d2(3)+d2(2)-d2(1))
            tmp2 = d2(1)*d2(2)*d2(3)
            ang = (0.375d0*t1*t2*t3/tmp2+1.0d0)/tmp2**1.5d0
            eabc = eabc+damp*ang*sqrt(c9)
          end do
        end do
      end do
    end if

!>-- LMO-LMO and LMO-atom electrostatics
    do i = 1,d%nlmo1
      iat = d%lmo1(i)
      do j = 1,d%nlmo2
        jat = d%lmo2(j)
        r2 = sum((d%rlmo1(1:3,i)-al2(1:3,j))**2)
        r = sqrt(r2)
        esl = esl+d%rlmo1(4,i)*al2(4,j)*r/(r2+r0tmp(jat,iat))
      end do
      do j = 1,n2
        r2 = sum((a2(1:3,j)-d%rlmo1(1:3,i))**2)
        r = sqrt(r2)
        esl = esl+d%rlmo1(4,i)*q2(j)*r/(r2+r0tmp(j,iat))
      end do
    end do
    do i = 1,d%nlmo2
      iat = d%lmo2(i)
      do j = 1,n1
        r2 = sum((d%xyz1(1:3,j)-al2(1:3,i))**2)
        r = sqrt(r2)
        esl = esl+al2(4,i)*q1(j)*r/(r2+r0tmp(iat,j))
      end do
    end do

    call drudescf(n1,n2,npair,pair,d%at1,d%at2,d%xyz1,q1,d%qdr1,a2,q2,d%qdr2,ei)

!>-- wall potentials
    esph = 0.0_wp
    if (d%qcg.and.d%nwalls > 0) then
      rabc = -1.0_wp
      do iw = 1,d%nwalls
        if (.not.allocated(d%walls(iw)%list)) then
          rabc(1:3) = d%walls(iw)%radius(1:3)
          exit
        end if
      end do
      call cavitye(n1,d%xyz1,rabc,esph)
      call cavitye(n2,a2,rabc,esph)
    else if (d%nwalls > 0) then
      allocate (xyz(3,n))
      xyz(1:3,1:n1) = d%xyz1(1:3,1:n1)
      xyz(1:3,n1+1:n) = a2(1:3,1:n2)
      call cavity_egrad(n,xyz,d%nwalls,d%walls,esph)
    end if

    ep = ep*par_rep_scal+ep2*par_xh2*0.01
    e = -ed+es+ep+esl+ei+gsolv+eabc+ect+esph

    !> QCG: erf-shaped attraction by the closest solute-solvent contact,
    !> (a/2)*erf(-(b*r-c*b))+a/2 with a = val(1), b = expo(1), c = expo(2)
    if (qcgpot) then
      e = e-((d%directed_val(1)/2) &
      &   *erf(-(d%directed_expo(1)*minval(rab_solu_solv))+(d%directed_expo(2)*d%directed_expo(1))) &
      &   +(d%directed_val(1)/2))
    end if

    if (isnan(e)) e = 1.0e42_wp

    if (pr.and.e .lt. 10000.2) then
      call rcma(n1,d%xyz1,d%at1,n2,a2,d%at2,r,rmin)
      write (stdout,'(''R CMA (Angst) :'',F10.3)') r*autoang
      write (stdout,'(''R min (Angst) :'',F10.3)') rmin*autoang
      call diptot(n1,n2,d%nlmo1,d%nlmo2,d%xyz1,a2,d%rlmo1,al2,q1,q2,dip)
      write (stdout,'(''dipole moment :'',F10.3)') sqrt(dip(1)**2+dip(2)**2+dip(3)**2)
      write (stdout,'(''intermolecular energies in kcal/mol'')')
      write (stdout,'(''E Pauli       :'',F10.3)') ep*au
      write (stdout,'(''E disp ATM    :'',F10.3)') eabc*au
      write (stdout,'(''E disp 2B     :'',F10.3)')-ed*au
      write (stdout,'(''E disp total  :'',F10.3)')-ed*au+eabc*au
      write (stdout,'(''E ES atom     :'',F10.3)') es*au
      write (stdout,'(''E ES LMO      :'',F10.3)') esl*au
      write (stdout,'(''E ES total    :'',F10.3)') (es+esl)*au
      write (stdout,'(''E induction   :'',F10.3)') ei*au
      write (stdout,'(''E CT          :'',F10.3)') ect*au
      if (d%nwalls > 0) write (stdout,'(''E cavity      :'',F10.3)') esph*au
      write (stdout,'(''Eint total,gas:'',F10.3)') (e-gsolv)*au
      write (stdout,'(''              '',F14.8,''  <== Gint total'')') e*au
    end if
  end subroutine iff_e

!========================================================================================!

  subroutine drudescf(n1,n2,npair,pair,at1,at2,xyz1,qat1,qdr1,xyz2,qat2,qdr2,edr)
    !*****************************************************************
    !* Drude induction energy: one steepest-descent step (alpha = 2)
    !* of the Drude positions, energy relative to the start
    !*****************************************************************
    integer,intent(in)   :: n1,n2,npair,pair(2,n1*n2),at1(n1),at2(n2)
    real(wp),intent(in)  :: xyz1(3,n1),xyz2(3,n2)
    real(wp),intent(in)  :: qdr1(n1),qat1(n1),qdr2(n2),qat2(n2)
    real(wp),intent(out) :: edr
    real(wp),parameter :: alpha = 2.0_wp
    real(wp) :: e0,q1(n1),q2(n2),gdr1(3,n1),gdr2(3,n2),xyzdr1(3,n1),xyzdr2(3,n2)
    integer :: i,ii

    xyzdr1 = xyz1
    xyzdr2 = xyz2
    do i = 1,n1
      q1(i) = qat1(i)-qdr1(i)
    end do
    do i = 1,n2
      q2(i) = qat2(i)-qdr2(i)
    end do
    e0 = 0.0_wp
    do ii = 1,2
      call indeg(n1,n2,npair,pair,xyz1,q1,xyzdr1,qdr1,xyz2,q2,xyzdr2,qdr2, &
      &          at1,at2,ii,edr,gdr1,gdr2)
      if (ii .eq. 1) then
        e0 = edr
        do i = 1,n1
          xyzdr1(1:3,i) = xyzdr1(1:3,i)-alpha*gdr1(1:3,i)
        end do
        do i = 1,n2
          xyzdr2(1:3,i) = xyzdr2(1:3,i)-alpha*gdr2(1:3,i)
        end do
      end if
      edr = edr-e0
    end do
  end subroutine drudescf

  subroutine indeg(n1,n2,npair,pair,xyz1,q1,xyzdr1,qdr1,xyz2,q2,xyzdr2,qdr2, &
  &                at1,at2,kode,ec,gdr1,gdr2)
    !*****************************************************************
    !* Drude spring + Coulomb energy ec; kode 1 adds the gradient
    !* w.r.t. the Drude positions, kode 2 energy only
    !*****************************************************************
    integer,intent(in)   :: n1,n2,npair,pair(2,n1*n2),kode,at1(n1),at2(n2)
    real(wp),intent(in)  :: xyz1(3,n1),xyzdr1(3,n1),q1(n1),qdr1(n1)
    real(wp),intent(in)  :: xyz2(3,n2),xyzdr2(3,n2),q2(n2),qdr2(n2)
    real(wp),intent(out) :: gdr1(3,n1),gdr2(3,n2),ec
    integer  :: i,i1,i2,ij
    real(wp) :: rab,rabd,dum,dum2,oner,dx,dy,dz,edr,espring,qq,damp
    logical  :: grad

    grad = (kode .eq. 1)
    damp = par_drude_damp
    if (grad) then
      gdr1 = 0.0_wp
      gdr2 = 0.0_wp
    end if

    espring = 0
    do i = 1,n1
      dx = xyzdr1(1,i)-xyz1(1,i)
      dy = xyzdr1(2,i)-xyz1(2,i)
      dz = xyzdr1(3,i)-xyz1(3,i)
      rab = dx*dx+dy*dy+dz*dz
      dum = par_drude_fc*gam(at1(i))
      espring = espring+dum*rab
      if (grad) then
        dum2 = 2.*dum
        gdr1(1,i) = dum2*dx
        gdr1(2,i) = dum2*dy
        gdr1(3,i) = dum2*dz
      end if
    end do
    do i = 1,n2
      dx = xyzdr2(1,i)-xyz2(1,i)
      dy = xyzdr2(2,i)-xyz2(2,i)
      dz = xyzdr2(3,i)-xyz2(3,i)
      rab = dx*dx+dy*dy+dz*dz
      dum = par_drude_fc*gam(at2(i))
      espring = espring+dum*rab
      if (grad) then
        dum2 = 2.*dum
        gdr2(1,i) = dum2*dx
        gdr2(2,i) = dum2*dy
        gdr2(3,i) = dum2*dz
      end if
    end do

    edr = 0
    do ij = 1,npair
      i1 = pair(1,ij)
      i2 = pair(2,ij)
!>-- Drude A - Drude B
      dx = xyzdr1(1,i1)-xyzdr2(1,i2)
      dy = xyzdr1(2,i1)-xyzdr2(2,i2)
      dz = xyzdr1(3,i1)-xyzdr2(3,i2)
      rab = sqrt(dx*dx+dy*dy+dz*dz)
      rabd = rab+damp
      qq = qdr1(i1)*qdr2(i2)
      edr = edr+qq/rabd
      if (grad) then
        oner = qq/(rabd*rabd*rab)
        gdr1(1,i1) = gdr1(1,i1)-dx*oner
        gdr1(2,i1) = gdr1(2,i1)-dy*oner
        gdr1(3,i1) = gdr1(3,i1)-dz*oner
        gdr2(1,i2) = gdr2(1,i2)+dx*oner
        gdr2(2,i2) = gdr2(2,i2)+dy*oner
        gdr2(3,i2) = gdr2(3,i2)+dz*oner
      end if
!>-- atom A - Drude B
      dx = xyz1(1,i1)-xyzdr2(1,i2)
      dy = xyz1(2,i1)-xyzdr2(2,i2)
      dz = xyz1(3,i1)-xyzdr2(3,i2)
      rab = sqrt(dx*dx+dy*dy+dz*dz)
      rabd = rab+damp
      qq = q1(i1)*qdr2(i2)
      edr = edr+qq/rabd
      if (grad) then
        oner = qq/(rabd*rabd*rab)
        gdr2(1,i2) = gdr2(1,i2)+dx*oner
        gdr2(2,i2) = gdr2(2,i2)+dy*oner
        gdr2(3,i2) = gdr2(3,i2)+dz*oner
      end if
!>-- atom B - Drude A
      dx = xyz2(1,i2)-xyzdr1(1,i1)
      dy = xyz2(2,i2)-xyzdr1(2,i1)
      dz = xyz2(3,i2)-xyzdr1(3,i1)
      rab = sqrt(dx*dx+dy*dy+dz*dz)
      rabd = rab+damp
      qq = q2(i2)*qdr1(i1)
      edr = edr+qq/rabd
      if (grad) then
        oner = qq/(rabd*rabd*rab)
        gdr1(1,i1) = gdr1(1,i1)+dx*oner
        gdr1(2,i1) = gdr1(2,i1)+dy*oner
        gdr1(3,i1) = gdr1(3,i1)+dz*oner
      end if
    end do

    ec = edr+espring
  end subroutine indeg

!========================================================================================!

  subroutine truees(m1,m2,i,j,r,r2,den1,den2,z1,z2,n1,n2,es)
    !*****************************************************************
    !* Electrostatics of atoms i (A) and j (B) at distance r with
    !* spherical 4-primitive Gaussian valence densities; z = valence
    !* charges, n = valence electron numbers (q = z-n)
    !*****************************************************************
    integer,intent(in)   :: i,j,m1,m2
    real(wp),intent(in)  :: r,r2,z1,z2,n1,n2
    real(wp),intent(in)  :: den1(2,4,m1),den2(2,4,m2)
    real(wp),intent(out) :: es
    real(wp),parameter :: tpi = 6.28318530717958_wp
    real(wp),parameter :: tpi25 = 34.98683665524963_wp
    real(wp) :: cc1,cc2,apb1,apb2,t,f02,vv,f1,v1,v2,abcd,eabcd
    integer  :: ip,jp,kp,lp
    integer,parameter :: np = 4

!>-- nucleus-nucleus and electron-nucleus
    es = z1*z2/r
    v1 = 0
    v2 = 0
    do ip = 1,np
      do jp = 1,ip-1
        apb1 = den1(1,ip,i)+den1(1,jp,i)
        cc1 = den1(2,ip,i)*den1(2,jp,i)
        apb2 = den2(1,ip,j)+den2(1,jp,j)
        cc2 = den2(2,ip,j)*den2(2,jp,j)
        call get_f02(apb1*r2,f02)
        v1 = v1+f02*cc1*tpi/apb1
        call get_f02(apb2*r2,f02)
        v2 = v2+f02*cc2*tpi/apb2
      end do
      apb1 = den1(1,ip,i)+den1(1,ip,i)
      cc1 = den1(2,ip,i)*den1(2,ip,i)
      apb2 = den2(1,ip,j)+den2(1,ip,j)
      cc2 = den2(2,ip,j)*den2(2,ip,j)
      call get_f02(apb1*r2,f02)
      v1 = v1+0.5*f02*cc1*tpi/apb1
      call get_f02(apb2*r2,f02)
      v2 = v2+0.5*f02*cc2*tpi/apb2
    end do
    v1 = 2.0d0*v1
    v2 = 2.0d0*v2
    es = es-v1*z2*n1-v2*z1*n2

!>-- electron-electron
    vv = 0
    do ip = 1,np
      do jp = 1,ip
        f1 = 2.0d0
        if (ip .eq. jp) f1 = 1.0d0
        apb1 = den1(1,ip,i)+den1(1,jp,i)
        cc1 = den1(2,ip,i)*den1(2,jp,i)*f1
        do kp = 1,np
          do lp = 1,kp-1
            apb2 = den2(1,kp,j)+den2(1,lp,j)
            cc2 = den2(2,kp,j)*den2(2,lp,j)*2.0d0
            eabcd = apb1*apb2
            abcd = apb1+apb2
            t = r2*eabcd/abcd
            call get_f02(t,f02)
            vv = vv+tpi25/(eabcd*dsqrt(abcd))*f02*cc1*cc2
          end do
          apb2 = den2(1,kp,j)+den2(1,kp,j)
          cc2 = den2(2,kp,j)*den2(2,kp,j)
          eabcd = apb1*apb2
          abcd = apb1+apb2
          t = r2*eabcd/abcd
          call get_f02(t,f02)
          vv = vv+tpi25/(eabcd*dsqrt(abcd))*f02*cc1*cc2
        end do
      end do
    end do
    es = es+n1*n2*vv
  end subroutine truees

  subroutine ovlp(n1,n2,r2ab,den1,den2,sab)
    !> overlap of the Gaussian atomic densities (zero beyond r2 = 500)
    integer,intent(in)   :: n1,n2
    real(wp),intent(in)  :: r2ab(n2,n1),den1(2,4,n1),den2(2,4,n2)
    real(wp),intent(out) :: sab(n2,n1)
    real(wp) :: apb,s,e1,e2,c1,c2
    integer  :: i,j,ip,jp
    integer,parameter :: np = 4
    sab = 0
    do i = 1,n1
      do j = 1,n2
        if (r2ab(j,i) .gt. 500.) cycle
        s = 0
        do ip = 1,np
          e1 = den1(1,ip,i)
          c1 = den1(2,ip,i)
          do jp = 1,np
            e2 = den2(1,jp,j)
            c2 = den2(2,jp,j)
            apb = 1.0d0/(e1+e2)
            s = s+c1*c2*(pi*apb)**1.50d0*dexp(-r2ab(j,i)*e1*e2*apb)
          end do
        end do
        sab(j,i) = s
      end do
    end do
  end subroutine ovlp

!========================================================================================!

  subroutine chrgtransfer(d,n1,n2,at1,at2,q1,q2,rab,pr,ect)
    !*****************************************************************
    !* Charge-transfer correction from a 3x3 CI of the ground state
    !* and the two HOMO->LUMO CT configurations (A->B, B->A); q1/q2
    !* are replaced by the CT-corrected charges
    !*****************************************************************
    type(iff_data),intent(in) :: d
    integer,intent(in)     :: n1,n2,at1(n1),at2(n2)
    real(wp),intent(inout) :: q1(n1),q2(n2)
    real(wp),intent(in)    :: rab(n2,n1)
    logical,intent(in)     :: pr
    real(wp),intent(out)   :: ect
    !> orbital gap shift a la sTDA-xTB, R and gap coupling decay, off-diagonal scaling
    real(wp),parameter :: shift = 0.37_wp,expo = 2.9_wp,scal = 0.02_wp
    real(wp),parameter :: expo2 = 100.0_wp,gfac = 1.0_wp
    real(wp) :: r,a(6),u(3,3),ee(3),coup1,coup2,u0,u12,u21,ex,ex1,ex2
    real(wp) :: q01(n1),q02(n2),dq1,dq2,e,e0,emo
    integer  :: i,j

    q01 = q1
    q02 = q2
    coup1 = 0
    coup2 = 0
    do i = 1,n1
      do j = 1,n2
        r = rab(j,i)
        if (r .gt. 20.0d0) cycle
        ex = exp(-expo*(r/d%alpab(j,i))**0.5d0)
        coup1 = coup1+d%qct1(i,1)*d%qct2(j,2)*ex
        coup2 = coup2+d%qct2(j,1)*d%qct1(i,2)*ex
      end do
    end do

    a(1) = 0.0d0
    ex1 = exp(-expo2*(d%elumo(2)-d%ehomo(1))**2)
    if (d%elumo(2)-d%ehomo(1) .lt. 0) ex1 = 1.
    ex2 = exp(-expo2*(d%elumo(1)-d%ehomo(2))**2)
    if (d%elumo(1)-d%ehomo(2) .lt. 0) ex2 = 1.
    a(2) = -coup1*scal*100.*ex1
    a(4) = -coup2*scal*100.*ex2
    a(3) = d%elumo(2)-d%ehomo(1)+shift
    a(5) = 0.0d0
    a(6) = d%elumo(1)-d%ehomo(2)+shift
    call rsp_packed(a,3,ee,u)
    i = 1
    u0 = u(1,1)**2
    if (u(1,2)**2 .gt. u0) i = 2
    if (u(1,3)**2 .gt. u0) i = 3
    u0 = u(1,i)**2
    u12 = u(2,i)**2
    u21 = u(3,i)**2

    dq1 = 0
    do i = 1,n1
      q1(i) = u0*q1(i)+u12*(q1(i)+d%qct1(i,1))+u21*(q1(i)-d%qct1(i,2))
      dq1 = dq1+q1(i)-q01(i)
    end do
    call esxtb(n1,at1,d%gab1,q1,e)
    call esxtb(n1,at1,d%gab1,q01,e0)
    ect = (e-e0)*gfac
    if (pr) then
      write (stdout,*) ' dq system 1       ',dq1
      write (stdout,*) ' De system 1 /kcal ',627.51*(e-e0)
    end if

    dq2 = 0
    do i = 1,n2
      q2(i) = u0*q2(i)+u12*(q2(i)-d%qct2(i,2))+u21*(q2(i)+d%qct2(i,1))
      dq2 = dq2+q2(i)-q02(i)
    end do
    call esxtb(n2,at2,d%gab2,q2,e)
    call esxtb(n2,at2,d%gab2,q02,e0)
    ect = ect+(e-e0)*gfac

    if (dq1 .gt. 0) then
      emo = -dq1*d%ehomo(1)+dq1*d%elumo(2)
    else
      emo = dq1*d%ehomo(2)-dq1*d%elumo(1)
    end if
    ect = ect+emo

    if (pr) then
      write (stdout,*) ' dq system 2       ',dq2
      write (stdout,*) ' De system 2 /kcal ',627.51*(e-e0)
      write (stdout,*) ' De MO       /kcal ',627.51*emo
    end if
  end subroutine chrgtransfer

  subroutine esxtb(n,at,gab,q,es)
    !> xTB-like intramolecular second- and third-order electrostatics
    integer,intent(in)   :: n,at(n)
    real(wp),intent(in)  :: q(n),gab(n,n)
    real(wp),intent(out) :: es
    !> PBE values
    real(wp),parameter :: gam3(94) = [ &
    &-0.02448_wp,0.178614_wp,0.194034_wp,0.154068_wp,0.173892_wp,0.16716_wp,0.156306_wp, &
    &0.321466_wp,0.163314_wp,0.170862_wp,0.256128_wp,0.18906_wp,0.146310_wp,0.136686_wp, &
    &0.123558_wp,0.122070_wp,0.119424_wp,0.115368_wp,0.175938_wp,0.152214_wp,0.24030_wp, &
    &0.204654_wp,0.186552_wp,0.178824_wp,0.174474_wp,0.205038_wp,0.188772_wp,0.180462_wp, &
    &0.180354_wp,0.176508_wp,0.158250_wp,0.139212_wp,0.131868_wp,0.119712_wp,0.115476_wp, &
    &0.108444_wp,0.091032_wp,0.07698_wp,0.102642_wp,0.095346_wp,0.088266_wp,0.086364_wp, &
    &0.085254_wp,0.088242_wp,0.087774_wp,0.088470_wp,0.091314_wp,0.090372_wp,0.110862_wp, &
    &0.093588_wp,0.079908_wp,0.074082_wp,0.069342_wp,0.064638_wp,0.077826_wp,0.0590214_wp, &
    &0.073614_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp, &
    &0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp, &
    &0.0_wp,0.085236_wp,0.07782_wp,0.074112_wp,0.07206_wp,0.070188_wp,0.069246_wp, &
    &0.072042_wp,0.07038_wp,0.072570_wp,0.108096_wp,0.084150_wp,0.070350_wp,0.06954_wp, &
    &0.064866_wp,0.0596874_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp, &
    &0.0_wp,0.0_wp,0.0_wp]
    real(wp) :: t,qq
    integer  :: i,j
    es = 0
    do i = 1,n
      do j = i+1,n
        es = es+q(i)*q(j)*gab(j,i)
      end do
    end do
    es = es*2.0d0
    t = 0
    do i = 1,n
      qq = q(i)*q(i)
      t = t+gam3(at(i))*q(i)*qq
      es = es+qq*gab(i,i)
    end do
    es = es*0.5d0+t/3.0d0
  end subroutine esxtb

  pure subroutine get_f02(arg,f02)
    !> Boys function F0
    real(wp),intent(in)  :: arg
    real(wp),intent(out) :: f02
    if (arg .ge. 1.0d-10) then
      f02 = 0.50d0*dsqrt(pi/arg)*derf(dsqrt(arg))
    else
      f02 = 1.0d0-0.33333333333333d0*arg
    end if
  end subroutine get_f02

!========================================================================================!

  subroutine intermole_probe(d,c2,e,eqp)
    !*****************************************************************
    !* Rare-gas probe atom at c2 around fragment A: neutral probe
    !* energy e (Pauli + scaled dispersion) and the electrostatic
    !* energy eqp of a 0.1 point charge at that position
    !*****************************************************************
    type(iff_data),intent(in) :: d
    real(wp),intent(in)  :: c2(3)
    real(wp),intent(out) :: e,eqp
    real(wp),allocatable :: rab_solu_solv(:)
    integer  :: i,iat,counter
    real(wp) :: ed,ep,r,r2,rr,r6,r42,r06,r08,dc6
    logical  :: qcgpot

    ed = 0
    ep = 0
    eqp = 0
    counter = 0
    qcgpot = (d%ndirected > 0.and.d%directed_type == p_atom_qcg)
    if (qcgpot) allocate (rab_solu_solv(d%ndirected),source=huge(1.0_wp))

    do i = 1,d%n1
      iat = d%at1(i)
      r2 = sum((d%xyz1(1:3,i)-c2(1:3))**2)
      if (r2 .gt. 5000.0d0) cycle
      r = sqrt(r2)
      eqp = eqp+d%q1(i)/(r+0.1)
      r6 = r2*r2*r2
      r42 = rrab(iat,probe_atom_type)
      r06 = r0ab6(iat,probe_atom_type)
      r08 = r0ab8(iat,probe_atom_type)
      if (qcgpot) then
        if (any(i == d%directed_atoms)) then
          counter = counter+1
          rab_solu_solv(counter) = r*autoang
        end if
      end if
      dc6 = 1./(r6+r06)+r42/(r6*r2+r08)
      ed = ed+dc6*d%cprob(i)
      if (d%ndirected > 0.and.d%directed_type == p_atom_att) then
        ed = ed-d%directed_val(i)*(exp(-0.1*((r-2)**2)))
      end if
      !> Pauli fitted to the full IFF potential
      rr = -3.85*sqrt(r2/r42)
      ep = ep+val_e(iat)*(erf(rr)+1.)/(r+0.6)
      if (d%ndirected > 0.and.d%directed_type == p_atom_pot) then
        ep = ep+d%directed_val(i)*(1/r**pot_expo)*pot_scal
      end if
    end do

    do i = 1,d%nlmo1
      r2 = sum((c2(1:3)-d%rlmo1(1:3,i))**2)
      if (r2 .gt. 1000.0d0) cycle
      r = sqrt(r2)
      eqp = eqp+d%rlmo1(4,i)/(r+0.1)
    end do

    e = 8.0d0*ep-ed*0.70d0
    !> same erf-shaped QCG attraction as in iff_e
    if (qcgpot) then
      e = e-((d%directed_val(1)/2) &
      &   *erf(-((d%directed_expo(1)*minval(rab_solu_solv))-(d%directed_expo(2)*d%directed_expo(1)))) &
      &   +(d%directed_val(1)/2))
    end if

    eqp = eqp*0.1
  end subroutine intermole_probe

!========================================================================================!

  subroutine alignmol(n1,n2,at1,at2,xyz1,xyz2,coord)
    !*****************************************************************
    !* Place molecule B (xyz2) into a grid-point cluster (xyz1) by a
    !* genetic minimization of a crude dispersion-like potential;
    !* coord returns the best translation/rotation
    !*****************************************************************
    integer,intent(in)   :: n1,n2,at1(n1),at2(n2)
    real(wp),intent(in)  :: xyz1(3,n1),xyz2(3,n2)
    real(wp),intent(out) :: coord(6)
    integer,parameter :: maxgen = 20,maxparent = 50,ndim = maxparent**2
    real(wp),allocatable :: found(:,:),found2(:,:)
    real(wp) :: rr,av,sig,displ(6),e,r1(3),c2(3,n2),c6i
    real :: f(6)
    integer  :: i,j,ii,jj,kk,icycle

    allocate (found(7,ndim),found2(7,ndim),source=0.0_wp)
    coord = 0

    call cmadock(n1,n1,at1,xyz1,r1)
    displ = 0
    do i = 1,maxparent
      displ(1:3) = r1(1:3)
      call rand6(0.8d0,0.10d0,6.28d0,displ)
      found(1:6,i) = displ(1:6)
    end do

    do icycle = 1,maxgen
      kk = 0
      do i = 1,maxparent
        do j = 1,maxparent
          kk = kk+1
          call random_number(f)
          displ(1:6) = found(1:6,i)*f(1:6)+found(1:6,j)*(1.0d0-f(1:6))
          if (i .ne. j) call rand6(0.8d0,1.00d0,0.10d0,displ)
          call move2(n2,xyz2,c2,displ)
          e = 0
          c6i = float(at1(1))
          do ii = 1,n1
            do jj = 1,n2
              rr = sum((xyz1(1:3,ii)-c2(1:3,jj))**2)
              e = e-c6i*float(at2(jj))/(rr*rr*rr+1.d+5)
            end do
          end do
          found2(1:6,kk) = displ(1:6)
          found2(7,kk) = e
        end do
      end do
      call sort6(kk,found2)
      found = found2
      av = sum(found(7,1:maxparent))/float(maxparent)
      sig = 0
      do i = 1,maxparent
        sig = sig+(found(7,i)-av)**2
      end do
      sig = sqrt(sig/float(maxparent-1))
      if (sig .lt. 0.001d0) exit
    end do

    coord(1:6) = found(1:6,1)
    deallocate (found,found2)
  end subroutine alignmol

!========================================================================================!
end module aiss_iff_energy
