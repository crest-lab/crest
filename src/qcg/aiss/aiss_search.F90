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
! Ported from the xtb docking submodule (github.com/grimme-lab/xtb,
! src/docking/search_nci.f90), LGPL-3.0-or-later.
!================================================================================!

!> aISS screening (RG grid, pocket, stack, angular, genetic) and final optimizations
module aiss_search
  use crest_parameters,only:wp,stdout
  use strucrd,only:coord,i2e
  use crest_calculator,only:calcdata,engrad
  use optimize_module,only:optimize_geometry
  use aiss_param
  use aiss_iff_energy,only:iff_e,intermole_probe,alignmol
  implicit none
  private

  public :: aiss_docking_search

!========================================================================================!
!========================================================================================!
contains  !> MODULE PROCEDURES START HERE
!========================================================================================!
!========================================================================================!

  subroutine aiss_docking_search(cfg,d,molA,molB,calc,molA_e,molB_e,iunit, &
  &                              best,eint,io)
    !******************************************************************
    !* Global search for the best placement of molecule B around A.
    !* Genes are (R_x,R_y,R_z,alpha,beta,gamma) of B relative to A.
    !* The n_opt best genes (plus pocket structures) are optimized
    !* with calc, and the lowest structure is returned.
    !*
    !* cfg            : run settings (local copy, modified)
    !* d              : initialized IFF data (CMA frames of A and B)
    !* molA, molB     : CMA-shifted fragments
    !* calc           : calculator for the complex (with constraints)
    !* molA_e,molB_e  : fragment reference energies (Eh)
    !* iunit          : output unit, 0 for silent runs
    !* best           : lowest optimized complex, best%energy in Eh
    !* eint           : lowest interaction energy (kcal/mol)
    !* io             : 0 on success
    !******************************************************************
    type(aiss_config),intent(inout) :: cfg
    type(iff_data),intent(in)       :: d
    type(coord),intent(in)          :: molA,molB
    type(calcdata),intent(inout)    :: calc
    real(wp),intent(in)             :: molA_e,molB_e
    integer,intent(in)              :: iunit
    type(coord),intent(out)         :: best
    real(wp),intent(out)            :: eint
    integer,intent(out)             :: io

    integer :: n,n1,n2
    real(wp),allocatable :: found(:,:),found2(:,:)
    real(wp),allocatable :: xyzprobc(:,:),xyzprob3(:,:),xyzprob4(:,:)
    real(wp),allocatable :: xyztmp(:,:),rprobc(:,:)
    real(wp),allocatable :: xyzprob0(:),xyzprob1(:),xyzprob2(:),rprob(:)
    integer,allocatable  :: ind0(:),ind1(:),ind2(:),rind(:),attmp(:)
    integer,allocatable  :: atfrag(:),molvec(:),iatnum(:)
    real(wp),allocatable :: xyzfrag(:,:)
    real(wp),allocatable :: xyz_pocket(:,:,:),pocket_e(:)
    real(wp),allocatable :: xyz_opt(:,:,:),final_e(:)
    real(wp),allocatable :: grad(:,:)
    real(wp) :: found_tmp(6),displ(6),f(6)
    real(wp) :: r,e,epq,dx,dy,dz,xx1,yy1,zz1,bord,bord_rg,stepr3(3),stepr4
    real(wp) :: dx_rg,dy_rg,dz_rg,xx1_rg,yy1_rg,zz1_rg,dx2,dy2,dz2
    real(wp) :: r1,r2,av,sig,etot,pmass,rmin
    real(wp) :: dum3(3,3),dum2(3),dum(3),mvec(3,14)
    real(wp) :: xyz_tmp(3,molB%nat)
    integer :: i,j,k,l,ii,ll,m1,m2,m3,mm1,mm2,mm3,nn,ndim,nclust,icycle
    integer :: kk,mfrag,maxlow,maxpock,icssym,sig_count,iostat
    logical :: pr
    type(coord) :: comb,molopt
    data mvec/1.,0.,0.,0.,1.,0.,0.,0.,1., &
    & -1.,0.,0.,0,-1.,0.,0.,0.,-1., &
    & 1.,1.,1.,-1.,-1.,-1., &
    & -1.,1.,1.,1.,-1.,1.,1.,1.,-1., &
    & 1.,-1.,-1.,-1.,1.,-1.,-1.,-1.,1./

    io = 0
    eint = 0.0_wp
    pr = (iunit /= 0)
    n = d%n
    n1 = d%n1
    n2 = d%n2
    kk = 0
    nclust = 0
    ll = 0

    if (pr) then
      write (iunit,*)
      write (iunit,*) '=============================================='
      write (iunit,*) '|         Starting Energy Screening          |'
      write (iunit,*) '=============================================='
      write (iunit,*)
      if (.not.cfg%fulle) then
        write (iunit,*) ' No ATM term employed (recommended)'
        write (iunit,*)
      end if
    end if

    call axis_docking(n1,d%at1,d%xyz1,dum,pmass,dum2,dum3)
    r1 = sqrt(dum2(1)+dum2(2)+dum2(3))
    call axis_docking(n2,d%at2,d%xyz2,dum,pmass,dum2,dum3)
    r2 = sqrt(dum2(1)+dum2(2)+dum2(3))
    call rcma(n1,d%xyz1,d%at1,n2,d%xyz2,d%at2,r,rmin)
    if (pr) then
      write (iunit,'('' Method for final opts.    :'',1x,a  )') cfg%optlvl
      write (iunit,'('' # of genetic optimizations:'',1x,i0  )') cfg%maxgen
      write (iunit,'('' # of parents              :'',1x,i0  )') cfg%maxparent
      write (iunit,'('' # of final geo. opts.     :'',1x,i0  )') cfg%n_opt
      write (iunit,'('' Rare gas grid step size   :'',F8.2)') cfg%stepr
      write (iunit,'('' ang step size /deg        :'',F8.2)') cfg%stepa
      write (iunit,'('' # angular grid points     :'',1x,i0  )') (360/int(cfg%stepa))**3
      write (iunit,*)
      if (cfg%pocket_grid) write (iunit,'('' Performing pocket search'')')
      if (cfg%stack_grid) write (iunit,'('' Performing stack search'')')
      if (cfg%angular_grid) write (iunit,'('' Performing angular search'')')
      write (iunit,*)
      write (iunit,'(''  Total '',a,'' energy molecule 1: '',F16.10)') cfg%optlvl,molA_e
      write (iunit,'(''  Total '',a,'' energy molecule 2: '',F16.10)') cfg%optlvl,molB_e
      write (iunit,*)
      write (iunit,*) '-----------------------------'
      write (iunit,*) ' Grid based energy screening '
      write (iunit,*) '-----------------------------'
      write (iunit,*)
    end if

    comb = molA
    comb%nat = n
    deallocate (comb%at,comb%xyz)
    allocate (comb%at(n),comb%xyz(3,n),grad(3,n))
    comb%at(1:n1) = molA%at
    comb%at(n1+1:n) = molB%at
    comb%chrg = molA%chrg+molB%chrg
    comb%uhf = molA%uhf+molB%uhf

!>-- grid borders
    xx1 = maxval(d%xyz1(1,1:n1))-minval(d%xyz1(1,1:n1))
    yy1 = maxval(d%xyz1(2,1:n1))-minval(d%xyz1(2,1:n1))
    zz1 = maxval(d%xyz1(3,1:n1))-minval(d%xyz1(3,1:n1))
    dum = [xx1,yy1,zz1]
    dx2 = maxval(d%xyz2(1,1:n2))-minval(d%xyz2(1,1:n2))
    dy2 = maxval(d%xyz2(2,1:n2))-minval(d%xyz2(2,1:n2))
    dz2 = maxval(d%xyz2(3,1:n2))-minval(d%xyz2(3,1:n2))
    dum2 = [dx2,dy2,dz2]
    bord = maxval(abs(dum))+maxval(dum2)+2
    bord_rg = 30

!>-- crude R grid (4x4x4) for the angular search
    dx = maxval(d%xyz1(1,1:n1))-minval(d%xyz1(1,1:n1))+bord
    dy = maxval(d%xyz1(2,1:n1))-minval(d%xyz1(2,1:n1))+bord
    dz = maxval(d%xyz1(3,1:n1))-minval(d%xyz1(3,1:n1))+bord
    xx1 = minval(d%xyz1(1,1:n1))-bord*0.5
    yy1 = minval(d%xyz1(2,1:n1))-bord*0.5
    zz1 = minval(d%xyz1(3,1:n1))-bord*0.5
    stepr3(1) = dx/3
    stepr3(2) = dy/3
    stepr3(3) = dz/3
    mm1 = 4
    mm2 = 4
    mm3 = 4

!>-- enlarged rare-gas grid (cheap, gives enough points for the pocket search)
    dx_rg = dx+bord_rg
    dy_rg = dy+bord_rg
    dz_rg = dz+bord_rg
    xx1_rg = xx1-bord_rg*0.5
    yy1_rg = yy1-bord_rg*0.5
    zz1_rg = zz1-bord_rg*0.5
    m1 = 1+int((dx_rg)/cfg%stepr)
    m2 = 1+int((dy_rg)/cfg%stepr)
    m3 = 1+int((dz_rg)/cfg%stepr)

    icssym = 0
    if (cfg%cssym) then
      dx2 = dx-bord
      dy2 = dy-bord
      dz2 = dz-bord
      if (dx2 .lt. dy2.and.dx2 .lt. dz2) then
        m1 = m1/2
        mm1 = mm1/2
        xx1 = xx1+bord*0.5
        icssym = 1
      end if
      if (dy2 .lt. dx2.and.dy2 .lt. dz2) then
        m2 = m2/2
        mm2 = mm2/2
        yy1 = yy1+bord*0.5
        icssym = 2
      end if
      if (dz2 .lt. dx2.and.dz2 .lt. dy2) then
        m3 = m3/2
        mm3 = mm3/2
        zz1 = zz1+bord*0.5
        icssym = 3
      end if
      if (pr) write (iunit,'(''  # Cs symmetry red. :'',i1  )') icssym
    end if

    nn = m1*m2*m3
    if (cfg%maxparent < 28) then
      !> fewer parents than genes would crash the genetic step
      cfg%maxparent = 28
      if (pr) write (iunit,*) 'Warning: Too few maxparent, taking 28'
    end if
    maxpock = 10
    ndim = cfg%maxparent**2

    allocate (found(7,ndim),found2(7,ndim),source=0.0_wp)

    if (pr) write (iunit,'(''  # probe RG points   :'',i0  )') nn
    allocate (xyzprob0(nn),xyzprob1(nn),xyzprob2(nn),xyzprobc(3,nn),source=0.0_wp)
    allocate (ind0(nn),ind1(nn),ind2(nn))

!>-- rare-gas probe grid around A
    dx_rg = xx1_rg
    l = 0
    do i = 1,m1
      dum(1) = dx_rg
      dy_rg = yy1_rg
      do j = 1,m2
        dum(2) = dy_rg
        dz_rg = zz1_rg
        do k = 1,m3
          dum(3) = dz_rg
          l = l+1
          xyzprobc(1:3,l) = dum(1:3)
          dz_rg = dz_rg+cfg%stepr
        end do
        dy_rg = dy_rg+cfg%stepr
      end do
      dx_rg = dx_rg+cfg%stepr
    end do
    if (cfg%debug) call wr_grid('RG_grid.xyz',molA%nat,nn,molA%at,36,molA%xyz,xyzprobc)

    !$omp parallel do default(none) &
    !$omp shared(xyzprob0,xyzprob1,xyzprob2,nn,xyzprobc,d) &
    !$omp private(l,e,epq,dum)
    do l = 1,nn
      dum(1:3) = xyzprobc(1:3,l)
      call intermole_probe(d,dum,e,epq)
      xyzprob0(l) = e
      xyzprob1(l) = e+epq
      xyzprob2(l) = e-epq
    end do
    !$omp end parallel do

    do i = 1,nn
      ind0(i) = i
      ind1(i) = i
      ind2(i) = i
    end do
    call qsort(xyzprob0,1,nn,ind0)
    call qsort(xyzprob1,1,nn,ind1)
    call qsort(xyzprob2,1,nn,ind2)

    if (pr) then
      write (iunit,'(''  Best rare gas probe energy/kcal   :'',F8.2)') xyzprob0(1)*au
      write (iunit,'(''  +0.1 charged probe energy/kcal:'',F8.2)') xyzprob1(1)*au
      write (iunit,'(''  -0.1 charged probe energy/kcal:'',F8.2)') xyzprob2(1)*au
    end if

!>-- pocket search: cluster the best RG points and place B into each pocket
    if (cfg%pocket_grid) then
      if (pr) write (iunit,'(/,a,/)') '  Starting pocket search'
      if (nn/4 < cfg%mxcent_clust) cfg%mxcent_clust = nn/4

      allocate (atfrag(cfg%mxcent_clust),molvec(cfg%mxcent_clust), &
      &         iatnum(cfg%mxcent_clust),xyzfrag(3,cfg%mxcent_clust))
      do i = 1,cfg%mxcent_clust
        xyzfrag(1:3,i) = xyzprobc(1:3,ind0(i))
        atfrag(i) = probe_atom_type
      end do
      call mrec(mfrag,xyzfrag,cfg%mxcent_clust,atfrag,molvec)
      iatnum = 0
      do i = 1,cfg%mxcent_clust
        k = molvec(i)
        iatnum(k) = iatnum(k)+1
      end do

      do i = 1,mfrag
        k = iatnum(i)
        allocate (xyztmp(3,k),attmp(k))
        l = 0
        do j = 1,cfg%mxcent_clust
          if (molvec(j) .eq. i) then
            l = l+1
            xyztmp(1:3,l) = xyzfrag(1:3,j)
            attmp(l) = atfrag(j)
          end if
        end do
        call axis_docking(k,attmp,xyztmp,dum,pmass,dum2,dum3)
        r = sqrt(dum2(1)+dum2(2)+dum2(3))
        !> only pockets large enough to host B
        if (r/r2 .gt. 0.5.and.k .gt. 10) then
          displ = 0
          call alignmol(k,n2,attmp,d%at2,xyztmp,d%xyz2,displ)
          kk = kk+1
          found2(1:6,kk) = displ(1:6)
        end if
        deallocate (xyztmp,attmp)
      end do
      if (pr) write (iunit,'(''# pocket optimizations  :'',2i5  )') kk
      allocate (xyz_pocket(3,n,kk),pocket_e(kk),source=0.0_wp)

      do i = 1,kk
        found_tmp = found2(1:6,i)
        call place_b(found_tmp)
        call optimize_complex(etot,iostat)
        xyz_pocket(:,:,i) = comb%xyz(:,:)
        pocket_e(i) = etot
        if (pr) write (iunit,'(''pocket '',i3,'' E /au :'',f9.2)') i,etot
      end do

      if (.not.cfg%stack_grid.and..not.cfg%angular_grid) then
        !> pocket-only search: return the lowest pocket structure
        if (kk == 0) then
          if (pr) write (iunit,*) 'Error: no pocket found'
          io = 1
          return
        end if
        call sortxyz_e(n,kk,pocket_e,xyz_pocket)
        final_e = pocket_e
        xyz_opt = xyz_pocket
        if (cfg%wrfiles) call write_final(kk)
        best = comb
        best%xyz = xyz_opt(:,:,1)
        best%energy = final_e(1)
        eint = (final_e(1)-molA_e-molB_e)*au
        call print_eint(kk)
        return
      end if
    end if

!>-- stack search along 14 directions, 4 orientations each
    if (cfg%stack_grid) then
      if (pr) write (iunit,'(/,a)') '  Starting stack search'
      ii = cfg%mxcma
      nn = 4*14*ii
      if (cfg%cssym) nn = 4*9*ii
      if (pr) write (iunit,'(''   Grid points: '',i0  )') nn
      allocate (rprobc(6,nn),rprob(nn),source=0.0_wp)
      allocate (rind(nn),source=0)
      stepr4 = 0.2*float(n1+n2)**0.10d0

      l = 0
      do j = 1,14
        if (cfg%cssym) then
          if (mvec(icssym,j) .lt. 0) cycle
        end if
        r = stepr4
        if (j .gt. 6) r = stepr4/sqrt(3.)
        dum2(1:3) = mvec(1:3,j)*r
        do i = 1,ii
          displ = 0
          displ(1:3) = float(ii+1-i)*dum2(1:3)
          l = l+1
          rprobc(1:6,l) = displ(1:6)
          do k = 4,6
            displ(k) = pi
            l = l+1
            rprobc(1:6,l) = displ(1:6)
          end do
        end do
      end do

      !$omp parallel do default(none) reduction(+:rprob) &
      !$omp shared(nn,rprobc,d) private(l,e,displ)
      do l = 1,nn
        displ(1:6) = rprobc(1:6,l)
        call iff_e(d,.false.,-1,e,displ)
        rprob(l) = e
      end do
      !$omp end parallel do

      do i = 1,nn
        rind(i) = i
      end do
      call qsort(rprob,1,nn,rind)
      if (pr) write (iunit,'(''   lowest found /kcal :'',F8.2)') rprob(1)*au
      do i = 1,nn/100
        if (nclust >= ndim) exit
        if (rprob(i) .lt. 0.1) then
          nclust = nclust+1
          found(1:6,nclust) = rprobc(1:6,rind(i))
          found(7,nclust) = rprob(i)*au
        end if
      end do
      ll = nclust
    end if

!>-- angular search: full rotation grid on the crude R grid
    if (cfg%angular_grid) then
      if (pr) write (iunit,'(/,a)') '  Starting angular search'
      m1 = 360/int(cfg%stepa)
      maxlow = mm1*mm2*mm3+2
      nn = maxlow*m1*m1*m1
      if (pr) write (iunit,'(''   Grid points:'',i0  )') nn
      allocate (xyzprob3(3,maxlow),xyzprob4(7,nn),source=0.0_wp)

      !> the best +/- polarized probe points come first
      xyzprob3(1:3,1) = xyzprobc(1:3,ind1(1))
      xyzprob3(1:3,2) = xyzprobc(1:3,ind2(1))
      deallocate (xyzprobc,xyzprob0,xyzprob1,xyzprob2,ind0,ind1,ind2)

      dx = xx1
      l = 2
      do i = 1,mm1
        dum(1) = dx
        dy = yy1
        do j = 1,mm2
          dum(2) = dy
          dz = zz1
          do k = 1,mm3
            dum(3) = dz
            l = l+1
            xyzprob3(1:3,l) = dum(1:3)
            dz = dz+stepr3(3)
          end do
          dy = dy+stepr3(2)
        end do
        dx = dx+stepr3(1)
      end do
      !> m1 distinct angles in [0,2pi)
      cfg%stepa = pi2/float(m1)

      !$omp parallel do default(none) &
      !$omp shared(xyzprob4,maxlow,m1,cfg,xyzprob3,d) &
      !$omp private(icycle,i,j,k,e,displ,l)
      do icycle = 1,maxlow
        do i = 1,m1
          displ(4) = (i-1)*cfg%stepa
          do j = 1,m1
            displ(5) = (j-1)*cfg%stepa
            do k = 1,m1
              displ(6) = (k-1)*cfg%stepa
              displ(1:3) = xyzprob3(1:3,icycle)
              call iff_e(d,.false.,-1,e,displ)
              l = k+m1*(j-1)+m1**2*(i-1)+(icycle-1)*m1**3
              xyzprob4(1:6,l) = displ(1:6)
              xyzprob4(7,l) = e
            end do
          end do
        end do
      end do
      !$omp end parallel do
      call sort6(nn,xyzprob4)

      !> merge stack and angular candidates and rank them
      ll = min(5*cfg%maxparent+nclust,ndim,nclust+nn)
      !$omp parallel do default(none) &
      !$omp shared(found,ll,xyzprob4,nclust,d) private(i,displ,e)
      do i = 1,ll
        if (i .le. nclust) then
          displ(1:6) = found(1:6,i)
        else
          displ(1:6) = xyzprob4(1:6,i-nclust)
        end if
        call iff_e(d,.false.,-1,e,displ)
        found(1:6,i) = displ(1:6)
        found(7,i) = e*au
      end do
      !$omp end parallel do
    end if

    call doubles(ll,ndim,3.0_wp,0.1_wp,found)
    call sort6(ll,found)
    if (cfg%debug) call wrc3('best_before_gen.xyz',n1,n2,d%at1,d%at2,d%xyz1,d%xyz2,1,found)

    if (pr) then
      write (iunit,'(/,a)') '  Interaction energy of lowest structures so far in kcal/mol:'
      do i = 1,10
        write (iunit,'(   F14.2)') found(7,i)
      end do
      write (iunit,*)
    end if

!>-- genetic optimization
    if (cfg%n_opt == 0) then
      sig_count = 15
    else
      sig_count = cfg%n_opt
    end if
    if (pr) then
      write (iunit,*) '------------------------------'
      write (iunit,*) 'genetic optimization algorithm'
      write (iunit,*) '------------------------------'
      write (iunit,*) '  cycle  Eint/kcal/mol  average Eint'
    end if
    ii = cfg%maxparent**2
    do icycle = 1,cfg%maxgen
      !$omp parallel do default(none) &
      !$omp shared(found2,cfg,found,d) private(i,j,ii,displ,e,f)
      do i = 1,cfg%maxparent
        do j = 1,cfg%maxparent
          ii = (i-1)*cfg%maxparent+j
          call crossover(0.00_wp,f)
          displ(1:6) = found(1:6,i)*f(1:6)+found(1:6,j)*(1.0_wp-f(1:6))
          !> mutations only on children
          if (i .ne. j) call rand6(0.5_wp,1.0_wp,1.0_wp,displ)
          call iff_e(d,.false.,-1,e,displ)
          found2(1:6,ii) = displ(1:6)
          found2(7,ii) = e*au
        end do
      end do
      !$omp end parallel do

      ii = cfg%maxparent**2
      call doubles(ii,ndim,3.0_wp,0.10_wp,found2)
      call sort6(ii,found2)

      found = found2
      if (cfg%n_opt == 0) then
        call wrc2('structures_after_gen.xyz',0,n1,n2,d%at1,d%at2,d%xyz1,d%xyz2,15,found)
      else if (cfg%debug) then
        call wrc2('structures_after_gen.xyz',0,n1,n2,d%at1,d%at2,d%xyz1,d%xyz2,cfg%n_opt,found)
      end if
      if (cfg%debug.or.cfg%n_opt == 0) &
      &  call wrc3('best_after_gen.xyz',n1,n2,d%at1,d%at2,d%xyz1,d%xyz2,1,found)

      av = sum(found(7,1:cfg%maxparent))/float(cfg%maxparent)
      sig = 0
      do i = 1,sig_count
        sig = sig+(found(7,i)-av)**2
      end do
      sig = sqrt(sig/float(max(sig_count-1,1)))
      if (pr) write (iunit,'(4x,i0,7x,F7.1,5x,F8.1,5x)') icycle,found(7,1),av
      if (sig .lt. 0.3_wp) exit
    end do

    if (cfg%n_opt == 0) then
      if (pr) then
        write (iunit,*)
        write (iunit,*) "  Skipping geometry optimizations due to user request."
        write (iunit,*) "  The best structure after screening can be found on best_after_gen.xyz."
      end if
      call place_b(found(1:6,1))
      best = comb
      best%energy = found(7,1)/au
      eint = found(7,1)
      return
    end if

    if (cfg%docking_ens) then
      !> ensemble mode: optimize every gene with Eint < -0.1 kcal/mol
      cfg%n_opt = 0
      do i = 1,ii
        if (found(7,i) < -0.1_wp) then
          cfg%n_opt = cfg%n_opt+1
        else
          exit
        end if
      end do
      cfg%n_opt = max(cfg%n_opt,1)
    end if

!>-- final optimizations
    allocate (xyz_opt(3,n,cfg%n_opt+kk),final_e(cfg%n_opt+kk),source=0.0_wp)
    if (pr) write (iunit,'(/,''Optimizing '',i0,'' best structures with '',a )') cfg%n_opt,cfg%optlvl
    do icycle = 1,cfg%n_opt
      found_tmp = found(1:6,icycle)
      call place_b(found_tmp)
      if (pr) write (iunit,*) icycle
      call optimize_complex(etot,iostat)
      xyz_opt(:,:,icycle) = comb%xyz(:,:)
      final_e(icycle) = etot
    end do

    if (cfg%pocket_grid) then
      do k = 1,kk
        final_e(cfg%n_opt+k) = pocket_e(k)
        xyz_opt(:,:,cfg%n_opt+k) = xyz_pocket(:,:,k)
      end do
    end if

    call sortxyz_e(n,cfg%n_opt+kk,final_e,xyz_opt)
    if (cfg%wrfiles) call write_final(cfg%n_opt+kk)

    best = comb
    best%xyz(:,:) = xyz_opt(:,:,1)
    best%energy = final_e(1)
    eint = (final_e(1)-molA_e-molB_e)*au

    call print_eint(cfg%n_opt+kk)

  contains

    subroutine print_eint(nall)
      !> interaction energies of all final structures
      integer,intent(in) :: nall
      integer :: ic
      if (.not.pr) return
      write (iunit,*)
      write (iunit,*) '  ---------------------------'
      write (iunit,*) '     Interaction energies'
      write (iunit,*) '  ---------------------------'
      write (iunit,*) 'Attention: monomers are not optimized'
      write (iunit,*) 'Interaction energies are not physical'
      write (iunit,*)
      write (iunit,'(2x,''Lowest Interaction Energy:'',1x, F10.2, 1x, ''kcal/mol'' )') eint
      write (iunit,*) ' #   E_int (kcal/mol)'
      do ic = 1,nall
        write (iunit,'(2x,i0,3x,F10.2)') ic,(final_e(ic)-molA_e-molB_e)*au
      end do
    end subroutine print_eint

    subroutine place_b(gene)
      !> comb = A (fixed) + B moved by gene
      real(wp),intent(in) :: gene(6)
      call move2(molB%nat,molB%xyz,xyz_tmp,gene)
      comb%xyz(1:3,1:n1) = molA%xyz(1:3,1:n1)
      comb%xyz(1:3,n1+1:n) = xyz_tmp(1:3,1:n2)
    end subroutine place_b

    subroutine optimize_complex(energy,stat)
      !> optimize comb in place; like xtb, a failed run keeps its last geometry
      real(wp),intent(out) :: energy
      integer,intent(out)  :: stat
      energy = 0.0_wp
      grad = 0.0_wp
      !> GFN-FF: fresh topology at the start geometry (xtb restart_gff); the
      !> fragment split keeps A and B unbonded
      if (allocated(calc%calcs(1)%ff_dat)) deallocate (calc%calcs(1)%ff_dat)
      call optimize_geometry(comb,molopt,calc,energy,grad,.false.,.false.,stat)
      if (allocated(molopt%xyz)) comb%xyz = molopt%xyz
      if (stat < 0.or.isnan(energy)) then
        !> failed runs must not compete for the best structure (stat > 0
        !> only means not converged within maxcycle, kept as in xtb)
        energy = 1.0e42_wp
        if (pr) write (iunit,*) 'Warning: geometry optimization failed'
      else if (stat > 0.and.pr) then
        write (iunit,*) 'Warning: geometry optimization did not converge'
      end if
    end subroutine optimize_complex

    subroutine write_final(nall)
      !> sorted, optimized structures as final_structures.xyz and best.xyz
      integer,intent(in) :: nall
      integer :: ich,ic,jat
      open (newunit=ich,file='final_structures.xyz')
      do ic = 1,nall
        write (ich,'(i0)') n
        write (ich,'(f20.14)') final_e(ic)
        do jat = 1,n
          write (ich,'(a4,2x,3f20.14)') i2e(comb%at(jat),'nc'),xyz_opt(:,jat,ic)*autoang
        end do
      end do
      close (ich)
      open (newunit=ich,file='best.xyz')
      write (ich,'(i0)') n
      write (ich,'(f20.14)') final_e(1)
      do jat = 1,n
        write (ich,'(a4,2x,3f20.14)') i2e(comb%at(jat),'nc'),xyz_opt(:,jat,1)*autoang
      end do
      close (ich)
    end subroutine write_final

  end subroutine aiss_docking_search

!========================================================================================!

  subroutine sortxyz_e(n,nall,e,xyz)
    !> sort structures by energy (minimum first)
    integer,intent(in)     :: n,nall
    real(wp),intent(inout) :: e(nall),xyz(3,n,nall)
    integer :: ind(nall),i
    real(wp) :: xyz_tmp(3,n,nall)
    do i = 1,nall
      ind(i) = i
    end do
    call qsort(e,1,nall,ind)
    xyz_tmp = xyz
    do i = 1,nall
      xyz(:,:,i) = xyz_tmp(:,:,ind(i))
    end do
  end subroutine sortxyz_e

!========================================================================================!

  subroutine doubles(nlow,ndim,rthr,ethr,found)
    !******************************************************
    !* flag duplicate genes (similar geometry and energy,
    !* or identical energy) with E = 1e42
    !******************************************************
    integer,intent(in)     :: nlow,ndim
    real(wp),intent(in)    :: rthr,ethr
    real(wp),intent(inout) :: found(7,ndim)
    real(wp),allocatable :: scr(:,:)
    integer :: i,j
    real(wp) :: r,de,dem

    scr = found
    do i = 1,nlow-1
      do j = i+1,nlow
        call rnorm(i,j,found,r)
        de = abs(found(7,i)-found(7,j))
        dem = abs(found(7,i)+found(7,j))*0.5
        if (r .lt. rthr.and.de/dem .lt. ethr) scr(7,j) = 1.d+42
        !> chiral or symmetry-equivalent case
        if (abs(found(7,i)-found(7,j)) .lt. 1.d-3) scr(7,j) = 1.d+42
      end do
    end do
    found = scr
  end subroutine doubles

!========================================================================================!
end module aiss_search
