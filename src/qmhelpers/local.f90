! Copyright (C) 2017-2020 Stefan Grimme, Sebastian Ehlert (xtb)
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
!
! Routines were adapted from the xtb code (github.com/grimme-lab/xtb)
! under the Open-source software LGPL-3.0 Licencse.
!================================================================================!
module mo_localize
  use crest_parameters
  use lopt_mod,only:lopt
  use strucrd,only:i2e
  use wiberg_mayer
  implicit none
  private

  public :: local

  logical,parameter :: pr_local = .true.

!========================================================================================!
!========================================================================================!
contains  !> MODULE PROCEDURES START HERE
!========================================================================================!
!========================================================================================!

  subroutine local(nat,at,nbf,nao,ihomo,xyz,z, &
  &                focc,s,p,cmo,eig, &
  &                aoat2,aoat,nprim,alp,lao,cont, &
  &                nprot,protxyz )
    implicit none
    !> INPUT
    integer,intent(in) :: nat     !> number of atoms
    integer,intent(in) :: at(nat) !> atomic numbers for each atom
    integer,intent(in) :: nbf     !> number of basis functions
    integer,intent(in) :: nao     !> number of AOs
    integer,intent(in) :: ihomo   !> number of highest occ MO
    real(wp),intent(in) :: xyz(3,nat)   !> molecular coordinates (Bohr)
    real(wp),intent(in) :: z(nat)       !> nuclear charge
    real(wp),intent(in) :: focc(nao)    !> occupations
    real(wp),intent(in) :: s(nao,nao)   !> overlap S
    real(wp),intent(in) :: p(nao,nao)   !> density matrix P
    real(wp),intent(in) :: cmo(nao,nao) !> MO coefficients C
    real(wp),intent(in) :: eig(nao)     !> orbital energies
    !real(wp),intent(in) :: q(nat)       !> atomic charges

    integer,intent(inout)  :: aoat2(nao) !> mapping AO -> atom number
    integer,intent(in)  :: aoat(nbf)  !> mapping BF -> atom
    integer,intent(in)  :: nprim(nbf) !> mappinf BF -> primitive number
    real(wp),intent(in) :: alp(9*nbf) !> mapping primitive -> primitive exponent
    integer,intent(in)  :: lao(nbf)   !> mapping BF -> azimudal quantum number l
    real(wp),intent(in) :: cont(9*nbf)  !> mapping primitive -> contraction coeffient
   
    !> OUTPUT
    integer,intent(out) :: nprot
    real(wp),intent(out),allocatable :: protxyz(:,:)

    !> LOCAL
    real(wp),allocatable :: op(:,:)
    real(wp),allocatable :: oc(:,:,:)
    real(wp),allocatable :: cca(:)
    real(wp),allocatable :: dip2(:)
    real(wp),allocatable :: qua(:,:)
    real(wp),allocatable :: dip(:,:)
    real(wp),allocatable :: ecent(:,:),qcent(:,:),ecent2(:,:)
    real(wp),allocatable :: d(:,:)
    real(wp),allocatable :: f(:)
    real(wp),allocatable :: eiga(:)
    real(wp),allocatable :: xcen(:)
    real(wp),allocatable :: qmo(:,:)
    real(wp),allocatable :: tmpq(:,:)
    real(wp),allocatable :: rr(:)
    real(wp),allocatable :: wbo(:,:)
    real(wp),allocatable :: xyztmp(:,:)
    real(sp),allocatable :: rklmo(:,:)
    integer,allocatable :: ind(:)
    integer,allocatable :: lneigh(:,:)
    integer,allocatable :: aneigh(:,:)

    integer :: mo,n,i,j,k,ii,jj,imem(nat),idum,isc,iso,nop,pilist(ihomo)
    integer :: lamdu1,lamdu2,imo1,imo2,iso1,iso2,ij,jdum,maxlp,maxpi,is1
    integer :: ilumo,klev,nlev,nn,m,ldum(ihomo),sigrel(3,ihomo),npi,is2
    integer :: i1,i2,i3,is3,ipi,piset(nat),ncyc,j1,j2,smo,pmo,nl,nm,nr
    integer :: imo,new,nao2
    real(wp) :: dd,dum,det,pp,dtot(3),diptot,thr,t0,w0,t1,w1,vec1(3),v(6)
    real(wp) :: enlumo,enhomo,qhl(nat,2),ps,efh,efl,r1,r2,pithr,vec2(3)
    real(wp) :: praxis(3,3),aa,bb,cc
    character(len=80) :: atmp
    character(len=5) :: lmostring(4)
    data lmostring/'sigma','LP   ','pi   ','delpi'/
    !data lmostring/'σ    ','LP   ','π    ','del.π'/
    logical l1,l2,l3,flip
    integer LWORK,IWORK,LIWORK,INFO
    integer :: iscreen,icoord,ilmoi,icent ! file handles
    integer :: ntmpsave
    real(wp),allocatable :: tmpsave(:,:)
    

    !> BLAS
    external dgemm

!    if (pr_local) then
!      write (*,*)
!      write (*,*) 'localization/xTB-IFF output generation'
!    end if
    n = ihomo
    ilumo = n+1


!    if (ilumo .gt. nao) then
!      enlumo = 1000.d0
!    else
!      enlumo = eig(ilumo)
!    end if
!    efh = eig(n)
!    efl = enlumo
!
!    !> HOMO/LUMO populations for xTB FF (CT terms)
!    qhl = 0
!    thr = 0.01
!
!    if (ihomo .eq. 0) then
!      !> the HOMO is non-existing, place at unrealistically low energy
!      enhomo = -9999.999999990d0
!      qhl(1:nat,1) = 0.0d0
!    else
!      klev = 0
!      enhomo = 0
!      do nlev = ihomo,1,-1 ! loop over highest levels
!        if (efh-eig(nlev) .gt. thr) exit
!        klev = klev+1
!        enhomo = enhomo+eig(nlev)
!        do i = 1,nao
!          ii = aoat2(i)
!          do j = 1,i-1
!            jj = aoat2(j)
!            ps = s(j,i)*cmo(j,nlev)*cmo(i,nlev)
!            qhl(ii,1) = qhl(ii,1)+ps
!            qhl(jj,1) = qhl(jj,1)+ps
!          end do
!          ps = s(i,i)*cmo(i,nlev)*cmo(i,nlev)
!          qhl(ii,1) = qhl(ii,1)+ps
!        end do
!      end do
!      dum = 1.0_wp/float(klev)
!      enhomo = enhomo*dum
!      qhl(1:nat,1) = qhl(1:nat,1)*dum
!      if (pr_local) write (*,*) 'averaging CT terms over ',klev,' occ. levels'
!    end if
!
!    klev = 0
!    enlumo = 0
!    do nlev = ilumo,nao !> loop over highest levels
!      if (eig(nlev)-efl .gt. thr) exit
!      klev = klev+1
!      enlumo = enlumo+eig(nlev)
!      do i = 1,nao
!        ii = aoat2(i)
!        do j = 1,i-1
!          jj = aoat2(j)
!          ps = s(j,i)*cmo(j,nlev)*cmo(i,nlev)
!          qhl(ii,2) = qhl(ii,2)+ps
!          qhl(jj,2) = qhl(jj,2)+ps
!        end do
!        ps = s(i,i)*cmo(i,nlev)*cmo(i,nlev)
!        qhl(ii,2) = qhl(ii,2)+ps
!      end do
!    end do
!    dum = 1.0d0/float(klev)
!    enlumo = enlumo*dum
!    qhl(1:nat,2) = qhl(1:nat,2)*dum
!    if (pr_local) write (*,*) 'averaging CT terms over ',klev,' virt. levels'
!    if (ilumo .gt. nao) then
!      enlumo = 1000.0d0
!      qhl(1:nat,2) = 0
!    end if

    allocate (cca(nao*nao),xcen(n),lneigh(4,n),aneigh(2,n))
    allocate (d(n,n),ecent(n,4),eiga(n),qcent(n,3),ecent2(n,4))

    !> do only occ. ones
    cca = 0
    k = 0
    do i = 1,n
      k = k+1
      eiga(k) = eig(i)
      do j = 1,nao
        cca(j+(k-1)*nao) = cmo(j,i)
      end do
    end do

    nop = 3
!>--- dipole integrals
    allocate (dip2(nao*nao),dip(nbf*(nbf+1)/2,3),op(n*(n+1)/2,nop))
    call Dints(nat,nbf,xyz,dip(:,1),dip(:,2),dip(:,3), &
    &          aoat,nprim,alp,lao,cont)

    nao2 = nao
    call cao2saop(nbf,nao2,dip(:,1),lao)
    call cao2saop(nbf,nao2,dip(:,2),lao)
    call cao2saop(nbf,nao2,dip(:,3),lao)
    do mo = 1,3
      call onetri(1,dip(:,mo),dip2,cca,nao2,n)
      k = 0
      do i = 1,n
        do j = 1,i
          k = k+1
          op(k,mo) = dip2(j+(i-1)*n)
        end do
      end do
    end do

!>--- compute dipole moment core part
    dtot = 0
    do i = 1,nat
      dtot(1) = dtot(1)+xyz(1,i)*z(i)
      dtot(2) = dtot(2)+xyz(2,i)*z(i)
      dtot(3) = dtot(3)+xyz(3,i)*z(i)
    end do
    k = 0
    do i = 1,n
      k = k+i
      dtot(1) = dtot(1)+op(k,1)*focc(i)
      dtot(2) = dtot(2)+op(k,2)*focc(i)
      dtot(3) = dtot(3)+op(k,3)*focc(i)
    end do
    diptot = sqrt(dtot(1)**2+dtot(2)**2+dtot(3)**2)
    if (pr_local) then
      write (*,*) 'dipole moment from electron density (au)'
      write (*,*) '    X       Y       Z   '
      write (*,'(3f9.4,''  total (Debye): '',f8.3)') &
      &      dtot(1),dtot(2),dtot(3),diptot*2.5418_wp
    end if

!>--- do optimization (i.e. the actual localization)
    if (pr_local) write (*,*) 'doing rotations ...'
    call lopt(.true.,n,3,1.d-6,op,d)

    if (pr_local) write (*,*) 'doing transformations ...'
!>--- MO charge center
    k = 0
    do i = 1,n
      k = k+i
      ecent(i,1) = -op(k,1)
      ecent(i,2) = -op(k,2)
      ecent(i,3) = -op(k,3)
    end do

!>--- LMO fock matrix
    allocate (f(n))
    do i = 1,n
      dum = 0.0d0
      do k = 1,n
        dum = dum+eiga(k)*d(k,i)*d(k,i)
      end do
      f(i) = dum
    end do
    call lmosort2(n,f,d,ecent)

!>--- get the lmos
    CALL dgemm('N','N',nao2,n,n,1.D0,cca,nao2,d,n,0.D0,cmo,nao2) ! non-std BLAS

    pithr = 2.20d0  ! below is pi, large is pi deloc (typ=4)

!>--- number of centers for each mo
    allocate (qmo(nat,n))
    call mocent(nat,nao2,n,cmo,s,qmo,xcen,aoat2)

    allocate (rklmo(5,2*n))
    allocate (tmpsave(4,2*n), source=0.0_wp)
    ntmpsave=0
    !if (pr_local) write (*,*) 'lmo centers(Z=2) and atoms on file <lmocent.coord>'
    if (pr_local) write (*,*) ' LMO type  Fii/eV   ncent   charge center         contributions...'
!    if (pr_local) open (newunit=iscreen,file='xtbscreen.xyz')
    allocate (tmpq(nat,n))
    tmpq(1:nat,1:n) = qmo(1:nat,1:n)
    maxlp = 0
    maxpi = 0
    do i = 1,n
      do j = 1,nat
        imem(j) = j
      end do
      call lmosort(nat,n,i,imem,qmo)
      idum = 1
      do j = 1,nat
        if (qmo(j,i) .gt. 0.07) idum = j
      end do
      if (nat .eq. 1) then
        jdum = 2
      else
        call lmotype(nat,at,xyz,ecent(i,1),ecent(i,2),ecent(i,3), &
        &                imem(1),imem(2),xcen(i),.false.,pithr,jdum)
      end if
      rklmo(1:3,i) = ecent(i,1:3)
      rklmo(4,i) = ecent(i,4)
      rklmo(5,i) = real(jdum)
      if (jdum .eq. 2) maxlp = i
      if (jdum .eq. 3) maxpi = i
    end do
    qmo(1:nat,1:n) = tmpq(1:nat,1:n)
    deallocate (tmpq)

    do i = 1,n
      do j = 1,nat
        imem(j) = j
      end do
      call lmosort(nat,n,i,imem,qmo)
      idum = 1
      do j = 1,nat
        if (qmo(j,i) .gt. 0.07) idum = j
      end do
      if (nat .eq. 1) then
        jdum = 1
      else
        call lmotype(nat,at,xyz,ecent(i,1),ecent(i,2),ecent(i,3), &
        &                imem(1),imem(2),xcen(i),.true.,pithr,jdum)
      end if
      if (pr_local) then
        write (*,'(i5,1x,a,1x,2f7.2,3f10.5,12(i5,2x,a2,'':'',f6.2))')  &
        &   i,lmostring(jdum),autoev*f(i),xcen(i),ecent(i,1:3), &
        &   (imem(j),i2e(at(imem(j))),qmo(j,i),j=1,idum)
      end if

!>--- write + LP/pi as H for protonation search
      if (pr_local) then
        if (jdum .gt. 1) then
!          if (i .eq. maxlp.or.(i .eq. maxpi.and.maxlp .eq. 0)) then
!            open (newunit=icoord,file='coordprot.0')
!            write (icoord,'(''$coord'')')
!            do ii = 1,nat
!              write (icoord,'(3F24.10,5x,a2)') xyz(1:3,ii),i2e(at(ii))
!            end do
!            write (icoord,'(3F24.10,5x,a2)') ecent(i,1:3),i2e(1)
!            write (icoord,'(''$end'')')
!            close (icoord)
!          else
!            write (iscreen,*) nat+1
!            write (iscreen,*)
!            do ii = 1,nat
!              write (iscreen,'(a2,3F24.10)') i2e(at(ii)),xyz(1:3,ii)*autoaa
!            end do
!            write (iscreen,'(a2,3F24.10)') i2e(1),ecent(i,1:3)*autoaa
!          end if
        end if
      end if

!>--- tmpsave
      if(jdum.gt.1)then
        ntmpsave = ntmpsave + 1
        tmpsave(1:3, ntmpsave) = ecent(i,1:3)
        tmpsave(4,ntmpsave) = float(jdum)
      endif

    end do

    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    new = 0
    if (maxval(rklmo(5,1:n)) .ge. 3.0d0) then
      if (pr_local) write (*,*) 'starting deloc pi regularization ...'

      call atomneigh(n,nat,xyz,ecent,aneigh) ! the nearest  atom neighbors for each LMO
      call lmoneigh(n,rklmo,ecent,aneigh,lneigh)   ! the nearest LMO neighbor but not LP

      ldum = -1
      do i = 1,n
        if (int(rklmo(5,i)) .ge. 3.and.ldum(i) .ne. 1) then
          flip = .false.
          !> check for same alignment of pi orbitals
          call piorient(ecent(i,:),ecent(lneigh(1,i),:),flip)
          if (flip) then
            ldum(i) = 1
            !> i itself is the pi DOWN counterpart and lneigh(1,i) contains the pi UP part
            ldum(lneigh(1,i)) = 0
          else
            ldum(i) = 0
            ldum(lneigh(1,i)) = 1  !> lneigh(1,i) contains the pi DOWN counterpart
          end if
        end if
      end do
      k = 0
      piset = 0
      do i = 1,n
        if (ldum(i) .eq. 0) then
          k = k+1
          pilist(k) = i            !> UP plane unique pi list
          piset(aneigh(1,i)) = 1
          piset(aneigh(2,i)) = 1
        end if
      end do
      npi = k

!>--- determine for each deloc pi bond sigma bonds which are in the pi system
      j = 0
      ldum = 0
      nn = 0
      sigrel = 0
      do i = 1,npi
        ipi = pilist(i)
        call irand3(i1,i2,i3)    !> change order of lookup
        is1 = lneigh(i1+1,ipi)
        is2 = lneigh(i2+1,ipi)
        is3 = lneigh(i3+1,ipi)
        !> do bonds on ipi belong to the pi system?
        l1 = bndcheck(nat,piset,aneigh(1,is1),aneigh(2,is1))
        l2 = bndcheck(nat,piset,aneigh(1,is2),aneigh(2,is2))
        l3 = bndcheck(nat,piset,aneigh(1,is3),aneigh(2,is3))
        if (l1.and.ldum(is1) .eq. 0) then
          ldum(is1) = 1
          sigrel(1,i) = is1
        end if
        if (l2.and.ldum(is2) .eq. 0) then
          ldum(is2) = 1
          sigrel(2,i) = is2
        end if
        if (l3.and.ldum(is3) .eq. 0) then
          ldum(is3) = 1
          sigrel(3,i) = is3
        end if
      end do

      ldum = 0
      do i = 1,npi
        ipi = pilist(i)
        do j = 1,3
          is1 = sigrel(j,i)
          if (is1 .gt. 0) ldum(is1) = ldum(is1)+1
        end do
      end do

      if (pr_local) write (*,'(1x,a,f6.2,a,i0)') 'ncent threshold ',pithr,' --> # pi deloc LMO ',npi

      allocate (wbo(nat,nat))
      wbo = 0.0d0
      call get_wbo_rhf(nat,nao2,P,S,aoat2,wbo)

!>--- now create new LMO
      k = 0
      m = 0
      do i = 1,npi
        pmo = pilist(i)
        i1 = aneigh(1,pmo)
        i2 = aneigh(2,pmo)
        dd = wbo(i2,i1)
        do j = 1,3
          smo = sigrel(j,i)
          if (rklmo(5,smo) .ne. 1) cycle
          j1 = aneigh(1,smo)
          j2 = aneigh(2,smo)
          dum = wbo(j2,j1)
          !> take difference of pi-bond order (i.e., WBO > 1) relative to its sum (i.e., number of pi-electrons distributed among the 3 centers)
          pp = abs((dum-dd)/(dum+dd-2.0d0))
          if (pp .ge. 0.50d0) cycle
          !> pi atom  center  sigma (nl,nm,nr)
          !> from 3 atoms A-B=C, find central atom (B)
          call threeoutfour(i1,i2,j1,j2,nl,nm,nr)
          if (nl .eq. 0) cycle !> fall back =do nothing if assignment fails
          vec1(1:3) = xyz(1:3,nm)
          vec2(1:3) = (xyz(1:3,nl)+xyz(1:3,nr))*0.5
          dtot(1:3) = rklmo(1:3,pmo)
          call calcrotation(dtot,vec2,vec1-vec2,pi)   !> project upper pi center to other bond
          k = k+1
          rklmo(1:3,n+k) = dtot(1:3)
          rklmo(4:5,n+k) = rklmo(4:5,pmo)

          !> add to screen file, protomer search
          if (pr_local) then
!            write (iscreen,*) nat+1
!            write (iscreen,*)
!            do ii = 1,nat
!              write (iscreen,'(a2,3F24.10)') i2e(at(ii)),xyz(1:3,ii)*autoaa
!            end do
!            write (iscreen,'(a2,3F24.10)') i2e(1),dtot(1:3)*autoaa
          end if
          !>--- tmpsave
          ntmpsave = ntmpsave + 1
          tmpsave(1:3, ntmpsave) = dtot(1:3)
          tmpsave(4,ntmpsave) = rklmo(5,pmo)


          imo = lneigh(1,pmo)
          vec1(1:3) = xyz(1:3,nm)
          vec2(1:3) = (xyz(1:3,nl)+xyz(1:3,nr))*0.5_wp
          dtot(1:3) = rklmo(1:3,imo)
          call calcrotation(dtot,vec2,vec1-vec2,pi)  !> project bottom pi center to other bond
          k = k+1
          rklmo(1:3,n+k) = dtot(1:3)
          rklmo(4:5,n+k) = rklmo(4:5,imo)
          rklmo(5,smo) = 0 ! remove sigma

          !> add to screen file, protomer search
          if (pr_local) then
!            write (iscreen,*) nat+1
!            write (iscreen,*)
!            do ii = 1,nat
!              write (iscreen,'(a2,3F24.10)') i2e(at(ii)),xyz(1:3,ii)*autoaa
!            end do
!            write (iscreen,'(a2,3F24.10)') i2e(1),dtot(1:3)*autoaa
          end if
          !>--- tmpsave
          ntmpsave = ntmpsave + 1
          tmpsave(1:3, ntmpsave) = dtot(1:3)
          tmpsave(4,ntmpsave) = rklmo(5,imo)

          m = m+1
        end do
      end do
      new = k

      if (pr_local) then
!        close (iscreen)
      end if
      deallocate (wbo)

    end if
 !ccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    !> collect possible protonation sites for output
!    nprot = 0
!    k = new+n
!    write(*,*) k,new,n
!    do i=1,k
!      if(rklmo(5,i) > 1) nprot=nprot+1
!    enddo
    nprot = ntmpsave
    if(nprot > 0)then
      allocate(protxyz(4,nprot), source=0.0_wp)
!      m=0
!      do i=1,k
!        if(nint(rklmo(5,i)) > 1)then
!          m=m+1
!          protxyz(1:3,m) = rklmo(1:3,i)  
!          protxyz(4,m) = rklmo(5,i)
!        endif
!      enddo  
      do i=1,nprot
         protxyz(:,i) = tmpsave(:,i)
      enddo
    endif
    deallocate(tmpsave)

 !ccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    !> If the normal xtb mode is used with --lmo, pr_local is true and
    !  a file with all information is written.
    !  If the docking mode is used, information are stored intrenally

!    if (pr_local) then
!      call open_file(ilmoi,set%lmoinfo_fname,'w')
!      write (ilmoi,*) nat
!      if (pr_local) call open_file(icent,'lmocent.coord','w')
!      if (pr_local) write (icent,'(''$coord'')')
!      do i = 1,nat
!        write (ilmoi,'(i2,3F20.10,E18.10)') at(i),xyz(1:3,i),q(i)
!        if (pr_local) write (icent,'(3F24.10,5x,a2)') xyz(1:3,i),i2e(at(i))
!      end do
!      k = 0
!      do i = 1,n+new
!        if (int(rklmo(5,i)) .gt. 0) k = k+1
!      end do
!      write (ilmoi,'(i5,5f14.8)') k,diptot,enlumo,enhomo
!      do i = 1,n+new
!        if (int(rklmo(5,i)) .gt. 0) then
!          write (ilmoi,'(i2,3F20.10,f14.8)') int(rklmo(5,i)),rklmo(1:3,i)
!          if (pr_local) write (icent,'(3F24.10,5x,a2)') rklmo(1:3,i),i2e(2)
!        end if
!      end do
!      write (ilmoi,'(10(F10.6))') (qhl(i,1),i=1,nat)  ! HOMO atom pop
!      write (ilmoi,'(10(F10.6))') (qhl(i,2),i=1,nat)  ! LUMO atom pop
!      if (pr_local) write (icent,'(''$end'')')
!      call close_file(ilmoi)
!      if (pr_local) call close_file(icent)
!    end if

    if (pr_local) then
!      write (*,*) 'files:'
!      write (*,*) 'coordprot.0/xtbscreen.xyz/xtblmoinfo/lmocent.coord'
!      write (*,*) 'with protonation site input, xtbdock and'
!      write (*,*) 'LMO center info written'
!      write (*,*)
    end if

    deallocate (xcen,cca,d,f,qmo,ecent,rklmo)

  end subroutine local

!========================================================================================!
  subroutine lmotype(n,at,xyz,ex,ey,ez,ia1,ia2,xcen,modi,pithr,typ)
!******************************************************
!* lmotype
!* determine type of LMO based on empirical thresholds
!******************************************************
    implicit none
    integer,intent(in)    :: n,ia1,ia2,at(n)
    integer,intent(out)   :: typ
    logical,intent(in)    :: modi
    real(wp),intent(in)    :: xyz(3,n)
    real(wp),intent(inout) :: ex,ey,ez
    real(wp),intent(in)    :: xcen,pithr
    real(wp) :: r1,r2,r,r0,f,norm

    if (xcen .lt. 1.3333333_wp) then
      typ = 2 ! LP
      f = -2.2_wp
      if (modi) then
        norm = sqrt((xyz(1,ia1)-ex)**2+(xyz(2,ia1)-ey)**2 &
        &              +(xyz(3,ia1)-ez)**2)
        ! the LP is on the atom e.g. in TM compounds so just shift it away
        if (norm .lt. 0.2) then
          call shiftlp(n,at,ia1,xyz,ex,ey,ez)
        else
          ! put it on the line along atom...LP center
          ex = ex+f*(xyz(1,ia1)-ex)
          ey = ey+f*(xyz(2,ia1)-ey)
          ez = ez+f*(xyz(3,ia1)-ez)
        end if
      end if
    else
      r1 = sqrt((xyz(1,ia1)-ex)**2+(xyz(2,ia1)-ey)**2+(xyz(3,ia1)-ez)**2)
      r2 = sqrt((xyz(1,ia2)-ex)**2+(xyz(2,ia2)-ey)**2+(xyz(3,ia2)-ez)**2)
      r = r1+r2
      r0 = sqrt(sum((xyz(:,ia1)-xyz(:,ia2))**2))
      if (r/r0 .gt. 1.04_wp) then
        typ = 3  ! pi
        if (xcen .gt. pithr) typ = 4
      else
        typ = 1  ! sigma
      end if
    end if

  end subroutine lmotype

!========================================================================================!
  subroutine mocent(n,ndim,ihomo,x,s,qmo,xcen,aoat2)
    implicit none
    integer :: n,ndim,ihomo
    integer,intent(in) :: aoat2(ndim)
    real(wp) :: xcen(ihomo),x(ndim,ndim),s(ndim,ndim)
    real(wp) :: qmo(n,ihomo)

    integer  :: ia,ib,m,i,j,k
    real(wp) :: dum,tmp

    qmo = 0
    !$OMP  PARALLEL PRIVATE (m,j,ia,k,ib,tmp,dum) &
    !$OMP& SHARED(qmo,xcen,X,S) &
    !$OMP& DEFAULT(SHARED)
    !$OMP DO
    do m = 1,ihomo
      do j = 1,ndim
        ia = aoat2(j)
        do k = 1,j-1
          ib = aoat2(k)
          tmp = X(j,m)*X(k,m)*s(k,j)
          qmo(ia,m) = qmo(ia,m)+tmp
          qmo(ib,m) = qmo(ib,m)+tmp
        end do
        qmo(ia,m) = qmo(ia,m)+X(j,m)*X(j,m)*s(j,j)
      end do
      dum = 0
      do i = 1,n
        dum = dum+qmo(i,m)**2
      end do
      xcen(m) = 1.0/(1.0d-8+dum)
    end do
    !$OMP END DO
    !$OMP END PARALLEL

  end subroutine mocent

  SUBROUTINE lmosort(ncent,ihomo,imo,imem,qmo)
    integer :: ncent,ihomo,imo
    integer :: imem(ncent)
    real(wp) :: qmo(ncent,ihomo)
    integer :: ii,i,k,j,ihilf
    real(wp) :: pp
    do ii = 2,ncent
      i = ii-1
      k = i
      pp = qmo(i,imo)
      do j = ii,ncent
        if (qmo(j,imo) .lt. pp) cycle
        k = j
        pp = qmo(j,imo)
      end do
      if (k .eq. i) cycle
      qmo(k,imo) = qmo(i,imo)
      qmo(i,imo) = pp

      ihilf = imem(i)
      imem(i) = imem(k)
      imem(k) = ihilf
    end do

  end subroutine lmosort

!========================================================================================!
  SUBROUTINE lmosort2(n,eps,d,ecent)
    real(wp) :: d(n,n),eps(n),ecent(n,3)
    integer :: n
    integer :: ii,j,i,k
    real(wp) :: hilf,pp
    do ii = 2,n
      i = ii-1
      k = i
      pp = eps(i)
      do j = ii,n
        if (eps(j) .gt. pp) cycle
        k = j
        pp = eps(j)
      end do
      if (k .eq. i) cycle
      eps(k) = eps(i)
      eps(i) = pp

      do j = 1,n
        hilf = d(j,i)
        d(j,i) = d(j,k)
        d(j,k) = hilf
      end do

      do j = 1,3
        hilf = ecent(i,j)
        ecent(i,j) = ecent(k,j)
        ecent(k,j) = hilf
      end do

    end do

  end subroutine lmosort2

!========================================================================================!
  subroutine lmoneigh(n,rk,ecent,aneigh,neigh)
    implicit none
    integer :: n,neigh(4,n),aneigh(2,*)
    real(wp) :: ecent(n,3)
    real(sp) :: rk(5,2*n)

    real(wp) rr(n)
    integer :: i,j,ind(n),i1,i2,j1,j2
    external qsort

    do i = 1,n
      i1 = aneigh(1,i)
      i2 = aneigh(2,i)
      do j = 1,n
        j1 = aneigh(1,j)
        j2 = aneigh(2,j)
        rr(j) = (ecent(i,1)-ecent(j,1))**2 &
           &            +(ecent(i,2)-ecent(j,2))**2 &
           &            +(ecent(i,3)-ecent(j,3))**2
        if (int(rk(5,j)) .eq. 2) rr(j) = 1.d+42  ! LP
        if (i1 .eq. j1.and.i2 .eq. j2) rr(j) = -1.0 ! same bond?
        if (i1 .eq. j2.and.i2 .eq. j1) rr(j) = -1.0
        if (i .eq. j) rr(j) = 2.d+42
        ind(j) = j
      end do
      call qsort(rr,1,n,ind)
      neigh(1,i) = ind(1)
      neigh(2,i) = ind(2)
      neigh(3,i) = ind(3)
      neigh(4,i) = ind(4)
    end do

  end subroutine lmoneigh

!========================================================================================!
  subroutine atomneigh(n,nat,xyz,ecent,neigh)
    implicit none
    integer,intent(in)    :: n,nat
    integer,intent(inout) :: neigh(2,n)
    real(wp),intent(in)    :: ecent(n,3)
    real(wp),intent(in)    :: xyz(3,nat)

    real(wp) :: rr(nat)
    integer  :: i,j,ind(nat)
    external qsort
    do i = 1,n
      do j = 1,nat
        rr(j) = (ecent(i,1)-xyz(1,j))**2 &
           &  +(ecent(i,2)-xyz(2,j))**2 &
           &  +(ecent(i,3)-xyz(3,j))**2
        ind(j) = j
      end do
      call qsort(rr,1,nat,ind)
      neigh(1,i) = ind(1)
      neigh(2,i) = ind(2)
    end do

  end subroutine atomneigh

!========================================================================================!
  pure function bndcheck(nat,list,i1,i2) result(check)
    integer,intent(in) :: nat,list(nat),i1,i2
    logical :: check
    check = .false.
    if (list(i1) .eq. 1.and.list(i2) .eq. 1) check = .true.
  end function bndcheck

!========================================================================================!
  subroutine irand3(n1,n2,n3)
    integer,intent(out) :: n1,n2,n3
    integer  :: irand
    real(sp) :: x

    call random_number(x)
    irand = int(3.1*x)
    if (irand .lt. 1) irand = 1
    if (irand .gt. 3) irand = 3
    n1 = irand
10  call random_number(x)
    irand = int(3.1*x)
    if (irand .lt. 1) irand = 1
    if (irand .gt. 3) irand = 3
    if (irand .ne. n1) then
      n2 = irand
    else
      goto 10
    end if
20  call random_number(x)
    irand = int(3.1*x)
    if (irand .lt. 1) irand = 1
    if (irand .gt. 3) irand = 3
    if (irand .ne. n1.and.irand .ne. n2) then
      n3 = irand
    else
      goto 20
    end if

  end subroutine irand3

!========================================================================================!
  pure subroutine threeoutfour(i1,i2,j1,j2,n1,n2,n3)
    implicit none
    integer,intent(in)  :: i1,i2,j1,j2
    integer,intent(out) :: n1,n2,n3

    n1 = 0
    if (i1 .eq. j1) then
      n1 = i2
      n2 = i1
      n3 = j2
      return
    end if
    if (i1 .eq. j2) then
      n1 = i2
      n2 = i1
      n3 = j1
      return
    end if
    if (i2 .eq. j1) then
      n1 = i1
      n2 = i2
      n3 = j2
      return
    end if
    if (i2 .eq. j2) then
      n1 = i1
      n2 = i2
      n3 = j1
      return
    end if

  end subroutine threeoutfour

!========================================================================================!
  pure subroutine calcrotation(x,ori,vec,phi)
    implicit none
    integer :: i,j
    real(wp),intent(inout) :: x(3)
    real(wp),intent(in)    :: ori(3)
    real(wp),intent(in)    :: vec(3)
    real(wp),intent(in)    :: phi
    real(wp) :: d(3)
    real(wp) :: xtmp(3)
    real(wp) :: absd

    d(1:3) = vec(1:3)

    xtmp(1:3) = x(1:3)-ori(1:3)

    absd = sqrt(d(1)**2+d(2)**2+d(3)**2)

    d(1:3) = d(1:3)/absd

    x(1) = ((d(2)**2+d(3)**2)*cos(phi)+d(1)**2)*xtmp(1)+  &
    &      (d(1)*d(2)*(1-cos(phi))-d(3)*sin(phi))*xtmp(2)+  &
    &      (d(1)*d(3)*(1-cos(phi))+d(2)*sin(phi))*xtmp(3)
    x(2) = (d(1)*d(2)*(1-cos(phi))+d(3)*sin(phi))*xtmp(1)+  &
    &      ((d(1)**2+d(3)**2)*cos(phi)+d(2)**2)*xtmp(2)+  &
    &    (d(2)*d(3)*(1-cos(phi))-d(1)*sin(phi))*xtmp(3)
    x(3) = (d(1)*d(3)*(1-cos(phi))-d(2)*sin(phi))*xtmp(1)+  &
    &      (d(2)*d(3)*(1-cos(phi))+d(1)*sin(phi))*xtmp(2)+  &
    &      ((d(1)**2+d(2)**2)*cos(phi)+d(3)**2)*xtmp(3)
    x(1:3) = x(1:3)+ori(1:3)

  end subroutine calcrotation

!========================================================================================!
  pure subroutine piorient(a,b,flip)
    implicit none
    real(wp),intent(in)  :: a(3),b(3)
    logical,intent(out) :: flip
    integer  :: i,j
    real(wp) :: d(3),d2(3)
    flip = .false.
    d = 0.0d0
    j = 0

    do i = 1,3
      d(i) = a(i)-b(i)
      d2(i) = d(i)*d(i)
    end do
    j = maxloc(d2,1)
    if (d(j) .lt. 0.0d0) flip = .true.
    return
  end subroutine piorient

!========================================================================================!

  subroutine Dints(n,nbf,xyz,S1,S2,S3, &
  &                aoat,nprim,alp,lao,cont)
    use xtb_intpack,only:opab1,propa
    implicit none
    !> INPUT
    integer,intent(in)  :: n
    integer,intent(in)  :: nbf
    real(wp),intent(in)  :: xyz(3,n)
    integer,intent(in)  :: aoat(nbf)   !> mapping BF -> atom
    integer,intent(in)  :: nprim(nbf)  !> mappinf BF -> primitive number
    real(wp),intent(in) :: alp(9*nbf)  !> mapping primitive -> primitive exponent
    integer,intent(in)  :: lao(nbf)    !> mapping BF -> azimudal quantum number l
    real(wp),intent(in) :: cont(9*nbf) !> mapping primitive -> contraction coeffient
    !> OUTPUT
    real(wp),intent(out) :: S1(nbf*(nbf+1)/2)
    real(wp),intent(out) :: S2(nbf*(nbf+1)/2)
    real(wp),intent(out) :: S3(nbf*(nbf+1)/2)
    !> LOCAL
    integer  :: i,j,k,l
    integer  :: iprimcount,jprimcount
    integer  :: npri,nprj
    integer  :: ii,iii,jj,jjj,ll,m,li,lj,mm,nn
    real(wp) :: xyza(3),xyzb(3),rab,est,ss,sss,gama,arg
    real(wp) :: point(3),gm2,ttt(3),tt1,tt2,tt3,intcut

    intcut = 20.0_wp

    point = 0.0_wp
    s1 = 0.0_wp
    s2 = 0.0_wp
    s3 = 0.0_wp

    k = 0
    iprimcount = 0
    do i = 1,nbf
      !> aufpunkt i
      xyza(1:3) = xyz(1:3,aoat(i))
      !> #prims
      npri = nprim(i)
      jprimcount = 0
      do j = 1,i
        k = k+1
        nprj = nprim(j)
        !> aufpunkt j
        xyzb(1:3) = xyz(1:3,aoat(j))
        rab = sum((xyza-xyzb)**2)
        if (rab .gt. 200) goto 42 !> cut-off gives crambin dipole accurate to 1d-3 Deb
        !> prim loop
        tt1 = 0.0_wp
        tt2 = 0.0_wp
        tt3 = 0.0_wp
        do ii = 1,npri
          iii = iprimcount+ii
          do jj = 1,nprj
            jjj = jprimcount+jj
            gama = 1.0_wp/(alp(iii)+alp(jjj))
            est = rab*alp(iii)*alp(jjj)*gama
            !> cutoff
            if (est .lt. intcut) then
              ttt = 0
              call propa(opab1,xyza,xyzb,point,alp(iii),alp(jjj),&
                 &       lao(i),lao(j),ttt,3)
              tt1 = tt1+ttt(1)*cont(iii)*cont(jjj)
              tt2 = tt2+ttt(2)*cont(iii)*cont(jjj)
              tt3 = tt3+ttt(3)*cont(iii)*cont(jjj)
            end if
          end do
        end do
        s1(k) = tt1
        s2(k) = tt2
        s3(k) = tt3
42      jprimcount = jprimcount+nprj
      end do
      iprimcount = iprimcount+npri
    end do

  end subroutine Dints

!========================================================================================!

  subroutine onetri(ity,s,s1,array,n,ival)
    !     ******designed for abelian groups only******
    !
    !     calling sequence:
    !     ity       =1 sym ao
    !               =-1 anti sym ao
    !     s         input property matrix over ao's
    !     s1        transformed integrals on output
    !     s2        scratch arrays large enough to hold a square matrix
    !     array     mo matrix over so's
    !     n         linear dimension of arrays
    !     bernd hess, university of bonn, january 1991
    integer,intent(in) :: n
    integer,intent(in) :: ity
    integer,intent(in) :: ival
    real(wp),intent(in) :: s(*)
    real(wp),intent(inout) :: s1(*)
    real(wp),intent(in) :: array(n,ival)
    real(wp) :: s2(n,n)
    external :: dsymm,dgemm ! can't use wrappers due to change of leading dim

    !
    !     determine if we have an antisymmetric integral
    if (ity /= -1) then
      !
      !     blow up symmetric matrix s
      call blowsy(ity,s,s1,n)

      !
      !     transformation of s
      call dsymm('l','l',n,ival,1.d0,s1,n,array,n,0.d0,s2,n)
      call dgemm('t','n',ival,ival,n,1.d0,array,n,s2,n,0.d0,s1,ival)
    else
      !
      !     blow up anti-symmetric matrix s
      call blowsy(ity,s,s1,n)
      !
      !     transformation of s
      call dgemm('n','n',n,ival,n,1.d0,s1,n,array,n,0.d0,s2,n)
      call dgemm('t','n',ival,ival,n,1.d0,array,n,s2,n,0.d0,s1,ival)
    end if
  end subroutine

!========================================================================================!

  pure subroutine blowsy(ity,a,b,n)
    implicit none
    ! blow up symmetric or antisymmetric matrix to full size
    integer,intent(in)  :: ity
    integer,intent(in)  :: n
    real(wp),intent(in)  :: a(*)
    real(wp),intent(out) :: b(n,n)
    integer :: ij,i,j
    ! determine if we have an antisymmetric integral

    ij = 0
    if (ity .eq. -1) then
      do i = 1,n
        do j = 1,i-1
          ij = ij+1
          b(j,i) = -a(ij)
          b(i,j) = a(ij)
        end do
        ij = ij+1
        b(i,i) = 0.d0
      end do
    else
      do i = 1,n
        do j = 1,i-1
          ij = ij+1
          b(j,i) = a(ij)
          b(i,j) = a(ij)
        end do
        ij = ij+1
        b(i,i) = a(ij)
      end do
    end if
  end subroutine blowsy

!========================================================================================!

  subroutine cao2saop(nbf,nao,s,lao)
    implicit none
    !> INPUT
    integer,intent(in)    :: nbf      !> number of BF
    integer,intent(inout) :: nao      !> number of MOs
    integer,intent(in)    :: lao(nbf) !> mapping BF -> azimudal quantum number l
    !> OUTPUT
    real(wp),intent(inout) :: s(nbf*(nbf+1)/2)
    !> LOCAL
    real(wp) :: sspher
    integer ::  lll(20),firstd(nbf)
    integer :: i,j,k,ii,jj,kk,iii,jjj,li,lj,mm,nn,m,n
    data lll/1,2,2,2,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4/

    real(wp) :: trafo(6,6)
    real(wp),allocatable :: sneu(:)

    trafo = 0.0d0
    ! dS
    trafo(1,1) = 1.d0/sqrt(3.d0)*sqrt(3.d0/5.d0)
    trafo(2,1) = 1.d0/sqrt(3.d0)*sqrt(3.d0/5.d0)
    trafo(3,1) = 1.d0/sqrt(3.d0)*sqrt(3.d0/5.d0)
    ! dx2-y2
    trafo(1,2) = 1.d0/sqrt(2.d0)*sqrt(3.d0/2.d0)
    trafo(2,2) = -1.d0/sqrt(2.d0)*sqrt(3.d0/2.d0)
    ! dz2
    trafo(1,3) = 0.50d0
    trafo(2,3) = 0.50d0
    trafo(3,3) = -1.0d0
    ! rest
    trafo(4,4) = 1.0d0
    trafo(5,5) = 1.0d0
    trafo(6,6) = 1.0d0

    nao = 0
    firstd = 0
    i = 1
42  if (lao(i) .gt. 4.and.lao(i) .le. 10) then
      nao = nao+1
      firstd(i:i+5) = i
      i = i+5
    end if
    i = i+1
    if (i .lt. nbf) goto 42

    if (nao .eq. 0) then
      nao = nbf
      return
    end if

    allocate (sneu(nbf*(nbf+1)/2))

    sneu = s

    k = 0
    do i = 1,nbf
      li = lll(lao(i))
      do j = 1,i
        lj = lll(lao(j))
        k = k+1
        ! d-d
        if (li .eq. 3.and.lj .eq. 3) then
          ii = lao(i)-4
          jj = lao(j)-4
          sspher = 0
          do m = 1,6
            mm = firstd(i)-1+m
            do n = 1,6
              nn = firstd(j)-1+n
              sspher = sspher+trafo(m,ii)*trafo(n,jj)*s(lin(mm,nn))
            end do
          end do
          sneu(k) = sspher
        end if
        ! d-sp
        if (li .eq. 3.and.lj .le. 2) then
          ii = lao(i)-4
          sspher = 0
          do m = 1,6
            mm = firstd(i)-1+m
            sspher = sspher+trafo(m,ii)*s(lin(mm,j))
          end do
          sneu(k) = sspher
        end if
        ! sp-d
        if (li .le. 2.and.lj .eq. 3) then
          jj = lao(j)-4
          sspher = 0
          do n = 1,6
            nn = firstd(j)-1+n
            sspher = sspher+trafo(n,jj)*s(lin(i,nn))
          end do
          sneu(k) = sspher
        end if

      end do
    end do

    s(1:nao*(nao+1)/2) = 0

    k = 0
    iii = 0
    do i = 1,nbf
      if (lao(i) .ne. 5) iii = iii+1
      jjj = 0
      do j = 1,i
        if (lao(j) .ne. 5) jjj = jjj+1
        k = k+1
        if (lao(i) .eq. 5.or.lao(j) .eq. 5) cycle
        s(lin(iii,jjj)) = sneu(k)
      end do
    end do

    nao = nbf-nao

    deallocate (sneu)
  end subroutine cao2saop

!========================================================================================!

  integer function lin(i1,i2)
    integer ::  i1,i2,idum1,idum2
    idum1 = max(i1,i2)
    idum2 = min(i1,i2)
    lin = idum2+idum1*(idum1-1)/2
    return
  end function lin

!========================================================================================!

  subroutine shiftlp(n,at,ia1,xyz,ex,ey,ez)
! shift the LP=protonation position to an "empty" region
    implicit none
    integer  :: n,at(n),ia1
    real(wp) :: xyz(3,n),ex,ey,ez

    real(wp) :: r0,r1,r2,exmin,eymin,ezmin,rmin,xsave,ysave,zsave,rmax
    real(wp) :: x,f
    integer :: icount,j

    r0 = 2.5d0
    rmin = 1.d+42
    rmax = 0
    icount = 0

10  f = 1.0
    call random_number(x)
    if (x .lt. 0.5) f = -1.0
    call random_number(x)
    ex = xyz(1,ia1)+x*f*r0
    call random_number(x)
    ey = xyz(2,ia1)+x*f*r0
    call random_number(x)
    ez = xyz(3,ia1)+x*f*r0
    r1 = sqrt((xyz(1,ia1)-ex)**2+(xyz(2,ia1)-ey)**2+(xyz(3,ia1)-ez)**2)
    if (abs(r0-r1) .gt. 0.1) goto 10
    rmin = 1.d+42
    do j = 1,n
      r2 = sqrt((xyz(1,j)-ex)**2+(xyz(2,j)-ey)**2+(xyz(3,j)-ez)**2)
      if (r2 .lt. rmin.and.j .ne. ia1) then
        rmin = r2
        exmin = ex
        eymin = ey
        ezmin = ez
      end if
    end do
    if (rmin .gt. rmax) then
      rmax = rmin
      xsave = exmin
      ysave = eymin
      zsave = ezmin
    end if
    icount = icount+1
    if (icount .lt. 100) goto 10

    ex = xsave
    ey = ysave
    ez = zsave

  end subroutine shiftlp

!========================================================================================!
!========================================================================================!
end module mo_localize
