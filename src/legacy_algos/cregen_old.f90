!===============================================================================!
! This file is part of crest.
!
! Copyright (C) 2018-2020 Stefan Grimme, Philipp Pracht
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
!===============================================================================!

!===============================================================================!
! THIS IS THE OLD VERSION OF CREGEN AND IS ONLY PRESENT FOR LEGACY SUPPORT
! SEE newcregen.f90 FOR THE ACTUAL ROUTINE USED PER DEFAULT
!===============================================================================!

subroutine cregen2(env)
      use iso_fortran_env, only : wp => real64, sp => real32, idp => int64
      use crest_data
      use crest_parameters, only: bohr,autokcal
      use ls_rmsd
      use iomod
      use strucrd, only: rdnat,rdcoord,wrc0,i2e,e2i
      use axis_module
      use miscdata, only: rcov
      use utilities
      use omp_lib
      implicit none
      type(systemdata) :: env    ! MAIN STORAGE OS SYSTEM DATA

      real(wp),allocatable :: xyz(:,:,:),e(:),p(:),er(:),c1(:,:),c2(:,:)
      real(wp),allocatable :: c1r(:,:),c2r(:,:)
      real(wp),allocatable :: xyznew(:,:,:),tmp(:,:),tmp3(:,:),enuc(:)
      real(wp),allocatable :: edum(:),dist(:,:,:),tmp2(:),cn(:),cn0(:) 
      real(wp),allocatable :: sd(:,:),rot(:,:),roth(:,:),metric(:,:),jfake(:)
      real(sp),allocatable :: rmat(:)
      integer,allocatable  :: at(:),ind(:),double(:),group(:),idum(:),ind2(:)
      integer,allocatable  :: atr(:)
      integer,allocatable  :: relat(:,:),equiv(:,:,:),molvec(:)
      integer,allocatable  :: pair(:),pre(:),nb(:,:),elist(:,:),flist(:,:)
      integer,allocatable  :: jnd(:),sames(:)
      logical,allocatable  :: vis(:)
      real(wp),allocatable :: xx(:),pg(:)
      real(wp),allocatable :: gdum(:,:),Udum(:,:),xdum(:), ydum(:)  ! rmsd dummy stuff
      integer,allocatable:: iref(:),dgen(:),glist(:,:),nmract(:)
      character(len=40),allocatable :: origin(:),originnew(:)
      integer,allocatable :: timetag(:)

      real(wp) :: beta,eav,A,T,rthr,s,g
      real(wp) :: ss,ethr,de,dr,bthr
      real(wp) :: elow,ewin,pthr,eref,athr,r
      real(wp) :: dum
      real(wp) :: rij(3),cnorm,erj
      !real*8 :: gdum(3,3),Udum(3,3),xdum(3), ydum(3)  ! rmsd dummy stuff
      real(wp), external :: rotdiff,shortest_distance

      integer :: i,j,k,l,m,nall,n,ig,ng,nall2
      integer :: m1,m2,s1,s2,maxg,current,ndoub
      integer :: molcount0,j1,iat,k2,m1end
      integer :: kk,memb,irr,ir,nr,jj,nall_old
      integer :: TID, nproc
      integer :: ich,ich2,io,ich3,rednat
      integer(idp) :: klong
      integer(idp),allocatable :: rmatmap1(:)
      integer,allocatable :: rmatmap2(:),includeRMSD(:)

      character(len=80) :: atmp,btmp,oname,fname,cname
      character(len=3) :: a3            
      character(len=512) :: outfile

      logical :: debug,heavy,l1,l2,ex,newfile,anal,rmsdchk
      logical :: fail
      logical :: substruc
      logical :: ttag
      logical, external :: distcheck,equalrot,equalrotall

      settingNames: associate( crestver => env%crestver, confgo => env%confgo, &
      &             methautocorr => env%methautocorr, printscoords => env%printscoords,  &
      &             ensemblename => env%ensemblename, compareens => env%compareens,      &
      &             elowest => env%elowest, ENSO => env%ENSO, subRMSD => env%subRMSD,    &
      &             trackorigin => env%trackorigin, autothreads => env%autothreads,      &
      &             omp => env%omp, MAXRUN => env%MAXRUN, threads => env%threads)
      rednat = env%rednat
      allocate(includeRMSD(env%nat))
      includeRMSD = env%includeRMSD

      thresholdNames: associate( cgf => env%cgf, thresholds => env%thresholds)

      outfile='cregen.out.tmp'
      if(cgf(6))outfile='tmp'

!------ the entire cregen output is printed to a seperate file
      call remove(outfile)
      if(confgo.and. .not.(env%properties.eq. -2))then
       ich=6
      else
        open(newunit=ich, file=outfile)
      endif      

!---- small header
      write(ich,'(''-------------------------------------'')')
      write(ich,'(''CREGEN - CONFORMER SYMMETRY ANALYSIS'')')
      write(ich,'(''-------------------------------------'')')


      call remove('LOWER_FOUND')
      
      allocate(gdum(3,3),Udum(3,3),xdum(3),ydum(3),nmract(100))

!---- setting the threads for OMP parallel usage
      if(autothreads)then
         call ompautoset(threads,4,omp,MAXRUN,0) !mode=4 --> Program intern Threads max
!$OMP PARALLEL PRIVATE(TID)
      TID = OMP_GET_THREAD_NUM()
      IF (TID .EQ. 0) THEN
         nproc = OMP_GET_NUM_THREADS()
         write(ich,*) '============================='
         write(ich,*) ' # threads =', nproc
         write(ich,*) '============================='
      END IF
!$OMP END PARALLEL 
      endif
!=====================================================================================!
      ewin = env%ewin      !ensemble energy window in kcal
      rthr = env%rthr      ! RMSD thr in Ang
      ethr = env%ethr      ! E threshold in kcal
      bthr = thresholds(4) ! rot const  thr  !legacy 
      athr = env%athr      ! to det. int. rotation. equal atoms for NMR, CRITICAL!
      pthr = env%pthr      ! population thr

      T    = env%tboltz  ! Temperature   

      debug  =cgf(1)
      newfile=cgf(2)
      anal   =cgf(3)
      heavy  =cgf(4) !compare only heavy atoms + OH with RMSD
      rmsdchk=cgf(5)

      end associate thresholdNames

      nmract = 0

      substruc=.false.

!======================================================================================!

      if(.not.rmsdchk) rthr = bthr  ! use rot thr for doubles check

      beta  =1./(T*8.314510/4.184/1000.+1.d-14)

      fname=trim(ensemblename)
      if(confgo.and.(index(trim(fname),'none selected').eq.0))then
        fname=trim(ensemblename)
        oname=trim(ensemblename)//'.sorted'
        write(ich,*) 'input  file name : ',trim(fname)
        write(ich,*) 'output file name : ',trim(oname)
        cname=ensemblefile !'crest_ensemble.xyz'
        if(env%fullcre)then
           ensemblename=trim(oname)
        endif
      else
        call checkname_xyz(crefile,fname,oname)
        write(ich,*) 'input  file name : ',trim(fname)
        write(ich,*) 'output file name : ',trim(oname)
        cname=conformerfile
      endif
      

      nall=1
      open(unit=1,file=fname)
 10   read(1,*,end=100) n
      read(1,'(a)') btmp
      do i=1,n
      read(1,'(a)') btmp
      enddo
      nall=nall+1
      goto 10
100   nall=nall-1
      if(nall.lt.1)then
         !stop 'read error'
         write(*,*) 'no new structures. continue.'
         return
      endif

!-----------------------------------------------------------------------------------------

      if(rednat.ne.n .and. subRMSD) substruc=.true.


      write(ich,'('' number of atoms                :'',i6)')n
      if(substruc)then
      write(ich,'('' atoms included in RMSD         :'',i6)')rednat
      endif
      write(ich,'('' number of points on xyz files  :'',i6)')nall
      write(ich,'('' RMSD threshold                 :'',f9.4)')rthr
      write(ich,'('' Bconst threshold               :'',f9.4)')bthr
      write(ich,'('' population threshold           :'',f9.4)')pthr
      write(ich,'('' conformer energy window  /kcal :'',f9.4)')ewin


      allocate(at(n),c1(3,n),c2(3,n),cn0(n),cn(n),molvec(n))
      allocate(c1r(3,rednat),c2r(3,rednat),atr(rednat))

      call rdcoord('coord',n,at,c1)
      call mrec(molcount0,c1,n,at,molvec)
      write(ich,'('' # fragment in coord            :'',i6)')molcount0

      if(substruc)then
        call cpincluded(n,rednat,c1,c1r,includeRMSD)
        k=1
        do i=1,n
          if(includeRMSD(i).gt.0)then
           atr(k)=at(i)
           k=k+1
          endif
        enddo 
        call ncoord(rednat,rcov,atr,c1r,cn0,500.0d0)
      else
        call ncoord(n,rcov,at,c1,cn0,500.0d0)
      endif

!c read file
      allocate(xx(10),iref(nall))

!--- determine the reliable points on the ensemble file
!    coord-file is the reference. Hence best use an optimized 
      rewind 1
      iref=0
!===================================================================!
      if(.not.substruc)then
!===================================================================!          
      do j=1,nall
         read(1,*) k
         read(1,'(a)') btmp
         call readl(btmp,xx,k)
         erj=xx(1)            ! from GMD ensemble i.e. in Eh
         do i=1,n
            read(1,*)a3,c2(1:3,i)
            if(j.eq.1) at(i) = e2i(a3)
         enddo
         c1(1:3,1:n)=c2(1:3,1:n)/bohr
         l1=distcheck(n,c1)  ! distance check
         cnorm=sum(abs(c1))
         if(abs(erj).gt.1.d-6.and.cnorm.gt.1.d-6.and.l1) &
      &  call ncoord(n,rcov,at,c1,cn,500.0d0)    ! further check for reactions based on CN
         k2=0        
         do i=1,n
            if(.not.(cn0(i).gt.0))cycle
            !if(abs((cn(i)-cn0(i)))/cn0(i).gt.0.3) k2=1       
         enddo
         if(k2.eq.1.or. &
         &  abs(erj).lt.1.d-6.or. &
         &  cnorm.lt.1.d-6.or. &
         &  (.not.l1) ) &
         &  then 
            write(ich,*) 'removing structure',j !,k2,er(j),cnorm,l1
            iref(j)=1
         endif
      enddo
!===================================================================!
      else  !  n /= rednat (only selected atoms contribute
!===================================================================!
      do j=1,nall
         read(1,*) k
         read(1,'(a)') btmp
         call readl(btmp,xx,k)
         erj=xx(1)            ! from GMD ensemble i.e. in Eh
         l=1
         do i=1,n
            if(includeRMSD(i).gt.0)then  
              read(1,*)a3,c2r(1:3,l)
              if(j.eq.1) atr(l) = e2i(a3)
              l=l+1
            else
              read(1,'(a)') btmp
            endif
         enddo
         c1r(1:3,1:rednat)=c2r(1:3,1:rednat)/bohr
         
         l1=distcheck(rednat,c1r)  ! distance check
         cnorm=sum(abs(c1r))
         if(abs(erj).gt.1.d-6.and.cnorm.gt.1.d-6.and.l1) &
      &  call ncoord(rednat,rcov,atr,c1r,cn,500.0d0)    ! further check for reactions based on CN
         k2=0
         do i=1,rednat
            !if(abs((cn(i)-cn0(i)))/cn0(i).gt.0.3) k2=1
         enddo
         if(k2.eq.1.or. &
         &  abs(erj).lt.1.d-6.or. &
         &  cnorm.lt.1.d-6.or. &
         &  (.not.l1) ) &
         &  then
            write(ich,*) 'removing structure',j !,k2,er(j),cnorm,l1
            iref(j)=1
         endif
       enddo
!==================================================================!
       endif
!==================================================================!       
    !---------
    ! iref is 0 for structures that are ok 
    ! and 1 for 'broken' structures.
    !--------
      nall_old=nall
      nall=nall-sum(iref(1:nall))
      nall2=nall
      write(ich,'('' number of reliable points      :'',i6)')nall    

      allocate(xyz(3,n,nall),er(nall),            &
     &         p(nall),e(nall),ind(nall),double(nall),tmp(3,nall),      &
     &         xyznew(3,n,nall),tmp2(n),jnd(n),sames(nall),ind2(nall),  &
     &         group(0:nall),idum(nall),edum(nall),                     &
     &         relat(0:n,n),enuc(nall),tmp3(3,nall),                    &
     &         jfake(n*(n+1)/2),sd(n,n),rot(3,nall),roth(3,nall),metric(n,n), &
     &         vis(n),pre(n),pair(n*(n+1)/2),nb(200,n),elist(n,n),flist(n,n))


      !array for tracking the origin
      if(trackorigin)then      
         allocate(origin(nall),originnew(nall))
         if(env%entropic)then
           allocate(timetag(nall))
         endif
      endif

!---- only use the reliable points
!!!
! THE ENSEMBLE FILE IS READ HERE
!!!
      jj=0
      rewind 1
      do j=1,nall_old
         if(iref(j).eq.0) jj=jj+1
         read(1,*) k
         read(1,'(a)') btmp
         call readl(btmp,xx,k)
         if(iref(j).eq.0)then
          er(jj)=xx(1)    
          if(trackorigin)then
             call getorigin(trim(btmp),origin(jj))
          endif
          do i=1,n
            read(1,*)a3,xyz(1:3,i,jj)
          enddo
         else
          do i=1,n
            read(1,*)a3,rij(1:3)
          enddo
         endif
      enddo
      close(1)

!---sort
!!!
! ind array is used to track the order
!!!
      do i=1,nall
         ind(i)=i
      enddo
      call Qsort(er,1,nall,ind) ! sort for E 
      elow=er(1)

      do i=1,nall
         xyznew(1:3,1:n,i)=xyz(1:3,1:n,ind(i))
      enddo

      do i=1,nall
         xyz(1:3,1:n,i)=xyznew(1:3,1:n,i)
      enddo
    
      !--- for the origin tracking. The arrays conatin a lable form which step the structure was obtained
      if(trackorigin)then
        do i=1,nall
         originnew(i)=origin(ind(i))
        enddo
        origin=originnew
      endif
      !if(trackorigin .and. env%entropic)then !check if we have valid time tags for entropy mode
      !  atmp=origin(1)
      !  if(atmp(1:1)=='t')then
      !      ttag=.true.
      !  else
            ttag=.false.
      !  endif
      !endif

!---to CMA
      do i=1,nall
         c1(1:3,1:n)=xyz(1:3,1:n,i)
         call axis(n,at,c1,c2,rot(1:3,i))   !-- this does the CMA trafo, and gives rot.const. MHz    
         !write(*,'(i5,2x,3f10.2)') i,rot(1:3,i) !-- rotational constant printout
         xyz(1:3,1:n,i)=c2(1:3,1:n)
      enddo
      if(substruc)then  !-- we have to do it again for the reduced mode, i.e., if atoms are excluded from the comparison
         do i=1,nall
            c1(1:3,1:n)=xyz(1:3,1:n,i)
            call cpincluded(n,rednat,c1,c1r,includeRMSD)
            call axis(rednat,atr,c1r,c2r,rot(1:3,i))
         enddo
      endif

      if(.not.env%allrot .and. .not.rmsdchk)then  !-- not needed by default
      rot = rot /dble(n)**2  !-- normalize such that it can be used instead of the RMSD
      endif

!--------
      if(env%esort) then    !-- just resort the input and make a cut (if inquired by cmd)
        open(unit=2,file=oname)                        
        do i=1,nall
           if((er(i)-elow)*autokcal.gt.ewin) goto 999
           c2(1:3,1:n)=xyz(1:3,1:n,i)
           write(2,'(2x,i0)') n
           if(.not.trackorigin)then
             write(2,'(2x,f18.8)') er(i)
           else
             write(2,'(2x,f18.8,4x,a)') er(i),'!'//trim(origin(i)) 
           endif
           do j=1,n
              write(2,'(1x,a2,1x,3f18.8)')i2e(at(j),'nc'),c2(1:3,j)
           enddo
        enddo
 999    close(2)
        write(ich,*) i-1,' out of ',nall,' taken'
        write(ich,*) 'normal resort termination of confg'
        goto 666
      endif
!--------
      inquire(file='.tmpxtbmodef',exist=ex)
      if(ex)then
       open(unit=66,file='.tmpxtbmodef')
       read(66,*) i
       read(66,*) eref
       close(66)
      else
       eref=er(1)
      endif
      if(crestver.eq.1)then
        write(ich,*)'reference state (Hessian) Etot :',eref
      else
        write(ich,*)'reference state Etot :',eref
      endif
!---get RMSD
  !========================================================!
  !--- crucial point: rmat is huge. VERY huge.
      !klong=nall
      !klong=klong*(nall+1)
      !klong=klong/2
      !write(*,*) nall,klong

  !--- for large ensembles the size of rmat can be largley reduced
  !    but this requires additional tracking and counting.
  !    It will be worth it, however.
      allocate(rmatmap1(nall),rmatmap2(nall))
      klong=0
      do i=1,nall
         rmatmap1(i)=klong
         l1=.true.
         do j=1,i-1
            !er(j) should always be smaller than er(i) because of presorting
            de=(er(i)-er(j))*autokcal

            if(de.lt.ethr)then
                klong=klong+1
                if(l1)then
                    rmatmap2(i)=j
                    l1=.false.
                endif
            endif
          enddo
      enddo   

      allocate(rmat(klong))
      rmat=0
  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
      if(rmsdchk)then
  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
      write(*,'(1x,a)') 'running RMSDs...'
      if(.not.substruc)then !regular case
      klong=0    
      do i=1,nall
         c1(1:3,1:n)=xyz(1:3,1:n,i)
!$OMP PARALLEL PRIVATE ( j,klong,c2,xdum,ydum,Udum,gdum,dum,de) &
!$OMP SHARED ( i,c1,rmat,n,xyz,rmatmap1,rmatmap2,er,ethr)
!$OMP DO 
         do j=1,i-1
           de=(er(i)-er(j))*autokcal
           if(de.lt.ethr)then
                c2(1:3,1:n)=xyz(1:3,1:n,j)
                call rmsd(n,c1,c2,0,Udum,xdum,ydum,dum,.false.,gdum) ! all atoms
               ! klong=lina(j,i)
                klong=linr(rmatmap1(i),rmatmap2(i),j)
                !$omp critical
                rmat(klong) = real(dum, 4)
                !$omp end critical
            endif
         enddo
!$OMP END DO
!$OMP END PARALLEL
      enddo
!--- substruc
      else  ! n /= rednat
      do i=1,nall
         c1(1:3,1:n)=xyz(1:3,1:n,i)
         call cpincluded(n,rednat,c1,c1r,includeRMSD)
!$OMP PARALLEL PRIVATE ( j,klong,c2,c2r,xdum,ydum,Udum,gdum,dum,de) &
!$OMP SHARED ( i,c1r,rmat,n,xyz,rednat,includeRMSD,rmatmap1,rmatmap2,er,ethr)
!$OMP DO 
         do j=1,i-1
           de=(er(i)-er(j))*autokcal
           if(de.lt.ethr)then
               c2(1:3,1:n)=xyz(1:3,1:n,j)
               call cpincluded(n,rednat,c2,c2r,includeRMSD)
               call rmsd(rednat,c1r,c2r,0,Udum,xdum,ydum,dum,.false.,gdum) ! all atoms
               ! klong=lina(j,i)
               klong=linr(rmatmap1(i),rmatmap2(i),j)
               !$omp critical
               rmat(klong) = real(dum, 4)
               !$omp end critical
           endif
         enddo
!$OMP END DO
!$OMP END PARALLEL
      enddo
      endif
      write(*,'(1x,a)') 'done.'
   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
      else   !if rotational constants are used instead of rmsd
   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
       do i=1,nall
         do j=1,i-1
            de=(er(i)-er(j))*autokcal
            if(de.lt.ethr)then
            !klong=lina(j,i)
            klong=linr(rmatmap1(i),rmatmap2(i),j)
            rmat(klong)=real(rotdiff(i,j,nall,rot), 4)
            endif
         enddo
       enddo
  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
      endif

!---apply RMSD threshold
! "sames" is used to track conformers (INCLUDING doubles).
!  IMPORTANT: we only have to check for duplicates within the boundaries of ETHR.
!  This allows for a much smaller RMAT
      double=0
      sames=0
      do i=1,nall
         do j=1,i-1
            ! klong=lina(j,i)
            de=(er(i)-er(j))*autokcal
            if(de.lt.ethr)then
              klong = linr(rmatmap1(i),rmatmap2(i),j) 
              dr = rmat(klong)
            else
              cycle  
            endif
            if(dr.lt.rthr) then
               double(i)=1  ! i > j, i.e., E(i) > E(j)
               sames(i)=j
               !write(*,*) i,j,rmat(k)
            elseif(dr.lt.2.*rthr) then
               if(.not.env%allrot)then
               l1=equalrot(i,j,nall,0.5*bthr,rot)
               else
               l1=equalrotall(i,j,nall,bthr,rot)
               endif
               if(l1)then !.and.l3) then
                  double(i)=1  
                  sames(i)=j
               endif
            endif
         enddo
      enddo
      deallocate(rmatmap2,rmatmap1)
      deallocate(rmat)

!---kick out doubles
      ndoub = 0
      do i=1,nall 
         if(double(i).eq.1) then
            er(i)=1.d+6           !with this high energy the conformer won't be considered later
            ndoub = ndoub + 1
         endif
      enddo
      write(ich,*)'number of doubles removed by rot/RMSD         :',ndoub

!---rel. energies    
      do i=2,nall
         er(i)=(er(i)-er(1))*autokcal
      enddo
      er(1)=0

!---keep only the true conformers in final ensemble
      edum=0
      tmp=0
      if(trackorigin) originnew=''
      m=0
      do i=1,nall
         if(er(i).gt.ewin) cycle   !--- discard everything >6kcal/mol
         m=m+1
         edum(m)=er(i)
         xyznew(1:3,1:n,m)=xyz(1:3,1:n,i)
         tmp   (1:3,m)    =rot(1:3,i)
         !tmp3  (1:3,m)    =roth(1:3,i)
         if(trackorigin)then
            originnew(m)=origin(i)
         endif
         ind2(m)=i
      enddo
      nall  = m
      xyz   = xyznew
      er    = edum
      rot   = tmp
      !roth  = tmp3
      write(ich,*)'total number unique points considered further :',nall 

      do k=1,nall
         c1(1:3,1:n)=xyz(1:3,1:n,k)
         enuc(k)=0                          
         do i=1,n-1
            do j=i+1,n
               r=(c1(1,i)-c1(1,j))**2 &
           &    +(c1(2,i)-c1(2,j))**2 &
           &    +(c1(3,i)-c1(3,j))**2+1.d-12
               enuc(k)=enuc(k)+at(i)*at(j)/r
            enddo
         enddo
      enddo

!--- IMPORTANT SORTING PROCEDURE
!    Loop over pre-sorted structure pairs, and check several things:
      double=0  !--- grouping tensor                                 
      SORTI : do i=1,nall
          SORTJ : do j=1,i 
            de=er(i)-er(j)                        !-- energy difference
            if(.not.env%allrot)then
              l1=equalrot(i,j,nall,bthr,rot)
            else
              l1=equalrotall(i,j,nall,bthr,rot)     !-- rotational constant difference
            endif

            l2=2.*abs(enuc(i)-enuc(j))/(enuc(i)+enuc(j)).lt.1.d-3  !-- nuclear permutation
!           k=lin(i,j)

            if(l1.and.l2.and.i.ne.j  &
           &   .and.double(i).eq.0   &
           &   .and.abs(de).lt.ethr) then
               double(i)=j   !-- this means "structure i is the same conformer as structure j"
                             !   for j=i it is 0 because the if statement is skipped.
               cycle SORTI
               goto 99
            endif
          enddo SORTJ
 99       continue
      enddo SORTI

       do i=1,nall
          call backtrack(double,i,j)
       enddo

!---assign structures to identical conformers
      ig=0      !-- total number of "groups", i.e. unique conformers
      group=0   !-- grouping tensor -- which structure belongs to which group
   !-- first, assign a representative for each new group, and determine how many groups there are
      do i=1,nall
         if(double(i).eq.0)then
            ig=ig+1                            !-- ascending order 1-to-"total number of groups"
            iref(ig)=i                         !-- i=the unique ref.structure in each group
            group(i)=ig                        !-- assign struc. i to group ig
         endif
      enddo
      ng=ig      !-- number of conformer groups, i.e., unique conformers
    !-- then, assign all remaining structures to their respective groups
      !!THE FOLLOWING PART IS WRONG!
      !k=group(1)
      !do i=1,nall
      !   if(group(i).ne.0)then
      !      k=group(i)
      !   else
      !      group(i)=k
      !   endif
      !enddo
      !!THIS SHOULD BE CORRECT!
      do i=1,nall
         if(group(i).eq.0)then
           j=double(i)
           group(i)=group(j)
         endif 
      enddo

      !-- repair the group order (can be importent when conformers with very similar energies exist)
      if(trackorigin .and. env%entropic .and. ttag)then
        call repairentropic(nall,n,ng,group,xyz,er,trackorigin,originnew)
        originnew(1) = 't1' !dirty hack, but conceptionally this is rigth -> lowest conformer must be present at the start
      else
        call repairorder(nall,n,ng,group,xyz,er,trackorigin,originnew)
      endif
     
      !-- determine # of members in each group (=rotamers)
      allocate(dgen(nall))
      dgen=0
      do i=1,nall
         dgen(group(i))=dgen(group(i))+1
      enddo

      allocate(glist(0:nall,ng))
      glist=0
      do i=1,ng                                ! det members of indents in each group
         do j=1,nall
            if(group(j).eq.i)then
               glist(0,i)=glist(0,i)+1
               glist(glist(0,i),i)=j
            endif
         enddo
      enddo

!----
! ADVANCED TRACKING FOR ENSO
!---
      if(ENSO.and.confgo)then
       do i=1,nall
         j=ind2(i)
         !if(i.eq.1.or.group(i-1).ne.group(i))then
         if(i.ne.1 .and. group(i-1).eq.group(i))then         
         k=ind2(i-1)
         sames(i)=k
         if(k.eq.i)sames(i)=0
         endif
       enddo
       do i=1,nall2
          j=0
          call backtrack(sames,i,j)
       enddo
       j=sum(sames)
       open(newunit=ich3,file='cregen.enso')
       if(j.gt.0)then  !doubles found
         write(ich3,*) "DUPLICATES FOUND"
         do i=1,nall2
            if(sames(i).gt.0)then
            !write(*,*) ind(i),ind(sames(i))
             write(ich3,*)ind(i),ind(sames(i))
            endif
         enddo
       else
         write(ich3,*) "ALL UNIQUE"
       endif
       close(ich3)
      endif


!---distance neighbor list
      c2(1:3,1:n)=xyz(1:3,1:n,1)/bohr
      call neighdist(n,at,c2,nb,metric)
      k=0
      pair=0
      do i=1,n-1
         do j=i+1,n
!---the shortest bond path
            current=j
            dum=shortest_distance(n, i, j, nb, metric, vis, pre)
            k=0
            do while (pre(current) /= 0)
               current = pre(current)
               k = k + 1
            end do ! End loop: while precessor(current) /= 0
            pair(lin(j,i))=k  ! # of bonds between i and j
         enddo
      enddo

      allocate(equiv(0:n,n,0:nall))
      equiv=0
! (costly) symmetry analyis of all rotamers for NMR. this is complicated stuff also
! and the end of the program where this is completed...
      if(anal)then
      allocate(dist(n,n,nall))
      do i=1,nall
         call distance(n,xyz(:,:,i),dist(:,:,i))   ! distance matrix
         do j=1,n
            do k=1,n
               tmp2(k)=dist(k,j,i)*dble(at(k))  ! the distance of j to all atoms * Z to distinguish
            enddo
            call qqsort(tmp2,1,n)
            dist(1:n,j,i)=tmp2(1:n)
         enddo
      enddo
      write(ich,*) 'compare ...'
      do i=1,ng 
         m=glist(0,i)
         if(m.lt.2) cycle                      ! det equivalent atoms in each group
         relat=0
         do m1=1,m
!$OMP PARALLEL PRIVATE ( m2, s1, s2 ) SHARED ( relat )
!$OMP DO 
            do m2=1,m1-1                       ! compare all members 
               s1=glist(m1,i)                  ! struc 1
               s2=glist(m2,i)                  ! struc 2
               call compare(n,nall,s1,s2,dist,athr,relat) ! athr is distance vector equivalence threshold
            enddo
!$OMP END DO
!$OMP END PARALLEL
         enddo
         equiv(0:n,1:n,i)=relat(0:n,1:n)
      enddo
      deallocate(dist)
! symmetrize result i.e. if iat is in list of jat, jat must be in list of iat
! done again at the end of this part
      do i=1,ng 
         do j1=1,n
            m1=equiv(0,j1,i)
            do k=1,m1
               iat=equiv(k,j1,i)
!              is atom j1 in the list of atom iat?
               ex=.false.
               m2=equiv(0,iat,i)
               do k2=1,m2
                  if(j1.eq.equiv(k2,iat,i)) ex=.true.
               enddo
               if(.not.ex)then
                  equiv(0,iat,i)=equiv(0,iat,i)+1
                  equiv(equiv(0,iat,i),iat,i)=j1
               endif
            enddo
         enddo
      enddo

      do i=1,ng 
        if(sum(equiv(0,1:n,i)).eq.0) cycle     
        if(debug)write(ich,*)'conformer set ',i
        memb=glist(0,i)
        do j=1,n
           m=equiv(0,j,i)
           if(m.eq.0) cycle         
           if(debug)write(ich,'(''atom :'',i4,'' equivalent to '',40i3)') &
        &  j,(equiv(k,j,i),k=1,m)
        enddo
      enddo
      endif ! symmetry analysis

!---populations
      call boltz(nall,T,er,p)
      allocate(pg(ng))
      do i=1,ng 
         m=glist(0,i)
         pg(i)=0
         do j=1,m
            pg(i)=pg(i)+p(glist(j,i))  !-- sum of pops in a group
         enddo
         if(pg(i).gt.pthr) maxg=i
      enddo

!---output and write
      c1(1:3,1:n)=xyz(1:3,1:n,1)
!      write(ich,*) 'full ensemble on file crest_rotamers_*.xyz'
!      write(ich,*) 
!     .'unique structures in ensemble on crest_conformers.xyz and scoord.*'
!      write(ich,*)'note: enantiomers are included.'
      if(.not.trackorigin)then
        write(ich,*)'  Erel/kcal     Etot      weight/tot conformer  set degen'
        do i=1,nall
           if(i.eq.1.or.group(i-1).ne.group(i))then
           write(ich,'(i5,f8.3,1x,3f11.5,2i5)') &
     &     i,er(i),er(i)/autokcal+elow,p(i),pg(group(i)), &
     &     group(i),dgen(group(i))
           else
           write(ich,'(i5,f8.3,1x,2f11.5,i5)')i,er(i),er(i)/autokcal+elow,p(i)
           endif
        enddo
      else !trackorigin
        write(ich,*)'  Erel/kcal     Etot      weight/tot conformer  set degen    origin'
        do i=1,nall
           if(i.eq.1.or.group(i-1).ne.group(i))then
           write(ich,'(i5,f8.3,1x,3f11.5,2i5,5x,a)') &
     &     i,er(i),er(i)/autokcal+elow,p(i),pg(group(i)), &
     &     group(i),dgen(group(i)),trim(originnew(i))
           else
           write(ich,'(i5,f8.3,1x,2f11.5,26x,a)') &
     &     i,er(i),er(i)/autokcal+elow,p(i),trim(originnew(i))
           endif
        enddo
      endif
!---- write a mapping file file for "-compare" mode
      if(compareens)then
         open(newunit=ich2,file='.cretrack')
         write(ich2,'(5x,i0)')maxval(group)
         kk=1
         do i=2,nall
            if(group(i-1).ne.group(i))then
              j=i-1
              write(ich2,'(1x,i5,1x,i7,1x,i7)')group(j),kk,j    
              kk=i
            endif
            if(i.eq.nall)then
              write(ich2,'(1x,i5,1x,i7,1x,i7)')group(nall),kk,nall
            endif
         enddo
         close(ich2)
      endif
!---write part
      if(newfile) then
        open(unit=2,file=oname)
        do i=1,nall
           c2(1:3,1:n)=xyz(1:3,1:n,i)
           call axis(n,at,xyz(:,:,i),c2,xdum)           
           !--- try to align all conformers the same way
           call xyzalign(n,c2)
           write(2,*) n
           if(.not.trackorigin)then
           write(2,*) er(i)/autokcal+elow,p(i)
           else
           write(2,*) er(i)/autokcal+elow,p(i),'!'//trim(originnew(i))
           endif
           do j=1,n
              write(2,'(1x,a2,1x,3f18.8)')i2e(at(j),'nc'),c2(1:3,j)
           enddo
        enddo
        close(2)
      endif

!---entropy and G
      A=0
      eav=0
      do i=1,nall
         A=A+p(i)*log(p(i)+1.d-12)   
         eav=eav+p(i)*er(i)     
      enddo

      g = (1./beta)*A
      s =-1000.*4.184*g/T
      ss=-1000.*      g/T

      write(ich,'(''T /K                                  :'', F9.2)')T 
      write(ich,'(''E lowest                              :'',f12.5)')elow
      elowest=elow  ! <--- global access to elow
      !---- elow printout in between routines
      if(.not.confgo)then
      write(*  ,'(1x,''E lowest :'',f12.5)')elow
      endif
      write(ich,'(''ensemble average energy (kcal)        :'', F9.3)')eav
      write(ich,'(''ensemble entropy (J/mol K, cal/mol K) :'',2F9.3)')s,ss
      write(ich,'(''ensemble free energy (kcal/mol)       : '',F8.3)')g
      write(ich,'(''population of lowest in %             : '',F8.3)')pg(1)*100.d0

      if((crestver.eq.1).and.(.not.confgo))then
      if((elow-eref)*autokcal.lt.-0.2)then
         write(ich,*) '...............................................'
         write(ich,*) 'WARNING: new (best) energy less than that from '
         write(ich,*) 'WARNING: preceding Hessian calculation:  '
         write(ich,*) 'Improved by ', eref-elow, ' Eh or ',(eref-elow)*autokcal, 'kcal'
         write(ich,*) '...............................................'
         call touch('LOWER_FOUND')
      endif
      endif

! now write the unique structures
! and the symmetry info for NMR
      maxg=0
      do i=1,ng
         k=glist(1,i) ! take first 
         if(er(k).lt.ewin) maxg=maxg+1
      enddo
      write(ich,*)'number of unique conformers for further calc ',maxg
      write(ich,*)'list of relative energies saved as "crest.energies"'
      open(newunit=ich2,file='crest.energies')
      do i=1,maxg
         write(ich2,'(i4,f12.3)') i,er(glist(1,i))
      enddo
      close(ich2)




!---NMR part
      if(anal)then
   !--- get NMR active nuclei
      atmp='.anmrrc'  ! <--- name of the .anmrrc written by ENSO
      call getanmrrc(atmp,fail)
      if(fail)then  !--- there is no .anmrrc from ENSO
        write(ich,*)'NMR mode.'
        nmract=0       ! all nuclei inactive
        nmract(1) = 1  ! H active
       !nmract(6) = 1  ! C active
        nmract(9) = 1  ! F active
        nmract(15)= 1  ! P active
      else          !--- there IS a .anmrrc, and it is used.
         write(ich,*)'NMR mode. Reading <',trim(atmp),'> for atomic NMR data'
         open(newunit=ich2,file=atmp)
           read(ich2,'(a)')atmp
           read(ich2,'(a)')atmp
           if(index(atmp,'ENSO').ne.0)then
              read(ich2,'(a)')atmp
           endif
           do
             read(ich2,*,iostat=io)i,xx(1:2),nmract(i)
             if(io<0) exit
           enddo
         close(ich2)
      endif


! inlcude equivalence info from the other conformers as well i.e.
! assume that all conformers have the same chemical equivalencies
! the result is put into equiv(:,:,0)
      equiv(0:n,1:n,0)=equiv(0:n,1:n,1)
      do i=2,ng
         do j=1,n
            m1end=equiv(0,j,0)   ! end of list of lowest
            do m=1,equiv(0,j,i)  ! list of higher
               k=equiv(m,j,i)    ! in the one in the higher list
               do m1=1,m1end
                  if(equiv(m1,j,0).eq.k) goto 19 !already there?
               enddo
               equiv(0,j,0)=equiv(0,j,0)+1 ! append
               equiv(equiv(0,j,0),j,0)=k
 19            continue
            enddo
         enddo
      enddo

!cc    this cc commented out code makes equivalencies for all confs seperately
!cc    do ig=0,ng     ! all conf groups
!cc    call mkd(ig)   ! make the NMR dir 
!cc    call getname(ig,'anmr_nucinfo',fname)  ! open the file 
!cc    open(unit=3,file=fname)
      ig=0
      atmp='anmr_nucinfo'
      open(unit=3,file=atmp)
!cc    write(*,'(''::::::::::: conformer group '',i3,'':::::::::::'')')ig
      write(ich,'(''::::::::::: conformer group all :::::::::::'')')
      write(3,*) n
!cccccccccccccccccc
!c chem eq. first
!cccccccccccccccccc
      elist=0
      !write(*,*) equiv(0,:,ig)
      if(methautocorr)then
        call methyl_autocomplete(n,c1/bohr,at,equiv(:,:,ig))
      endif
      !write(*,*) equiv(0,:,ig)
      do i=1,n
         m=equiv(0,i,ig)   
         do k=1,m              
            l=equiv(k,i,ig)
            elist(l,i)=1
         enddo
         elist(i,i)=1
      enddo
      do i=1,n
         do j=1,equiv(0,i,ig)
            k=equiv(j,i,ig)
            elist(1:n,k)=elist(1:n,k)+elist(1:n,i)
         enddo
      enddo
!---  prepare write out
      do i=1,n
         k=1
         equiv(1,i,ig)=i
         elist(i,i)=0
         do j=1,n
            if(elist(j,i).ne.0)then
               k=k+1
               equiv(k,i,ig)=j               
            endif
         enddo
         equiv(0,i,ig)=k
      enddo
      write(ich,*)'chemical equivalencies (mag.active nuclei):'
      jnd=1
      do j=1,n
         m=equiv(0,j,ig)
         write(3,'(3x,i0,3x,i0)') j, m
         do l=1,m
          if(l.ne.m)then
           write(3,'(1x,i0)',advance='no') equiv(l,j,ig)  ! include the atom ie if there are no equiv.
          else
           write(3,'(1x,i0)',advance='yes') equiv(l,j,ig)
          endif
         enddo
         if(nmract(at(j)).eq.0) cycle
         if(m.gt.1.and.jnd(j).eq.1)then  ! just print
         write(ich,'(''reference atom'',i4,'' # :'',i2)') equiv(1,j,ig),m
         do k=1,m
            jnd(equiv(k,j,ig))=0
         enddo
         endif
      enddo
!cccccccccccccccccc
!c mag eq. 
!cccccccccccccccccc
!c make a check list of atoms for the mag. eq.
      elist=0
      flist=1
      do i=1,n
         m=equiv(0,i,ig)        ! the following lines fill the equiv list  
         do k=1,m              
            l=equiv(k,i,ig)
            elist(l,i)=1
         enddo
         elist(i,i)=1
      enddo
      flist=elist

      do i=1,n    
         m=equiv(0,i,ig)
         do k=1,m
            l=equiv(k,i,ig)
            if(l.eq.i) cycle
            do j=1,n
               if(flist(j,i).eq.1.or.nmract(at(j)).eq.0) cycle   ! don't check non-magnetic nuclei
!c              write(*,*) l,j,pair(lin(i,j)),pair(lin(l,j))      ! and chem. equiv. ones (ie in the same
               if(pair(lin(i,j)).ne.pair(lin(l,j))) elist(l,i)=0 ! group
            enddo
         enddo
      enddo
!---  symmetrize 
      do i=1,n
         k=1
         equiv(1,i,ig)=i
         elist(i,i)=0
         do j=1,n
            if(elist(j,i).ne.0)then
               k=k+1
               equiv(k,i,ig)=j               
            endif
         enddo
         equiv(0,i,ig)=k
      enddo
      do i=1,n
         do j=1,equiv(0,i,ig)
            k=equiv(j,i,ig)
            elist(1:n,k)=elist(1:n,k)+elist(1:n,i)
         enddo
      enddo
!---  prepare write out
      do i=1,n
         k=1
         equiv(1,i,ig)=i
         elist(i,i)=0
         do j=1,n
            if(elist(j,i).ne.0)then
               k=k+1
               equiv(k,i,ig)=j               
            endif
         enddo
!old     equiv(0,i,ig)=k
         if(k.gt.2) then
            equiv(0,i,ig)=k    ! CH3 etc
         else
            equiv(0,i,ig)=1    ! this makes CH2-CH2 not mag. equiv. !
         endif
      enddo
      jnd=1
      write(ich,*)'magnetic equivalencies:'
      do j=1,n
         m=equiv(0,j,ig)
         write(3,*) j, m
         write(3,'(20i5)') (equiv(l,j,ig),l=1,m)  ! include the atom ie if there are no equiv.
         if(nmract(at(j)).eq.0) cycle
         if(m.gt.1.and.jnd(j).eq.1)then  ! just print
         write(ich,'(''reference atom'',i4,'' # :'',i2)') equiv(1,j,ig),m
         do k=1,m
            jnd(equiv(k,j,ig))=0
         enddo
         endif
      enddo
      close(3)

!cc    enddo  ! next conformer group (ig)


!ccccccccccccccccccccc
!c J averaging matrix
!ccccccccccccccccccccc
 
      fname='anmr_rotamer'
      open(unit=112,file=fname,form='unformatted')
      !open(unit=112,file=fname)
      write(112) ng 

      jfake=0
      do ig=1,ng     ! all conf groups
      nr=glist(0,ig) ! how many rotamers?
      write(112) nr         
      do ir=1,nr
         irr=glist(ir,ig)
         call distance(n,xyz(:,:,irr),sd)   ! distance matrix
         c1(1:3,1:n)=xyz(1:3,1:n,irr)
         call ncoord(n,rcov,at,c1,cn,500.0d0)
         do i=1,n-1
            do j=i+1,n
               jfake(lin(j,i))=cn(i)*cn(j)*sqrt(dble(at(i)*at(j))) &
          &    /(dble(pair(lin(j,i)))*sd(j,i)**5) ! the approx. "J" is topologically equivalent to J
                                                  ! R^3 was wrong in one case because Hs were artificially paired
                                                  ! R^5 seems to be save
            enddo
         enddo
!c        call prmat(6,dd,n,n,'J fake') 
         write(112) jfake(1:n*(n+1)/2)        ! read by anmr 
      enddo
      enddo
      close(112)

   !--- how many rotamers per conformer   
      open(newunit=ich3,file='cre_degen')
      write(ich3,'(3x,i0)') ng
      do ir=1,ng
        write(ich3,'(3x,i0,2x,i0)') ir,glist(0,ir)
      enddo
      close(ich3) 

      !if(env%entropic)then
      ! if(ttag)then    
      !  call origin2time(nall,originnew,timetag)
      ! else
      !  timetag=1
      ! endif   
      !open(newunit=ich3,file='cre_tag')
      !if(env%mdtime.lt.0.0_wp)then
      !  write(ich3,'(3x,i0,1x,f10.2)') ng,1.0_wp
      !else
      !  write(ich3,'(3x,i0,1x,f10.2)') ng,env%mdtime
      !endif
      !do ir=1,ng
      !  k=glist(1,ir)
      !  write(ich3,'(3x,i0,2x,i0)') ir,timetag(k)
      !enddo
      !close(ich3)
      !endif    

      endif ! equivalence anal
!===============================================================!
      call remove('cre_members')
      open(newunit=ich3,file='cre_members')
      write(ich3,'(3x,i0)') ng
      do i=1,ng
        k=glist(0,i)
        write(ich3,'(3x,i5,1x,i10,1x,i10)') &
      &   k,glist(1,i),glist(k,i)
      enddo
      close(ich3)

      open(unit=2,file=trim(cname))
      open(unit=4,file='crest_best.xyz')
      do i=1,ng  
         k=glist(1,i) ! take first 
         if(er(k).ge.ewin) cycle
         if(i.eq.1.or.printscoords)then   !write a scoord.* for each conformer? ,scoord.1 is always written  
           call getname1(i,atmp)
           c2(1:3,1:n)=xyz(1:3,1:n,k)/bohr
           call wrc0(atmp,n,at,c2)
         endif
         write(2,'(2x,i0)') n
         write(2,'(2x,f18.8)') er(k)/autokcal+elow
         do j=1,n
            write(2,'(1x,a2,1x,3f20.10)')i2e(at(j),'nc'),xyz(1:3,j,k)
         enddo
         if(i.eq.1)then
         write(4,'(2x,i0)') n
         write(4,'(2x,f18.8)') er(k)/autokcal+elow
         do j=1,n
            write(4,'(1x,a2,1x,3f20.10)')i2e(at(j),'nc'),xyz(1:3,j,k)
         enddo
         endif
      enddo
      close(2)
      close(4)

      !because of "esort"-flag these have to be deallocated first
      deallocate(pg,glist,dgen)

 666  write(ich,*) 'Normal termination.'

      if(.not.confgo) close(ich) !close the output file
!deallocate all the rest
      if(trackorigin)then
         deallocate(originnew,origin)
      endif
      deallocate(flist,elist,nb,pair,pre,vis,metric,rot)
      deallocate(sd,jfake,enuc,equiv,relat,edum)
      deallocate(idum,group,jnd,tmp2,xyznew,tmp,double)
      deallocate(ind,e,p,er,xyz,ind2,sames)

      deallocate(iref,xx)

      deallocate(molvec,cn,cn0,c2,c1,at)
      deallocate(nmract,ydum,xdum,Udum,gdum)

      end associate settingNames
   
end subroutine cregen2

!================================================
! cpincluded is a small routine to map only
! the selected atoms from structure c into
! structure cr, via the array inc (value 1 or 0)
!================================================
subroutine cpincluded(n,nr,c,cr,inc)
      use iso_fortran_env, only : wp => real64
      implicit none
      integer,intent(in) :: n
      integer,intent(in) :: nr
      real(wp),intent(in) :: c(3,n)
      real(wp),intent(out) :: cr(3,nr)
      integer,intent(in) :: inc(n)
      integer :: i,k
      k=1
      do i=1,n
        if(inc(i).gt.0)then
        cr(1:3,k)=c(1:3,i)
        k=k+1
        endif
      enddo
      return
end subroutine cpincluded

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


subroutine compare(n,nall,s1,s2,dist,athr,relat)
      implicit none
      integer :: n,nall,s1,s2
      integer :: relat(0:n,n)
      real*8 :: dist(n,n,nall),athr

      integer :: i,j,k
      real*8 :: rsum,athr2
      logical :: ok

      athr2=athr*dble(n)
      do i=1,n
         do j=1,n
            if(i.eq.j) cycle
            rsum=0
!c           bp=.true.
            do k=1,n
               rsum=rsum+abs(dist(k,i,s1)-dist(k,j,s2))
!c              if(k.ne.i.and.k.ne.j.and.
!c    .            pair(lin(k,i)).ne.pair(lin(k,j))) then
!c                 bp=.false. ! check bond path length
!c                 exit
!c              endif
               if(rsum.gt.athr2) goto 99
            enddo
!c           if(rsum/n.lt.athr.and.bp) then
            if(rsum/n.lt.athr       ) then
               ok=.true.
               do k=1,relat(0,i)
                  if(relat(k,i).eq.j) ok=.false.   ! already existing?
               enddo
               if(ok)then
                  relat(0,i)=relat(0,i)+1
                  relat(relat(0,i),i)=j
               endif
            endif
 99         continue
         enddo
      enddo

end subroutine compare

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

logical function distcheck(n,xyz)
      implicit none
      real*8,allocatable :: rij(:)
      integer :: n
      real*8 :: xyz(3,n)
      integer i,j
      distcheck=.true.
      allocate(rij(3))
      do i=1,n-1
         do j=i+1,n
         rij=xyz(:,j)-xyz(:,i)
         if(sum(rij*rij).lt.1.d-3) distcheck=.false.
         enddo
      enddo
      deallocate(rij)
      return
end function distcheck


!==========================================================!
! bactrack is used to find the element "val" of array "arr"
! that contains 0 when the content of each field of arr
! references another field of it, starting from element i
!==========================================================!
recursive subroutine backtrack(arr,i,val)
      implicit none
      integer,intent(in) :: i
      integer,intent(inout) :: arr(*)
      integer,intent(out) :: val

      if(arr(i).gt.0)then
         call backtrack(arr,arr(i),val)
         arr(i)=val
      else
         val=i
      endif

      return
end subroutine backtrack


!==========================================================!
! align a structure in space according to the position
! of the first atom
!==========================================================!
subroutine xyzalign(nat,xyz)
      use iso_fortran_env, wp => real64
      use geo, only: rotRz180,rotRx180
      implicit none
      integer :: nat
      real(wp) :: xyz(3,nat)
      integer :: i
      if(xyz(1,1).lt.0.0_wp)then
       do i=1,nat
         call rotRz180(xyz(1:3,i))
       enddo
      endif
      if(xyz(2,1).lt.0.0_wp)then
       do i=1,nat
         call rotRx180(xyz(1:3,i))
       enddo
      endif
      return
end subroutine xyzalign

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine repairorder(nall,n,ng,order,xyz,er,trackorigin,originnew)
      use iso_fortran_env, wp => real64
      implicit none
      integer :: nall
      integer :: n,ng
      integer :: order(0:nall)
      real(wp) :: xyz(3,n,nall)
      real(wp) :: er(nall)
      logical :: trackorigin
      character(len=40) :: originnew(nall)
      integer :: i,j,k

      integer,allocatable :: dumorder(:)
      real(wp),allocatable :: dumxyz(:,:,:)
      real(wp),allocatable :: dumer(:)
      character(len=40),allocatable :: dumorig(:)

      allocate(dumorder(nall),dumxyz(3,n,nall),dumer(nall))
      if(trackorigin)then
      allocate(dumorig(nall))
      endif
      k=0
      do i=1,ng
         do j=1,nall  
            if(order(j).eq.i)then
            k=k+1     
            dumorder(k)=order(j)
            dumxyz(:,:,k)=xyz(:,:,j)
            dumer(k)=er(j)
            if(trackorigin)then
            dumorig(k)=originnew(j)
            endif
            endif
         enddo
      enddo
      order(1:nall)=dumorder(1:nall)
      xyz=dumxyz
      er=dumer
      if(trackorigin)then
      originnew=dumorig
      endif 

      if(trackorigin)then
      deallocate(dumorig)
      endif
      deallocate(dumer,dumxyz,dumorder)

      return
end subroutine repairorder

!===============================================================================
subroutine repairentropic(nall,n,ng,order,xyz,er,trackorigin,originnew)
      use iso_fortran_env, wp => real64
      implicit none
      integer :: nall
      integer :: n,ng
      integer :: order(0:nall)
      real(wp) :: xyz(3,n,nall)
      real(wp) :: er(nall)
      logical :: trackorigin
      character(len=40) :: originnew(nall)
      integer :: i,j,k,t
      integer :: tmax

      integer,allocatable :: dumorder(:)
      real(wp),allocatable :: dumxyz(:,:,:)
      real(wp),allocatable :: dumer(:)
      character(len=40),allocatable :: dumorig(:)
      integer,allocatable :: timetag(:)

      allocate(dumorder(nall),dumxyz(3,n,nall),dumer(nall))
      allocate(dumorig(nall))
      allocate(timetag(nall), source=0)
      if(trackorigin)then
      call origin2time(nall,originnew,timetag)
      endif
      tmax=maxval(timetag,1)
      k=0
      ILOOP : do i=1,ng
        TLOOP : do t=1,tmax
         JLOOP : do j=1,nall
            if(order(j).eq.i .and. timetag(j).eq.t)then
              k=k+1
              dumorder(k)=order(j)
              dumxyz(:,:,k)=xyz(:,:,j)
              dumer(k)=er(j)
              if(trackorigin)then
              dumorig(k)=originnew(j)
              endif
            else
              cycle JLOOP
            endif
         enddo JLOOP
        enddo TLOOP
      enddo ILOOP
      order(1:nall)=dumorder(1:nall)
      xyz=dumxyz
      er=dumer
      if(trackorigin)then
      originnew=dumorig
      endif

      if(trackorigin)then
      deallocate(timetag)
      deallocate(dumorig)
      endif
      deallocate(dumer,dumxyz,dumorder)

      return
end subroutine repairentropic



