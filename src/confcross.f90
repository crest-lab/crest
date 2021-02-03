!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2018-2020 Stefan Grimme
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

subroutine confcross(env,maxgen,kk)

      use iso_fortran_env, only : wp => real64
      use crest_data
      use ls_rmsd
      use iomod
      use strucrd, only: rdnat,rdcoord,wrxyz,e2i
      implicit none

      type(systemdata) :: env
      !type(options)    :: opt

      real(wp),allocatable::xyz(:,:,:),new(:,:),er(:),c0(:,:),c1(:,:),tmp3(:,:)
      real(wp),allocatable::coord(:,:),tmp(:,:,:),cn(:),c2(:,:),cnref(:)
      integer,allocatable::at(:),ind(:),na(:),nb(:),nc(:)
      integer,allocatable::at2(:),clist(:),cstat(:)
      real(wp),allocatable :: xx(:),rcov(:),rms(:)
      real(wp),allocatable :: gdum(:,:),Udum(:,:),xdum(:),ydum(:)  !rmsd dummy stuff
      real*4  ,allocatable :: csave(:,:,:),cint(:,:,:) ! save some mem because list can be long

      real(wp) :: autokcal,cthr,rthr,ewin,de,dr,bthr,ddum,rmsdav
      real(wp) :: elow,dum,aa,bb,cc,eref,athr,pi,sd,rthrbohr
      real(wp) :: cnorm
      real(wp) :: thresholds(8)
      real(wp) :: percent,pcount

      parameter (autokcal=627.509541d0)
      parameter (pi =  3.14159265358979D0)

      integer :: i,j,k,l,m,nall,n,ig,ng,m1,ierr,nmax,nl,kl,kk
      integer :: maxgen,r,n2,nk,ident,nmaxmax,ii
      integer :: io,ich
      integer TID, OMP_GET_NUM_THREADS, OMP_GET_THREAD_NUM, nproc
      integer, external :: lin

      character(len=512) :: thispath
      character(len=80) :: atmp,btmp,ctmp,fname,oname
      character(len=52) :: bar
      character(len=3) :: a3            

      logical debug,error,ldum,ex,fail
      
      associate( thresholds => env%thresholds, cgf => env%cgf)

!---- setting the threads for OMP parallel usage
      if(env%autothreads)then
         call ompautoset(env%threads,4,env%omp,env%MAXRUN,0) !mode=4 --> Program intern Threads max
         call ompprint_intern
      endif

      
      allocate(xx(10),rcov(94))
      allocate(gdum(3,3),Udum(3,3),xdum(3),ydum(3))

      call setrcov(rcov) ! for CN calc

      call checkname_xyz(crefile,fname,oname)
      write(*,*) 'input  file name : ',trim(fname)
!------------------------------------------------------------------------
      !thresholds
      !ewin=thresholds(1)        ! E threshold in kcal
      ewin=env%ewin
      rthr=env%rthr*2.0d0  ! standard RMSD threshold

!SG   cthr=0.5  ! CN clash threshold
      cthr=0.3  ! CN clash threshold SG change, was 0.2
!-----------------------------------------------------------------------

      nall=1
      open(unit=1,file=fname)
      do
        read(1,*,iostat=io) n
        if(io<0) exit
        read(1,'(a)') btmp
        cnorm=0
        do i=1,n
           read(1,'(a)') btmp
           call readl(btmp,xx,k)
           cnorm=cnorm+sum(abs(xx(1:3)))
        enddo
        if(cnorm.gt.1.d-6) nall=nall+1
      enddo
      nall=nall-1
      if(nall.lt.1) error stop 'read error'
      close(1)

      write(*,'(a,i6)')'number of atoms                :',n   
      write(*,'(a,i6)')'number of points on xyz files  :',nall     
      write(*,'(a,f8.2)')'conformer energy window  /kcal :',ewin     
      write(*,'(a,f8.4)')'CN per atom difference cut-off :',cthr     
      write(*,'(a,f8.4)')'RMSD threshold                 :',rthr     
      write(*,'(a,1x,i0)')'max. # of generated structures :',maxgen   
      rthrbohr=rthr/bohr

      if(nall.lt.2)then
         write(*,*) 'Not enough structures to perform GC!'
         kk=0
         return
         goto 666
      endif

      allocate(xyz(3,n,nall),er(nall),at(n),                            &
      &         coord(3,n),cnref(n),cn(n),tmp3(3,n),                    &
      &         clist(nall*nall),                                       &
      &         tmp(3,n,nall*nall)) 

      call rdcoord('coord',n,at,xyz(1:3,1:n,1))

! read file
      write(*,*) 'reading xyz file ...'
      open(unit=1,file=fname)
      i=1
 20   read(1,*,end=200) n
      read(1,'(a)') btmp
      call readl(btmp,xx,k)
      er(i)=xx(1)            ! Etot in au  
      cnorm=0
      do l=1,n
         read(1,*)a3,xyz(1:3,l,i)
         xyz(1:3,l,i)=xyz(1:3,l,i)/bohr
         if(i.eq.1) at(l) = e2i(a3)
         cnorm=cnorm+sum(abs(xyz(1:3,l,i)))
      enddo
      if(cnorm.gt.1.d-6) i=i+1
      goto 20
200   close(1)

      tmp3(1:3,1:n)=xyz(1:3,1:n,1)
      call ycoord(n,rcov,at,tmp3,cnref,100.0d0)   ! best = ref CN !SG loose threshold

      allocate(na(n),nb(n),nc(n),c0(3,n),         &
      &        c1(3,n),c2(3,n),new(3,n),at2(n))

      nmax=nall
      do i=1,nall
         if((er(i)-er(1))*autokcal.ge.ewin) then  ! E window
           nmax=i-1
           exit
         endif
      enddo

      call XYZINT(tmp3,n,na,nb,nc,1.0d0,c0)
!     call zmatpr(n,at,c0,na,nb,nc,1)

      nmaxmax=lin(nmax-1,nmax)
      write(*,*) '# in E window       ',nmax
      write(*,*) 'generating pairs ...',nmaxmax
      allocate(cstat(nmaxmax),csave(3,n,nmaxmax),rms(nmaxmax))
      allocate(cint(3,n,nmax),ind(nmaxmax))

      if(env%niceprint)then
        call printemptybar()
      endif

! all internals 
      do i=1,nmax
         call XYZGEO(XYZ(1,1,i),N,NA,NB,NC,1.0d0,c1)
         cint(1:3,1:n,i)=c1(1:3,1:n)              
      enddo
! first built them all
      cstat=0
      pcount=0
!$OMP PARALLEL PRIVATE(m,k,new,coord,l,fail,ldum) &
!$OMP SHARED(cnref,rcov,cthr,c0,cstat,csave,n,na,nb,nc,xyz,cint,rms) &
!$omp shared(percent,pcount,bar)
!$OMP DO 
      do i=1,nmax
         do j=1,i-1
            l=lin(j,i) ! common array pointer for all threads
            do m=1,n 
               do k=1,3
                  new(k,m)=c0(k,m)+cint(k,m,j)-cint(k,m,i) ! new internal coords
               enddo
               if(pi-new(2,m).lt.0.001) new(2,m)=pi-0.005 ! prevent linear bends
               if(   new(2,m).lt.0.001) new(2,m)=0.005
            enddo
            call GMETRY(n, new, coord, na,nb,nc,fail)
            if(fail) cstat(l)=99                    ! problem with Cart gen in GMETRY
            call ycoord2(n,rcov,at,coord,cnref,100.d0,cthr,ldum)  !SG check clash and skip directly
            if(ldum) cstat(l)=1 
            if(cstat(l).eq.0) then
               csave(1:3,1:n,l)=coord(1:3,1:n) ! take it
               call rmsd(n,tmp3,coord,0,Udum,xdum,ydum,dum,.false.,gdum) ! w.r.t. first
               rms(l)=dum
               cstat(l)=-1
            endif
         enddo

         call initsignal()
         !$omp critical
          pcount=pcount+1.0d0
          percent=pcount/float(nmax)*100d0
          if(env%niceprint)then
            call  progbar(percent,bar)
            call printprogbar(percent,bar)
          else
             if(mod(i,100).eq.0) write(*,'(f6.1,'' % done'')')percent
          endif
         !$omp end critical

      enddo
!$OMP END DO
!$OMP END PARALLEL      


      l=lin(nmax-1,nmax)
      ierr=0
      n2=0
      rmsdav=0
      do i=1,l
         if(cstat(i).ge.0)then
            if(cstat(i).eq.1) ierr=ierr+1 ! clash
         else
            n2=n2+1                       ! struc ok
            rmsdav=rmsdav+rms(i)
         endif
      enddo
      rmsdav=rmsdav/dble(n2)
      sd=0
      do i=1,l
         ind(i)=i 
         if(cstat(i).eq.-1)then
            rms(i)=-(rms(i)-rmsdav)**2    ! sort criterion
            sd=sd-rms(i)
         else
            rms(i)=1.d+9
         endif
      enddo
      call qsort(rms,1,l,ind)             ! sort

      if(env%niceprint) write(*,*)

      if(n2.lt.1)then
        write(*,*)'number of clash discarded :',ierr
        write(*,*) 'No new structures from GC.'
         kk=0
         return
         goto 666
      else
        write(*,*)'generated pairs           :',n2  
        write(*,*)'number of clash discarded :',ierr
      endif     

      write(*,'(1x,a,f8.5)')'average rmsd w.r.t input  :',rmsdav   
      write(*,'(1x,a,f8.5)')'sd of ensemble            :',sqrt(sd/dble(n2)) 
      kk=min(2*maxgen,n2) ! check faktor 2 more than needed

      l=0
      ident=0
      outer: do ii=1,kk                  
            i=ind(ii)  ! sort index
            coord(1:3,1:n)=csave(1:3,1:n,i)
!$OMP PARALLEL PRIVATE(c2,Udum,xdum,ydum,gdum,dum) &
!$OMP SHARED(tmp,rms,coord)
!$OMP DO     
            do kl=1,l  
               c2(1:3,1:n)=tmp(1:3,1:n,kl)
               call rmsd(n,c2,coord,0,Udum,xdum,ydum,dum,.false.,gdum) 
               rms(kl)=dum
            enddo
!$OMP END DO
!$OMP END PARALLEL      
            do kl=1,l    
               if(rms(kl).lt.rthrbohr) then
                  ident=ident+1
                  cycle outer                          ! remove identical
               endif
            enddo

            l=l+1                                      ! save and write new one
            tmp(1:3,1:n,l)=coord(1:3,1:n)              
!           next pair
      enddo outer

      deallocate(cstat,csave,rms)

      write(*,*)'number of new structures      :',l     
      if(l.lt.1)then
         write(*,*) 'No new structures from GC.'
         kk=0
         return
         goto 666
      end if
      write(*,*)'removed identical structures  :',ident

      kk=min(maxgen,l)

!---- don't pollute the working directory, so make an dir OPTIM, for the optimizations
      call getcwd(thispath)
      call rmrf('OPTIM')
      io = makedir('OPTIM')
      call copysub('coord','OPTIM')
      !call copysub('.CHRG','OPTIM')
      !call copysub('.UHF','OPTIM')
      call env%wrtCHRG('OPTIM')   
      call copysub(env%fixfile,'OPTIM')
      !call copysub(env%constraints,'OPTIM')
      if(env%useqmdff)then
      call copysub('solvent','OPTIM')
      endif
      if(env%gfnver=='--gff')then
         io = sylnk(trim(thispath)//'/'//'gfnff_topo','OPTIM'//'/'//'gfnff_topo')
      endif


      call chdir('OPTIM')
!---------- 
!---- write TMPCONF dirs   
      write(*,*)'writing ',kk,' TMPCONF* dirs ...'     
      do i=1,kk
         j=i        
         write(ctmp,'(''TMPCONF'',i0)')i
         io = makedir(trim(ctmp))
         btmp=trim(ctmp)//'/'//'confcross.xyz'
         open(newunit=ich,file=btmp)
         tmp(1:3,1:n,j)=tmp(1:3,1:n,j)*bohr
         call wrxyz(ich,n,at,tmp(:,:,j))
         write(ich,'(a)')'$end'
         write(ich,'(a)')'$opt'
         if(env%crestver.eq.2)then
         write(ich,'(1x,a,f14.4)')'hlow=',env%hlowopt
         write(ich,'(1x,a,i0)')'microcycle=',nint(env%microopt)
         write(ich,'(1x,a,f14.4)')'s6=',env%s6opt
         endif
         if(env%gcmultiopt)then
           write(ich,'(1x,a,i6)')'optlevel=',-1
         else
           write(ich,'(1x,a,i6)')'optlevel=',nint(env%optlev)
         endif
         write(ich,'(a)')'$end'
         call write_cts(ich,env%cts)
         call write_cts_biasext(ich,env%cts)
         close(ich)

         !call copysub('.CHRG',trim(ctmp))
         !call copysub('.UHF',trim(ctmp))
         call env%wrtCHRG(trim(ctmp))   
         call copysub(env%fixfile,trim(ctmp))
         if(env%useqmdff)then
         call copysub('solvent',trim(ctmp))
         endif
         if(env%gfnver=='--gff')then
           io = sylnk(trim(thispath)//'/'//'OPTIM'//'/'//'gfnff_topo',trim(ctmp)//'/'//'gfnff_topo')
         endif

      enddo
!---- change back to working directory
      call chdir(trim(thispath))

666   continue

      deallocate(ind)
      deallocate(at2,new,c2,c1,c0,nc,nb,na)
      deallocate(tmp3,cn,cnref,coord,at,er,xyz)
      deallocate(tmp,clist)
      deallocate(ydum,xdum,Udum,gdum)
      deallocate(rcov,xx,cint)


      end associate

end subroutine confcross

