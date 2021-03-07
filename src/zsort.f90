!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2018-2020 Stefan Grimme, Jakob Seibert
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
subroutine zsort  
      use iomod
      use strucrd, only: rdnat,rdcoord,wrc0
      implicit none
      integer n
      real*8,allocatable::xyz(:,:)                             ! comment out if subroutine
      real*8,allocatable::xyznew(:,:)
      real*8,allocatable::cn(:),tmp(:,:),tmpgeo(:,:)
      real*8,allocatable::geo(:,:)                             ! comment out if subroutine
      integer,allocatable::nat_mols(:),mvec(:),itmp(:),atnew(:)        
      integer,allocatable::na(:),nb(:),nc(:),at(:)             ! comment out if subroutine
      integer,allocatable::tmpna(:),tmpnb(:),tmpnc(:)          
      real*8,allocatable :: rcov(:),xx(:)     

      integer :: i,j,k,l,nn,nmol,refmol
      integer :: i1,i2,fail,try,zmatcnt,zmatstart

      character(len=80) ::fname

      logical ex

!      write(*,*)
!      write(*,'(7x,''========================================'')')
!      write(*,'(7x,''|             Z S O R T                |'')')
!      write(*,'(7x,''|JS, Universitaet Bonn, MCTC, 05/2017  |'')')
!      write(*,'(7x,''========================================'')')
!      write(*,*)'input on file coord, Z matrix on file zmatrix'
!      write(*,*)'sorted coord on file zcoord'

      allocate(xx(10),rcov(94))

      call setrcov(rcov) ! for CN calc
      inquire(file='coord', exist=ex)                        ! comment out if subroutine
      if(.not.ex) stop 'ERROR: no coord file in folder!'     ! comment out if subroutine

      call rdnat('coord',n)
      if(n < 4)then !without many atoms there is no point in sorting the zmat
          call copy('coord','zcoord')
          return
      endif

      allocate(xyz(3,n),at(n),mvec(n),xyznew(3,n),atnew(n), &! comment out if subroutine
      &        geo(3,n),na(n),nb(n),nc(n))

                
      write(*,*)'total number of atoms :',n            ! comment out if subroutine

      call rdcoord('coord',n,at,xyz)                ! comment out if subroutine

      call mrec(nmol,xyz,n,at,mvec)
      write(*,*)'total number of frags :',nmol         ! comment out if subroutine
      allocate(nat_mols(nmol))
      nat_mols=0
      zmatcnt=1
      zmatstart=0
! check for largest fragment
       do i=1,nmol
         nn=0
         do j=1,n
            if(mvec(j).eq.i) nn=nn+1
         enddo
         nat_mols(i)=nn
       enddo
       refmol=maxloc(nat_mols,1)
! loop over fragments
      do k=1,nmol
           
      i=maxloc(nat_mols,1)   
      write(*,*) "Fragment ",i                         ! comment out if subroutine
         allocate(itmp(nat_mols(i)),tmp(3,nat_mols(i)), &
     &     tmpna(nat_mols(i)) &
     &   ,tmpnb(nat_mols(i)),tmpnc(nat_mols(i)),tmpgeo(3,nat_mols(i)))
         nn=0
         do j=1,n
            if(mvec(j).eq.i)then
               nn=nn+1
               tmp(1:3,nn)=xyz(1:3,j)
               itmp(nn)=at(j)
            endif
         enddo
         call recursive_zmat(tmp,itmp,nat_mols(i),tmpna,tmpnb,tmpnc,tmpgeo)
         call zmatpr(nat_mols(i),itmp,tmpgeo,tmpna,tmpnb,tmpnc,i)      ! comment out if subroutine  ! write zmat file
                 
        if (i .eq. refmol) then
           do j=1,nat_mols(i)
               na(j)=tmpna(j)
               nb(j)=tmpnb(j)
               nc(j)=tmpnc(j)
               atnew(j)=itmp(j)
               xyznew(:,j)=tmp(:,j)
               zmatcnt=zmatcnt+1
           enddo
        else
            na(1+zmatstart)=1
            nb(1+zmatstart)=2
            nc(1+zmatstart)=3
            atnew(1+zmatstart)=itmp(1)
            xyznew(:,1+zmatstart)=tmp(:,1)
            zmatcnt=zmatcnt+1
            if(nat_mols(i).gt.1)then
            na(2+zmatstart)=tmpna(2)+zmatstart
            nb(2+zmatstart)=1
            nc(2+zmatstart)=2
            atnew(2+zmatstart)=itmp(2)
            xyznew(:,2+zmatstart)=tmp(:,2)
            zmatcnt=zmatcnt+1
            endif
            if(nat_mols(i).gt.2)then
            na(3+zmatstart)=tmpna(3)+zmatstart
            nb(3+zmatstart)=tmpnb(3)+zmatstart
            nc(3+zmatstart)=1
            atnew(3+zmatstart)=itmp(3)
            xyznew(:,3+zmatstart)=tmp(:,3)
            zmatcnt=zmatcnt+1
            endif
            do j=4,nat_mols(i)
               na(j+zmatstart)=tmpna(j)+zmatstart
               nb(j+zmatstart)=tmpnb(j)+zmatstart
               nc(j+zmatstart)=tmpnc(j)+zmatstart
               atnew(j+zmatstart)=itmp(j)
               xyznew(:,j+zmatstart)=tmp(:,j)
               zmatcnt=zmatcnt+1
            enddo

        endif
         deallocate(itmp,tmp,tmpna,tmpnb,tmpnc,tmpgeo)
         zmatstart=zmatstart+nat_mols(i)
         nat_mols(i)=0 
      enddo 
      if (nmol .ne. 1) then
      !write(*,*) "Merged Z matrix"                           ! comment out if subroutine
      call xyzgeo(xyznew,n,na,nb,nc,1.0d0,geo)
      !call zmatpr(n,atnew,geo,na,nb,nc,0)                    ! comment out if subroutine
      endif

!     if you want to overwrite old xyz and at
!      xyz=xyznew
!      at=atnew
 
      call wrc0('zcoord',n,atnew,xyznew)

      write(*,*)'terminated normally'

      deallocate(xyz,at,mvec,xyznew,atnew,geo,na,nb,nc,nat_mols)
!     deallocate(mvec,xyznew,atnew,nat_mols) if subroutine
      deallocate(rcov,xx)

end subroutine zsort

subroutine recursive_zmat(xyz,at,nat,na,nb,nc,geo)
      implicit none
!     recursive Z matrix generation
!     i/o block
      integer, intent(in) :: nat        ! number of atoms
      integer, intent(in) :: at(nat)    ! element number of atoms
      real*8, intent(in) :: xyz(3,nat)  ! coordinates
      integer, intent(out) :: na(nat)   ! column 1 of Z matrix (distance)
      integer, intent(out) :: nb(nat)   ! column 2 of Z matrix (angle)
      integer, intent(out) :: nc(nat)   ! column 3 of Z matrix (dihedral)
      real*8, intent(out)  :: geo(3,nat)! dist, ang and dihed in Z matrix
!     further variables
      real*8 cn(nat),bond(nat,nat)
      integer i,zmatcnt,nref(nat),zmatpos(nat)
      logical taken(nat)
!     initialization
      geo=0.0d0
      nref=0
      taken=.false.
      zmatpos=0
      zmatcnt=1
      na=0
      nb=0
      nc=0
!     code
      call xcoord(nat,at,xyz,cn,bond)
      do i=1,nat
!       j=maxloc(cn,1)
       if(.not.taken(i)) then
         
         taken(i)=.true.
         zmatpos(i)=zmatcnt
         call find_buddies(i,xyz,taken,nat,cn,bond,na,nb,nc, &
         &                 zmatcnt,nref,zmatpos)
      endif
      enddo
      call zmatsort(nat,at,xyz,na,nb,nc,zmatpos) ! sort Z matrix 
      call xyzgeo(xyz,nat,na,nb,nc,1.0d0,geo)    ! compute dist, ang, dihed
end subroutine recursive_zmat

recursive subroutine find_buddies(i,xyz,taken,nat,cn,bond,na,nb,nc,zmatcnt,nref,zmatpos)
      implicit none
      real*8 xyz(3,nat),cntmp,cn(nat),bond(nat,nat)
      integer i,nat,j,icn,k,na(nat),nb(nat),nc(nat)
      integer zmatcnt, nref(nat),zmatpos(nat)
      logical taken(nat)

      icn=nint(cn(i))

      do k=1,icn
         j=maxloc(bond(:,i),1)
         bond(j,i)=0.0d0
         
         if (i .eq. j) cycle
         if (.not.taken(j)) then
            taken(j)=.true.
            zmatcnt=zmatcnt+1
            if (zmatcnt .gt. nat ) exit
            zmatpos(j)=zmatcnt
            if (nref(i) .ne. 0) then
               if (nref(nref(i)) .eq. 0 ) then
                  na(j)=na(nref(i))
                  if (zmatcnt .le. 3 .or. nb(nref(i)).eq. 0) then!
                     nb(j)=nref(i)!
                     nref(i)=j!
                  else!
                     nb(j)=nb(nref(i))
                  endif! 
                  if (zmatcnt .le. 3 ) then !
                     nc(j)=0!
                  else!
                     nc(j)=nref(i)
                  endif
               else if (nref(nref(i)) .ne. 0 ) then
                  na(j)=i
                  nb(j)=nref(i)
                  if (nb(nref(i)) .eq. 0) then
                     nc(j)=nref(nref(i))
                  else
                     nc(j)=nb(nref(i))
                  endif
                   nref(i)=j
               else
                  na(j)=i                   
                  nb(j)=na(i)                   
                  nc(j)=nb(i)              
                  nref(i)=j
               endif
           else
               na(j)=i
               if (na(i) .eq. 0 .and. nref(i) .ne. 0) then !
               nb(j)=nref(i) !
               else!
               nb(j)=na(i)
               endif!
               nc(j)=nb(i)
               nref(i)=j
           endif    
         call find_buddies(j,xyz,taken,nat,cn,bond,na,nb,nc,zmatcnt,nref,zmatpos)
         endif
      enddo
end subroutine find_buddies

subroutine zmatsort(nat,at,xyz,na,nb,nc,zmatpos)
      implicit none
!     sort zmat according to zmat position from previous subroutine     
      integer nat, at(nat), na(nat),nb(nat),nc(nat),zmatpos(nat)
      integer tmpat(nat),tmpna(nat),tmpnb(nat),tmpnc(nat)
      integer i,j,k,l,oldpos,mode
      real*8 xyz(3,nat),tmpxyz(3,nat)
      logical considered(nat)
        tmpat=0
        tmpna=0
        tmpnb=0
        tmpnc=0
        tmpxyz=0.0d0
         do i=nat,1,-1
      
            oldpos=maxloc(zmatpos,1)  
            tmpat(i)=at(oldpos)
            if (na(oldpos) .ne. 0 ) then
               tmpna(i)=zmatpos(na(oldpos))
            else 
               tmpna(i)=0
            endif
            if (nb(oldpos) .ne. 0 ) then
               tmpnb(i)=zmatpos(nb(oldpos))
            else
               tmpnb(i)=0
            endif
            if (nc(oldpos) .ne. 0 ) then
               tmpnc(i)=zmatpos(nc(oldpos))
            else
               tmpnc(i)=0
            endif

            tmpxyz(:,i)=xyz(:,oldpos)
            
            zmatpos(oldpos)=0

         enddo

          at=tmpat
          na=tmpna
          nb=tmpnb
          nc=tmpnc
          xyz=tmpxyz
!     check if Z matrix reasonable
      i=1
      if ( i .eq. 1) then 
            if ( na(1) .ne. 0  .or. nb(1) .ne. 0  &
     &      .or. nc(1) .ne. 0 ) stop "ERROR: first Z matrix  &
     &     entry is erroneous"
            if (nat .gt. 1) then
               if ( na(2) .eq. 0 .or. nb(2) .ne. 0  &
     &         .or. nc(2) .ne. 0 ) stop "ERROR: second Z matrix &
     &         entry is erroneous"
               if (nat .gt. 2) then
                  if ( na(3) .eq. 0 .or. nb(3) .eq. 0  &
     &            .or. nc(3) .ne. 0 ) stop "ERROR: third Z matrix &
     &            entry is erroneous"
                if (nat .gt. 3) then
                  do i=4,nat
                     if(na(i).eq.0.or.nb(i).eq.0.or.nc(i).eq.0) then 
                      write(*,*) "ERROR: Erroneous entry in Z matrix at position ",i
                      stop
                     endif

                  enddo
               endif
              endif
           endif
      endif
end subroutine zmatsort

!------------------------------------------------------------------------
subroutine zmatpr(nat,at,geo,na,nb,nc,molnum)
      use strucrd, only: i2e
      implicit none
      integer :: nat,na(nat),nb(nat),nc(nat),at(nat)
      real*8 :: geo(3,nat)
      character(len=20) ::  filename
      logical ex
      integer i,l,m,n,molnum
      real*8 bl,ang,dihed,pi,au2ang
      parameter (pi =  3.14159265358979D0)
      parameter (au2ang = 0.529177210920d0 )
         do  i=1,nat             
            l=1                  
            m=1                  
            n=1                  
            if(i.eq.1)then       
               l=0               
               m=0               
               n=0               
            endif                
            if(i.eq.2)then       
               m=0               
               n=0               
            endif                
            if(i.eq.3)n=0        
         bl=geo(1,i)*au2ang
         ang=geo(2,i)*180./pi        
         dihed=geo(3,i)*180./pi          
         if(dihed.gt.180.0d0)dihed=dihed-360.0d0                       
!        write(*,'(2x,a2,f14.8,i4,f14.8,i5,f14.8,i5,i6,2i5)') &
!        & i2e(at(i)),bl,l,ang,m,dihed,n,na(i),nb(i),nc(i)                                                      
         write(*,'(i4,2x,a2,f12.6,2x,f10.4,2x,f10.4,i6,2i5)') &
         &   i,i2e(at(i)),bl,ang,dihed,na(i),nb(i),nc(i)   
         enddo
         
         write(filename,'("zmatrix",i0,".zmat")'),molnum
         open(unit=42, file=filename)

         write(42,'(a2)') i2e(at(1))
         if(nat.gt.1)then
         write(42,'(a2,x,i0,x,f8.3)') i2e(at(2)), na(2), geo(1,2)*au2ang
         endif
         if(nat.gt.2)then
         write(42,'(a2,x,i0,x,f8.3,x,i0,x,f8.3)') i2e(at(3)), na(3) &
         &                  ,geo(1,3)*au2ang,nb(3), geo(2,3)*180./pi     
         endif
 
         do i=4,nat
            bl=geo(1,i)*au2ang
            ang=geo(2,i)*180./pi         
            dihed=geo(3,i)*180./pi         
            if(dihed.gt.180.0d0)dihed=dihed-360.0d0
            write(42,'(a2,x,i0,x,f8.3,x,i0,x,f8.3,x,i0,x,f8.3)') &
            &         i2e(at(i)),na(i),bl,nb(i),ang,nc(i),dihed 
         enddo
         write(42,*) 
         close(42)                                             
end subroutine zmatpr
