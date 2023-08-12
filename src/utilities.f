!================================================================================!
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
!================================================================================!

!TODO replace with utilmod.f90


!--- formerly in "lin.f"
***********************************************************************
* address in packed array
***********************************************************************

      integer function lin(i1,i2)
      integer i1,i2,idum1,idum2
      idum1=max(i1,i2)
      idum2=min(i1,i2)
      lin=idum2+idum1*(idum1-1)/2        
      return
      end

      !--- determine i and j for a given index k
      !    (i.e., the reverse of the lin function)
      subroutine revlin(k,i,j)
          implicit none
          integer*8 :: k,b,x
          integer :: i,j
          real*8 :: kf
          real*8 :: idum
          kf = real(k, 8)
          idum = 1.0 + 8.0*kf
          !idum = sqrt(idum)
          idum = idum ** 0.5
          idum = (1.0 + idum)/2.0
          i =  floor(idum)
          x = i-1
          x = x*i      
          b = x/2
          b = k - b    
          j = int(b)
          return
      end subroutine revlin

      !--- determine i and j for a given index k
      !    by counting until k is reached (only for testing, because too
      !    expensive)
      subroutine revlin_count(k,i,j,dimij)
          implicit none
          integer*8 :: k
          integer :: i,j
          integer :: a,b
          integer :: dimij
          integer*8 :: kdum
          integer*8 :: lina

          OUTER : do a=1,dimij
              do b = 1,a
                 kdum = lina(a,b)
                 if(kdum == k)then
                    i=a
                    j=b
                    exit OUTER
                 endif
              enddo
          enddo OUTER

          return                                        
      end subroutine revlin_count

      integer*8 function lina(i1,i2)
      integer i1,i2,idum1,idum2
      idum1=max(i1,i2)
      idum2=min(i1,i2)
      !lina=idum2+idum1*(idum1-1)/2        
      lina=idum1
      lina=lina*(idum1-1)
      lina=lina/2
      lina=lina+idum2
      return
      end

      integer*8 function linr(o1,o2,i)
       integer*8 :: o1
       integer   :: o2,i
       linr = (o1 + 1)
       linr = linr + i
       linr = linr - o2
       return
      end

!--- formerly in "boltz.f"
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine boltz(n,t,e,p)
      real*8 :: e(*),p(*)          
      real*8,allocatable :: e2(:)
      real*8 :: t,f,hsum,esum
      allocate(e2(n))
      !f=8.314*t/4.184d+3
      f = 0.593d0 / 298.15d0
      f = f*t
      esum=0
      do i=1,n
         !e2(i) = (e(i) -e(1))* 627.5095
          e2(i) = e(i)
      enddo

      do i=1,n
         esum=esum+exp(-e2(i)/f)
      enddo
      hsum=0
      do i=1,n
         p(i)=exp(-e2(i)/f)/esum
      enddo
      deallocate(e2)
      end

      subroutine boltz2(n,e,p)
      use iso_fortran_env, only : wp => real64
      implicit none

      integer,intent(in)   :: n     ! Number of molecules
      real(wp),intent(in)  :: e(n)  ! Molecule energies
      real(wp),intent(out) :: p(n)  ! Population

      !> Boltzman constant k in Eh/K
      real(wp),parameter :: kh = 3.1668114d-6 
      !> Room temperature 
      real(wp),parameter :: T = 298.15_wp     

      real(wp) :: val
      real(wp) :: denom
      real(wp) :: emin
      integer  :: i

      p = 0.0_wp
      emin = minval(e)

      do i=1,n
         val  = -(e(i)-emin) / (kh * T)
         p(i) = exp(val)
      enddo
      denom = sum(p)
      p = p / denom
      end subroutine boltz2


!--- formerly in "rmsd.f"
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c but including H in OH      

      subroutine heavyrmsd(n,nall,k,l,xyz,at,rmsdval)
      use ls_rmsd
      use crest_parameters, only: bohr
      implicit none
      integer n,at(n),j,nall,k,l,nn
      real*8 xyz(3,n,nall),rmsdval

      logical oh,ohbonded2

      real*8 ,allocatable :: xyz1(:,:), xyz2(:,:)
      !Dummys:
      real*8  g(3,3), U(3,3), x_center(3), y_center(3)
      integer i
      
      nn=0
      do j=1,n
      oh=ohbonded2(n,j,xyz(1,1,k),at)
      if(at(j).gt.2.or.oh) nn=nn+1
      enddo
      allocate(xyz1(3,nn),xyz2(3,nn))

      i=0
      do j=1,n
      oh=ohbonded2(n,j,xyz(1,1,k),at)
      if(at(j).gt.2.or.oh) then
         i=i+1
         xyz1(1:3,i)=xyz(1:3,j,k)*bohr
         xyz2(1:3,i)=xyz(1:3,j,l)*bohr
      endif
      enddo

      call rmsd(i,xyz1,xyz2,0,U,x_center,y_center,rmsdval,.false.,g)

      deallocate(xyz1,xyz2)
      end

      subroutine heavyrmsdfile(k,l,rmsdval)
      use ls_rmsd
      use strucrd, only: rdnat,rdcoord
      use crest_parameters, only: bohr
      implicit none
      integer k,l
      real*8 rmsdval
      character(len=80) atmp1,atmp2
      !Dummys:
      real*8  g(3,3), U(3,3), x_center(3), y_center(3)
      integer n
      real*8 ,allocatable :: xyz1(:,:), xyz2(:,:)
      integer,allocatable :: at(:)

      call getname1(k,atmp1)
      call getname1(l,atmp2)
      call rdnat(atmp1,n)
      allocate(xyz1(3,n),xyz2(3,n),at(n))

      call rdcoord(atmp1,n,at,xyz1)
      call rdcoord(atmp2,n,at,xyz2)
      xyz1 = xyz1 * bohr
      xyz2 = xyz2 * bohr

      call rmsd(n,xyz1,xyz2,0,U,x_center,y_center,rmsdval,.false.,g)

      deallocate(xyz1,xyz2)
      end

      logical function ohbonded2(n,m,xyz,at)
      integer n,at(n),m             
      real*8 xyz(3,n)
      real*8 :: r

      ohbonded2=.false.
      if(at(m).ne.1) return

      do i=1,n
         if(i.eq.m.or.at(i).ne.8) cycle
         r  =sqrt((xyz(1,i)-xyz(1,m))**2
     .           +(xyz(2,i)-xyz(2,m))**2
     .           +(xyz(3,i)-xyz(3,m))**2)
         if(r*0.52917726d0.lt.1.1) then
          ohbonded2=.true.
          exit
         endif
      enddo

      end


!--- formerly in "ohbonded.f"
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      logical function ohbonded(n,m,xyz,at,acid)
      integer n,at(n),m,acid(86)       
      real*8 xyz(3,n)
      real*8 :: r

      ohbonded=.false.
      if(at(m).ne.1) return

      do i=1,n
         if(i.eq.m) cycle
         if(acid(at(i)).eq.0) cycle 
         r  =sqrt((xyz(1,i)-xyz(1,m))**2
     .           +(xyz(2,i)-xyz(2,m))**2
     .           +(xyz(3,i)-xyz(3,m))**2)
         if(r*0.52917726d0.lt.1.2) then
            ohbonded=.true.
            goto 99
         endif
      enddo
 99   continue
c     write(*,*) r*0.52917726,at(m),at(i),ohbonded

      end

!--- formerly in "dist.f"
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine distance(n,xyz,r)
      implicit none  
      integer n
      real*8 xyz(3,n),r(n,n)

      real*8 dx,dy,dz
      integer i,j

      do i=1,n
         do j=1,n       
            dx=xyz(1,j)-xyz(1,i)
            dy=xyz(2,j)-xyz(2,i)
            dz=xyz(3,j)-xyz(3,i)
            r(j,i)=sqrt(dx*dx+dy*dy+dz*dz)   
!c           r(i,j)=r(j,i)                    
         enddo
         r(i,i)=0
      enddo

      end 

      real*8 function distcma(n,j,xyz)
      implicit none  
      integer n,j
      real*8 xyz(3,n)

      real*8 dx,dy,dz

      dx=xyz(1,j)
      dy=xyz(2,j)
      dz=xyz(3,j)
      distcma=sqrt(dx*dx+dy*dy+dz*dz)   

      end 
