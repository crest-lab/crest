!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2018-2020 Philipp Pracht
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

!=================================================================!
! subroutine mrec
! molcount: number of total fragments (increased during search)
! xyz: overall Cart. coordinates
! nat: overall number of atoms
! at: atomic number array
! molvec: assignment vector of atom to fragment
!=================================================================!
subroutine mrec(molcount,xyz,nat,at,molvec)
      implicit none
      real*8 xyz(3,nat),cn(nat),bond(nat,nat)
      integer nat,molvec(nat),i,molcount,at(nat)
      logical taken(nat)
      molvec=0
      molcount=1
      taken=.false.
      cn=0.0d0
      bond=0.0d0
      call xcoord(nat,at,xyz,cn,bond)
      do i=1,nat
       if(.not.taken(i)) then
         molvec(i)=molcount
         taken(i)=.true.
         call neighbours(i,xyz,at,taken,nat,cn,bond,molvec,molcount)
         molcount=molcount+1
      endif
      enddo
      molcount=molcount-1
end subroutine mrec

!=================================================================!
! subroutine mreclm
! a variant of the mrec routine with less allocate statements
! should be faster than the old version if several thousand
! calls are done.
!=================================================================!
subroutine mreclm(molcount,nat,at,xyz,molvec,bond,rcov,cn)
      use iso_fortran_env, only: wp => real64
      implicit none
      integer :: molcount,nat
      integer :: at(nat),molvec(nat)
      real(wp) :: xyz(3,nat)
      real(wp) :: cn(nat),bond(nat,nat)
      real(wp) :: bref(nat,nat)
      real(wp) :: rcov(94)
      integer :: i
      logical :: taken(nat)
      molvec=0
      molcount=1
      taken=.false.
      cn=0.0d0
      bond=0.0d0
      !call xcoord(nat,at,xyz,cn,bond)
      call xcoord2(nat,at,xyz,rcov,cn,400.0_wp,bond)
      bref = bond
      do i=1,nat
       if(.not.taken(i)) then
         molvec(i)=molcount
         taken(i)=.true.
         call neighbours(i,xyz,at,taken,nat,cn,bond,molvec,molcount)
         molcount=molcount+1
      endif
      enddo
      molcount=molcount-1
      bond = bref
end subroutine mreclm

!==================================================================================!
recursive subroutine neighbours(i,xyz,iat,taken,nat,cn,bond,molvec,molcnt)
      implicit none
      real*8 xyz(3,nat),tr,xi(3),cn(nat),bond(nat,nat)
      integer i,nat, molcnt,molvec(nat),j,iat(nat),icn,k
      logical taken(nat)
      tr=2.0d0
      
      xi(1:3)=xyz(1:3,i) 
      icn=nint(cn(i))
      do k=1,icn
         j=maxloc(bond(:,i),1)
         bond(j,i)=0.0d0
         if (i .eq. j) cycle
         if (.not.taken(j)) then
            molvec(j)=molcnt
            taken(j)=.true.
            call neighbours(j,xyz,iat,taken,nat,cn,bond,molvec,molcnt)
         endif
      enddo
end subroutine neighbours

!================================================================================!
! compute coordination numbers by adding an inverse damping function
!================================================================================!
subroutine xcoord(natoms,iz,xyz,cn,bond)
      implicit none
      integer iz(natoms),natoms,i,k1
      real*8 xyz(3,natoms),cn(natoms)
      real*8 cn_thr,bond(natoms,natoms)
      integer iat
      real*8 dx,dy,dz,r,damp,xn,rr,rco,r2,rcovi,rcovj
      real*8,allocatable :: rcov(:)
      allocate(rcov(94))
      call setrcov(rcov)

      cn_thr=400.0d0
      k1=16
      bond=0.0d0
      cn=0.0d0
      do i=1,natoms
        xn=0.0d0
        rcovi=rcov(iz(i))
        do iat=1,natoms
         if(iat.ne.i)then
            dx=xyz(1,iat)-xyz(1,i)
            dy=xyz(2,iat)-xyz(2,i)
            dz=xyz(3,iat)-xyz(3,i)
            r2=dx*dx+dy*dy+dz*dz
            r=sqrt(r2)
            if (r2.gt.cn_thr) cycle
            rcovj=rcov(iz(iat))
! covalent distance in Bohr
            rco=(rcovi+rcovj)*1.0  ! this scaling reduces the size of the clusters
            rr=rco/r
! counting function exponential has a better long-range behavior than MHGs inverse damping
            damp=1.d0/(1.d0+exp(-k1*(rr-1.0d0)))
            bond(iat,i)=damp
            xn=xn+damp
         endif
        enddo
        cn(i)=xn
      enddo
      deallocate(rcov)
end subroutine xcoord


subroutine ycoord2(natoms,rcov,iz,xyz,cn,cn_thr,cthr,clash)
      implicit none  
      real*8 k1,k3
      parameter (k1     =16)
      parameter (k3     =-4)
      integer iz(*),natoms,i
      real*8 xyz(3,*),cn(*),rcov(94)
      real*8 cn_thr,cthr
      logical clash

      integer iat    
      real*8 dx,dy,dz,r,damp,xn,rr,rco,r2

      clash=.false.
      do i=1,natoms
      xn=0.0d0
      do iat=1,natoms
         if(iat.ne.i)then
            dx=xyz(1,iat)-xyz(1,i)
            dy=xyz(2,iat)-xyz(2,i)
            dz=xyz(3,iat)-xyz(3,i)
            r2=dx*dx+dy*dy+dz*dz 
            if (r2.gt.cn_thr) cycle 
            r=sqrt(r2)                  
!c covalent distance in Bohr
            rco=rcov(iz(i))+rcov(iz(iat))
            rr=rco/r
!c counting function exponential has a better long-range behavior than MHGs inverse damping
            damp=1.d0/(1.d0+exp(-k1*(rr-1.0d0)))
            xn=xn+damp*iz(iat)
         endif
      enddo
      if(abs(cn(i)-xn).gt.cthr) then                  !   clash check
         clash=.true.
         return
      endif
      enddo
end subroutine ycoord2


subroutine ycoord(natoms,rcov,iz,xyz,cn,cn_thr)
      implicit none  
      real*8 k1,k3
      parameter (k1     =16)
      parameter (k3     =-4)
      integer iz(*),natoms,i
      real*8 xyz(3,*),cn(*),rcov(94)
      real*8 cn_thr

      integer iat    
      real*8 dx,dy,dz,r,damp,xn,rr,rco,r2

      do i=1,natoms
      xn=0.0d0
      do iat=1,natoms
         if(iat.ne.i)then
            dx=xyz(1,iat)-xyz(1,i)
            dy=xyz(2,iat)-xyz(2,i)
            dz=xyz(3,iat)-xyz(3,i)
            r2=dx*dx+dy*dy+dz*dz 
            if (r2.gt.cn_thr) cycle 
            r=sqrt(r2)                  
!c covalent distance in Bohr
            rco=rcov(iz(i))+rcov(iz(iat))
            rr=rco/r
!c counting function exponential has a better long-range behavior than MHGs inverse damping
            damp=1.d0/(1.d0+exp(-k1*(rr-1.0d0)))
            xn=xn+damp*iz(iat)
         endif
      enddo
      cn(i)=xn
      enddo

end subroutine ycoord

subroutine ycoord3(natoms,iz,xyz,cn,cn_thr,bond)
      implicit none
      integer :: natoms      ! = n
      integer :: iz(natoms)  ! = at
      integer :: i,k1
      real*8 xyz(3,natoms),cn(natoms)
      real*8 cn_thr,bond(natoms,natoms)
      integer iat
      real*8 dx,dy,dz,r,damp,xn,rr,rco,r2,rcovi,rcovj
      real*8,allocatable :: rcov(:)
     
      allocate(rcov(94))
      !cn_thr=1600.0d0
      call setrcov(rcov)
      k1=16
      bond=0.0d0
      cn=0.0d0
      do i=1,natoms
        xn=0.0d0
        !call setrcov(iz(i),rcovi)
        rcovi=rcov(iz(i))
        do iat=1,natoms
         if(iat.ne.i)then
            dx=xyz(1,iat)-xyz(1,i)
            dy=xyz(2,iat)-xyz(2,i)
            dz=xyz(3,iat)-xyz(3,i)
            r2=dx*dx+dy*dy+dz*dz
            r=sqrt(r2)
            if (r2.gt.cn_thr) cycle
            !call setrcov(iz(iat),rcovj)
            rcovj=rcov(iz(iat))
!c covalent distance in Bohr
            rco=rcovi+rcovj
            rr=rco/r
!            rr=rr*0.90d0
!c counting function exponential has a better long-range behavior than MHGs inverse damping
            damp=1.d0/(1.d0+exp(-k1*(rr-1.0d0)))
            bond(iat,i)=damp
!            bond(i,iat)=damp
            xn=xn+damp
         endif
        enddo
        cn(i)=xn
      enddo
      deallocate(rcov)
end subroutine ycoord3

subroutine ncoord(natoms,rcov,iz,xyz,cn,cn_thr)
      implicit none  
      real*8 :: k1,k3
      parameter (k1     =16)
      parameter (k3     =-4)
      integer :: iz(*),natoms,i
      real*8 :: xyz(3,*),cn(*),rcov(94)
      real*8 :: cn_thr

      integer :: iat    
      real*8 :: dx,dy,dz,r,damp,xn,rr,rco,r2

      do i=1,natoms
      xn=0.0d0
      do iat=1,natoms
         if(iat.ne.i)then
            dx=xyz(1,iat)-xyz(1,i)
            dy=xyz(2,iat)-xyz(2,i)
            dz=xyz(3,iat)-xyz(3,i)
            r2=dx*dx+dy*dy+dz*dz 
            if (r2.gt.cn_thr) cycle 
            r=sqrt(r2)                  
!c covalent distance in Bohr
            rco=rcov(iz(i))+rcov(iz(iat))
            rr=rco/r
!c counting function exponential has a better long-range behavior than MHGs inverse damping
            damp=1.d0/(1.d0+exp(-k1*(rr-1.0d0)))
            xn=xn+damp
         endif
      enddo
      cn(i)=xn
      enddo
end subroutine ncoord


