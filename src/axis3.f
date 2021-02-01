!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2017 Stefan Grimme
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

c*****************************************************************         
c transform to CMA for output on trj xyz file
c input coords remain unchanged
c*****************************************************************         
                                               
      subroutine axis3(mode,numat,nat,coord,coordout,eig)            

      use atmasses
      implicit none                                                
      integer numat,nat(numat),mode
      real*8  coord(3,numat),coordout(3,numat),eig(3)

      real*8 sumw,sumwx,sumwy,sumwz,atmass,xsum,eps
      real*8 t(6), evec(3,3)         
      real*8 x(numat),y(numat),z(numat),coordtmp(3,numat)
      real*8 coord1(3,numat)
      integer i,j,k
      data t /6*0.d0/                                              
                                                                   
      sumw=1.d-20                                                  
      sumwx=0.d0                                                   
      sumwy=0.d0                                                   
      sumwz=0.d0                                                   
      do 10 i=1,numat                                           
         if(mode.eq.0)then
            atmass=ams(nat(i))                                     
         else
            atmass=1./ams(nat(i))                                     
         endif
         sumw=sumw+atmass                                       
         sumwx=sumwx+atmass*coord(1,i)                          
         sumwy=sumwy+atmass*coord(2,i)                          
         sumwz=sumwz+atmass*coord(3,i)                          
   10 continue                                                  
                                                                   
      eps=1.d-3
      sumwx=sumwx/sumw                                             
      sumwy=sumwy/sumw                                             
      sumwz=sumwz/sumw                                             
      do i=1,numat                                              
         x(i)=coord(1,i)-sumwx 
         y(i)=coord(2,i)-sumwy 
         z(i)=coord(3,i)-sumwz 
         coordtmp(1,i)=x(i)
         coordtmp(2,i)=y(i)
         coordtmp(3,i)=z(i)
      enddo
      
      do 40 i=1,6                                                  
   40 t(i)=dble(i)*1.0d-10                                         
                                                                   
      do 50 i=1,numat                                           
         atmass=ams(nat(i))                                     
         t(1)=t(1)+atmass*(y(i)**2+z(i)**2)+eps                 
         t(2)=t(2)-atmass*x(i)*y(i)
         t(3)=t(3)+atmass*(z(i)**2+x(i)**2)+eps                 
         t(4)=t(4)-atmass*z(i)*x(i)
         t(5)=t(5)-atmass*y(i)*z(i)
         t(6)=t(6)+atmass*(x(i)**2+y(i)**2)+eps                 
   50 continue                                                  
                                                                   
      call rsp(t,3,3,eig,evec)                                     

c   now to orient the molecule so the chirality is preserved       
      xsum=evec(1,1)*(evec(2,2)*evec(3,3)-evec(3,2)*evec(2,3)) +    
     1     evec(1,2)*(evec(2,3)*evec(3,1)-evec(2,1)*evec(3,3)) +    
     2     evec(1,3)*(evec(2,1)*evec(3,2)-evec(2,2)*evec(3,1))      
      if( xsum .lt. 0) then                                         
         do 80 j=1,3                                               
   80    evec(j,1)=-evec(j,1)                                      
      endif                                                        
c     call prmat(6,evec,3,3,'Rmat')
                                                                   
      do 120 i=1,numat                                             
         do 120 j=1,3                                              
            xsum=0.d0                                               
            do 110 k=1,3                                           
  110       xsum=xsum+coordtmp(k,i)*evec(k,j)                           
  120 coord1(j,i)=xsum                                              
                                                                   
      do i=1,numat                                                 
         coordout(1,i)=coord1(1,i)                                    
         coordout(2,i)=coord1(2,i)                                    
         coordout(3,i)=coord1(3,i)                                    
      enddo                                                        
                                                                   
      return               
      end                 



c*****************************************************************         
c transform to CMA for output on trj xyz file
c input coords remain unchanged
c*****************************************************************         
                                               
      subroutine axis3heavy(mode,numat,nat,nm,coord,coordout,eig)

      use atmasses
      implicit none                                                
      integer numat,nat(numat),mode
      real*8  coord(3,numat),coordout(3,numat),eig(3)

      real*8 sumw,sumwx,sumwy,sumwz,atmass,xsum,eps
      real*8 t(6), evec(3,3)         
      real*8 x(numat),y(numat),z(numat),coordtmp(3,numat)
      real*8 coord1(3,numat)
      real*8 nm
      integer i,j,k
      data t /6*0.d0/                                              
                                                                   
      sumw=1.d-20                                                  
      sumwx=0.d0                                                   
      sumwy=0.d0                                                   
      sumwz=0.d0                                                   

      do 10 i=1,numat                                           
         if(mode.eq.0)then
            atmass=ams(nat(i))+nm                                     
!            atmass=ams(nat(i))
         else
            atmass=1./(ams(nat(i))+nm)                                  
!            atmass=1./ams(nat(i))
         endif
         sumw=sumw+atmass                                       
         sumwx=sumwx+atmass*coord(1,i)                          
         sumwy=sumwy+atmass*coord(2,i)                          
         sumwz=sumwz+atmass*coord(3,i)                          
   10 continue                                                  
                                                                   
      eps=1.d-3
      sumwx=sumwx/sumw                                             
      sumwy=sumwy/sumw                                             
      sumwz=sumwz/sumw                                             
      do i=1,numat                                              
         x(i)=coord(1,i)-sumwx 
         y(i)=coord(2,i)-sumwy 
         z(i)=coord(3,i)-sumwz 
         coordtmp(1,i)=x(i)
         coordtmp(2,i)=y(i)
         coordtmp(3,i)=z(i)
      enddo
      
      do 40 i=1,6                                                  
   40 t(i)=dble(i)*1.0d-10                                         
                                                                   
      do 50 i=1,numat                                           
         atmass=ams(nat(i))+nm                                  
!         atmass=ams(nat(i))
         t(1)=t(1)+atmass*(y(i)**2+z(i)**2)+eps                 
         t(2)=t(2)-atmass*x(i)*y(i)
         t(3)=t(3)+atmass*(z(i)**2+x(i)**2)+eps                 
         t(4)=t(4)-atmass*z(i)*x(i)
         t(5)=t(5)-atmass*y(i)*z(i)
         t(6)=t(6)+atmass*(x(i)**2+y(i)**2)+eps                 
   50 continue                                                  
                                                                   
      call rsp(t,3,3,eig,evec)                                     

c   now to orient the molecule so the chirality is preserved       
      xsum=evec(1,1)*(evec(2,2)*evec(3,3)-evec(3,2)*evec(2,3)) +    
     1     evec(1,2)*(evec(2,3)*evec(3,1)-evec(2,1)*evec(3,3)) +    
     2     evec(1,3)*(evec(2,1)*evec(3,2)-evec(2,2)*evec(3,1))      
      if( xsum .lt. 0) then                                         
         do 80 j=1,3                                               
   80    evec(j,1)=-evec(j,1)                                      
      endif                                                        
c     call prmat(6,evec,3,3,'Rmat')
                                                                   
      do 120 i=1,numat                                             
         do 120 j=1,3                                              
            xsum=0.d0                                               
            do 110 k=1,3                                           
  110       xsum=xsum+coordtmp(k,i)*evec(k,j)                           
  120 coord1(j,i)=xsum                                              
                                                                   
      do i=1,numat                                                 
         coordout(1,i)=coord1(1,i)                                    
         coordout(2,i)=coord1(2,i)                                    
         coordout(3,i)=coord1(3,i)                                    
      enddo                                                        
                                                                   
      return               
      end                 




c*****************************************************************         
c transform to CMA of the first n0 atoms
c input coords (all atoms) replaced
c*****************************************************************         
                                               
      subroutine axistrf(numat,n0,nat,coord)             

      use atmasses
      implicit none                                                
      integer numat,n0,nat(numat)
      real*8  coord(3,numat)
      real*8 sumw,sumwx,sumwy,sumwz,atmass,xsum,eps
      real*8 t(6), evec(3,3),eig(3)  
      real*8 x(numat),y(numat),z(numat),coordtmp(3,numat)
      real*8 coord1(3,numat)
      integer i,j,k
      data t /6*0.d0/                                              

      sumw=1.d-20                                                  
      sumwx=0.d0                                                   
      sumwy=0.d0                                                   
      sumwz=0.d0                                                   
      do 10 i=1,n0                                              
         atmass=ams(nat(i))                                     
         sumw=sumw+atmass                                       
         sumwx=sumwx+atmass*coord(1,i)                          
         sumwy=sumwy+atmass*coord(2,i)                          
         sumwz=sumwz+atmass*coord(3,i)                          
   10 continue                                                  
                                                                   
      eps=1.d-3
      sumwx=sumwx/sumw                                             
      sumwy=sumwy/sumw                                             
      sumwz=sumwz/sumw                                             
      do i=1,numat                                              
         x(i)=coord(1,i)-sumwx 
         y(i)=coord(2,i)-sumwy 
         z(i)=coord(3,i)-sumwz 
         coordtmp(1,i)=x(i)
         coordtmp(2,i)=y(i)
         coordtmp(3,i)=z(i)
      enddo
      
      do 40 i=1,6                                                  
   40 t(i)=dble(i)*1.0d-10                                         
                                                                   
      do 50 i=1,n0                                              
         atmass=ams(nat(i))                                     
         t(1)=t(1)+atmass*(y(i)**2+z(i)**2)+eps                 
         t(2)=t(2)-atmass*x(i)*y(i)
         t(3)=t(3)+atmass*(z(i)**2+x(i)**2)+eps                 
         t(4)=t(4)-atmass*z(i)*x(i)
         t(5)=t(5)-atmass*y(i)*z(i)
         t(6)=t(6)+atmass*(x(i)**2+y(i)**2)+eps                 
   50 continue                                                  
                                                                   
      call rsp(t,3,3,eig,evec)                                     

c   now to orient the molecule so the chirality is preserved       
      xsum=evec(1,1)*(evec(2,2)*evec(3,3)-evec(3,2)*evec(2,3)) +    
     1     evec(1,2)*(evec(2,3)*evec(3,1)-evec(2,1)*evec(3,3)) +    
     2     evec(1,3)*(evec(2,1)*evec(3,2)-evec(2,2)*evec(3,1))      
      if( xsum .lt. 0) then                                         
         do 80 j=1,3                                               
   80    evec(j,1)=-evec(j,1)                                      
      endif                                                        
                                                                   
      do 120 i=1,numat                                             
         do 120 j=1,3                                              
            xsum=0.d0                                               
            do 110 k=1,3                                           
  110       xsum=xsum+coordtmp(k,i)*evec(k,j)                           
  120 coord1(j,i)=xsum                                              
                                                                   
      do i=1,numat                                                 
         coord(1,i)=coord1(1,i)                                    
         coord(2,i)=coord1(2,i)                                    
         coord(3,i)=coord1(3,i)                                    
      enddo                                                        
                                                                   
      return               
      end                 

      subroutine axis2(pr,numat,nat,coord,eax)         

      use atmasses
      implicit none                                                
      integer numat,nat(numat)
      real*8  coord(3,numat),eax(3)
      real*8 sumw,sumwx,sumwy,sumwz,atmass,xsum,eps
      real*8 t(6), evec(3,3), eig(3)
      real*8 x(numat),y(numat),z(numat)
      integer i,j,k
      data t /6*0.d0/
      logical pr
      common / atomradii / rad(94)
      real*8 rad
                                                                   
      sumw=1.d-20                                                  
      sumwx=0.d0                                                   
      sumwy=0.d0                                                   
      sumwz=0.d0                                                   
      do 10 i=1,numat                                           
         atmass=ams(nat(i))
         sumw=sumw+atmass                                       
         sumwx=sumwx+atmass*coord(1,i)                          
         sumwy=sumwy+atmass*coord(2,i)                          
         sumwz=sumwz+atmass*coord(3,i)                          
   10 continue                                                  
                                                                   
      eps=1.d-9
      sumwx=sumwx/sumw                                             
      sumwy=sumwy/sumw                                             
      sumwz=sumwz/sumw                                             
      do i=1,numat                                              
         x(i)=coord(1,i)-sumwx 
         y(i)=coord(2,i)-sumwy 
         z(i)=coord(3,i)-sumwz 
      enddo
      
      do 40 i=1,6                                                  
   40 t(i)=dble(i)*1.0d-10                                         
                                                                   
      do 50 i=1,numat                                           
         atmass=ams(nat(i))
         t(1)=t(1)+atmass*(y(i)**2+z(i)**2)+eps                 
         t(2)=t(2)-atmass*x(i)*y(i)
         t(3)=t(3)+atmass*(z(i)**2+x(i)**2)+eps                 
         t(4)=t(4)-atmass*z(i)*x(i)
         t(5)=t(5)-atmass*y(i)*z(i)
         t(6)=t(6)+atmass*(x(i)**2+y(i)**2)+eps                 
   50 continue                                                  
                                                                   
      call rsp(t,3,3,eig,evec)                                     

      eig=1./(eig+eps)**0.25   !0.1666666666  ! 0.25 ???

      sumw=sum(eig)
      eax=eig/sumw
      if(pr)then
        if(eax(2).lt.0.1.and.eax(3).lt.0.1) then
          eax(2:3)=0.4
          write(*,'(7x,''adjusting axis to finite length'')')
        endif
      endif
      sumw=sum(eax)
      eax=eax/sumw
      if(pr)write(*,'(7x,''unit ellipsoid axis a,b,c     :'',3f8.3)')eax

      end

!=============================================================================================================
! Original axis routine adapted from RANAL tool
!=============================================================================================================


!      subroutine axis(numat,nat,coord,aa,bb,cc,avmom,sumw)             
      subroutine axis(numat,nat,coord,rot,avmom)

      use atmasses

      implicit double precision (a-h,o-z)                               
c ams size is wrong  
      dimension coord(3,numat)
      integer nat(numat)
c atomic masses . Products in data statements may cause problems 
                                                                        
      PARAMETER (BOHR=0.52917726)                                       
      dimension t(6), rot(3), xyzmom(3), eig(3), evec(3,3)              
      dimension x(10000),y(10000),z(10000)
      data t /6*0.d0/                                                   
************************************************************************
*     const1 =  10**40/(n*a*a)                                          
*               n = avergadro's number                                  
*               a = cm in an angstrom                                   
*               10**40 is to allow units to be 10**(-40)gram-cm**2      
*                                                                       
************************************************************************
      const1 = 1.66053d0                                                
************************************************************************
*                                                                       
*     const2 = conversion factor from angstrom-amu to cm**(-1)          
*                                                                       
*            = (planck's constant*n*10**16)/(8*pi*pi*c)                 
*            = 6.62618*10**(-27)[erg-sec]*6.02205*10**23*10**16/        
*              (8*(3.1415926535)**2*2.997925*10**10[cm/sec])            
*                                                                     
************************************************************************
      const2=16.8576522d0                                               
c    first we centre the molecule about the centre of gravity,          
c    this depends on the isotopic masses, and the cartesian geometry.   
c                                                                       
      sumw=1.d-20                                                       
      sumwx=0.d0                                                        
      sumwy=0.d0                                                        
      sumwz=0.d0                                                        

      !coord(1:3,1:numat)=coord(1:3,1:numat)*bohr
                                                                        
         do 10 i=1,numat                                                
            atmass=ams(nat(i)) 
            sumw=sumw+atmass                                            
            sumwx=sumwx+atmass*coord(1,i)                               
            sumwy=sumwy+atmass*coord(2,i)
            sumwz=sumwz+atmass*coord(3,i)
   10    continue                                                       
c                                                                       
c     write(*,'(/10x,''molecular weight ='',f8.2,/)')                   
c    2min(99999.99d0,sumw)                                              
      sumwx=sumwx/sumw                                                  
      sumwy=sumwy/sumw                                                  
      sumwz=sumwz/sumw                                                  
      f=1.0d0/bohr
c     write(*,*)'center of mass (ang)',sumwx,sumwy,sumwz
c     write(*,*)'center of mass (au) ',f*sumwx,f*sumwy,f*sumwz
      do 30 i=1,numat                                                   
         x(i)=coord(1,i)-sumwx                                          
         y(i)=coord(2,i)-sumwy                                          
   30 z(i)=coord(3,i)-sumwz                                             

************************************************************************
*    matrix for moments of inertia is of form                           
*                                                                       
*           |   y**2+z**2                         |                     
*           |    -y*x       z**2+x**2             | -i =0               
*           |    -z*x        -z*y       x**2+y**2 |                     
*                                                                       
************************************************************************
      do 40 i=1,6                                                       
   40 t(i)=dble(i)*1.0d-10                                              
         do 50 i=1,numat                                                
            atmass=ams(nat(i))                                          
            t(1)=t(1)+atmass*(y(i)**2+z(i)**2)                          
            t(2)=t(2)-atmass*x(i)*y(i)                                  
            t(3)=t(3)+atmass*(z(i)**2+x(i)**2)                          
            t(4)=t(4)-atmass*z(i)*x(i)                                  
            t(5)=t(5)-atmass*y(i)*z(i)                                  
            t(6)=t(6)+atmass*(x(i)**2+y(i)**2)                          
   50    continue                                                       
      call rsp(t,3,3,eig,evec)                                          
         do 70 i=1,3                                                    
            if(eig(i).lt.3.d-4) then                                    
               eig(i)=0.d0                                              
               rot(i)=0.d0                                              
            else                                                        
               rot(i)=2.9979245d+4*const2/eig(i)
            endif                                                       
   70    xyzmom(i)=eig(i)*const1                                        
c        write(*,'(10x,''(10-47 kg m^2)  a ='',f14.6,''   b ='',f14.6,  
c    1''   c ='',f14.6)')(xyzmom(i),i=1,3)                           
c        write(*,'(10x,''(Mhz         )  a ='',f14.6,''   b ='',f14.6,  
c    1''   c ='',f14.6)')(rot(i),i=1,3)                           
c        write(*,'(10x,''(cm-1        )  a ='',f14.6,''   b ='',f14.6,  
c    1''   c ='',f14.6)')(rot(i)/2.9979245d+4,i=1,3)                    

         aa=rot(1)
         bb=rot(2)
         cc=rot(3)
         avmom=1.d-47*(xyzmom(1)+xyzmom(2)+xyzmom(3))/3.

c        write(*,*) sqrt(avmom/(sumw*1.66053886E-27))
      return                                                            
      end                                                               


