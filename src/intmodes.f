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
c project cartesian on internal mode to determine str, bend, tors
      subroutine intmodestep(n,bmat,u,step,geo,na,nb,nc,coord)
      implicit none
      integer n
      integer na(n),nb(n),nc(n)
      real*8 bmat(3*n-6,3*n),u(3*n),geo(3,n),coord(3,n),step

      real*8 dit(3*n-6),norm,geo2(3,n)
      integer k,j,kl,l,lend,n36
      logical fail

      dit = 0
      do k=1,3*n    
         do j=1,3*n-6
            dit(j)=dit(j)+bmat(j,k)*u(k)
         enddo
      enddo

      kl=0
      geo2=0
      do k=2,n
         lend=3
         if(k.eq.2) lend=1
         if(k.eq.3) lend=2
         do l=1,lend
            kl=kl+1
            geo(l,k)=geo(l,k)+dit(kl)*step
         enddo
      enddo
c     new Cartesians
      call GMETRY(n,geo,coord,na,nb,nc,fail)

      end

c project cartesian on internal mode to determine str, bend, tors
      subroutine modetyp(n,bmat,u,root,vtyp)  
      implicit none
      integer n,root
      real*8 bmat(3*n-6,3*n),u(3*n,3*n),vtyp(3)

      real*8 dit(3*n-6),norm
      integer k,j,kl,l,lend

      dit = 0
      do k=1,3*n    
         do j=1,3*n-6
            dit(j)=dit(j)+bmat(j,k)*u(k,root)
         enddo
      enddo

      vtyp = 0
      kl=0
      do k=2,n
         lend=3
         if(k.eq.2) lend=1
         if(k.eq.3) lend=2
         do l=1,lend
            kl=kl+1
            vtyp(l)=vtyp(l)+dit(kl)**2
         enddo
      enddo

      norm=sum(vtyp)
      vtyp=vtyp/norm

      end

c Bmatrix dZmat/dxyz
      subroutine bzmat(n,at,xyzin,bmat)  
      implicit none
      integer n,at(n)
      real*8 xyzin(3,n)
      real*8 bmat(3*n-6,3*n)

      real*8 xyz(3,n),geo(3,n),br(3*n-6),bl(3*n-6),step,one
      parameter (one=1.0d0)
      integer i,k,n3,ia,ii,ic,na(n),nb(n),nc(n),n36,j,kl
      integer lend,l,ij

      n36=3*n-6
      n3=3*n
      xyz=xyzin*0.52916790d0
      step=1.d-5

c Eq. zmat
      geo = 0
      call XYZINT(XYZ,n,NA,NB,NC,one,GEO)

c make Bmat
      ij=0
      do i=1,n
         do j=1,3
            ij=ij+1
            xyz(j,i)=xyz(j,i)+step
            call XYZGEO(XYZ,n,NA,NB,NC,one,GEO)
            kl=0
            do k=2,n
               lend=3
               if(k.eq.2) lend=1
               if(k.eq.3) lend=2
               do l=1,lend
                  kl=kl+1
                  br(kl)=geo(l,k)
               enddo
            enddo
            xyz(j,i)=xyz(j,i)-step*2.
            call XYZGEO(XYZ,n,NA,NB,NC,one,GEO)
            xyz(j,i)=xyz(j,i)+step
            kl=0
            do k=2,n
               lend=3
               if(k.eq.2) lend=1
               if(k.eq.3) lend=2
               do l=1,lend
                  kl=kl+1
                  bl(kl)=geo(l,k)
               enddo
            enddo
            do k=1,n36
               bmat(k,ij)=0.5*(br(k)-bl(k))/step
            enddo
         enddo
      enddo

      end


C     *****************************************************************

      SUBROUTINE makenabc(XYZin,molvec,at,nat,N,n2,nmol,
     .                    fragind,NA,NB,NC)
      IMPLICIT none                       
      integer n,n2,NA(n2), NB(n2), NC(n2), at(n), nmol, fragind(100)
      integer molvec(n),nat(n2)
      real*8 xyzin(3,n)

      integer ind(n2,n2),i,j,k,m,molv(n2),nm,ntot,mm,idum(n)
      real*8 xyz(3,n2),R,rr(n2),sum3(3),tmp(3,n2),one
      parameter (one=1.0d0)

      if(nmol.gt.1)then
      k=1
      ntot=0
      do nm=1,nmol  ! loop over fragments
         m=0
         do i=k,fragind(nm)
            m=m+1
            tmp(1:3,m)=xyzin(1:3,i)
            idum(m)=at(i)
         enddo
         call cma(m,idum,tmp,sum3)   
         ntot=ntot+1
         nat(ntot)=107
         molv(ntot)=nm            
         xyz(1:3,ntot)=sum3(1:3)
         do i=k,fragind(nm)
            ntot=ntot+1
            xyz(1:3,ntot)=xyzin(1:3,i)
            molv(ntot)=molvec(i)
            nat(ntot)=at(i)
         enddo
         k=k+fragind(nm)
      enddo
      else
      nat=at
      molv=molvec
      xyz=xyzin
      endif

      do i=1,n2
         do j=1,n2     
            R=(XYZ(1,I)-XYZ(1,J))**2+
     1        (XYZ(2,I)-XYZ(2,J))**2+
     2        (XYZ(3,I)-XYZ(3,J))**2
            rr(j)=R+abs(molv(i)-molv(j))*100000.
            if(i.eq.j) rr(j)=1.d+42
            ind(j,i)=j
         enddo
         call qsort(rr,1,n,ind(1,i))
      enddo
      na(2)=1
      call bonded(3,na,nb,nc,ind,n2,molv)
      NA(1)=0
      NB(1)=0
      NC(1)=0
      NB(2)=0
      NC(2)=0
      NC(3)=0

      write(*,*) 'Z-matrix connectivity:'
      mm=1
      do i=1,n2
         na(i)=0
         if(i.gt.1)then
            if(molv(i).ne.molv(i-1)) then
               na(i)=mm
               mm=na(i)
            endif
         endif
         call bonded(i,na,nb,nc,ind,n2,molv)
c        write(*,*)'atom ',i,'  NA,NB,NC :',na(i),nb(i),nc(i)
      enddo

      call XYZINT(XYZ,n2,NA,NB,NC,1.d0,tmp)
      call xyzgeo(xyz,n2,na,nb,nc,one,tmp)
      call zmatpr(n2,nat,tmp,na,nb,nc,1)

      RETURN
      END

C     *****************************************************************

      SUBROUTINE cart2zmat(XYZin,molvec,at,nat,N,n2,nmol,
     .                    fragind,NA,NB,NC,geo)
      IMPLICIT none                       
      integer n,n2,NA(n2), NB(n2), NC(n2), at(n), nmol, fragind(100)
      integer molvec(n),nat(n2)
      real*8 xyzin(3,n),geo(3,n2)

      integer ind(n2,n2),i,j,k,m,molv(n2),nm,ntot,mm,idum(n)
      real*8 xyz(3,n2),R,rr(n2),sum3(3),tmp(3,n),one
      parameter (one=1.0d0)

      if(nmol.gt.1)then
      k=1
      ntot=0
      do nm=1,nmol  ! loop over fragments
         m=0
         do i=k,fragind(nm)
            m=m+1
            tmp(1:3,m)=xyzin(1:3,i)
            idum(m)=at(i)
         enddo
         call cma(m,idum,tmp,sum3)   
         ntot=ntot+1
         xyz(1:3,ntot)=sum3(1:3)
         do i=k,fragind(nm)
            ntot=ntot+1
            xyz(1:3,ntot)=xyzin(1:3,i)
         enddo
         k=k+fragind(nm)
      enddo
      else
      xyz=xyzin
      endif

      call xyzgeo(xyz,n2,na,nb,nc,one,geo)

      RETURN
      END

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine bonded(i,na,nb,nc,ind,n,molvec)
      implicit none
      integer i,n,ind(n,n),j,is,molvec(n),icon,na(n),nb(n),nc(n)
 
      nb(i)=0
      nc(i)=0
      is=0
      if(na(i).ne.0) goto 10
      do j=1,n
         icon=ind(j,i)
         if(icon.lt.i.and.molvec(icon).eq.molvec(i)) then
            is=j
            na(i)=icon
            goto 10
         endif
      enddo
 10   continue
      if(na(i).eq.0)then
      do j=1,n
         icon=ind(j,i)
         if(icon.lt.i) then
            is=j
            na(i)=icon
            goto 20
         endif
      enddo
 20   continue
      endif
      do j=is+1,n
         icon=ind(j,i)
         if(icon.lt.i.and.molvec(icon).eq.molvec(i)
     .   .and.icon.ne.na(i)) then
            is=j
            nb(i)=icon
            goto 30
         endif
      enddo
 30   continue
      if(nb(i).eq.0)then
      do j=is+1,n
         icon=ind(j,i)
         if(icon.lt.i
     .   .and.icon.ne.na(i)) then
            is=j
            nb(i)=icon
            goto 40
         endif
      enddo
      endif                 
 40   continue
      do j=is+1,n
         icon=ind(j,i)
         if(icon.lt.i.and.molvec(icon).eq.molvec(i)
     .   .and.icon.ne.na(i).and.icon.ne.nb(i)) then
            is=j
            nc(i)=icon
            goto 50
         endif
      enddo
 50   continue
      if(nc(i).eq.0)then
      do j=is+1,n
         icon=ind(j,i)
         if(icon.lt.i
     .   .and.icon.ne.na(i).and.icon.ne.nb(i)) then
            is=j
            nc(i)=icon
            return
         endif
      enddo
      endif

      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine zmat2cart(n,n2,at,geo,xyz,na,nb,nc,fail) 
      implicit none
      integer n,n2,na(n2),nb(n2),nc(n2),at(n2)
      real*8 xyz(3,n),geo(3,n2)
      logical fail

      integer k,m
      real*8 tmp(3,n2)

      call GMETRY(n2,geo,tmp,na,nb,nc,fail)

      m=0
      do k=1,n2
         if(at(k).lt.100) then
            m=m+1
            xyz(1:3,m)=tmp(1:3,k)
         endif
      enddo

      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c     subroutine int2cart(nall,s1,s2,xyzref,xyz,n,nmol,fragind) 
c     use ls_rmsd
c     implicit none
c     integer n,nall,s1,s2,fragind(100),nmol
c     real*8 xyzref(3,n,nall),xyz(3,n)

c     integer nm,i,j,k,m,mm
c     real*8,allocatable::tmp(:,:),tmp2(:,:)
c     real*8 gdum(3,3), U(3,3), x(3), y(3), dum
c
c     k=1
c     mm=0
c     do nm=1,nmol  ! loop over fragments
c        m=0
c        do i=k,fragind(nm)
c           m=m+1
c        enddo
c        mm=m
c        allocate(tmp(3,mm),tmp2(3,mm))
c        m=0
c        do i=k,fragind(nm)
c           m=m+1
c           tmp (1:3,m)=xyzref(1:3,i,1)
c           tmp2(1:3,m)=xyzref(1:3,i,s1)                 ! copy bad Zmatrix
c        enddo
c        call rmsd(mm,tmp,tmp2,1,U1,x,y1,dum,.false.,gdum) ! det orient
c        m=0
c        do i=k,fragind(nm)
c           m=m+1
c           tmp (1:3,m)=xyzref(1:3,i,1)
c           tmp2(1:3,m)=xyzref(1:3,i,s2)                 ! copy bad Zmatrix
c        enddo
c        call rmsd(mm,tmp,tmp2,1,U2,x,y2,dum,.false.,gdum) ! det orient
c        write(*,*) nm,dum
c        call prmat(6,U,3,3,'U')
c        write(*,*) x
c        write(*,*) y
c        do i=1,mm                                       ! shift
c           tmp2(1:3,i)=tmp2(1:3,i)-y(1:3)
c        enddo
c        tmp2=matmul(U,tmp2)                             ! rotate
c        m=0
c        do i=k,fragind(nm)
c           m=m+1
c           xyz(1:3,i)=tmp2(1:3,m)+x(1:3)                ! shift back
c        enddo
c        deallocate(tmp,tmp2)
c        k=k+fragind(nm)
c     enddo

c     end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE XYZGEO(XYZ,NUMAT,NA,NB,NC,DEGREE,GEO)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION XYZ(3,*), NA(*), NB(*), NC(*), GEO(3,*)
***********************************************************************
*
*   XYZGEO CONVERTS COORDINATES FROM CARTESIAN TO INTERNAL.
*
*     ON INPUT XYZ  = ARRAY OF CARTESIAN COORDINATES
*              NUMAT= NUMBER OF ATOMS
*              NA   = NUMBERS OF ATOM TO WHICH ATOMS ARE RELATED
*                     BY DISTANCE
*              NB   = NUMBERS OF ATOM TO WHICH ATOMS ARE RELATED
*                     BY ANGLE
*              NC   = NUMBERS OF ATOM TO WHICH ATOMS ARE RELATED
*                     BY DIHEDRAL
*
*    ON OUTPUT GEO  = INTERNAL COORDINATES IN ANGSTROMS, RADIANS,
*                     AND RADIANS
*
***********************************************************************
      DO 10 I=2,NUMAT
         J=NA(I)
         K=NB(I)
         L=NC(I)
         IF(I.LT.3) GOTO 10
         II=I
         CALL BANGLE(XYZ,II,J,K,GEO(2,I))
         GEO(2,I)=GEO(2,I)*DEGREE
         IF(I.LT.4) GOTO 10
C
C   MAKE SURE DIHEDRAL IS MEANINGLFUL
C
c        CALL BANGLE(XYZ,J,K,L,ANGL)
c        TOL=0.2617994D0
c        IF(ANGL.GT.3.1415926D0-TOL.OR.ANGL.LT.TOL)THEN
C
C  ANGLE IS UNSATISFACTORY, LET'S SEARCH FOR ANOTHER ATOM FOR
C  DEFINING THE DIHEDRAL.
c 88       SUM=100.D0
c        DO 74 I1=1,II-1
c           R=(XYZ(1,I1)-XYZ(1,K))**2+
c    1          (XYZ(2,I1)-XYZ(2,K))**2+
c    2          (XYZ(3,I1)-XYZ(3,K))**2
c           IF(R.LT.SUM.AND.I1.NE.J.AND.I1.NE.K) THEN
c         CALL BANGLE(XYZ,J,K,I1,ANGL)
c        IF(ANGL.LT.3.1415926D0-TOL.OR.ANGL.GT.TOL)THEN
c              SUM=R
c              L=I1
c              NC(II)=L
c           ENDIF
c           ENDIF
c 74  CONTINUE
c     IF(SUM.GT.99.D0.AND.TOL.GT.0.1D0)THEN
c
c ANYTHING WITHIN 5 DEGREES?
c
c     TOL=0.087266D0
c     GOTO 88
c     ENDIF
c           ENDIF
         CALL DIHED(XYZ,II,J,K,L,GEO(3,I))
         GEO(3,I)=GEO(3,I)*DEGREE
10    GEO(1,I)= SQRT((XYZ(1,I)-XYZ(1,J))**2+
     1               (XYZ(2,I)-XYZ(2,J))**2+
     2               (XYZ(3,I)-XYZ(3,J))**2)
      GEO(1,1)=0.D0
      GEO(2,1)=0.D0
      GEO(3,1)=0.D0
      GEO(2,2)=0.D0
      GEO(3,2)=0.D0
      GEO(3,3)=0.D0
      RETURN
      END

C     *****************************************************************
      SUBROUTINE XYZINT(XYZ,NUMAT,NA,NB,NC,DEGREE,GEO)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION XYZ(3,*), NA(*), NB(*), NC(*), GEO(3,*)
***********************************************************************
*
* XYZINT WORKS OUT THE INTERNAL COORDINATES OF A MOLECULE.
*        THE "RULES" FOR THE CONNECTIVITY ARE AS FOLLOWS:
*        ATOM I IS DEFINED AS BEING AT A DISTANCE FROM THE NEAREST
*        ATOM J, ATOM J ALREADY HAVING BEEN DEFINED.
*        ATOM I MAKES AN ANGLE WITH ATOM J AND THE ATOM K, WHICH HAS
*        ALREADY BEEN DEFINED, AND IS THE NEAREST ATOM TO J
*        ATOM I MAKES A DIHEDRAL ANGLE WITH ATOMS J, K, AND L. L HAVING
*        BEEN DEFINED AND IS THE NEAREST ATOM TO K, AND J, K AND L
*        HAVE A CONTAINED ANGLE IN THE RANGE 15 TO 165 DEGREES,
*        IF POSSIBLE.
*
*        NOTE THAT GEO AND XYZ MUST NOT BE THE SAME IN THE CALL.
*
*   ON INPUT XYZ    = CARTESIAN ARRAY OF NUMAT ATOMS
*            DEGREE = 1 IF ANGLES ARE TO BE IN RADIANS
*            DEGREE = 57.29578 IF ANGLES ARE TO BE IN RADIANS
*
***********************************************************************
      NAI1=0
      NAI2=0
      DO 20 I=1,NUMAT
         NA(I)=2
         NB(I)=3
         NC(I)=4
         IM1=I-1
         IF(IM1.EQ.0)GOTO 20
         SUM=100.D0
         DO 10 J=1,IM1
            R=(XYZ(1,I)-XYZ(1,J))**2+
     1          (XYZ(2,I)-XYZ(2,J))**2+
     2          (XYZ(3,I)-XYZ(3,J))**2
            IF(R.LT.SUM.AND.NA(J).NE.J.AND.NB(J).NE.J) THEN
               SUM=R
               K=J
            ENDIF
   10    CONTINUE
C
C   ATOM I IS NEAREST TO ATOM K
C
         NA(I)=K
         IF(I.GT.2)NB(I)=NA(K)
         IF(I.GT.3)NC(I)=NB(K)
C
C   FIND ANY ATOM TO RELATE TO NA(I)
C
   20 CONTINUE
      NA(1)=0
      NB(1)=0
      NC(1)=0
      NB(2)=0
      NC(2)=0
      NC(3)=0
      
      CALL XYZGEO(XYZ,NUMAT,NA,NB,NC,DEGREE,GEO)
      RETURN
      END

      SUBROUTINE BANGLE(XYZ,I,J,K,ANGLE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION XYZ(3,*)
*********************************************************************
*
* BANGLE CALCULATES THE ANGLE BETWEEN ATOMS I,J, AND K. THE
*        CARTESIAN COORDINATES ARE IN XYZ.
*
*********************************************************************
      D2IJ = (XYZ(1,I)-XYZ(1,J))**2+
     1       (XYZ(2,I)-XYZ(2,J))**2+
     2       (XYZ(3,I)-XYZ(3,J))**2
      D2JK = (XYZ(1,J)-XYZ(1,K))**2+
     1       (XYZ(2,J)-XYZ(2,K))**2+
     2       (XYZ(3,J)-XYZ(3,K))**2
      D2IK = (XYZ(1,I)-XYZ(1,K))**2+
     1       (XYZ(2,I)-XYZ(2,K))**2+
     2       (XYZ(3,I)-XYZ(3,K))**2
      XY = SQRT(D2IJ*D2JK+1.d-14)
      TEMP = 0.5D0 * (D2IJ+D2JK-D2IK) / XY
      IF (TEMP .GT. 1.0D0) TEMP=1.0D0
      IF (TEMP .LT. -1.0D0) TEMP=-1.0D0
      ANGLE = ACOS( TEMP )
      RETURN
      END
      SUBROUTINE DIHED(XYZ,I,J,K,L,ANGLE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION XYZ(3,*)
*********************************************************************
*
*      DIHED CALCULATES THE DIHEDRAL ANGLE BETWEEN ATOMS I, J, K,
*            AND L.  THE CARTESIAN COORDINATES OF THESE ATOMS
*            ARE IN ARRAY XYZ.
*
*     DIHED IS A MODIFIED VERSION OF A SUBROUTINE OF THE SAME NAME
*           WHICH WAS WRITTEN BY DR. W. THEIL IN 1973.
*
*********************************************************************
      DATA PI/3.14159265358979D0/
      XI1=XYZ(1,I)-XYZ(1,K)
      XJ1=XYZ(1,J)-XYZ(1,K)
      XL1=XYZ(1,L)-XYZ(1,K)
      YI1=XYZ(2,I)-XYZ(2,K)
      YJ1=XYZ(2,J)-XYZ(2,K)
      YL1=XYZ(2,L)-XYZ(2,K)
      ZI1=XYZ(3,I)-XYZ(3,K)
      ZJ1=XYZ(3,J)-XYZ(3,K)
      ZL1=XYZ(3,L)-XYZ(3,K)
C      ROTATE AROUND Z AXIS TO PUT KJ ALONG Y AXIS
      DIST= SQRT(XJ1**2+YJ1**2+ZJ1**2+1.d-14)
      COSA=ZJ1/DIST
      IF(COSA.GT.1.0D0) COSA=1.0D0
      IF(COSA.LT.-1.0D0) COSA=-1.0D0
      DDD=1.0D0-COSA**2
      IF(DDD.LE.0.0) GO TO 10
      YXDIST=DIST* SQRT(DDD)
      IF(YXDIST.GT.1.0D-12) GO TO 20
   10 CONTINUE
      XI2=XI1
      XL2=XL1
      YI2=YI1
      YL2=YL1
      COSTH=COSA
      SINTH=0.D0
      GO TO 30
   20 COSPH=YJ1/YXDIST
      SINPH=XJ1/YXDIST
      XI2=XI1*COSPH-YI1*SINPH
      XJ2=XJ1*COSPH-YJ1*SINPH
      XL2=XL1*COSPH-YL1*SINPH
      YI2=XI1*SINPH+YI1*COSPH
      YJ2=XJ1*SINPH+YJ1*COSPH
      YL2=XL1*SINPH+YL1*COSPH
C      ROTATE KJ AROUND THE X AXIS SO KJ LIES ALONG THE Z AXIS
      COSTH=COSA
      SINTH=YJ2/DIST
   30 CONTINUE
      YI3=YI2*COSTH-ZI1*SINTH
      YL3=YL2*COSTH-ZL1*SINTH
      CALL DANG(XL2,YL3,XI2,YI3,ANGLE)
c     IF (ANGLE .LT. 0.) ANGLE=2.0D0*PI+ANGLE
csg   IF (ANGLE .GE. 2.0d0*PI    ) ANGLE=0.D0
      RETURN
      END
      SUBROUTINE DANG(A1,A2,B1,B2,RCOS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
**********************************************************************
*
*    DANG  DETERMINES THE ANGLE BETWEEN THE POINTS (A1,A2), (0,0),
*          AND (B1,B2).  THE RESULT IS PUT IN RCOS.
*
**********************************************************************
      DATA PI/3.14159265358979D0/
      ZERO=1.0D-10
      IF( ABS(A1).LT.ZERO.AND. ABS(A2).LT.ZERO) GO TO 10
      IF( ABS(B1).LT.ZERO.AND. ABS(B2).LT.ZERO) GO TO 10
      ANORM=1.0D0/ SQRT(A1**2+A2**2)
      BNORM=1.0D0/ SQRT(B1**2+B2**2)
      A1=A1*ANORM
      A2=A2*ANORM
      B1=B1*BNORM
      B2=B2*BNORM
      SINTH=(A1*B2)-(A2*B1)
      COSTH=A1*B1+A2*B2
      IF(COSTH.GT.1.0D0) COSTH=1.0D0
      IF(COSTH.LT.-1.0D0) COSTH=-1.0D0
      RCOS= ACOS(COSTH)
      IF( ABS(RCOS).LT.ZERO) GO TO 10
      IF(SINTH.GT.0.D0) RCOS=2.0D0*PI-RCOS
      RCOS=-RCOS
      RETURN
   10 RCOS=0.0D0
      RETURN
      END

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE GMETRY(NATOMS, geo, coord, na,nb,nc,fail)
      IMPLICIT REAL*8 (A-H,O-Z)

      DIMENSION GEO(3,*), coord(3,*), na(*), nb(*), nc(*)
      logical fail

      fail=.false.

      COORD(1,1)=0.0D00
      COORD(2,1)=0.0D00
      COORD(3,1)=0.0D00
      COORD(1,2)=GEO(1,2)
      COORD(2,2)=0.0D00
      COORD(3,2)=0.0D00
      IF(NATOMS.EQ.2) GOTO 110
      CCOS=COS(GEO(2,3))
      IF(NA(3).EQ.1)THEN
         COORD(1,3)=COORD(1,1)+GEO(1,3)*CCOS
      ELSE
         COORD(1,3)=COORD(1,2)-GEO(1,3)*CCOS
      ENDIF
      COORD(2,3)=GEO(1,3)*SIN(GEO(2,3))
      COORD(3,3)=0.0D00
      DO 100 I=4,NATOMS
         COSA=COS(GEO(2,I))
         MB=NB(I)
         MC=NA(I)
         XB=COORD(1,MB)-COORD(1,MC)
         YB=COORD(2,MB)-COORD(2,MC)
         ZB=COORD(3,MB)-COORD(3,MC)
         RBC=1.0D00/DSQRT(XB*XB+YB*YB+ZB*ZB)
         IF (ABS(COSA).LT.0.9999999999D0) GO TO 40
C
C     ATOMS MC, MB, AND (I) ARE COLLINEAR
C
         RBC=GEO(1,I)*RBC*COSA
         COORD(1,I)=COORD(1,MC)+XB*RBC
         COORD(2,I)=COORD(2,MC)+YB*RBC
         COORD(3,I)=COORD(3,MC)+ZB*RBC
         GO TO 100
C
C     THE ATOMS ARE NOT COLLINEAR
C
   40    MA=NC(I)
         XA=COORD(1,MA)-COORD(1,MC)
         YA=COORD(2,MA)-COORD(2,MC)
         ZA=COORD(3,MA)-COORD(3,MC)
C
C     ROTATE ABOUT THE Z-AXIS TO MAKE YB=0, AND XB POSITIVE.  IF XYB IS
C     TOO SMALL, FIRST ROTATE THE Y-AXIS BY 90 DEGREES.
C
         XYB=DSQRT(XB*XB+YB*YB)
         K=-1
         IF (XYB.GT.0.1D00) GO TO 50
         XPA=ZA
         ZA=-XA
         XA=XPA
         XPB=ZB
         ZB=-XB
         XB=XPB
         XYB=DSQRT(XB*XB+YB*YB)
         K=1
C
C     ROTATE ABOUT THE Y-AXIS TO MAKE ZB VANISH
C
   50    COSTH=XB/XYB
         SINTH=YB/XYB
         XPA=XA*COSTH+YA*SINTH
         YPA=YA*COSTH-XA*SINTH
         SINPH=ZB*RBC
         COSPH=DSQRT(ABS(1.D00-SINPH*SINPH))
         XQA=XPA*COSPH+ZA*SINPH
         ZQA=ZA*COSPH-XPA*SINPH
C
C     ROTATE ABOUT THE X-AXIS TO MAKE ZA=0, AND YA POSITIVE.
C
         YZA=DSQRT(YPA*YPA+ZQA*ZQA)
         IF(YZA.LT.1.D-4 )THEN
            IF(YZA.LT.1.D-4)GOTO 70
c           WRITE(*,'(/9X,'' ATOMS'',I3,'','',I3,'', AND'',I3,
c    1'' ARE WITHIN'',F7.4,'' ANGSTROMS OF A STRAIGHT LINE'')')
c    2MC,MB,MA,YZA
            fail=.true.
            return
         ENDIF
         COSKH=YPA/YZA
         SINKH=ZQA/YZA
         GOTO 80
   70    CONTINUE
C
C   ANGLE TOO SMALL TO BE IMPORTANT
C
         COSKH=1.D0
         SINKH=0.D0
   80    CONTINUE
C
C     COORDINATES :-   A=(XQA,YZA,0),   B=(RBC,0,0),  C=(0,0,0)
C     NONE ARE NEGATIVE.
C     THE COORDINATES OF I ARE EVALUATED IN THE NEW FRAME.
C
         SINA=SIN(GEO(2,I))
         SIND=-SIN(GEO(3,I))
         COSD=COS(GEO(3,I))
         XD=GEO(1,I)*COSA
         YD=GEO(1,I)*SINA*COSD
         ZD=GEO(1,I)*SINA*SIND
C
C     TRANSFORM THE COORDINATES BACK TO THE ORIGINAL SYSTEM.
C
         YPD=YD*COSKH-ZD*SINKH
         ZPD=ZD*COSKH+YD*SINKH
         XPD=XD*COSPH-ZPD*SINPH
         ZQD=ZPD*COSPH+XD*SINPH
         XQD=XPD*COSTH-YPD*SINTH
         YQD=YPD*COSTH+XPD*SINTH
         IF (K.LT.1) GO TO 90
         XRD=-ZQD
         ZQD=XQD
         XQD=XRD
   90    COORD(1,I)=XQD+COORD(1,MC)
         COORD(2,I)=YQD+COORD(2,MC)
         COORD(3,I)=ZQD+COORD(3,MC)
  100 CONTINUE
  110 CONTINUE

      RETURN
      END

      subroutine cma(numat,nat,coord,sum3)   

c atomic masses
      use atmasses

      implicit double precision (a-h,o-z)                      
      dimension coord(3,*),nat(*), sum3(3)                     
                                                                 
      sumw=0                                                    
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

      sum3(1)=sumwx/sumw                                           
      sum3(2)=sumwy/sumw                                            
      sum3(3)=sumwz/sumw                                             

      end
