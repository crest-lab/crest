CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C non-linear Fit by Levenberg-Marquart Algorithm
C based on a Pascal program by S.G from the early 90s
C S.G., 10/08
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

c the function f(x):
      real*8 function f_marq(p,x)
      real*8 x,p(*)

      f_marq=p(1)*(1.0d0-exp(-p(2)*x**p(3)))

      end

      subroutine Marqfit(pr,pr2,n_in,n_use,y_in,nconf_exp,val,rmsd,ok)
      implicit none
      logical :: pr,pr2   ! print flag
      integer n_in        ! # points
      integer n_use       ! # points used in fit
      real*8 y_in(*)      ! data
      real*8 nconf_exp    ! fixed parameter
      real*8 val          ! extrapolated value
      real*8 rmsd         ! of fit
      logical ok          ! fit ok ?

      real*8, allocatable :: Nmat(:,:),fx(:),x(:),y(:)
      real*8, allocatable :: rr(:),ll(:),t(:,:),a(:,:),dp(:),p(:)
      real*8 s0,s1,dsq,inkre,dp1,eps,sumx,det,f_marq,r2
      integer np
      integer n 
      integer i,j,kk,L,iter

c # parameters
      np=2
      if(n_in-1 .gt. n_use) then
        n = n_use
        j = n_in - n_use + 1    ! first value used for fit
      else
        n = n_in - 1 
        j = 2
      endif

      allocate(Nmat(np,np),fx(n),p(np+1),x(n),y(n),
     .         rr(n),ll(n),t(np,n),a(n,np),dp(np))

      do i=1,n
         x(i) = float(i)   ! taker iter blocks
      enddo
      y(1:n) = y_in(j:n_in) - y_in(1) ! subtract base (MTD=iter0) value

c initial set of parameters 
      p(1)=y(n)+y(n)*0.1 !just start values
      p(2)=0.3            ! "     "     "
      p(3)=nconf_exp      ! nconf_exp is a fixed parameter

      if(pr)then
      write(*,*) 'initial P:',p(1:np+1)
      do i=1,n
         write(*,*) x(i),y(i),f_marq(p,x(i))
      enddo
      endif
      
      eps=1.d-4  ! exit thr
      inkre=0.5  ! damping

      s1 =1.d+42
      do iter=1,100

         do i=1,n
            fx(i)=f_marq(p,x(i))
         enddo
         s0=s1
         s1=0
         do i=1,n
            s1=s1+(fx(i)-y(i))**2
         enddo
         dsq=s1-s0

c        if(pr)write(*,*) iter,s1,dsq/s1
         if(abs(dsq/s1).lt.eps.or.s1.lt.eps) exit

         do kk=1,np
            dp1=abs(0.001*p(kk))
            if(dp1.lt.1.d-8) dp1=1.d-8
            p(kk)=p(kk)+dp1
            do i=1,n
               rr(i)=f_marq(p,x(i))
            enddo
            p(kk)=p(kk)-2.0*dp1
            do i=1,n
               ll(i)=f_marq(p,x(i))
            enddo
            p(kk)=p(kk)+dp1
            do i=1,n
               a(i,kk)=(rr(i)-ll(i))/(2.*dp1)
            enddo
         enddo

         do i=1,n
            do j=1,np
               t(j,i)=a(i,j)
            enddo
         enddo
c                                               (*Matrizenmultiplikation *)
c                                               (* T^*A^ = Nmat          *)
         do i = 1, np 
            do j = 1, np 
            sumx = 0
            do L = 1,   n 
               sumx = sumx+t(i, L)*a(L, j)
            enddo
            Nmat(i, j) = sumx
            Enddo
         Enddo

c damping of diagonal
         do i=1,np
            Nmat(i,i)=Nmat(i,i)+inkre*Nmat(i,i)
         enddo
         inkre=inkre/1.1
c                                            (*  T^ * residuen = dp*)
         do i = 1, np 
            sumx = 0
            do L = 1, n 
               sumx = sumx+T(i, L)*(y(L)-fx(L))
            Enddo
            dp(i) = sumx
         Enddo
c                                                (* Matrix invertieren *)
         call minv(Nmat,np,det)
         if(det.lt.1.d-12) then
            if(pr)write(*,*)'no matrix inversion possible. exit.'
            goto 999
         endif
c                                                (* N^-1 * (transp.*residuen) *)

         do i = 1, np 
            sumx = 0
            do L = 1, np 
               sumx = sumx+Nmat(i, L)*dp(L)
            Enddo
c                                                 (* neue Parameter *)
            p(i)=p(i)+sumx    
         Enddo


      enddo

999   rmsd = sqrt(s1/float(n))
      s0  =0
      sumx=0
      dsq =0
      do i=1,n
         sumx=sumx+(y(i)-f_marq(p,x(i)))**2
         s0  =s0  +abs(y(i)-f_marq(p,x(i)))/(y(i)+1.d-6)
         dsq =dsq +(y(i)-s0)**2
      enddo
      r2 = 1.0d0-sumx/(dsq+1.d-9)
      s0 = 100.*s0 / float(n)
      s1 = 100.*abs(p(1)-y(n)) / y(n)

!     if(pr)then
!     do i=1,n
!        write(42,*) i,y(i)
!     enddo
!     write(42,*) 
!     do i=1,n
!        write(42,*) i,f_marq(p,x(i))
!     enddo
!     write(42,*) 
!     write(42,*) '0 ',p(1)
!     write(42,*) i+1,p(1)
!     endif

      ok = .true. 
      if( s0 .gt. 4.00 ) ok = .false.! result trustworthy if av. rel. dev of fit < 4 % 
      if( s1 .gt.100.0 ) ok = .false.! extrapolated value very different from last  
      if( r2 .lt. 0.95 ) ok = .false.! bad correlation

      if(pr)then
      write(*,'(''P         :'',6F12.6)')(p(kk),kk=1,np),nconf_exp
      write(*,'(''RMSD      :'', F12.6)') rmsd
      write(*,'(''av % dev  :'', F12.6)') s0
      write(*,'(''% dev last:'', F12.6)') s1
      write(*,'(''ok fit    :'', L    )') ok     
      endif

      if(pr2) then
        write(*,'(''   block      Siter   Siter (fit) '')') 
        do i=1,n
           write(*,'(i9,2F9.3)') int(x(i)),y(i),f_marq(p,x(i))
        enddo
        write(*,'(''R²       :'', f12.6)') r2             
      endif

      val = p(1) + y_in(1) ! add base value

      end

      SUBROUTINE MINV(A,N,D)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(*)
**********************************************************************
*
*     INVERT A MATRIX USING GAUSS-JORDAN METHOD.  PART OF DIIS
*     A - INPUT MATRIX (MUST BE A GENERAL MATRIX), DESTROYED IN
*        COMPUTATION AND REPLACED BY RESULTANT INVERSE.
*     N - ORDER OF MATRIX A
*     D - RESULTANT DETERMINANT
*
**********************************************************************
      integer M(N), L(N)
C
C     SEARCH FOR LARGEST ELEMENT
C
      D=1.0D0
      NK=-N
      DO 180 K=1,N
         NK=NK+N
         L(K)=K
         M(K)=K
         KK=NK+K
         BIGA=A(KK)
         DO 20 J=K,N
            IZ=N*(J-1)
            DO 20 I=K,N
               IJ=IZ+I
   10          IF (ABS(BIGA).LT.ABS(A(IJ)))THEN
                  BIGA=A(IJ)
                  L(K)=I
                  M(K)=J
               ENDIF
   20    CONTINUE
C
C     INTERCHANGE ROWS
C
         J=L(K)
         IF (J-K) 50,50,30
   30    KI=K-N
         DO 40 I=1,N
            KI=KI+N
            HOLD=-A(KI)
            JI=KI-K+J
            A(KI)=A(JI)
   40    A(JI)=HOLD
C
C     INTERCHANGE COLUMNS
C
   50    I=M(K)
         IF (I-K) 80,80,60
   60    JP=N*(I-1)
         DO 70 J=1,N
            JK=NK+J
            JI=JP+J
            HOLD=-A(JK)
            A(JK)=A(JI)
   70    A(JI)=HOLD
C
C     DIVIDE COLUMN BY MINUS PIVOT (VALUE OF PIVOT ELEMENT IS
C     CONTAINED IN BIGA)
C
   80    IF (BIGA) 100,90,100
   90    D=0.0
         RETURN
  100    DO 120 I=1,N
            IF (I-K) 110,120,110
  110       IK=NK+I
            A(IK)=A(IK)/(-BIGA)
  120    CONTINUE
C  REDUCE MATRIX
         DO 150 I=1,N
            IK=NK+I
            HOLD=A(IK)
            IJ=I-N
            DO 150 J=1,N
               IJ=IJ+N
               IF (I-K) 130,150,130
  130          IF (J-K) 140,150,140
  140          KJ=IJ-I+K
               A(IJ)=HOLD*A(KJ)+A(IJ)
  150    CONTINUE
C
C     DIVIDE ROW BY PIVOT
C
         KJ=K-N
         DO 170 J=1,N
            KJ=KJ+N
            IF (J-K) 160,170,160
  160       A(KJ)=A(KJ)/BIGA
  170    CONTINUE
C
C     PRODUCT OF PIVOTS
C
         D=MAX(-1.D25,MIN(1.D25,D))
         D=D*BIGA
C
C     REPLACE PIVOT BY RECIPROCAL
C
         A(KK)=1.0/BIGA
  180 CONTINUE
C
C     FINAL ROW AND COLUMN INTERCHANGE
C
      K=N
  190 K=(K-1)
      IF (K) 260,260,200
  200 I=L(K)
      IF (I-K) 230,230,210
  210 JQ=N*(K-1)
      JR=N*(I-1)
      DO 220 J=1,N
         JK=JQ+J
         HOLD=A(JK)
         JI=JR+J
         A(JK)=-A(JI)
  220 A(JI)=HOLD
  230 J=M(K)
      IF (J-K) 190,190,240
  240 KI=K-N
      DO 250 I=1,N
         KI=KI+N
         HOLD=A(KI)
         JI=KI-K+J
         A(KI)=-A(JI)
  250 A(JI) =HOLD
      GO TO 190
  260 RETURN
      END
