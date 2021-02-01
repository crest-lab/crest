!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2020 Stefan Grimme
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

      SUBROUTINE PREIG(IO,OCC,F,E,istart,NORBS)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION E(*), OCC(*)

      write(io,'(/,10x,''eigenvalues'')')
      N=8 
      NTIMES=(NORBS-istart+1)/N
      NREST=MOD(NORBS-istart+1,N)
      IF(NTIMES.EQ.0) NREST=NORBS-istart+1

      J=istart
      N2=N+istart-1

      DO K=1,NTIMES
      WRITE(IO,100)(I,I=J,N2)
      WRITE(IO,200)(occ(i),I=J,N2)
      WRITE(IO,300)(F*E(i),I=J,N2)
      J =J +N
      N2=N2+N
      ENDDO

      IF(NREST.GT.0.OR.NTIMES.EQ.0) THEN
      WRITE(IO,100)(I,I=J,J+NREST-1)
      WRITE(IO,200)(occ(i),I=J,J+NREST-1)
      WRITE(IO,300)(F*E(i),I=J,J+NREST-1)
      ENDIF

      RETURN

 100  FORMAT(' #    : ',2X,8(3X,I6,2X))
 200  FORMAT(' occ. : ',2X,8(4X,F6.3,1X))
 300  FORMAT(' eps  : ',2X,8F11.3)
      END

      SUBROUTINE PREIG2(IO,OCC,E,NORBS)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION E(*), OCC(*)

      N=6
      NTIMES=NORBS/N
      NREST=MOD(NORBS,N)
      IF(NTIMES.EQ.0) NREST=NORBS

      J=1
      N2=N

      DO K=1,NTIMES
      WRITE(IO,100)(I,I=J,N2)
      WRITE(IO,200)(occ(i),I=J,N2)
      WRITE(IO,300)(E(i),I=J,N2)
      J =J +N
      N2=N2+N
      ENDDO

      IF(NREST.GT.0.OR.NTIMES.EQ.0) THEN
      WRITE(IO,100)(I,I=J,J+NREST-1)
      WRITE(IO,200)(occ(i),I=J,J+NREST-1)
      WRITE(IO,300)(E(i),I=J,J+NREST-1)
      ENDIF

      RETURN

 100  FORMAT('#       :',2X,6(4X,I4,2X))
 200  FORMAT('# atoms :',2X,6(4X,F5.3,1X))
 300  FORMAT('shift ev:',2X,6F10.5)
      END

      SUBROUTINE PREIG3(IO,E,NORBS)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION E(*)
      N=10
      NTIMES=NORBS/N
      NREST=MOD(NORBS,N)
      IF(NTIMES.EQ.0) NREST=NORBS
      J=1
      N2=N
      DO K=1,NTIMES
c     WRITE(IO,100)(I,I=J,N2)
      WRITE(IO,300)J,N2,(E(i),I=J,N2)
      J =J +N
      N2=N2+N
      ENDDO
      IF(NREST.GT.0.OR.NTIMES.EQ.0) THEN
      WRITE(IO,300)J,J+NREST-1,(E(i),I=J,J+NREST-1)
      ENDIF

 100  FORMAT('atoms : ',2X,6(4X,I4,2X))
 300  FORMAT(' value',i5,'-',i5,':',2X,12F6.2)
      END

      SUBROUTINE PRMAT(IUOUT,R,N,M,HEAD)
      REAL*8 R
      CHARACTER*(*) HEAD
      DIMENSION R(*)
C     SUBROUTINE PRINTS MATRIX R,WHICH IS SUPPOSED
C     TO HAVE DIMENSION N,M  WHEN M IS NONZERO AND
C     ((N+1)*N)/2 WHEN M IS ZERO

      WRITE(IUOUT,1001) HEAD
      NKPB=6
      IF(M)10,10,80
C
   10 CONTINUE
      IBL=N/NKPB
      IR=N-IBL*NKPB
      J1=1
      K1S=1
      KD=0
      IF(IBL.EQ.0) GO TO 50
      J2=NKPB
      DO 40 I=1,IBL
      WRITE(IUOUT,1002)(J,J=J1,J2)
      K1=K1S
      K2=K1
      KK=0
      DO 20 J=J1,J2
      WRITE(IUOUT,1003)J,(R(K),K=K1,K2)
      KK=KK+1
      K1=K1+KD+KK
   20 K2=K1+KK
      J1=J1+NKPB
      IF(J1.GT.N) RETURN
      J2=J2+NKPB
      K2=K1-1
      K1=K2+1
      K2=K1+(NKPB-1)
      K1S=K2+1
      KK=KD+NKPB
      DO 30 J=J1,N
      WRITE(IUOUT,1003)J,(R(K),K=K1,K2)
      KK=KK+1
      K1=K1+KK
   30 K2=K2+KK
   40 KD=KD+NKPB
   50 IF(IR.EQ.0) GO TO 70
      K1=K1S
      J2=J1+IR-1
      KK=0
      K2=K1
      WRITE(IUOUT,1002)(J,J=J1,J2)
      WRITE(IUOUT,1003)
      DO 60 J=J1,J2
      WRITE(IUOUT,1003)J,(R(K),K=K1,K2)
      KK=KK+1
      K1=K1+KD+KK
   60 K2=K1+KK
   70 RETURN
   80 IBL=M/NKPB
      IR=M-IBL*NKPB
      I2=0
      K2=0
      IF(IBL.EQ.0) GO TO 100
      DO 90 I=1,IBL
      I1=(I-1)*N*NKPB+1
      I2=I1+(NKPB-1)*N
      K1=K2+1
      K2=K1+(NKPB-1)
      WRITE(IUOUT,1002)(K,K=K1,K2)
      DO 90 J=1,N
      WRITE(IUOUT,1003)J,(R(IJ),IJ=I1,I2,N)
      I1=I1+1
   90 I2=I1+(NKPB-1)*N
  100 IF(IR.EQ.0) GO TO 120
      I1=IBL*N*NKPB+1
      I2=I1+(IR-1)*N
      K1=K2+1
      K2=M
      WRITE(IUOUT,1002)(K,K=K1,K2)
      WRITE(IUOUT,1003)
      DO 110 J=1,N
      WRITE(IUOUT,1003)J,(R(IJ),IJ=I1,I2,N)
      I1=I1+1
      I2=I1+(IR-1)*N
  110 CONTINUE
  120 WRITE(IUOUT,1003)
      RETURN
 1001 FORMAT(/' MATRIX PRINTED:',2X,A)
 1002 FORMAT(/,' ',4X,6(3X,I4,3X),/)
 1003 FORMAT(' ',I4,6F10.5)
      END

      SUBROUTINE PRMAT4(IUOUT,R,N,M,HEAD)
      REAL*4 R
      CHARACTER*(*) HEAD
      DIMENSION R(*)
C     SUBROUTINE PRINTS MATRIX R,WHICH IS SUPPOSED
C     TO HAVE DIMENSION N,M  WHEN M IS NONZERO AND
C     ((N+1)*N)/2 WHEN M IS ZERO

      WRITE(IUOUT,1001) HEAD
      NKPB=6
      IF(M)10,10,80
C
   10 CONTINUE
      IBL=N/NKPB
      IR=N-IBL*NKPB
      J1=1
      K1S=1
      KD=0
      IF(IBL.EQ.0) GO TO 50
      J2=NKPB
      DO 40 I=1,IBL
      WRITE(IUOUT,1002)(J,J=J1,J2)
      K1=K1S
      K2=K1
      KK=0
      DO 20 J=J1,J2
      WRITE(IUOUT,1003)J,(R(K),K=K1,K2)
      KK=KK+1
      K1=K1+KD+KK
   20 K2=K1+KK
      J1=J1+NKPB
      IF(J1.GT.N) RETURN
      J2=J2+NKPB
      K2=K1-1
      K1=K2+1
      K2=K1+(NKPB-1)
      K1S=K2+1
      KK=KD+NKPB
      DO 30 J=J1,N
      WRITE(IUOUT,1003)J,(R(K),K=K1,K2)
      KK=KK+1
      K1=K1+KK
   30 K2=K2+KK
   40 KD=KD+NKPB
   50 IF(IR.EQ.0) GO TO 70
      K1=K1S
      J2=J1+IR-1
      KK=0
      K2=K1
      WRITE(IUOUT,1002)(J,J=J1,J2)
      WRITE(IUOUT,1003)
      DO 60 J=J1,J2
      WRITE(IUOUT,1003)J,(R(K),K=K1,K2)
      KK=KK+1
      K1=K1+KD+KK
   60 K2=K1+KK
   70 RETURN
   80 IBL=M/NKPB
      IR=M-IBL*NKPB
      I2=0
      K2=0
      IF(IBL.EQ.0) GO TO 100
      DO 90 I=1,IBL
      I1=(I-1)*N*NKPB+1
      I2=I1+(NKPB-1)*N
      K1=K2+1
      K2=K1+(NKPB-1)
      WRITE(IUOUT,1002)(K,K=K1,K2)
      DO 90 J=1,N
      WRITE(IUOUT,1003)J,(R(IJ),IJ=I1,I2,N)
      I1=I1+1
   90 I2=I1+(NKPB-1)*N
  100 IF(IR.EQ.0) GO TO 120
      I1=IBL*N*NKPB+1
      I2=I1+(IR-1)*N
      K1=K2+1
      K2=M
      WRITE(IUOUT,1002)(K,K=K1,K2)
      WRITE(IUOUT,1003)
      DO 110 J=1,N
      WRITE(IUOUT,1003)J,(R(IJ),IJ=I1,I2,N)
      I1=I1+1
      I2=I1+(IR-1)*N
  110 CONTINUE
  120 WRITE(IUOUT,1003)
      RETURN
 1001 FORMAT(/' MATRIX PRINTED:',2X,A)
 1002 FORMAT(/,' ',4X,6(3X,I4,3X),/)
 1003 FORMAT(' ',I4,6F10.5)
      END

      SUBROUTINE PRMATI(IUOUT,RR,N,M,HEAD)
      CHARACTER*(*) HEAD
      integer RR(*)
      DIMENSION R(n*n)
C     SUBROUTINE PRINTS MATRIX R,WHICH IS SUPPOSED
C     TO HAVE DIMENSION N,M  WHEN M IS NONZERO AND
C     ((N+1)*N)/2 WHEN M IS ZERO

      R(1:n*n)=float(RR(1:n*n))
      WRITE(IUOUT,1001) HEAD
      NKPB=6
      IF(M)10,10,80
C
   10 CONTINUE
      IBL=N/NKPB
      IR=N-IBL*NKPB
      J1=1
      K1S=1
      KD=0
      IF(IBL.EQ.0) GO TO 50
      J2=NKPB
      DO 40 I=1,IBL
      WRITE(IUOUT,1002)(J,J=J1,J2)
      K1=K1S
      K2=K1
      KK=0
      DO 20 J=J1,J2
      WRITE(IUOUT,1003)J,(R(K),K=K1,K2)
      KK=KK+1
      K1=K1+KD+KK
   20 K2=K1+KK
      J1=J1+NKPB
      IF(J1.GT.N) RETURN
      J2=J2+NKPB
      K2=K1-1
      K1=K2+1
      K2=K1+(NKPB-1)
      K1S=K2+1
      KK=KD+NKPB
      DO 30 J=J1,N
      WRITE(IUOUT,1003)J,(R(K),K=K1,K2)
      KK=KK+1
      K1=K1+KK
   30 K2=K2+KK
   40 KD=KD+NKPB
   50 IF(IR.EQ.0) GO TO 70
      K1=K1S
      J2=J1+IR-1
      KK=0
      K2=K1
      WRITE(IUOUT,1002)(J,J=J1,J2)
      WRITE(IUOUT,1003)
      DO 60 J=J1,J2
      WRITE(IUOUT,1003)J,(R(K),K=K1,K2)
      KK=KK+1
      K1=K1+KD+KK
   60 K2=K1+KK
   70 RETURN
   80 IBL=M/NKPB
      IR=M-IBL*NKPB
      I2=0
      K2=0
      IF(IBL.EQ.0) GO TO 100
      DO 90 I=1,IBL
      I1=(I-1)*N*NKPB+1
      I2=I1+(NKPB-1)*N
      K1=K2+1
      K2=K1+(NKPB-1)
      WRITE(IUOUT,1002)(K,K=K1,K2)
      DO 90 J=1,N
      WRITE(IUOUT,1003)J,(R(IJ),IJ=I1,I2,N)
      I1=I1+1
   90 I2=I1+(NKPB-1)*N
  100 IF(IR.EQ.0) GO TO 120
      I1=IBL*N*NKPB+1
      I2=I1+(IR-1)*N
      K1=K2+1
      K2=M
      WRITE(IUOUT,1002)(K,K=K1,K2)
      WRITE(IUOUT,1003)
      DO 110 J=1,N
      WRITE(IUOUT,1003)J,(R(IJ),IJ=I1,I2,N)
      I1=I1+1
      I2=I1+(IR-1)*N
  110 CONTINUE
  120 WRITE(IUOUT,1003)
      RETURN
 1001 FORMAT(/' MATRIX PRINTED:',2X,A)
 1002 FORMAT(/,' ',4X,6(3X,I4,3X),/)
 1003 FORMAT(' ',I4,6F10.5)
      END

      SUBROUTINE PRMATS(IUOUT,R,N,M,HEAD)
      REAL*8 R
      CHARACTER*(*) HEAD
      DIMENSION R(*)
C     SUBROUTINE PRINTS MATRIX R,WHICH IS SUPPOSED
C     TO HAVE DIMENSION N,M  WHEN M IS NONZERO AND
C     ((N+1)*N)/2 WHEN M IS ZERO

      WRITE(IUOUT,1001) HEAD
      NKPB=12
      IF(M)10,10,80
C
   10 CONTINUE
      IBL=N/NKPB
      IR=N-IBL*NKPB
      J1=1
      K1S=1
      KD=0
      IF(IBL.EQ.0) GO TO 50
      J2=NKPB
      DO 40 I=1,IBL
      WRITE(IUOUT,1002)(J,J=J1,J2)
      K1=K1S
      K2=K1
      KK=0
      DO 20 J=J1,J2
      WRITE(IUOUT,1003)J,(R(K),K=K1,K2)
      KK=KK+1
      K1=K1+KD+KK
   20 K2=K1+KK
      J1=J1+NKPB
      IF(J1.GT.N) RETURN
      J2=J2+NKPB
      K2=K1-1
      K1=K2+1
      K2=K1+(NKPB-1)
      K1S=K2+1
      KK=KD+NKPB
      DO 30 J=J1,N
      WRITE(IUOUT,1003)J,(R(K),K=K1,K2)
      KK=KK+1
      K1=K1+KK
   30 K2=K2+KK
   40 KD=KD+NKPB
   50 IF(IR.EQ.0) GO TO 70
      K1=K1S
      J2=J1+IR-1
      KK=0
      K2=K1
      WRITE(IUOUT,1002)(J,J=J1,J2)
      WRITE(IUOUT,1003)
      DO 60 J=J1,J2
      WRITE(IUOUT,1003)J,(R(K),K=K1,K2)
      KK=KK+1
      K1=K1+KD+KK
   60 K2=K1+KK
   70 RETURN
   80 IBL=M/NKPB
      IR=M-IBL*NKPB
      I2=0
      K2=0
      IF(IBL.EQ.0) GO TO 100
      DO 90 I=1,IBL
      I1=(I-1)*N*NKPB+1
      I2=I1+(NKPB-1)*N
      K1=K2+1
      K2=K1+(NKPB-1)
      WRITE(IUOUT,1002)(K,K=K1,K2)
      DO 90 J=1,N
      WRITE(IUOUT,1003)J,(R(IJ),IJ=I1,I2,N)
      I1=I1+1
   90 I2=I1+(NKPB-1)*N
  100 IF(IR.EQ.0) GO TO 120
      I1=IBL*N*NKPB+1
      I2=I1+(IR-1)*N
      K1=K2+1
      K2=M
      WRITE(IUOUT,1002)(K,K=K1,K2)
      WRITE(IUOUT,1003)
      DO 110 J=1,N
      WRITE(IUOUT,1003)J,(R(IJ),IJ=I1,I2,N)
      I1=I1+1
      I2=I1+(IR-1)*N
  110 CONTINUE
  120 WRITE(IUOUT,1003)
      RETURN
 1001 FORMAT(/' MATRIX PRINTED:',2X,A)
 1002 FORMAT(/,' ',2X,12(2X,I3,1X),/)
 1003 FORMAT(' ',I4,12F6.2)
      END

      SUBROUTINE PREIGF(IO,E,NORBS)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION E(*)
      N=6
      NTIMES=NORBS/N
      NREST=MOD(NORBS,N)
      IF(NTIMES.EQ.0) NREST=NORBS
      J=1
      N2=N
      DO K=1,NTIMES
c     WRITE(IO,100)(I,I=J,N2)
c     WRITE(IO,300)J,N2,(E(i),I=J,N2)
      WRITE(IO,300)     (E(i),I=J,N2)
      J =J +N
      N2=N2+N
      ENDDO
      IF(NREST.GT.0.OR.NTIMES.EQ.0) THEN
c     WRITE(IO,100)(I,I=J,J+NREST-1)
c     WRITE(IO,300)J,J+NREST-1,(E(i),I=J,J+NREST-1)
      WRITE(IO,300)            (E(i),I=J,J+NREST-1)
      ENDIF

 100  FORMAT(' #     : ',2X,6(4X,I4,2X))
!300  FORMAT('eigval',i4,'-',i4,':',2X,10F9.2)
 300  FORMAT('eigval : '           ,2X,10F9.2)
      END

      SUBROUTINE PREIGF0(IO,E,NORBS)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION E(*)
      N=6
      NTIMES=NORBS/N
      NREST=MOD(NORBS,N)
      IF(NTIMES.EQ.0) NREST=NORBS
      J=1
      N2=N
      DO K=1,NTIMES
      WRITE(IO,300)     (E(i),I=J,N2)
      J =J +N
      N2=N2+N
      ENDDO
      IF(NREST.GT.0.OR.NTIMES.EQ.0) THEN
      WRITE(IO,300)            (E(i),I=J,J+NREST-1)
      ENDIF

 100  FORMAT(' #     : ',2X,6(4X,I4,2X))
 300  FORMAT('eig.   : '           ,2X,10F9.2)
      END
