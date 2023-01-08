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

      SUBROUTINE XYZGEO(XYZ,NUMAT,NA,NB,NC,DEGREE,GEO)
        IMPLICIT DOUBLE PRECISION(A-H,O-Z)
        DIMENSION XYZ(3,*),NA(*),NB(*),NC(*),GEO(3,*)
!***********************************************************************
!*
!*   XYZGEO CONVERTS COORDINATES FROM CARTESIAN TO INTERNAL.
!*
!*     ON INPUT XYZ  = ARRAY OF CARTESIAN COORDINATES
!*              NUMAT= NUMBER OF ATOMS
!*              NA   = NUMBERS OF ATOM TO WHICH ATOMS ARE RELATED
!*                     BY DISTANCE
!*              NB   = NUMBERS OF ATOM TO WHICH ATOMS ARE RELATED
!*                     BY ANGLE
!*              NC   = NUMBERS OF ATOM TO WHICH ATOMS ARE RELATED
!*                     BY DIHEDRAL
!*
!*    ON OUTPUT GEO  = INTERNAL COORDINATES IN ANGSTROMS, RADIANS,
!*                     AND RADIANS
!*
!***********************************************************************
        DO I = 2,NUMAT
          J = NA(I)
          K = NB(I)
          L = NC(I)
          GEO(1,I) = SQRT((XYZ(1,I)-XYZ(1,J))**2+ &
      &                  (XYZ(2,I)-XYZ(2,J))**2+ &
      &                  (XYZ(3,I)-XYZ(3,J))**2)

          IF (I .LT. 3) CYCLE
          II = I
          CALL BANGLE(XYZ,II,J,K,GEO(2,I))
          GEO(2,I) = GEO(2,I)*DEGREE
          IF (I .LT. 4) cycle
          CALL DIHED(XYZ,II,J,K,L,GEO(3,I))
          GEO(3,I) = GEO(3,I)*DEGREE
        END DO
        GEO(1,1) = 0.D0
        GEO(2,1) = 0.D0
        GEO(3,1) = 0.D0
        GEO(2,2) = 0.D0
        GEO(3,2) = 0.D0
        GEO(3,3) = 0.D0
        RETURN
      END SUBROUTINE XYZGEO

!C     *****************************************************************
      SUBROUTINE XYZINT(XYZ,NUMAT,NA,NB,NC,DEGREE,GEO)
        IMPLICIT DOUBLE PRECISION(A-H,O-Z)
        DIMENSION XYZ(3,*),NA(*),NB(*),NC(*),GEO(3,*)
!***********************************************************************
!*
!* XYZINT WORKS OUT THE INTERNAL COORDINATES OF A MOLECULE.
!*        THE "RULES" FOR THE CONNECTIVITY ARE AS FOLLOWS:
!*        ATOM I IS DEFINED AS BEING AT A DISTANCE FROM THE NEAREST
!*        ATOM J, ATOM J ALREADY HAVING BEEN DEFINED.
!*        ATOM I MAKES AN ANGLE WITH ATOM J AND THE ATOM K, WHICH HAS
!*        ALREADY BEEN DEFINED, AND IS THE NEAREST ATOM TO J
!*        ATOM I MAKES A DIHEDRAL ANGLE WITH ATOMS J, K, AND L. L HAVING
!*        BEEN DEFINED AND IS THE NEAREST ATOM TO K, AND J, K AND L
!*        HAVE A CONTAINED ANGLE IN THE RANGE 15 TO 165 DEGREES,
!*        IF POSSIBLE.
!*
!*        NOTE THAT GEO AND XYZ MUST NOT BE THE SAME IN THE CALL.
!*
!*   ON INPUT XYZ    = CARTESIAN ARRAY OF NUMAT ATOMS
!*            DEGREE = 1 IF ANGLES ARE TO BE IN RADIANS
!*            DEGREE = 57.29578 IF ANGLES ARE TO BE IN RADIANS
!*
!***********************************************************************
        NAI1 = 0
        NAI2 = 0
        K = 1
        DO I = 1,NUMAT
          NA(I) = 2
          NB(I) = 3
          NC(I) = 4
          IM1 = I-1
          IF (IM1 .EQ. 0) CYCLE
          SUM = 100.D0
          DO J = 1,IM1
            R = (XYZ(1,I)-XYZ(1,J))**2+ &
      &       (XYZ(2,I)-XYZ(2,J))**2+ &
      &       (XYZ(3,I)-XYZ(3,J))**2
            IF (R .LT. SUM.AND.NA(J) .NE. J.AND.NB(J) .NE. J) THEN
              SUM = R
              K = J
            END IF
          END DO
!C
!C   ATOM I IS NEAREST TO ATOM K
!C
          NA(I) = K
          IF (I .GT. 2) NB(I) = NA(K)
          IF (I .GT. 3) NC(I) = NB(K)
!C
!C   FIND ANY ATOM TO RELATE TO NA(I)
!C
        END DO
        NA(1) = 0
        NB(1) = 0
        NC(1) = 0
        NB(2) = 0
        NC(2) = 0
        NC(3) = 0

        CALL XYZGEO(XYZ,NUMAT,NA,NB,NC,DEGREE,GEO)
        RETURN
      END SUBROUTINE XYZINT

      SUBROUTINE BANGLE(XYZ,I,J,K,ANGLE)
        IMPLICIT DOUBLE PRECISION(A-H,O-Z)
        DIMENSION XYZ(3,*)
!*********************************************************************
!*
!* BANGLE CALCULATES THE ANGLE BETWEEN ATOMS I,J, AND K. THE
!*        CARTESIAN COORDINATES ARE IN XYZ.
!*
!*********************************************************************
        D2IJ = (XYZ(1,I)-XYZ(1,J))**2+ &
       &       (XYZ(2,I)-XYZ(2,J))**2+ &
       &       (XYZ(3,I)-XYZ(3,J))**2
        D2JK = (XYZ(1,J)-XYZ(1,K))**2+ &
       &       (XYZ(2,J)-XYZ(2,K))**2+ &
       &       (XYZ(3,J)-XYZ(3,K))**2
        D2IK = (XYZ(1,I)-XYZ(1,K))**2+ &
       &       (XYZ(2,I)-XYZ(2,K))**2+ &
       &       (XYZ(3,I)-XYZ(3,K))**2
        XY = SQRT(D2IJ*D2JK+1.d-14)
        TEMP = 0.5D0*(D2IJ+D2JK-D2IK)/XY
        IF (TEMP .GT. 1.0D0) TEMP = 1.0D0
        IF (TEMP .LT. -1.0D0) TEMP = -1.0D0
        ANGLE = ACOS(TEMP)
        RETURN
      END SUBROUTINE BANGLE
      SUBROUTINE DIHED(XYZ,I,J,K,L,ANGLE)
        IMPLICIT DOUBLE PRECISION(A-H,O-Z)
        DIMENSION XYZ(3,*)
!*********************************************************************
!*
!*      DIHED CALCULATES THE DIHEDRAL ANGLE BETWEEN ATOMS I, J, K,
!*            AND L.  THE CARTESIAN COORDINATES OF THESE ATOMS
!*            ARE IN ARRAY XYZ.
!*
!*     DIHED IS A MODIFIED VERSION OF A SUBROUTINE OF THE SAME NAME
!*           WHICH WAS WRITTEN BY DR. W. THEIL IN 1973.
!*
!*********************************************************************
        DATA PI/3.14159265358979D0/
        XI1 = XYZ(1,I)-XYZ(1,K)
        XJ1 = XYZ(1,J)-XYZ(1,K)
        XL1 = XYZ(1,L)-XYZ(1,K)
        YI1 = XYZ(2,I)-XYZ(2,K)
        YJ1 = XYZ(2,J)-XYZ(2,K)
        YL1 = XYZ(2,L)-XYZ(2,K)
        ZI1 = XYZ(3,I)-XYZ(3,K)
        ZJ1 = XYZ(3,J)-XYZ(3,K)
        ZL1 = XYZ(3,L)-XYZ(3,K)
!C      ROTATE AROUND Z AXIS TO PUT KJ ALONG Y AXIS
        DIST = SQRT(XJ1**2+YJ1**2+ZJ1**2+1.d-14)
        COSA = ZJ1/DIST
        IF (COSA .GT. 1.0D0) COSA = 1.0D0
        IF (COSA .LT. -1.0D0) COSA = -1.0D0
        DDD = 1.0D0-COSA**2
        IF (DDD .LE. 0.0) GO TO 10
        YXDIST = DIST*SQRT(DDD)
        IF (YXDIST .GT. 1.0D-12) GO TO 20
10      CONTINUE
        XI2 = XI1
        XL2 = XL1
        YI2 = YI1
        YL2 = YL1
        COSTH = COSA
        SINTH = 0.D0
        GO TO 30
20      COSPH = YJ1/YXDIST
        SINPH = XJ1/YXDIST
        XI2 = XI1*COSPH-YI1*SINPH
        XJ2 = XJ1*COSPH-YJ1*SINPH
        XL2 = XL1*COSPH-YL1*SINPH
        YI2 = XI1*SINPH+YI1*COSPH
        YJ2 = XJ1*SINPH+YJ1*COSPH
        YL2 = XL1*SINPH+YL1*COSPH
!C      ROTATE KJ AROUND THE X AXIS SO KJ LIES ALONG THE Z AXIS
        COSTH = COSA
        SINTH = YJ2/DIST
30      CONTINUE
        YI3 = YI2*COSTH-ZI1*SINTH
        YL3 = YL2*COSTH-ZL1*SINTH
        CALL DANG(XL2,YL3,XI2,YI3,ANGLE)
!c     IF (ANGLE .LT. 0.) ANGLE=2.0D0*PI+ANGLE
!csg   IF (ANGLE .GE. 2.0d0*PI    ) ANGLE=0.D0
        RETURN
      END SUBROUTINE DIHED
      SUBROUTINE DANG(A1,A2,B1,B2,RCOS)
        IMPLICIT DOUBLE PRECISION(A-H,O-Z)
!**********************************************************************
!*
!*    DANG  DETERMINES THE ANGLE BETWEEN THE POINTS (A1,A2), (0,0),
!*          AND (B1,B2).  THE RESULT IS PUT IN RCOS.
!*
!**********************************************************************
        DATA PI/3.14159265358979D0/
        ZERO = 1.0D-10
        IF (ABS(A1) .LT. ZERO.AND.ABS(A2) .LT. ZERO) GO TO 10
        IF (ABS(B1) .LT. ZERO.AND.ABS(B2) .LT. ZERO) GO TO 10
        ANORM = 1.0D0/SQRT(A1**2+A2**2)
        BNORM = 1.0D0/SQRT(B1**2+B2**2)
        A1 = A1*ANORM
        A2 = A2*ANORM
        B1 = B1*BNORM
        B2 = B2*BNORM
        SINTH = (A1*B2)-(A2*B1)
        COSTH = A1*B1+A2*B2
        IF (COSTH .GT. 1.0D0) COSTH = 1.0D0
        IF (COSTH .LT. -1.0D0) COSTH = -1.0D0
        RCOS = ACOS(COSTH)
        IF (ABS(RCOS) .LT. ZERO) GO TO 10
        IF (SINTH .GT. 0.D0) RCOS = 2.0D0*PI-RCOS
        RCOS = -RCOS
        RETURN
10      RCOS = 0.0D0
        RETURN
      END SUBROUTINE DANG

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE GMETRY(NATOMS,geo,coord,na,nb,nc,fail)
        IMPLICIT REAL*8(A-H,O-Z)

        DIMENSION GEO(3,*),coord(3,*),na(*),nb(*),nc(*)
        logical fail

        fail = .false.

        COORD(1,1) = 0.0D00
        COORD(2,1) = 0.0D00
        COORD(3,1) = 0.0D00
        COORD(1,2) = GEO(1,2)
        COORD(2,2) = 0.0D00
        COORD(3,2) = 0.0D00
        IF (NATOMS .EQ. 2) GOTO 110
        CCOS = COS(GEO(2,3))
        IF (NA(3) .EQ. 1) THEN
          COORD(1,3) = COORD(1,1)+GEO(1,3)*CCOS
        ELSE
          COORD(1,3) = COORD(1,2)-GEO(1,3)*CCOS
        END IF
        COORD(2,3) = GEO(1,3)*SIN(GEO(2,3))
        COORD(3,3) = 0.0D00
        DO I = 4,NATOMS
          COSA = COS(GEO(2,I))
          MB = NB(I)
          MC = NA(I)
          XB = COORD(1,MB)-COORD(1,MC)
          YB = COORD(2,MB)-COORD(2,MC)
          ZB = COORD(3,MB)-COORD(3,MC)
          RBC = 1.0D00/DSQRT(XB*XB+YB*YB+ZB*ZB)
          !IF (ABS(COSA).LT.0.9999999999D0) GO TO 40
          IF (ABS(COSA) .GE. 0.9999999999D0) THEN
!C
!C     ATOMS MC, MB, AND (I) ARE COLLINEAR
!C
            RBC = GEO(1,I)*RBC*COSA
            COORD(1,I) = COORD(1,MC)+XB*RBC
            COORD(2,I) = COORD(2,MC)+YB*RBC
            COORD(3,I) = COORD(3,MC)+ZB*RBC
            !GO TO 100
            CYCLE
          END IF
!C
!C     THE ATOMS ARE NOT COLLINEAR
!C
40        MA = NC(I)
          XA = COORD(1,MA)-COORD(1,MC)
          YA = COORD(2,MA)-COORD(2,MC)
          ZA = COORD(3,MA)-COORD(3,MC)
!C
!C     ROTATE ABOUT THE Z-AXIS TO MAKE YB=0, AND XB POSITIVE.  IF XYB IS
!C     TOO SMALL, FIRST ROTATE THE Y-AXIS BY 90 DEGREES.
!C
          XYB = DSQRT(XB*XB+YB*YB)
          K = -1
          !IF (XYB.GT.0.1D00) GO TO 50
          IF (XYB .LE. 0.1D00) THEN
            XPA = ZA
            ZA = -XA
            XA = XPA
            XPB = ZB
            ZB = -XB
            XB = XPB
            XYB = DSQRT(XB*XB+YB*YB)
            K = 1
          END IF
!C
!C     ROTATE ABOUT THE Y-AXIS TO MAKE ZB VANISH
!C
50        COSTH = XB/XYB
          SINTH = YB/XYB
          XPA = XA*COSTH+YA*SINTH
          YPA = YA*COSTH-XA*SINTH
          SINPH = ZB*RBC
          COSPH = DSQRT(ABS(1.D00-SINPH*SINPH))
          XQA = XPA*COSPH+ZA*SINPH
          ZQA = ZA*COSPH-XPA*SINPH
!C
!C     ROTATE ABOUT THE X-AXIS TO MAKE ZA=0, AND YA POSITIVE.
!C
          YZA = DSQRT(YPA*YPA+ZQA*ZQA)
          IF (YZA .LT. 1.D-4) THEN
!            IF(YZA.LT.1.D-4)GOTO 70
!c           WRITE(*,'(/9X,'' ATOMS'',I3,'','',I3,'', AND'',I3,
!c    1'' ARE WITHIN'',F7.4,'' ANGSTROMS OF A STRAIGHT LINE'')')
!c    2MC,MB,MA,YZA
!            fail=.true.
!            return
            COSKH = 1.D0
            SINKH = 0.D0
          ELSE
            COSKH = YPA/YZA
            SINKH = ZQA/YZA
          END IF
!C
!C     COORDINATES :-   A=(XQA,YZA,0),   B=(RBC,0,0),  C=(0,0,0)
!C     NONE ARE NEGATIVE.
!C     THE COORDINATES OF I ARE EVALUATED IN THE NEW FRAME.
!C
          SINA = SIN(GEO(2,I))
          SIND = -SIN(GEO(3,I))
          COSD = COS(GEO(3,I))
          XD = GEO(1,I)*COSA
          YD = GEO(1,I)*SINA*COSD
          ZD = GEO(1,I)*SINA*SIND
!C
!C     TRANSFORM THE COORDINATES BACK TO THE ORIGINAL SYSTEM.
!C
          YPD = YD*COSKH-ZD*SINKH
          ZPD = ZD*COSKH+YD*SINKH
          XPD = XD*COSPH-ZPD*SINPH
          ZQD = ZPD*COSPH+XD*SINPH
          XQD = XPD*COSTH-YPD*SINTH
          YQD = YPD*COSTH+XPD*SINTH
          !IF (K.LT.1) GO TO 90
          IF (K .GE. 1) THEN
            XRD = -ZQD
            ZQD = XQD
            XQD = XRD
          END IF
90        COORD(1,I) = XQD+COORD(1,MC)
          COORD(2,I) = YQD+COORD(2,MC)
          COORD(3,I) = ZQD+COORD(3,MC)
        END DO
100     CONTINUE
110     CONTINUE

        RETURN
      END SUBROUTINE GMETRY

