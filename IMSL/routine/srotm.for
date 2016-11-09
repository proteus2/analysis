C   IMSL ROUTINE NAME   - VBLA=SROTM
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - APPLY A MODIFIED GIVENS PLANE ROTATION
C                           (SINGLE PRECISION)
C
C   USAGE               - CALL SROTM (N,SX,INCX,SY,INCY,SPARAM)
C
C   ARGUMENTS    N      - LENGTH OF VECTORS X AND Y. (INPUT)
C                SX     - REAL VECTOR OF LENGTH MAX(N*IABS(INCX),1).
C                           (INPUT/OUTPUT)
C                           SROTM REPLACES X(I) WITH H11*X(I) +
C                           H12*Y(I) FOR I=1,...,N.
C                           X(I) AND Y(I) REFER TO SPECIFIC ELEMENTS
C                           OF SX AND SY, RESPECTIVELY. SEE INCX AND
C                           INCY ARGUMENT DESCRIPTIONS. THE H
C                           COMPONENTS REFER TO THE TRANSFORMATION
C                           DEFINED BY SPARAM.
C                INCX   - DISPLACEMENT BETWEEN ELEMENTS OF SX. (INPUT)
C                           X(I) IS DEFINED TO BE..
C                           SX(1+(I-1)*INCX) IF INCX.GE.0 OR
C                           SX(1+(I-N)*INCX) IF INCX.LT.0.
C                SY     - REAL VECTOR OF LENGTH MAX(N*IABS(INCY),1).
C                           (INPUT/OUTPUT)
C                           SROTM REPLACES Y(I) WITH H21*X(I) +
C                           H22*Y(I) FOR I=1,...,N.
C                           X(I) AND Y(I) REFER TO SPECIFIC ELEMENTS
C                           OF SX AND SY, RESPECTIVELY. SEE INCX AND
C                           INCY ARGUMENT DESCRIPTIONS. THE H
C                           COMPONENTS REFER TO THE TRANSFORMATION
C                           DEFINED BY SPARAM.
C                INCY   - DISPLACEMENT BETWEEN ELEMENTS OF SY. (INPUT)
C                           Y(I) IS DEFINED TO BE..
C                           SY(1+(I-1)*INCY) IF INCY.GE.0 OR
C                           SY(1+(I-N)*INCY) IF INCY.LT.0.
C                SPARAM - REAL VECTOR OF LENGTH 5 WHICH DEFINES THE
C                           TRANSFORMATION MATRIX H. SEE REMARKS.
C                           (INPUT)
C
C   REQD. IMSL ROUTINES - NONE REQUIRED
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS      SROTM APPLIES THE MODIFIED GIVENS TRANSFORMATION H
C                TO THE 2 BY N MATRIX     (X(1)   X(N)).
C                                         (    ...    )
C                                         (Y(1)   Y(N))
C                H TAKES ONE OF THE FOLLOWING FORMS,
C                SPARAM(1)=-2.0
C                  H11 = 1.0          H12 = 0.0
C                  H21 = 0.0          H22 = 1.0
C                SPARAM(1)=-1.0
C                  H11 = SPARAM(2)    H12 = SPARAM(4)
C                  H21 = SPARAM(3)    H22 = SPARAM(5)
C                SPARAM(1)=0.0
C                  H11 = 1.0          H12 = SPARAM(4)
C                  H21 = SPARAM(3)    H22 = 1.0
C                SPARAM(1)=1.0
C                  H11 = SPARAM(2)    H12 = 1.0
C                  H21 = -1.0         H22 = SPARAM(5)
C
C-----------------------------------------------------------------------
C
      SUBROUTINE SROTM  (N,SX,INCX,SY,INCY,SPARAM)
C
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N,INCX,INCY
      REAL               SX(1),SY(1),SPARAM(5)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,NSTEPS,KX,KY
      REAL               SFLAG,SH12,SH21,W,Z,SH11,SH2,TWO,ZERO
      DATA               ZERO,TWO/0.E0,2.E0/
C                                  FIRST EXECUTABLE STATEMENT
      SFLAG = SPARAM(1)
      IF (N.LE.0.OR.(SFLAG+TWO.EQ.ZERO)) GO TO 70
      IF (.NOT.(INCX.EQ.INCY.AND.INCX.GT.0)) GO TO 35
C
      NSTEPS = N*INCX
      IF (SFLAG) 25,5,15
    5 CONTINUE
      SH12 = SPARAM(4)
      SH21 = SPARAM(3)
      DO 10 I=1,NSTEPS,INCX
         W = SX(I)
         Z = SY(I)
         SX(I) = W+Z*SH12
         SY(I) = W*SH21+Z
   10 CONTINUE
      GO TO 70
   15 CONTINUE
      SH11 = SPARAM(2)
      SH22 = SPARAM(5)
      DO 20 I=1,NSTEPS,INCX
         W = SX(I)
         Z = SY(I)
         SX(I) = W*SH11+Z
         SY(I) = -W+SH22*Z
   20 CONTINUE
      GO TO 70
   25 CONTINUE
      SH11 = SPARAM(2)
      SH12 = SPARAM(4)
      SH21 = SPARAM(3)
      SH22 = SPARAM(5)
      DO 30 I=1,NSTEPS,INCX
         W = SX(I)
         Z = SY(I)
         SX(I) = W*SH11+Z*SH12
         SY(I) = W*SH21+Z*SH22
   30 CONTINUE
      GO TO 70
   35 CONTINUE
      KX = 1
      KY = 1
      IF (INCX.LT.0) KX = 1+(1-N)*INCX
      IF (INCY.LT.0) KY = 1+(1-N)*INCY
C
      IF (SFLAG) 60,40,50
   40 CONTINUE
      SH12 = SPARAM(4)
      SH21 = SPARAM(3)
      DO 45 I=1,N
         W = SX(KX)
         Z = SY(KY)
         SX(KX) = W+Z*SH12
         SY(KY) = W*SH21+Z
         KX = KX+INCX
         KY = KY+INCY
   45 CONTINUE
      GO TO 70
   50 CONTINUE
      SH11 = SPARAM(2)
      SH22 = SPARAM(5)
      DO 55 I=1,N
         W = SX(KX)
         Z = SY(KY)
         SX(KX) = W*SH11+Z
         SY(KY) = -W+SH22*Z
         KX = KX+INCX
         KY = KY+INCY
   55 CONTINUE
      GO TO 70
   60 CONTINUE
      SH11 = SPARAM(2)
      SH12 = SPARAM(4)
      SH21 = SPARAM(3)
      SH22 = SPARAM(5)
      DO 65 I=1,N
         W = SX(KX)
         Z = SY(KY)
         SX(KX) = W*SH11+Z*SH12
         SY(KY) = W*SH21+Z*SH22
         KX = KX+INCX
         KY = KY+INCY
   65 CONTINUE
   70 CONTINUE
      RETURN
      END
