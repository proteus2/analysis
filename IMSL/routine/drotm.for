C   IMSL ROUTINE NAME   - VBLA=DROTM
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - APPLY A MODIFIED GIVENS PLANE ROTATION
C                           (DOUBLE PRECISION)
C
C   USAGE               - CALL DROTM (N,DX,INCX,DY,INCY,DPARAM)
C
C   ARGUMENTS    N      - LENGTH OF VECTORS X AND Y. (INPUT)
C                DX     - DOUBLE PRECISION VECTOR OF LENGTH
C                           MAX(N*IABS(INCX),1). (INPUT/OUTPUT)
C                           DROTM REPLACES X(I) WITH H11*X(I) +
C                           H12*Y(I) FOR I=1,...,N.
C                           X(I) AND Y(I) REFER TO SPECIFIC ELEMENTS
C                           OF DX AND DY, RESPECTIVELY. SEE INCX AND
C                           INCY ARGUMENT DESCRIPTIONS. THE H
C                           COMPONENTS REFER TO THE TRANSFORMATION
C                           DEFINED BY DPARAM.
C                INCX   - DISPLACEMENT BETWEEN ELEMENTS OF DX. (INPUT)
C                           X(I) IS DEFINED TO BE..
C                           DX(1+(I-1)*INCX) IF INCX.GE.0 OR
C                           DX(1+(I-N)*INCX) IF INCX.LT.0.
C                DY     - DOUBLE PRECISION VECTOR OF LENGTH
C                           MAX(N*IABS(INCY),1). (INPUT/OUTPUT)
C                           DROTM REPLACES Y(I) WITH H21*X(I) +
C                           H22*Y(I) FOR I=1,...,N.
C                           X(I) AND Y(I) REFER TO SPECIFIC ELEMENTS
C                           OF DX AND DY, RESPECTIVELY. SEE INCX AND
C                           INCY ARGUMENT DESCRIPTIONS. THE H
C                           COMPONENTS REFER TO THE TRANSFORMATION
C                           DEFINED BY DPARAM.
C                INCY   - DISPLACEMENT BETWEEN ELEMENTS OF DY. (INPUT)
C                           Y(I) IS DEFINED TO BE..
C                           DY(1+(I-1)*INCY) IF INCY.GE.0 OR
C                           DY(1+(I-N)*INCY) IF INCY.LT.0.
C                DPARAM - DOUBLE PRECISION VECTOR OF LENGTH 5 WHICH
C                           DEFINES THE TRANSFORMATION MATRIX H.
C                           SEE REMARKS. (INPUT)
C
C   REQD. IMSL ROUTINES - NONE REQUIRED
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS      DROTM APPLIES THE MODIFIED GIVENS TRANSFORMATION H
C                TO THE 2 BY N MATRIX     (X(1)   X(N)).
C                                         (    ...    )
C                                         (Y(1)   Y(N))
C                H TAKES ONE OF THE FOLLOWING FORMS,
C                DPARAM(1)=-2.0
C                  H11 = 1.0          H12 = 0.0
C                  H21 = 0.0          H22 = 1.0
C                DPARAM(1)=-1.0
C                  H11 = DPARAM(2)    H12 = DPARAM(4)
C                  H21 = DPARAM(3)    H22 = DPARAM(5)
C                DPARAM(1)=0.0
C                  H11 = 1.0          H12 = DPARAM(4)
C                  H21 = DPARAM(3)    H22 = 1.0
C                DPARAM(1)=1.0
C                  H11 = DPARAM(2)    H12 = 1.0
C                  H21 = -1.0         H22 = DPARAM(5)
C
C-----------------------------------------------------------------------
C
      SUBROUTINE DROTM  (N,DX,INCX,DY,INCY,DPARAM)
C
C                                  SPECIFICATIONS FOR ARGUMENTS
      DOUBLE PRECISION   DX(1),DY(1),DPARAM(5)
      INTEGER            N,INCX,INCY
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      DOUBLE PRECISION   DFLAG,DH12,DH22,TWO,Z,DH11,DH21,W,ZERO
      INTEGER            I,NSTEPS,KX,KY
      DATA               ZERO,TWO/0.D0,2.D0/
C                                  FIRST EXECUTABLE STATEMENT
C
      DFLAG = DPARAM(1)
      IF (N.LE.0.OR.(DFLAG+TWO.EQ.ZERO)) GO TO 70
      IF (.NOT.(INCX.EQ.INCY.AND.INCX.GT.0)) GO TO 35
C
      NSTEPS = N*INCX
      IF (DFLAG) 25,5,15
    5 CONTINUE
      DH12 = DPARAM(4)
      DH21 = DPARAM(3)
      DO 10 I=1,NSTEPS,INCX
         W = DX(I)
         Z = DY(I)
         DX(I) = W+Z*DH12
         DY(I) = W*DH21+Z
   10 CONTINUE
      GO TO 70
   15 CONTINUE
      DH11 = DPARAM(2)
      DH22 = DPARAM(5)
      DO 20 I=1,NSTEPS,INCX
         W = DX(I)
         Z = DY(I)
         DX(I) = W*DH11+Z
         DY(I) = -W+DH22*Z
   20 CONTINUE
      GO TO 70
   25 CONTINUE
      DH11 = DPARAM(2)
      DH12 = DPARAM(4)
      DH21 = DPARAM(3)
      DH22 = DPARAM(5)
      DO 30 I=1,NSTEPS,INCX
         W = DX(I)
         Z = DY(I)
         DX(I) = W*DH11+Z*DH12
         DY(I) = W*DH21+Z*DH22
   30 CONTINUE
      GO TO 70
   35 CONTINUE
      KX = 1
      KY = 1
      IF (INCX.LT.0) KX = 1+(1-N)*INCX
      IF (INCY.LT.0) KY = 1+(1-N)*INCY
C
      IF (DFLAG) 60,40,50
   40 CONTINUE
      DH12 = DPARAM(4)
      DH21 = DPARAM(3)
      DO 45 I=1,N
         W = DX(KX)
         Z = DY(KY)
         DX(KX) = W+Z*DH12
         DY(KY) = W*DH21+Z
         KX = KX+INCX
         KY = KY+INCY
   45 CONTINUE
      GO TO 70
   50 CONTINUE
      DH11 = DPARAM(2)
      DH22 = DPARAM(5)
      DO 55 I=1,N
         W = DX(KX)
         Z = DY(KY)
         DX(KX) = W*DH11+Z
         DY(KY) = -W+DH22*Z
         KX = KX+INCX
         KY = KY+INCY
   55 CONTINUE
      GO TO 70
   60 CONTINUE
      DH11 = DPARAM(2)
      DH12 = DPARAM(4)
      DH21 = DPARAM(3)
      DH22 = DPARAM(5)
      DO 65 I=1,N
         W = DX(KX)
         Z = DY(KY)
         DX(KX) = W*DH11+Z*DH12
         DY(KY) = W*DH21+Z*DH22
         KX = KX+INCX
         KY = KY+INCY
   65 CONTINUE
   70 CONTINUE
      RETURN
      END
