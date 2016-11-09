C   IMSL ROUTINE NAME   - VBLA=SROT
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - APPLY GIVENS PLANE ROTATION (SINGLE
C                           PRECISION)
C
C   USAGE               - CALL SROT (N,SX,INCX,SY,INCY,SC,SS)
C
C   ARGUMENTS    N      - LENGTH OF VECTORS X AND Y. (INPUT)
C                SX     - REAL VECTOR OF LENGTH MAX(N*IABS(INCX),1).
C                           (INPUT/OUTPUT)
C                           SROT REPLACES X(I) WITH SC*X(I) +
C                           SS*Y(I) FOR I=1,...,N.
C                           X(I) AND Y(I) REFER TO SPECIFIC ELEMENTS
C                           OF SX AND SY, RESPECTIVELY. SEE INCX AND
C                           INCY ARGUMENT DESCRIPTIONS.
C                INCX   - DISPLACEMENT BETWEEN ELEMENTS OF SX. (INPUT)
C                           X(I) IS DEFINED TO BE..
C                           SX(1+(I-1)*INCX) IF INCX.GE.0 OR
C                           SX(1+(I-N)*INCX) IF INCX.LT.0.
C                SY     - REAL VECTOR OF LENGTH MAX(N*IABS(INCY),1).
C                           (INPUT/OUTPUT)
C                           SROT REPLACES Y(I) WITH -SS*X(I) +
C                           SC*Y(I).
C                           X(I) AND Y(I) REFER TO SPECIFIC ELEMENTS
C                           OF SX AND SY, RESPECTIVELY. SEE INCX AND
C                           INCY ARGUMENT DESCRIPTIONS.
C                INCY   - DISPLACEMENT BETWEEN ELEMENTS OF SY. (INPUT)
C                           Y(I) IS DEFINED TO BE..
C                           SY(1+(I-1)*INCY) IF INCY.GE.0 OR
C                           SY(1+(I-N)*INCY) IF INCY.LT.0.
C                SC     - ELEMENT OF TRANSFORMATION MATRIX. (INPUT)
C                SS     - ELEMENT OF TRANSFORMATION MATRIX. (INPUT)
C
C   REQD. IMSL ROUTINES - NONE REQUIRED
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS      SROT APPLIES ( SC SS ) TO (X(1)     X(N)).
C                             (       )    (             )
C                             (-SS SC )    (Y(1) ... Y(N))
C
C-----------------------------------------------------------------------
C
      SUBROUTINE SROT   (N,SX,INCX,SY,INCY,SC,SS)
C
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N,INCX,INCY
      REAL               SX(1),SY(1),SC,SS
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,NSTEPS,KX,KY
      REAL               ONE,W,Z,ZERO
      DATA               ZERO,ONE/0.E0,1.E0/
C                                  FIRST EXECUTABLE STATEMENT
      IF (N.LE.0.OR.(SS.EQ.ZERO.AND.SC.EQ.ONE)) GO TO 20
      IF (.NOT.(INCX.EQ.INCY.AND.INCX.GT.0)) GO TO 10
C
      NSTEPS = INCX*N
      DO 5 I=1,NSTEPS,INCX
         W = SX(I)
         Z = SY(I)
         SX(I) = SC*W+SS*Z
         SY(I) = -SS*W+SC*Z
    5 CONTINUE
      GO TO 20
C
   10 CONTINUE
      KX = 1
      KY = 1
C
      IF (INCX.LT.0) KX = 1-(N-1)*INCX
      IF (INCY.LT.0) KY = 1-(N-1)*INCY
C
      DO 15 I=1,N
         W = SX(KX)
         Z = SY(KY)
         SX(KX) = SC*W+SS*Z
         SY(KY) = -SS*W+SC*Z
         KX = KX+INCX
         KY = KY+INCY
   15 CONTINUE
   20 CONTINUE
C
      RETURN
      END
