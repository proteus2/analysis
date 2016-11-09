C   IMSL ROUTINE NAME   - VBLA=DROT
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - APPLY GIVENS PLANE ROTATION (DOUBLE
C                           PRECISION)
C
C   USAGE               - CALL DROT (N,DX,INCX,DY,INCY,DC,DS)
C
C   ARGUMENTS    N      - LENGTH OF VECTORS X AND Y. (INPUT)
C                DX     - DOUBLE PRECISION VECTOR OF LENGTH
C                           MAX(N*IABS(INCX),1). (INPUT/OUTPUT)
C                           DROT REPLACES X(I) WITH DC*X(I)+DS*Y(I)
C                           FOR I=1,...,N.
C                           X(I) AND Y(I) REFER TO SPECIFIC ELEMENTS
C                           OF DX AND DY, RESPECTIVELY. SEE INCX AND
C                           INCY ARGUMENT DESCRIPTIONS.
C                INCX   - DISPLACEMENT BETWEEN ELEMENTS OF DX. (INPUT)
C                           X(I) IS DEFINED TO BE..
C                           DX(1+(I-1)*INCX) IF INCX.GE.0 OR
C                           DX(1+(I-N)*INCX) IF INCX.LT.0.
C                DY     - DOUBLE PRECISION VECTOR OF LENGTH
C                           MAX(N*IABS(INCY),1). (INPUT/OUTPUT)
C                           DROT REPLACES Y(I) WITH -DS*X(I)+DC*Y(I)
C                           FOR I=1,...,N.
C                           X(I) AND Y(I) REFER TO SPECIFIC ELEMENTS
C                           OF DX AND DY, RESPECTIVELY. SEE INCX AND
C                           INCY ARGUMENT DESCRIPTIONS.
C                INCY   - DISPLACEMENT BETWEEN ELEMENTS OF DY. (INPUT)
C                           Y(I) IS DEFINED TO BE..
C                           DY(1+(I-1)*INCY) IF INCY.GE.0 OR
C                           DY(1+(I-N)*INCY) IF INCY.LT.0.
C                DC     - DOUBLE PRECISION ELEMENT OF TRANSFORMATION
C                           MATRIX. (INPUT)
C                DS     - DOUBLE PRECISION ELEMENT OF TRANSFORMATION
C                           MATRIX. (INPUT)
C
C   REQD. IMSL ROUTINES - NONE REQUIRED
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS      DROT APPLIES ( DC DS ) TO (X(1)     X(N)).
C                             (       )    (             )
C                             (-DS DC )    (Y(1) ... Y(N))
C
C-----------------------------------------------------------------------
C
      SUBROUTINE DROT   (N,DX,INCX,DY,INCY,DC,DS)
C
C                                  SPECIFICATIONS FOR ARGUMENTS
      DOUBLE PRECISION   DX(1),DY(1),DC,DS
      INTEGER            N,INCX,INCY
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      DOUBLE PRECISION   ZERO,ONE,W,Z
      INTEGER            I,NSTEPS,KX,KY
      DATA               ZERO,ONE/0.D0,1.D0/
C                                  FIRST EXECUTABLE STATEMENT
      IF (N.LE.0.OR.(DS.EQ.ZERO.AND.DC.EQ.ONE)) GO TO 20
      IF (.NOT.(INCX.EQ.INCY.AND.INCX.GT.0)) GO TO 10
C
      NSTEPS = INCX*N
      DO 5 I=1,NSTEPS,INCX
         W = DX(I)
         Z = DY(I)
         DX(I) = DC*W+DS*Z
         DY(I) = -DS*W+DC*Z
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
         W = DX(KX)
         Z = DY(KY)
         DX(KX) = DC*W+DS*Z
         DY(KY) = -DS*W+DC*Z
         KX = KX+INCX
         KY = KY+INCY
   15 CONTINUE
   20 CONTINUE
C
      RETURN
      END
