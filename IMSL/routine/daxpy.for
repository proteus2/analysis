C   IMSL ROUTINE NAME   - VBLA=DAXPY
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - COMPUTE A CONSTANT TIMES A VECTOR PLUS
C                           A VECTOR, ALL DOUBLE PRECISION
C
C   USAGE               - CALL DAXPY (N,DA,DX,INCX,DY,INCY)
C
C   ARGUMENTS    N      - LENGTH OF VECTORS X AND Y. (INPUT)
C                DA     - DOUBLE PRECISION SCALAR. (INPUT)
C                DX     - DOUBLE PRECISION VECTOR OF LENGTH
C                           MAX(N*IABS(INCX),1). (INPUT)
C                INCX   - DISPLACEMENT BETWEEN ELEMENTS OF DX. (INPUT)
C                           X(I) IS DEFINED TO BE..
C                           DX(1+(I-1)*INCX) IF INCX.GE.0 OR
C                           DX(1+(I-N)*INCX) IF INCX.LT.0.
C                DY     - DOUBLE PRECISION VECTOR OF LENGTH
C                           MAX(N*IABS(INCY),1). (INPUT/OUTPUT)
C                           DAXPY REPLACES Y(I) WITH DA*X(I)+Y(I) FOR
C                           I=1,...,N.
C                           X(I) AND Y(I) REFER TO SPECIFIC ELEMENTS
C                           OF DX AND DY, RESPECTIVELY. SEE INCX AND
C                           INCY ARGUMENT DESCRIPTIONS.
C                INCY   - DISPLACEMENT BETWEEN ELEMENTS OF DY. (INPUT)
C                           Y(I) IS DEFINED TO BE..
C                           DY(1+(I-1)*INCY) IF INCY.GE.0 OR
C                           DY(1+(I-N)*INCY) IF INCY.LT.0.
C
C   REQD. IMSL ROUTINES - NONE REQUIRED
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE DAXPY  (N,DA,DX,INCX,DY,INCY)
C
C                                  SPECIFICATIONS FOR ARGUMENTS
      DOUBLE PRECISION   DX(1),DY(1),DA
      INTEGER            N,INCX,INCY
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,IY,M,MP1,NS,IX
C                                  FIRST EXECUTABLE STATEMENT
      IF (N.LE.0.OR.DA.EQ.0.D0) RETURN
      IF (INCX.EQ.INCY) IF (INCX-1) 5,15,35
    5 CONTINUE
C                                  CODE FOR NONEQUAL OR NONPOSITIVE
C                                    INCREMENTS.
      IX = 1
      IY = 1
      IF (INCX.LT.0) IX = (-N+1)*INCX+1
      IF (INCY.LT.0) IY = (-N+1)*INCY+1
      DO 10 I=1,N
         DY(IY) = DY(IY)+DA*DX(IX)
         IX = IX+INCX
         IY = IY+INCY
   10 CONTINUE
      RETURN
C                                  CODE FOR BOTH INCREMENTS EQUAL TO 1
C                                    CLEAN-UP LOOP SO REMAINING VECTOR
C                                    LENGTH IS A MULTIPLE OF 4.
   15 M = N-(N/4)*4
      IF (M.EQ.0) GO TO 25
      DO 20 I=1,M
         DY(I) = DY(I)+DA*DX(I)
   20 CONTINUE
      IF (N.LT.4) RETURN
   25 MP1 = M+1
      DO 30 I=MP1,N,4
         DY(I) = DY(I)+DA*DX(I)
         DY(I+1) = DY(I+1)+DA*DX(I+1)
         DY(I+2) = DY(I+2)+DA*DX(I+2)
         DY(I+3) = DY(I+3)+DA*DX(I+3)
   30 CONTINUE
      RETURN
C                                  CODE FOR EQUAL, POSITIVE, NONUNIT
C                                    INCREMENTS.
   35 CONTINUE
      NS = N*INCX
      DO 40 I=1,NS,INCX
         DY(I) = DA*DX(I)+DY(I)
   40 CONTINUE
      RETURN
      END
