C   IMSL ROUTINE NAME   - VBLA=DSWAP
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - INTERCHANGE VECTORS X AND Y, BOTH
C                           DOUBLE PRECISION
C
C   USAGE               - CALL DSWAP (N,DX,INCX,DY,INCY)
C
C   ARGUMENTS    N      - LENGTH OF VECTORS X AND Y. (INPUT)
C                DX     - DOUBLE PRECISION VECTOR OF LENGTH
C                           MAX(N*IABS(INCX),1). (INPUT/OUTPUT)
C                           DSWAP INTERCHANGES X(I) WITH Y(I) FOR
C                           I=1,...,N.
C                           X(I) AND Y(I) REFER TO SPECIFIC ELEMENTS
C                           OF DX AND DY, RESPECTIVELY. SEE INCX AND
C                           INCY ARGUMENT DESCRIPTIONS.
C                INCX   - DISPLACEMENT BETWEEN ELEMENTS OF DX. (INPUT)
C                           X(I) IS DEFINED TO BE..
C                           DX(1+(I-1)*INCX) IF INCX.GE.0 OR
C                           DX(1+(I-N)*INCX) IF INCX.LT.0.
C                DY     - DOUBLE PRECISION VECTOR OF LENGTH
C                           MAX(N*IABS(INCY),1). (INPUT/OUTPUT)
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
      SUBROUTINE DSWAP  (N,DX,INCX,DY,INCY)
C
C                                  SPECIFICATIONS FOR ARGUMENTS
      DOUBLE PRECISION   DX(1),DY(1)
      INTEGER            N,INCX,INCY
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      DOUBLE PRECISION   DTEMP1,DTEMP2,DTEMP3
      INTEGER            I,IX,IY,M,MP1,NS
C                                  FIRST EXECUTABLE STATEMENT
      IF (N.LE.0) RETURN
      IF (INCX.EQ.INCY) IF (INCX-1) 5,15,35
    5 CONTINUE
C                                  CODE FOR UNEQUAL OR NONPOSITIVE
C                                    INCREMENTS.
      IX = 1
      IY = 1
      IF (INCX.LT.0) IX = (-N+1)*INCX+1
      IF (INCY.LT.0) IY = (-N+1)*INCY+1
      DO 10 I=1,N
         DTEMP1 = DX(IX)
         DX(IX) = DY(IY)
         DY(IY) = DTEMP1
         IX = IX+INCX
         IY = IY+INCY
   10 CONTINUE
      RETURN
C                                  CODE FOR BOTH INCREMENTS EQUAL TO 1
C                                    CLEAN-UP LOOP SO REMAINING VECTOR
C                                    LENGTH IS A MULTIPLE OF 3.
   15 M = N-(N/3)*3
      IF (M.EQ.0) GO TO 25
      DO 20 I=1,M
         DTEMP1 = DX(I)
         DX(I) = DY(I)
         DY(I) = DTEMP1
   20 CONTINUE
      IF (N.LT.3) RETURN
   25 MP1 = M+1
      DO 30 I=MP1,N,3
         DTEMP1 = DX(I)
         DTEMP2 = DX(I+1)
         DTEMP3 = DX(I+2)
         DX(I) = DY(I)
         DX(I+1) = DY(I+1)
         DX(I+2) = DY(I+2)
         DY(I) = DTEMP1
         DY(I+1) = DTEMP2
         DY(I+2) = DTEMP3
   30 CONTINUE
      RETURN
   35 CONTINUE
C                                  CODE FOR EQUAL, POSITIVE, NONUNIT
C                                    INCREMENTS.
      NS = N*INCX
      DO 40 I=1,NS,INCX
         DTEMP1 = DX(I)
         DX(I) = DY(I)
         DY(I) = DTEMP1
   40 CONTINUE
      RETURN
      END
