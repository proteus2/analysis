C   IMSL ROUTINE NAME   - VBLA=SSWAP
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - INTERCHANGE VECTORS X AND Y, BOTH
C                           SINGLE PRECISION
C
C   USAGE               - CALL SSWAP (N,SX,INCX,SY,INCY)
C
C   ARGUMENTS    N      - LENGTH OF VECTORS X AND Y. (INPUT)
C                SX     - REAL VECTOR OF LENGTH MAX(N*IABS(INCX),1).
C                           (INPUT/OUTPUT)
C                           SSWAP INTERCHANGES X(I) AND Y(I) FOR
C                           I=1,...,N.
C                           X(I) AND Y(I) REFER TO SPECIFIC ELEMENTS
C                           OF SX AND SY, RESPECTIVELY. SEE INCX AND
C                           INCY ARGUMENT DESCRIPTIONS.
C                INCX   - DISPLACEMENT BETWEEN ELEMENTS OF SX. (INPUT)
C                           X(I) IS DEFINED TO BE..
C                           SX(1+(I-1)*INCX) IF INCX.GE.0 OR
C                           SX(1+(I-N)*INCX) IF INCX.LT.0.
C                SY     - REAL VECTOR OF LENGTH MAX(N*IABS(INCY),1).
C                           (INPUT/OUTPUT)
C                INCY   - DISPLACEMENT BETWEEN ELEMENTS OF SY. (INPUT)
C                           Y(I) IS DEFINED TO BE..
C                           SY(1+(I-1)*INCY) IF INCY.GE.0 OR
C                           SY(1+(I-N)*INCY) IF INCY.LT.0.
C
C   REQD. IMSL ROUTINES - NONE REQUIRED
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE SSWAP  (N,SX,INCX,SY,INCY)
C
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N,INCX,INCY
      REAL               SX(1),SY(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,IX,IY,M,MP1,NS
      REAL               STEMP1,STEMP2,STEMP3
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
         STEMP1 = SX(IX)
         SX(IX) = SY(IY)
         SY(IY) = STEMP1
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
         STEMP1 = SX(I)
         SX(I) = SY(I)
         SY(I) = STEMP1
   20 CONTINUE
      IF (N.LT.3) RETURN
   25 MP1 = M+1
      DO 30 I=MP1,N,3
         STEMP1 = SX(I)
         STEMP2 = SX(I+1)
         STEMP3 = SX(I+2)
         SX(I) = SY(I)
         SX(I+1) = SY(I+1)
         SX(I+2) = SY(I+2)
         SY(I) = STEMP1
         SY(I+1) = STEMP2
         SY(I+2) = STEMP3
   30 CONTINUE
      RETURN
   35 CONTINUE
C                                  CODE FOR EQUAL, POSITIVE, NONUNIT
C                                    INCREMENTS.
      NS = N*INCX
      DO 40 I=1,NS,INCX
         STEMP1 = SX(I)
         SX(I) = SY(I)
         SY(I) = STEMP1
   40 CONTINUE
      RETURN
      END
