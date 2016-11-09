C   IMSL ROUTINE NAME   - VBLA=SDOT
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - COMPUTE SINGLE PRECISION DOT PRODUCT
C
C   USAGE               - FUNCTION SDOT (N,SX,INCX,SY,INCY)
C
C   ARGUMENTS    SDOT   - SUM FROM I=1 TO N OF X(I)*Y(I). (OUTPUT)
C                           X(I) AND Y(I) REFER TO SPECIFIC ELEMENTS
C                           OF SX AND SY, RESPECTIVELY. SEE INCX AND
C                           INCY ARGUMENT DESCRIPTIONS.
C                N      - LENGTH OF VECTORS X AND Y. (INPUT)
C                SX     - REAL VECTOR OF LENGTH MAX(N*IABS(INCX),1).
C                           (INPUT)
C                INCX   - DISPLACEMENT BETWEEN ELEMENTS OF SX. (INPUT)
C                           X(I) IS DEFINED TO BE..
C                           SX(1+(I-1)*INCX) IF INCX.GE.0 OR
C                           SX(1+(I-N)*INCX) IF INCX.LT.0.
C                SY     - REAL VECTOR OF LENGTH MAX(N*IABS(INCY),1).
C                           (INPUT)
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
      REAL FUNCTION SDOT (N,SX,INCX,SY,INCY)
C
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N,INCX,INCY
      REAL               SX(1),SY(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,M,MP1,NS,IX,IY
C                                  FIRST EXECUTABLE STATEMENT
      SDOT = 0.0E0
      IF (N.LE.0) RETURN
      IF (INCX.EQ.INCY) IF (INCX-1) 5,15,35
    5 CONTINUE
C                                  CODE FOR UNEQUAL INCREMENTS OR
C                                    NONPOSITIVE INCREMENTS.
      IX = 1
      IY = 1
      IF (INCX.LT.0) IX = (-N+1)*INCX+1
      IF (INCY.LT.0) IY = (-N+1)*INCY+1
      DO 10 I=1,N
         SDOT = SDOT+SX(IX)*SY(IY)
         IX = IX+INCX
         IY = IY+INCY
   10 CONTINUE
      RETURN
C                                  CODE FOR BOTH INCREMENTS EQUAL TO 1
C                                    CLEAN-UP LOOP SO REMAINING VECTOR
C                                    LENGTH IS A MULTIPLE OF 5.
   15 M = N-(N/5)*5
      IF (M.EQ.0) GO TO 25
      DO 20 I=1,M
         SDOT = SDOT+SX(I)*SY(I)
   20 CONTINUE
      IF (N.LT.5) RETURN
   25 MP1 = M+1
      DO 30 I=MP1,N,5
         SDOT = SDOT+SX(I)*SY(I)+SX(I+1)*SY(I+1)+SX(I+2)*SY(I+2)+SX(I
     1   +3)*SY(I+3)+SX(I+4)*SY(I+4)
   30 CONTINUE
      RETURN
C                                  CODE FOR POSITIVE EQUAL INCREMENTS
C                                    .NE.1.
   35 CONTINUE
      NS = N*INCX
      DO 40 I=1,NS,INCX
         SDOT = SDOT+SX(I)*SY(I)
   40 CONTINUE
      RETURN
      END
