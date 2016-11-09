C   IMSL ROUTINE NAME   - VBLA=SDSDOT
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - COMPUTE SINGLE PRECISION DOT PRODUCT AND
C                           ADD A CONSTANT USING DOUBLE PRECISION
C                           ACCUMULATION
C
C   USAGE               - FUNCTION SDSDOT (N,SB,SX,INCX,SY,INCY)
C
C   ARGUMENTS    SDSDOT - SUM FROM I=1 TO N OF X(I)*Y(I)+SB. (OUTPUT)
C                           X(I) AND Y(I) REFER TO SPECIFIC ELEMENTS
C                           OF SX AND SY, RESPECTIVELY. SEE INCX AND
C                           INCY ARGUMENT DESCRIPTIONS.
C                N      - LENGTH OF VECTORS X AND Y. (INPUT)
C                SB     - REAL SCALAR. (INPUT)
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
      REAL FUNCTION SDSDOT (N,SB,SX,INCX,SY,INCY)
C
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N,INCX,INCY
      REAL               SX(1),SY(1),SB
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      DOUBLE PRECISION   DSDOT
      INTEGER            I,NS,KX,KY
C                                  FIRST EXECUTABLE STATEMENT
      DSDOT = DBLE(SB)
      IF (N.LE.0) GO TO 10
      IF (INCX.EQ.INCY.AND.INCX.GT.0) GO TO 15
      KX = 1
      KY = 1
      IF (INCX.LT.0) KX = 1+(1-N)*INCX
      IF (INCY.LT.0) KY = 1+(1-N)*INCY
      DO 5 I=1,N
         DSDOT = DSDOT+DBLE(SX(KX))*DBLE(SY(KY))
         KX = KX+INCX
         KY = KY+INCY
    5 CONTINUE
   10 SDSDOT = SNGL(DSDOT)
      RETURN
   15 CONTINUE
      NS = N*INCX
      DO 20 I=1,NS,INCX
         DSDOT = DSDOT+DBLE(SX(I))*DBLE(SY(I))
   20 CONTINUE
      SDSDOT = SNGL(DSDOT)
      RETURN
      END
 
R; T=0.02/0.19 22:12:28
