C   IMSL ROUTINE NAME   - VBLA=DSDOT
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - COMPUTE SINGLE PRECISION DOT PRODUCT USING
C                           DOUBLE PRECISION ACCUMULATION
C
C   USAGE               - FUNCTION DSDOT (N,SX,INCX,SY,INCY)
C
C   ARGUMENTS    DSDOT  - DOUBLE PRECISION SUM FROM I=1 TO N OF
C                           X(I)*Y(I). (OUTPUT)
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
      DOUBLE PRECISION FUNCTION DSDOT (N,SX,INCX,SY,INCY)
C
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N,INCX,INCY
      REAL               SX(1),SY(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,NS,KX,KY
C                                  FIRST EXECUTABLE STATEMENT
      DSDOT = 0.D0
      IF (N.LE.0) RETURN
      IF (INCX.EQ.INCY.AND.INCX.GT.0) GO TO 10
      KX = 1
      KY = 1
      IF (INCX.LT.0) KX = 1+(1-N)*INCX
      IF (INCY.LT.0) KY = 1+(1-N)*INCY
      DO 5 I=1,N
         DSDOT = DSDOT+DBLE(SX(KX))*DBLE(SY(KY))
         KX = KX+INCX
         KY = KY+INCY
    5 CONTINUE
      RETURN
   10 CONTINUE
      NS = N*INCX
      DO 15 I=1,NS,INCX
         DSDOT = DSDOT+DBLE(SX(I))*DBLE(SY(I))
   15 CONTINUE
      RETURN
      END
