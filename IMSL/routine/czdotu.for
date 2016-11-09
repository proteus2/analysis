C   IMSL ROUTINE NAME   - VBLA=CZDOTU
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - COMPUTE COMPLEX DOT PRODUCT USING
C                           UNCONJUGATED VECTOR COMPONENTS (AND
C                           DOUBLE PRECISION ACCUMULATION)
C
C   USAGE               - FUNCTION CZDOTU (N,CX,INCX,CY,INCY)
C
C   ARGUMENTS    CZDOTU - COMPLEX SUM FROM I=1 TO N OF X(I)*Y(I).
C                           (OUTPUT)
C                           X(I) AND Y(I) REFER TO SPECIFIC ELEMENTS
C                           OF CX AND CY, RESPECTIVELY. SEE INCX AND
C                           INCY ARGUMENT DESCRIPTIONS.
C                N      - LENGTH OF VECTORS X AND Y. (INPUT)
C                CX     - COMPLEX VECTOR OF LENGTH MAX(N*IABS(INCX),1).
C                           (INPUT)
C                INCX   - DISPLACEMENT BETWEEN ELEMENTS OF CX. (INPUT)
C                           X(I) IS DEFINED TO BE..
C                           CX(1+(I-1)*INCX) IF INCX.GE.0 OR
C                           CX(1+(I-N)*INCX) IF INCX.LT.0.
C                CY     - COMPLEX VECTOR OF LENGTH MAX(N*IABS(INCY),1).
C                           (INPUT)
C                INCY   - DISPLACEMENT BETWEEN ELEMENTS OF CY. (INPUT)
C                           Y(I) IS DEFINED TO BE..
C                           CY(1+(I-1)*INCY) IF INCY.GE.0 OR
C                           CY(1+(I-N)*INCY) IF INCY.LT.0.
C
C   REQD. IMSL ROUTINES - NONE REQUIRED
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      COMPLEX FUNCTION CZDOTU (N,CX,INCX,CY,INCY)
C
C                                  SPECIFICATIONS FOR ARGUMENTS
      COMPLEX            CX(1),CY(1)
      INTEGER            N,INCX,INCY
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      DOUBLE PRECISION   A,B,C,D,SUMR,SUMI
      INTEGER            I,NS,KX,KY
      REAL               ZERO
      DATA               ZERO/0.0/
C                                  FIRST EXECUTABLE STATEMENT
      CZDOTU = CMPLX(ZERO,ZERO)
      IF (N.LE.0) RETURN
      SUMR = DBLE(ZERO)
      SUMI = DBLE(ZERO)
      IF (INCX.EQ.INCY.AND.INCX.GT.0) GO TO 10
      KX = 1
      KY = 1
      IF (INCX.LT.0) KX = 1+(1-N)*INCX
      IF (INCY.LT.0) KY = 1+(1-N)*INCY
      DO 5 I=1,N
         A = DBLE(REAL(CX(KX)))
         B = DBLE(AIMAG(CX(KX)))
         C = DBLE(REAL(CY(KY)))
         D = DBLE(AIMAG(CY(KY)))
         SUMR = SUMR+A*C-B*D
         SUMI = SUMI+A*D+B*C
         KX = KX+INCX
         KY = KY+INCY
    5 CONTINUE
      CZDOTU = CMPLX(SNGL(SUMR),SNGL(SUMI))
      RETURN
   10 CONTINUE
      NS = N*INCX
      DO 15 I=1,NS,INCX
         A = DBLE(REAL(CX(I)))
         B = DBLE(AIMAG(CX(I)))
         C = DBLE(REAL(CY(I)))
         D = DBLE(AIMAG(CY(I)))
         SUMR = SUMR+A*C-B*D
         SUMI = SUMI+A*D+B*C
   15 CONTINUE
      CZDOTU = CMPLX(SNGL(SUMR),SNGL(SUMI))
      RETURN
      END
