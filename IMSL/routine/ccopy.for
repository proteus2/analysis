C   IMSL ROUTINE NAME   - VBLA=CCOPY
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - COPY A VECTOR X TO A VECTOR Y, BOTH
C                           COMPLEX
C
C   USAGE               - CALL CCOPY (N,CX,INCX,CY,INCY)
C
C   ARGUMENTS    N      - LENGTH OF VECTORS X AND Y. (INPUT)
C                CX     - COMPLEX VECTOR OF LENGTH MAX(N*IABS(INCX),1).
C                           (INPUT)
C                INCX   - DISPLACEMENT BETWEEN ELEMENTS OF CX. (INPUT)
C                           X(I) IS DEFINED TO BE..
C                           CX(1+(I-1)*INCX) IF INCX.GE.0 OR
C                           CX(1+(I-N)*INCX) IF INCX.LT.0.
C                CY     - COMPLEX VECTOR OF LENGTH MAX(N*IABS(INCY),1).
C                           (OUTPUT)
C                           CCOPY COPIES X(I) TO Y(I) FOR I=1,...,N.
C                           X(I) AND Y(I) REFER TO SPECIFIC ELEMENTS
C                           OF CX AND CY, RESPECTIVELY. SEE INCX AND
C                           INCY ARGUMENT DESCRIPTIONS.
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
      SUBROUTINE CCOPY  (N,CX,INCX,CY,INCY)
C
C                                  SPECIFICATIONS FOR ARGUMENTS
      COMPLEX            CX(1),CY(1)
      INTEGER            N,INCX,INCY
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,KX,KY,NS
C                                  FIRST EXECUTABLE STATEMENT
      IF (N.LE.0) RETURN
      IF (INCX.EQ.INCY.AND.INCX.GT.0) GO TO 10
      KX = 1
      KY = 1
      IF (INCX.LT.0) KX = 1+(1-N)*INCX
      IF (INCY.LT.0) KY = 1+(1-N)*INCY
      DO 5 I=1,N
         CY(KY) = CX(KX)
         KX = KX+INCX
         KY = KY+INCY
    5 CONTINUE
      RETURN
   10 CONTINUE
      NS = N*INCX
      DO 15 I=1,NS,INCX
         CY(I) = CX(I)
   15 CONTINUE
      RETURN
      END
