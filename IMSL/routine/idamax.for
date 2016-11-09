C   IMSL ROUTINE NAME   - VBLA=IDAMAX
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - FIND THE SMALLEST INDEX OF THE MAXIMUM
C                           MAGNITUDE OF A DOUBLE PRECISION VECTOR
C
C   USAGE               - FUNCTION IDAMAX (N,DX,INCX)
C
C   ARGUMENTS    IDAMAX - THE SMALLEST INDEX I SUCH THAT DABS(X(I))
C                           IS THE MAXIMUM OF DABS(X(J)) FOR J=1 TO N.
C                           (OUTPUT)
C                           X(I) REFERS TO A SPECIFIC ELEMENT OF DX.
C                           SEE INCX ARGUMENT DESCRIPTION.
C                N      - LENGTH OF VECTOR X. (INPUT)
C                DX     - DOUBLE PRECISION VECTOR OF LENGTH N*INCX.
C                           (INPUT)
C                INCX   - DISPLACEMENT BETWEEN ELEMENTS OF DX. (INPUT)
C                           X(I) IS DEFINED TO BE DX(1+(I-1)*INCX).
C                           INCX MUST BE GREATER THAN ZERO.
C
C   REQD. IMSL ROUTINES - NONE REQUIRED
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      INTEGER FUNCTION IDAMAX (N,DX,INCX)
C                                    SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N,INCX
      DOUBLE PRECISION   DX(1)
C                                    SPECIFICATIONS FOR LOCAL VARIABLES
      DOUBLE PRECISION   DMAX, XMAG
      INTEGER            I,II,NS
C                                  FIRST EXECUTABLE STATEMENT
      IDAMAX = 0
      IF (N.LE.0) RETURN
      IDAMAX = 1
      IF (N.LE.1) RETURN
      IF (INCX.EQ.1) GO TO 15
C                                  CODE FOR INCREMENTS NOT EQUAL TO 1.
      DMAX = DABS(DX(1))
      NS = N*INCX
      II = 1
      DO 10 I=1,NS,INCX
         XMAG = DABS(DX(I))
         IF (XMAG.LE.DMAX) GO TO 5
         IDAMAX = II
         DMAX = XMAG
    5    II = II+1
   10 CONTINUE
      RETURN
C                                  CODE FOR INCREMENTS EQUAL TO 1.
   15 DMAX = DABS(DX(1))
      DO 20 I=2,N
         XMAG = DABS(DX(I))
         IF (XMAG.LE.DMAX) GO TO 20
         IDAMAX = I
         DMAX = XMAG
   20 CONTINUE
      RETURN
      END
