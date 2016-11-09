C   IMSL ROUTINE NAME   - VBLA=ISAMAX
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - FIND THE SMALLEST INDEX OF THE MAXIMUM
C                           MAGNITUDE OF A SINGLE PRECISION VECTOR
C
C   USAGE               - FUNCTION ISAMAX (N,SX,INCX)
C
C   ARGUMENTS    ISAMAX - THE SMALLEST INDEX I SUCH THAT ABS(X(I))
C                           IS THE MAXIMUM OF ABS(X(J)) FOR J=1 TO N.
C                           (OUTPUT)
C                           X(I) REFERS TO A SPECIFIC ELEMENT OF SX.
C                           SEE INCX ARGUMENT DESCRIPTION.
C                N      - LENGTH OF VECTOR X. (INPUT)
C                SX     - REAL VECTOR OF LENGTH N*INCX. (INPUT)
C                INCX   - DISPLACEMENT BETWEEN ELEMENTS OF SX. (INPUT)
C                           X(I) IS DEFINED TO BE SX(1+(I-1)*INCX).
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
      INTEGER FUNCTION ISAMAX (N,SX,INCX)
C
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N,INCX
      REAL               SX(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,II,NS
      REAL               SMAX,XMAG
C                                  FIRST EXECUTABLE STATEMENT
      ISAMAX = 0
      IF (N.LE.0) RETURN
      ISAMAX = 1
      IF (N.LE.1) RETURN
      IF (INCX.EQ.1) GO TO 15
C                                  CODE FOR INCREMENTS NOT EQUAL TO 1.
      SMAX = ABS(SX(1))
      NS = N*INCX
      II = 1
      DO 10 I=1,NS,INCX
         XMAG = ABS(SX(I))
         IF (XMAG.LE.SMAX) GO TO 5
         ISAMAX = II
         SMAX = XMAG
    5    II = II+1
   10 CONTINUE
      RETURN
C                                  CODE FOR INCREMENTS EQUAL TO 1.
   15 SMAX = ABS(SX(1))
      DO 20 I=2,N
         XMAG = ABS(SX(I))
         IF (XMAG.LE.SMAX) GO TO 20
         ISAMAX = I
         SMAX = XMAG
   20 CONTINUE
      RETURN
      END
