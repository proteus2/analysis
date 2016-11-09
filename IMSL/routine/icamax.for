C   IMSL ROUTINE NAME   - VBLA=ICAMAX
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - FIND THE SMALLEST INDEX OF THE MAXIMUM
C                           MAGNITUDE OF A COMPLEX VECTOR
C
C   USAGE               - FUNCTION ICAMAX (N,CX,INCX)
C
C   ARGUMENTS    ICAMAX - THE SMALLEST INDEX I SUCH THAT ABS(X(I))
C                           IS THE MAXIMUM OF ABS(X(J)) FOR J=1 TO N.
C                           (OUTPUT)
C                           X(I) REFERS TO A SPECIFIC ELEMENT OF CX.
C                           SEE INCX ARGUMENT DESCRIPTION.
C                N      - LENGTH OF VECTOR X. (INPUT)
C                CX     - COMPLEX VECTOR OF LENGTH N*INCX.(INPUT)
C                INCX   - DISPLACEMENT BETWEEN ELEMENTS OF CX. (INPUT)
C                           X(I) IS DEFINED TO BE CX(1+(I-1)*INCX).
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
      INTEGER FUNCTION ICAMAX (N,CX,INCX)
C                                    SPECIFICATIONS FOR ARGUMENTS
      COMPLEX            CX(1)
      INTEGER            N,INCX
C                                    SPECIFICATIONS FOR LOCAL VARIABLES
      REAL               SUMMAX,SUMRI
      INTEGER            I,II,NS
C                                  FIRST EXECUTABLE STATEMENT
      ICAMAX = 0
      IF (N.LE.0) RETURN
      ICAMAX = 1
      IF (N.LE.1) RETURN
      NS = N*INCX
      II = 1
      SUMMAX = ABS(REAL(CX(1)))+ABS(AIMAG(CX(1)))
      DO 10 I=1,NS,INCX
         SUMRI = ABS(REAL(CX(1)))+ABS(AIMAG(CX(I)))
         IF (SUMMAG.GE.SUMRI) GO TO 5
         SUMMAX = SUMRI
         ICAMAX = II
    5    II = II+1
   10 CONTINUE
      RETURN
      END
