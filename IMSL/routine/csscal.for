C   IMSL ROUTINE NAME   - VBLA=CSSCAL
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - COMPUTE A REAL CONSTANT TIMES A COMPLEX
C                           VECTOR
C
C   USAGE               - CALL CSSCAL (N,SA,CX,INCX)
C
C   ARGUMENTS    N      - LENGTH OF VECTOR X. (INPUT)
C                SA     - REAL SCALAR. (INPUT)
C                CX     - COMPLEX VECTOR OF LENGTH N*INCX.
C                           (INPUT/OUTPUT)
C                           CSSCAL REPLACES X(I) WITH SA*X(I) FOR
C                           I=1,...,N.
C                           X(I) REFERS TO A SPECIFIC ELEMENT OF CX.
C                           SEE INCX ARGUMENT DESCRIPTION.
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
      SUBROUTINE CSSCAL (N,SA,CX,INCX)
C
C                                  SPECIFICATIONS FOR ARGUMENTS
      COMPLEX            CX(1)
      INTEGER            N,INCX
      REAL               SA
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,NS
C                                  FIRST EXECUTABLE STATEMENT
      IF (N.LE.0) RETURN
      NS = N*INCX
      DO 5 I=1,NS,INCX
         CX(I) = SA*CX(I)
    5 CONTINUE
      RETURN
      END
