C   IMSL ROUTINE NAME   - VBLA=CSCAL
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - COMPUTE A COMPLEX CONSTANT TIMES A
C                           COMPLEX VECTOR
C
C   USAGE               - CALL CSCAL (N,CA,CX,INCX)
C
C   ARGUMENTS    N      - LENGTH OF VECTOR X. (INPUT)
C                CA     - COMPLEX SCALAR. (INPUT)
C                CX     - COMPLEX VECTOR OF LENGTH N*INCX.
C                           (INPUT/OUTPUT)
C                           CSCAL REPLACES X(I) WITH CA*X(I) FOR
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
      SUBROUTINE CSCAL  (N,CA,CX,INCX)
C
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N,INCX
      COMPLEX            CA,CX(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,NS
C                                  FIRST EXECUTABLE STATEMENT
      IF (N.LE.0) RETURN
      NS = N*INCX
      DO 5 I=1,NS,INCX
         CX(I) = CA*CX(I)
    5 CONTINUE
      RETURN
      END
