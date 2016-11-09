C   IMSL ROUTINE NAME   - VBLA=SCASUM
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - COMPUTE COMPLEX SUM OF ABSOLUTE VALUES
C
C   USAGE               - FUNCTION SCASUM (N,CX,INCX)
C
C   ARGUMENTS    SCASUM - SUM FROM I=1 TO N OF ABS(REAL X(I))*
C                           ABS(AIMAG(X(I))). (OUTPUT)
C                           X(I) REFERS TO A SPECIFIC ELEMENT OF CX.
C                           SEE INCX ARGUMENT DESCRIPTION.
C                N      - LENGTH OF VECTOR X. (INPUT)
C                CX     - COMPLEX VECTOR OF LENGTH N*INCX. (INPUT)
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
      REAL FUNCTION SCASUM (N,CX,INCX)
C
C                                  SPECIFICATIONS FOR ARGUMENTS
      COMPLEX            CX(1)
      INTEGER            N,INCX
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,NS
C                                  FIRST EXECUTABLE STATEMENT
      SCASUM = 0.
      IF (N.LE.0) RETURN
      NS = N*INCX
      DO 5 I=1,NS,INCX
         SCASUM = SCASUM+ABS(REAL(CX(I)))+ABS(AIMAG(CX(I)))
    5 CONTINUE
      RETURN
      END
