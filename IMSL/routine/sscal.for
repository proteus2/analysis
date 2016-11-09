C   IMSL ROUTINE NAME   - VBLA=SSCAL
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - COMPUTE A SINGLE PRECISION CONSTANT
C                           TIMES A SINGLE PRECISION VECTOR
C
C   USAGE               - CALL SSCAL (N,SA,SX,INCX)
C
C   ARGUMENTS    N      - LENGTH OF VECTOR X. (INPUT)
C                SA     - REAL SCALAR. (INPUT)
C                SX     - REAL VECTOR OF LENGTH N*INCX. (INPUT/OUTPUT)
C                           SSCAL REPLACES X(I) WITH SA*X(I) FOR
C                           I=1,...,N.
C                           X(I) REFERS TO A SPECIFIC ELEMENT OF SX.
C                           SEE INCX ARGUMENT DESCRIPTION.
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
      SUBROUTINE SSCAL  (N,SA,SX,INCX)
C
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            INCX,N
      REAL               SA,SX(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,M,MP1,NS
C                                  FIRST EXECUTABLE STATEMENT
      IF (N.LE.0) RETURN
      IF (INCX.EQ.1) GO TO 10
C                                  CODE FOR INCREMENTS NOT EQUAL TO 1.
      NS = N*INCX
      DO 5 I=1,NS,INCX
         SX(I) = SA*SX(I)
    5 CONTINUE
      RETURN
C                                  CODE FOR INCREMENTS EQUAL TO 1.
C                                    CLEAN-UP LOOP SO REMAINING VECTOR
C                                    LENGTH IS A MULTIPLE OF 5.
   10 M = N-(N/5)*5
      IF (M.EQ.0) GO TO 20
      DO 15 I=1,M
         SX(I) = SA*SX(I)
   15 CONTINUE
      IF (N.LT.5) RETURN
   20 MP1 = M+1
      DO 25 I=MP1,N,5
         SX(I) = SA*SX(I)
         SX(I+1) = SA*SX(I+1)
         SX(I+2) = SA*SX(I+2)
         SX(I+3) = SA*SX(I+3)
         SX(I+4) = SA*SX(I+4)
   25 CONTINUE
      RETURN
      END
