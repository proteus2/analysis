C   IMSL ROUTINE NAME   - VBLA=SASUM
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - COMPUTE SINGLE PRECISION SUM OF ABSOLUTE
C                           VALUES
C
C   USAGE               - FUNCTION SASUM (N,SX,INCX)
C
C   ARGUMENTS    SASUM  - SUM FROM I=1 TO N OF ABS(X(I)). (OUTPUT)
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
      REAL FUNCTION SASUM (N,SX,INCX)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N,INCX
      REAL               SX(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,M,MP1,NS
C                                  FIRST EXECUTABLE STATEMENT
      SASUM = 0.0E0
      IF (N.LE.0) RETURN
      IF (INCX.EQ.1) GO TO 10
C                                  CODE FOR INCREMENTS NOT EQUAL TO 1.
      NS = N*INCX
      DO 5 I=1,NS,INCX
         SASUM = SASUM+ABS(SX(I))
    5 CONTINUE
      RETURN
C                                  CODE FOR INCREMENTS EQUAL TO 1.
C                                    CLEAN-UP LOOP SO REMAINING VECTOR
C                                    LENGTH IS A MULTIPLE OF 6.
   10 M = N-(N/6)*6
      IF (M.EQ.0) GO TO 20
      DO 15 I=1,M
         SASUM = SASUM+ABS(SX(I))
   15 CONTINUE
      IF (N.LT.6) RETURN
   20 MP1 = M+1
      DO 25 I=MP1,N,6
         SASUM = SASUM+ABS(SX(I))+ABS(SX(I+1))+ABS(SX(I+2))+ABS(SX(I
     1   +3))+ABS(SX(I+4))+ABS(SX(I+5))
   25 CONTINUE
      RETURN
      END
