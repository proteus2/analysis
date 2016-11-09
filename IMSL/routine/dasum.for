C   IMSL ROUTINE NAME   - VBLA=DASUM
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - COMPUTE DOUBLE PRECISION SUM OF ABSOLUTE
C                           VALUES
C
C   USAGE               - FUNCTION DASUM (N,DX,INCX)
C
C   ARGUMENTS    DASUM  - DOUBLE PRECISION SUM FROM I=1 TO N OF
C                           DABS(X(I)). (OUTPUT)
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
      DOUBLE PRECISION FUNCTION DASUM (N,DX,INCX)
C
C                                  SPECIFICATIONS FOR ARGUMENTS
      DOUBLE PRECISION   DX(1)
      INTEGER            N,INCX
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,M,MP1,NS
C                                  FIRST EXECUTABLE STATEMENT
      DASUM = 0.D0
      IF (N.LE.0) RETURN
      IF (INCX.EQ.1) GO TO 10
C                                  CODE FOR INCREMENTS NOT EQUAL TO 1.
      NS = N*INCX
      DO 5 I=1,NS,INCX
         DASUM = DASUM+DABS(DX(I))
    5 CONTINUE
      RETURN
C                                  CODE FOR INCREMENTS EQUAL TO 1.
C                                    CLEAN-UP LOOP SO REMAINING VECTOR
C                                    LENGTH IS A MULTIPLE OF 6.
   10 M = N-(N/6)*6
      IF (M.EQ.0) GO TO 20
      DO 15 I=1,M
         DASUM = DASUM+DABS(DX(I))
   15 CONTINUE
      IF (N.LT.6) RETURN
   20 MP1 = M+1
      DO 25 I=MP1,N,6
         DASUM = DASUM+DABS(DX(I))+DABS(DX(I+1))+DABS(DX(I+2))
     1   +DABS(DX(I+3))+DABS(DX(I+4))+DABS(DX(I+5))
   25 CONTINUE
      RETURN
      END
