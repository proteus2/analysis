C   IMSL ROUTINE NAME   - VBLA=DSCAL
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - COMPUTE A DOUBLE PRECISION CONSTANT
C                           TIMES A DOUBLE PRECISION VECTOR
C
C   USAGE               - CALL DSCAL (N,DA,DX,INCX)
C
C   ARGUMENTS    N      - LENGTH OF VECTOR X. (INPUT)
C                DA     - DOUBLE PRECISION SCALAR. (INPUT)
C                DX     - DOUBLE PRECISION VECTOR OF LENGTH N*INCX.
C                           (INPUT/OUTPUT)
C                           DSCAL REPLACES X(I) WITH DA*X(I) FOR
C                           I=1,...,N.
C                           X(I) REFERS TO A SPECIFIC ELEMENT OF DX.
C                           SEE INCX ARGUMENT DESCRIPTION.
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
      SUBROUTINE DSCAL  (N,DA,DX,INCX)
C
C                                  SPECIFICATIONS FOR ARGUMENTS
      DOUBLE PRECISION   DA,DX(1)
      INTEGER            N,INCX
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,M,MP1,NS
C                                  FIRST EXECUTABLE STATEMENT
      IF (N.LE.0) RETURN
      IF (INCX.EQ.1) GO TO 10
C                                  CODE FOR INCREMENTS NOT EQUAL TO 1.
      NS = N*INCX
      DO 5 I=1,NS,INCX
         DX(I) = DA*DX(I)
    5 CONTINUE
      RETURN
C                                  CODE FOR INCREMENTS EQUAL TO 1.
C                                    CLEAN-UP LOOP SO REMAINING VECTOR
C                                    LENGTH IS A MULTIPLE OF 5.
   10 M = N-(N/5)*5
      IF (M.EQ.0) GO TO 20
      DO 15 I=1,M
         DX(I) = DA*DX(I)
   15 CONTINUE
      IF (N.LT.5) RETURN
   20 MP1 = M+1
      DO 25 I=MP1,N,5
         DX(I) = DA*DX(I)
         DX(I+1) = DA*DX(I+1)
         DX(I+2) = DA*DX(I+2)
         DX(I+3) = DA*DX(I+3)
         DX(I+4) = DA*DX(I+4)
   25 CONTINUE
      RETURN
      END
