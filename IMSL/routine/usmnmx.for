C   IMSL ROUTINE NAME   - USMNMX
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - DETERMINATION OF THE MINIMUM AND MAXIMUM
C                           VALUES OF A VECTOR
C
C   USAGE               - CALL USMNMX (X,N,INC,XMIN,XMAX)
C
C   ARGUMENTS    X      - INPUT VECTOR OF LENGTH N FROM WHICH MINIMUM,
C                           MAXIMUM VALUES ARE TO BE TAKEN.
C                N      - LENGTH OF THE INPUT VECTOR X. (INPUT)
C                INC    - DISPLACEMENT BETWEEN CONSECUTIVE VALUES OF X
C                           TO BE CONSIDERED.
C                XMIN   - OUTPUT SCALAR CONTAINING MINIMUM VALUE OF X.
C                XMAX   - OUTPUT SCALAR CONTAINING MAXIMUM VALUE OF X.
C
C   REQD. IMSL ROUTINES - NONE REQUIRED
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE USMNMX (X,N,INC,XMIN,XMAX)
C
      DIMENSION          X(N)
C                                  FIRST EXECUTABLE STATEMENT
      XMIN = X(1)
      XMAX = X(1)
      DO 10 I=1,N,INC
         IF (X(I) .GE. XMIN) GO TO 5
         XMIN = X(I)
         GO TO 10
   5     IF (X(I) .GT. XMAX) XMAX = X(I)
  10  CONTINUE
      RETURN
      END
