C   IMSL ROUTINE NAME   - VSRTAD
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - SORTING OF DOUBLE PRECISION ARRAYS BY
C                           ALGEBRAIC VALUE
C
C   USAGE               - CALL VSRTAD (A,LA)
C
C   ARGUMENTS    A      - ON INPUT, A CONTAINS THE DOUBLE PRECISION
C                           ARRAY TO BE SORTED.
C                         ON OUTPUT, A CONTAINS THE SORTED DOUBLE
C                           PRECISION ARRAY.
C                LA     - INPUT VARIABLE CONTAINING THE NUMBER OF
C                           ELEMENTS IN THE ARRAY TO BE SORTED.
C
C   REQD. IMSL ROUTINES - NONE REQUIRED
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE VSRTAD (A,LA)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            LA
      DOUBLE PRECISION   A(LA)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      REAL               R
      DOUBLE PRECISION   T,TT
      INTEGER            IU(21),IL(21),I,M,J,K,IJ,L
C                                  FIRST EXECUTABLE STATEMENT
      M=1
      I=1
      J=LA
      R=.375
      IF (LA.LE.0) RETURN
   10 IF (I .EQ. J) GO TO 55
   15 IF (R .GT. .5898437) GO TO 20
      R=R+3.90625E-2
      GO TO 25
   20 R=R-.21875
   25 K=I
C                                  SELECT A CENTRAL ELEMENT OF THE
C                                  ARRAY AND SAVE IT IN LOCATION T
      IJ=I+(J-I)*R
      T=A(IJ)
C                                  IF FIRST ELEMENT OF ARRAY IS GREATER
C                                  THAN T, INTERCHANGE WITH T
      IF (A(I) .LE. T) GO TO 30
      A(IJ)=A(I)
      A(I)=T
      T=A(IJ)
   30 L=J
C                                  IF LAST ELEMENT OF ARRAY IS LESS THAN
C                                  T, INTERCHANGE WITH T
      IF (A(J) .GE. T) GO TO 40
      A(IJ)=A(J)
      A(J)=T
      T=A(IJ)
C                                  IF FIRST ELEMENT OF ARRAY IS GREATER
C                                  THAN T, INTERCHANGE WITH T
      IF (A(I) .LE. T) GO TO 40
      A(IJ)=A(I)
      A(I)=T
      T=A(IJ)
      GO TO 40
   35 IF(A(L).EQ.A(K)) GO TO 40
      TT=A(L)
      A(L)=A(K)
      A(K)=TT
C                                  FIND AN ELEMENT IN THE SECOND HALF OF
C                                  THE ARRAY WHICH IS SMALLER THAN T
   40 L=L-1
      IF (A(L) .GT. T) GO TO 40
C                                  FIND AN ELEMENT IN THE FIRST HALF OF
C                                  THE ARRAY WHICH IS GREATER THAN T
   45 K=K+1
      IF (A(K) .LT. T) GO TO 45
C                                  INTERCHANGE THESE ELEMENTS
      IF (K .LE. L) GO TO 35
C                                  SAVE UPPER AND LOWER SUBSCRIPTS OF
C                                  THE ARRAY YET TO BE SORTED
      IF (L-I .LE. J-K) GO TO 50
      IL(M)=I
      IU(M)=L
      I=K
      M=M+1
      GO TO 60
   50 IL(M)=K
      IU(M)=J
      J=L
      M=M+1
      GO TO 60
C                                  BEGIN AGAIN ON ANOTHER PORTION OF
C                                  THE UNSORTED ARRAY
   55 M=M-1
      IF (M .EQ. 0) RETURN
      I=IL(M)
      J=IU(M)
   60 IF (J-I .GE. 11) GO TO 25
      IF (I .EQ. 1) GO TO 10
      I=I-1
   65 I=I+1
      IF (I .EQ. J) GO TO 55
      T=A(I+1)
      IF (A(I) .LE. T) GO TO 65
      K=I
   70 A(K+1)=A(K)
      K=K-1
      IF (T .LT. A(K)) GO TO 70
      A(K+1)=T
      GO TO 65
      END
