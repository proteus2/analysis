C   IMSL ROUTINE NAME   - VSENA
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1982
C
C   PURPOSE             - NUCLEUS CALLED BY IMSL ROUTINES WHICH SORT
C                           A MATRIX WITH RESPECT TO VECTORS.
C
C   REQD. IMSL ROUTINES - NONE REQUIRED
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE VSENA (KA,LA)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            LA,KA(LA)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,IJ,IL(21),IT,ITT,IU(21),J,K,L,M
      REAL               R
C                                  FIRST EXECUTABLE STATEMENT
      M = 1
      I = 1
      J = LA
      R = .375
      IF (LA.LE.0) RETURN
    5 IF (I.EQ.J) GO TO 45
      IF (R.GT..5898437) GO TO 10
      R = R+3.90625E-2
      GO TO 15
   10 R = R-.21875
   15 K = I
C                                  SELECT A CENTRAL ELEMENT OF THE
C                                  ARRAY AND SAVE IT IN LOCATION T
      IJ = I+(J-I)*R
      IT = KA(IJ)
C                                  IF FIRST ELEMENT OF ARRAY IS GREATER
C                                  THAN T, INTERCHANGE WITH T
      IF (KA(I).LE.IT) GO TO 20
      KA(IJ) = KA(I)
      KA(I) = IT
      IT = KA(IJ)
   20 L = J
C                                  IF LAST ELEMENT OF ARRAY IS LESS THAN
C                                  T, INTERCHANGE WITH T
      IF (KA(J).GE.IT) GO TO 30
      KA(IJ) = KA(J)
      KA(J) = IT
      IT = KA(IJ)
C                                  IF FIRST ELEMENT OF ARRAY IS GREATER
C                                  THAN T, INTERCHANGE WITH T
      IF (KA(I).LE.IT) GO TO 30
      KA(IJ) = KA(I)
      KA(I) = IT
      IT = KA(IJ)
      GO TO 30
   25 IF (KA(L).EQ.KA(K)) GO TO 30
      ITT = KA(L)
      KA(L) = KA(K)
      KA(K) = ITT
C                                  FIND AN ELEMENT IN THE SECOND HALF OF
C                                  THE ARRAY WHICH IS SMALLER THAN T
   30 L = L-1
      IF (KA(L).GT.IT) GO TO 30
C                                  FIND AN ELEMENT IN THE FIRST HALF OF
C                                  THE ARRAY WHICH IS GREATER THAN T
   35 K = K+1
      IF (KA(K).LT.IT) GO TO 35
C                                  INTERCHANGE THESE ELEMENTS
      IF (K.LE.L) GO TO 25
C                                  SAVE UPPER AND LOWER SUBSCRIPTS OF
C                                  THE ARRAY YET TO BE SORTED
      IF (L-I.LE.J-K) GO TO 40
      IL(M) = I
      IU(M) = L
      I = K
      M = M+1
      GO TO 50
   40 IL(M) = K
      IU(M) = J
      J = L
      M = M+1
      GO TO 50
C                                  BEGIN AGAIN ON ANOTHER PORTION OF
C                                  THE UNSORTED ARRAY
   45 M = M-1
      IF (M.EQ.0) RETURN
      I = IL(M)
      J = IU(M)
   50 IF (J-I.GE.11) GO TO 15
      IF (I.EQ.1) GO TO 5
      I = I-1
   55 I = I+1
      IF (I.EQ.J) GO TO 45
      IT = KA(I+1)
      IF (KA(I).LE.IT) GO TO 55
      K = I
   60 KA(K+1) = KA(K)
      K = K-1
      IF (IT.LT.KA(K)) GO TO 60
      KA(K+1) = IT
      GO TO 55
      END
