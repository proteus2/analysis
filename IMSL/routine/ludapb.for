C   IMSL ROUTINE NAME   - LUDAPB
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - DECOMPOSITION OF A POSITIVE DEFINITE BAND
C                           SYMMETRIC MATRIX - BAND SYMMETRIC STORAGE
C                           MODE
C
C   USAGE               - CALL LUDAPB (A,N,NC,IA,UL,IU,D1,D2,IER)
C
C   ARGUMENTS    A      - N BY N POSITIVE DEFINITE BAND SYMMETRIC
C                           MATRIX STORED IN BAND SYMMETRIC STORAGE
C                           MODE. A SHOULD BE DIMENSIONED AT LEAST
C                           N BY NC+1. (INPUT)
C                N      - ORDER OF A. (INPUT)
C                NC     - NUMBER OF UPPER OR LOWER CODIAGONALS OF A.
C                           (INPUT)
C                IA     - ROW DIMENSION OF MATRIX A EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM. (INPUT)
C                UL     - OUTPUT MATRIX L WHERE A = L*L-TRANSPOSE. L IS
C                           STORED IN BAND STORAGE MODE. UL SHOULD BE
C                           DIMENSIONED AT LEAST N BY NC+1. NOTE - THE
C                           DIAGONAL OF UL CONTAINS THE RECIPROCALS
C                           OF THE ACTUAL DIAGONAL ELEMENTS.
C                IU     - ROW DIMENSION OF MATRIX UL EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM. (INPUT)
C                D1     - COMPONENTS OF THE DETERMINANT OF A.
C                D2         DETERMINANT(A) = D1*2.**D2. (OUTPUT)
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           N = 129 INDICATES THAT THE MATRIX A IS
C                             ALGORITHMICALLY NOT POSITIVE DEFINITE.
C
C   REQD. IMSL ROUTINES - UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE LUDAPB (A,N,NC,IA,UL,IU,D1,D2,IER)
C
      REAL               ZERO,A(IA,1),UL(IU,1),D1,D2,ONE,HALF,XINF,SUM,
     *                   RN,FOUR,SIXTN,SIXTH
      DATA               ZERO/0.0/,FOUR/4.0/,SIXTN/16./,
     *                   SIXTH/.0625/,ONE/1.0/
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
      RN = ONE/(N*SIXTN)
      D1 = ONE
      D2 = ZERO
      NCP1 = NC+1
      IF (NC .EQ. 0) GO TO 15
C                                  INITIALIZE ZERO ELEMENTS
      DO 10 I = 1,NC
         DO 5 J = I,NC
            K = NCP1-J
            UL(I,K) = ZERO
    5    CONTINUE
   10 CONTINUE
C                                  I IS ROW INDEX OF ELEMENT BEING
C                                  COMPUTED
   15 DO 60 I = 1,N
         IMNCP1 = I-NCP1
         I1 = MAX0(1,1-IMNCP1)
C                                  J IS COLUMN INDEX OF ELEMENT BEING
C                                  COMPUTED
         DO 60 J = I1,NCP1
C                                  L IS ROW INDEX OF PREVIOUSLY COMPUTED
C                                  VECTOR BEING USED TO COMPUTE INNER
C                                  PRODUCT
            L = IMNCP1+J
            I2 = NCP1-J
            SUM = A(I,J)
            JM1 = J-1
            IF (JM1) 30,30,20
   20       DO 25 K = 1,JM1
C                                  M IS COLUMN INDEX
               M = I2+K
               SUM = SUM-UL(I,K)*UL(L,M)
   25       CONTINUE
   30       IF (J .NE. NCP1) GO TO 55
            IF(A(I,J)+SUM*RN .LE. A(I,J))GO TO 65
            UL(I,J) = ONE/SQRT(SUM)
C                                  UPDATE THE DETERMINANT
            D1 = D1*SUM
   35       IF (ABS(D1)-ONE) 45,45,40
   40       D1 = D1*SIXTH
            D2 = D2+FOUR
            GO TO 35
   45       IF (ABS(D1)-SIXTH) 50,50,60
   50       D1 = D1*SIXTN
            D2 = D2-FOUR
            GO TO 45
   55       UL(I,J) = SUM*UL(L,NCP1)
   60 CONTINUE
      GO TO 9005
   65 IER = 129
 9000 CONTINUE
      CALL UERTST(IER,6HLUDAPB)
 9005 RETURN
      END
