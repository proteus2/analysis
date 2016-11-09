C   IMSL ROUTINE NAME   - LUDECP
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - DECOMPOSITION OF A POSITIVE DEFINITE MATRIX -
C                           SYMMETRIC STORAGE MODE
C
C   USAGE               - CALL LUDECP (A,UL,N,D1,D2,IER)
C
C   ARGUMENTS    A      - INPUT VECTOR OF LENGTH N(N+1)/2 CONTAINING
C                           THE N BY N POSITIVE DEFINITE SYMMETRIC
C                           MATRIX STORED IN SYMMETRIC STORAGE MODE.
C                UL     - OUTPUT VECTOR OF LENGTH N(N+1)/2 CONTAINING
C                           THE DECOMPOSED MATRIX L SUCH THAT A = L*
C                           L-TRANSPOSE. L IS STORED IN SYMMETRIC
C                           STORAGE MODE. THE DIAGONAL OF L CONTAINS THE
C                           RECIPROCALS OF THE ACTUAL DIAGONAL ELEMENTS.
C                N      - ORDER OF A. (INPUT)
C                D1,D2  - COMPONENTS OF THE DETERMINANT OF A.
C                           DETERMINANT(A) = D1*2.**D2. (OUTPUT)
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER = 129 INDICATES THAT MATRIX A IS
C                             ALGORITHMICALLY NOT POSITIVE DEFINITE.
C                             (SEE THE CHAPTER L PRELUDE).
C
C   REQD. IMSL ROUTINES - UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE LUDECP (A,UL,N,D1,D2,IER)
C
      DIMENSION          A(1),UL(1)
      DATA               ZERO,ONE,FOUR,SIXTN,SIXTH/
     *                   0.0,1.,4.,16.,.0625/
C                                  FIRST EXECUTABLE STATEMENT
      D1=ONE
      D2=ZERO
      RN = ONE/(N*SIXTN)
      IP = 1
      IER=0
      DO 45 I = 1,N
         IQ = IP
         IR = 1
         DO 40 J = 1,I
            X = A(IP)
            IF (J .EQ. 1) GO TO 10
            DO 5  K=IQ,IP1
               X = X - UL(K) * UL(IR)
               IR = IR+1
    5       CONTINUE
   10       IF (I.NE.J) GO TO 30
            D1 = D1*X
            IF (A(IP) + X*RN .LE. A(IP)) GO TO 50
   15       IF (ABS(D1).LE.ONE) GO TO 20
            D1 = D1 * SIXTH
            D2 = D2 + FOUR
            GO TO 15
   20       IF (ABS(D1) .GE. SIXTH) GO TO 25
            D1 = D1 * SIXTN
            D2 = D2 - FOUR
            GO TO 20
   25       UL(IP) = ONE/SQRT(X)
            GO TO 35
   30       UL(IP) = X * UL(IR)
   35       IP1 = IP
            IP = IP+1
            IR = IR+1
   40    CONTINUE
   45 CONTINUE
      GO TO 9005
   50 IER = 129
 9000 CONTINUE
      CALL UERTST(IER,6HLUDECP)
 9005 RETURN
      END
