C   IMSL ROUTINE NAME   - VMULQB
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - MATRIX MULTIPLICATION (BAND SYMMETRIC BY BAND
C                           MATRICES)
C
C   USAGE               - CALL VMULQB (A,IA,B,IB,N,C,IC)
C
C   ARGUMENTS    A      - N(1) BY N(1) BAND SYMMETRIC MATRIX WITH N(2)
C                           UPPER OR LOWER CODIAGONALS. A IS STORED IN
C                           BAND SYMMETRIC STORAGE MODE AND THEREFORE
C                           HAS DIMENSION N(1) BY (N(2)+1). (INPUT)
C                IA     - ROW DIMENSION OF MATRIX A EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM. (INPUT)
C                B      - N(1) BY N(1) BAND MATRIX WITH N(3) LOWER
C                           CODIAGONALS AND N(4) UPPER CODIAGONALS. B
C                           IS STORED IN BAND STORAGE MODE AND
C                           THEREFORE HAS DIMENSION N(1) BY
C                           (N(3)+N(4)+1). (INPUT)
C                IB     - ROW DIMENSION OF MATRIX B EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM. (INPUT)
C                N      - INPUT VECTOR OF LENGTH 4 CONTAINING THE
C                           FOLLOWING INFORMATION:
C                         N(1) CONTAINS THE NUMBER OF ROWS AND COLUMNS
C                           IN MATRICES A,B,AND C.
C                         N(2) CONTAINS THE NUMBER OF UPPER OR LOWER
C                           CODIAGONALS IN MATRIX A.
C                         N(3) CONTAINS THE NUMBER OF LOWER CODIAGONALS
C                           IN MATRIX B.
C                         N(4) CONTAINS THE NUMBER OF UPPER CODIAGONALS
C                           IN MATRIX B.
C                C      - OUTPUT MATRIX CONTAINING THE PRODUCT C = A*B.
C                           C IS A BAND MATRIX WITH M LOWER CODIAGONALS,
C                           WHERE M = MIN(N(1)-1,N(2)+N(3)) AND L
C                           UPPER CODIAGONALS, WHERE L = MIN(N(1)-1,
C                           N(2)+N(4)). C IS STORED IN BAND STORAGE
C                           MODE AND THEREFORE HAS DIMENSION N(1) BY
C                           (L+M+1). (OUTPUT)
C                IC     - ROW DIMENSION OF MATRIX C EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM. (INPUT)
C
C   REQD. IMSL ROUTINES - NONE REQUIRED
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE VMULQB (A,IA,B,IB,N,C,IC)
C
      INTEGER            N(4)
      REAL               A(IA,1),B(IB,1),C(IC,1),ZERO
      DOUBLE PRECISION   SUM
      EQUIVALENCE        (NLC,M),(NUC,L)
      DATA               ZERO/0.0/
C                                  FIRST EXECUTABLE STATEMENT
      NN = N(1)
      NCA = N(2)
      NLB = N(3)
      NUB = N(4)
      NCP1 = NCA+1
      NCP2 = NCA+2
      NCB = NLB + 1 + NUB
      M = MIN0(NCA+NLB,NN-1)
      L = MIN0(NCA+NUB,NN-1)
      NCC = M + L + 1
      NLCP1 = NLC + 1
      NUBP1 = NUB + 1
C                                  COMPUTE 1 ROW OF MATRIX C
      DO 75 IRC = 1,NN
C                                  INITIALIZE C-COL AND B-VECTOR
         JCC = 1
         IVB = 1
         LDZROC = NLCP1 - IRC
         IF (LDZROC) 15,20,5
C                                  INSERT LEADING ZEROES IN ROW OF C
C                                  (AS COMPRESSED) AND INCREMENT
C                                    COL. INDEX
    5    DO 10 JCC = 1,LDZROC
            C(IRC,JCC) = ZERO
   10    CONTINUE
         JCC = LDZROC + 1
         GO TO 20
C                                  INCREMENT B-VECTOR FOR COMPRESSION
C                                    OF C
   15    IVB = 1 - LDZROC
C                                  COMPUTE NO. OF LEADING ZEROES IN
C                                    NEXT A-VECTOR
   20    LDZROA = MAX0(IRC-NCP1,0)
         NCP2X = NCP2 - IRC
C                                  COMPUTE 1 ELEMENT OF C, USING 1 ROW
C                                    OF A AND 1 COL OF B
         DO 65 NVB = IVB,NN
C                                  INITIALIZE C-ELEMENT SUM
            SUM = 0.0D0
C                                  SET INDICES FOR FIRST NON-ZERO
C                                    ELEMENT OF A-VECTOR
            IRA = IRC
            JCA = MAX0(NCP2X,1)
C                                  COMPUTE NO. OF LEADING ZEROES IN
C                                    NEXT B-VECTOR
            LDZROB = MAX0(NVB-NUBP1,0)
C                                  SET INDICES FOR FIRST NON-ZERO
C                                    ELEMENT OF NEXT B-VECTOR
            IRB = LDZROB + 1
            JCB = MIN0(NVB+NLB,NCB)
C                                  ADJUST A OR B INDICES FOR FIRST SET
C                                    OF FACTORS (IF UNEQUAL NO. OF
C                                    LEADING ZEROES)
            NLZBA = LDZROB-LDZROA
            IF (NLZBA) 40,45,25
   25       DO 35 INCA = 1,NLZBA
               IF ((IRA.NE.IRC).OR.(JCA.GE.NCP1)) GO TO 30
               JCA = JCA + 1
               GO TO 35
   30          JCA = JCA - 1
               IRA = IRA + 1
   35       CONTINUE
            GO TO 45
   40       JCB = JCB+NLZBA
            IRB = IRB-NLZBA
C                                  MULTIPLY 2 ELEMENTS (IF END OF
C                                    VECTORS NOT EXCEEDED) AND
C                                    ACCUMULATE INTO SUM
   45       IF ((JCA.LT.1).OR.(IRA.GT.NN).OR.(JCB.LT.1).OR.(IRB.GT.NN))
     *        GO TO 60
            SUM = SUM + DBLE(A(IRA,JCA)) * DBLE(B(IRB,JCB))
C                                  INCREMENT A AND B INDICES
            IF ((IRA.NE.IRC).OR.(JCA.GE.NCP1)) GO TO 50
            JCA = JCA + 1
            GO TO 55
   50       JCA = JCA - 1
            IRA = IRA + 1
   55       JCB = JCB - 1
            IRB = IRB + 1
            GO TO 45
C                                  STORE C-ELEMENT
   60       C(IRC,JCC) = SUM
C                                  INCREASE C-COL INDEX
            JCC = JCC + 1
         IF (JCC .GT. NCC) GO TO 75
   65    CONTINUE
C                                  INSERT TRAILING ZEROES IN C-ROW IF
C                                    REQUIRED
         DO 70 ICC = JCC,NCC
            C(IRC,ICC) = ZERO
   70    CONTINUE
   75 CONTINUE
      RETURN
      END
