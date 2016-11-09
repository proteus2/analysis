C   IMSL ROUTINE NAME   - VMULBB
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - MATRIX MULTIPLICATION (BAND STORAGE MODE)
C
C   USAGE               - CALL VMULBB (A,IA,B,IB,N,C,IC)
C
C   ARGUMENTS    A      - N(1) BY N(1) BAND MATRIX WITH N(2) LOWER
C                           CODIAGINALS AND N(3) UPPER CODIAGONALS.
C                           A IS STORED IN BAND STORAGE MODE AND
C                           THEREFORE HAS DIMENSION N(1) BY
C                           (N(2)+N(3)+1). (INPUT)
C                IA     - ROW DIMENSION OF MATRIX A EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM. (INPUT)
C                B      - N(1) BY N(1) BAND MATRIX WITH N(4) LOWER
C                           CODIAGONALS AND N(5) UPPER CODIAGONALS. B
C                           IS STORED IN BAND STORAGE MODE AND
C                           THEREFORE HAS DIMENSION N(1) BY
C                           (N(4)+N(5)+1). (INPUT)
C                IB     - ROW DIMENSION OF MATRIX B EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM. (INPUT)
C                N      - INPUT VECTOR OF LENGTH 5 CONTAINING THE
C                           FOLLOWING:
C                           N(1) CONTAINS THE NUMBER OF ROWS AND
C                             COLUMNS IN A,B,AND C.
C                           N(2) CONTAINS THE NUMBER OF LOWER
C                             CODIAGONALS IN MATRIX A.
C                           N(3) CONTAINS THE NUMBER OF UPPER
C                             CODIAGONALS IN A.
C                           N(4) CONTAINS THE NUMBER OF LOWER
C                             CODIAGONALS IN B.
C                           N(5) CONTAINS THE NUMBER OF UPPER
C                             CODIAGONALS IN B.
C                C      - OUTPUT MATRIX CONTAINING THE PRODUCT C=A*B.
C                           C IS A BAND MATRIX WITH M LOWER CODIAGONALS,
C                           WHERE M=MIN(N(1)-1,N(2)+N(4)) AND L UPPER
C                           CODIAGONALS, WHERE L=MIN(N(1)-1,N(3)+N(5)).
C                           C IS STORED IN BAND STORAGE MODE AND
C                           THEREFORE HAS DIMENSION N(1) BY (L+M+1).
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
      SUBROUTINE VMULBB (A,IA,B,IB,N,C,IC)
C
      INTEGER            N(5)
      REAL               A(IA,1),B(IB,1),C(IC,1),ZERO
      DOUBLE PRECISION   SUM
      DATA               ZERO/0.0/
C                                  FIRST EXECUTABLE STATEMENT
      NN = N(1)
      NLA = N(2)
      NLAP2 = NLA+2
      NUA = N(3)
      NLB = N(4)
      NUB = N(5)
      NCA = NLA+NUA+1
      NCB = NUB+NLB+1
      NMNUC = NN-NUA
      NLC = MIN0(NLB+NLA,NN-1)
      NLCP1 =NLC+1
      NUC = MIN0(NUB+NUA,NN-1)
      MM = NUC+NLC+1
      N1 = NN+1
C                                  ZERO OUT CORNERS OF THE
C                                  RESULTING MATRIX
      IF (NLC .LE. 0) GO TO 15
      NX = NLC
      DO 10 I = 1,NLC
         DO 5 J = 1,NX
            C(I,J) = ZERO
    5    CONTINUE
         NX = NX-1
   10 CONTINUE
   15 IF (NUC .LE. 0) GO TO 30
      NX = NLC+2
      DO 25 I = 1,NUC
         DO 20 J = NX,MM
            C(N1-I,J) = ZERO
   20    CONTINUE
         NX = NX+1
   25 CONTINUE
   30 NCB = NUB+NLB+1
      JBBEG = NLC
      JCBEG = NLC+1
      DO 50 IRA = 1,NN
         DO 45 JB = 1,NN
            JABEG = NLAP2-IRA
            IRB = 1
            JBBEG = NLB+JB
            KL = MAX0(IDIM(JBBEG,NCB),IDIM(1,JABEG))
            JABEG = JABEG+KL
            JBBEG = JBBEG-KL
            IRB = IRB+KL
            IF (JABEG .GT. NCA) GO TO 45
            IF (IRB .LT. 1) GO TO 45
            IF (JBBEG .LT. 1) GO TO 45
            SUM = 0.0D0
            DO 35 L = 1,NN
               SUM = SUM+DBLE(A(IRA,JABEG))*DBLE(B(IRB,JBBEG))
               JABEG = JABEG+1
               IRB = IRB+1
               JBBEG = JBBEG-1
               IF (JABEG .GT. NCA) GO TO 40
               IF (JBBEG .LT. 1) GO TO 40
               IF (IRB .GT. NN) GO TO 40
   35       CONTINUE
   40       C(IRA,JCBEG) = SUM
            JCBEG = JCBEG+1
   45    CONTINUE
         JCBEG = MAX0(NLCP1-IRA,1)
         IF (IRA .GT. NMNUC) NCA = NCA-1
   50 CONTINUE
      RETURN
      END
