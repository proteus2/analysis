C   IMSL ROUTINE NAME   - MDHYP
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - HYPERGEOMETRIC PROBABILITY DISTRIBUTION
C                           FUNCTION
C
C   USAGE               - CALL MDHYP (K,N,L,ND,PEQK,PLEK,IER)
C
C   ARGUMENTS    K      - INPUT NUMBER OF DEFECTIVES FOR WHICH PROBA-
C                           BILITIES ARE DESIRED
C                N      - INPUT SAMPLE SIZE
C                L      - INPUT LOT SIZE
C                ND     - INPUT NUMBER OF DEFECTIVES IN LOT L
C                PEQK   - OUTPUT PROBABILITY OF OBTAINING K DEFECTIVES
C                           IN THE N-SIZED SAMPLE (DOUBLE PRECISION)
C                PLEK   - OUTPUT PROBABILITY OF OBTAINING K OR LESS
C                           DEFECTIVES (DOUBLE PRECISION)
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER = 129 INDICATES K IS NOT IN THE INCLU-
C                             SIVE RANGE (0,N)
C                           IER = 130 INDICATES N IS NOT IN THE INCLU-
C                             SIVE RANGE (1,L)
C                           IER = 131 INDICATES ND IS NOT IN THE INCLU-
C                             SIVE RANGE (0,L)
C                           IER = 132 INDICATES N IS GREATER THAN 32767
C
C   REQD. IMSL ROUTINES - UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS  1.  MDHYP COMPUTES THE PROBABILITY THAT A RANDOM VARIABLE
C                FOLLOWING THE HYPERGEOMETRIC PROBABILITY DISTRIBUTION
C                IS LESS THAN OR EQUAL TO K. THE RANDOM VARIABLE REPRE-
C                SENTS THE NUMBER OF DEFECTIVES IN A SAMPLE OF SIZE N
C                DRAWN WITHOUT REPLACEMENT FROM A LOT OF SIZE L CON-
C                TAINING ND DEFECTIVES.
C            2.  PEQK AND PLEK ARE DOUBLE PRECISION VARIABLES. THE USER
C                MUST SPECIFY THE CORRESPONDING ACTUAL ARGUMENTS AS SUCH
C                IN THE CALLING PROGRAM.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE MDHYP  (K,N,L,ND,PEQK,PLEK,IER)
C
      DOUBLE PRECISION   EPS,P,SP,U,A,B,A1,B1,AA,BB,ANND,ANXD,AJ
      DOUBLE PRECISION   PEQK,PLEK
      DATA               EPS/1.0D-78/
C                                  CHECK RANGE OF K
C                                  FIRST EXECUTABLE STATEMENT
      IF ((K .GE. 0) .AND. (K .LE. N)) GO TO 5
      IER = 129
      GO TO 9000
C                                  CHECK RANGE OF N
    5 IF ((N .GE. 1) .AND. (N .LE. L)) GO TO 10
      IER = 130
      GO TO 9000
C                                  CHECK RANGE OF ND
   10 IF ((ND .GE. 0) .AND. (ND .LE. L)) GO TO 15
      IER = 131
      GO TO 9000
   15 IF (N .LE. 32767) GO TO 20
      IER = 132
      GO TO 9000
C                                  INITIALIZATION FOR RECURSIVE FORMULA
   20 IER = 0
      P = 0.0D0
      SP = 0.0D0
      IF (K .GT. ND) GO TO 65
      IF ((N-K) .GT. (L-ND)) GO TO 70
      P = 1.0D0
      NND = MIN0(N,ND)
      NXD = MAX0(N,ND)
      B = L
      ANND = NND
      ANXD = NXD
C                                  CHOOSE MORE EFFICIENT DIRECTION FOR
C                                  COMPUTATION OF PROBABILITIES
      IF ((K*(L+2)) .GT. ((ND+1)*(N+1))) GO TO 25
C                                  INITIALIZATION FOR CALCULATION OF
C                                  P(X .LE. K)
      IFLAG = 0
      AA = ANXD - ANND
      BB = B - ANXD - ANND
      K0 = N +(ND - L)
      IF (K0 .LT. 0) K0=0
      AJ = K0
      A1 = ANND-AJ
      B1 = AJ + 1.0D0
      MM1 = K - K0
      A = B - ANXD + AJ
      MM = NND - K0
      GO TO 30
C                                  INITIALIZATION FOR CALCULATION OF
C                                  P(X .GT. K)
   25 IFLAG = 1
      AA = B - ANXD - ANND
      BB = ANXD - ANND
      A1 = ANND
      B1 = 1.0D0
      MM1 = NND - K
      MM = L - NXD
      A = B - ANND
      IF (NND .GE. MM) GO TO 30
      MM = NND
      A = ANXD
   30 ICNT = 0
      IF (MM .EQ. 0) GO TO 45
C                                  LOOP FOR PEQK
      DO 40 M=1,MM
         U = A/B
C                                  DETECT AND CORRECT FOR POTENTIAL
C                                  UNDERFLOW
         IF (U .GE. EPS/P) GO TO 35
         P = P/EPS
         ICNT = ICNT + 1
   35    P = P*U
         A = A - 1.0D0
         B = B - 1.0D0
   40 CONTINUE
   45 IF (MM1 .EQ. 0) GO TO 60
C                                  LOOP FOR PLEK
      DO 55 M=1,MM1
         IF (ICNT .EQ. 0) SP=SP+P
         U = A1*(A1+AA)/(B1*(B1+BB))
         P = U*P
         IF (P .LT. 1.0D0) GO TO 50
C                                  INVERSE UNDERFLOW ADJUSTMENT
         P = P*EPS
         ICNT = ICNT - 1
   50    A1 = A1 - 1.0D0
         B1 = B1 + 1.0D0
   55 CONTINUE
   60 IF (ICNT .NE. 0) P=0.0D0
      IF (IFLAG .NE. 0) GO TO 65
C                                  P(X .LE. K)=P(X .EQ. K)+P(X .LT. K)
      SP = SP + P
      GO TO 70
   65 SP = 1.0D0 - SP
   70 PEQK = P
      PLEK = SP
      GO TO 9005
 9000 CONTINUE
      CALL UERTST (IER,6HMDHYP )
 9005 RETURN
      END
