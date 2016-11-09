C   IMSL ROUTINE NAME   - NRBHA
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - BHAPKAR V TEST
C
C   USAGE               - CALL NRBHA (X,NPS,IC,IR,W,V,Q,IER)
C
C   ARGUMENTS    X      - INPUT VECTOR CONTAINING THE OBSERVATIONS.
C                           X IS A VECTOR OF LENGTH
C                           NPS(1)+NPS(2)+...+NPS(IC).
C                           X(1),...,X(NPS(1)) CONTAINS THE SAMPLE 1
C                           OBSERVATIONS.  THE OBSERVATIONS FOR
C                           SAMPLE I FOLLOW THOSE FOR SAMPLE I-1 IN THE
C                           X VECTOR, FOR I=2,...,IC.
C                         ON OUTPUT, X IS DESTROYED.
C                NPS    - INPUT VECTOR OF LENGTH IC.
C                           NPS(I) CONTAINS THE NUMBER OF OBSERVATIONS
C                           IN SAMPLE I. THE ELEMENTS OF NPS MUST BE
C                           POSITIVE.
C                IC     - INPUT NUMBER OF SAMPLES.
C                IR     - WORK VECTOR OF LENGTH
C                           NPS(1)+NPS(2)+...+NPS(IC)+3.
C                W      - WORK VECTOR OF LENGTH IC.
C                V      - OUTPUT BHAPKAR STATISTIC.
C                Q      - OUTPUT PROBABILITY OF EXCEEDING V IF THE
C                           HYPOTHESIS OF EQUALITY IS TRUE.
C                           IT IS ASSUMED THAT V IS DISTRIBUTED AS
C                           A CHI-SQUARE VARIATE WITH IC-1 DEGREES OF
C                           FREEDOM.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER=129 INDICATES THAT AN ERROR OCCURRED IN
C                             IMSL ROUTINE MDCH.
C                           IER=130 INDICATES THAT SOME ELEMENT OF NPS
C                             WAS LESS THAN OR EQUAL TO ZERO.
C
C   REQD. IMSL ROUTINES - H32/MDCH,MDNOR,MGAMAD=DGAMMA,UERTST,UGETIO,
C                           VSRTR
C                       - H36,H48,H60/MDCH,MDNOR,MGAMA=GAMMA,
C                           UERTST,UGETIO,VSRTR
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE NRBHA  (X,NPS,IC,IR,W,V,Q,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            NPS(1),IC,IR(1),IER
      REAL               X(1),W(1),V,Q
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,IBEG,IEND,IEND1,IRJ,IRR,I1,J,K,M,MM
      REAL               CI,DQ,FM,P,PM,PN,PS,S,S1,Y,Y2
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
      IEND = 0
      M = 1
C                                  ASSIGN A SAMPLE NUMBER TO EACH
C                                  ELEMENT OF THE X VECTOR
      DO 10 I = 1,IC
         IF (NPS(I) .GT. 0) GO TO 5
C                                  TERMINAL ERROR - SOME SAMPLE SIZE IS
C                                  NON-POSITIVE
         IER = 130
         GO TO 9000
    5    IBEG = IEND+1
         IEND = IEND+NPS(I)
         W(I) = 0.
         M = M*NPS(I)
         DO 10 J = IBEG,IEND
   10    IR(J) = I
      DO 15 I = 1,IC
   15 IR(IEND+I) = NPS(I)
C                                  SORT ELEMENTS OF THE X VECTOR WITH
C                                  RETENTION OF THE SAMPLE NUMBERS
C                                  ASSOCIATED WITH EACH ELEMENT
      CALL VSRTR (X,IEND,IR)
C                                  CALCULATE THE NUMBER OF IC-PLETS
C                                  WHICH MAY BE FORMED SUCH THAT THE
C                                  OBSERVATION FROM SAMPLE I IS THE
C                                  LEAST
      IEND1 = IEND-IC+1
      DO 45 I = 1,IEND1
         I1 = I+1
         Y = X(I)
         K = I1
         DO 20 J = I1,IEND
            IF (Y .NE. X(J)) GO TO 25
            IF (IR(I) .NE. IR(J)) K = K+1
   20    CONTINUE
   25    DO 30 J = 1,IC
   30    NPS(J) = 0
         DO 35 J = K,IEND
            IRJ=IR(J)
   35    NPS(IRJ)=NPS(IRJ)+1
         IRR = IR(I)
         NPS(IRR) = 1
         MM = 1
         DO 40 J = 1,IC
   40    MM = MM*NPS(J)
   45 W(IRR) = MM + W(IRR)
C                                  CALCULATE BHAPHAR V STATISTIC
      CI = 1./IC
      PM = M
      MM = IC-1
      PS = IEND
      PN = IEND*(IC+MM)
      S = 0.
      S1 = 0.
      DO 50 I = 1,IC
         Y = W(I)/PM-CI
         Y2 = Y*Y
         NPS(I) = IR(IEND+I)
         P = NPS(I)/PS
         S = S+P*Y2
   50    S1 = S1+P*Y
      S1 = S1*S1
      V = PN*(S-S1)
      FM=MM
C                                  CALCULATE PROBABILITY OF THE COMPUTED
C                                  VALUE OR A MORE EXTREME VALUE OF V
      CALL MDCH (V,FM,DQ,IER)
      Q=1.0-DQ
      IF(IER .LE. 34) GO TO 9005
 9000 CONTINUE
      CALL UERTST (IER,6HNRBHA )
 9005 RETURN
      END
