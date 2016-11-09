C   IMSL ROUTINE NAME   - GGSTA
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1982
C
C   PURPOSE             - STABLE DISTRIBUTION RANDOM DEVIATE GENERATOR
C
C   USAGE               - CALL GGSTA (DSEED,ALPHA,BPRIM,NR,R)
C
C   ARGUMENTS    DSEED  - INPUT/OUTPUT DOUBLE PRECISION VARIABLE
C                           ASSIGNED AN INTEGER VALUE IN THE
C                           EXCLUSIVE RANGE (1.D0, 2147483647.D0).
C                           DSEED IS REPLACED BY A NEW VALUE TO BE
C                           USED IN A SUBSEQUENT CALL.
C                ALPHA  - INPUT CHARACTERISTIC EXPONENT, GREATER THAN
C                           ZERO, AND LESS THAN OR EQUAL TO 2.
C                BPRIM  - INPUT SKEWNESS PARAMETER IN REVISED
C                           PARAMETERIZATION.
C                NR     - INPUT NUMBER OF STABLE RANDOM VARIATES
C                           TO BE GENERATED.
C                R      - OUTPUT VECTOR OF LENGTH NR CONTAINING THE
C                           STABLE RANDOM DEVIATES GENERATED.
C
C   REQD. IMSL ROUTINES - GGEXN,GGSTA1,GGUBFS,GGUBS
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE GGSTA(DSEED,ALPHA,BPRIM,NR,R)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            NR
      REAL               ALPHA,BPRIM,R(NR)
      DOUBLE PRECISION   DSEED
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I
      REAL               AN1,PT5,ONE,TWO,THR1,W(1)
      DOUBLE PRECISION   DA,DB,PIBY2,PHIBY2,EPH,EPI,EPS
      DOUBLE PRECISION   DONE,DALPHA,D2S,GGSTA1
      DOUBLE PRECISION   P1,P2,Q1,Q2,Q3,PV
      DATA               P1/.8400668525364832D3/
      DATA               P2/.2000111415899646D2/
      DATA               Q1/.1680133705079266D4/
      DATA               Q2/.1801337040739002D3/
      DATA               Q3/1.D0/
      DATA               PIBY2/1.570796326794897D0/
      DATA               THR1/0.99/,DONE/1.0D0/
      DATA               AN1/-.99/,PT5/.5/,ONE/1.0/,TWO/2.0/
C                                  FIRST EXECUTABLE STATEMENT
      DALPHA = DBLE(ALPHA)
      EPS = DONE-DALPHA
      EPI = EPS*PIBY2
C                                  GENERATE VARIATES
      DO 25  I=1,NR
         U = GGUBFS(DSEED)
         PHIBY2 = PIBY2*(U-PT5)
         A = PHIBY2 * GGSTA1(PHIBY2)
         EPH = EPS*PHIBY2
         BB = GGSTA1(EPH)
         B = EPH*BB
         IF (EPS.GT.AN1) TAU = BPRIM/(GGSTA1(EPI)*PIBY2)
         IF (EPS.LE.AN1) TAU = BPRIM*EPI*DALPHA*GGSTA1(DALPHA*PIBY2)
C                                  COMPUTE NEEDED SUBEXPRESSIONS
C                                   IF PHI NEAR PI BY 2,USE DBL.PREC.
         IF (A.GT.THR1) GO TO 5
C                                   SINGLE PRECISION
         A2 = A**2
         A2P = ONE+A2
         A2 = ONE-A2
         B2 = B**2
         B2P = ONE+B2
         B2 = ONE-B2
         GO TO 10
C                                   DOUBLE PRECISION
    5    DA = DBLE(A)**2
         DB = DBLE(B)**2
         A2 = DONE-DA
         A2P = DONE+DA
         B2 = DONE-DB
         B2P = DONE+DB
C                                   COMPUTE COEFFICIENT
   10    CALL GGEXN(DSEED,ONE,1,W)
         Z = A2P*(B2+TWO*PHIBY2*BB*TAU)/(W(1)*A2*B2P)
C                                   COMPUTE EXPONENTIAL-TYPE EXPRESSION
         ALOGZ = ALOG(Z)
         D1 = EPS*ALOGZ/(ONE-EPS)
         IF(ABS(D1) .GT. 0.1) GO TO 15
         D2S=D1*D1
         PV=P1+D2S*P2
         D2=2.D0*PV/(Q1+D2S*(Q2+D2S*Q3)-D1*PV)
         GO TO 20
   15    D2=(EXP(D1)-1.0)/D1
   20    D = D2*(ALOGZ/(ONE-EPS))
C                                   COMPUTE STABLE
         R(I) = (ONE+EPS*D)*TWO
     1        * ((A-B)*(ONE+A*B) - PHIBY2*TAU*BB*(B*A2-TWO*A))
     2        / (A2*B2P) + TAU*D
   25 CONTINUE
      RETURN
      END
