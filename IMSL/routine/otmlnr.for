C   IMSL ROUTINE NAME   - OTMLNR
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - MAXIMUM LIKELIHOOD ESTIMATION FROM GROUPED
C                           AND/OR CENSORED NORMAL DATA
C
C   USAGE               - CALL OTMLNR (X,IXI,N,IP,RM,SIGMA,E1,E2,
C                           MAXITS,IDS,COV,NOBS,K,IER)
C   ARGUMENTS    X      - INPUT N BY 2 MATRIX CONTAINING THE OBSER-
C                           VATIONS. X MAY BE DESTROYED ON OUTPUT.
C                IXI    - INPUT ROW DIMENSION OF MATRIX X EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN
C                           THE CALLING PROGRAM.
C                N      - INPUT NUMBER OF OBSERVATIONS
C                IP     - INPUT VECTOR OF LENGTH N CONTAINING CODES
C                           DESCRIBING THE OBSERVATIONS IN X (SEE
C                           REMARKS).  IP MAY BE DESTROYED ON OUTPUT.
C                RM     - ON INPUT, OPTIONALLY CONTAINS THE STARTING
C                           VALUE FOR THE MEAN.
C                         ON EXIT, RM CONTAINS THE CALCULATED MAXIMUM
C                           LIKELIHOOD ESTIMATE FOR THE MEAN.
C                SIGMA  - ON INPUT, SIGMA CONTAINS THE STARTING VALUE
C                           FOR THE STANDARD DEVIATION. IF SIGMA IS
C                           NEGATIVE, THE PROGRAM OTMLNR WILL PROVIDE
C                           THE STARTING VALUE FOR THE MEAN (RM) AND THE
C                           STANDARD DEVIATION (SIGMA).
C                         ON EXIT, SIGMA CONTAINS THE CALCULATED
C                           MAXIMUM LIKELIHOOD ESTIMATE FOR THE
C                           STANDARD DEVIATION.
C                E1     - INPUT PARAMETER CONTAINING THE CONVERGENCE
C                           CRITERION FOR THE MEAN.
C                E2     - INPUT PARAMETER CONTAINING THE CONVERGENCE
C                           CRITERION FOR THE STANDARD DEVIATION. E1
C                           AND E2 SHOULD NOT BE SMALLER THAN MACHINE
C                           PRECISION.
C                MAXITS - INPUT PARAMETER CONTAINING THE MAXIMUM NUMBER
C                           OF ITERATIONS ALLOWED. MAXITS COULD BE SET
C                           AT 25.
C                IDS    - INPUT CONVERGENCE PARAMETER. IF IDS CONSECU-
C                           TIVE VALUES FOR DELTA SIGMA INCREASE, THE
C                           PROCESS IS ASSUMED TO BE DIVERGENT AND THE
C                           PROCESS IS TERMINATED WITH IER = 130. A
C                           SUGGESTED VALUE FOR IDS IS 3.
C                COV    - OUTPUT MATRIX OF ORDER 2 CONTAINING AN
C                           ESTIMATE OF THE COVARIANCE MATRIX OF THE
C                           ESTIMATES. IF NO CONVERGENCE WAS OBTAINED,
C                           (SEE PARAMETER - IER) COV(1,2) AND COV(2,1)
C                           CONTAIN THE LAST CORRECTIONS TO THE MEAN
C                           AND STANDARD DEVIATION ESTIMATES.
C                NOBS   - OUTPUT VECTOR CONTAINING THE NUMBERS OF EACH
C                           TYPE OF OBSERVATION
C                           NOBS(1) -- THOSE SPECIFIED BY EXACT
C                             OBSERVATIONS
C                           NOBS(2) -- THOSE SPECIFIED BY LOWER BOUND
C                           NOBS(3) -- THOSE SPECIFIED BY UPPER BOUND
C                           NOBS(4) -- THOSE SPECIFIED BY TWO BOUNDS
C                K      - OUTPUT PARAMETER CONTAINING THE NUMBER OF
C                           ITERATIONS TAKEN
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR        .
C                           N = 129 INDICATES THAT AN ERROR OCCURRED IN
C                             SUBROUTINE MSMRAT.
C                           N = 130 INDICATES THAT THE PROCESS FAILED
C                             TO CONVERGE.
C                           N = 131 INDICATES THAT THE NUMBER OF
C                             OBSERVATIONS WAS LESS THAN 2.
C                           N = 132 INDICATES THAT THE UPPER AND LOWER
C                             BOUNDS ON AN OBSERVATION, STANDARDIZED
C                             BY SUBTRACTION OF AND DIVISION BY THE
C                             CURRENT ESTIMATES OF THE MEAN AND
C                             STANDARD DEVIATION, BECAME EQUAL DURING
C                             THE ITERATION.
C                           N = 133 MEANS THAT ONLY LOWER BOUND
C                             SPECIFICATIONS WERE PROVIDED. THE PROCESS
C                             FAILS.
C                           N = 134 MEANS THAT ONLY UPPER BOUND
C                             SPECIFICATIONS WERE PROVIDED. THE PROCESS
C                             FAILS.
C
C   REQD. IMSL ROUTINES - MERRC=ERFC,MSMRAT,UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS      THE INPUT DATA FOR OBSERVATION I, CONTAINED IN MATRIX
C                X AND VECTOR IP, MUST BE SPECIFIED AS FOLLOWS.
C
C                  IF THE OBSERVATION IS EXACT, X(I,1) CONTAINS THE
C                  OBSERVATION AND IP(I) = 0.
C
C                  IF THE OBSERVATION IS KNOWN BY A LOWER BOUND, X(I,1)
C                  CONTAINS THE LOWER BOUND AND IP(I) = 1.
C
C                  IF THE OBSERVATION IS KNOWN BY AN UPPER BOUND,
C                  X(I,1) CONTAINS THE UPPER BOUND AND IP(I) = -1.
C
C                  IF THE OBSERVATION IS KNOWN BY AN UPPER AND LOWER
C                  BOUND, X(I,1) CONTAINS THE LOWER BOUND, X(I,2)
C                  CONTAINS THE UPPER BOUND, AND IP(I) = 2.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE OTMLNR (X,IXI,N,IP,RM,SIGMA,E1,E2,MAXITS,IDS,COV,NOBS,
     1                   K,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            IXI,N,IP(1),MAXITS,IDS,NOBS(4),K,IER
      REAL               X(IXI,2),COV(2,2),RM,SIGMA,E1,E2
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            IX,I2,IL,IU,I,M
      DOUBLE PRECISION   S,T,S1,T1,W,Y,SIG,XMEAN,A,RSIG,DS1,B,C,D,E,F,G,
     1                   RK,A0,RL11,RL22,RL12,Y1,A1,DS,ASIG,BSIG,AMU,DM,
     2                   RL1,RL2,XX
C                                  FIRST EXECUTABLE STATEMENT
      K = 0
      IER = 0
      IF (N.GE.2) GO TO 5
C                                  TERMINAL ERROR - LESS THAN TWO
C                                  OBSERVATIONS WERE SPECIFIED
      IER = 131
      GO TO 9000
    5 IX = 0
      I2 = 0
      IL = 0
      IU = 0
      S = 0.0D0
      S1 = 0.0D0
      T = 0.D0
      T1 = 0.0D0
      DO 25 I=1,N
         IF (IP(I).EQ.0) GO TO 20
         IF (IP(I).NE.2) GO TO 10
         IF (X(I,1).EQ.X(I,2)) GO TO 20
C                                  ELIMINATE EXACTLY SPECIFIED ELEMENTS
C                                  AND CALCULATE STARTING VALUES FOR
C                                  THE MEAN AND SIGMA
   10    IU = IU+1
         IP(IU) = IP(I)
         X(IU,1) = X(I,1)
         IF (IP(I).NE.2) GO TO 15
C                                  ELEMENTS SPECIFIED BY TWO BOUNDS ARE
C                                  COUNTED
         I2 = I2+1
         IF (IX.EQ.0.AND.I2.EQ.1) W = .5D0*(X(I,1)+X(I,2))
         X(IU,2) = X(I,2)
         Y = .5D0*(X(I,1)+X(I,2))-W
         T1 = T1+Y
         S1 = S1+Y*Y
         GO TO 25
C                                  ELEMENTS SPECIFIED BY LOWER BOUND
C                                  ARE COUNTED
   15    IF (IP(I).EQ.1) IL = IL+1
         GO TO 25
C                                  EXACTLY SPECIFIED ELEMENTS ARE
C                                  COUNTED
   20    IX = IX+1
         IF (I2.EQ.0.AND.IX.EQ.1) W = X(I,1)
         Y = X(I,1)-W
         T = T+Y
         S = S+Y*Y
   25 CONTINUE
C                                  SAVE COUNTS
      M = IX+I2
      NOBS(1) = IX
      XX = IX
      NOBS(2) = IL
      NOBS(3) = N-IL-M
      NOBS(4) = I2
      IF (SIGMA.GT.0..AND.IX.NE.N) GO TO 35
      IF (M.GT.1) GO TO 30
C                                  TOO FEW EXACTLY SPECIFIED ELEMENTS TO
C                                  CALCULATE STARTING VALUES
      XMEAN = 1.0D0
      SIG = 1.0D0
      GO TO 40
C                                  CALCULATE STARTING VALUES
   30 T1 = T1+T
      S1 = S1+S
      A = T1/M
      XMEAN = A+W
      SIG = DSQRT((S1-T1*A)/(M-1))
      IF (IX.NE.N) GO TO 40
C                                  ALL ELEMENTS WERE EXACTLY SPECIFIED
      RM = XMEAN
      SIGMA = SIG
      COV(1,1) = SIG*SIG/N
      COV(2,2) = .5*COV(1,1)
      COV(1,2) = 0.0
      COV(2,1) = 0.0
      GO TO 9005
C                                  STARTING VALUES WERE USER SUPPLIED
   35 SIG = SIGMA
      XMEAN = RM
   40 A = T+XX*W
      S = S+W*(T+A)
      T = A
      M = N-IX
      I2 = -1
      DS1 = 0.0D0
      IF (SIG.EQ.0.0D0) SIG = 1.D0
   45 RSIG = 1.D0/SIG
      A = 0.0D0
      B = 0.0D0
      C = 0.0D0
      D = 0.0D0
      E = 0.0D0
      F = 0.0D0
      G = 0.0D0
      K = K+1
      DO 70 I=1,M
         Y = X(I,1)
         U = (Y-XMEAN)*RSIG
         IF (IP(I).EQ.-1) GO TO 50
         CALL MSMRAT (U,TR,IER)
         GO TO 55
   50    CALL MSMRAT (-U,TR,IER)
         IF (IER.GT.128) GO TO 9000
         TR = -TR
   55    IF (IP(I).NE.2) GO TO 65
         Y1 = X(I,2)
         U1 = (Y1-XMEAN)*RSIG
         IF (U.NE.U1) GO TO 60
         IER = 132
         GO TO 9000
   60    CALL MSMRAT (U1,TR1,IER)
         IF (IER.GE.128) GO TO 9000
         W = DEXP(.5D0*(DBLE(U1)*DBLE(U1)-DBLE(U)*DBLE(U)))
         RK = DBLE(TR)*DBLE(TR1)/(TR1*W-TR)
         A0 = RK*(W-1.D0)
         A = A+A0
         RL11 = Y*W
         RL22 = Y*Y
         RL12 = Y1*Y1
         A1 = RK*(RL11-Y1)
         B = B+A1
         C = C+RK*(RL22*W-RL12)
         D = D+RK*(RL22*RL11-Y1*RL12)
         E = E+A0*A0
         F = F+A0*A1
         G = G+A1*A1
         GO TO 70
   65    W = Y*Y
         A1 = Y*TR
         A = A+TR
         B = B+A1
         C = C+W*TR
         D = D+W*A1
         E = E+DBLE(TR)*DBLE(TR)
         F = F+A1*TR
         G = G+A1*A1
   70 CONTINUE
C                                  CALCULATE FIRST AND SECOND
C                                  DERIVATIVES OF THE LIKELIHOOD
C                                  FUNCTION
      A0 = T-XX*XMEAN
      A1 = B-XMEAN*A
      DS = RSIG*RSIG
      RL1 = A0*DS+A*RSIG
      RL2 = (((S-XMEAN*(A0+T))*RSIG+A1)*RSIG-XX)*RSIG
      AMU = A1*DS-E*RSIG
      ASIG = (((-A1-B)*XMEAN+C)*RSIG+XMEAN*E-F)*DS
      BSIG = (((B*XMEAN-C-C)*XMEAN+D)*RSIG+XMEAN*F-G)*DS
      RL11 = (AMU-XX*RSIG)*RSIG
      RL12 = RSIG*(ASIG-RSIG*(A+2.D0*RSIG*A0))
      RL22 = DS*(XX+BSIG-XMEAN*ASIG-RSIG*(A1+A1)+3.D0*RSIG*(XMEAN*(-A0
     1-T)+S))
C                                  FIND CORRECTIONS BY SOLVING NEWTON-
C                                  RAPHSON EQUATIONS
      A = RL12*RL12-RL11*RL22
      DS = (RL2*RL11-RL1*RL12)/A
      DM = -(RL1+DS*RL12)/RL11
C                                  FIND NEW APPROXIMATIONS TO MAXIMUM
C                                  LIKELIHOOD ESTIMATES
      XMEAN = XMEAN+DM
      SIG = SIG+DS
      RL1 = DABS(DS)
      IF (RL1.GT.DS1) GO TO 75
      I2 = 0
      GO TO 80
   75 I2 = I2+1
   80 DS1 = RL1
C                                  IF SIG WOULD BECOME NEGATIVE, KEEP
C                                  SIG POSITIVE AND MOVE IT HALF WAY TO
C                                  ZERO
      IF (SIG.LE.0.D0) SIG = .5D0*(SIG-DS)
      IF (DABS(DM).LE.E1.OR.RL1.LE.E2) GO TO 85
      IF (I2.LT.IDS.AND.K.LT.MAXITS) GO TO 45
C                                  TERMINAL ERROR - THE PROCESS IS
C                                  ASSUMED TO BE DIVERGING SINCE THREE
C                                  CONSECUTIVE VALUES OF DS WERE
C                                  INCREASING
      IER = 130
C                                  ALL LOWER BOUNDS SPECIFIED
      IF (N.EQ.NOBS(2)) IER = 133
C                                  ALL UPPER BOUNDS SPECIFIED
      IF (N.EQ.NOBS(3)) IER = 134
   85 RM = XMEAN
      SIGMA = SIG
      COV(1,1) = RL22/A
      COV(2,2) = RL11/A
      IF (IER.EQ.0) GO TO 90
      COV(1,2) = DM
      COV(2,1) = DS
      GO TO 9000
   90 COV(1,2) = -RL12/A
      COV(2,1) = COV(1,2)
      GO TO 9005
 9000 CONTINUE
      CALL UERTST (IER,6HOTMLNR)
 9005 RETURN
      END
