C   IMSL ROUTINE NAME   - GGNSM
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1981
C
C   PURPOSE             - MULTIVARIATE NORMAL RANDOM DEVIATE GENERATOR
C                           WITH GIVEN COVARIANCE MATRIX
C
C   USAGE               - CALL GGNSM (DSEED,NR,K,SIGMA,IR,RVEC,WKVEC,IER
C
C   ARGUMENTS    DSEED  - INPUT/OUTPUT DOUBLE PRECISION VARIABLE
C                           ASSIGNED AN INTEGER VALUE IN THE
C                           EXCLUSIVE RANGE (1.D0, 2147483647.D0).
C                           DSEED IS REPLACED BY A NEW VALUE TO BE
C                           USED IN A SUBSEQUENT CALL.
C                NR     - INPUT.  NUMBER OF K-DEVIATE VECTORS TO BE
C                           GENERATED.
C                K      - INPUT.  NUMBER OF RANDOM DEVIATES PER VECTOR.
C                SIGMA  - INPUT/OUTPUT VECTOR OF LENGTH K(K+1)/2. ON
C                           INPUT SIGMA CONTAINS THE VARIANCE-COVARIANCE
C                           VALUES. SIGMA IS A POSITIVE DEFINITE MATRIX
C                           STORED IN SYMMETRIC STORAGE MODE. AFTER THE
C                           FIRST CALL TO GGNSM, SIGMA IS REPLACED BY
C                           ITS FACTOR (SQUARE ROOT) ON OUTPUT. (SEE
C                           REMARKS)
C                IR     - INPUT. ROW DIMENSION OF MATRIX RVEC EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM.
C                RVEC   - OUTPUT.  NR BY K MATRIX OF MULTIVARIATE NORMAL
C                           DEVIATES.
C                WKVEC  - INPUT WORK VECTOR OF LENGTH K. WKVEC(1)
C                           SHOULD BE SET TO 0.0 ON THE FIRST OF A
C                           SERIES OF CALLS TO GGNSM. FOR ALL SUB-
C                           SEQUENT CALLS, WKVEC(1) SHOULD BE NONZERO.
C                           IF ONLY ONE CALL IS REQUIRED, SET WKVEC(1)
C                           TO 0.0. (SEE REMARKS) THE REMAINDER OF
C                           WKVEC IS USED AS WORK AREA FOR THE NORMAL
C                           DEVIATE GENERATION.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                           TERMINAL ERROR
C                             IER = 129 INDICATES THAT INPUT MATRIX
C                               SIGMA IS ALGORITHMICALLY NOT POSITIVE
C                               DEFINITE.
C
C   REQD. IMSL ROUTINES - GGNML,GGUBS,MDNRIS,MERFI,UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS      IF THE USER WISHES TO CONTINUE GENERATING MULTIVARIATE
C                NORMAL DEVIATE VECTORS DISTRIBUTED WITH THE SAME SIGMA,
C                THEN MULTIPLE CALLS MAY BE MADE TO GGNSM WITH WKVEC(1)
C                NONZERO ON INPUT. WKVEC(1) SET TO 0.0 ON INPUT TRIGGERS
C                THE CALCULATION OF THE FACTOR (SQUARE ROOT) OF SIGMA.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE GGNSM   (DSEED,NR,K,SIGMA,IR,RVEC,WKVEC,IER)
C
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            NR,K,IR,IER
      REAL               SIGMA(1),RVEC(IR,K),WKVEC(K)
      DOUBLE PRECISION   DSEED
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            L,I,II,J,KKK,IQ,IP1,IP,IRR
      REAL               ZERO,ONE,FOUR,SIXTN,SIXTH,X,RN,D1,D2,Q
      DOUBLE PRECISION   TEMP
      DATA               ZERO/0.0E0/,ONE/1.0E0/,FOUR/4.0E0/,
     1                   SIXTN/16.0E0/,SIXTH/.0625E0/
C                                  FIRST EXECUTABLE STATEMENT
      IF (WKVEC(1).NE.0.0) GO TO 65
C                                  DECOMPOSE SIGMA MATRIX
      D1 = ONE
      D2 = ZERO
      RN = ONE/(K*SIXTN)
      IP = 1
      IER = 0
      DO 45 I=1,K
         IQ = IP
         IRR = 1
         DO 40 J=1,I
            X = SIGMA(IP)
            IF (J.EQ.1) GO TO 10
            DO 5 KKK=IQ,IP1
               X = X-SIGMA(KKK)*SIGMA(IRR)
               IRR = IRR+1
    5       CONTINUE
   10       IF (I.NE.J) GO TO 30
            D1 = D1*X
            Q = SIGMA(IP)+X*RN
            IF (Q.LE.SIGMA(IP)) GO TO 50
   15       IF (ABS(D1).LE.ONE) GO TO 20
            D1 = D1*SIXTH
            D2 = D2+FOUR
            GO TO 15
   20       IF (ABS(D1).GE.SIXTH) GO TO 25
            D1 = D1*SIXTN
            D2 = D2-FOUR
            GO TO 20
   25       SIGMA(IP) = ONE/SQRT(X)
            GO TO 35
   30       SIGMA(IP) = X*SIGMA(IRR)
   35       IP1 = IP
            IP = IP+1
            IRR = IRR+1
   40    CONTINUE
   45 CONTINUE
      GO TO 55
   50 IER = 129
      GO TO 9000
C                                  RECALCULATE DIAGONAL OF SIGMA
   55 L = 0
      DO 60 I=1,K
         L = L+I
   60 SIGMA(L) = 1.0/SIGMA(L)
C                                  GENERATE NR X K NORMAL RANDOM DEVIATE
   65 DO 80 J=1,NR
C                                  GENERATE K NORMAL DEVIATES
         CALL GGNML (DSEED,K,WKVEC)
         L = 1
C                                  CONVERT THE K UNIVARIATE NORMAL
C                                    RANDOM DEVIATES TO MULTIVARIATE
C                                    NORMAL DEVIATES.
         DO 75 II=1,K
            TEMP = 0.D0
            DO 70 I=1,II
               TEMP = TEMP + DBLE(WKVEC(I))*DBLE(SIGMA(L))
   70       L = L+1
            RVEC(J,II) = TEMP
   75    CONTINUE
   80 CONTINUE
      GO TO 9005
 9000 CONTINUE
      CALL UERTST (IER,'GGNSM ')
 9005 RETURN
      END
