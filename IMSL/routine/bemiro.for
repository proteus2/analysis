C   IMSL ROUTINE NAME   - BEMIRO
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1980
C
C   PURPOSE             - ESTIMATES OF MEANS, SIMPLE REGRESSION
C                           COEFFICIENTS, THEIR INTERCEPTS, STANDARD
C                           ERRORS OF THE REGRESSION COEFFICIENTS, AND
C                           STANDARD DEVIATIONS FOR ARRAYS WHICH
C                           CONTAIN MISSING VALUES. (OUT-OF-CORE
C                           VERSION)
C
C   USAGE               - CALL BEMIRO (X,N,M,IND,XMEAN,B,A,S,IBAS,INCD,
C                           IER)
C
C   ARGUMENTS    X      - INPUT VECTOR OF LENGTH M CONTAINING ONE
C                           OBSERVATION FOR EACH OF M VARIABLES.
C                         ON OUTPUT, X IS EQUAL TO X MINUS THE
C                           TEMPORARY MEAN.
C                N      - INPUT.  NUMBER OF OBSERVATIONS PER VARIABLE.
C                           N MUST BE GREATER THAN 1.
C                M      - INPUT.  NUMBER OF VARIABLES.  M MUST BE
C                           GREATER THAN 1.
C                IND    - INPUT INDICATOR FOR MEAN VALUES.
C                           ON INITIAL ENTRY IND MUST BE LESS THAN OR
C                             EQUAL TO ZERO.
C                           IF IND EQUALS ZERO THE FIRST OBSERVATION
C                             ON THE M VARIABLES WILL BE USED AS A
C                             TEMPORARY MEAN VECTOR.
C                           IF IND IS LESS THAN ZERO THE FIRST COLUMN
C                             OF MATRIX B MUST CONTAIN M TEMPORARY
C                             MEAN VALUES.
C                         OUTPUT.  COUNTER FOR THE NUMBER OF
C                           OBSERVATIONS THAT HAVE BEEN ENTERED.
C                           IND IS SET TO ONE ON INITIAL ENTRY AND
C                           INCREMENTED EACH ADDITIONAL TIME BEMIRO
C                           IS ENTERED.
C                XMEAN  - OUTPUT VECTOR OF LENGTH M CONTAINING MEANS
C                           OF THE OBSERVATIONS ON M VARIABLES.
C                           XMEAN MUST BE A DOUBLE PRECISION VECTOR.
C                B      - INPUT.  WHEN IND IS LESS THAN ZERO ON INITIAL
C                           ENTRY, THE FIRST COLUMN OF B MUST CONTAIN
C                           M TEMPORARY MEAN VALUES.
C                         OUTPUT MATRIX OF DIMENSION M BY M CONTAINING
C                           REGRESSION COEFFICIENTS. B(I,J) CONTAINS
C                           THE COEFFICIENT FOR THE REGRESSION OF
C                           VARIABLE J ON VARIABLE I. THE DIAGONAL OF
C                           B IS NOT USED.
C                A      - OUTPUT MATRIX OF DIMENSION M BY 2M CONTAINING
C                           THE INTERCEPTS IN THE FIRST M COLUMNS.
C                           A(I,J), FOR J+1,2,...,M, CONTAINS THE
C                           INTERCEPT FOR THE REGRESSION OF VARIABLE J
C                           ON VARIABLE I. THE A(I,I) ELEMENTS ARE NOT
C                           USED. THE REMAINING LOCATIONS ARE WORK
C                           STORAGE. A MUST BE A DOUBLE PRECISION
C                           MATRIX.
C                S      - OUTPUT MATRIX OF DIMENSION M BY M CONTAINING
C                           THE STANDARD ERRORS OF THE REGRESSION
C                           COEFFICIENTS AND STANDARD DEVIATIONS.
C                           THE STANDARD DEVIATIONS LIE ON THE DIAGONAL
C                           OF THE MATRIX.
C                IBAS   - INPUT. ROW DIMENSION OF MATRICES B, A, AND S
C                           EXACTLY AS SPECIFIED IN THE DIMENSION
C                           STATEMENT IN THE CALLING PROGRAM.
C                INCD   - OUTPUT VECTOR CONTAINING THE INCIDENCE VALUES.
C                           INDICATES THE NUMBER OF PAIRS OF
C                           OBSERVATIONS USED IN THE CALCULATIONS OF
C                           B, A, AND S.  INCD IS STORED IN SYMMETRIC
C                           MATRIX MODE AND MUST BE OF LENGTH AT LEAST
C                           M(M+1)/2.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         WARNING ERROR
C                           IER=33 INDICATES THAT FEWER THAN 3 PAIRS OF
C                             OBSERVATIONS WERE PRESENT.
C                             THE STANDARD ERRORS ARE SET TO
C                             NEGATIVE MACHINE INFINITY.
C                           IER=34 INDICATES THAT FEWER THAN 2 PAIRS OF
C                             OBSERVATIONS WERE PRESENT.
C                             THE STANDARD DEVIATIONS, STANDARD ERRORS,
C                             REGRESSION COEFFICIENTS, AND INTERCEPTS
C                             ARE SET TO NEGATIVE MACHINE INFINITY.
C                           IER=35 INDICATES THAT THE OBESERVATIONS ON
C                             SOME VARIABLE WERE CONSTANT. THE PERTINENT
C                             REGRESSION STATISTICS ARE SET TO NEGATIVE
C                             INFINITY.
C
C   REQD. IMSL ROUTINES - UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE BEMIRO (X,N,M,IND,XMEAN,B,A,S,IBAS,INCD,IER)
C
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N,M,IND,IBAS,INCD(1),IER
      REAL               X(M),B(IBAS,M),S(IBAS,M)
      DOUBLE PRECISION   XMEAN(M)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      DOUBLE PRECISION   A(IBAS,1),Z,Z1
      DATA               XMISS/-999999./
      DATA               SINFM/ZFFFFFFFF/
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
C                                  TEST FOR FIRST CALL OF A SEQUENCE
      IF (IND.GT.0) GO TO 30
C                                  TEST FOR LOCATION OF TEMP MEANS
      IF (IND.NE.0) GO TO 10
C                                  PUT THE TEMP MEANS IN B
      DO 5 I=1,M
    5 B(I,1) = X(I)
C                                  INITIALIZE ACCUMULATORS
   10 IND = 0
      DO 20 I=1,M
         DO 15 J=1,M
            JM = J+M
            A(I,J) = 0.0D0
   15    A(I,JM) = 0.0D0
   20 XMEAN(I) = 0.0D0
      K = (M*(M+1))/2
      DO 25 I=1,K
   25 INCD(I) = 0
   30 IND = IND+1
      DO 35 I=1,M
         IF (B(I,1).EQ.XMISS) B(I,1) = X(I)
   35 CONTINUE
      L = 1
C                                  ACCUMULATE OBSERVATIONS AND A VALUES
C                                  FOR CALCULATIONS ON LAST CALL
      DO 55 I=1,M
         IF (X(I).EQ.XMISS) GO TO 50
         XMEAN(I) = XMEAN(I)+DBLE(X(I))
         X(I) = X(I)-B(I,1)
         Z1 = X(I)
         I2 = M+I
         DO 45 K=1,I
            Z = X(K)
            IF (X(K).EQ.XMISS) GO TO 40
            INCD(L) = INCD(L)+1
            A(K,I) = A(K,I)+(Z*Z1)
            A(K,I2) = A(K,I2)+Z1
            IF (I.EQ.K) GO TO 40
            I3 = M+K
            A(I,I3) = A(I,I3)+X(K)
   40       L = L+1
   45    CONTINUE
         GO TO 55
   50    L = L+I
   55 CONTINUE
      IF (IND.LT.N) GO TO 9005
C                                  N-TH ENTRY.
C                                  COMPUTE MEAN VALUES AND ADJUST A
      L = 0
      DO 65 I=1,M
         KK = L
         L = L+I
         IF (INCD(L).EQ.0) GO TO 65
         XMEAN(I) = XMEAN(I)/INCD(L)
         IF (A(I,I).EQ.0.D0) GO TO 65
         TEMP = B(I,1)-XMEAN(I)
         I2 = M+I
         DO 60 K=1,I
            IF (A(K,I).EQ.0.D0) GO TO 60
            I3 = M+K
            XT = B(K,1)-XMEAN(K)
            A(K,I) = A(K,I)+XT*A(K,I2)+TEMP*A(I,I3)+XT*TEMP*INCD(KK+K)
   60    CONTINUE
   65 CONTINUE
C                                  COMPUTE REGRESSION COEFFICIENTS
C                                  STANDARD ERRORS AND INTERCEPTS
      L1 = 1
      L2 = 1
      K = 1
      DO 105 J=2,M
         L1 = L1+J
         L = 0
         ZL1M1 = INCD(L1)-1.
         IF (ZL1M1.GT.0.) ZJJ = A(J,J)/ZL1M1
         DO 100 I=1,K
            L = L+I
            L2 = L2+1
            ZL2M1 = INCD(L2)-1.
            ZL2M2 = INCD(L2)-2.
            IF (ZL2M1.GT.0.) ZIJ = A(I,J)/ZL2M1
            ZLM1 = INCD(L)-1.
            IF (ZLM1.GT.0.) ZII = A(I,I)/ZLM1
            S(I,J) = 0.
            S(J,I) = 0.
            IF (ZL2M1.GT.1.) GO TO 70
            S(I,J) = SINFM
            S(J,I) = SINFM
            IER = 33
   70       IF (ZL2M1.GE.1.) GO TO 75
            B(I,J) = SINFM
            A(I,J) = SINFM
            B(J,I) = SINFM
            A(J,I) = SINFM
            IER = 34
            GO TO 85
   75       IF (ZII.LE.0) GO TO 90
            B(I,J) = ZIJ/ZII
            A(I,J) = XMEAN(J)-B(I,J)*XMEAN(I)
   80       IF (ZJJ.LE.0) GO TO 95
            B(J,I) = ZIJ/ZJJ
            A(J,I) = XMEAN(I)-B(J,I)*XMEAN(J)
   85       IF (S(I,J).NE.0.) GO TO 100
            DNUM = ZII*ZJJ-ZIJ*ZIJ
            IF (DNUM.LE.0.) GO TO 100
            S(J,I) = SQRT(DNUM/(ZJJ*ZJJ*ZL2M2))
            S(I,J) = SQRT(DNUM/(ZII*ZII*ZL2M2))
            GO TO 100
   90       B(I,J) = SINFM
            A(I,J) = SINFM
            S(I,J) = SINFM
            IER = 35
            GO TO 80
   95       B(J,I) = SINFM
            A(J,I) = SINFM
            S(J,I) = SINFM
            IER = 35
  100    CONTINUE
         K = J
         B(J,J) = 0.
         L2 = L2+1
  105 CONTINUE
      B(1,1) = 0.
C                                  COMPUTE STANDARD DEVIATIONS.
      L = 0
      DO 115 I=1,M
         L = L+I
         IF (INCD(L).LT.2) GO TO 110
         S(I,I) = DSQRT(A(I,I)/(INCD(L)-1.0D0))
         A(I,I) = 0.0D0
         GO TO 115
  110    S(I,I) = SINFM
  115 CONTINUE
      IF (IER.EQ.0) GO TO 9005
 9000 CONTINUE
      CALL UERTST (IER,'BEMIRO')
 9005 RETURN
      END
