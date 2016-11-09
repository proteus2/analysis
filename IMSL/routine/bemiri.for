C   IMSL ROUTINE NAME   - BEMIRI
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1980
C
C   PURPOSE             - ESTIMATES OF MEANS, SIMPLE REGRESSION
C                           COEFFICIENTS, THEIR INTERCEPTS, STANDARD
C                           ERRORS OF THE REGRESSION COEFFICIENTS, AND
C                           STANDARD DEVIATIONS FOR ARRAYS WHICH
C                           CONTAIN MISSING VALUES. (IN-CORE VERSION)
C
C   USAGE               - CALL BEMIRI (X,N,M,IX,XMEAN,B,A,S,IBAS,INCD,
C                           IER)
C
C   ARGUMENTS    X      - INPUT MATRIX OF DIMENSION N BY M CONTAINING
C                           N OBSERVATIONS ON M VARIABLES.
C                         ON OUTPUT, X EQUALS X ADJUSTED (X-XMEAN).
C                N      - INPUT.  NUMBER OF OBSERVATIONS PER VARIABLE.
C                           N MUST BE GREATER THAN 1.
C                M      - INPUT.  NUMBER OF VARIABLES.  M MUST BE
C                           GREATER THAN 1.
C                IX     - ROW DIMENSION OF MATRIX X EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM. (INPUT)
C                XMEAN  - OUTPUT VECTOR OF LENGTH M CONTAINING MEANS
C                           OF THE M COLUMNS OF X.
C                B      - OUTPUT MATRIX CONTAINING THE REGRESSION
C                           COEFFICIENTS.  B MUST BE OF LENGTH AT LEAST
C                           M BY M. B(I,J) CONTAINS THE COEFFICIENT
C                           FOR THE REGRESSION OF VARIABLE J ON
C                           VARIABLE I. THE DIAGONAL OF B IS NOT USED.
C                A      - OUTPUT MATRIX CONTAINING THE INTERCEPTS.
C                           A MUST BE OF LENGTH AT LEAST M BY M.
C                           A(I,J) CONTAINS THE INTERCEPT FOR THE
C                           REGRESSION OF VARIABLE J ON VARIABLE I. THE
C                S      - OUTPUT MATRIX CONTAINING THE STANDARD ERRORS
C                           OF THE REGRESSION COEFFICIENTS AND THE
C                           STANDARD DEVIATIONS.  THE STANDARD
C                           DEVIATIONS LIE ON THE DIAGONAL OF THE
C                           MATRIX.  S MUST BE OF LENGTH AT LEAST
C                           M BY M.  SEE DOCUMENT.
C                IBAS   - ROW DIMENSION OF MATRICES B, A, AND S EXACTLY
C                           AS SPECIFIED IN  DIMENSION STATEMENT IN THE
C                           THE CALLING PROGRAM. (INPUT)
C                INCD   - OUTPUT VECTOR CONTAINING THE INCIDENCE VALUES.
C                           INDICATES THE NUMBER OF PAIRS OF
C                           OBSERVATIONS USED IN THE CALCULATIONS OF B,
C                           A, AND S.  INCD IS STORED IN SYMMETRIC
C                           MATRIX MODE AND MUST BE OF LENGTH AT LEAST
C                           M(M+1)/2.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         WARNING ERROR
C                           IER=33 INDICATES THAT FEWER THAN 3 PAIRS OF
C                             OBSERVATIONS WERE PRESENT. THE STANDARD
C                             ERRORS ARE SET TO NEGATIVE MACHINE
C                             INFINITY.
C                           IER=34 INDICATES THAT FEWER THAN 2 PAIRS OF
C                             OBSERVATIONS WERE PRESENT. THE STANDARD
C                             DEVIATIONS, STANDARD ERRORS, REGRESSION
C                             COEFFICIENTS, AND INTERCEPTS ARE SET TO
C                             NEGATIVE MACHINE INFINITY.
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
      SUBROUTINE BEMIRI (X,N,M,IX,XMEAN,B,A,S,IBAS,INCD,IER)
C
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N,M,IX,IBAS,INCD(1),IER
      REAL               X(IX,M),XMEAN(M),B(IBAS,M),A(IBAS,M),S(IBAS,M)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      DOUBLE PRECISION   TEMP
      DATA               XMISS/-999999./
      DATA               SINFM/ZFFFFFFFF/
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
C                                  COMPUTE THE MEAN VALUES
      DO 10 I=1,M
         TEMP = 0.0D0
         XN = 0.0
         XMEAN(I) = 0.0
         DO 5 J=1,N
            IF (X(J,I).EQ.XMISS) GO TO 5
            XN = XN+1.
            TEMP = TEMP+DBLE(X(J,I))
    5    CONTINUE
         IF (XN.EQ.0.) GO TO 10
         XMEAN(I) = TEMP/XN
   10 CONTINUE
      DO 20 J=1,M
         DO 15 I=1,N
            IF (X(I,J).EQ.XMISS) GO TO 15
            X(I,J) = X(I,J)-XMEAN(J)
   15    CONTINUE
   20 CONTINUE
      L = 1
      DO 35 J=1,M
         DO 30 I=1,J
            INCD(L) = 0
            TEMP = 0.0D0
            DO 25 K=1,N
               IF (X(K,J).EQ.XMISS.OR.X(K,I).EQ.XMISS) GO TO 25
               INCD(L) = INCD(L)+1
               TEMP = TEMP+(DBLE(X(K,J))*DBLE(X(K,I)))
   25       CONTINUE
            A(I,J) = TEMP
   30    L = L+1
   35 CONTINUE
C                                  COMPUTE REGRESSION COEFFICIENTS
C                                  STANDARD ERRORS AND INTERCEPTS
      L1 = 1
      L2 = 1
      K = 1
      DO 75 J=2,M
         L1 = L1+J
         L = 0
         ZL1M1 = INCD(L1)-1.
         IF (ZL1M1.GT.0.) ZJJ = A(J,J)/ZL1M1
         DO 70 I=1,K
            L = L+I
            L2 = L2+1
            ZL2M1 = INCD(L2)-1.
            ZL2M2 = INCD(L2)-2.
            IF (ZL2M1.GT.0.) ZIJ = A(I,J)/ZL2M1
            ZLM1 = INCD(L)-1.
            IF (ZLM1.GT.0.) ZII = A(I,I)/ZLM1
            S(I,J) = 0.
            S(J,I) = 0.
            IF (ZL2M1.GT.1) GO TO 40
            S(I,J) = SINFM
            S(J,I) = SINFM
            IER = 33
   40       IF (ZL2M1.GE.1) GO TO 45
            B(I,J) = SINFM
            A(I,J) = SINFM
            B(J,I) = SINFM
            A(J,I) = SINFM
            IER = 34
            GO TO 55
   45       IF (ZII.LE.0) GO TO 60
            B(I,J) = ZIJ/ZII
            A(I,J) = XMEAN(J)-B(I,J)*XMEAN(I)
   50       IF (ZJJ.LE.0) GO TO 65
            B(J,I) = ZIJ/ZJJ
            A(J,I) = XMEAN(I)-B(J,I)*XMEAN(J)
   55       IF (S(I,J).NE.0.) GO TO 70
            DNUM = ZII*ZJJ-ZIJ*ZIJ
            IF (DNUM.LE.0.) GO TO 70
            S(J,I) = SQRT(DNUM/(ZJJ*ZJJ*ZL2M2))
            S(I,J) = SQRT(DNUM/(ZII*ZII*ZL2M2))
            GO TO 70
   60       B(I,J) = SINFM
            A(I,J) = SINFM
            S(I,J) = SINFM
            IER = 35
            GO TO 50
   65       B(J,I) = SINFM
            A(J,I) = SINFM
            S(J,I) = SINFM
            IER = 35
   70    CONTINUE
         K = J
         B(J,J) = 0.
         L2 = L2+1
   75 CONTINUE
      B(1,1) = 0.
C                                  COMPUTE STANDARD DEVIATIONS.
      L = 0
      DO 85 I=1,M
         L = L+I
         IF (INCD(L).LT.2) GO TO 80
         S(I,I) = SQRT(A(I,I)/(INCD(L)-1.))
         A(I,I) = 0.0
         GO TO 85
   80    S(I,I) = SINFM
   85 CONTINUE
      IF (IER.EQ.0) GO TO 9005
 9000 CONTINUE
      CALL UERTST (IER,'BEMIRI')
 9005 RETURN
      END
