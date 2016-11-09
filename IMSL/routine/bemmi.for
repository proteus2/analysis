C   IMSL ROUTINE NAME   - BEMMI
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1980
C
C   PURPOSE             - ESTIMATES OF MEANS, STANDARD DEVIATIONS,
C                           CORRELATION COEFFICIENTS, AND COEFFICIENTS
C                           OF SKEWNESS AND KURTOSIS FROM A DATA MATRIX
C                           CONTAINING MISSING OBSERVATIONS (IN-CORE
C                           VERSION)
C
C   USAGE               - CALL BEMMI (X,N,M,IX,XMEAN,U,A,INCD,IER)
C
C   ARGUMENTS    X      - INPUT MATRIX OF DIMENSION N BY M CONTAINING
C                           N OBSERVATIONS ON M VARIABLES.
C                         ON OUTPUT, X IS SET TO X ADJUSTED (X-XMEAN),
C                           EXCEPTING THE MISSING X READINGS.
C                N      - INPUT.  NUMBER OF OBSERVATIONS PER VARIABLE.
C                           N MUST BE GREATER THAN 1.
C                M      - INPUT.  NUMBER OF VARIABLES.
C                IX     - INPUT. ROW DIMENSION OF MATRIX X EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM.
C                XMEAN  - OUTPUT VECTOR OF LENGTH M CONTAINING MEAN
C                           VALUES.
C                U      - OUTPUT MATRIX OF DIMENSION 3 BY M CONTAINING
C                           FOR I=1,...,M
C                             U(1,I)=STANDARD DEVIATIONS.
C                             U(2,I)=COEFFICIENTS OF SKEWNESS.
C                             U(3,I)=COEFFICIENTS OF EXCESS (KURTOSIS).
C                A      - OUTPUT VECTOR CONTAINING THE CORRELATION
C                           COEFFICIENTS.  A IS STORED IN SYMMETRIC
C                           STORAGE MODE AND MUST BE OF LENGTH AT LEAST
C                           M(M+1)/2.
C                INCD   - OUTPUT VECTOR CONTAINING INCIDENCE VALUES.
C                           INDICATES THE NUMBER OF PAIRS OF
C                           OBSERVATIONS USED IN THE CALCULATIONS OF U
C                           AND A.  INCD IS STORED IN SYMMETRIC STORAGE
C                           MODE AND MUST BE OF LENGTH AT LEAST
C                           M(M+1)/2.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         WARNING ERROR
C                           IER=34 INDICATES THAT THE OBSERVATIONS ON
C                             SOME VARIABLE WERE CONSTANT. THE PERTINENT
C                             CORRELATION COEFFICIENTS AND COEFFICIENTS
C                             OF SKEWNESS AND EXCESS ARE SET TO NEGATIVE
C                             MACHINE INFINITY.
C                           IER=35 INDICATES THAT FEWER THAN TWO
C                             OBSERVATIONS WERE PRESENT.  THE STANDARD
C                             DEVIATIONS, CORRELATION COEFFICIENTS, AND
C                             COEFFICIENTS OF SKEWNESS AND EXCESS ARE
C                             SET TO NEGATIVE MACHINE INFINITY.
C
C   REQD. IMSL ROUTINES - UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE BEMMI  (X,N,M,IX,XMEAN,U,A,INCD,IER)
C
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            INCD(1),N,M,IX,IER
      REAL               X(IX,M),XMEAN(M),U(3,M),A(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      DOUBLE PRECISION   TEMP,TEMP1,TEMP2
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
            IF(X(J,I) .EQ. XMISS) GO TO 5
            XN = XN+1.
            TEMP = TEMP+DBLE(X(J,I))
    5    CONTINUE
         IF(XN .EQ. 0.) GO TO 10
         XMEAN(I) = TEMP/XN
   10 CONTINUE
      DO 20 J=1,M
         DO 15 I=1,N
            IF(X(I,J) .EQ. XMISS) GO TO 15
            X(I,J) = X(I,J)-XMEAN(J)
   15    CONTINUE
   20 CONTINUE
      L = 1
      DO 30 J=1,M
         DO 30 I=1,J
            INCD(L) = 0
            TEMP = 0.0D0
            DO 25 K=1,N
               IF(X(K,J) .EQ. XMISS .OR. X(K,I) .EQ. XMISS) GO TO 25
               INCD(L) = INCD(L)+1
               TEMP = TEMP+(DBLE(X(K,J))*DBLE(X(K,I)))
   25       CONTINUE
            A(L) = TEMP
            L = L+1
   30 CONTINUE
C                                  COMPUTE STANDARD DEVIATIONS,
C                                  COEFFICIENTS OF SKEWNESS, KURTOSIS.
      L = 0
      DO 50 I=1,M
         L = L+I
         IF(INCD(L) .LT. 2 ) GO TO 40
         TEMP1 = 0.0D0
         TEMP2 = 0.0D0
         DO 35 J=1,N
            IF(X(J,I) .EQ. XMISS) GO TO 35
            TEMP = DBLE(X(J,I))**3
            TEMP1 = TEMP1+(TEMP*DBLE(X(J,I)))
            TEMP2 = TEMP2+TEMP
   35    CONTINUE
         XT = INCD(L)
         XTT = A(L)/XT
         U(1,I) = SQRT(A(L)/(XT-1.))
         IF(A(L) .EQ. 0.0) GO TO 45
         U(2,I) = TEMP2 / ((SQRT(XTT)**3)*XT)
         U(3,I) = (TEMP1/((XTT**2)*XT)) - 3.
         GO TO 50
   40    U(1,I) = SINFM
         IER = 34
   45    U(2,I) = SINFM
         U(3,I) = SINFM
   50 CONTINUE
      IF (M.LT.2) GO TO 75
C                                  COMPUTE THE CORRELATION COEFFICIENTS
      L1 = 1
      L2 = 2
      DO 70 J=2,M
         L1 = L1+J
         L = 0
         K = J-1
         DO 65 I=1,K
            L = L+I
            IF(A(L) .EQ. 0. .OR. A(L1) .EQ. 0.) GO TO 55
            IF(INCD(L2) .LT. 2) GO TO 55
            A(L2) = (A(L2)*SQRT((INCD(L)-1.)*(INCD(L1)-1.))) /
     *        ((INCD(L2)-1)*SQRT(A(L1)*A(L)))
            IF(A(L2) .GT. 1.) A(L2) = 1.
            IF(A(L2) .LT. -1.) A(L2) = -1.
            GO TO 60
   55       A(L2) = SINFM
            IF(IER .EQ. 0) IER = 33
   60       L2 = L2+1
   65    CONTINUE
         L2 = L2+1
   70 CONTINUE
   75 CONTINUE
      L = 0
      DO 80 I=1,M
         L = L+I
         A(L) = 1.0
   80 CONTINUE
      IF(IER .EQ. 0) GO TO 9005
 9000 CONTINUE
      CALL UERTST(IER,'BEMMI ')
 9005 RETURN
      END
