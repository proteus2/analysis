C   IMSL ROUTINE NAME   - BEMMO
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1980
C
C   PURPOSE             - ESTIMATES OF MEANS, STANDARD DEVIATIONS,
C                           CORRELATION COEFFICIENTS, AND COEFFICIENTS
C                           OF SKEWNESS AND KURTOSIS FROM A DATA MATRIX
C                           CONTAINING MISSING OBSERVATIONS (OUT OF
C                           CORE VERSION)
C
C   USAGE               - CALL BEMMO (X,N,M,IND,XMEAN,U,A,INCD,WK,IER)
C
C   ARGUMENTS    X      - INPUT VECTOR OF LENGTH M CONTAINING ONE
C                           OBSERVATION FOR EACH OF M VARIABLES.
C                         ON OUTPUT, X IS EQUAL TO X MINUS THE
C                           TEMPORARY MEAN.
C                N      - INPUT.  NUMBER OF OBSERVATIONS PER VARIABLE.
C                           N MUST BE GREATER THAN 1.
C                M      - INPUT.  NUMBER OF VARIABLES.
C                IND    - INPUT INDICATOR FOR MEAN VALUES-INITIAL ENTRY.
C                           ON INITIAL ENTRY IND MUST BE LESS THAN OR
C                             EQUAL TO ZERO.
C                           IF IND EQUALS ZERO THE FIRST OBSERVATION ON
C                             THE M VARIABLES WILL BE USED AS A
C                             TEMPORARY MEAN VECTOR.
C                           IF IND IS LESS THAN ZERO THE FIRST ROW OF
C                             U MUST CONTAIN THE TEMPORARY MEANS.
C                         OUTPUT COUNT FOR THE NUMBER OF
C                           OBSERVATIONS THAT HAVE BEEN ENTERED.  IND
C                           IS INCREMENTED EACH TIME BEMMO IS CALLED.
C                XMEAN  - OUTPUT VECTOR OF LENGTH M CONTAINING MEAN
C                           VALUES. XMEAN MUST BE DOUBLE PRECISION.
C                U      - IF IND.LT.0 ON INITIAL ENTRY, THEN ON ALL
C                           ENTRIES TO BEMMO U(1,I) FOR I=1,M MUST
C                           CONTAIN THE TEMPORARY MEANS.
C                         OUTPUT MATRIX OF DIMENSION 3 BY M CONTAINING
C                           FOR I=1,...,M
C                             U(1,I)=STANDARD DEVIATIONS.
C                             U(2,I)=COEFFICIENTS OF SKEWNESS.
C                             U(3,I)=COEFFICIENTS OF EXCESS (KURTOSIS).
C                A      - OUTPUT VECTOR (DOUBLE PRECISION) CONTAINING
C                           THE CORRELATION COEFFICIENTS.  A IS STORED
C                           IN SYMMETRIC STORAGE MODE AND MUST BE
C                           OF LENGTH AT LEAST 3(M(M+1)/2).
C                           THE FIRST (M(M+1))/2 LOCATIONS CONTAIN THE
C                           CORRELATION COEFFICIENTS.  THE REMAINING
C                           LOCATIONS ARE WORK STORAGE.
C                INCD   - OUTPUT VECTOR CONTAINING THE INCIDENCE VALUES.
C                           INDICATES THE NUMBER OF PAIRS OF
C                           OBSERVATIONS USED IN THE CALCULATIONS OF
C                           U AND A.  INCD IS STORED IN SYMMETRIC
C                           STORAGE MODE AND MUST BE OF LENGTH AT LEAST
C                           M(M+1)/2.
C                WK     - WORK AREA MATRIX WHICH MUST BE DIMENSIONED
C                           AT LEAST 2 X M IN LENGTH.  IT IS USED TO
C                           ACCUMULATE VALUES FOR 3RD AND 4TH MOMENTS.
C                           WK MUST BE DOUBLE PRECISION.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         WARNING ERROR
C                           IER=33 INDICATES THAT THE OBSERVATIONS ON
C                             SOME VARIABLE WERE CONSTANT. THE PERTINENT
C                             CORRELATION COEFFICIENTS, AND COEFFICIENTS
C                             OF SKEWNESS AND EXCESS ARE SET TO NEGATIVE
C                             MACHINE INFINITY.
C                           IER=34 INDICATES THAT FEWER THAN TWO
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
      SUBROUTINE BEMMO  (X,N,M,IND,XMEAN,U,A,INCD,WK,IER)
C
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N,M,IND,INCD(1),IER
      DOUBLE PRECISION   A(1),WK(2,M),XMEAN(M)
      REAL               X(M),U(3,M)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      DOUBLE PRECISION   TEMP,XT,XTT
      DATA               XMISS/-999999./
      DATA               SINFM/ZFFFFFFFF/
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
C                                  TEST FOR FIRST CALL OF A SEQUENCE
      IF (IND.GT.0) GO TO 25
C                                  TEST FOR LOCATION OF TEMPORARY MEANS
      IF (IND.NE.0) GO TO 10
C                                  GET TEMPORARY MEANS FROM X ENTRIES
      DO 5 I=1,M
    5 U(1,I) = X(I)
C                                  INITIALIZE FOR ACCUMULATION PURPOSES
   10 K = (M*(M+1))/2
      KK = K+K
      IND = 0
      DO 15 I=1,K
         A(I) = 0.0D0
         A(K+I) = 0.0D0
         A(KK+I) = 0.0D0
   15 INCD(I) = 0
      DO 20 I=1,M
         U(2,I) = 0.0
         U(3,I) = 0.0
         WK(1,I) = 0.0D0
         WK(2,I) = 0.0D0
   20 XMEAN(I) = 0.0D0
C                                  ACCUMULATE OBSERVATIONS AND A VALUES
C                                  FOR CALCULATIONS ON LAST CALL
   25 IND = IND+1
      DO 30 I=1,M
         IF (U(1,I).EQ.XMISS) U(1,I) = X(I)
   30 CONTINUE
      K = (M*(M+1))/2
      KK = K+K
      L = 1
      DO 45 I=1,M
         IF (X(I).EQ.XMISS) GO TO 35
         XMEAN(I) = XMEAN(I)+DBLE(X(I))
         X(I) = X(I)-U(1,I)
   35    DO 40 J=1,I
            K = K+1
            KK = KK+1
            IF (X(J).EQ.XMISS.OR.X(I).EQ.XMISS) GO TO 40
            INCD(L) = INCD(L)+1
            A(L) = A(L)+(DBLE(X(I))*DBLE(X(J)))
            A(K) = A(K)+DBLE(X(J))
            A(KK) = A(KK)+DBLE(X(I))
   40    L = L+1
   45 CONTINUE
C                                  TEST FOR CONSTANT OBSERVATION SET
C                                  AND ACCUMULATE FOR 3RD AND 4TH MOMENT
      IF (IND.EQ.1) GO TO 55
      DO 50 I=1,M
         IF (X(I).EQ.XMISS) GO TO 50
         IF (X(I).NE.U(2,I)) U(3,I) = 1.0
   50 CONTINUE
   55 DO 60 I=1,M
         IF (X(I).EQ.XMISS) GO TO 60
         TEMP = DBLE(X(I))**3
         WK(1,I) = WK(1,I)+TEMP
         WK(2,I) = WK(2,I)+(TEMP*DBLE(X(I)))
         U(2,I) = X(I)
   60 CONTINUE
      IF (IND.LT.N) GO TO 9005
C                                  COMPUTE STANDARD DEVIATIONS, 3RD
C                                  AND 4TH MOMENTS
      L = 0
      DO 65 I=1,M
         L = L+I
         IF (U(3,I).EQ.0.0) A(L) = 0.0D0
   65 CONTINUE
      K0 = KK-K
      K1 = 0
      K2 = K0
      K3 = K
      L = 0
      DO 85 I=1,M
         L = L+I
         IF (INCD(L).LT.2) GO TO 75
         XT = INCD(L)
         XMEAN(I) = XMEAN(I)/XT
         IF (A(L).EQ.0.D0) GO TO 80
         TEMP = U(1,I)-XMEAN(I)
         XTT = TEMP*TEMP
         WK(2,I) = WK(2,I)+4.D0*TEMP*WK(1,I)+6.D0*XTT*A(L)-3.D0*XT*XTT
     1   *XTT
         WK(1,I) = WK(1,I)+3.D0*TEMP*A(L)-2.D0*XT*(TEMP*XTT)
         XTT = (A(L)/XT)-XTT
         U(2,I) = WK(1,I)/((DSQRT(XTT)**3)*XT)
         U(3,I) = (WK(2,I)/((XTT**2)*XT))-3.D0
         DO 70 J=1,I
            K1 = K1+1
            K2 = K2+1
            K3 = K3+1
            IF (A(K1).EQ.0.D0) GO TO 70
            XT = U(1,J)-XMEAN(J)
            A(K1) = A(K1)+XT*A(K3)+TEMP*A(K2)+XT*TEMP*INCD(K1)
   70    CONTINUE
         GO TO 85
   75    IER = 34
   80    U(2,I) = SINFM
         U(3,I) = SINFM
         K1 = L
         K2 = K0+K1
         K3 = K+K1
   85 CONTINUE
C                                  COMPUTE CORRELATION COEFFICIENTS
      L1 = 1
      L2 = 2
      U(1,1) = SINFM
      IF (INCD(1).GT.1) U(1,1) = DSQRT(A(1)/(INCD(1)-1))
      IF (M.LT.2) GO TO 110
      DO 105 J=2,M
         L1 = L1+J
         U(1,J) = SINFM
         IF (INCD(L1).GT.1) U(1,J) = DSQRT(A(L1)/(INCD(L1)-1))
         L = 0
         K = J-1
         DO 100 I=1,K
            L = L+I
            IF (A(L).EQ.0.D0.OR.A(L1).EQ.0.D0) GO TO 90
            IF (INCD(L2).LT.2) GO TO 90
            A(L2) = (A(L2)*DSQRT((INCD(L)-1.D0)*(INCD(L1)-1.D0)))
     1      /((INCD(L2)-1.D0)*DSQRT(A(L1)*A(L)))
            IF (A(L2).GT.1.D0) A(L2) = 1.D0
            IF (A(L2).LT.-1.D0) A(L2) = -1.D0
            GO TO 95
   90       A(L2) = SINFM
            IF (IER.EQ.0) IER = 33
   95       L2 = L2+1
  100    CONTINUE
         L2 = L2+1
  105 CONTINUE
  110 CONTINUE
C                                  SET DIAGONAL OF A TO ONE
      L = 0
      DO 115 I=1,M
         L = L+I
         A(L) = 1.0D0
  115 CONTINUE
      IF (IER.EQ.0) GO TO 9005
 9000 CONTINUE
      CALL UERTST (IER,'BEMMO ')
 9005 RETURN
      END
