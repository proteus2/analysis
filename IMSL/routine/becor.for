C   IMSL ROUTINE NAME   - BECOR
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - ESTIMATES OF MEANS, STANDARD DEVIATIONS, AND
C                           CORRELATION COEFFICIENTS (OUT-OF-CORE
C                           VERSION)
C
C   USAGE               - CALL BECOR (X,N,M,I,IND,TEMP,XMD,SD,RD,IER)
C
C   ARGUMENTS    X      - INPUT/OUTPUT VECTOR OF LENGTH M.
C                         ON INPUT, X IS THE I-TH ROW OF AN N BY M
C                           MATRIX OF OBSERVATIONS.
C                         ON OUTPUT, X IS THE I-TH ROW OF THE N BY M
C                           MATRIX OF OBSERVATIONS ADJUSTED BY THE
C                           TEMPORARY MEANS.
C                N      - NUMBER OF OBSERVATIONS PER VARIABLE. N MUST
C                           BE GREATER THAN 1. (INPUT)
C                M      - NUMBER OF VARIABLES. (INPUT)
C                I      - NUMBER OF THE ROW ENTERED TO BECOR.  (INPUT)
C                IND    - TEMPORARY MEAN INDICATOR. (INPUT)
C                           IF IND=1, TEMPORARY MEANS ARE USER SUPPLIED.
C                           OTHERWISE, THEY ARE NOT.
C                TEMP   - DOUBLE PRECISION WORKING VECTOR OF LENGTH 2M
C                           USED AS WORKING STORAGE. IF IND = 1, TEMP
C                           CONTAINS THE USER SUPPLIED TEMPORARY MEANS
C                           ON THE FIRST ENTRY IN THE FIRST M LOCATIONS.
C                XMD    - OUTPUT VECTOR OF LENGTH M CONTAINING THE MEANS
C                           OF THE M VARIABLES. XMD MUST BE TYPED
C                           DOUBLE PRECISION.
C                SD     - OUTPUT VECTOR OF LENGTH M CONTAINING THE
C                           STANDARD DEVIATIONS OF THE M VARIABLES. SD
C                           MUST BE DOUBLE PRECISION.
C                RD     - OUTPUT MATRIX OF DIMENSION M BY M CONTAINING
C                           THE SIMPLE CORRELATION COEFFICIENTS. RD IS
C                           STORED IN SYMMETRIC STORAGE MODE AND
C                           THEREFORE REQUIRES AT LEAST M*(M+1)/2
C                           STORAGE LOCATIONS. RD MUST BE DOUBLE
C                           PRECISION.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER=129 INDICATES THAT I IS GREATER THAN N
C                             OR LESS THAN 1.
C                         WARNING ERROR
C                           IER=34 INDICATES THAT AT LEAST ONE COLUMN
C                             OF X IS CONSTANT. THE CORRELATION
C                             COEFFICIENTS INVOLVING THE CONSTANT
C                             COLUMN HAVE BEEN SET TO MACHINE INFINITY.
C
C   REQD. IMSL ROUTINES - UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS      WHEN SINGLE PRECISION RESULTS ARE REQUIRED, THE
C                DOUBLE PRECISION OUTPUT OF BECOR MUST BE CONVERTED
C                IN THE CALLING PROGRAM.  SEE PROGRAMMING NOTES.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE BECOR  (X,N,M,I,IND,TEMP,XMD,SD,RD,IER)
C
      REAL               X(M),EPS
      DOUBLE PRECISION   XMD(M),SD(M),RD(1),TEMP(1),Z,RN,ONEDN,D,ONE,
     *                   ZZ,BIGNO,ZERO
      DATA               BIGNO/Z7FFFFFFFFFFFFFFF/
      DATA               ONE/1.0D0/,ZERO/0.0D0/
      DATA               EPS/Z3C100000/
C                                  FIRST EXECUTABLE STATEMENT
      IER=0
      IF (I .LE. N .AND. I .GE. 1) GO TO 5
C                                  TERMINAL ERROR - I IS GREATER THAN
C                                  N OR LESS THAN 1
      IER=129
      GO TO 9000
C                                  FIRST ENTRY
    5 IF (I .GT. 1) GO TO 30
      IF (IND .EQ. 1) GO TO 15
C                                  USE FIRST ROW AS TEMPORARY MEANS
      DO 10 J = 1,M
         TEMP(J) = X(J)
   10 CONTINUE
C                                  INITIALIZE MEANS, STANDARD DEVIATIONS
C                                  AND CORRELATION COEFFICIENTS
   15 DO 20 J = 1,M
         TEMP(M+J) = ZERO
         XMD(J)=ZERO
         SD(J)=ZERO
   20 CONTINUE
      JJ = (M*(M+1))/2
      DO 25 J = 1,JJ
         RD(J) = ZERO
   25 CONTINUE
   30 L = 1
      DO 45 J = 1,M
         XMD(J)=X(J)+XMD(J)
         Z=X(J)-TEMP(J)
         ZZ = DABS(Z)
         IF (ZZ .GT. TEMP(M+J)) TEMP(J+M) = ZZ
         X(J)=Z
         SD(J)=SD(J)+Z*Z
         IF (J .EQ. 1) GO TO 40
         DO 35 K = 1,LL
            RD(L)=RD(L)+Z*X(K)
            L = L+1
   35    CONTINUE
   40    L = L+1
         LL = J
   45 CONTINUE
      IF (I .LT. N) GO TO 9005
C                                  N-TH ENTRY
      L=1
      RN = N
      ONEDN = ONE/N
      DO 80 J = 1,M
C                                  FIND MEAN
         XMD(J) = XMD(J)*ONEDN
         IF (TEMP(M+J) .GT. DABS(EPS*XMD(J))) GO TO 55
         IER = 34
         SD(J) = ZERO
         IF (J .EQ. 1) GO TO 75
         DO 50 K = 1,II
            RD(L) = BIGNO
            L = L+1
   50    CONTINUE
         GO TO 75
   55    Z = XMD(J)-TEMP(J)
         D = RN*Z
         SD(J) = DSQRT(SD(J)-D*Z)
C                                  FIND CORRELATION COEFFICIENTS
         IF (J .EQ. 1) GO TO 75
         DO 70 K = 1,II
            IF (SD(K) .EQ. ZERO) GO TO 60
            RD(L) = (RD(L)-D*(XMD(K)-TEMP(K)))/(SD(J)*SD(K))
            GO TO 65
   60       RD(L) = BIGNO
   65       L = L+1
   70    CONTINUE
   75    RD(L) = ONE
         II = J
         L = L+1
   80 CONTINUE
C                                  FIND STANDARD DEVIATIONS
      Z = ONE/DSQRT(RN-ONE)
      DO 85 J = 1,M
         SD(J) = SD(J)*Z
   85 CONTINUE
      IF (IER .EQ. 0) GO TO 9005
 9000 CONTINUE
      CALL UERTST (IER,'BECOR ')
 9005 RETURN
      END
