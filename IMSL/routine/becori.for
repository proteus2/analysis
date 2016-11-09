C   IMSL ROUTINE NAME   - BECORI
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - ESTIMATES OF MEANS, STANDARD DEVIATIONS, AND
C                           CORRELATION COEFFICIENTS (IN-CORE VERSION)
C
C   USAGE               - CALL BECORI (X,N,M,IX,XM,S,R,IER)
C
C   ARGUMENTS    X      - INPUT/OUTPUT MATRIX OF DIMENSION N BY M.
C                         ON INPUT, X IS A MATRIX OF N OBSERVATIONS
C                           ON M VARIABLES.
C                         ON OUTPUT, THE ELEMENTS OF MATRIX X, HAVE
C                           BEEN ADJUSTED BY THE MEANS OF THE VARIABLES.
C                N      - NUMBER OF OBSERVATIONS. (INPUT)
C                M      - NUMBER OF VARIABLES. (INPUT)
C                IX     - ROW DIMENSION OF MATRIX X EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM. (INPUT)
C                XM     - OUTPUT VECTOR OF LENGTH M CONTAINING THE
C                           MEANS OF THE M VARIABLES.
C                S      - OUTPUT VECTOR OF LENGTH M CONTAINING THE
C                           STANDARD DEVIATIONS OF THE M VARIABLES.
C                R      - OUTPUT MATRIX OF DIMENSION M BY M CONTAINING
C                           THE SIMPLE CORRELATION COEFFICIENTS. R IS
C                           STORED IN SYMMETRIC STORAGE MODE AND
C                           THEREFORE IS A VECTOR OF LENGTH M*(M+1)/2.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         WARNING ERROR
C                           IER=33 INDICATES THAT AT LEAST ONE COLUMN
C                             OF X IS CONSTANT. THE CORRELATION
C                             COEFFICIENTS INVOLVING THE CONSTANT
C                             COLUMN HAVE BEEN SET TO MACHINE INFINITY.
C
C   REQD. IMSL ROUTINES - SINGLE/UERTST,UGETIO
C                       - DOUBLE/UERTST,UGETIO,VXADD,VXMUL,VXSTO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE BECORI (X,N,M,IX,XM,S,R,IER)
C
      DOUBLE PRECISION   TEMP,DM,D,V,ONED,ONEDN
      REAL               X(IX,M),XM(M),S(M),R(1),Y,EPS,
     *                   ONE,BIGNO,DD,ZERO
      DATA               EPS/Z3C100000/,BIGNO/Z7FFFFFFF/
      DATA               ONE/1.0/,ZERO/0.0/
      DATA               ONED/1.D0/
C                                  FIRST EXECUTABLE STATEMENT
C                                  CALCULATE MEANS
      IER=0
      R(1) = ONE
      L=2
      ONEDN = ONED/N
      DO 60 I = 1,M
         TEMP=0.D0
         DO 5 J = 1,N
            TEMP = TEMP+DBLE(X(J,I))
    5    CONTINUE
         XM(I) = TEMP*ONEDN
C                                  ADJUST THE MATRIX BY THE MEANS
         Y=XM(I)
         DD = ZERO
         DO 10 J = 1,N
            X(J,I)=X(J,I)-Y
            D=ABS(X(J,I))
            IF (D .GT. DD) DD = D
   10    CONTINUE
         IF(DD .GT.  ABS(EPS*Y)) GO TO 20
C                                  WARNING ERROR - COLUMN IS CONSTANT
         IER = 33
         S(I) = ZERO
         IF (I .EQ. 1) GO TO 55
         DO 15 J = 1,K
            R(L) = BIGNO
            L = L+1
   15    CONTINUE
         GO TO 50
C                                  CALCULATE STANDARD DEVIATIONS
   20    DM = 0.D0
         DO 25 J = 1,N
            V=DBLE(X(J,I))
            DM = DM+V*V
   25    CONTINUE
         DM=DSQRT(DM)
         S(I)=DM
C                                  CALCULATE CORRELATION COEFFICIENTS
         IF (I .EQ. 1) GO TO 55
         DO 45 J = 1,K
            D=0.D0
            DO 30 JJ = 1,N
               D = D+DBLE(X(JJ,I))*DBLE(X(JJ,J))
   30       CONTINUE
            IF (S(J) .EQ. ZERO) GO TO 35
            R(L) = D/(DM*S(J))
            GO TO 40
   35       R(L) = BIGNO
   40       L = L+1
   45    CONTINUE
   50    R(L) = ONE
         L=L+1
   55    K = I
   60 CONTINUE
      D = ONE/SQRT(N-ONE)
      DO 65 I = 1,M
         S(I) = S(I)*D
   65 CONTINUE
      IF (IER .EQ. 0) GO TO 9005
 9000 CONTINUE
      CALL UERTST (IER,'BECORI')
 9005 RETURN
      END
