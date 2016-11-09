C   IMSL ROUTINE NAME   - RLDCVA
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - VARIANCE ESTIMATES FOR DECODED ORTHOGONAL
C                           POLYNOMIAL REGRESSION COEFFICIENTS
C
C   USAGE               - CALL RLDCVA (V,ID,A,B,SC,T,IT,IER)
C
C   ARGUMENTS    V      - INPUT VECTOR OF LENGTH ID+1 CONTAINING THE
C                           VARIANCES OF THE CODED INTERCEPT AND
C                           REGRESSION COEFFICIENTS. SEE REMARKS.
C                         ON OUTPUT, V IS A VECTOR OF LENGTH ID+1
C                           CONTAINING THE VARIANCES OF THE DECODED
C                           INTERCEPT AND REGRESSION COEFFICIENTS.
C                ID     - INPUT DEGREE OF THE FITTED MODEL.
C                A      - INPUT VECTOR OF LENGTH ID CONTAINING
C                           CONSTANTS USED IN GENERATING THE ORTHOGONAL
C                           POLYNOMIALS. SEE REMARKS.
C                B      - INPUT VECTOR OF LENGTH ID CONTAINING
C                           ADDITIONAL CONSTANTS USED IN GENERATING THE
C                           ORTHOGONAL POLYNOMIALS. SEE REMARKS.
C                SC     - INPUT VECTOR OF LENGTH 2 CONTAINING THE
C                           SCALING CONSTANTS USED IN THE LINEAR
C                           TRANSFORMATION SC(1)*X+SC(2) WHICH MAPS THE
C                           VALUES OF X, THE INDEPENDENT VARIABLE, TO
C                           THE INTERVAL (-2,2), INCLUSIVELY.
C                T      - DOUBLE PRECISION WORK MATRIX OF DIMENSION
C                           (ID+1) BY (ID+3).
C                IT     - INPUT ROW DIMENSION OF MATRIX T EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT
C                           IN THE CALLING PROGRAM.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER = 129 INDICATES THAT ONE OR MORE
C                             ELEMENTS OF V WERE LESS THAN OR EQUAL TO
C                             ZERO ON INPUT.
C
C   REQD. IMSL ROUTINES - UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS  1.  RLDCVA IS DESIGNED TO BE USED IN CONJUNCTION WITH
C                IMSL ROUTINES RLFOTH AND RLFOTW. THE INPUT A AND B
C                USED BY RLDCVA ARE COMPUTED BY BOTH RLFOTH AND RLFOTW.
C                THE INPUT V MAY BE COMPUTED BY RLDCW.
C            2.  RLDCVA ASSUMES THE INDEPENDENT VARIABLE WAS CODED TO
C                THE INTERVAL (-2,2), INCLUSIVELY, BEFORE THE MODEL
C                WAS FITTED BY THE ORTHOGONAL POLYNOMIAL METHOD.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE RLDCVA (V,ID,A,B,SC,T,IT,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            ID,IT,IER
      REAL               V(1),A(1),B(1),SC(1)
      DOUBLE PRECISION   T(IT,1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            ID1,ID2,ID3,I,J,K,I1,J1,K1,L,IONE,ITWO
      REAL               ZERO
      DOUBLE PRECISION   BETA,SSS,SS,S,SL,SKL,SUM,DZERO,ONE
      DATA               DZERO/0.0D0/,ONE/1.0D0/
      DATA               ZERO/0.0/
      DATA               IONE/1/,ITWO/2/
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
      ID1 = ID+1
      ID2 = ID1+1
      ID3 = ID2+1
      T(IONE,ITWO) = -A(IONE)
      DO 5 I = 1,ID1
C                                  TERMINAL ERROR - SOME INPUT VARIANCE
C                                  WAS LESS THAN OR EQUAL TO ZERO
         IF (V(I) .LE. ZERO) GO TO 90
C                                  CALCULATE COEFFICIENTS OF THE CODED
C                                  REGRESSION COEFFICIENT ESTIMATES,
C                                  USED TO OBTAIN PARTIALLY DECODED
C                                  ESTIMATES
         T(I,I) = ONE
         IF (I .GT. 2) T(IONE,I) = -T(IONE,K)*A(K)-T(IONE,I-2)*B(K)
         K = I
    5 CONTINUE
      IF (ID .LT. 2) GO TO 17
      DO 15 J = 2,ID
         J1 = J-1
         I1 = J+2
         K = J+1
         T(J,K) = T(J1,J)-A(J)
         IF (I1 .GT. ID1) GO TO 15
         DO 10 I = I1,ID1
            T(J,I) = T(J1,K)-A(K)*T(J,K)-B(K)*T(J,K-1)
            K = I
   10    CONTINUE
   15 CONTINUE
C                                  CALCULATE VARIANCES AND COVARIANCES
C                                  OF PARTIALLY DECODED ESTIMATES
   17 CONTINUE
      DO 30 I = 1,ID1
         J1 = I+1
         SUM = V(I)
         IF (J1 .GT. ID1) GO TO 25
         DO 20 J = J1,ID1
            SUM = SUM+V(J)*T(I,J)*T(I,J)
   20    CONTINUE
   25    T(I,ID2) = SUM
   30 CONTINUE
      DO 45 I = 1,ID
         J1 = I+1
         DO 40 J = J1,ID1
            SUM = DZERO
            DO 35 K = J,ID1
               SUM = SUM+V(K)*T(I,K)*T(J,K)
   35       CONTINUE
            T(J,I) = SUM
   40    CONTINUE
   45 CONTINUE
C                                  CALCULATE COEFFICIENTS OF THE
C                                  PARTIALLY DECODED REGRESSION
C                                  COEFFICIENT ESTIMATES, USED TO
C                                  OBTAIN COMPLETELY DECODED ESTIMATES
      BETA = SC(IONE)
      SSS = ONE
      SS = SC(ITWO)
      S = ONE
      SL = ONE
      SKL = ONE
      DO 65 J = 1,ID1
         K1 = J
         K = J+1
         T(J,K1) = SSS
         IF (K .GT. ID1) GO TO 60
         IF (J .EQ. 1) GO TO 50
         SKL = J
         SL = ONE
   50    DO 55 I = K,ID1
            S = SS*S*SKL/SL
            K1 = K1+1
            T(J,K1) = S
            IF (J .EQ. 1) GO TO 55
            SKL = SKL+ONE
            SL = SL+ONE
   55    CONTINUE
   60    SSS = SSS*BETA
         S = SSS
   65 CONTINUE
C                                  CALCULATE VARIANCES OF COMPLETELY
C                                  DECODED REGRESSION COEFFICIENTS
      DO 75 I = 1,ID1
         SUM = DZERO
         DO 70 J = I,ID1
            SUM = SUM+T(J,ID2)*T(I,J)*T(I,J)
   70    CONTINUE
         T(I,ID3) = SUM
   75 CONTINUE
      DO 85 I = 1,ID
         SUM = DZERO
         DO  80 L = I,ID
            K = L + 1
            DO  80 J = K,ID1
               SUM = SUM+T(I,L)*T(I,J)*T(J,L)
   80    CONTINUE
         V(I) = T(I,ID3)+SUM+SUM
   85 CONTINUE
      V(ID1) = T(ID1,ID3)
      GO TO 9005
   90 IER = 129
 9000 CONTINUE
      CALL UERTST (IER,6HRLDCVA)
 9005 RETURN
      END
