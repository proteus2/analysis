C   IMSL ROUTINE NAME   - BECVL
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - VARIANCES AND COVARIANCES OF LINEAR FUNCTIONS
C                           (OUT-OF-CORE VERSION)
C
C   USAGE               - CALL BECVL (X,M,Y,C,R,V)
C
C   ARGUMENTS    X      - INPUT VECTOR OF LENGTH M CONTAINING
C                           COEFFICIENTS OF THE FIRST LINEAR FUNCTION.
C                M      - NUMBER OF COEFFICIENTS. (INPUT)
C                Y      - INPUT VECTOR OF LENGTH M CONTAINING
C                           COEFFICIENTS OF THE SECOND LINEAR FUNCTION.
C                C      - INPUT VARIANCE-COVARIANCE MATRIX OF THE
C                           RANDOM VARIABLES IN THE LINEAR FUNCTION(S).
C                           C IS A SYMMETRIC M X M MATRIX STORED IN
C                           SYMMETRIC MODE. C REQUIRES M*(M+1)/2
C                           STORAGE LOCATIONS.
C                R      - WORK STORAGE VECTOR OF LENGTH M. R MUST BE
C                           TYPED DOUBLE PRECISION.
C                V      - VARIANCE (OR COVARIANCE) OF INPUT LINEAR
C                           FUNCTION (OR FUNCTIONS). (OUTPUT)
C
C   REQD. IMSL ROUTINES - SINGLE/NONE REQUIRED
C                       - DOUBLE/VXADD,VXMUL,VXSTO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE BECVL  (X,M,Y,C,R,V)
C
      REAL               X(1),Y(1),C(1),V
      DOUBLE PRECISION   TEMP,R(1)
C                                  FIRST EXECUTABLE STATEMENT
C                                  MULTIPLY X TIMES C
      LI=1
         DO 15 J=1,M
         LS=LI
         TEMP=0.
            DO 10 K=1,M
            TEMP=TEMP+DBLE(X(K))*DBLE(C(LS))
            IF (K .GE. J) GO TO 5
            LS=LS+1
            GO TO 10
    5       LS=LS+K
   10       CONTINUE
         R(J)=TEMP
   15    LI=LI+J
C                                  MULTIPLY R TIMES Y
      TEMP=0.
         DO 20 K=1,M
   20    TEMP=TEMP+R(K)*Y(K)
      V=TEMP
      RETURN
      END
