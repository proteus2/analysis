C   IMSL ROUTINE NAME   - BECVLI
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - VARIANCES AND COVARIANCES OF LINEAR FUNCTIONS
C                           (IN-CORE VERSION)
C
C   USAGE               - CALL BECVLI (X,N,M,IX,C,R,IOPT,V)
C
C   ARGUMENTS    X      - N BY M MATRIX OF LINEAR FUNCTION COEFFICIENTS.
C                           (INPUT)
C                N      - NUMBER OF ROWS IN X. (INPUT)
C                M      - NUMBER OF COLUMNS IN X. (INPUT)
C                IX     - ROW DIMENSION OF MATRIX X EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM. (INPUT)
C                C      - INPUT SYMMETRIC VARIANCE-COVARIANCE MATRIX OF
C                           THE RANDOM VARIABLES IN THE LINEAR
C                           FUNCTIONS. C MUST BE AN M X M SYMMETRIC
C                           MATRIX STORED IN SYMMETRIC STORAGE MODE.
C                           THEREFORE, A VECTOR WITH M*(M+1)/2
C                           LOCATIONS IS NEEDED.
C                R      - MATRIX OR VECTOR CONTAINING N*M LOCATIONS
C                           USED FOR WORK STORAGE. R MUST BE DOUBLE
C                           PRECISION.
C                IOPT   - INPUT COVARIANCE OPTION PARAMETER.
C                           IF IOPT=1, COVARIANCES ARE COMPUTED.
C                           OTHERWISE, THEY ARE NOT.
C                V      - OUTPUT MATRIX OF VARIANCES AND COVARIANCES.
C                           V IS AN N X N SYMMETRIC MATRIX AND IS
C                           STORED IN SYMMETRIC MODE. THEREFORE, A
C                           VECTOR V WITH N*(N+1)/2 LOCATIONS IS NEEDED.
C                           HOWEVER, IF ONLY THE VARIANCES ARE DESIRED
C                           (IOPT .NE. 1), THEN ONLY A VECTOR WITH N
C                           LOCATIONS IS REQUIRED.
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
      SUBROUTINE BECVLI (X,N,M,IX,C,R,IOPT,V)
C
      REAL               X(IX,1),C(1),V(1)
      DOUBLE PRECISION   TEMP,R(N,1)
C                                  FIRST EXECUTABLE STATEMENT
C                                  COMPUTE X TIMES C
         DO 15 I=1,N
         LI=1
            DO 15 J=1,M
            LS=LI
            TEMP=0.
               DO 10 K=1,M
               TEMP=TEMP+DBLE(X(I,K))*DBLE(C(LS))
               IF (K .GE. J) GO TO 5
               LS=LS+1
               GO TO 10
    5          LS=LS+K
   10          CONTINUE
            R(I,J)=TEMP
   15       LI=LI+J
C                                  MULTIPLY R TIMES ALL X
      L=1
      JJ=1
         DO 25 I=1,N
         IF (IOPT .NE. 1) JJ=I
            DO 25 J=JJ,I
            TEMP=0.
               DO 20 K=1,M
   20          TEMP=TEMP+R(I,K)*X(J,K)
            V(L)=TEMP
   25       L=L+1
      RETURN
      END
