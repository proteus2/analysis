C   IMSL ROUTINE NAME   - RLGQMI
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1980
C
C   PURPOSE             - CENTERING OF INDEPENDENT VARIABLE SETTINGS AND
C                           GENERATION OF CENTERED SQUARE AND CROSS
C                           PRODUCT TERMS - IN CORE VERSION
C
C   USAGE               - CALL RLGQMI (X,N,M,IX,XBAR)
C
C   ARGUMENTS    X      - INPUT AND OUTPUT N BY M*(M+3)/2 MATRIX.
C                         ON INPUT, X IS AN N BY M MATRIX OF N SETTINGS
C                           FOR EACH OF M INDEPENDENT VARIABLES.
C                         ON OUTPUT, X IS AN N BY M*(M+3)/2 MATRIX OF
C                           N SETTINGS FOR EACH OF THE CENTERED M INPUT
C                           VARIABLES AND THE M*(M+1)/2 GENERATED AND
C                           CENTERED SQUARE AND CROSS PRODUCT VARIABLES.
C                N      - INPUT NUMBER OF ROWS IN THE INPUT MATRIX X.
C                M      - INPUT NUMBER OF COLUMNS IN THE INPUT MATRIX X.
C                IX     - INPUT ROW DIMENSION OF THE MATRIX X EXACTLY
C                          AS SPECIFIED IN THE DIMENSION STATEMENT IN
C                          THE CALLING PROGRAM.
C                XBAR   - INPUT AND OUTPUT VECTOR OF LENGTH M*(M+3)/2.
C                         ON INPUT, XBAR IS A VECTOR OF LENGTH M
C                           CONTAINING THE MEANS OF THE M INDEPENDENT
C                           VARIABLES.
C                         ON OUTPUT, XBAR IS A VECTOR OF LENGTH
C                           M*(M+3)/2 WHOSE FIRST M COMPONENTS ARE
C                           UNCHANGED AND WHOSE LAST M*(M+1)/2
C                           COMPONENTS CONTAIN THE MEANS OF THE
C                           GENERATED SQUARE AND CROSS PRODUCT
C                           VARIABLES.
C
C   REQD. IMSL ROUTINES - SINGLE/NONE REQUIRED
C                       - DOUBLE/VXADD,VXMUL,VXSTO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS      THE INPUT VARIABLES MAY BE CENTERED OR UNCENTERED.
C                THE INPUT MEANS MUST BE THE MEANS OF THESE VARIABLES.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE RLGQMI (X,N,M,IX,XBAR)
C                                 SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N,M,IX
      REAL               X(IX,1),XBAR(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,J,JJ,L,L1
      REAL               VMAX
      DOUBLE PRECISION   TEMP,T
C                                  FIRST EXECUTABLE STATEMENT
      T=1.D0/N
      DO 5 J=1,M
         VMAX=XBAR(J)
         IF (VMAX .EQ. 0.) GO TO 10
         DO 5 I=1,N
            X(I,J)=X(I,J)-VMAX
    5 CONTINUE
   10 CONTINUE
      L1=M+M
      DO 30 J=1,M
         DO 30 JJ=J,M
            TEMP=0.0D0
            DO 15 I=1,N
               TEMP=TEMP+DBLE(X(I,J))*DBLE(X(I,JJ))
   15       CONTINUE
            TEMP=TEMP*T
            IF (J .EQ. JJ) GO TO 20
            L1=L1+1
            L=L1
            GO TO 25
   20       L=M+J
   25       XBAR(L)=TEMP
            DO 30 I=1,N
   30          X(I,L)=DBLE(X(I,J))*DBLE(X(I,JJ))-TEMP
      RETURN
      END
