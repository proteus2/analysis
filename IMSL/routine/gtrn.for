C   IMSL ROUTINE NAME   - GTRN
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1980
C
C   PURPOSE             - RUNS TEST
C
C   USAGE               - CALL GTRN (RUNS,N,E,CS,Q,STD,IER)
C
C   ARGUMENTS    RUNS   - INPUT VECTOR OF LENGTH 8 CONTAINING THE TAL-
C                           LIES OF THE OBSERVED RUNS (CALCULATED BY
C                           IMSL ROUTINE GTRTN)
C                N      - INPUT TOTAL LENGTH OF THE SEQUENCE OF RANDOM
C                           NUMBERS
C                E      - OUTPUT VECTOR OF LENGTH 8. E(M) CONTAINS THE
C                           EXPECTED NUMBER OF RUNS OF LENGTH M WHERE
C                           M IS IN THE INCLUSIVE RANGE (1,7). E(8)
C                           CONTAINS THE EXPECTED NUMBER OF RUNS OF
C                           LENGTH MORE THAN 7.
C                CS     - OUTPUT APPROXIMATE CHI-SQUARE STATISTIC
C                Q      - OUTPUT PROBABILITY OF EXCEEDING CS IF THE
C                           UNIFORMITY NULL HYPOTHESIS IS TRUE
C                STD    - OUTPUT STANDARDIZED APPROXIMATE CHI-SQUARE
C                           STATISTIC.  STD=(CS-DF)/SQRT(DF+DF) WHERE
C                           DF IS THE DEGREES OF FREEDOM
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER = 129 INDICATES AN ERROR OCCURRED IN
C                             IMSL ROUTINE MDCH.
C                         WARNING (WITH FIX)
C                         NOTE- IN THE FOLLOWING, X IS THE VECTOR OF
C                         ADJUSTED EXPECTED VALUES OF RUNS.
C                           IER = 66 INDICATES THAT TWO X(I) WERE IN
C                             THE RANGE (1,5). THOSE TWO CELLS WERE
C                             POOLED AND 6 DEGREES OF FREEDOM WERE USED
C                             IN THE CHI-SQUARE TEST. SET RUNS(I)=-1,
C                             WHERE I IS SET FOR THE ONE WHICH WAS
C                             POOLED.
C                           IER = 67 INDICATES THAT THREE X(I) WERE
C                             IN THE RANGE (1,5). THE TWO SMALLEST CELLS
C                             WERE POOLED AND 6 DEGREES OF FREEDOM WERE
C                             USED IN THE CHI-SQUARE TEST. SET RUNS(I)
C                             =-1, WHERE I IS SET FOR THE ONE WHICH WAS
C                             POOLED. IT IS STILL POSSIBLE THAT TWO
C                             CELLS WERE IN THE RANGE (1,5).
C                         WARNING
C                           IER = 36 INDICATES THAT FOUR OR MORE X(I)
C                             ARE IN THE RANGE (1,5). CHI-SQUARE IS
C                             CALCULATED. Q AND THE STANDARDIZED CHI-
C                             SQUARE ARE SET TO MACHINE INFINITY.
C                           IER = 37 INDICATES THAT X(I) LESS THAN 1
C                             APPEARED. CHI-SQUARE IS CALCULATED. Q AND
C                             THE STANDARDIZED CHI-SQUARE ARE SET TO
C                             MACHINE INFINITY.
C                           IER = 38 INDICATES THAT ONE X(I) IS IN
C                             (1,5). PROCESSING CONTINUES.
C
C   REQD. IMSL ROUTINES - H32/MDCH,MDNOR,MERRC=ERFC,MGAMAD=DGAMMA
C                           UERTST,UGETIO
C                       - H36,H48,H60/MDCH,MDNOR,MERRC=ERFC,MGAMA=GAMMA,
C                           UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE GTRN   (RUNS,N,E,CS,Q,STD,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N,IER
      REAL               RUNS(8),E(8),CS,Q,STD
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            II(3),N2,I,K,J,L,JER
      REAL               X(8),RINFP,XX,XN,SR,Z,DF,DQ,CON1,THIRD
      DATA               CON1/.4166667/
      DATA               THIRD/.3333333/
      DATA               RINFP/Z7FFFFFFF/
C                                  FIRST EXECUTABLE STATEMENT
      IER=0
      STD=CON1*(N-3)
      E(1)=STD
      XX=11.
      XN=8.
      SR=RUNS(1)
      CS=60.
      Q=6.
C                                  CALCULATE EXPECTED NUMBER OF RUNS OF
C                                  LENGTH I
      N2=N-2
         DO 5 I=2,7
         Z=XX/CS*(N2-I)
         XX=XX+XN
         CS=CS*Q
         XN=XN+2
         Q=Q+1
         STD=STD+Z
         SR=SR+RUNS(I)
    5    E(I)=Z
      DF=(N+N-7.)*THIRD
      E(8)=DF-STD
      SR=(SR+RUNS(8))/DF
      DF=7
      K=0
         DO 15 I=1,8
         Z=E(I)*SR
         IF (Z .GE. 5) GO TO 15
         IF (Z .LT. 1) GO TO 10
         K=K+1
         IF (K .GT. 4) K = 4
         IF (K .LE. 3) II(K) = I
         GO TO 15
   10    IER=37
         Q=RINFP
         STD=RINFP
   15    X(I)=Z
      IF (K .EQ. 0 .OR. IER .EQ. 37) GO TO 35
      IF (K .EQ. 1) K = 6
      IER = 32+K
      Q=RINFP
      STD=RINFP
      IF (K .NE. 2) GO TO 25
      I=II(1)
      J=II(2)
   20 RUNS(I)=RUNS(I)+RUNS(J)
      RUNS(J)=-1
      X(I) = X(I)+X(J)
      DF=6.
      IER = IER+32
      GO TO 35
   25 IF (K .NE. 3) GO TO 35
      I=II(1)
      J=II(2)
      L=II(3)
      IF (X(I) .GE. X(L)) GO TO 30
      I=II(1)
      GO TO 20
   30 I=II(3)
      GO TO 20
C                                  CALCULATE CHI-SQUARE STATISTIC
   35 CS=0.
         DO 40 I=1,8
         IF (RUNS(I) .LT. 0) GO TO 40
         Z=RUNS(I)-X(I)
         CS=CS+Z*Z/X(I)
   40    CONTINUE
C                                  CALCULATE PROBABILITY
      IF (IER .EQ. 36 .OR. IER .EQ. 37) GO TO 9000
      CALL MDCH (CS,DF,DQ,JER)
      Q=1.0-DQ
      IF (JER .GE. 128) IER = JER
C                                  CALCULATE STANDARDIZED CHI-SQUARE
C                                  STATISTIC
      STD=(CS-DF)/SQRT(DF+DF)
      IF (IER .EQ. 0) GO TO 9005
 9000 CONTINUE
      CALL UERTST (IER,'GTRN  ')
 9005 RETURN
      END
