C   IMSL ROUTINE NAME   - NMCC
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - CALCULATE AND TEST THE SIGNIFICANCE OF THE
C                           KENDALL COEFFICIENT OF CONCORDANCE
C
C   USAGE               - CALL NMCC (X,M,N,IX,EPS,IR,R,STAT,IER)
C
C   ARGUMENTS    X      - INPUT/OUTPUT M BY N MATRIX.
C                         ON INPUT, X IS A MATRIX OF N SETS OF M
C                           OBSERVATIONS OR X CAN BE N SETS OF M
C                           RANKS OF OBSERVATIONS.
C                         ON OUTPUT, X CONTAINS THE RANKS IN EITHER
C                           CASE.
C                M      - INPUT NUMBER OF OBSERVATIONS OR RANKS PER SET.
C                           M MUST BE GREATER THAN OR EQUAL TO 2.
C                N      - INPUT NUMBER OF SETS OF OBSERVATIONS OR
C                           RANKINGS. M MUST BE GREATER THAN OR EQUAL
C                           TO 2.
C                IX     - INPUT ROW DIMENSION OF THE MATRIX X EXACTLY
C                           AS SPECIFIED IN THE DIMENSION STATEMENT IN
C                           THE CALLING PROGRAM.
C                EPS    - INPUT NUMBER TO BE USED TO DETERMINE TIES
C                           IN RANKING. IF, WITHIN A SET OF OBSERVATIONS
C                           (COLUMN OF X), THE DIFFERENCE BETWEEN TWO
C                           ELEMENTS IS LESS THAN OR EQUAL TO EPS,
C                           A TIE IS COUNTED.
C                IR     - WORK VECTOR OF LENGTH M.
C                R      - WORK VECTOR OF LENGTH M. ON OUTPUT, R CONTAINS
C                           THE SUMS OF THE RANKS WHICH ARE IN X.
C                STAT   - OUTPUT VECTOR OF LENGTH 4.
C                         STAT(1) CONTAINS W, THE COEFFICIENT OF
C                           CONCORDANCE
C                         STAT(2) CONTAINS THE RESULTING CHI-SQUARED
C                           STATISTIC, WITH M-1 DEGREES OF FREEDOM.
C                         STAT(3) CONTAINS THE PROBABILITY OF EXCEEDING
C                           STAT(2) IF THE NULL HYPOTHESIS OF
C                           INDEPENDENCE IS CORRECT.
C                         STAT(4) CONTAINS THE KENDALL S - SUMS OF
C                           SQUARES OF DEVIATIONS FROM EXPECTED SUMS OF
C                           RANKS.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER=129 MEANS THAT X WAS NOT AT LEAST 2 BY 2
C                         WARNING ERROR
C                           IER=34 MEANS THE CHI-SQUARED DEGREES OF
C                             FREEDOM IS LESS THAN 7 AND STAT(3) IS TO
C                             BE REGARDED WITH SUSPICION.
C
C   REQD. IMSL ROUTINES - H32/MDCH,MDNOR,MERRC=ERFC,MGAMAD=DGAMMA,NNUC,
C                           UERTST,UGETIO,VSRTR
C                       - H36,H48,H60/MDCH,MDNOR,MERRC=ERFC,MGAMA=GAMMA,
C                           NNUC,UERTST,UGETIO,VSRTR
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE NMCC   (X,M,N,IX,EPS,IR,R,STAT,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            M,N,IX,IR(1),IER
      REAL               X(IX,1),EPS,R(1),STAT(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,IET,J
      REAL               FM1,P,Q,S,SB,SUM,T,T1,W,XSQRD,PZET
      DATA               PZET/.8333333E-1/
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
      T1 = 0.0
C                                  RANK THE COLUMNS OF X.
      DO 5 I=1,N
         CALL NNUC (X(1,I),M,EPS,IR,R,X(1,I),S,T)
C                                  SUM THE TIE STATISTICS
    5 T1 = T1+T
      T1 = T1*PZET
      DO 10 I=1,M
   10 R(I) = 0.0
C                                  SUM THE RANKS INTO R
      S = 0.0
      DO 20  J=1,M
         DO 15 I=1,N
   15    R(J) = R(J) + X(J,I)
   20 S = S+R(J)
      SB = S/M
      S = 0.0
C                                  OBTAIN THE STATISTIC S
      DO 25 I=1,M
         SUM = R(I)-SB
   25 S = SUM*SUM+S
      P = M*N*(M+1)*PZET
C                                  AND THE KENDALL W
      FM1 = M-1.
      W = S/(P*FM1*N-N*T1)
C                                  TRANSFORM S TO A CHI-SQUARED DEVIATE
C                                    XSQRD
      XSQRD = S/(P-T1/FM1)
      IF (FM1 .LT. 7.) IER = 34
      STAT(1) = W
      STAT(2) = XSQRD
      STAT(4) = S
C                                  OBTAIN THE PROBABILITY OF EXCEEDING
C                                    XSQRD
      CALL MDCH(XSQRD,FM1,Q,IET)
      IF (IET.GT.127) GO TO 30
      STAT(3) = 1.-Q
      IF (IER-34) 9005,9000,9005
   30 IER = IET
 9000 CONTINUE
       CALL UERTST(IER,6HNMCC  )
 9005 RETURN
      END
