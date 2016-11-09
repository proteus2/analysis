C   IMSL ROUTINE NAME   - NAK1
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - KRUSKAL-WALLIS TEST FOR IDENTICAL POPULATIONS
C
C   USAGE               - CALL NAK1 (X,NI,M,EPS,IR,R,STAT,IER)
C
C   ARGUMENTS    X      - INPUT/OUTPUT VECTOR OF LENGTH NI(1)+...+NI(M).
C                         ON INPUT, X CONTAINS THE GROUPED OBSERVATIONS.
C                           X(1),...,X(NI(1)) CONTAINS THE SAMPLE FROM
C                           POPULATION 1. THE OBSERVATIONS FROM
C                           POPULATION I FOLLOW THOSE FOR POPULATION I-1
C                           IN THE X VECTOR, FOR I=2,...,M.
C                         ON OUTPUT, X CONTAINS THE RANKS OF THE
C                           ORIGINAL OBSERVATIONS.
C                NI     - INPUT VECTOR OF LENGTH M CONTAINING THE NUMBER
C                           OF OBSERVATIONS IN EACH GROUP.
C                           NI(I) CONTAINS THE NUMBER OF OBSERVATIONS IN
C                           GROUP I, FOR I=1,...,M.
C                           EACH ELEMENT OF NI MUST BE GREATER THAN
C                           ZERO.
C                M      - INPUT NUMBER OF GROUPS.
C                           M MUST BE GREATER THAN OR EQUAL TO 2.
C                EPS    - INPUT EPSILON USED TO DETERMINE TIES IN X. IF
C                           (AFTER SORTING) ABS(X(I)-X(I+1))-EPS IS
C                           NEGATIVE, A TIE IS COUNTED.
C                IR     - WORK VECTOR OF LENGTH NI(1)+...+NI(M).
C                R      - WORK VECTOR OF LENGTH NI(1)+...+NI(M).
C                STAT   - OUTPUT VECTOR OF LENGTH 4 CONTAINING THE
C                           RESULTING STATISTICS.
C                         STAT(1) CONTAINS THE KRUSKAL-WALLIS H
C                           STATISTIC
C                         STAT(2) CONTAINS THE PROBABILITY OF EXCEEDING
C                           H IF THE NULL HYPOTHESIS OF IDENTICAL
C                           POPULATIONS IS TRUE
C                         STAT(3) CONTAINS H CORRECTED FOR TIES, H*
C                         STAT(4) CONTAINS THE PROBABILITY OF EXCEEDING
C                           H* IF THE NULL HYPOTHESIS OF IDENTICAL
C                           POPULATIONS IS TRUE
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER=129 INDICATES SOME ELEMENT OF NI IS LESS
C                             THAN OR EQUAL TO ZERO
C                           IER=130 INDICATES M IS LESS THAN 2
C                         WARNING ERROR
C                           IER=35 INDICATES THE CHI-SQUARE DEGREES OF
C                             FREEDOM ARE LESS THAN 5 (AND THE BETA
C                             APPROXIMATION IS USED)
C                           IER=36 INDICATES AT LEAST 1 TIE WAS DETECTED
C
C   REQD. IMSL ROUTINES - H32/MDBETA,MDCH,MDNOR,MERRC=ERFC,
C                           MGAMAD=DGAMMA,MLGAMD=DLGAMA,NNUC,UERTST,
C                           UGETIO,VSRTR
C                       - H36,H48,H60/MDBETA,MDCH,MDNOR,MERRC=ERFC,
C                           MGAMA=GAMMA,MLGAMA=ALGAMA,NNUC,UERTST,
C                           UGETIO,VSRTR
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE NAK1   (X,NI,M,EPS,IR,R,STAT,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            NI(1),M,IR(1),IER
      REAL               X(1),EPS,R(1),STAT(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,IER2,J,N,NEL,NHI,NSUM
      REAL               ACH,D,EM,FM1,F1,F2,Q,RSUM,S,SUMRI,T,XN,XSUMRI
C                                  FIRST EXECUTABLE STATEMENT
      IF (M .GE. 2) GO TO 5
      IER = 130
      GO TO 9000
    5 IER = 0
C                                  ACCUMULATE TOTAL LENGTH OF X-VECTOR
      N = 0
      DO 10 I=1,M
         IF (NI(I) .LT. 1) GO TO 70
         N = N + NI(I)
   10 CONTINUE
C                                  INVOKE NNUC FOR RANKING OF X
      CALL NNUC (X,N,EPS,IR,R,X,S,T)
      IF (S .GT. 0.0) IER = 36
C                                  FIND SQUARE OF RANK SUM OF EACH GROUP
C                                  AND SUM THEM
      SUMRI = 0.0
      NSUM = 0
      DO 20 I=1,M
         NHI = NI(I)
         RSUM = 0.0
         DO 15 J=1,NHI
            NEL = NSUM + J
C                                  RSUM=SUM OF RANKS IN I-TH GROUP
            RSUM = RSUM + X(NEL)
   15    CONTINUE
C                                  RSUM=RSUM X AVG RANK OF I-TH GROUP
         RSUM = RSUM * RSUM / NHI
C                                  SUMRI = TOTAL RSUM OVER ALL M GROUPS
         SUMRI = SUMRI + RSUM
         NSUM = NSUM + NHI
   20 CONTINUE
C                                  COMPUTE H=STAT(1)
      XN = N
      XSUMRI = SUMRI
      STAT(1) = 12.0*XSUMRI/(XN*(XN+1.0)) - 3.0*(XN+1.0)
      IF (S .EQ. 0.0) GO TO 25
C                                  COMPUTE H*=STAT(3)
      T = 1.0 - T/(XN*(XN*XN-1.0))
      STAT(3) = STAT(1)/T
      GO TO 30
C                                  H*=H WHEN TIES NONEXISTENT
   25 STAT(3) = STAT(1)
   30 I = 2
      ACH = STAT(1)
   35 IF (M .GE. 6) GO TO 50
C                                  BETA APPROX. FOR M LESS THAN 6
      IF (I .EQ. 4) GO TO 45
C                                  TRANSFORM VARIABLES FOR BETA PDF
      IER = 35
      XSUMRI = 0.0
      DO 40 J=1,M
         XSUMRI = XSUMRI + 1.0/NI(J)
   40 CONTINUE
      EM = M
      XSUMRI = XSUMRI - EM**2/XN
      Q = 0.5 * XN * (XN+1.0)
      Q = Q/((EM-1.0) * (XN-EM))
      Q = Q * XSUMRI
      XSUMRI = 1.2 + XN/(1.0-Q)
      XSUMRI = 1.0/XSUMRI
      Q = (XN+1.0)/(XN-1.0)
      D = 0.5 - 0.6*Q*XSUMRI
      F1 = D * (EM-1.0)
      F2 = D * (XN-EM)
   45 D = ACH/(XN-1.0)
C                                  INVOKE BETA PDF
      CALL MDBETA(D,F1,F2,Q,IER2)
      STAT(I) = 1.0 - Q
      GO TO 55
C                                  CHI-SQ PDF FOR  M NOT LESS THAN 6
   50 FM1 = M-1
      CALL MDCH(ACH,FM1,RSUM,IER2)
      STAT(I) = 1. - RSUM
C                                  CHECK FOR FURTHER CALL TO PDF
   55 IF (I .EQ. 4) GO TO 65
      IF (T .NE. 0.0) GO TO 60
      STAT(4) = STAT(2)
      GO TO 65
C                                  SET PARAMETERS FOR H* CALL TO PDF
   60 I = 4
      ACH = STAT(3)
      GO TO 35
   65 IF (IER .EQ. 0) GO TO 9005
      GO TO 9000
   70 IER = 129
 9000 CONTINUE
      CALL UERTST(IER,6HNAK1  )
 9005 RETURN
      END
