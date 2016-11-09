C   IMSL ROUTINE NAME   - OIND
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - NOVEMBER 1, 1984
C
C   PURPOSE             - WILKS TEST FOR THE INDEPENDENCE OF K SETS
C                           OF MULTI-NORMAL VARIATES
C
C   USAGE               - CALL OIND (S,N,IP,K,STAT,WKAREA,IER)
C
C   ARGUMENTS    S      - INPUT SYMMETRIC VARIANCE-COVARIANCE MATRIX
C                           WHICH HAS BEEN PARTITIONED INTO SUBMATRICES
C                           OF ORDER IP(I). S IS M BY M AND IS STORED IN
C                           SYMMETRIC STORAGE MODE IN A VECTOR OF LENGTH
C                           M(M+1)/2 WHERE M IS THE NUMBER OF VARIATES.
C                N      - THE NUMBER OF INDEPENDENT OBSERVATION VECTORS
C                           USED IN CALCULATING S. N MUST BE GREATER
C                           THAN OR EQUAL TO M + 1.
C                IP     - INPUT VECTOR OF LENGTH K IN WHICH THE ITH
C                           ELEMENT GIVES THE ORDER OF MATRIX S(I,I),
C                           THE ITH PARTITION OF S. M, THE ORDER OF S,
C                           IS THE SUM OF THE ELEMENTS OF VECTOR IP.
C                K      - INPUT NUMBER OF PARTITIONS OF S.
C                STAT   - OUTPUT VECTOR OF LENGTH 4 CONTAINING
C                           STATISTICS. STAT(I) CONTAINS, FOR I =
C                           1--THE STATISTIC V.
C                           2--THE CHI-SQUARED STATISTIC ASSOCIATED
C                             WITH V.
C                           3--THE DEGREES OF FREEDOM OF THE
C                             CHI-SQUARED STATISTIC.
C                           4--THE PROBABILITY OF EXCEEDING THE
C                             STATISTIC STAT(2), IF THE HYPOTHESIS
C                             OF INDEPENDENCE IS CORRECT.
C                WKAREA - WORKING VECTOR OF LENGTH L(L+1)/2 WHERE
C                           L=MAXIMUM, OVER I, OF IP(I).
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER = 129 INDICATES INPUT S IS SINGULAR.
C                           IER = 130 INDICATES THAT ONE OF THE SUB-
C                             MATRICES OF S IS SINGULAR.
C                           IER = 131 INDICATES THAT K IS LESS THAN 2.
C                           IER = 132 INDICATES THAT N IS NOT GREATER
C                             THAN OR EQUAL TO M + 1.
C                           IER = 133 INDICATES THAT AN ELEMENT OF IP
C                             IS NOT GREATER THAN OR EQUAL TO 1.
C                         WARNING ERROR
C                           IER = 35 INDICATES THAT AN UNDERFLOW WAS
C                             PRODUCED DURING THE COMPUTATION OF THE
C                             PROBABILITY P.
C
C   REQD. IMSL ROUTINES - SINGLE(H32)/LUDECP,MDCH,MDNOR,MERRC=ERFC,
C                           MGAMAD=DGAMMA,UERTST,UGETIO
C                       - SINGLE(H36,H48,H60)/LUDECP,MDCH,MDNOR,
C                           MERRC=ERFC,MGAMA=GAMMA,UERTST,UGETIO
C                       - DOUBLE/LUDECP,MDCH,MDNOR,MERRC=ERFC,
C                           MGAMAD=DGAMMA,UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE OIND    (S,N,IP,K,STAT,WKAREA,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N,K,IP(K),IER
      REAL               S(1),STAT(1),WKAREA(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            JER,M,I,ISUM,IND,KI,II,L,KK,J,INDX,NM1,MM
      REAL               STAT2,STAT3,P
      REAL               D1,D2,ALN2,ZERO,ONE,DS,SUMD,E2,E3,CINVS
      DATA               ALN2/.6931472/
      DATA               ZERO/0.0/,ONE/1.0/
C                                  INITIALIZE IER
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
      JER = 0
      IF (K.GE.2) GO TO 5
      IER=131
      GO TO 9000
    5 E2 = ZERO
      E3 = ZERO
      M = 0
C                                  SUM THE IP ARRAY TO GET THE ORDER OF
C                                  MATRIX S.  CHECK FOR IP(I) GE 1
      DO 15 I=1,K
         IF (IP(I) .GE. 1) GO TO 10
            IER = 133
            GO TO 9000
   10    CONTINUE
         M = M + IP(I)
         V = IP(I) * IP(I)
         E2 = E2 + V
         E3 = E3 + V * IP(I)
   15 CONTINUE
C                                  CHECK IF N GE M + 1
      IF (N .GE. M+1) GO TO 20
      IER = 132
      GO TO 9000
   20 CONTINUE
C                                  EVALUATE THE DETERMINANTS
C                                  OF SUB MATRICES OF S
      ISUM = 0
      SUMD = 0
      DO 35 I=1,K
         IND=ISUM
         KI=0
         II=IP(I)
         DO 30 L=1,II
            IND = IND + 1
            KK = IND *(IND+1)/2
            DO 25 J=1,L
               INDX=KK+J-L
               KI=KI+1
               WKAREA(KI)=S(INDX)
   25       CONTINUE
   30    CONTINUE
         CALL LUDECP (WKAREA(1),WKAREA(1),II,D1,D2,JER)
         IF(JER.NE.0) GO TO 50
         SUMD = SUMD + ALOG(D1)+D2*ALN2
         ISUM=ISUM+II
   35 CONTINUE
C                                  EVALUATE THE DETERMINANT
C                                  OF S
      CALL LUDECP (S,S,M,D1,D2,JER)
      IF (JER.NE.0) GO TO 40
      DS = ALOG(D1)+D2*ALN2
      DS = DS-SUMD
      GO TO 45
C                                  S IS SINGULAR, ROUTINE EXITS
   40 IER = 129
      GO TO 9000
C                                  COMPUTE THE STATISTIC V
   45 STAT(1) = EXP(DS)
      NM1 = N-1
      MM = M*M
      E2 = MM-E2
      E3 = MM*M-E3
      CINVS = ONE-ONE/(6.0*E2*NM1)*(E3+E3+E2+E2+E2)
C                                  COMPUTE THE CHI-SQUARED
C                                  STATISTIC ASSOCIATED WITH V
      STAT(2) = -(NM1)*CINVS*DS
C                                  COMPUTE THE DEGREE OF
C                                  FREEDOM OF THE CHI-SQUARED STATISTIC
      STAT(3) = .5*E2
      STAT2 = STAT(2)
      STAT3 = STAT(3)
      CALL MDCH (STAT2,STAT3,P,JER)
      IF(JER.EQ.34) IER = 35
      STAT(4) = ONE - P
      IF(IER.EQ.0) GO TO 9005
      GO TO 9000
   50 IER=130
 9000 CONTINUE
      CALL UERTST (IER,6HOIND  )
 9005 RETURN
      END
