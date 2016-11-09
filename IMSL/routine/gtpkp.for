C   IMSL ROUTINE NAME   - GTPKP
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - PROBABILITY DISTRIBUTION OF N ELEMENTS INTO
C                           TWO EQUI-PROBABLE STATES
C
C   USAGE               - CALL GTPKP (N,P,IER)
C
C   ARGUMENTS    N      - INPUT NUMBER OF ELEMENTS. N MUST EXCEED 1.
C                P      - OUTPUT VECTOR OF LENGTH K+1 WHERE K IS THE
C                           GREATEST INTEGER IN N/2.  THE ITH ELEMENT
C                           OF P CONTAINS THE PROBABILITY THAT I-1
C                           ELEMENTS OF THE N ELEMENTS ARE IN THE SAME
C                           STATE, AND N-(I-1) ARE IN THE OTHER STATE.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER = 129 INDICATES THAT N IS LESS THAN 2.
C
C   REQD. IMSL ROUTINES - UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS      GTPKP IS USED IN PREPARATION OF EXPECTED VALUES FOR
C                THE POKER TEST.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE GTPKP(N,P,IER)
C
      DIMENSION          P(1)
      DOUBLE PRECISION   FN,FI,ONE,TEMP
      DATA               ONE/1.0D0/
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
      IF(N.GE.2)GO TO 5
C                                  TERMINAL ERROR
      IER=129
      GO TO 9000
    5 FN = N
C                                  P(1) = (.5)**N
      TEMP=.6931471805599453D0
      TEMP=DEXP(-FN*TEMP)
      P(1)=TEMP
      K = .5*N + .001
C                                  P(I), USING COMBINATORIAL FORMULA
      DO 10 I=1,K
         FI = I
         TEMP=(TEMP*(FN-FI+ONE))/FI
         P(I+1)=TEMP
   10 CONTINUE
C                                  IS N EVEN
      IF(K+K.EQ.N) GO TO 15
      K=K+1
C                                  DOUBLE THE PROBABILITIES
   15 DO 20 I = 1,K
   20 P(I)=P(I)+P(I)
      GO TO 9005
 9000 CONTINUE
      CALL UERTST (IER,'GTPKP ')
 9005 RETURN
      END
