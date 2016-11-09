C   IMSL ROUTINE NAME   - GTCN
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1980
C
C   PURPOSE             - SAMPLE SIZE OR NUMBER OF CLASS INTERVALS
C                           DETERMINATION FOR CHI-SQUARED TEST
C                           APPLICATIONS
C
C   USAGE               - CALL GTCN (Q,IOPT,B,K,N,IER)
C
C   ARGUMENTS    Q      - INPUT PROBABILITY ASSOCIATED WITH THE
C                           CRITICAL REGION UNDER THE NULL HYPOTHESIS.
C                           Q SHOULD BE GREATER THAN .005 AND LESS
C                           THAN 0.1 . GENERALLY Q=.05 OR Q=.01 IS USED.
C                IOPT   - INPUT OPTION PARAMETER.
C                           IOPT=1 INDICATES N IS INPUT
C                           IOPT=0 INDICATES K IS GIVEN
C                B      - INPUT PARAMETER, ORDINARILY BETWEEN 2 AND 4
C                           INCLUSIVE.  B=4 WILL YEILD AN OPTIMAL VALUE
C                           FOR N OR K FOR THE POINT AT WHICH THE POWER
C                           IS .5.  (SEE THE REFERENCES.)
C                K      - BEST NUMBER OF CLASS INTERVALS, CORRESPOND-
C                           ING TO N, FOR USE IN SUBDIVIDING THE RANGE
C                           OF THE RANDOM NUMBERS TO BE TESTED. K IS
C                           OUTPUT IF IOPT=1, INPUT IF IOPT=0.
C                N      - LENGTH OF THE SEQUENCE OF RANDOM NUMBERS TO
C                           BE TESTED. N IS INPUT IF IOPT=1, OUTPUT IF
C                           IOPT=0.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER = 129 INDICATES Q IS SPECIFIED INCOR-
C                             RECTLY.
C                         WARNING ERROR
C                           IER = 34 INDICATES THAT EITHER N IS LESS
C                             THAN OR EQUAL TO 200 OR K IS LESS THAN
C                             OR EQUAL TO 26.
C
C   REQD. IMSL ROUTINES - MDNRIS,MERFI,UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE GTCN   (Q,IOPT,B,K,N,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            IOPT,K,N,IER
      REAL               Q,B
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      REAL               P,C,D,SSQ1H,CON1
      DATA               SSQ1H/.7071068/
      DATA               CON1/1.148698/
C                                  FIRST EXECUTABLE STATEMENT
      IER=0
      IF (Q .GT. .5E-2 .AND. Q .LT. .1) GO TO 5
C                                  SET TERMINAL ERROR
      IER=129
      GO TO 9000
    5 P=1.-Q
C                                  OBTAIN CRITICAL POINT C
      CALL MDNRIS(P,C,IER)
      IF (IOPT .EQ. 1) GO TO 10
C                                  WARNING ERROR
      IF (K .LE. 26) IER=34
      D=K/B
      N = SQRT(D)*D*D*C*SSQ1H+1
      GO TO 15
C                                  WARNING ERROR
   10 IF (N .LE. 200) IER=34
      D=((N-1)/C)**.4
      K = D*B*CON1
   15 IF(IER.EQ.0) GO TO 9005
 9000 CONTINUE
      CALL UERTST(IER,'HGTCN  ')
 9005 RETURN
      END
