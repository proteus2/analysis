C   IMSL ROUTINE NAME   - BEMNON
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1980
C
C   PURPOSE             - LOCATION (MEAN) INFERENCES USING A SAMPLE
C                           FROM A NORMAL POPULATION WITH KNOWN
C                           VARIANCE
C
C   USAGE               - CALL BEMNON (Y,N,IOP,CRIT,M,YMN,STAT,IER)
C
C   ARGUMENTS    Y      - INPUT VECTOR OF LENGTH N CONTAINING THE
C                           SAMPLE.
C                N      - SAMPLE SIZE. (INPUT)
C                IOP    - INPUT VECTOR OF LENGTH 2 INDICATING, FOR THE
C                           MEAN, THE HYPOTHESIS TEST AND INTERVAL
C                           ESTIMATE TO BE COMPUTED.
C                           IOP(1)=I IMPLIES THE FOLLOWING TEST, WHEN
C                           I=1, TWO-TAILED.
C                           I=2, UPPER ONE-TAILED.
C                           I=3, LOWER ONE-TAILED.
C                           I NOT EQUAL TO 1, 2, OR 3 IMPLIES TEST
C                           IS NOT DESIRED.  SEE REMARKS.
C                           IOP(2)=I IMPLIES THE FOLLOWING INTERVAL
C                           ESTIMATE, WHEN
C                           I=1, TWO-SIDED.
C                           I=2, UPPER ONE-SIDED.
C                           I=3, LOWER ONE-SIDED.
C                           I NOT EQUAL TO 1, 2, OR 3 IMPLIES ESTIMATE
C                           IS NOT DESIRED.
C                CRIT   - INPUT VECTOR OF LENGTH 3 CONTAINING CONSTANTS
C                           REQUIRED FOR INTERVAL ESTIMATE AND HYPOTHE-
C                           SIS TEST. THE I-TH ELEMENT OF CRIT CONTAINS,
C                           WHEN
C                           I=1, CONFIDENCE COEFFICIENT FOR INTERVAL
C                             ESTIMATE OF MEAN. THE CHOICE 0.95 IS A
C                             COMMON ONE. (REQUIRED ONLY WHEN IOP(2)=1,
C                             2, OR 3).
C                           I=2, HYPOTHESIZED VALUE OF THE MEAN (RE-
C                             QUIRED ONLY WHEN IOP(1).= 1, 2, OR 3)
C                           I=3, KNOWN VARIANCE OF POPULATION.
C                M      - NUMBER OF RESPONSES IN EARLIER SAMPLE YIELDING
C                           ADDITIONAL ESTIMATE OF MEAN. M LESS THAN
C                           2 IMPLIES NO ADDITIONAL ESTIMATE IS
C                           AVAILABLE. (INPUT).
C                YMN    - IF M IS GREATER THAN OR EQUAL TO 2, YMN
C                           CONTAINS ADDITIONAL ESTIMATE OF MEAN ON
C                           INPUT AND COMBINED ESTIMATE ON OUTPUT. IF M
C                           IS LESS THAN 2, YMN IS UNDEFINED ON INPUT
C                           AND CONTAINS ESTIMATE OF MEAN ON OUTPUT.
C                STAT   - OUTPUT VECTOR OF LENGTH 3 CONTAINING RESULTS
C                           OF HYPOTHESIS TEST AND INTERVAL ESTIMATE FOR
C                           MEAN. THE I-TH ELEMENT OF STAT CONTAINS,
C                           WHEN
C                           I=1, PROBABILITY OF COMPUTED OR MORE EXTREME
C                             VALUE OF TEST STATISTIC FOR HYPOTHESIS
C                             TEST ON MEAN. (DEFINED ONLY WHEN IOP(1)=1,
C                             2, OR 3).
C                           I=2, LOWER LIMIT FOR MEAN (DEFINED ONLY WHEN
C                             IOP(2)=1 OR 3).
C                           I=3, UPPER LIMIT FOR MEAN (DEFINED ONLY WHEN
C                             IOP(2)=1 OR 2).
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER=129 MEANS NUMBER OF RESPONSES SPECIFIED
C                             IS LESS THAN 1.
C                           IER=130 MEANS AN ERROR OCCURRED IN MDNRIS.
C
C
C   REQD. IMSL ROUTINES - SINGLE/MDNOR,MDNRIS,MERFI,MERRC=ERFC,UERTST,
C                           UGETIO
C                       - DOUBLE/MDNORD,MDNRIS,MERFI,MERRCD=DERFC,
C                           UERTST,UGETIO,VXADD,VXSTO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS      LETTING (M,S) BE THE TRUE MEAN AND KNOWN VARIANCE OF
C                THE POPULATION, RESPECTIVELY, AND DEFINING
C                           M0 = CRIT(2)
C                THE TABLE BELOW CLARIFIES HYPOTHESIS TESTING OPTIONS.
C
C                      IOP(1)            NULL            ALTERNATIVE
C                                     HYPOTHESIS         HYPOTHESIS
C                      ------         ----------         -----------
C                        1             M.EQ.M0             M.NE.M0
C                        2             M.LE.M0             M.GT.M0
C                        3             M.GE.M0             M.LT.M0
C
C-----------------------------------------------------------------------
C
      SUBROUTINE BEMNON (Y,N,IOP,CRIT,M,YMN,STAT,IER)
C
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N,IOP(2),M,IER
      REAL               Y(N),CRIT(3),YMN,STAT(3)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      DOUBLE PRECISION   YC,YCOM
      REAL               SDN,XN,XM,C,Z,ZERO,HALF
      DATA               C/.7071068/
      DATA               ZERO/0.0/,HALF/0.5/
C                                  FIRST EXECUTABLE STATEMENT
      IF (N .GE. 1) GO TO 5
C                                  TERMINAL ERROR - NUMBER OF RESPONSES
C                                  (N) SPECIFIED .LT. 1
      IER = 129
      GO TO 9000
    5 IER = 0
      STAT(1) = ZERO
      STAT(2) = ZERO
      STAT(3) = ZERO
      IOP1 = IOP(1)
      IOP2 = IOP(2)
      XN = N
      XM = M
      SDN = SQRT(CRIT(3)/XN)
C                                  COMPUTE THE MEAN
      YC = 0.D0
      DO 10 I = 1,N
         YC = YC + Y(I)
   10 CONTINUE
      YCOM = YC/XN
      IF (M .LT. 2) GO TO 15
      YMN = (XN*YCOM + XM*YMN)/(XN + XM)
      SDN = SQRT(CRIT(3)/(XN+XM))
      GO TO 20
   15 YMN = YCOM
   20 IF (IOP1 .LT. 1 .OR. IOP1 .GT. 3) GO TO 40
C                                  COMPUTE PROBABILITY FOR HYPOTHESIS
C                                    TEST FOR THE MEAN
      Z  = (YMN - CRIT(2))/SDN
      GO TO (25,30,35), IOP1
   25 IF (Z .LT. ZERO) Z = -Z
      STAT(1) = ERFC (C*Z)
      GO TO 40
   30 CALL MDNOR (-Z,STAT(1))
      GO TO 40
   35 CALL MDNOR (Z,STAT(1))
   40 IF (IOP2 .GT. 3 .OR. IOP2 .LT. 1) GO TO 9005
C                                  COMPUTE INTERVAL FOR MEAN
      P = CRIT(1)
      IF (IOP2 .EQ. 1) P = HALF + HALF * CRIT(1)
      CALL MDNRIS (P,A,IEM)
      IF (IEM .EQ. 0) GO TO 45
C                                  ERROR CONDITION RETURNED FROM MDNRIS
      IER = 130
      GO TO 9000
   45 IF (IOP2 .NE. 2) STAT(2) = YMN - A*SDN
      IF (IOP2 .NE. 3) STAT(3) = YMN + A*SDN
      GO TO 9005
 9000 CONTINUE
      CALL UERTST (IER,'BEMNON')
 9005 RETURN
      END
