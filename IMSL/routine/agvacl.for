C   IMSL ROUTINE NAME   - AGVACL
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - ONE OR TWO-SIDED INTERVAL ESTIMATE OF A
C                           VARIANCE COMPONENT
C
C   USAGE               - CALL AGVACL (V,FDF,S,CF,IOP,STAT,IER)
C
C   ARGUMENTS    V      - INPUT VECTOR OF LENGTH 2 CONTAINING THE MEAN
C                           SQUARES INVOLVED IN ESTIMATING THE VARIANCE
C                           COMPONENT. V(1) MUST BE GREATER THAN V(2).
C                FDF    - INPUT VECTOR OF LENGTH 2 CONTAINING THE
C                           DEGREES OF FREEDOM CORRESPONDING TO THE MEAN
C                           SQUARES IN V. THE ORDER OF THE ELEMENTS OF
C                           FDF CORRESPONDS TO THE ORDERING OF V.
C                S      - INPUT ESTIMATE OF THE VARIANCE COMPONENT.
C                           S IS DEFINED TO BE (V(1)-V(2))/C, WHERE C
C                           IS SOME POSITIVE CONSTANT.
C                CF     - INPUT CONFIDENCE COEFFICIENT FOR THE INTERVAL
C                           ESTIMATE (0.95 IS A COMMON CHOICE).
C                IOP    - INPUT OPTION PARAMETER INDICATING THE TYPE OF
C                           INTERVAL ESTIMATE TO BE COMPUTED.
C                           IOP EQUAL TO 1 IMPLIES A TWO-SIDED ESTIMATE
C                           IOP EQUAL TO 2 IMPLIES A UPPER ONE-SIDED
C                             ESTIMATE
C                           IOP NOT EQUAL TO 1 OR 2 IMPLIES A LOWER
C                             ONE-SIDED ESTIMATE
C                STAT   - OUTPUT VECTOR OF LENGTH 2 CONTAINING THE
C                           INTERVAL ESTIMATE.
C                           STAT(1) CONTAINS LOWER LIMIT (NOT DEFINED
C                             WHEN IOP IS EQUAL TO 2)
C                           STAT(2) CONTAINS UPPER LIMIT (DEFINED ONLY
C                             WHEN IOP IS EQUAL TO 1 OR 2)
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         WARNING ERROR WITH FIX
C                           IER=65 MEANS THE COMPUTED LOWER LIMIT WAS
C                             NEGATIVE AND WAS SET TO ZERO.
C                         TERMINAL ERROR
C                           IER=130 MEANS V(1) IS LESS THAN OR EQUAL TO
C                             V(2).  THE INTERVAL ESTIMATE IS NOT
C                             RETURNED.
C                           IER=131 MEANS AN ERROR OCCURRED IN IMSL
C                             ROUTINE MDFI.
C
C   REQD. IMSL ROUTINES - H32/MDBETA,MDBETI,MDFI,MLGAMD=DLGAMA,
C                           UERTST,UGETIO
C                       - H36,H48,H60/MDBETA,MDBETI,MDFI,
C                           MLGAMA=ALGAMA,UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE AGVACL (V,FDF,S,CF,IOP,STAT,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            IOP,IER
      REAL               V(2),FDF(2),S,CF,STAT(2)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            JER
      REAL               TEMP,F1,F2
      REAL               F0,ZERO,HALF,ONE
      DATA               ZERO,HALF,ONE/0.0,0.5,1.0/
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
    5 STAT(1) = ZERO
      STAT(2) = ZERO
      IF (V(1) .GT. V(2)) GO TO 10
C                                  TERMINAL ERROR - V(1).LE.V(2)
      IER = 130
      GO TO 9000
C                                  COMPUTE F0
   10 F0 = V(1)/V(2)
C                                  COMPUTE F1
      IF (IOP .LT. 1 .OR. IOP .GT. 3) IOP = 3
      TEMP = ONE-CF
      IF (IOP .EQ. 1) TEMP = TEMP*HALF
      IF (IOP .EQ. 3) GO TO 15
      CALL MDFI (TEMP,FDF(1),FDF(2),F1,JER)
      IF (JER .NE. 0) GO TO 25
C                                  COMPUTE UPPER LIMIT
      STAT(2) = ((F0-F1)/((F0-ONE)*F1)) * S
      IF (IOP .EQ. 2) GO TO 9005
C                                  COMPUTE F2
      TEMP = (ONE+CF)*HALF
      GO TO 20
   15 TEMP = CF
   20 CALL MDFI (TEMP,FDF(1),FDF(2),F2,JER)
      IF (JER .NE. 0) GO TO 25
C                                  COMPUTE LOWER LIMIT
      STAT(1) = ((F0-F2)/((F0-ONE)*F2)) * S
      IF (STAT(1) .GE. ZERO) GO TO 9005
      IER = 65
      STAT(1) = ZERO
      GO TO 9000
C                                  TERMINAL ERROR RETURNED FROM MDFI
   25 IER = 131
 9000 CONTINUE
      CALL UERTST (IER,'AGVACL')
 9005 RETURN
      END
