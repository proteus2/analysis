C   IMSL ROUTINE NAME   - GTNOR
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - TEST FOR NORMALITY OF RANDOM DEVIATES
C
C   USAGE               - CALL GTNOR (R,N,K,STAT,OBSC,CSOBS,IER)
C
C   ARGUMENTS    R      - INPUT VECTOR OF NORMAL (0,1) RANDOM DEVIATES.
C                           R IS N LONG.
C                N      - INPUT LENGTH OF VECTOR R
C                K      - INPUT NUMBER OF EQUIPROBABLE CATEGORIES INTO
C                           WHICH THE ELEMENTS OF R ARE TO BE TALLIED.
C                           IF THE TALLY OF THE TOTAL SEQUENCE REQUIRES
C                           SEVERAL INVOCATIONS OF GTNOR, K MUST NOT BE
C                           CHANGED BETWEEN CALLS.
C                STAT   - OUTPUT VECTOR OF LENGTH 3. THE I-TH ELEMENT OF
C                           STAT CONTAINS, WHEN
C                           I=1, THE CHI-SQUARE STATISTIC RESULTING FROM
C                             THE TEST
C                           I=2, THE CHI-SQUARE STATISTIC, STANDARDIZED
C                           I=3, THE PROBABILITY OF EXCEEDING THAT
C                             STATISTIC, GIVEN THE HYPOTHESIS OF
C                             NORMALITY IS TRUE. ON ALL ENTRIES EXCEPT
C                             THE FINAL ENTRY, STAT(3) MUST BE ZERO.
C                             ON THE FINAL ENTRY, STAT(3) MUST BE NON-
C                             ZERO. (SEE REMARKS)
C                OBSC   - OUTPUT VECTOR WHICH CONTAINS TALLY OF OCCUR-
C                           RENCES OF ELEMENTS OF R IN EACH OF K CATEGO-
C                           RIES. OBSC IS K LONG AND MUST BE NULLED BY
C                           THE CALLING PROGRAM PRIOR TO ITS FIRST INVO-
C                           CATION OF GTNOR.
C                CSOBS  - OUTPUT VECTOR WHICH CONTAINS THE K INDIVIDUAL
C                           COMPONENTS OF THE CHI-SQUARE STATISTIC.
C                           CSOBS IS K LONG.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER = 129 INDICATES CHI-SQUARE STATISTIC IS
C                             LESS THAN ZERO
C                           IER = 130 INDICATES K IS LESS THAN 2
C                             (INSUFFICIENT NUMBER OF CATEGORIES)
C                         WARNING ERROR (WITH FIX)
C                           IER = 67 INDICATES EXPECTED VALUE OF OBSC
C                             CELLS IS LESS THAN 1.  STAT(3) WILL BE SET
C                             AT MACHINE INFINITY AND MDCH WILL NOT BE
C                             INVOKED.
C                         WARNING ERROR
C                           IER = 36 INDICATES EXPECTED VALUE OF OBSC
C                             CELLS IS LESS THAN 5.  (SMALLER THAN
C                             ADVISABLE)
C
C   REQD. IMSL ROUTINES - H32/MDCH,MDNOR,MERRC=ERFC,MGAMAD=DGAMMA,
C                           UERTST,UGETIO
C                       - H36,H48,H60/MDCH,MDNOR,MERRC=ERFC,MGAMA=GAMMA,
C                           UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS      THE VECTOR OBSC MUST BE NULLED BY THE CALLING PROGRAM
C                PRIOR TO THE FIRST CALL TO GTNOR. IF THE TOTAL
C                SEQUENCE OF NUMBERS TO BE TESTED CAN BE HELD IN
C                MEMORY, ONE CALL TO GTNOR WILL PERFORM THE TEST AND
C                STAT(3) MUST BE NON-ZERO ON THAT CALL. IF THE TOTAL
C                SEQUENCE MUST BE PARTITIONED AND READ INTO MEMORY IN
C                SECTIONS, THEN GTNOR MUST BE CALLED MORE THAN ONCE.
C                IN THE LATTER CASE, STAT(3) MUST BE ZERO ON ALL
C                CALLS EXCEPT THE FINAL CALL. NEITHER K NOR OBSC
C                MAY BE MODIFIED BETWEEN CALLS.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE GTNOR (R,N,K,STAT,OBSC,CSOBS,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N,K,IER
      REAL               R(1),STAT(1),OBSC(1),CSOBS(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,J,IER2
      REAL               P,E,DF,RINFP,QUE,SSQ1H
C                                  RINFP=LGST POS REAL -SINGLE PRECISION
      DATA               RINFP/Z7FFFFFFF/
      DATA               SSQ1H/.7071068/
C                                  CHECK FOR K LESS THAN 2
C                                  FIRST EXECUTABLE STATEMENT
      IF (K .GE. 2) GO TO 5
      IER = 130
      GO TO 9000
    5 IER = 0
C                                  DETERMINE PROPER CELL FOR EACH R(I)
      DO 10 I=1,N
         CALL MDNOR (R(I),P)
         J = K *P + 1.0
         IF (P .EQ. 1.0) J=K
C                                  INCREMENT TALLY FOR THAT CELL
         OBSC(J) = OBSC(J) + 1.0
   10 CONTINUE
C                                  CHECK FOR LAST ENTRY
      IF (STAT(3) .EQ. 0.0) GO TO 9005
C                                  DETERMINE TOTAL OF ALL TALLIES
      E = 0.0
      DO 15 I=1,K
         E = E + OBSC(I)
   15 CONTINUE
C                                  COMPUTE EXPECTED TALLY
      E = E/K
C                                  COMPUTE CHI-SQ COMPONENTS
      P = 0.0
      DO 20 I=1,K
         CSOBS(I) = (OBSC(I)-E)**2/E
C                                  ACCUMULATE OVERALL CHI-SQ
      P = P + CSOBS(I)
   20 CONTINUE
      STAT(1) = P
      DF = K-1
C                                  COMPUTE STANDARDIZED CHI-SQ
      STAT(2) = (P-DF) * SSQ1H/SQRT(DF)
C                                  CHECK VALUE OF E BEFORE CALL
      IF (E .LT. 5.0) IER=36
      IF (E .GE. 1.0) GO TO 25
      IER = 67
      STAT(3) = RINFP
      GO TO 9000
C                                  INVOKE CHI-SQ ROUTINE
   25 CALL MDCH (P,DF,QUE,IER2)
C                                  CHECK FOR ERROR IN MDCH
      IF (IER2 .LT. 129) GO TO 30
      IER = 129
      GO TO 9000
   30 STAT(3) = 1.0-QUE
      IF (IER .EQ. 0) GO TO 9005
 9000 CONTINUE
      CALL UERTST (IER,'GTNOR ')
 9005 RETURN
      END
