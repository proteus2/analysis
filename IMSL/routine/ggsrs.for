C   IMSL ROUTINE NAME   - GGSRS
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1980
C
C   PURPOSE             - GENERATE A SIMPLE RANDOM SAMPLE
C                           FROM A FINITE POPULATION.
C
C   USAGE               - CALL GGSRS (DSEED,IOPT,NPOP,IP,MPOP,POP,
C                           NSAMP,MSAMP,SAMP,IX,IER)
C
C   ARGUMENTS    DSEED  - INPUT/OUTPUT DOUBLE PRECISION VARIABLE
C                           ASSIGNED AN INTEGER VALUE IN THE
C                           EXCLUSIVE RANGE (1.D0, 2147483647.D0).
C                           DSEED IS REPLACED BY A NEW VALUE TO BE
C                           USED IN A SUBSEQUENT CALL.
C                IOPT   - INPUT PARAMETER.
C                           IOPT = 0 INDICATES ONLY A SIMPLE RANDOM
C                           SAMPLE OF INDICES IS TO BE PRODUCED.  IF
C                           IOPT .NE. 0, THEN A SIMPLE RANDOM SAMPLE
C                           OF THE VALUES IN POP WILL BE OUTPUT IN SAMP.
C                NPOP   - INPUT. THE SIZE OF THE POPULATION.
C                IP     - INPUT/OUTPUT. ON THE FIRST CALL TO GGSRS
C                           IP=0.  THE USER DOES NOT CHANGE IP ON
C                           SUBSEQUENT CALLS. FOLLOWING EACH CALL
C                           IP IS THE CUMULATED POPULATION INDEX.
C                MPOP   - INPUT.  IF IOPT .NE. 0, MPOP IS THE SIZE OF
C                           THE SUBPOPULATION INCLUDED IN EACH
C                           CALL (EXCEPT POSSIBLY THE LAST) TO GGSRS.
C                           IF THE FULL POPULATION RESIDES IN MEMORY,
C                           MPOP=NPOP AND ONLY ONE CALL TO GGSRS IS
C                           NECESSARY.  IF IOPT=0 MPOP IS NOT NEEDED.
C                POP    - INPUT VECTOR.  IF IOPT .NE. 0, POP IS OF
C                           LENGTH MPOP AND CONTAINS THE SUBPOPULATION
C                           (OR FULL POPULATION IF MPOP=NPOP).
C                           IF IOPT=0, POP IS NOT USED AND MAY BE
C                           DIMENSIONED AS POP(1).
C                NSAMP  - INPUT. THE SAMPLE SIZE DESIRED.
C                MSAMP  - OUTPUT.  THE SIZE OF THE SAMPLE COLLECTED
C                           UP TO THE POINT OF RETURN FROM GGSRS.
C                           AFTER THE FINAL CALL TO GGSRS, MSAMP=NSAMP.
C                SAMP   - OUTPUT VECTOR.  IF IOPT=0 SAMP IS NOT USED
C                           AND MAY BE DIMENSIONED AS SAMP(1).
C                           IF IOPT .NE. 0 SAMP IS OF LENGTH NSAMP AND
C                           CONTAINS THE SAMPLE VALUES.
C                IX     - OUTPUT VECTOR OF LENGTH NSAMP CONTAINING
C                           THE INDICES OF THE POPULATION WHICH
C                           ARE INCLUDED IN THE SAMPLE.
C                IER    - ERROR PARAMETER.  (OUTPUT)
C                         WARNING ERROR
C                           IER = 33 MEANS THE POPULATION SIZE IS
C                           EQUAL TO OR LESS THAN THE SAMPLE SIZE
C                           DESIRED.  THE SAMPLE WILL CONTAIN THE
C                           ENTIRE POPULATION.
C
C   REQD. IMSL ROUTINES - GGUBS,UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C
C-----------------------------------------------------------------------
C
      SUBROUTINE GGSRS (DSEED,IOPT,NPOP,IP,MPOP,POP,
     1                   NSAMP,MSAMP,SAMP,IX,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            IOPT,NPOP,IP,MPOP,NSAMP,MSAMP,IX(NSAMP),IER,ICP
      REAL               POP(1),SAMP(1),R(1)
      DOUBLE PRECISION   DSEED
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            ILIM
      REAL               TEMP
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
      IF (NPOP .LE. NSAMP) IER = 33
      IF (IOPT .EQ. 0) MPOP = NPOP
      IF (IP .EQ. 0) MSAMP = 0
      ILIM = MIN0(MPOP,NPOP-IP)
      ICP = IP
      DO 5 I=1,ILIM
         IF (MSAMP .EQ. NSAMP) GO TO 10
         CALL GGUBS (DSEED,1,R)
         TEMP = FLOAT(NSAMP-MSAMP)/FLOAT(NPOP-IP)
         IP = IP + 1
         IF (R(1) .GT. TEMP) GO TO 5
         MSAMP = MSAMP + 1
         IX(MSAMP) = IP
         IF (IOPT .EQ. 0) GO TO 5
         SAMP(MSAMP) = POP(I)
    5 CONTINUE
      GO TO 9000
   10 IP = ICP + ILIM
 9000 CONTINUE
      IF (IER .GT. 0) CALL UERTST (IER,'GGSRS ')
 9005 RETURN
      END
