C   IMSL ROUTINE NAME   - GGBN
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1980
C
C   PURPOSE             - BINOMIAL RANDOM DEVIATE GENERATOR
C
C   USAGE               - CALL GGBN (DSEED,NR,NIND,P,IR)
C
C   ARGUMENTS    DSEED  - INPUT/OUTPUT DOUBLE PRECISION VARIABLE
C                           ASSIGNED AN INTEGER VALUE IN THE
C                           EXCLUSIVE RANGE (1.D0, 2147483647.D0).
C                           DSEED IS REPLACED BY A NEW VALUE TO BE
C                           USED IN A SUBSEQUENT CALL.
C                NR     - INPUT NUMBER OF BINOMIAL RANDOM DEVIATES
C                           TO BE GENERATED.
C                NIND   - INPUT NON-NEGATIVE INTEGER BINOMIAL PARA-
C                           METER SPECIFYING THE NUMBER OF TRIALS.
C                P      - INPUT BINOMIAL PARAMETER (PROBABILITY OF
C                           SUCCESS ON ANY TRIAL) IN THE INCLUSIVE
C                           (0,1) RANGE.
C                IR     - OUTPUT VECTOR OF NR RANDOM DEVIATES.
C
C   REQD. IMSL ROUTINES - GGBTR,GGUBFS,GGUBS
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE GGBN (DSEED,NR,NIND,P,IR)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            NR,NIND,IR(1)
      REAL               P
      DOUBLE PRECISION   DSEED
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            NS,IRI
      REAL               A,AA(1),PS,PS2,PS3,PS4,U
      LOGICAL            SMALL
      DATA               ONE/1.0/
C
C                                  FIRST EXECUTABLE STATEMENT
      IF(NIND .LT. 35) GO TO 40
      DO 35 I=1,NR
         IRI = 0
         NS = NIND
         PS = P
    5    IF (NS.LT.15) GO TO 20
C                                  IS N ODD
         IF (MOD(NS,2).NE.0) GO TO 10
         NS = NS-1
         IF (GGUBFS(DSEED).LE.PS) IRI = IRI+1
C                                  GET THE MEDIAN OF NS UNIFORM DEVIATES
C                                    DEVIATES
   10    NS = NS/2
C                                  GET BETA DEVIATE
         Q = FLOAT(NS+1)
         CALL GGBTR (DSEED,Q,Q,1,AA)
         A = AA(1)
         IF (A.LE.PS) GO TO 15
         PS = PS/A
         GO TO 5
C                                  UPDATE COUNT OF POINTS LESS THAN P.
   15    IRI = IRI+NS+1
         PS = (PS-A)/(1.0-A)
         GO TO 5
C                                  COUNT UNIFORM DEVIATES LESS THAN P.
   20    IF (NS.LE.0) GO TO 30
         DO 25 J=1,NS
            IF (GGUBFS(DSEED).LE.PS) IRI = IRI+1
   25    CONTINUE
   30    IR(I) = IRI
   35 CONTINUE
      GO TO 9000
C                                  ALTERNATE LOOP FOR NIND LT 35.
C                                    COUNT UNIFORM DEVIATES.  REUSE
C                                    UNIFORMS IF POSSIBLE.
   40 PS = P
      SMALL = .FALSE.
      IF(PS .LT. .5) SMALL = .TRUE.
      IF(SMALL) PS = ONE - PS
      PS2 = PS*PS
      PS3 = PS*PS2
      PS4 = PS*PS3
      DO 60   I=1,NR
         IRI = 0
         NS = NIND
   45    IF(NS .LE. 0) GO TO 55
         U = GGUBFS(DSEED)
         IF(U .GT. PS) GO TO 50
         IRI = IRI + 1
         NS = NS -1
         IF(NS .LE. 0) GO TO 55
         IF(U .GT. PS2) GO TO 50
         IRI = IRI + 1
         NS = NS -1
         IF(NS .LE. 0) GO TO 55
         IF(U .GT. PS3) GO TO 50
         IRI = IRI + 1
         NS = NS -1
         IF(NS .LE. 0) GO TO 55
         IF(U .LE. PS4) IRI = IRI + 1
   50    NS = NS -1
         GO TO 45
   55    IR(I) = IRI
         IF(SMALL) IR(I) = NIND - IRI
   60 CONTINUE
 9000 CONTINUE
      RETURN
      END
