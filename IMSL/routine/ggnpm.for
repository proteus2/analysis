C   IMSL ROUTINE NAME   - GGNPM
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - NOVEMBER 1, 1984
C
C   PURPOSE             - NORMAL RANDOM DEVIATE GENERATOR VIA THE POLAR
C                           METHOD
C
C   USAGE               - CALL GGNPM (DSEED,NR,R)
C
C   ARGUMENTS    DSEED  - INPUT/OUTPUT DOUBLE PRECISION VARIABLE
C                           ASSIGNED AN INTEGER VALUE IN THE
C                           EXCLUSIVE RANGE (1.D0, 2147483647.D0).
C                           DSEED IS REPLACED BY A NEW VALUE TO BE
C                           USED IN A SUBSEQUENT CALL.
C                NR     - INPUT.  THE NUMBER OF NORMAL DEVIATES TO BE
C                           GENERATED.
C                R      - OUTPUT VECTOR OF LENGTH NR CONTAINING THE
C                           NORMAL DEVIATES.
C
C   REQD. IMSL ROUTINES - NONE REQUIRED
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   WARRANTY            - IMSL WARRANTS ONLY THAT IMSL TESTING HAS BEEN
C                           APPLIED TO THIS CODE.  NO OTHER WARRANTY,
C                           EXPRESSED OR IMPLIED, IS APPLICABLE.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE GGNPM (DSEED,NR,R)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            NR
      REAL               R(NR)
      DOUBLE PRECISION   DSEED
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,NN,M
      REAL               U,V,SUM,SLN
      DOUBLE PRECISION   D2P31M,D2PN31
C                                  D2P31M = (2**31)-1
C                                  D2PN31 = (2**31)
      DATA               D2P31M/2147483647.D0/
      DATA               D2PN31/2147483648.D0/
C                                  FIRST EXECUTABLE STATEMENT
      NN = NR
C                                  TEST FOR NR ODD
      M = MOD(NR,2)
      IF (M .NE. 0) NN = NN-1
      IF (NR .EQ. 1) GO TO 15
      DO 10 I = 1,NN,2
    5    DSEED = DMOD(16807.D0*DSEED,D2P31M)
         U = DSEED/D2PN31
         DSEED = DMOD(16807.D0*DSEED,D2P31M)
         V = DSEED/D2PN31
         U = U+U-1.0
         V = V+V-1.0
         SUM = U*U+V*V
         IF(SUM .GE. 1.0) GO TO 5
         SLN = ALOG(SUM)
         SLN = SQRT((-SLN-SLN)/SUM)
         R(I) = U*SLN
         R(I+1) = V*SLN
   10 CONTINUE
      IF (NR .EQ. NN) GO TO 20
   15    DSEED = DMOD(16807.D0*DSEED,D2P31M)
         U = DSEED/D2PN31
         DSEED = DMOD(16807.D0*DSEED,D2P31M)
         V = DSEED/D2PN31
         U = U+U-1.0
         V = V+V-1.0
         SUM = U*U+V*V
         IF (SUM .GE. 1.0) GO TO 15
         SLN = ALOG(SUM)
         SLN = SQRT((-SLN-SLN)/SUM)
         R(NR) = U*SLN
   20 CONTINUE
      RETURN
      END
