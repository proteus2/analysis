C   IMSL ROUTINE NAME   - NMTIE
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - TIE STATISTICS, GIVEN A SAMPLE OF OBSERVATIONS
C
C   USAGE               - CALL NMTIE (X,M,EPS,TIES)
C
C   ARGUMENTS    X      - INPUT VECTOR OF LENGTH M CONTAINING THE
C                           SAMPLE OF OBSERVATIONS. X MUST BE
C                           ORDERED MONOTONICALLY INCREASING.
C                M      - INPUT NUMBER OF OBSERVATIONS IN THE X VECTOR
C                EPS    - INPUT EPSILON USED TO DETERMINE IF TIES
C                           EXIST.  IF TWO OR MORE CONSECUTIVE ELEMENTS
C                           OF X ARE WITHIN EPS OF EACH OTHER, A TIE
C                           IS COUNTED.  IF THE ELEMENTS OF X WERE
C                           EQUALLY SPACED AND WITHIN EPS OF EACH
C                           OTHER, ALL ELEMENTS WOULD BE CONSIDERED
C                           TIED.
C                TIES   - OUTPUT VECTOR OF LENGTH 4.
C                           IN THE DESCRIPTION BELOW, T REFERS TO
C                           THE NUMBER OF OBSERVATIONS TIED FOR A
C                           GIVEN RANK, AND THE SUM IS OVER ALL RANKS.
C                         TIES(1) CONTAINS THE (SUM OF T*(T-1))/2
C                         TIES(2) CONTAINS THE (SUM OF T*(T-1)*(T+1))/12
C                         TIES(3) CONTAINS THE SUM OF T*(T-1)*(2*T+5)
C                         TIES(4) CONTAINS THE SUM OF T*(T-1)*(T-2)
C
C   REQD. IMSL ROUTINES - NONE REQUIRED
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS      ROUTINES IN CHAPTER V OF THE IMSL LIBRARY MIGHT
C                BE USED TO ORDER X MONOTONICALLY INCREASING.
C
C----------------------------------------------------------
------------
C
      SUBROUTINE NMTIE  (X,M,EPS,TIES)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            M
      REAL               X(1),EPS,TIES(4)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,NT
      REAL               PZET,ST,T
      DATA               PZET/.8333333E-1/
C                                  FIRST EXECUTABLE STATEMENT
      DO 5 I=1,4
         TIES(I) = 0.0
    5 CONTINUE
      NT = 0
      ST = 0.0
C                                  ACCUMULATE TIE STATISTICS
      DO 25 I=2,M
         IF ((X(I)-X(I-1)) .LE. EPS) GO TO 10
         NT = 0
         GO TO 15
   10    ST = ST + 1.0
         NT = NT + 1
         IF (I .EQ. M) GO TO 20
         GO TO 25
   15    IF (ST .EQ. 0.0) GO TO 25
   20    T = ST * (ST+1.0)
         TIES(1) = T + TIES(1)
         TIES(2) = T*(ST+2.0) + TIES(2)
         TIES(3) = T*(ST+ST+7.0) + TIES(3)
         TIES(4) = T*(ST-1.0) + TIES(4)
         ST = 0.0
   25 CONTINUE
      TIES(1) = 0.5 * TIES(1)
      TIES(2) = PZET * TIES(2)
      RETURN
      END
