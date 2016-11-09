C   IMSL ROUTINE NAME   - GTPBC
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - COUNT OF THE NUMBER OF ZERO BITS IN A GIVEN
C                           SUBSET OF A REAL WORD
C
C   USAGE               - CALL GTPBC (I1,I2,R,KZERO)
C
C   ARGUMENTS    I1     - INPUT VALUE OF THE BIT POSITION WITHIN WORD R
C                           WHERE COUNT WILL BEGIN. (I1 MUST BE GREATER
C                           THAN 0 AND LESS THAN OR EQUAL TO RBITS WHERE
C                           RBITS REPRESENTS THE NUMBER OF BITS IN A
C                           REAL WORD. NO TEST IS MADE FOR I1 BEING
C                           WITHIN THIS RANGE.)
C                I2     - INPUT VALUE OF BIT POSITION WITHIN WORD R
C                           WHERE COUNT WILL STOP. (I2 MUST BE GREATER
C                           THAN OR EQUAL TO I1 AND LESS THAN OR EQUAL
C                           TO RBITS. NO TEST IS MADE FOR I2 BEING
C                           WITHIN THIS RANGE.)
C                R      - INPUT WORD UPON WHICH THE COUNT OF BITS
C                           IS TO BE MADE
C                KZERO  - OUTPUT COUNT OF NUMBER OF ZERO BITS IN R
C                           BETWEEN BITS I1 AND I2, INCLUSIVELY.
C
C   REQD. IMSL ROUTINES - NONE REQUIRED
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE GTPBC  (I1,I2,R,KZERO)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            I1,I2,KZERO
      REAL               R
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            IR,ISTART,K,IBITS,J1,J2,IDJ
      REAL               S
      EQUIVALENCE        (IR,S)
      DATA               IBITS/32/
C                                  FIRST EXECUTABLE STATEMENT
      S = R
      J1 = I1
      J2 = I2
      K = J2-J1+1
      IDJ = 2**(IBITS-J2)
      ISTART = IR/IDJ
      KZERO = 0
      DO 10 I=1,K
         IF(MOD(ISTART,2).EQ.1) GO TO 10
         KZERO = KZERO+1
   10 ISTART = ISTART/2
      RETURN
      END
