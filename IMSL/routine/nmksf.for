C   IMSL ROUTINE NAME   - NMKSF
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1982
C
C   PURPOSE             - FREQUENCY DISTRIBUTION OF K AND THE
C                           PROBABILITY OF EQUALLING OR EXCEEDING K,
C                           WHERE K, THE TOTAL SCORE FROM THE KENDALL
C                           RANK CORRELATION COEFFICIENT CALCULATIONS,
C                           AND N, THE SAMPLE SIZE, ARE GIVEN
C
C   USAGE               - CALL NMKSF (K,N,Z,ZW,P)
C
C   ARGUMENTS    K      - INPUT SCORE FOR WHICH THE PROBABILITY (OF
C                           BEING NOT LESS THAN K) IS CALCULATED.
C                           K MUST BE IN THE RANGE
C                           (-N*(N-1)/2, N*(N-1)/2), INCLUSIVELY.
C                N      - INPUT SIZE OF THE SAMPLE.
C                           N MUST BE GREATER THAN 1 AND LESS THAN 34,
C                           OR LESS THAN 56, ON SOME COMPUTERS.  (SEE
C                           PROGRAMMING NOTES.)
C                Z      - OUTPUT VECTOR OF LENGTH N*(N-1)/2+1 CONTAINING
C                           THE FREQUENCY DISTRIBUTION OF POSSIBLE
C                           VALUES OF K. K WILL RANGE BETWEEN
C                           -N*(N-1)/2 TO N*(N-1)/2, INCLUSIVELY, IN
C                           INCREMENTS OF 2, WITH FREQUENCY Z(I) FOR
C                           A K VALUE OF -N*(N-1)/2+2*(I-1),
C                           FOR I=1,...,N*(N-1)/2+1.
C                ZW     - WORKING AREA VECTOR OF LENGTH (N-1)*(N-2)/2+1.
C                           ON OUTPUT, ZW(1) = I, WHERE I IS THE
C                           ELEMENT OF Z CORRESPONDING TO THE SCORE K.
C                P      - OUTPUT PROBABILITY OF EQUALLING OR EXCEEDING
C                           K IF THE SAMPLES ON WHICH K IS BASED ARE
C                           UNCORRELATED.
C
C   REQD. IMSL ROUTINES - NONE REQUIRED
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE NMKSF  (K,N,Z,ZW,P)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            K,N
      REAL               Z(1),ZW(1),P
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,IC,I1,I2,JC,JST,KC,L,M
      REAL               SUM,X,Y
C                                  FIRST EXECUTABLE STATEMENT
      IC = 2
      KC = 1
      Z(1) = 1
      Z(2) = 1
C                                  DO ROW SUMMATION CALCULATION
    5 IF(IC.EQ.N) GO TO 25
      IC = IC+1
      KC = KC+IC
      JC = (KC+1)/2
      DO 15  I=1,JC
         SUM = 0.0
C                                  MOVE A MEMBER TO WORKAREA
         ZW(I) = Z(I)
         JST = I-IC+1
         IF(JST.LT.1) JST=1
C                                  SUM THE ROW TO NEXT Z
         DO 10 I1=JST,I
   10    SUM = SUM+ZW(I1)
   15 Z(I) = SUM
      IF(IC.GT.3) KC=KC-1
      I2 = KC
      I1 = KC-JC
      DO 20 I=1,I1
         Z(I2) = Z(I)
C                                  SET LAST HALF OF ROW
   20 I2 = I2-1
      GO TO 5
   25 IC = 1
      M = -N*(N-1)*.5
      L = -M+1
C                                  DO THE PROBABILITY CALCULATION
   30 IF(K.LE.M) GO TO 35
         M = M+2
         IC = IC+1
      GO TO 30
   35 X = 0.0
      Y = 0.0
      DO 40 I=1,L
   40 Y = Y + Z(I)
      DO 45 I=IC,L
   45 X = X+Z(I)
      P = X/Y
      ZW(1) = IC
      RETURN
      END
