C   IMSL ROUTINE NAME   - GGCHS
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1980
C
C   PURPOSE             - CHI-SQUARED RANDOM DEVIATE GENERATOR
C
C   USAGE               - CALL GGCHS (DSEED,N,R,CHI2)
C
C   ARGUMENTS    DSEED  - INPUT/OUTPUT DOUBLE PRECISION VARIABLE
C                           ASSIGNED AN INTEGER VALUE IN THE
C                           EXCLUSIVE RANGE (1.D0, 2147483647.D0).
C                           DSEED IS REPLACED BY A NEW VALUE TO BE
C                           USED IN A SUBSEQUENT CALL.
C                N      - INPUT DEGREES OF FREEDOM OF THE CHI-SQUARED
C                           DISTRIBUTION.
C                R      - WORK VECTOR OF LENGTH THE LARGEST INTEGER
C                           NOT GREATER THAN N/2. R CONTAINS UNIFORM
C                           RANDOM DEVIATES.
C                CHI2   - OUTPUT CHI-SQUARED DEVIATE.
C
C   REQD. IMSL ROUTINES - GGNQF,GGUBS,MDNRIS,MERFI,UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE GGCHS  (DSEED,N,R,CHI2)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N
      REAL               R(1),CHI2
      DOUBLE PRECISION   DSEED
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            M,L,I
      REAL               D1,D2,B,SE,SMALL,ALRAD,TWOM,S
C                        ALRAD = LN(RADIX)
C                        SE = RADIX**S WHERE S IS THE SMALLEST NUMBER
C                        SUCH THAT RADIX**S IS GREATER THAN OR EQUAL TO
C                        2**32
C
C                        SMALL = 4*16**-8
      DATA               SE/.4294967E10/
      DATA               SMALL/.9313226E-09/
      DATA               ALRAD/2.772589/
      DATA               S/8./
      DATA               TWOM/-2.0/
C                                  FIRST EXECUTABLE STATEMENT
      IF(N .GT. 1) GO TO 5
      CHI2 = 0.0
      GO TO 20
    5 M = N/2
      L = N-M
      CALL GGUBS(DSEED,M,R)
      D1 = 1.0
      D2 = 0.0
      DO 15 I=1,M
         D1 = D1*R(I)
         IF(D1 .GT. SMALL) GO TO 15
   10    D1 = D1*SE
         D2 = D2-S
         IF(D1 .LT. 1.0) GO TO 10
   15 CONTINUE
C                                  CHI2 = -2*LN(D1*RADIX**D2)
      CHI2 = TWOM*(ALOG(D1)+D2*ALRAD)
C                                  TEST FOR N EVEN
      IF(L .EQ. M) GO TO 25
   20 B = GGNQF(DSEED)
      CHI2 = CHI2+B*B
   25 RETURN
      END
