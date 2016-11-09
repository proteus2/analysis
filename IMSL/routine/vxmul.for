C   IMSL ROUTINE NAME   - VXMUL
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - EXTENDED PRECISION MULTIPLY
C
C   USAGE               - CALL VXMUL (A,B,ACC)
C
C   ARGUMENTS    A      - INPUT DOUBLE PRECISION NUMBER
C                B      - INPUT DOUBLE PRECISION NUMBER
C                ACC    - ACCUMULATOR. (INPUT AND OUTPUT)
C                           ACC IS A DOUBLE PRECISION VECTOR OF LENGTH
C                           2.  ON OUTPUT, ACC CONTAINS THE SUM OF
C                           INPUT ACC AND A*B.
C
C   REQD. IMSL ROUTINES - VXADD
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS      VXMUL ADDS THE PRODUCT A*B TO THE EXTENDED PRECISION
C                ACCUMULATOR, ACC. THE SUBROUTINE ASSUMES THAT AN
C                EXTENDED PRECISION NUMBER IS ALREADY IN THE
C                ACCUMULATOR.  THEREFORE, BEFORE THE FIRST CALL TO
C                VXMUL, ACC(1) AND ACC(2) MUST BE SET TO ZERO.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE VXMUL (A,B,ACC)
C                                  SPECIFICATIONS FOR ARGUMENTS
      DOUBLE PRECISION   A,B,ACC(2)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      DOUBLE PRECISION   X,HA,TA,HB,TB
      INTEGER            IX(2),I
      LOGICAL*1          LX(8),LI(4)
      EQUIVALENCE        (X,LX(1),IX(1)),(I,LI(1))
      DATA               I/0/
C                                  SPLIT A = HA+TA
C                                        B = HB+TB
C                                  FIRST EXECUTABLE STATEMENT
      X = A
      LI(4) = LX(5)
      IX(2) = 0
      I = (I/16)*16
      LX(5) = LI(4)
      HA=X
      TA=A-HA
      X = B
      LI(4) = LX(5)
      IX(2) = 0
      I = (I/16)*16
      LX(5) = LI(4)
      HB = X
      TB = B-HB
C                                  COMPUTE HA*HB,HA*TB,TA*HB, AND TA*TB
C                                    AND CALL VXADD TO ACCUMULATE THE
C                                    SUM
      X = TA*TB
      CALL VXADD(X,ACC)
      X = HA*TB
      CALL VXADD(X,ACC)
      X = TA*HB
      CALL VXADD(X,ACC)
      X = HA*HB
      CALL VXADD(X,ACC)
      RETURN
      END
