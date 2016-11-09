C   IMSL ROUTINE NAME   - VXADD
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - EXTENDED PRECISION ADD
C
C   USAGE               - CALL VXADD (A,ACC)
C
C   ARGUMENTS    A      - DOUBLE PRECISION NUMBER TO BE ADDED TO THE
C                           ACCUMULATOR. (INPUT)
C                ACC    - ACCUMULATOR. (INPUT AND OUTPUT)
C                           ACC IS A DOUBLE PRECISION VECTOR OF LENGTH
C                           2. ON OUTPUT, ACC CONTAINS THE SUM OF
C                           INPUT ACC AND A.
C
C   REQD. IMSL ROUTINES - NONE REQUIRED
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS      VXADD ADDS THE DOUBLE PRECISION NUMBER A TO THE
C                EXTENDED PRECISION ACCUMULATOR, ACC. THE SUBROUTINE
C                ASSUMES THAT AN EXTENDED PRECISION NUMBER IS ALREADY IN
C                THE ACCUMULATOR. THEREFORE, BEFORE THE FIRST CALL TO
C                VXADD, ACC(1) AND ACC(2) MUST BE SET TO ZERO.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE VXADD(A,ACC)
C
C                                  SPECIFICATIONS FOR ARGUMENTS
      DOUBLE PRECISION   A,ACC(2)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      DOUBLE PRECISION   X,Y,Z,ZZ
C                                  FIRST EXECUTABLE STATEMENT
      X = ACC(1)
      Y = A
      IF (DABS(ACC(1)).GE.DABS(A)) GO TO 1
      X = A
      Y = ACC(1)
C                                  COMPUTE Z+ZZ = ACC(1)+A EXACTLY
    1 Z = X+Y
      ZZ = (X-Z)+Y
C                                  COMPUTE ZZ+ACC(2) USING DOUBLE
C                                    PRECISION ARITHMETIC
      ZZ = ZZ+ACC(2)
C                                  COMPUTE ACC(1)+ACC(2) = Z+ZZ EXACTLY
      ACC(1) = Z+ZZ
      ACC(2) = (Z-ACC(1))+ZZ
      RETURN
      END
