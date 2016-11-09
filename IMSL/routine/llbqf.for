C   IMSL ROUTINE NAME   - LLBQF
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - NOVEMBER 1, 1984
C
C   PURPOSE             - SOLUTION OF LINEAR LEAST SQUARES PROBLEM -
C                           HIGH ACCURACY SOLUTION
C
C   USAGE               - CALL LLBQF (A,IA,M,N,B,IB,NB,IND,C,X,IX,IWK,
C                           WK,IER)
C
C   ARGUMENTS    A      - M BY N COEFFICIENT MATRIX. (INPUT)
C                IA     - ROW DIMENSION OF MATRIX A EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN
C                           THE CALLING PROGRAM. (INPUT)
C                M      - NUMBER OF ROWS IN A AND B. (INPUT)
C                N      - NUMBER OF COLUMNS IN A AND ROWS IN X. (INPUT)
C                B      - M BY NB MATRIX OF RIGHT-HAND-SIDES. (INPUT)
C                IB     - ROW DIMENSION OF MATRIX B EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN
C                           THE CALLING PROGRAM. (INPUT)
C                NB     - NUMBER OF RIGHT-HAND-SIDES (COLUMNS IN B).
C                           (INPUT)
C                IND    - OPTION SELECTION INDICATOR. (INPUT)
C                         IF IND=0, ALL DEFAULT OPTIONS ARE SELECTED
C                           AND THERE IS NO NEED TO INITIALIZE ELEMENTS
C                           OF THE COMMUNICATION VECTOR C. IN THIS CASE,
C                           LLBQF SOLVES THE LEAST SQUARES PROBLEM GIVEN
C                           BY A AND B AND RETURNS THE SOLUTION IN X.
C                         IF IND=1, CERTAIN OPTIONS MAY BE SELECTED BY
C                           INITIALIZING C(I),I=1,...,3 TO ZERO AND
C                           SPECIFYING NON ZERO VALUES WHERE OPTIONS
C                           ARE DESIRED.
C                C      - COMMUNICATION VECTOR OF LENGTH 4. (INPUT IF
C                           IND.EQ.1, AND OUTPUT)
C                         C(1) - EQUALITY CONSTRAINTS.
C                           IF C(1)=M1 IS GREATER THAN ZERO AND LESS
C                           THAN OR EQUAL TO N, THE FIRST M1 EQUATIONS
C                           OF AX=B ARE INTERPRETED AS EQUALITY
C                           CONSTRAINTS WHICH MUST BE SATISFIED EXACTLY.
C                           THE FIRST M1 ROWS OF A MUST BE LINEARLY
C                           INDEPENDENT IN THIS CASE. FOR ALL OTHER
C                           VALUES OF C(1), INCLUDING C(1)=0, THERE ARE
C                           NO EQUALITY CONSTRAINTS.
C                         C(2) - TOLERANCE SPECIFICATION.
C                           COLUMN PIVOTING IS PERFORMED TO INTRODUCE
C                           COLUMNS OF A, ONE AT A TIME, INTO THE BASIS.
C                           AT EACH STEP, THE COLUMN THAT PRODUCES THE
C                           LARGEST REDUCTION OF THE RESIDUAL SUM OF
C                           SQUARES IS SELECTED. THE PROCESS IS
C                           TERMINATED IF, AFTER N1 STEPS, THE REMAINING
C                           N-N1 COLUMNS OF A ARE NEARLY LINEAR
C                           COMBINATIONS OF THE N1 BASIS COLUMNS. LET SS
C                           BE THE SUM OF SQUARES OF ELEMENTS OF ONE OF
C                           THE REMAINING COLUMNS. LET RS BE THE
C                           RESIDUAL SUM OF SQUARES AFTER APPROXIMATING
C                           THAT REMAINING COLUMN OF A BY THE N1
C                           LINEARLY INDEPENDENT COLUMNS. THEN THE
C                           CONDITION OF BEING NEARLY A LINEAR
C                           COMBINATION IS SATISFIED IF
C                                  SQRT(RS).LE.C(2)*SQRT(SS).
C                           WHEN THIS CONDITION HOLDS FOR ALL REMAINING
C                           COLUMNS, PIVOTING IS TERMINATED.
C                           THE VALUE OF N1 DETERMINED BY THE PROCESS
C                           IS RETURNED IN C(4).
C                         C(3) - BEST LEAST SQUARES SOLUTION OPTION.
C                           IF C(3)=1 THE BEST LEAST SQUARES SOLUTION
C                           IS COMPUTED. THAT IS, ALL COLUMNS OF A ARE
C                           USED TO COMPUTE THE LEAST SQUARES SOLUTION
C                           WITH MINIMUM NORM.
C                           FOR ALL OTHER VALUES OF C(3), INCLUDING
C                           C(3)=0, A BASIC LEAST SQUARES SOLUTION IS
C                           COMPUTED USING ONLY THE LINEARLY
C                           INDEPENDENT COLUMNS OF A AS DESCRIBED ABOVE
C                           (TOLERANCE SPECIFICATION). IF THE J-TH
C                           COLUMN OF A IS OMITTED FROM THE BASIS, X(J)
C                           IS SET TO ZERO.
C                         C(4) - NUMERICAL RANK OF A AS DETERMINED BY
C                           LLBQF. IF N1=C(4) IS LESS THAN N, THE BASIS
C                           COLUMNS ARE GIVEN BY THE FIRST N1 ELEMENTS
C                           OF IWK. IWK(I)=J MEANS THAT COLUMN J OF A
C                           WAS SELECTED FOR INCLUSION IN THE BASIS ON
C                           THE I-TH STEP.
C                X      - N BY NB SOLUTION MATRIX. (OUTPUT)
C                IX     - ROW DIMENSION OF MATRIX X EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN
C                           THE CALLING PROGRAM. (INPUT)
C                IWK    - WORK VECTOR OF LENGTH N.
C                WK     - WORK VECTOR OF LENGTH (M+N)*(N+3)+N. THE FINAL
C                           RESIDUAL VECTOR IS LOCATED IN WK(I), WHERE
C                           I = 1,...,M. (OUTPUT)
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                            IER=129 INDICATES INCORRECT ARGUMENTS N
C                              OR M. EITHER N.LE.0 OR M.LE.0.
C                            IER=130 INDICATES UNABLE TO COMPUTE A
C                              SOLUTION BECAUSE THE RANK OF A IS LESS
C                              THE NUMBER OF EQUALITY CONSTRAINTS.
C                            IER=131 INDICATES ITERATIVE REFINEMENT
C                              FAILED TO CONVERGE. X CONTAINS THE
C                              LAST APPROXIMATE SOLUTION COMPUTED, BUT
C                              THIS MAY BE INACCURATE.
C
C   REQD. IMSL ROUTINES - SINGLE/LLBQG,LLBQH,LLBQI,UERTST,UGETIO,
C                           VBLA=SDSDOT,VBLA=SDOT
C                       - DOUBLE/LLBQG,LLBQH,LLBQI,UERTST,UGETIO,
C                           VBLA=DDOT,VXADD,VXMUL,VXSTO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE LLBQF  (A,IA,M,N,B,IB,NB,IND,C,X,IX,IWK,WK,IER)
C
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            IA,M,N,IB,NB,IND,IX,IER,IWK(1)
      REAL               A(IA,N),B(IB,NB),C(4),X(IX,NB),WK(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            ID,IFAIL,IFF,IQRIJ,IQR,IRES,IV,IY,I,J,K,MPN,
     1                   NP1,M1,N1
      LOGICAL            BASIC
      REAL               TOL
C                                  FIRST EXECUTABLE STATEMENT
      IER = 129
      IF ((N.LE.0).OR.(M.LE.0)) GO TO 9000
C                                  SETUP
      IF (IND.EQ.1) GO TO 10
      DO 5 I=1,4
    5 C(I) = 0.0
   10 CONTINUE
      M1 = C(1)
      IF (M1.LT.0.OR.M1.GT.N) M1 = 0
      TOL = ABS(C(2))
      BASIC = C(3).NE.1.0
      MPN = M+N
      NP1 = N+1
      IRES = 1
      ID = IRES+M
      IQR = ID+N
      IFF = IQR+MPN*NP1
      IY = IFF+MPN
      IQRIJ = IQR
      DO 25 J=1,N
         DO 15 K=1,NB
   15    X(J,K) = 0.0
         DO 20 I=1,M
            WK(IQRIJ) = A(I,J)
            IQRIJ = IQRIJ+1
   20    CONTINUE
         IQRIJ = IQRIJ+N
   25 CONTINUE
      DO 30 I=1,M
         WK(IQRIJ) = B(I,1)
         IQRIJ = IQRIJ+1
   30 CONTINUE
C                                  QR DECOMPOSITION OF A
      CALL LLBQG (WK(IQR),MPN,NP1,M,N,M1,N1,BASIC,TOL,IWK,WK(ID))
      C(4) = N1
      IF (N1.EQ.0) GO TO 9005
      IER = 130
      IF (N1.LT.M1) GO TO 9000
      IER = 0
      IFAIL = 0
C                                  SOLVE NB RIGHT-HAND-SIDES
      DO 35 IV=1,NB
   35 CALL LLBQH (A,IA,M,N,M1,N1,B(1,IV),WK(IQR),MPN,BASIC,X(1,IV),IWK,W
     1K(IRES),WK(ID),WK(IFF),WK(IY),IFAIL)
      IF (IFAIL.EQ.0) GO TO 9005
      IER = 131
 9000 CONTINUE
      CALL UERTST (IER,6HLLBQF )
 9005 RETURN
      END
