C   IMSL ROUTINE NAME   - ZSCNT
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1982
C
C   PURPOSE             - SOLVE A SYSTEM OF NONLINEAR EQUATIONS
C
C   USAGE               - CALL ZSCNT (FCN,NSIG,N,ITMAX,PAR,X,FNORM,
C                           WK,IER)
C
C   ARGUMENTS    FCN    - THE NAME OF A USER-SUPPLIED SUBROUTINE WHICH
C                           EVALUATES THE SYSTEM OF EQUATIONS TO BE
C                           SOLVED. FCN MUST BE DECLARED EXTERNAL IN
C                           THE CALLING PROGRAM AND MUST HAVE THE
C                           FOLLOWING FORM,
C                             SUBROUTINE FCN(X,F,N,PAR)
C                             DIMENSION X(N),F(N),PAR(1)
C                             F(1)=
C                              .
C                             F(N)=
C                             RETURN
C                             END
C                           GIVEN X(1)...X(N), FCN MUST EVALUATE THE
C                           FUNCTIONS F(1)...F(N) WHICH ARE TO BE MADE
C                           ZERO. X SHOULD NOT BE ALTERED BY FCN. THE
C                           PARAMETERS IN VECTOR PAR (SEE ARGUMENT
C                           PAR BELOW) MAY ALSO BE USED IN THE
C                           CALCULATION OF F(1)...F(N).
C                NSIG   - THE NUMBER OF DIGITS OF ACCURACY DESIRED
C                           IN THE COMPUTED ROOT (INPUT).
C                N      - THE NUMBER OF EQUATIONS TO BE SOLVED AND
C                           THE NUMBER OF UNKNOWNS (INPUT).
C                ITMAX  - THE MAXIMUM ALLOWABLE NUMBER OF ITERATIONS
C                           (INPUT).
C                PAR    - PAR CONTAINS A PARAMETER SET WHICH IS
C                           PASSED TO THE USER SUPPLIED FUNCTION FCN.
C                           PAR MAY BE USED TO PASS ANY AUXILIARY
C                           PARAMETERS NECESSARY FOR COMPUTATION OF
C                           THE FUNCTION FCN. (INPUT)
C                X      - A VECTOR OF LENGTH N. (INPUT/OUTPUT) ON INPUT,
C                           X IS THE INITIAL APPROXIMATION TO THE ROOT.
C                           ON OUTPUT, X IS THE BEST APPROXIMATION TO
C                           THE ROOT FOUND BY ZSCNT.
C                FNORM  - ON OUTPUT, FNORM IS EQUAL TO
C                           F(1)**2+...F(N)**2 AT THE POINT X.
C                WK     - WORK VECTOR OF LENGTH (N+1)*(3*N+8)
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER = 129 INDICATES THAT ZSCNT FAILED TO
C                             CONVERGE WITHIN ITMAX ITERATIONS. THE
C                             USER MAY INCREASE ITMAX OR TRY A NEW
C                             INITIAL GUESS.
C                           IER = 130 INDICATES THE ALGORITHM WAS
C                             UNABLE TO IMPROVE ON THE RETURNED VALUE
C                             OF X. THIS SITUATION ARISES WHEN THE
C                             SOLUTION CANNOT BE DETERMINED TO NSIG
C                             DIGITS DUE TO ERRORS IN THE FUNCTION
C                             VALUES. IT MAY ALSO INDICATE THAT THE
C                             ROUTINE IS TRAPPED IN THE AREA OF A
C                             LOCAL MINIMUM. THE USER MAY TRY A NEW
C                             INITIAL GUESS.
C
C   REQD. IMSL ROUTINES - SINGLE/GGUBFS,LEQT2F,LUDATN,LUELMN,LUREFN,
C                           UERSET,UERTST,UGETIO,ZSCNU
C                       - DOUBLE/GGUBFS,LEQT2F,LUDATN,LUELMN,LUREFN,
C                           UERSET,UERTST,UGETIO,VXADD,VXMUL,VXSTO,
C                           ZSCNU
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE ZSCNT  (FCN,NSIG,N,ITMAX,PAR,X,FNORM,WK,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            IER,ITMAX,N,NSIG
      REAL               FNORM,PAR(1),WK(1),X(N)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I1,I2,I3,I4,I5,LNEW,LOLD,N1
      EXTERNAL           FCN
C                                  FIRST EXECUTABLE STATEMENT
      N1=N+1
      I1 = N1*N1+1
      I2 = I1+N*N1
      I3 = I2+N1
      I4 = I3+N1
      I5 = I4+N1
      CALL UERSET(0,LOLD)
      CALL ZSCNU(X,N,FCN,NSIG,N1,WK(1),WK(I1),WK(I2),WK(I3),WK(I4)
     * ,WK(I5),ITMAX,PAR,IER)
      CALL FCN(X,WK,N,PAR)
      FNORM = 0.0
      DO 5 I = 1,N
         FNORM = FNORM + WK(I)*WK(I)
    5 CONTINUE
      CALL UERSET(LOLD,LNEW)
      IF (IER.EQ.0) GO TO 9005
 9000 CONTINUE
      CALL UERTST(IER,6HZSCNT )
 9005 RETURN
      END
