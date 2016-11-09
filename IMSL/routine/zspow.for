C   IMSL ROUTINE NAME   - ZSPOW
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JULY 1, 1983
C
C   PURPOSE             - SOLVE A SYSTEM OF NONLINEAR EQUATIONS
C
C   USAGE               - CALL ZSPOW (FCN,NSIG,N,ITMAX,PAR,X,FNORM,
C                           WK,IER)
C
C   ARGUMENTS    FCN    - THE NAME OF A USER-SUPPLIED SUBROUTINE WHICH
C                           EVALUATES THE SYSTEM OF EQUATIONS TO BE
C                           SOLVED. FCN MUST BE DECLARED EXTERNAL IN
C                           THE CALLING PROGRAM AND MUST HAVE THE
C                           FOLLOWING FORM,
C                             SUBROUTINE FCN(X,F,N,PAR)
C                             REAL X(N),F(N),PAR(1)
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
C                           IN THE COMPUTED ROOT. (INPUT)
C                N      - THE NUMBER OF EQUATIONS TO BE SOLVED AND
C                           THE NUMBER OF UNKNOWNS. (INPUT)
C                ITMAX  - THE MAXIMUM ALLOWABLE NUMBER OF ITERATIONS.
C                           (INPUT) THE MAXIMUM NUMBER OF CALLS TO FCN
C                           IS ITMAX*(N+1). SUGGESTED VALUE = 200.
C                PAR    - PAR CONTAINS A PARAMETER SET WHICH IS
C                           PASSED TO THE USER-SUPPLIED FUNCTION FCN.
C                           PAR MAY BE USED TO PASS ANY AUXILIARY
C                           PARAMETERS NECESSARY FOR COMPUTATION OF
C                           THE FUNCTION FCN. (INPUT)
C                X      - A VECTOR OF LENGTH N. (INPUT/OUTPUT) ON INPUT,
C                           X IS THE INITIAL APPROXIMATION TO THE ROOT.
C                           ON OUTPUT, X IS THE BEST APPROXIMATION TO
C                           THE ROOT FOUND BY ZSPOW.
C                FNORM  - ON OUTPUT, FNORM IS EQUAL TO
C                           F(1)**2+...F(N)**2 AT THE POINT X.
C                WK     - WORK VECTOR OF LENGTH N*(3*N+15)/2
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER = 129 INDICATES THAT THE NUMBER OF
C                             CALLS TO FCN HAS EXCEEDED ITMAX*(N+1).
C                             THE USER MAY TRY A NEW INITIAL GUESS.
C                           IER = 130 INDICATES THAT NSIG IS TOO
C                             LARGE.  NO FURTHER IMPROVEMENT IN THE
C                             APPROXIMATE SOLUTION X IS POSSIBLE.
C                             THE USER SHOULD DECREASE NSIG.
C                           IER = 131 INDICATES THAT THE ITERATION
C                             HAS NOT MADE GOOD PROGRESS.  THE USER
C                             MAY TRY A NEW INITIAL GUESS.
C
C   REQD. IMSL ROUTINES - SINGLE/UERTST,UGETIO,VBLA=SNRM2,ZSPWA,
C                           ZSPWB,ZSPWC,ZSPWD,ZSPWE,ZSPWF,ZSPWG
C                       - DOUBLE/UERTST,UGETIO,VBLA=DNRM2,ZSPWA,
C                           ZSPWB,ZSPWC,ZSPWD,ZSPWE,ZSPWF,ZSPWG
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE ZSPOW (FCN,NSIG,N,ITMAX,PAR,X,FNORM,WK,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            NSIG,N,ITMAX,IER
      REAL               PAR(1),X(N),FNORM,WK(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            INDEX2,INDEX,INFO,I,J,LR,MAXFEV,ML,MODE,MU,
     1                   NFEV,NPRINT
      REAL               EPSFCN,FACTOR,ONE,XTOL,ZERO
      EXTERNAL           FCN
      DATA               FACTOR,ONE,ZERO /1.0E2,1.0E0,0.0E0/
C                                  FIRST EXECUTABLE STATEMENT
      INFO = 0
C                                  CALL ZSPWA
      MAXFEV = ITMAX*(N + 1)
      XTOL = 0.1**NSIG
      ML = N - 1
      MU = N - 1
      EPSFCN = ZERO
      MODE = 2
      DO 5 J = 1, N
         WK(J) = ONE
    5 CONTINUE
      NPRINT = 0
      LR = (N*(N + 1))/2
      INDEX = 7*N + LR
      CALL ZSPWA(FCN,N,X,WK(6*N+1),XTOL,MAXFEV,ML,MU,EPSFCN,WK(1),
     * MODE,FACTOR,NPRINT,INFO,NFEV,WK(INDEX+1),N,WK(7*N+1),LR,
     * WK(N+1),WK(2*N+1),WK(3*N+1),WK(4*N+1),WK(5*N+1),PAR)
      IF (INFO .EQ. 5) INFO = 4
      FNORM = 0.0
      DO 10 I=1,N
         INDEX2 = 6*N+I
         FNORM = FNORM+WK(INDEX2)*WK(INDEX2)
   10 CONTINUE
      IER = 0
      IF (INFO .EQ. 2) IER = 129
      IF (INFO .EQ. 3) IER = 130
      IF (INFO .EQ. 4) IER = 131
      IF (IER .GT. 0) CALL UERTST(IER,6HZSPOW )
      RETURN
      END
