C   IMSL ROUTINE NAME   - ZREAL1
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - THE REAL ZEROS OF A REAL FUNCTION - TO BE
C                           USED WHEN INITIAL GUESSES ARE POOR
C
C   USAGE               - CALL ZREAL1 (F,EPS,EPS2,ETA,NSIG,N,X,ITMAX,
C                           IER)
C
C   ARGUMENTS    F      - A SINGLE-ARGUMENT REAL FUNCTION SUBPROGRAM
C                           SUPPLIED BY THE USER. (INPUT)
C                           F MUST BE DECLARED EXTERNAL IN THE CALLING
C                           PROGRAM. F DEFINES THE FUNCTION FOR WHICH
C                           THE ROOTS ARE TO BE FOUND.
C                EPS    - CONVERGENCE CRITERION. (INPUT)
C                           A ROOT, X(I), IS ACCEPTED IF
C                           ABS(F(X(I)) .LE. EPS.
C                EPS2   - SPREAD CRITERIA FOR MULTIPLE ROOTS. (INPUT)
C                ETA        IF THE ROOT X(I) HAS BEEN COMPUTED AND IT IS
C                           FOUND THAT
C                           ABS(X(I)-X(J)) .LT. EPS2
C                           WHERE X(J) IS A PREVIOUSLY COMPUTED ROOT,
C                           THEN THE COMPUTATION IS RESTARTED WITH A
C                           GUESS EQUAL TO X(I) + ETA.
C                NSIG   - CONVERGENCE CRITERION. (INPUT)
C                           A ROOT IS ACCEPTED IF TWO SUCCESSIVE
C                           APPROXIMATIONS TO A GIVEN ROOT AGREE
C                           IN THE FIRST NSIG DIGITS.
C                         NOTE THAT IF EITHER CONVERGENCE CRITERION
C                         IS SATISFIED, THE ROOT IS ACCEPTED.
C                N      - THE NUMBER OF ROOTS TO BE FOUND. (INPUT)
C                X      - VECTOR OF LENGTH N. (INPUT/OUTPUT)
C                           ON INPUT, X CONTAINS THE INITIAL GUESSES
C                           FOR THE ROOTS.
C                           ON OUTPUT, X CONTAINS THE COMPUTED ROOTS.
C                ITMAX  - ITERATION INDICATOR. (INPUT/OUTPUT)
C                           ON INPUT, ITMAX IS THE MAXIMUM NUMBER OF
C                           ITERATIONS TO BE TAKEN PER ROOT.
C                           ON OUTPUT, ITMAX IS THE NUMBER OF ITERATIONS
C                           USED IN FINDING THE LAST ROOT.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         WARNING ERROR
C                           IER=33 INDICATES THAT CONVERGENCE WAS NOT
C                             OBTAINED WITHIN ITMAX ITERATIONS
C                             FOR AT LEAST ONE INITIAL GUESS,
C                             X(I), I=1,...,N. X(I) IS SET TO 111111.
C                             NOTE THAT THE ROUTINE IS DESIGNED SO THAT
C                             A MULTIPLE ROOT WILL NOT APPEAR IN THE
C                             OUTPUT VECTOR X.
C
C   REQD. IMSL ROUTINES - UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS  1.  ZREAL1 ASSUMES THAT THERE EXIST N DISTINCT REAL ROOTS
C                FOR THE FUNCTION F AND THAT THEY CAN BE REACHED FROM
C                THE INITIAL GUESSES SUPPLIED. THE ROUTINE IS DESIGNED
C                SO THAT CONVERGENCE TO ANY SINGLE ROOT CANNOT BE
C                OBTAINED FROM TWO DIFFERENT INITIAL GUESSES.
C            2.  SCALING THE X VECTOR IN THE FUNCTION F MAY BE REQUIRED
C                IF ANY OF THE ROOTS ARE KNOWN TO BE LESS THAN ONE.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE ZREAL1 (F,EPS,EPS2,ETA,NSIG,N,X,ITMAX,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            NSIG,N,ITMAX,IER
      REAL               F,EPS,EPS2,ETA,X(N)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            L,I,JK
      REAL               P,P1,P2,X0,X1,X2,RT,FRT,FPRT,D,DD,DI,H,BI,DEN,
     1                   DN,DM,TEM,DIGT,TEN,ONE,ZERO,P9,P11,HALF,PP1,F4
      DATA               TEN,ONE,ZERO,P9/10.0,1.0,0.0,.9/
      DATA               P11,HALF,PP1,F4/1.1,.5,.1,4.0/
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
      DIGT = TEN**(-NSIG)
      P = -ONE
      P1 = ONE
      P2 = ZERO
      H = ZERO
      DO 95 L = 1,N
         JK = 0
         IF (X(L) .EQ. ZERO) GO TO 5
         P = P9*X(L)
         P1 = P11*X(L)
         P2 = X(L)
    5    RT = P
         GO TO 65
   10    IF (JK .NE. 1) GO TO 15
         RT = P1
         X0 = FPRT
         GO TO 65
   15    IF (JK .NE. 2) GO TO 20
         RT = P2
         X1 = FPRT
         GO TO 65
   20    IF (JK .NE. 3) GO TO 55
         X2 = FPRT
         D = -HALF
         IF (X(L) .EQ. ZERO) GO TO 25
         H =-PP1*X(L)
         GO TO 30
   25    H = -ONE
   30    DD = ONE+D
         BI = X0*D**2-X1*DD**2+X2*(DD+D)
         DEN = BI**2 -F4*X2*D*DD*(X0*D-(X1*DD)+X2)
         IF (DEN .LE. ZERO) GO TO 35
         DEN = SQRT(DEN)
         GO TO 40
   35    DEN = ZERO
   40    DN = BI + DEN
         DM = BI - DEN
         IF (ABS(DN) .LE. ABS(DM)) GO TO 45
         DEN = DN
         GO TO 50
   45    DEN = DM
   50    IF (DEN .EQ. ZERO) DEN = ONE
         DI=-DD*(X2+X2)/DEN
         H = DI * H
         RT = RT + H
C                                  TEST FOR CONVERGENCE
         IF (ABS(H) .LT. ABS(RT)*DIGT) GO TO 90
         GO TO 65
   55    IF (ABS(FPRT) .GE. ABS(X2*10.0))  GO TO 60
         X0 = X1
         X1 = X2
         X2 = FPRT
         D = DI
         GO TO 30
   60    DI = DI * HALF
         H = H * HALF
         RT = RT - H
   65    JK = JK + 1
         IF (JK .LT. ITMAX)  GO TO 75
C                                  WARNING  ERROR ITERATIONS = MAXIMUM
         IER=33
         X(L)=111111.
         GO TO 95
   75    FRT = F(RT)
         FPRT = FRT
         IF (L .LT. 2) GO TO 81
         DO 80 I = 2,L
            TEM = RT - X(I-1)
            IF (ABS(TEM) .LT. EPS2)  GO TO 85
            FPRT = FPRT/TEM
   80    CONTINUE
C                                  TEST FOR CONVERGENCE
   81    IF ((ABS(FRT) .LT. EPS) .AND. (ABS(FPRT) .LT. EPS))  GO TO 90
         GO TO 10
   85    RT = RT + ETA
         JK = JK - 1
         GO TO 65
   90    X(L) = RT
   95 CONTINUE
      ITMAX = JK
      IF(IER.EQ.0) GO TO 9005
 9000 CONTINUE
      CALL UERTST (IER,6HZREAL1)
 9005 RETURN
      END
