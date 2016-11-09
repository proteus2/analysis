C   IMSL ROUTINE NAME   - ZREAL2
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - THE REAL ZEROS OF A REAL FUNCTION - TO BE
C                           USED WHEN INITIAL GUESSES ARE GOOD
C
C   USAGE               - CALL ZREAL2 (F,EPS,EPS2,ETA,NSIG,N,X,ITMAX,
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
C                           IER=33 INDICATES THAT FOR ONE ROOT,
C                             CONVERGENCE WAS NOT OBTAINED WITHIN ITMAX
C                             ITERATIONS. THAT ROOT IS SET TO 111111.
C                           IER=34 INDICATES THAT FOR ONE ROOT,
C                             THE DERIVATIVE OF THE FUNCTION AT THAT
C                             ROOT WAS TOO SMALL. THAT ROOT IS SET TO
C                             222222.
C                           IER=35 INDICATES THAT THE ERROR CONDITIONS
C                             DESCRIBED FOR IER=33 AND IER=34 ABOVE,
C                             OCCURRED MORE THAN ONCE. THE ROOTS FOR
C                             WHICH THE ERROR OCCURRED ARE SET TO
C                             111111. OR 222222., DEPENDING ON THE TYPE
C                             OF ERROR.
C
C   REQD. IMSL ROUTINES - UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS  1.  ZREAL2 ASSUMES THAT THERE EXIST N DISTINCT REAL ROOTS
C                FOR THE FUNCTION F AND THAT THE INITIAL GUESSES
C                SUPPLIED BY THE USER ARE SUFFICIENTLY CLOSE TO THE
C                ROOTS TO OBTAIN CONVERGENCE BY NEWTONS METHOD.
C                THE ROUTINE IS DESIGNED SO THAT CONVERGENCE TO ANY
C                SINGLE ROOT CANNOT BE OBTAINED FROM TWO DIFFERENT
C                INITIAL GUESSES. THIS ROUTINE IS INTENDED PRIMARILY
C                FOR THE REFINEMENT OF N KNOWN ROUGH APPROXIMATIONS
C                OF THE ROOTS OF F.
C            2.  SCALING THE X VECTOR IN THE FUNCTION F MAY BE REQUIRED
C                IF ANY OF THE ROOTS ARE KNOWN TO BE LESS THAN ONE.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE ZREAL2 (F,EPS,EPS2,ETA,NSIG,N,X,ITMAX,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            NSIG,N,ITMAX,IER
      REAL               F,EPS,EPS2,ETA,X(N)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,J,NM1,IR,IC
      REAL               XI,AXI,FXI,AFXI,DI,HI,FXIPHI,SEPS,SINF,
     1                   SQREPS,DER,XIPI,ADER,ERR1,CRIT1,P1,ONE,TEN
      DATA               P1,ONE,TEN/.1,1.0,10.0/
      DATA               SEPS/Z3C100000/
      DATA               SINF/Z7FFFFFFF/
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
      IR=0
      CRIT1 = TEN**(-NSIG)
      SQREPS = SQRT(SEPS)
      DO 30 I=1,N
         IC = 1
         XI = X(I)
    5    AXI = ABS(XI)
         IF (I .EQ. 1) GO TO 15
         NM1=I-1
         DO 10 J = 1,NM1
            IF (ABS(XI - X(J)) .LT. EPS2) XI = XI + ETA
   10    CONTINUE
   15    FXI = F(XI)
         AFXI = ABS(FXI)
C                                  TEST FOR CONVERGENCE
         IF (AFXI .LE. EPS) GO TO 25
         DI = SQREPS
         IF (AXI .GE. P1) DI = SQREPS*AXI
         HI = AMIN1(AFXI,DI)
         FXIPHI = F(XI + HI)
         DER = (FXIPHI - FXI)/HI
         ADER = ABS(DER)
         IF (ADER .GE. ONE) GO TO 16
         IF (AFXI .GE. SINF*ADER) GO TO 20
   16    XIPI=FXI/DER
         XI=XI-XIPI
C                                  TEST FOR CONVERGENCE
         ERR1 = ABS(XIPI)/AMAX1(P1,AXI)
         IF(ERR1.LE.CRIT1) GO TO 25
         IC = IC + 1
         IF (IC .LE. ITMAX)  GO TO 5
C                                  ROOT NOT FOUND, NO CONVERGENCE
         X(I) = 111111.
         IR=IR+1
         IER=33
         GO TO 30
C                                  ROOT NOT FOUND, DERIVATIVE = 0.
   20    X(I) = 222222.
         IR=IR+1
         IER=34
         GO TO 30
   25    X(I)=XI
   30 CONTINUE
      ITMAX = IC
      IF(IER.EQ.0) GO TO 9005
      IF(IR.LE.1) GO TO 9000
      IER=35
 9000 CONTINUE
      CALL UERTST(IER,6HZREAL2)
 9005 RETURN
      END
