C   IMSL ROUTINE NAME   - DVERK
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1982
C
C   PURPOSE             - DIFFERENTIAL EQUATION SOLVER - RUNGE
C                           KUTTA-VERNER FIFTH AND SIXTH ORDER METHOD
C
C   USAGE               - CALL DVERK (N,FCN,X,Y,XEND,TOL,IND,C,NW,W,IER)
C
C   ARGUMENTS    N      - NUMBER OF EQUATIONS. (INPUT)
C                FCN    - NAME OF SUBROUTINE FOR EVALUATING FUNCTIONS.
C                           (INPUT)
C                           THE SUBROUTINE ITSELF MUST ALSO BE PROVIDED
C                             BY THE USER AND IT SHOULD BE OF THE
C                             FOLLOWING FORM
C                               SUBROUTINE FCN(N,X,Y,YPRIME)
C                               REAL Y(N),YPRIME(N)
C                                    .
C                                    .
C                                    .
C                           FCN SHOULD EVALUATE YPRIME(1),...,YPRIME(N)
C                             GIVEN N,X, AND Y(1),...,Y(N).  YPRIME(I)
C                             IS THE FIRST DERIVATIVE OF Y(I) WITH
C                             RESPECT TO X.
C                           FCN MUST APPEAR IN AN EXTERNAL STATEMENT IN
C                             THE CALLING PROGRAM AND N,X,Y(1),...,Y(N)
C                             MUST NOT BE ALTERED BY FCN.
C                X      - INDEPENDENT VARIABLE. (INPUT AND OUTPUT)
C                           ON INPUT, X SUPPLIES THE INITIAL VALUE.
C                           ON OUTPUT, X IS REPLACED WITH XEND UNLESS
C                             ERROR CONDITIONS ARISE.  SEE THE DES-
C                             CRIPTION OF PARAMETER IND.
C                Y      - DEPENDENT VARIABLES, VECTOR OF LENGTH N.
C                           (INPUT AND OUTPUT)
C                           ON INPUT, Y(1),...,Y(N) SUPPLY INITIAL
C                             VALUES.
C                           ON OUTPUT, Y(1),...,Y(N) ARE REPLACED WITH
C                             AN APPROXIMATE SOLUTION AT XEND UNLESS
C                             ERROR CONDITIONS ARISE.  SEE THE DES-
C                             CRIPTION OF PARAMETER IND.
C                XEND   - VALUE OF X AT WHICH SOLUTION IS DESIRED.
C                           (INPUT)
C                           XEND MAY BE LESS THAN THE INITIAL VALUE OF
C                             X.
C                TOL    - TOLERANCE FOR ERROR CONTROL. (INPUT)
C                           THE SUBROUTINE ATTEMPTS TO CONTROL A NORM
C                             OF THE LOCAL ERROR IN SUCH A WAY THAT THE
C                             GLOBAL ERROR IS PROPORTIONAL TO TOL.
C                             MAKING TOL SMALLER IMPROVES ACCURACY AND
C                             MORE THAN ONE RUN, WITH DIFFERENT VALUES
C                             OF TOL, CAN BE USED IN AN ATTEMPT TO
C                             ESTIMATE THE GLOBAL ERROR.
C                           IN THE DEFAULT CASE (IND=1), THE GLOBAL
C                             ERROR IS
C                                 MAX(ABS(E(1)),...,ABS(E(N)))
C                             WHERE E(K)=(Y(K)-YT(K))/MAX(1,ABS(Y(K)))
C                             YT(K) IS THE TRUE SOLUTION, AND
C                             Y(K) IS THE COMPUTED SOLUTION AT XEND,
C                             FOR K=1,2,...,N.
C                             OTHER ERROR CONTROL OPTIONS ARE AVAILABLE.
C                             SEE THE DESCRIPTION OF PARAMETERS IND AND
C                             C BELOW.
C                IND    - INDICATOR. (INPUT AND OUTPUT)
C                           ON INITIAL ENTRY IND MUST BE SET EQUAL TO
C                           EITHER 1 OR 2.
C                           IND = 1 CAUSES ALL DEFAULT OPTIONS TO BE
C                             USED AND ELIMINATES THE NEED TO SET
C                             SPECIFIC VALUES IN THE COMMUNICATIONS
C                             VECTOR C.
C                           IND = 2 ALLOWS OPTIONS TO BE SELECTED.  IN
C                             THIS CASE, THE FIRST 9 COMPONENTS OF C
C                             MUST BE INITIALIZED TO SELECT OPTIONS AS
C                             DESCRIBED BELOW.
C                         THE SUBROUTINE WILL NORMALLY RETURN WITH
C                           IND = 3, HAVING REPLACED THE INITIAL VALUES
C                           OF X AND Y WITH, RESPECTIVELY, THE VALUE
C                           XEND AND AN APPROXIMATION TO Y AT XEND.
C                         THE SUBROUTINE CAN BE CALLED REPEATEDLY WITH
C                           NEW VALUES OF XEND WITHOUT CHANGING ANY
C                           OF THE OTHER PARAMETERS.
C                         THREE ERROR RETURNS ARE ALSO POSSIBLE, IN
C                           WHICH CASE X AND Y WILL BE THE MOST
C                           RECENTLY ACCEPTED VALUES.
C                           IND = -3 INDICATES THAT THE SUBROUTINE WAS
C                             UNABLE TO SATISFY THE ERROR REQUIREMENT.
C                             THIS MAY MEAN THAT TOL IS TOO SMALL.
C                           IND = -2 INDICATES THAT THE VALUE OF HMIN
C                             (MINIMUM STEP-SIZE) IS GREATER THAN HMAX
C                             (MAXIMUM STEP-SIZE), WHICH PROBABLY MEANS
C                             THAT THE REQUESTED TOL (WHICH IS USED IN
C                             THE CALCULATION OF HMIN) IS TOO SMALL.
C                           IND = -1 INDICATES THAT THE ALLOWED MAXIMUM
C                             NUMBER OF FCN EVALUATIONS HAS BEEN
C                             EXCEEDED.  THIS CAN ONLY OCCUR IF OPTION
C                             C(7), AS DESCRIBED BELOW, HAS BEEN USED.
C                C      - COMMUNICATIONS VECTOR OF LENGTH 24. (INPUT IF
C                           IND.NE.1, AND OUTPUT).
C                           C IS USED TO SELECT OPTIONS AND TO RETAIN
C                             INFORMATION BETWEEN CALLS.  THE USER NEED
C                             NOT BE CONCERNED WITH THE FOLLOWING
C                             DESCRIPTION OF THE ELEMENTS OF C WHEN
C                             DEFAULT OPTIONS ARE USED (IND=1).
C                             HOWEVER, WHEN IT IS DESIRED TO USE IND=2
C                             AND SELECT OPTIONS, A BASIC UNDERSTANDING
C                             OF DVERK IS REQUIRED.  THE FOLLOWING
C                             PARAGRAPH DESCRIBES, BRIEFLY, THE BASIC
C                             TERMS.  FOR MORE DETAILS, SEE THE DOCUMENT
C                             REFERENCE.
C                             DVERK ADVANCES THE INDEPENDENT VARIABLE
C                             X ONE STEP AT A TIME UNTIL XEND IS
C                             REACHED.  THE SOLUTION IS COMPUTED AT
C                             XTRIAL = X+HTRIAL ALONG WITH AN ERROR
C                             ESTIMATE EST.  IF EST IS LESS THAN OR
C                             EQUAL TO TOL (SUCCESSFUL STEP), THE STEP
C                             IS ACCEPTED AND X IS ADVANCED TO XTRIAL.
C                             IF EST IS GREATER THAN TOL (FAILURE)
C                             HTRIAL IS ADJUSTED AND THE SOLUTION IS
C                             RECOMPUTED.  HMAG = ABS(HTRIAL) IS NEVER
C                             ALLOWED TO EXCEED HMAX NOR IS IT ALLOWED
C                             TO BECOME SMALLER THAN HMIN.  THE FIRST
C                             TRIAL STEP IS HSTART.  DURING THE
C                             COMPUTATION, A COUNTER (C(23)) IS
C                             INCREMENTED EACH TIME A TRIAL STEP FAILS
C                             TO PROVIDE A SOLUTION SATISFYING THE ERROR
C                             TOLERANCE.  ANOTHER COUNTER (C(22)) IS
C                             USED TO RECORD THE NUMBER OF SUCCESSFUL
C                             STEPS.  AFTER A SUCCESSFUL STEP, C(23) IS
C                             SET TO ZERO.
C                         OPTIONS.  IF THE SUBROUTINE IS ENTERED WITH
C                           IND=2, THE FIRST 9 COMPONENTS OF THE
C                           COMMUNICATIONS VECTOR MUST BE INITIALIZED
C                           BY THE USER.  NORMALLY THIS IS DONE BY
C                           FIRST SETTING THEM ALL TO ZERO, AND THEN
C                           THOSE CORRESPONDING TO PARTICULAR OPTIONS
C                           ARE MADE NON-ZERO.
C                C(1)   - ERROR CONTROL INDICATOR.
C                           THE SUBROUTINE ATTEMPTS TO CONTROL A NORM
C                           OF THE LOCAL ERROR IN SUCH A WAY THAT THE
C                           GLOBAL ERROR IS PROPORTIONAL TO TOL.
C                           THE DEFINITION OF GLOBAL ERROR FOR THE
C                           DEFAULT CASE (IND=1) IS GIVEN IN THE
C                           DESCRIPTION OF PARAMETER TOL.  THE DEFAULT
C                           WEIGHTS ARE 1/MAX(1,ABS(Y(K))).  WHEN IND=2
C                           IS USED, THE WEIGHTS ARE DETERMINED
C                           ACCORDING TO THE VALUE OF C(1).
C                           IF C(1)=1 THE WEIGHTS ARE 1
C                                     (ABSOLUTE ERROR CONTROL)
C                           IF C(1)=2 THE WEIGHTS ARE 1/ABS(Y(K))
C                                     FOR K=1,2,...,N.
C                                     (RELATIVE ERROR CONTROL)
C                           IF C(1)=3 THE WEIGHTS ARE
C                                     1/MAX(ABS(C(2)),ABS(Y(K)))
C                                     FOR K=1,2,...,N.
C                                     (RELATIVE ERROR CONTROL, UNLESS
C                                     ABS(Y(K)) IS LESS THAN THE FLOOR
C                                     VALUE,ABS(C(2)))
C                           IF C(1)=4 THE WEIGHTS ARE
C                                     1/MAX(ABS(C(K+30)),ABS(Y(K)))
C                                     FOR K=1,2,...,N.
C                                     (HERE INDIVIDUAL FLOOR VALUES
C                                     ARE USED)
C                                     IN THIS CASE,  THE DIMENSION OF C
C                                     MUST BE GREATER THAN OR EQUAL TO
C                                     N+30 AND C(31), C(32),...,C(N+30)
C                                     MUST BE INITIALIZED BY THE USER.
C                           IF C(1)=5 THE WEIGHTS ARE 1/ABS(C(K+30))
C                                     FOR K=1,2,...,N.
C                                     IN THIS CASE,  THE DIMENSION OF C
C                                     MUST BE GREATER THAN OR EQUAL TO
C                                     N+30 AND C(31), C(32),...,C(N+30)
C                                     MUST BE INITIALIZED BY THE USER.
C                           FOR ALL OTHER VALUES OF C(1), INCLUDING
C                              C(1)=0 THE DEFAULT VALUES OF
C                                     THE WEIGHTS ARE TAKEN TO BE
C                                     1/MAX(1,ABS(Y(K)))
C                                     FOR K=1,2,...,N.
C                C(2)   - FLOOR VALUE.  USED WHEN THE INDICATOR C(1)
C                           HAS THE VALUE 3.
C                C(3)   - HMIN SPECIFICATION.  IF NOT ZERO, THE SUB-
C                           ROUTINE CHOOSES HMIN TO BE ABS(C(3)).
C                           OTHERWISE IT USES THE DEFAULT VALUE
C                           10*MAX(DWARF,RREB*MAX(NORM(Y)/TOL,ABS(X)))
C                           WHERE DWARF IS A VERY SMALL POSITIVE MACHINE
C                           NUMBER AND RREB IS THE RELATIVE ROUNDOFF
C                           ERROR BOUND.
C                C(4)   - HSTART SPECIFICATION.  IF NOT ZERO, THE SUB-
C                           ROUTINE WILL USE AN INITIAL HMAG EQUAL TO
C                           ABS(C(4)), EXCEPT OF COURSE FOR THE RE-
C                           STRICTIONS IMPOSED BY HMIN AND HMAX.
C                           OTHERWISE IT USES THE DEFAULT VALUE
C                             HMAX*(TOL)**(1/6).
C                C(5)   - SCALE SPECIFICATION.  THIS IS INTENDED TO BE
C                           A MEASURE OF THE SCALE OF THE PROBLEM.
C                           LARGER VALUES OF SCALE TEND TO MAKE THE
C                           METHOD MORE RELIABLE, FIRST BY POSSIBLY RE-
C                           STRICTING HMAX (AS DESCRIBED BELOW) AND
C                           SECOND, BY TIGHTENING THE ACCEPTANCE
C                           REQUIREMENT.  IF C(5) IS ZERO, A DEFAULT
C                           VALUE OF 1 IS USED. FOR LINEAR HOMOGENEOUS
C                           PROBLEMS WITH CONSTANT COEFFICIENTS, AN
C                           APPROPRIATE VALUE FOR SCALE IS A NORM OF
C                           THE ASSOCIATED MATRIX. FOR OTHER PROBLEMS,
C                           AN APPROXIMATION TO AN AVERAGE VALUE OF A
C                           NORM OF THE JACOBIAN ALONG THE TRAJEC-
C                           TORY MAY BE APPROPRIATE.
C                C(6)   - HMAX SPECIFICATION.  FOUR CASES ARE POSSIBLE,
C                           IF C(6).NE.0 AND C(5).NE.0, HMAX IS TAKEN
C                             TO BE MIN(ABS(C(6)),2/ABS(C(5))).
C                           IF C(6).NE.0 AND C(5).EQ.0, HMAX IS TAKEN
C                             TO BE ABS(C(6)).
C                           IF C(6).EQ.0 AND C(5).NE.0, HMAX IS TAKEN
C                             TO BE 2/ABS(C(5)).
C                           IF C(6).EQ.0 AND C(5).EQ.0, HMAX IS GIVEN
C                             A DEFAULT VALUE OF 2.
C                C(7)   - MAXIMUM NUMBER OF FUNCTION EVALUATIONS.  IF
C                           NOT ZERO, AN ERROR RETURN WITH IND = -1
C                           WILL BE CAUSED WHEN THE NUMBER OF FUNCTION
C                           EVALUATIONS EXCEEDS ABS(C(7)).
C                C(8)   - INTERRUPT NUMBER 1 . IF NOT ZERO, THE SUB-
C                           ROUTINE WILL INTERRUPT THE CALCULATIONS
C                           AFTER IT HAS CHOSEN ITS PRELIMINARY VALUE
C                           OF HMAG, AND JUST BEFORE CHOOSING HTRIAL
C                           AND XTRIAL IN PREPARATION FOR TAKING A STEP
C                           (HTRIAL MAY DIFFER FROM HMAG IN SIGN, AND
C                           MAY REQUIRE ADJUSTMENT IF XEND IS NEAR).
C                           THE SUBROUTINE RETURNS WITH IND = 4, AND
C                           WILL RESUME CALCULATION AT THE POINT OF
C                           INTERRUPTION IF  RE-ENTERED WITH IND = 4.
C                C(9)   - INTERRUPT NUMBER 2.  IF NOT ZERO, THE SUB-
C                           ROUTINE WILL INTERRUPT THE CALCULATIONS
C                           IMMEDIATELY AFTER IT HAS DECIDED WHETHER OR
C                           NOT TO ACCEPT THE RESULT OF THE MOST RECENT
C                           TRIAL STEP, WITH IND = 5 IF IT PLANS TO
C                           ACCEPT, OR IND = 6 IF IT PLANS TO REJECT.
C                           Y(*) IS THE PREVIOUSLY ACCEPTED RESULT,
C                           WHILE W(*,9) IS THE NEWLY COMPUTED TRIAL
C                           VALUE, AND W(*,2) IS THE UNWEIGHTED ERROR
C                           ESTIMATE VECTOR. THE SUBROUTINE WILL RESUME
C                           CALCULATIONS AT THE POINT OF INTERRUPTION
C                           ON RE-ENTRY WITH IND = 5 OR 6.
C                           IND MAY BE CHANGED BY THE USER IN ORDER TO
C                           FORCE ACCEPTANCE OF A STEP (BY CHANGING IND
C                           FROM 6 TO 5) THAT WOULD OTHERWISE BE
C                           REJECTED, OR VICE VERSA.
C                NW     - ROW DIMENSION OF THE MATRIX W EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT
C                           IN THE CALLING PROGRAM. (INPUT)
C                           NW MUST BE GREATER THAN OR EQUAL TO N.
C                W      - WORKSPACE MATRIX.
C                           THE FIRST DIMENSION OF W MUST BE NW AND THE
C                           SECOND MUST BE GREATER THAN OR EQUAL TO 9.
C                           W MUST REMAIN UNCHANGED BETWEEN SUCCESSIVE
C                           CALLS DURING INTEGRATION.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER = 129, NW IS LESS THAN N OR TOL IS LESS
C                             THAN OR EQUAL TO ZERO.
C                           IER = 130, IND IS NOT IN THE RANGE 1 TO 6.
C                           IER = 131, XEND HAS NOT BEEN CHANGED FROM
C                             PREVIOUS CALL OR X IS NOT SET TO
C                             THE PREVIOUS XEND VALUE.
C                           IER = 132, THE RELATIVE ERROR CONTROL
C                             OPTION (C(1)=2) WAS SELECTED AND
C                             ONE OF THE SOLUTION COMPONENTS
C                             IS ZERO.
C
C   REQD. IMSL ROUTINES - UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS  1.  IN A TYPICAL SITUATION, DVERK IS CALLED
C                REPEATEDLY WITH A SEQUENCE OF VALUES FOR XEND.
C                AFTER EACH SUCH CALL, THE USER SHOULD INTERROGATE
C                IND AND IER. ERROR CONDITIONS ARE SIGNALED WHEN
C                IND IS LESS THAN ZERO AND/OR IER IS GREATER THAN
C                ZERO. CORRECTIVE ACTION (SUCH AS CHANGING CERTAIN
C                PARAMETER VALUES) MUST BE TAKEN PRIOR TO RE-ENTRY.
C            2.  WHEN ERROR CONDITIONS ARISE, IT IS OFTEN HELPFUL
C                TO EXAMINE COMPONENTS OF THE COMMUNICATIONS VECTOR
C                C. A SUMMARY FOLLOWS-
C
C                PRESCRIBED AT THE OPTION OF THE USER
C
C                  C(1) ERROR CONTROL INDICATOR
C                  C(2) FLOOR VALUE
C                  C(3) HMIN SPECIFICATION
C                  C(4) HSTART SPECIFICATION
C                  C(5) SCALE SPECIFICATION
C                  C(6) HMAX SPECIFICATION
C                  C(7) MAXIMUM NUMBER OF FCN EVALUATIONS
C                  C(8) INTERRUPT NUMBER 1
C                  C(9) INTERRUPT NUMBER 2
C
C                DETERMINED BY THE PROGRAM
C
C                  C(10) RREB (RELATIVE ROUNDOFF ERROR BOUND)
C                  C(11) DWARF (VERY SMALL MACHINE NUMBER)
C                  C(12) WEIGHTED NORM OF Y
C                  C(13) HMIN
C                  C(14) HMAG
C                  C(15) SCALE
C                  C(16) HMAX
C                  C(17) XTRIAL
C                  C(18) HTRIAL
C                  C(19) EST
C                  C(20) PREVIOUS XEND
C                  C(21) FLAG FOR XEND
C                  C(22) NUMBER OF SUCCESSFUL STEPS
C                  C(23) NUMBER OF SUCCESSIVE FAILURES
C                  C(24) NUMBER OF FCN EVALUATIONS
C
C                IF C(1) = 4 OR 5, C(31),C(32),...,C(N+30) ARE FLOOR
C                VALUES.
C            3.  PARAMETER NW GIVES THE ROW DIMENSION OF W EXACTLY AS
C                IT APPEARS IN THE DIMENSION STATEMENT IN THE CALLING
C                PROGRAM. IF ONLY ONE SYSTEM OF EQUATIONS IS BEING
C                SOLVED, NW NORMALLY WILL HAVE THE SAME VALUE AS N.
C                HOWEVER, IF MORE THAN ONE SYSTEM IS BEING HANDLED,
C                AND THEY ARE TO USE A COMMON WORKSPACE, W, ONE AFTER
C                THE OTHER, THE VALUE OF NW (AND HENCE, THE ROW
C                DIMENSION OF W IN THE CALLING PROGRAM) MUST BE AS
C                LARGE AS THE MAXIMUM VALUE OF THE INDIVIDUAL N VALUES.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE DVERK (N,FCN,X,Y,XEND,TOL,IND,C,NW,W,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N,IND,NW,IER
      REAL               X,Y(N),XEND,TOL,C(1),W(NW,9)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            K
      REAL               ZERO,ONE,TWO,THREE,FOUR,FIVE,SEVEN,TEN,HALF,P9
      REAL               C4D15,C2D3,C5D6,C1D6,C1D15,C2D96,TEMP
      REAL               RK(39),REPS,RTOL
      DATA               ZERO/0.0/,ONE/1.0/,TWO/2.0/,THREE/3.0/
      DATA               FOUR/4.0/,FIVE/5.0/,SEVEN/7.0/
      DATA               TEN/10.0/,HALF/0.5/,P9/0.9/
      DATA               C4D15/.2666667/
      DATA               C2D3/.6666667/
      DATA               C5D6/.8333333/
      DATA               C1D6/.1666667/
      DATA               C1D15/.6666667E-1/
      DATA               C2D96/120.4273/
      DATA               REPS/Z3C100000/
      DATA               RTOL/Z05100000/
      DATA               RK( 1)/.1666667E+00/
      DATA               RK( 2)/.5333333E-01/
      DATA               RK( 3)/.2133333E+00/
      DATA               RK( 4)/.8333333E+00/
      DATA               RK( 5)/.2666667E+01/
      DATA               RK( 6)/.2500000E+01/
      DATA               RK( 7)/.2578125E+01/
      DATA               RK( 8)/.9166667E+01/
      DATA               RK( 9)/.6640625E+01/
      DATA               RK(10)/.8854167E+00/
      DATA               RK(11)/.2400000E+01/
      DATA               RK(12)/.8000000E+01/
      DATA               RK(13)/.6560458E+01/
      DATA               RK(14)/.3055556E+00/
      DATA               RK(15)/.3450980E+00/
      DATA               RK(16)/.5508667E+00/
      DATA               RK(17)/.1653333E+01/
      DATA               RK(18)/.9455882E+00/
      DATA               RK(19)/.3240000E+00/
      DATA               RK(20)/.2337882E+00/
      DATA               RK(21)/.2035465E+01/
      DATA               RK(22)/.6976744E+01/
      DATA               RK(23)/.5648180E+01/
      DATA               RK(24)/.1373816E+00/
      DATA               RK(25)/.2863023E+00/
      DATA               RK(26)/.1441786E+00/
      DATA               RK(27)/.7500000E-01/
      DATA               RK(28)/.3899287E+00/
      DATA               RK(29)/.3194444E+00/
      DATA               RK(30)/.1350384E+00/
      DATA               RK(31)/.1078330E-01/
      DATA               RK(32)/.6980519E-01/
      DATA               RK(33)/.6250000E-02/
      DATA               RK(34)/.6963012E-02/
      DATA               RK(35)/.6944444E-02/
      DATA               RK(36)/.6138107E-02/
      DATA               RK(37)/.6818182E-01/
      DATA               RK(38)/.1078330E-01/
      DATA               RK(39)/.6980519E-01/
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
C                                  BEGIN INITIALIZATION, PARAMETER
C                                    CHECKING, INTERRUPT RE-ENTRIES
C                                  ABORT IF IND OUT OF RANGE 1 TO 6
      IF (IND.LT.1.OR.IND.GT.6) GO TO 290
C                                  CASES - INITIAL ENTRY, NORMAL
C                                    RE-ENTRY, INTERRUPT RE-ENTRIES
      GO TO (5,5,40,145,265,265), IND
C                                  CASE 1 - INITIAL ENTRY (IND .EQ. 1
C                                    OR 2) ABORT IF N.GT.NW OR TOL.LE.0
C
    5 IF (N.GT.NW.OR.TOL.LE.ZERO) GO TO 295
      IF (IND.EQ.2) GO TO 15
C                                  INITIAL ENTRY WITHOUT OPTIONS (IND
C                                    .EQ. 1) SET C(1) TO C(9) EQUAL TO
C                                    0
      DO 10 K=1,9
         C(K) = ZERO
   10 CONTINUE
      GO TO 30
C                                  SUMMARY OF THE COMPONENTS OF THE
C                                    COMMUNICATIONS VECTOR
C                                    PRESCRIBED AT THE OPTION
C                                        OF THE USER
C
C                                    C(1) ERROR CONTROL INDICATOR
C                                    C(2) FLOOR VALUE
C                                    C(3) HMIN SPECIFICATION
C                                    C(4) HSTART SPECIFICATION
C                                    C(5) SCALE SPECIFICATION
C                                    C(6) HMAX SPECIFICATION
C                                    C(7) MAX NO OF FCN EVALS
C                                    C(8) INTERRUPT NO 1
C                                    C(9) INTERRUPT NO 2
C
C                                    DETERMINED BY THE PROGRAM
C
C                                    C(10) RREB(REL ROUNDOFF ERROR BND)
C                                    C(11) DWARF (VERY SMALL MACH NO)
C                                    C(12) WEIGHTED NORM Y
C                                    C(13) HMIN
C                                    C(14) HMAG
C                                    C(15) SCALE
C                                    C(16) HMAX
C                                    C(17) XTRIAL
C                                    C(18) HTRIAL
C                                    C(19) EST
C                                    C(20) PREVIOUS XEND
C                                    C(21) FLAG FOR XEND
C                                    C(22) NO OF SUCCESSFUL STEPS
C                                    C(23) NO OF SUCCESSIVE FAILURES
C                                    C(24) NO OF FCN EVALS
C                                    IF C(1) = 4 OR 5, C(31),C(32),...
C                                      C(N+30) ARE FLOOR VALUES
   15 CONTINUE
C                                  INITIAL ENTRY WITH OPTIONS (IND .EQ.
C                                    2) MAKE C(1) TO C(9) NON-NEGATIVE
      DO 20 K=1,9
         C(K) = ABS(C(K))
   20 CONTINUE
C                                  MAKE FLOOR VALUES NON-NEGATIVE IF
C                                    THEY ARE TO BE USED
C
      IF (C(1).NE.FOUR.AND.C(1).NE.FIVE) GO TO 30
      DO 25 K=1,N
         C(K+30) = ABS(C(K+30))
   25 CONTINUE
   30 CONTINUE
C                                  INITIALIZE RREB, DWARF, PREV XEND,
C                                    FLAG, COUNTS
      C(10) = REPS
      C(11) = RTOL
C                                  SET PREVIOUS XEND INITIALLY TO
C                                    INITIAL VALUE OF X
      C(20) = X
      DO 35 K=21,24
         C(K) = ZERO
   35 CONTINUE
      GO TO 45
C                                  CASE 2 - NORMAL RE-ENTRY (IND .EQ.
C                                    3) ABORT IF XEND REACHED, AND
C                                    EITHER X CHANGED OR XEND NOT
   40 IF (C(21).NE.ZERO.AND.(X.NE.C(20).OR.XEND.EQ.C(20))) GO TO 285
C
C                                  RE-INITIALIZE FLAG
      C(21) = ZERO
C                                  CASE 3 - RE-ENTRY FOLLOWING AN
C                                    INTERRUPT (IND .EQ. 4 TO 6)
C                                    TRANSFER CONTROL TO THE
C                                    APPROPRIATE RE-ENTRY POINT. THIS
C                                    HAS ALREADY BEEN HANDLED BY THE
C                                    COMPUTED GO TO.
   45 CONTINUE
C                                  LOOP THROUGH THE FOLLOWING 4 STAGES,
C                                    ONCE FOR EACH TRIAL STEP UNTIL THE
C                                    OCCURRENCE OF ONE OF THE FOLLOWING
C                                    (A) THE NORMAL RETURN (WITH IND
C                                    .EQ. 3) ON REACHING XEND IN STAGE
C                                    4 (B) AN ERROR RETURN (WITH IND
C                                    .LT. 0) IN STAGE 1 OR STAGE 4 (C)
C                                    AN INTERRUPT RETURN (WITH IND .EQ.
C                                    4, 5 OR 6), IF REQUESTED, IN STAGE
C                                    1 OR STAGE 4
   50 CONTINUE
C                                  STAGE 1 - PREPARE - DO CALCULATIONS
C                                    OF HMIN, HMAX, ETC., AND SOME
C                                    PARAMETER CHECKING, AND END UP
C                                    WITH SUITABLE VALUES OF HMAG,
C                                    XTRIAL AND HTRIAL IN PREPARATION
C                                    FOR TAKING AN INTEGRATION STEP.
C                                    ERROR RETURN (WITH IND=-1) IF NO
C                                    OF FCN EVALS TOO GREAT
C
      IF (C(7).EQ.ZERO.OR.C(24).LT.C(7)) GO TO 55
      IND = -1
      GO TO 9005
   55 CONTINUE
C                                  CALCULATE SLOPE (ADDING 1 TO NO OF
C                                    FCN EVALS) IF IND .NE. 6
      IF (IND.EQ.6) GO TO 60
      CALL FCN (N,X,Y,W(1,1))
      C(24) = C(24)+ONE
   60 CONTINUE
C                                  CALCULATE HMIN - USE DEFAULT UNLESS
C                                    VALUE PRESCRIBED
      C(13) = C(3)
      IF (C(3).NE.ZERO) GO TO 120
C                                  CALCULATE DEFAULT VALUE OF HMIN
C                                    FIRST CALCULATE WEIGHTED NORM Y -
C                                    C(12) - AS SPECIFIED BY THE ERROR
C                                    CONTROL INDICATOR C(1)
      TEMP = ZERO
      IF (C(1).NE.ONE) GO TO 70
C                                  ABSOLUTE ERROR CONTROL - WEIGHTS ARE
C                                    1
      DO 65 K=1,N
         TEMP = AMAX1(TEMP,ABS(Y(K)))
   65 CONTINUE
      C(12) = TEMP
      GO TO 115
   70 IF (C(1).NE.TWO) GO TO 75
C                                  RELATIVE ERROR CONTROL - WEIGHTS ARE
C                                    1/DABS(Y(K)) SO WEIGHTED NORM Y IS
C                                    1
      C(12) = ONE
      GO TO 115
   75 IF (C(1).NE.THREE) GO TO 85
C                                  WEIGHTS ARE 1/MAX(C(2),ABS(Y(K)))
      DO 80 K=1,N
         TEMP = AMAX1(TEMP,ABS(Y(K))/C(2))
   80 CONTINUE
      C(12) = AMIN1(TEMP,ONE)
      GO TO 115
   85 IF (C(1).NE.FOUR) GO TO 95
C                                  WEIGHTS ARE 1/MAX(C(K+30),ABS(Y(K)))
      DO 90 K=1,N
         TEMP = AMAX1(TEMP,ABS(Y(K))/C(K+30))
   90 CONTINUE
      C(12) = AMIN1(TEMP,ONE)
      GO TO 115
   95 IF (C(1).NE.FIVE) GO TO 105
C                                  WEIGHTS ARE 1/C(K+30)
      DO 100 K=1,N
         TEMP = AMAX1(TEMP,ABS(Y(K))/C(K+30))
  100 CONTINUE
      C(12) = TEMP
      GO TO 115
  105 CONTINUE
C                                  DEFAULT CASE - WEIGHTS ARE
C                                    1/MAX(1,ABS(Y(K)))
      DO 110 K=1,N
         TEMP = AMAX1(TEMP,ABS(Y(K)))
  110 CONTINUE
      C(12) = AMIN1(TEMP,ONE)
  115 CONTINUE
      C(13) = TEN*AMAX1(C(11),C(10)*AMAX1(C(12)/TOL,ABS(X)))
  120 CONTINUE
C                                  CALCULATE SCALE - USE DEFAULT UNLESS
C                                    VALUE PRESCRIBED
      C(15) = C(5)
      IF (C(5).EQ.ZERO) C(15) = ONE
C                                  CALCULATE HMAX - CONSIDER 4 CASES
C                                    CASE 1 - BOTH HMAX AND SCALE
C                                    PRESCRIBED
C
      IF (C(6).NE.ZERO.AND.C(5).NE.ZERO) C(16) = AMIN1(C(6),TWO/C(5))
C
C                                  CASE 2 - HMAX PRESCRIBED, BUT SCALE
C                                    NOT
C
      IF (C(6).NE.ZERO.AND.C(5).EQ.ZERO) C(16) = C(6)
C
C                                  CASE 3 - HMAX NOT PRESCRIBED, BUT
C                                    SCALE IS
C
      IF (C(6).EQ.ZERO.AND.C(5).NE.ZERO) C(16) = TWO/C(5)
C
C                                  CASE 4 - NEITHER HMAX NOR SCALE IS
C                                    PROVIDED
C
      IF (C(6).EQ.ZERO.AND.C(5).EQ.ZERO) C(16) = TWO
C
C                                  ERROR RETURN (WITH IND=-2) IF HMIN
C                                    .GT. HMAX
      IF (C(13).LE.C(16)) GO TO 125
      IND = -2
      GO TO 9005
  125 CONTINUE
C                                  CALCULATE PRELIMINARY HMAG -
C                                    CONSIDER 3 CASES
      IF (IND.GT.2) GO TO 130
C                                  CASE 1 - INITIAL ENTRY - USE
C                                    PRESCRIBED VALUE OF HSTART, IF
C                                    ANY, ELSE DEFAULT
      C(14) = C(4)
      IF (C(4).EQ.ZERO) C(14) = C(16)*TOL**C1D6
      GO TO 140
  130 IF (C(23).GT.ONE) GO TO 135
C                                  CASE 2 - AFTER A SUCCESSFUL STEP, OR
C                                    AT MOST ONE FAILURE, USE MIN(2,
C                                    .9*(TOL/EST)**(1/6))*HMAG, BUT
C                                    AVOID POSSIBLE OVERFLOW. THEN
C                                    AVOID REDUCTION BY MORE THAN HALF.
      TEMP = TWO*C(14)
      IF (TOL.LT.C2D96*C(19)) TEMP = P9*(TOL/C(19))**C1D6*C(14)
      C(14) = AMAX1(TEMP,HALF*C(14))
      GO TO 140
  135 CONTINUE
C                                  CASE 3 - AFTER TWO OR MORE
C                                    SUCCESSIVE FAILURES
      C(14) = HALF*C(14)
  140 CONTINUE
C                                  CHECK AGAINST HMAX
      C(14) = AMIN1(C(14),C(16))
C                                  CHECK AGAINST HMIN
      C(14) = AMAX1(C(14),C(13))
C                                  INTERRUPT NO 1 (WITH IND=4) IF
C                                    REQUESTED
      IF (C(8).EQ.ZERO) GO TO 145
      IND = 4
      GO TO 9005
C                                  RESUME HERE ON RE-ENTRY WITH IND
C                                    .EQ. 4
  145 CONTINUE
C                                  CALCULATE HMAG, XTRIAL - DEPENDING
C                                    ON PRELIMINARY HMAG, XEND
C
      IF (C(14).GE.ABS(XEND-X)) GO TO 150
C
C                                  DO NOT STEP MORE THAN HALF WAY TO
C                                    XEND
C
      C(14) = AMIN1(C(14),HALF*ABS(XEND-X))
      C(17) = X+SIGN(C(14),XEND-X)
      GO TO 155
  150 CONTINUE
C                                  HIT XEND EXACTLY
      C(14) = ABS(XEND-X)
      C(17) = XEND
  155 CONTINUE
C                                  CALCULATE HTRIAL
      C(18) = C(17)-X
C                                  STAGE 2 - CALCULATE YTRIAL (ADDING 7
C                                    TO NO OF FCN EVALS). W(*,2), ...
C                                    W(*,8) HOLD INTERMEDIATE RESULTS
C                                    NEEDED IN STAGE 3. W(*,9) IS
C                                    TEMPORARY STORAGE UNTIL FINALLY IT
C                                    HOLDS YTRIAL.
      DO 160 K=1,N
         W(K,9) = Y(K)+C(18)*W(K,1)*RK(1)
  160 CONTINUE
      CALL FCN (N,X+C(18)*C1D6,W(1,9),W(1,2))
C
      DO 165 K=1,N
         W(K,9) = Y(K)+C(18)*(W(K,1)*RK(2)+W(K,2)*RK(3))
  165 CONTINUE
      CALL FCN (N,X+C(18)*C4D15,W(1,9),W(1,3))
C
      DO 170 K=1,N
         W(K,9) = Y(K)+C(18)*(W(K,1)*RK(4)-W(K,2)*RK(5)+W(K,3)*RK(6))
  170 CONTINUE
      CALL FCN (N,X+C(18)*C2D3,W(1,9),W(1,4))
C
      DO 175 K=1,N
         W(K,9) = Y(K)+C(18)*(-W(K,1)*RK(7)+W(K,2)*RK(8)-W(K,3)*RK(9)
     1   +W(K,4)*RK(10))
  175 CONTINUE
      CALL FCN (N,X+C(18)*C5D6,W(1,9),W(1,5))
C
      DO 180 K=1,N
         W(K,9) = Y(K)+C(18)*(W(K,1)*RK(11)-W(K,2)*RK(12)+W(K,3)*RK(13)
     1   -W(K,4)*RK(14)+W(K,5)*RK(15))
  180 CONTINUE
      CALL FCN (N,X+C(18),W(1,9),W(1,6))
C
      DO 185 K=1,N
         W(K,9) = Y(K)+C(18)*(-W(K,1)*RK(16)+W(K,2)*RK(17)-W(K,3)
     1   *RK(18)-W(K,4)*RK(19)+W(K,5)*RK(20))
  185 CONTINUE
      CALL FCN (N,X+C(18)*C1D15,W(1,9),W(1,7))
C
      DO 190 K=1,N
         W(K,9) = Y(K)+C(18)*(W(K,1)*RK(21)-W(K,2)*RK(22)+W(K,3)*RK(23)
     1   -W(K,4)*RK(24)+W(K,5)*RK(25)+W(K,7)*RK(26))
  190 CONTINUE
      CALL FCN (N,X+C(18),W(1,9),W(1,8))
C                                  CALCULATE YTRIAL, THE EXTRAPOLATED
C                                    APPROXIMATION AND STORE IN W(*,9)
      DO 195 K=1,N
         W(K,9) = Y(K)+C(18)*(W(K,1)*RK(27)+W(K,3)*RK(28)+W(K,4)*RK(29)
     1   +W(K,5)*RK(30)+W(K,7)*RK(31)+W(K,8)*RK(32))
  195 CONTINUE
C                                  ADD 7 TO THE NO OF FCN EVALS
      C(24) = C(24)+SEVEN
C                                  STAGE 3 - CALCULATE THE ERROR
C                                    ESTIMATE EST. FIRST CALCULATE THE
C                                    UNWEIGHTED ABSOLUTE ERROR ESTIMATE
C                                    VECTOR (PER UNIT STEP) FOR THE
C                                    UNEXTRAPOLATED APPROXIMATION AND
C                                    STORE IT IN W(*,2). THEN CALCULATE
C                                    THE WEIGHTED MAX NORM OF W(*,2) AS
C                                    SPECIFIED BY THE ERROR CONTROL
C                                    INDICATOR C(1). FINALLY, MODIFY
C                                    THIS RESULT TO PRODUCE EST, THE
C                                    ERROR ESTIMATE (PER UNIT STEP) FOR
C                                    THE EXTRAPOLATED APPROXIMATION
C                                    YTRIAL.
C
C                                  CALCULATE THE UNWEIGHTED ABSOLUTE
C                                    ERROR ESTIMATE VECTOR.
      DO 200 K=1,N
         W(K,2) = W(K,1)*RK(33)+W(K,3)*RK(34)-W(K,4)*RK(35)+W(K,5)
     1   *RK(36)+W(K,6)*RK(37)-W(K,7)*RK(38)-W(K,8)*RK(39)
  200 CONTINUE
C                                  CALCULATE THE WEIGHTED MAX NORM OF
C                                    W(*,2) AS SPECIFIED BY THE ERROR
C                                    CONTROL INDICATOR C(1)
      TEMP = ZERO
      IF (C(1).NE.ONE) GO TO 210
C                                  ABSOLUTE ERROR CONTROL
      DO 205 K=1,N
         TEMP = AMAX1(TEMP,ABS(W(K,2)))
  205 CONTINUE
      GO TO 260
  210 IF (C(1).NE.TWO) GO TO 220
C                                  RELATIVE ERROR CONTROL
      DO 215 K=1,N
         IF (Y(K).EQ.ZERO) GO TO 280
         TEMP = AMAX1(TEMP,ABS(W(K,2)/Y(K)))
  215 CONTINUE
      GO TO 260
  220 IF (C(1).NE.THREE) GO TO 230
C                                  WEIGHTS ARE 1/MAX(C(2),ABS(Y(K)))
      DO 225 K=1,N
         TEMP = AMAX1(TEMP,ABS(W(K,2))/AMAX1(C(2),ABS(Y(K))))
  225 CONTINUE
      GO TO 260
  230 IF (C(1).NE.FOUR) GO TO 240
C                                  WEIGHTS ARE 1/MAX(C(K+30),ABS(Y(K)))
      DO 235 K=1,N
         TEMP = AMAX1(TEMP,ABS(W(K,2))/AMAX1(C(K+30),ABS(Y(K))))
  235 CONTINUE
      GO TO 260
  240 IF (C(1).NE.FIVE) GO TO 250
C                                  WEIGHTS ARE 1/C(K+30)
      DO 245 K=1,N
         TEMP = AMAX1(TEMP,ABS(W(K,2)/C(K+30)))
  245 CONTINUE
      GO TO 260
  250 CONTINUE
C                                  DEFAULT CASE - WEIGHTS ARE
C                                    1/MAX(1,ABS(Y(K)))
      DO 255 K=1,N
         TEMP = AMAX1(TEMP,ABS(W(K,2))/AMAX1(ONE,ABS(Y(K))))
  255 CONTINUE
  260 CONTINUE
C                                  CALCULATE EST - (THE WEIGHTED MAX
C                                    NORM OF W(*,2))*HMAG*SCALE - EST
C                                    IS INTENDED TO BE A MEASURE OF THE
C                                    ERROR PER UNIT STEP IN YTRIAL
      C(19) = TEMP*C(14)*C(15)
C                                  STAGE 4 - MAKE DECISIONS. SET IND=5
C                                    IF STEP ACCEPTABLE, ELSE SET IND=6
      IND = 5
      IF (C(19).GT.TOL) IND = 6
C                                  INTERRUPT NO 2 IF REQUESTED
      IF (C(9).NE.ZERO) GO TO 9005
C                                  RESUME HERE ON RE-ENTRY WITH IND
C                                    .EQ. 5 OR 6
  265 CONTINUE
C
      IF (IND.EQ.6) GO TO 275
C                                  STEP ACCEPTED (IND .EQ. 5), SO
C                                    UPDATE X, Y FROM XTRIAL, YTRIAL,
C                                    ADD 1 TO THE NO OF SUCCESSFUL
C                                    STEPS, AND SET THE NO OF
C                                    SUCCESSIVE FAILURES TO ZERO
      X = C(17)
      DO 270 K=1,N
         Y(K) = W(K,9)
  270 CONTINUE
      C(22) = C(22)+ONE
      C(23) = ZERO
C                                  RETURN(WITH IND=3, XEND SAVED, FLAG
C                                    SET) IF X .EQ. XEND
      IF (X.NE.XEND) GO TO 50
      IND = 3
      C(20) = XEND
      C(21) = ONE
      GO TO 9005
  275 CONTINUE
C                                  STEP NOT ACCEPTED (IND .EQ. 6), SO
C                                    ADD 1 TO THE NO OF SUCCESSIVE
C                                    FAILURES
      C(23) = C(23)+ONE
C                                  ERROR RETURN (WITH IND=-3) IF HMAG
C                                    .LE. HMIN
      IF (C(14).GT.C(13)) GO TO 50
      IND = -3
      GO TO 9005
C                                  END STAGE 4
C
C                                  END LOOP
  280 CONTINUE
C                                  RELATIVE ERROR OPTION SELECTED AND
C                                    Y(K) IS ZERO
      IER = 132
      GO TO 9000
  285 CONTINUE
C                                  X OR XEND WAS NOT CHANGED FROM
C                                    PREVIOUS CALL
      IER = 131
      GO TO 9000
C
  290 CONTINUE
C                                  IND OUT OF RANGE
      IER = 130
      GO TO 9000
  295 CONTINUE
C                                  N.GT.NW OR TOL.LE.0
      IER = 129
C                                  BEGIN ABORT ACTION
 9000 CONTINUE
      CALL UERTST (IER,6HDVERK )
 9005 CONTINUE
      RETURN
      END
