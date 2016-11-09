C   IMSL ROUTINE NAME   - DGEAR
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - NOVEMBER 1, 1984
C
C   PURPOSE             - DIFFERENTIAL EQUATION SOLVER - VARIABLE ORDER
C                           ADAMS PREDICTOR CORRECTOR METHOD OR
C                           GEARS METHOD
C
C   USAGE               - CALL DGEAR (N,FCN,FCNJ,X,H,Y,XEND,TOL,METH,
C                           MITER,INDEX,IWK,WK,IER)
C
C   ARGUMENTS    N      - INPUT NUMBER OF FIRST-ORDER DIFFERENTIAL
C                           EQUATIONS.
C                FCN    - NAME OF SUBROUTINE FOR EVALUATING FUNCTIONS.
C                           (INPUT)
C                           THE SUBROUTINE ITSELF MUST ALSO BE PROVIDED
C                             BY THE USER AND IT SHOULD BE OF THE
C                             FOLLOWING FORM
C                               SUBROUTINE FCN (N,X,Y,YPRIME)
C                               REAL X,Y(N),YPRIME(N)
C                                    .
C                                    .
C                                    .
C                           FCN SHOULD EVALUATE YPRIME(1),...,YPRIME(N)
C                             GIVEN N,X, AND Y(1),...,Y(N). YPRIME(I)
C                             IS THE FIRST DERIVATIVE OF Y(I) WITH
C                             RESPECT TO X.
C                           FCN MUST APPEAR IN AN EXTERNAL STATEMENT IN
C                             THE CALLING PROGRAM AND N,X,Y(1),...,Y(N)
C                             MUST NOT BE ALTERED BY FCN.
C                FCNJ   - NAME OF THE SUBROUTINE FOR COMPUTING THE
C                           JACOBIAN MATRIX OF PARTIAL DERIVATIVES.
C                           (INPUT)
C                           THE SUBROUTINE ITSELF MUST ALSO BE PROVIDED
C                             BY THE USER.
C                           IF MITER=1 IT SHOULD BE OF THE FOLLOWING
C                             FORM
C                               SUBROUTINE FCNJ (N,X,Y,PD)
C                               REAL X,Y(N),PD(N,N)
C                                    .
C                                    .
C                           FCNJ MUST EVALUATE PD(I,J), THE PARTIAL
C                             DERIVATIVE OF YPRIME(I) WITH RESPECT TO
C                             Y(J), FOR I=1,N AND J=1,N.
C                           IF MITER= -1 IT SHOULD BE OF THE FOLLOWING
C                             FORM
C                               SUBROUTINE FCNJ (N,X,Y,PD)
C                               REAL X,Y(N),PD(1)
C                                    .
C                                    .
C                           FCNJ MUST EVALUATE PD IN BAND STORAGE MODE.
C                             THAT IS, PD(N*(J-I+NLC)+I) IS THE PARTIAL
C                             DERIVATIVE OF YPRIME(I) WITH RESPECT TO
C                             Y(J).  NLC IS THE NUMBER OF LOWER
C                             CODIAGONALS FOR THE BAND MATRIX.
C                           FCNJ MUST APPEAR IN AN EXTERNAL STATEMENT IN
C                             THE CALLING PROGRAM AND N,X,Y(1),...,Y(N)
C                             MUST NOT BE ALTERED BY FCNJ.
C                           FCNJ IS USED ONLY IF MITER IS EQUAL TO
C                             1 OR -1. OTHERWISE A DUMMY ROUTINE CAN
C                             BE SUBSTITUTED. SEE REMARK 1.
C                X      - INDEPENDENT VARIABLE. (INPUT AND OUTPUT)
C                           ON INPUT, X SUPPLIES THE INITIAL VALUE
C                             AND IS USED ONLY ON THE FIRST CALL.
C                           ON OUTPUT, X IS REPLACED WITH THE CURRENT
C                             VALUE OF THE INDEPENDENT VARIABLE AT WHICH
C                             INTEGRATION HAS BEEN COMPLETED.
C                H      - INPUT/OUTPUT.
C                           ON INPUT, H CONTAINS THE NEXT STEP SIZE IN
C                             X. H IS USED ONLY ON THE FIRST CALL.
C                           ON OUTPUT, H CONTAINS THE STEP SIZE USED
C                             LAST, WHETHER SUCCESSFULLY OR NOT.
C                Y      - DEPENDENT VARIABLES, VECTOR OF LENGTH N.
C                           (INPUT AND OUTPUT)
C                           ON INPUT, Y(1),...,Y(N) SUPPLY INITIAL
C                             VALUES.
C                           ON OUTPUT, Y(1),...,Y(N) ARE REPLACED WITH
C                             A COMPUTED VALUE AT XEND.
C                XEND   - INPUT VALUE OF X AT WHICH SOLUTION IS DESIRED
C                           NEXT. INTEGRATION WILL NORMALLY GO
C                           BEYOND XEND AND THE ROUTINE WILL INTERPOLATE
C                           TO X = XEND.
C                         NOTE THAT (X-XEND)*H MUST BE LESS THAN
C                           ZERO (X AND H AS SPECIFIED ON INPUT).
C                TOL    - INPUT RELATIVE ERROR BOUND. TOL MUST BE
C                           GREATER THAN ZERO. TOL IS USED ONLY ON THE
C                           FIRST CALL UNLESS INDEX IS EQUAL TO -1.
C                           TOL SHOULD BE AT LEAST AN ORDER OF
C                           MAGNITUDE LARGER THAN THE UNIT ROUNDOFF
C                           BUT GENERALLY NOT LARGER THAN .001.
C                           SINGLE STEP ERROR ESTIMATES DIVIDED BY
C                           YMAX(I) WILL BE KEPT LESS THAN TOL IN
C                           ROOT-MEAN-SQUARE NORM (EUCLIDEAN NORM
C                           DIVIDED BY SQRT(N)). THE VECTOR YMAX OF
C                           WEIGHTS IS COMPUTED INTERNALLY AND STORED
C                           IN WORK VECTOR WK. INITIALLY YMAX(I) IS
C                           THE ABSOLUTE VALUE OF Y(I), WITH A DEFAULT
C                           VALUE OF ONE IF Y(I) IS EQUAL TO ZERO.
C                           THEREAFTER, YMAX(I) IS THE LARGEST VALUE
C                           OF THE ABSOLUTE VALUE OF Y(I) SEEN SO FAR,
C                           OR THE INITIAL VALUE OF YMAX(I) IF THAT IS
C                           LARGER.
C                METH   - INPUT BASIC METHOD INDICATOR.
C                           USED ONLY ON THE FIRST CALL UNLESS INDEX IS
C                           EQUAL TO -1.
C                         METH = 1, IMPLIES THAT THE ADAMS METHOD IS
C                           TO BE USED.
C                         METH = 2, IMPLIES THAT THE STIFF METHODS OF
C                           GEAR, OR THE BACKWARD DIFFERENTIATION
C                           FORMULAE ARE TO BE USED.
C                MITER  - INPUT ITERATION METHOD INDICATOR.
C                           MITER = 0, IMPLIES THAT FUNCTIONAL
C                             ITERATION IS USED. NO PARTIAL
C                             DERIVATIVES ARE NEEDED. A DUMMY FCNJ
C                             CAN BE USED.
C                           MITER = 1, IMPLIES THAT THE CHORD METHOD
C                             IS USED WITH AN ANALYTIC JACOBIAN. FOR
C                             THIS METHOD, THE USER SUPPLIES
C                             SUBROUTINE FCNJ.
C                           MITER = 2, IMPLIES THAT THE CHORD METHOD
C                             IS USED WITH THE JACOBIAN CALCULATED
C                             INTERNALLY BY FINITE DIFFERENCES.
C                             A DUMMY FCNJ CAN BE USED.
C                           MITER = 3, IMPLIES THAT THE CHORD METHOD
C                             IS USED WITH THE JACOBIAN REPLACED BY
C                             A DIAGONAL APPROXIMATION BASED ON A
C                             DIRECTIONAL DERIVATIVE.
C                             A DUMMY FCNJ CAN BE USED.
C                           MITER = -1 OR -2, IMPLIES USE THE SAME
C                             METHOD AS FOR MITER= 1 OR 2, RESPECTIVELY,
C                             BUT USING A BANDED JACOBIAN MATRIX.  IN
C                             THESE TWO CASES BANDWIDTH INFORMATION
C                             MUST BE PASSED TO DGEAR THROUGH THE
C                             COMMON BLOCK
C                                COMMON /DBAND/ NLC,NUC
C                             WHERE NLC=NUMBER OF LOWER CODIAGONALS
C                                   NUC=NUMBER OF UPPER CODIAGONALS
C                INDEX  - INPUT AND OUTPUT PARAMETER USED TO INDICATE
C                           THE TYPE OF CALL TO THE SUBROUTINE.  ON
C                           OUTPUT INDEX IS RESET TO 0 IF INTEGRATION
C                           WAS SUCCESSFUL.  OTHERWISE, THE VALUE OF
C                           INDEX IS UNCHANGED.
C                         ON INPUT, INDEX = 1, IMPLIES THAT THIS IS THE
C                           FIRST CALL FOR THIS PROBLEM.
C                         ON INPUT, INDEX = 0, IMPLIES THAT THIS IS NOT
C                           THE FIRST CALL FOR THIS PROBLEM.
C                         ON INPUT, INDEX = -1, IMPLIES THAT THIS IS NOT
C                           THE FIRST CALL FOR THIS PROBLEM, AND THE
C                           USER HAS RESET TOL.
C                         ON INPUT, INDEX = 2, IMPLIES THAT THIS IS NOT
C                           THE FIRST CALL FOR THIS PROBLEM. INTEGRATION
C                           IS TO CONTINUE AND XEND IS TO BE HIT EXACTLY
C                           (NO INTERPOLATION IS DONE). THIS VALUE OF
C                           INDEX ASSUMES THAT XEND IS BEYOND THE
C                           CURRENT VALUE OF X.
C                         ON INPUT, INDEX = 3, IMPLIES THAT THIS IS NOT
C                           THE FIRST CALL FOR THIS PROBLEM. INTEGRATION
C                           IS TO CONTINUE AND CONTROL IS TO BE RETURNED
C                           TO THE CALLING PROGRAM AFTER ONE STEP. XEND
C                           IS IGNORED.
C                IWK    - INTEGER WORK VECTOR OF LENGTH N. USED ONLY IF
C                           MITER = 1 OR 2
C                WK     - REAL WORK VECTOR OF LENGTH 4*N+NMETH+NMITER.
C                           THE VALUE OF NMETH DEPENDS ON THE VALUE OF
C                             METH.
C                             IF METH IS EQUAL TO 1,
C                               NMETH IS EQUAL TO N*13.
C                             IF METH IS EQUAL TO 2,
C                               NMETH IS EQUAL TO N*6.
C                           THE VALUE OF NMITER DEPENDS ON THE VALUE OF
C                             MITER.
C                             IF MITER IS EQUAL TO 1 OR 2,
C                               NMITER IS EQUAL TO N*(N+1)
C                             IF MITER IS EQUAL TO -1 OR -2,
C                               NMITER IS EQUAL TO (2*NLC+NUC+3)*N
C                                WHERE NLC=NUMBER OF LOWER CODIAGONALS
C                                      NUC=NUMBER OF UPPER CODIAGONALS
C                             IF MITER IS EQUAL TO 3,
C                               NMITER IS EQUAL TO N.
C                             IF MITER IS EQUAL TO 0,
C                               NMITER IS EQUAL TO 1.
C                           WK MUST REMAIN UNCHANGED BETWEEN SUCCESSIVE
C                           CALLS DURING INTEGRATION.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         WARNING ERROR
C                           IER = 33, IMPLIES THAT X+H WILL EQUAL X ON
C                             THE NEXT STEP. THIS CONDITION DOES NOT
C                             FORCE THE ROUTINE TO HALT. HOWEVER, IT
C                             DOES INDICATE ONE OF TWO CONDITIONS.
C                             THE USER MIGHT BE REQUIRING TOO MUCH
C                             ACCURACY VIA THE INPUT PARAMETER TOL.
C                             IN THIS CASE THE USER SHOULD CONSIDER
C                             INCREASING THE VALUE OF TOL. THE OTHER
C                             CONDITION WHICH MIGHT GIVE RISE TO THIS
C                             ERROR MESSAGE IS THAT THE SYSTEM OF
C                             DIFFERENTIAL EQUATIONS BEING SOLVED
C                             IS STIFF (EITHER IN GENERAL OR OVER
C                             THE SUBINTERVAL OF THE PROBLEM BEING
C                             SOLVED AT THE TIME OF THE ERROR). IN
C                             THIS CASE THE USER SHOULD CONSIDER
C                             USING A NONZERO VALUE FOR THE INPUT
C                             PARAMETER MITER.
C                         WARNING WITH FIX ERROR
C                           IER = 66, IMPLIES THAT THE ERROR TEST
C                             FAILED. H WAS REDUCED BY .1 ONE OR MORE
C                             TIMES AND THE STEP WAS TRIED AGAIN
C                             SUCCESSFULLY.
C                           IER = 67, IMPLIES THAT CORRECTOR
C                             CONVERGENCE COULD NOT BE ACHIEVED.
C                             H WAS REDUCED BY .1 ONE OR MORE TIMES AND
C                             THE STEP WAS TRIED AGAIN SUCCESSFULLY.
C                         TERMINAL ERROR
C                           IER = 132, IMPLIES THE INTEGRATION WAS
C                             HALTED AFTER FAILING TO PASS THE ERROR
C                             TEST EVEN AFTER REDUCING H BY A FACTOR
C                             OF 1.0E10 FROM ITS INITIAL VALUE.
C                             SEE REMARKS.
C                           IER = 133, IMPLIES THE INTEGRATION WAS
C                             HALTED AFTER FAILING TO ACHIEVE
C                             CORRECTOR CONVERGENCE EVEN AFTER
C                             REDUCING H BY A FACTOR OF 1.0E10 FROM
C                             ITS INITIAL VALUE. SEE REMARKS.
C                           IER = 134, IMPLIES THAT AFTER SOME INITIAL
C                             SUCCESS, THE INTEGRATION WAS HALTED EITHER
C                             BY REPEATED ERROR TEST FAILURES OR BY
C                             A TEST ON TOL. SEE REMARKS.
C                           IER = 135, IMPLIES THAT ONE OF THE INPUT
C                             PARAMETERS N,X,H,XEND,TOL,METH,MITER, OR
C                             INDEX WAS SPECIFIED INCORRECTLY.
C                           IER = 136, IMPLIES THAT INDEX HAD A VALUE
C                             OF -1 ON INPUT, BUT THE DESIRED CHANGES
C                             OF PARAMETERS WERE NOT IMPLEMENTED
C                             BECAUSE XEND WAS NOT BEYOND X.
C                             INTERPOLATION TO X = XEND WAS PERFORMED.
C                             TO TRY AGAIN, SIMPLY CALL AGAIN WITH
C                             INDEX EQUAL TO -1 AND A NEW VALUE FOR
C                             XEND.
C
C   REQD. IMSL ROUTINES - DGRCS,DGRIN,DGRPS,DGRST,LUDATF,LUELMF,LEQT1B,
C                           UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS  1.  THE EXTERNAL SUBROUTINE FCNJ IS USED ONLY WHEN
C                INPUT PARAMETER MITER IS EQUAL TO 1 OR -1. OTHERWISE,
C                A DUMMY FUNCTION CAN BE USED. THE DUMMY SUBROUTINE
C                SHOULD BE OF THE FOLLOWING FORM
C                  SUBROUTINE FCNJ (N,X,Y,PD)
C                  INTEGER N
C                  REAL Y(N),PD(N,N),X
C                  RETURN
C                  END
C            2.  AFTER THE INITIAL CALL, IF A NORMAL RETURN OCCURRED
C                (IER=0) AND A NORMAL CONTINUATION IS DESIRED, SIMPLY
C                RESET XEND AND CALL DGEAR AGAIN. ALL OTHER
C                PARAMETERS WILL BE READY FOR THE NEXT CALL. A CHANGE
C                OF PARAMETERS WITH INDEX EQUAL TO -1 CAN BE MADE
C                AFTER EITHER A SUCCESSFUL OR AN UNSUCCESSFUL RETURN.
C            3.  THE COMMON BLOCKS /DBAND/ AND /GEAR/ NEED TO BE
C                PRESERVED BETWEEN CALLS TO DGEAR. IF IT IS NECESSARY
C                FOR THE COMMON BLOCKS TO EXIST IN THE CALLING PROGRAM
C                THE FOLLOWING STATEMENTS SHOULD BE INCLUDED
C                  COMMON  /DBAND/ NLC,NUC
C                  COMMON  /GEAR/ DUMMY(48),SDUMMY(4),IDUMMY(38)
C                WHERE DUMMY, SDUMMY, AND IDUMMY ARE VARIABLE NAMES NOT
C                USED ELSEWHERE IN THE CALLING PROGRAM.  (FOR DOUBLE
C                PRECISION DUMMY IS TYPE DOUBLE AND SDUMMY IS TYPE REAL)
C            4.  THE CHOICE OF VALUES FOR METH AND MITER MAY REQUIRE
C                SOME EXPERIMENTATION, AND ALSO SOME CONSIDERATION OF
C                THE NATURE OF THE PROBLEM AND OF STORAGE REQUIREMENTS.
C                THE PRIME CONSIDERATION IS STIFFNESS. IF
C                THE PROBLEM IS NOT STIFF, THE BEST CHOICE IS PROBABLY
C                METH = 1 WITH MITER = 0. IF THE PROBLEM IS STIFF TO A
C                SIGNIFICANT DEGREE, THEN METH SHOULD BE 2 AND MITER
C                SHOULD BE 1,2,-1,-2 OR 3. IF THE USER HAS NO KNOWLEDGE
C                OF THE INHERENT TIME CONSTANTS OF THE PROBLEM, WITH
C                WHICH TO PREDICT ITS STIFFNESS, ONE WAY TO DETERMINE
C                THIS IS TO TRY METH = 1 AND MITER = 0 FIRST, AND LOOK
C                AT THE BEHAVIOR OF THE SOLUTION COMPUTED AND THE STEP
C                SIZES USED. IF THE TYPICAL VALUES OF H ARE MUCH
C                SMALLER THAN THE SOLUTION BEHAVIOR WOULD SEEM TO
C                REQUIRE (THAT IS, MORE THAN 100 STEPS ARE TAKEN OVER
C                AN INTERVAL IN WHICH THE SOLUTIONS CHANGE BY LESS
C                THAN ONE PERCENT), THEN THE PROBLEM IS PROBABLY STIFF
C                AND THE DEGREE OF STIFFNESS CAN BE ESTIMATED FROM THE
C                VALUES OF H USED AND THE SMOOTHNESS OF THE SOLUTION.
C                IF THE DEGREE OF STIFFNESS IS ONLY SLIGHT, IT MAY BE
C                THAT METH=1 IS MORE EFFICIENT THAN METH=2.
C                EXPERIMENTATION WOULD BE REQUIRED TO DETERMINE THIS.
C                REGARDLESS OF METH, THE LEAST EFFECTIVE VALUE OF
C                MITER IS 0, AND THE MOST EFFECTIVE IS 1,-1,2,OR -2.
C                MITER = 3 IS GENERALLY SOMEWHERE IN BETWEEN. SINCE
C                THE STORAGE REQUIREMENTS GO UP IN THE SAME ORDER AS
C                EFFECTIVENESS, TRADE-OFF CONSIDERATIONS ARE
C                NECESSARY. FOR REASONS OF ACCURACY AND SPEED, THE
C                CHOICE OF ABS(MITER)=1 IS GENERALLY PREFERRED TO
C                ABS(MITER)=2, UNLESS THE SYSTEM IS FAIRLY COMPLICATED
C                (AND FCNJ IS THUS NOT FEASIBLE TO CODE). THE
C                ACCURACY OF THE FCNJ CALCULATION CAN BE CHECKED BY
C                COMPARISON OF THE JACOBIAN WITH THAT GENERATED WITH
C                ABS(MITER)=2. IF THE JACOBIAN MATRIX IS SIGNIFICANTLY
C                DIAGONALLY DOMINANT, THEN THE OPTION MITER = 3 IS
C                LIKELY TO BE NEARLY AS EFFECTIVE AS ABS(MITER)=1 OR 2,
C                AND WILL SAVE CONSIDERABLE STORAGE AND RUN TIME.
C                IT IS POSSIBLE, AND POTENTIALLY QUITE DESIRABLE, TO
C                USE DIFFERENT VALUES OF METH AND MITER IN DIFFERENT
C                SUBINTERVALS OF THE PROBLEM. FOR EXAMPLE, IF THE
C                PROBLEM IS NON-STIFF INITIALLY AND STIFF LATER,
C                METH = 1 AND MITER = 0 MIGHT BE SET INITIALLY, AND
C                METH = 2 AND MITER = 1 LATER.
C            5.  THE INITIAL VALUE OF THE STEP SIZE, H, SHOULD BE
C                CHOSEN CONSIDERABLY SMALLER THAN THE AVERAGE VALUE
C                EXPECTED FOR THE PROBLEM, AS THE FIRST-ORDER METHOD
C                WITH WHICH DGEAR BEGINS IS NOT GENERALLY THE MOST
C                EFFICIENT ONE. HOWEVER, FOR THE FIRST STEP, AS FOR
C                EVERY STEP, DGEAR TESTS FOR THE POSSIBILITY THAT
C                THE STEP SIZE WAS TOO LARGE TO PASS THE ERROR TEST
C                (BASED ON TOL), AND IF SO ADJUSTS THE STEP SIZE
C                DOWN AUTOMATICALLY. THIS DOWNWARD ADJUSTMENT, IF
C                ANY, IS NOTED BY IER HAVING THE VALUES 66 OR 67,
C                AND SUBSEQUENT RUNS ON THE SAME OR SIMILAR PROBLEM
C                SHOULD BE STARTED WITH AN APPROPRIATELY SMALLER
C                VALUE OF H.
C            6.  SOME OF THE VALUES OF INTEREST LOCATED IN THE
C                COMMON BLOCK /GEAR/ ARE
C                A. HUSED, THE STEP SIZE H LAST USED SUCCESSFULLY
C                   (DUMMY(8))
C                B. NQUSED, THE ORDER LAST USED SUCCESSFULLY
C                   (IDUMMY(6))
C                C. NSTEP, THE CUMULATIVE NUMBER OF STEPS TAKEN
C                   (IDUMMY(7))
C                D. NFE, THE CUMULATIVE NUMBER OF FCN EVALUATIONS
C                   (IDUMMY(8))
C                E. NJE, THE CUMULATIVE NUMBER OF JACOBIAN
C                   EVALUATIONS, AND HENCE ALSO OF MATRIX LU
C                   DECOMPOSITIONS (IDUMMY(9))
C            7.  THE NORMAL USAGE OF DGEAR MAY BE SUMMARIZED AS FOLLOWS
C                A. SET THE INITIAL VALUES IN Y.
C                B. SET N, X, H, TOL, METH, AND MITER.
C                C. SET XEND TO THE FIRST OUTPUT POINT, AND INDEX TO 1.
C                D. CALL DGEAR
C                E. EXIT IF IER IS GREATER THAN 128.
C                F. OTHERWISE, DO DESIRED OUTPUT OF Y.
C                G. EXIT IF THE PROBLEM IS FINISHED.
C                H. OTHERWISE, RESET XEND TO THE NEXT OUTPUT POINT, AND
C                   RETURN TO STEP D.
C            8.  THE ERROR WHICH IS CONTROLLED BY WAY OF THE PARAMETER
C                TOL IS AN ESTIMATE OF THE LOCAL TRUNCATION ERROR, THAT
C                IS, THE ERROR COMMITTED ON TAKING A SINGLE STEP WITH
C                THE METHOD, STARTING WITH DATA REGARDED AS EXACT. THIS
C                IS TO BE DISTINGUISHED FROM THE GLOBAL TRUNCATION
C                ERROR, WHICH IS THE ERROR IN ANY GIVEN COMPUTED VALUE
C                OF Y(X) AS A RESULT OF THE LOCAL TRUNCATION ERRORS
C                FROM ALL STEPS TAKEN TO OBTAIN Y(X). THE LATTER ERROR
C                ACCUMULATES IN A NON-TRIVIAL WAY FROM THE LOCAL
C                ERRORS, AND IS NEITHER ESTIMATED NOR CONTROLLED BY
C                THE ROUTINE. SINCE IT IS USUALLY THE GLOBAL ERROR THAT
C                A USER WANTS TO HAVE UNDER CONTROL, SOME
C                EXPERIMENTATION MAY BE NECESSARY TO GET THE RIGHT
C                VALUE OF TOL TO ACHIEVE THE USERS NEEDS. IF THE
C                PROBLEM IS MATHEMATICALLY STABLE, AND THE METHOD USED
C                IS APPROPRIATELY STABLE, THEN THE GLOBAL ERROR AT A
C                GIVEN X SHOULD VARY SMOOTHLY WITH TOL IN A MONOTONE
C                INCREASING MANNER.
C            9.  IF THE ROUTINE RETURNS WITH IER VALUES OF 132, 133,
C                OR 134, THE USER SHOULD CHECK TO SEE IF TOO MUCH
C                ACCURACY IS BEING REQUIRED. THE USER MAY WISH TO
C                SET TOL TO A LARGER VALUE AND CONTINUE. ANOTHER
C                POSSIBLE CAUSE OF THESE ERROR CONDITIONS IS AN
C                ERROR IN THE CODING OF THE EXTERNAL FUNCTIONS FCN
C                OR FCNJ. IF NO ERRORS ARE FOUND, IT MAY BE NECESSARY
C                TO MONITOR INTERMEDIATE QUANTITIES GENERATED BY THE
C                ROUTINE. THESE QUANTITIES ARE STORED IN THE WORK VECTOR
C                WK AND INDEXED BY SPECIFIC ELEMENTS IN THE COMMON BLOCK
C                /GEAR/. IF IER IS 132 OR 134, THE COMPONENTS CAUSING
C                THE ERROR TEST FAILURE CAN BE IDENTIFIED FROM LARGE
C                VALUES OF THE QUANTITY
C                  WK(IDUMMY(11)+I)/WK(I), FOR I=1,...,N.
C                ONE CAUSE OF THIS MAY BE A VERY SMALL BUT NONZERO
C                INITIAL VALUE OF ABS(Y(I)).
C                IF IER IS 133, SEVERAL POSSIBILITIES EXIST.
C                IT MAY BE INSTRUCTIVE TO TRY DIFFERENT VALUES OF MITER.
C                ALTERNATIVELY, THE USER MIGHT MONITOR SUCCESSIVE
C                CORRECTOR ITERATES CONTAINED IN WK(IDUMMY(12)+I), FOR
C                I=1,...,N. ANOTHER POSSIBILITY MIGHT BE TO MONITOR
C                THE JACOBIAN MATRIX, IF ONE IS USED, STORED, BY
C                COLUMN, IN WK(IDUMMY(10)+I), FOR I=1,...,N*N IF
C                ABS(MITER) IS EQUAL TO 1 OR 2, OR FOR I=1,...,N IF
C                MITER IS EQUAL TO 3.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE DGEAR  (N,FCN,FCNJ,X,H,Y,XEND,TOL,METH,MITER,INDEX,
     1                   IWK,WK,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N,METH,MITER,INDEX,IWK(1),IER
      REAL               X,H,Y(N),XEND,TOL,WK(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            NERROR,NSAVE1,NSAVE2,NPW,NY,NC,MFC,KFLAG,
     1                   JSTART,NSQ,NQUSED,NSTEP,NFE,NJE,I,N0,NHCUT,KGO,
     2                   JER,KER,NN,NEQUIL,IDUMMY(21),NLC,NUC
      REAL               SDUMMY(4)
      REAL               T,HH,HMIN,HMAX,EPSC,UROUND,EPSJ,HUSED,TOUTP,
     1                   AYI,D,DN,SEPS,DUMMY(39)
      EXTERNAL           FCN,FCNJ
      COMMON /DBAND/     NLC,NUC
      COMMON /GEAR/      T,HH,HMIN,HMAX,EPSC,UROUND,EPSJ,HUSED,DUMMY,
     1                   TOUTP,SDUMMY,NC,MFC,KFLAG,JSTART,NSQ,NQUSED,
     2                   NSTEP,NFE,NJE,NPW,NERROR,NSAVE1,NSAVE2,NEQUIL,
     3                   NY,IDUMMY,N0,NHCUT
      DATA               SEPS/Z3C100000/
C                                  FIRST EXECUTABLE STATEMENT
      IF (MITER.GE.0) NLC = -1
      KER = 0
      JER = 0
      UROUND = SEPS
C                                  COMPUTE WORK VECTOR INDICIES
      NERROR = N
      NSAVE1 = NERROR+N
      NSAVE2 = NSAVE1+N
      NY = NSAVE2+N
      IF (METH.EQ.1) NEQUIL = NY+13*N
      IF (METH.EQ.2) NEQUIL = NY+6*N
      NPW = NEQUIL + N
      IF (MITER.EQ.0.OR.MITER.EQ.3) NPW = NEQUIL
      MFC = 10*METH+IABS(MITER)
C                                  CHECK FOR INCORRECT INPUT PARAMETERS
C
      IF (MITER.LT.-2.OR.MITER.GT.3) GO TO 85
      IF (METH.NE.1.AND.METH.NE.2) GO TO 85
      IF (TOL.LE.0.) GO TO 85
      IF (N.LE.0) GO TO 85
      IF ((X-XEND)*H.GE.0.) GO TO 85
      IF (INDEX.EQ.0) GO TO 10
      IF (INDEX.EQ.2) GO TO 15
      IF (INDEX.EQ.-1) GO TO 20
      IF (INDEX.EQ.3) GO TO 25
      IF (INDEX.NE.1) GO TO 85
C                                  IF INITIAL VALUES OF YMAX OTHER THAN
C                                    THOSE SET BELOW ARE DESIRED, THEY
C                                    SHOULD BE SET HERE. ALL YMAX(I)
C                                    MUST BE POSITIVE. IF VALUES FOR
C                                    HMIN OR HMAX, THE BOUNDS ON
C                                    DABS(HH), OTHER THAN THOSE BELOW
C                                    ARE DESIRED, THEY SHOULD BE SET
C                                    BELOW.
      DO 5 I=1,N
         WK(I) = ABS(Y(I))
         IF (WK(I).EQ.0.) WK(I) = 1.
         WK(NY+I) = Y(I)
    5 CONTINUE
      NC = N
      T = X
      HH = H
      IF ((T+HH).EQ.T) KER = 33
      HMIN = ABS(H)
      HMAX = ABS(X-XEND)*10.
      EPSC = TOL
      JSTART = 0
      N0 = N
      NSQ = N0*N0
      EPSJ = SQRT(UROUND)
      NHCUT = 0
      DUMMY(2) = 1.0
      DUMMY(14) = 1.0
      GO TO 30
C                                  TOUTP IS THE PREVIOUS VALUE OF XEND
C                                    FOR USE IN HMAX.
   10 HMAX = ABS(XEND-TOUTP)*10.
      GO TO 45
C
   15 HMAX = ABS(XEND-TOUTP)*10.
      IF ((T-XEND)*HH.GE.0.) GO TO 95
      GO TO 50
C
   20 IF ((T-XEND)*HH.GE.0.) GO TO 90
      JSTART = -1
      NC = N
      EPSC = TOL
C
   25 IF ((T+HH).EQ.T) KER = 33
C
   30 NN = N0
      CALL DGRST (FCN,FCNJ,WK(NY+1),WK,WK(NERROR+1),WK(NSAVE1+1),
     1 WK(NSAVE2+1),WK(NPW+1),WK(NEQUIL+1),IWK,NN)
C
      KGO = 1-KFLAG
      GO TO (35,55,70,80), KGO
C                                  KFLAG = 0, -1, -2, -3
   35 CONTINUE
C                                  NORMAL RETURN FROM INTEGRATOR. THE
C                                    WEIGHTS YMAX(I) ARE UPDATED. IF
C                                    DIFFERENT VALUES ARE DESIRED, THEY
C                                    SHOULD BE SET HERE. A TEST IS MADE
C                                    FOR TOL BEING TOO SMALL FOR THE
C                                    MACHINE PRECISION. ANY OTHER TESTS
C                                    OR CALCULATIONS THAT ARE REQUIRED
C                                    AFTER EVERY STEP SHOULD BE
C                                    INSERTED HERE. IF INDEX = 3, Y IS
C                                    SET TO THE CURRENT SOLUTION ON
C                                    RETURN. IF INDEX = 2, HH IS
C                                    CONTROLLED TO HIT XEND (WITHIN
C                                    ROUNDOFF ERROR), AND THEN THE
C                                    CURRENT SOLUTION IS PUT IN Y ON
C                                    RETURN. FOR ANY OTHER VALUE OF
C                                    INDEX, CONTROL RETURNS TO THE
C                                    INTEGRATOR UNLESS XEND HAS BEEN
C                                    REACHED. THEN INTERPOLATED VALUES
C                                    OF THE SOLUTION ARE COMPUTED AND
C                                    STORED IN Y ON RETURN.
C                                    IF INTERPOLATION IS NOT
C                                    DESIRED, THE CALL TO DGRIN SHOULD
C                                    BE REMOVED AND CONTROL TRANSFERRED
C                                    TO STATEMENT 95 INSTEAD OF 105.
      D = 0.
      DO 40 I=1,N
         AYI = ABS(WK(NY+I))
         WK(I) = AMAX1(WK(I),AYI)
   40 D = D+(AYI/WK(I))**2
      D = D*(UROUND/TOL)**2
      DN = N
      IF (D.GT.DN) GO TO 75
      IF (INDEX.EQ.3) GO TO 95
      IF (INDEX.EQ.2) GO TO 50
   45 IF ((T-XEND)*HH.LT.0.) GO TO 25
      NN = N0
      CALL DGRIN (XEND,WK(NY+1),NN,Y)
      X = XEND
      GO TO 105
   50 IF (((T+HH)-XEND)*HH.LE.0.) GO TO 25
      IF (ABS(T-XEND).LE.UROUND*AMAX1(10.*ABS(T),HMAX)) GO TO 95
      IF ((T-XEND)*HH.GE.0.) GO TO 95
      HH = (XEND-T)*(1.-4.*UROUND)
      JSTART = -1
      GO TO 25
C                                  ON AN ERROR RETURN FROM INTEGRATOR,
C                                    AN IMMEDIATE RETURN OCCURS IF
C                                    KFLAG = -2, AND RECOVERY ATTEMPTS
C                                    ARE MADE OTHERWISE. TO RECOVER, HH
C                                    AND HMIN ARE REDUCED BY A FACTOR
C                                    OF .1 UP TO 10 TIMES BEFORE GIVING
C                                    UP.
   55 JER = 66
   60 IF (NHCUT.EQ.10) GO TO 65
      NHCUT = NHCUT+1
      HMIN = HMIN*.1
      HH = HH*.1
      JSTART = -1
      GO TO 25
C
   65 IF (JER.EQ.66) JER = 132
      IF (JER.EQ.67) JER = 133
      GO TO 95
C
   70 JER = 134
      GO TO 95
C
   75 JER = 134
      KFLAG = -2
      GO TO 95
C
   80 JER = 67
      GO TO 60
C
   85 JER = 135
      GO TO 110
C
   90 JER = 136
      NN = N0
      CALL DGRIN (XEND,WK(NY+1),NN,Y)
      X = XEND
      GO TO 110
C
   95 X = T
      DO 100 I=1,N
  100 Y(I) = WK(NY+I)
  105 IF (JER.LT.128) INDEX = KFLAG
      TOUTP = X
      IF (KFLAG.EQ.0) H = HUSED
      IF (KFLAG.NE.0) H = HH
  110 IER = MAX0(KER,JER)
 9000 CONTINUE
      IF (KER.NE.0.AND.JER.LT.128) CALL UERTST (KER,6HDGEAR )
      IF (JER.NE.0) CALL UERTST (JER,6HDGEAR )
 9005 RETURN
      END
