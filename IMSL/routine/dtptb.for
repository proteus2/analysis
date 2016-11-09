C   IMSL ROUTINE NAME   - DTPTB
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1982
C
C   PURPOSE             - SOLVE A SYSTEM OF ORDINARY DIFFERENTIAL
C                           EQUATIONS WITH BOUNDARY CONDITIONS AT TWO
C                           POINTS, USING A MULTIPLE SHOOTING METHOD
C
C   USAGE               - CALL DTPTB (N,FCNI,FCNJ,FCNB,XA,XB,NITER,X,
C                           MAX,Y,IY,DTOL,BTOL,WORK,IER)
C
C   ARGUMENTS    N      - NUMBER OF DIFFERENTIAL EQUATIONS. (INPUT)
C                FCNI   - NAME OF SUBROUTINE FOR EVALUATING DERIVATIVES.
C                           (INPUT)
C                           THE SUBROUTINE ITSELF MUST ALSO BE PROVIDED
C                             BY THE USER AND IT SHOULD BE OF THE
C                             FOLLOWING FORM
C                               SUBROUTINE FCNI(N,X,Y,YPRIME)
C                               REAL Y(N),YPRIME(N)
C                                    .
C                                    .
C                                    .
C                           FCNI SHOULD EVALUATE YPRIME(1)...YPRIME(N)
C                             GIVEN N,X, AND Y(1)...Y(N).  YPRIME(I) IS
C                             THE DERIVATIVE OF Y(I) WITH RESPECT TO X.
C                           FCNI MUST APPEAR IN AN EXTERNAL STATEMENT IN
C                             THE CALLING PROGRAM.
C                FCNJ   - NAME OF SUBROUTINE FOR EVALUATING THE N BY N
C                           JACOBIAN MATRIX OF PARTIAL DERIVATIVES.
C                           (INPUT)
C                           THE SUBROUTINE ITSELF MUST ALSO BE PROVIDED
C                             BY THE USER AND IT SHOULD BE OF THE
C                             FOLLOWING FORM
C                               SUBROUTINE FCNJ(N,X,Y,PD)
C                               REAL Y(N),PD(N,N)
C                                    .
C                                    .
C                                    .
C                           FCNJ SHOULD EVALUATE PD(I,J) FOR I,J=1,N
C                             GIVEN N,X, AND Y(1)...Y(N). PD(I,J) IS
C                             THE PARTIAL DERIVATIVE OF YPRIME(I) WITH
C                             RESPECT TO Y(J).
C                           FCNJ MUST APPEAR IN AN EXTERNAL STATEMENT IN
C                             THE CALLING PROGRAM.
C                FCNB   - NAME OF THE SUBROUTINE FOR EVALUATING THE
C                           BOUNDARY CONDITIONS. (INPUT)
C                           THE SUBROUTINE ITSELF MUST ALSO BE PROVIDED
C                             BY THE USER AND IT SHOULD BE OF THE
C                             FOLLOWING FORM
C                               SUBROUTINE FCNB(N,YA,YB,F)
C                               REAL YA(N),YB(N),F(N)
C                                    .
C                                    .
C                                    .
C                           FCNB SHOULD EVALUATE F(1)...F(N) GIVEN
C                             YA(1)...YA(N),YB(1)...YB(N).  YA(I) AND
C                             YB(I) ARE THE VALUES OF Y(I) AT XA AND
C                             XB, RESPECTIVELY, AND THE BOUNDARY
C                             CONDITIONS ARE DEFINED BY F(I)=0.0, I=1,N.
C                           FCNB MUST APPEAR IN AN EXTERNAL STATEMENT
C                             IN THE CALLING PROGRAM.
C                XA,XB  - TWO POINTS WHERE BOUNDARY CONDITIONS ARE
C                           GIVEN. (INPUT) XA MUST BE LESS THAN XB.
C                NITER  - MAXIMUM NUMBER OF ITERATIONS OF NEWTONS METHOD
C                           TO BE PERMITTED. (INPUT) ITERATION WILL STOP
C                           IF CONVERGENCE OCCURS EARLY. SUGGESTED
C                           VALUES ARE NITER=2 FOR LINEAR AND NITER=9
C                           FOR NON-LINEAR PROBLEMS.
C                X      - ARRAY OF LENGTH IABS(MAX) OF SHOOTING
C                           POINTS.
C                           IF MAX IS POSITIVE, X IS DEFINED ON OUTPUT
C                             AND NEED NOT BE DEFINED BY THE USER.
C                           IF MAX IS NEGATIVE, X SHOULD CONTAIN
C                             THE SHOOTING POINTS ON INPUT. X SHOULD
C                             BE AN INCREASING SEQUENCE WITH X(1)=XA
C                             AND X(IABS(MAX))=XB.
C                MAX    - IF THE PROGRAM IS TO CHOOSE THE SHOOTING
C                           POINTS (NOTE- PROGRAM SELECTION OF SHOOTING
C                           POINTS MAY BE EXPENSIVE, AND SHOULD BE USED
C                           ONLY FOR LINEAR PROBLEMS)
C                             ON INPUT, MAX SHOULD BE AN UPPER BOUND ON
C                               THE NUMBER OF SHOOTING POINTS.
C                               MAX=10 WILL WORK FOR MANY PROBLEMS.
C                             ON OUTPUT, MAX WILL RETURN THE NUMBER OF
C                               POINTS ACTUALLY USED. THE FIRST MAX
C                               COMPONENTS OF X WILL CONTAIN THESE
C                               POINTS. IF MORE THAN MAX POINTS ARE
C                               NEEDED FOR STABILITY, AN ERROR MESSAGE
C                               WILL BE PRINTED AND MAX WILL RETURN A
C                               ROUGH ESTIMATE OF THE NUMBER NEEDED.
C                         IF THE USER IS TO SUPPLY THE SHOOTING POINTS
C                           (SUGGESTED AFTER FAILURE WITH PROGRAM-
C                           SELECTED POINTS)
C                             ON INPUT, MAX SHOULD BE A NEGATIVE INTEGER
C                               WHOSE ABSOLUTE VALUE GIVES THE NUMBER OF
C                               POINTS SUPPLIED BY THE USER. MAX=-10 IS
C                               A SUGGESTED VALUE FOR THE FIRST ATTEMPT.
C                             ON OUTPUT, MAX IS UNCHANGED.
C                Y      - N BY IABS(MAX) SOLUTION ARRAY. (OUTPUT AND
C                           -IF MAX IS NEGATIVE-INPUT).  Y(J,I) WILL
C                           RETURN AN APPROXIMATION TO THE JTH SOLUTION
C                           AT X(I).
C                         USING THIS OUTPUT ONE OF THE IMSL INITIAL
C                           VALUE ROUTINES CAN BE USED TO CALCULATE THE
C                           SOLUTION AT OTHER POINTS.
C                         IF MAX IS NEGATIVE, Y(J,I) IS EXPECTED TO
C                           CONTAIN, ON INPUT, A STARTING APPROXIMATION
C                           TO THE VALUE OF THE JTH SOLUTION AT X(I),
C                           FOR USE IN THE NEWTON ITERATION.
C                IY     - FIRST DIMENSION OF Y EXACTLY AS SPECIFIED IN
C                           THE DIMENSION STATEMENT IN THE CALLING
C                           PROGRAM. (INPUT)
C                DTOL   - DIFFERENTIAL EQUATION ERROR TOLERANCE. (INPUT)
C                           DTPTB ATTEMPTS TO CONTROL THE LOCAL ERROR IN
C                           SUCH A WAY THAT THE GLOBAL ERROR IS
C                           PROPORTIONAL TO DTOL. DTOL IS USED AS THE
C                           INPUT ERROR CONTROL PARAMETER FOR THE IMSL
C                           INITIAL VALUE ROUTINE DTPTD (A MODIFIED
C                           COPY OF IMSL ROUTINE DVERK). SEE THE
C                           DESCRIPTION OF DVERK ARGUMENT TOL (DEFAULT
C                           CASE EXPLANATION) FOR MORE DETAIL.
C                BTOL   - BOUNDARY CONDITION ERROR TOLERANCE. (INPUT)
C                           THE SOLUTION RETURNED WILL SATISFY ALL
C                           BOUNDARY CONDITIONS F(I)=0.0 TO WITHIN BTOL
C                           TOLERANCE.
C                WORK   - WORK ARRAY OF LENGTH
C                           N*(N+1)*(IABS(MAX)+12)+N+30
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         WARNING ERROR
C                           IER = 33, CAUTION, CONVERGENCE HAS OCCURRED,
C                             BUT TO GET ACCURATE APPROXIMATIONS TO Y AT
C                             A POINT X, IT MAY BE NECESSARY TO
C                             START THE INITIAL VALUE SOLVER (E.G.
C                             DVERK) AT THE NEAREST X(I) WITH X(I).LE.X
C                             (USING Y(J,I), J=1..N AS INITIAL VALUES)
C                             RATHER THAN START AT XA.
C                         TERMINAL ERROR
C                           IER = 129, INITIAL VALUE INTEGRATOR FAILED.
C                             RELAX ERROR TOLERANCE OR SEE REMARK 1.
C                           IER = 130, NUMBER OF SHOOTING POINTS NEEDED
C                             FOR STABILITY EXCEEDS MAX. ON OUTPUT
C                             MAX WILL GIVE A ROUGH ESTIMATE OF
C                             THE NUMBER NEEDED. IF THIS ESTIMATE IS
C                             GREATER THAN 100, THE USER MUST SUPPLY
C                             THE SHOOTING POINTS AS INPUT.
C                           IER = 131, NEWTONS METHOD FAILS TO CONVERGE
C                             IN NITER ITERATIONS.
C                             FOR LINEAR PROBLEMS, DO AN EXTRA
C                               ITERATION. IF IER = 131 STILL OCCURS,
C                               CHECK TO SEE THAT SUBROUTINE FCNJ IS
C                               GIVING CORRECT DERIVATIVES. IF THE
C                               DERIVATIVES ARE CORRECT, SEE REMARK 1.
C                             FOR NON-LINEAR PROBLEMS, SEE REMARK 1.
C                           IER = 132, LINEAR EQUATION SOLVER FAILED.
C                             THE PROBLEM MAY NOT HAVE A UNIQUE
C                             SOLUTION. OTHERWISE, SEE REMARK 1.
C
C   REQD. IMSL ROUTINES - SINGLE/DTPTC,DTPTD,DTPTE,LEQT2F,LUDATN,LUELMN,
C                           LUREFN,UERSET,UERTST,UGETIO
C                       - DOUBLE/DTPTC,DTPTD,DTPTE,LEQT2F,LUDATN,LUELMN,
C                           LUREFN,UERSET,UERTST,UGETIO,VXADD,VXMUL,
C                           VXSTO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS  1.  NON-LINEAR (OR DIFFICULT LINEAR) PROBLEMS
C
C                DTPTB SHOULD HANDLE MOST LINEAR PROBLEMS
C                AUTOMATICALLY USING PROGRAM SELECTED SHOOTING POINTS,
C                BUT NON-LINEAR PROBLEMS MAY REQUIRE MORE USER
C                EFFORT. IF IER=129, 131 OR 132 THE FOLLOWING
C                STEPS MAY BE TAKEN TO INCREASE THE PROBABILITY
C                OF CONVERGENCE (SUGGESTED ORDER AS INDICATED)
C                A. INCREASE THE NUMBER OF SHOOTING POINTS
C                   (I.E. INCREASE IABS(MAX)). REMEMBER TO
C                   MAKE THE CORRESPONDING DIMENSION
C                   STATEMENT INCREASES. WITH MANY SHOOTING
C                   POINTS THE PROGRAM ESSENTIALLY USES
C                   A FINITE DIFFERENCE METHOD, WHICH HAS
C                   LESS TROUBLE WITH NON-LINEARITIES THAN
C                   SHOOTING METHODS. AFTER A CERTAIN POINT
C                   HOWEVER, INCREASING THE NUMBER OF POINTS
C                   WILL NO LONGER HELP CONVERGENCE.
C                B. PARAMETERIZE THE PROBLEM WITH A VARIABLE
C                   ALPHA SO THAT WITH ALPHA=0.0 THE PROBLEM
C                   IS EASY TO SOLVE (E.G. LINEAR) AND WITH
C                   ALPHA=1.0 IT REDUCES TO THE ORIGINAL
C                   NON-LINEAR PROBLEM. IF NITER.GE.10, ALPHA,
C                   WHICH WILL BE THE ONLY VARIABLE IN COMMON
C                   BLOCK /DALPHA/, WILL BE SET EQUAL TO
C                     ALPHA=MIN((ITER-1)/(NITER-9),1.0)
C                   WHERE ITER=1,...,NITER IS THE ITERATION
C                   COUNTER. THUS THE FIRST ITERATION WILL
C                   BE DONE WITH ALPHA=0.0. PROGRESSIVELY
C                   MORE DIFFICULT PROBLEMS WILL BE SOLVED
C                   ON THE FOLLOWING ITERATIONS, AND THE
C                   LAST 9 (LESS IF CONVERGENCE OCCURS
C                   EARLIER) ITERATIONS ARE DONE WITH ALPHA=1.0
C                   (THE ORIGINAL PROBLEM). A GOOD
C                   PARAMETERIZATION SHOULD BE SUCH THAT THE
C                   SOLUTION OF THE (EASY) PROBLEM WITH ALPHA=0.0
C                   IS CLOSE TO THE SOLUTION OF THE (DIFFICULT)
C                   PROBLEM WITH ALPHA=1.0.
C            2.  DTPTB IS NOT RECOMMENDED FOR STIFF SYSTEMS OF
C                DIFFERENTIAL EQUATIONS.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE DTPTB  (N,FCNI,FCNJ,FCNB,XA,XB,NITER,X,MAX,Y,IY,DTOL,
     1                   BTOL,WORK,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N,NITER,MAX,IY,IER
      REAL               XA,XB,X(1),Y(IY,1),DTOL,BTOL,WORK(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            NS,NA,NF,NTA,NTB,NTC
      EXTERNAL           FCNI,FCNJ,FCNB
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
C                                  SET UP WORK ARRAY INDICIES
      NS=N+1
      NA=NS+9*N*NS
      NF=NA+N*N*IABS(MAX)
      NTA=NF+N*IABS(MAX)
      NTB=NTA+N*NS
      NTC=NTB+N*NS
C                                  CALL DTPTC TO PERFORM ACTUAL
C                                    CALCULATIONS
      CALL DTPTC(FCNI,FCNJ,FCNB,XA,XB,X,MAX,Y,IY,N,NITER,WORK(1),
     1 WORK(NS),DTOL,BTOL,WORK(NA),WORK(NF),WORK(NTA),WORK(NTB),
     2 WORK(NTC),IER)
C
 9000 CONTINUE
      IF (IER.NE.0) CALL UERTST (IER,6HDTPTB )
 9005 RETURN
      END
