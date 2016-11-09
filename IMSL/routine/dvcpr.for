C   IMSL ROUTINE NAME   - DVCPR
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1982
C
C   PURPOSE             - SOLVE A SYSTEM OF ORDINARY DIFFERENTIAL
C                           EQUATIONS WITH BOUNDARY CONDITIONS AT TWO
C                           POINTS, USING A VARIABLE ORDER, VARIABLE
C                           STEP SIZE FINITE DIFFERENCE METHOD WITH
C                           DEFERRED CORRECTIONS
C
C   USAGE               - CALL DVCPR(N,FCNI,FCNJ,FCNB,XA,XB,NGMAX,
C                           NGRID,IP,IR,TOL,X,Y,IY,ABT,PAR,
C                           WORK,IWORK,IER)
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
C                             THE INITIAL CONDITIONS MUST BE DEFINED
C                             FIRST, THEN THE COUPLED CONDITIONS, AND
C                             THEN THE FINAL CONDITIONS.
C                           FCNB MUST APPEAR IN AN EXTERNAL STATEMENT
C                             IN THE CALLING PROGRAM.
C                XA,XB  - TWO POINTS WHERE BOUNDARY CONDITIONS ARE
C                           GIVEN. (INPUT) XA MUST BE LESS THAN XB.
C                NGMAX  - MAXIMUM NUMBER OF GRID POINTS TO BE
C                           ALLOWED (INPUT).
C                NGRID  - NUMBER OF POINTS IN THE INPUT GRID
C                           (COUNTING THE ENDPOINTS). NGRID MUST
C                           BE .GT. 3.  ON OUTPUT, NGRID WILL
C                           CONTAIN THE FINAL NUMBER OF GRID POINTS.
C                           (INPUT/OUTPUT)
C                IP     - NUMBER OF INITIAL CONDITIONS (INPUT).
C                           0.LE.IP.LT.N
C                IR     - NUMBER OF COUPLED BOUNDARY CONDITIONS
C                           (INPUT). 0.LT.IP+IR.LE.N.
C                TOL    - RELATIVE ERROR CONTROL PARAMETER (INPUT).
C                           THE COMPUTATIONS STOP WHEN
C                             ABS(ERROR(J,I))/AMAX1(ABS(Y(J,I)),1.0)
C                           IS LESS THAN TOL FOR ALL J=1...N,
C                           I=1...NGRID, WHERE ERROR(J,I) IS THE
C                           ESTIMATED ERROR IN Y(J,I).
C                X      - VECTOR OF LENGTH NGMAX CONTAINING THE
C                           FINAL GRID. IF PAR(4)=0 THE PROGRAM
C                           INITIALIZES X TO A UNIFORM MESH OF NGRID
C                           POINTS.  OTHERWISE THE USER MUST SUPPLY
C                           THE INITIAL GRID ON INPUT. (INPUT/OUTPUT)
C                Y      - MATRIX OF DIMENSION N BY NGMAX CONTAINING
C                           THE COMPUTED SOLUTION ON THE FINAL GRID.
C                           Y(J,I) WILL RETURN AN APPROXIMATION TO THE
C                           JTH SOLUTION COMPONENT AT X(I).  IF PAR(4)=0
C                           THE PROGRAM INITIALIZES Y TO ZERO.
C                           OTHERWISE THE USER MUST SUPPLY INITIAL
C                           VALUES FOR Y. (INPUT/OUTPUT)
C                IY     - ROW DIMENSION OF MATRIX Y, EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT(INPUT).
C                           IY MUST BE .GE. N.
C                ABT    - VECTOR OF LENGTH N CONTAINING, IN ITS JTH
C                           COMPONENT, AN ESTIMATE OF THE MAXIMUM
C                           ABSOLUTE ERROR OVER THE GRID POINTS FOR
C                           THE JTH SOLUTION COMPONENT. (OUTPUT)
C                PAR    - OPTIONS VECTOR OF LENGTH 5 (INPUT).  IF
C                           PAR(1)=0 THE DEFAULT OPTIONS ARE USED AND
C                           THE REMAINING COMPONENTS ARE IGNORED.
C                           IF PAR(1)=1, ALL REMAINING COMPONENTS OF
C                           PAR MUST BE GIVEN A VALUE.
C                           THE DEFAULT VALUE OF PAR(I) IN EACH CASE
C                           IS ZERO.
C                         PAR(2).GT.0, IMPLIES THAT CONTINUATION IS
C                           TO BE DONE FOR THIS HIGHLY NONLINEAR
C                           PROBLEM.  IT IS ASSUMED THAT THE USER HAS
C                           EMBEDDED HIS PROBLEM IN A ONE PARAMETER
C                           FAMILY  DY/DX = YPRIME(X,Y,EPSNU)
C                                   F(Y(A),Y(B),EPSNU)=0
C                           SUCH THAT FOR EPSNU=0 THE PROBLEM IS
C                           SIMPLE (E.G. LINEAR), AND FOR EPSNU=1
C                           THE ORIGINAL PROBLEM IS RECOVERED.  THE
C                           PROGRAM WILL AUTOMATICALLY ATTEMPT TO GO
C                           FROM EPSNU=0 TO EPSNU=1.  PAR(2) IS THE
C                           STARTING STEP IN THE CONTINUATION.  THE
C                           STEP MAY BE VARIED BY DVCPR BUT A LOWER
C                           BOUND ON THE STEPSIZE OF 0.01 IS IMPOSED.
C                           THE FOLLOWING COMMON BLOCK SHOULD APPEAR
C                           IN SUBROUTINES FCNI,FCNJ AND FCNB-
C                             COMMON /C1/ EPSNU,CONT
C                             REAL EPSNU
C                             LOGICAL CONT
C                           IF CONT=.TRUE. VECTORS YPRIME IN SUBROUTINE
C                           FCNI AND F IN SUBROUTINE FCNB SHOULD BE
C                           DEFINED BY
C                             YPRIME(I) = D(YPRIME(I))/D(EPSNU)
C                             F(I) = D(F(I))/D(EPSNU)
C                           AND WHEN CONT=.FALSE., YPRIME AND F
C                           SHOULD HAVE THEIR NORMAL DEFINITIONS.
C                         PAR(3)=1, IMPLIES THAT INTERMEDIATE
C                           OUTPUT IS TO BE PRINTED (FOR DEBUGGING
C                           PURPOSES)
C                         PAR(4)=1, IMPLIES THAT INITIAL VALUES FOR
C                           X AND Y ARE SUPPLIED BY THE USER.
C                         PAR(5)=1, IMPLIES THAT THE DIFFERENTIAL
C                           EQUATIONS AND BOUNDARY CONDITIONS ARE
C                           LINEAR, AND THE ALGORITHM SHOULD TAKE
C                           ADVANTAGE OF THIS FACT.
C                WORK   - REAL WORK VECTOR OF LENGTH
C                           N*(3*N*NGMAX+4*N+1)+NGMAX*(7*N+2)
C                IWORK  - INTEGER WORK VECTOR OF LENGTH
C                           2*N*NGMAX+N+NGMAX
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER = 129 ILLEGAL VALUES FOR N,NGRID,
C                             IP OR IR
C                           IER = 130 MORE THAN NGMAX GRID POINTS ARE
C                             NEEDED TO SOLVE THE PROBLEM.
C                           IER = 131 NEWTON'S ITERATION DIVERGED.
C                           IER = 132 NEWTON'S ITERATION REACHED
C                             ROUNDOFF ERROR LEVEL.  IF
C                             REQUESTED PRECISION IS NOT
C                             ATTAINED, THIS MEANS THAT TOL
C                             IS TOO SMALL.
C
C   REQD. IMSL ROUTINES - DVCPS,DVCPT,DVCPU,DVCPV,DVCPW,DVCPX,DVCPY,
C                           UGETIO,UERTST
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE DVCPR (N,FCNI,FCNJ,FCNB,XA,XB,NGMAX,NGRID,IP,IR,TOL,X,
     *                   Y,IY,ABT,PAR,WORK,IWORK,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N,NGMAX,NGRID,IP,IR,IY,IER,IWORK(1)
      REAL               XA,XB,TOL,X(1),Y(IY,1),ABT(1),PAR(5),WORK(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            IM17,IM2,IM3,IM4,I,JERROR,J,M10,M11,M12,M14,
     *                   M15,M16,M17,M2,M3,M4,M5,M6,M7,M8,M9,MMAX2,MMAX,
     *                   MTNMAX,MT
      EXTERNAL           FCNI,FCNJ,FCNB
C                                  FIRST EXECUTABLE STATEMENT
      MMAX = N
      MTNMAX = N*NGMAX
      MMAX2 = 2*MMAX
      M2 = 1+MMAX
      M3 = M2+MMAX**2
      M4 = M3+MMAX**2
      M5 = M4+NGMAX
      MT = MTNMAX*MMAX
      M6 = M5+MT
      M7 = M6+MT
      M8 = M7+MT
      M9 = M8+MTNMAX
      M10 = M9+MTNMAX
      M11 = M10+MTNMAX
      M12 = M11+NGMAX
      M14 = M12+MTNMAX
      M15 = M14+MTNMAX
      M16 = M15+MMAX2*MMAX
      M17 = M16+MTNMAX
      IM2 = 1+MTNMAX
      IM3 = IM2+MTNMAX
      IM4 = IM3+NGMAX
      IF (PAR(1).EQ.0.0) PAR(4) = 0.0
      IF (PAR(4).EQ.0.0) GO TO 15
      IM17 = M17
      DO 10 I=1,NGRID
         DO 5 J=1,N
            WORK(IM17) = Y(J,I)
            IM17 = IM17+1
    5    CONTINUE
   10 CONTINUE
   15 CONTINUE
      CALL DVCPS(MMAX,N,NGMAX,NGRID,IP,IR,MTNMAX,MMAX2,XA,XB,PAR,TOL,X,
     *WORK(M17),ABT,WORK(1),WORK(M2),WORK(M3),WORK(M4),WORK(M5),WORK(M6)
     *,WORK(M7),WORK(M8),WORK(M9),WORK(M10),WORK(M11),WORK(M12),
     *WORK(M14),WORK(M15),WORK(M16),IWORK(1),IWORK(IM2),IWORK(IM3),
     *IWORK(IM4),FCNI,FCNJ,FCNB,JERROR)
      IM17 = M17
      DO 25 I=1,NGRID
         DO 20 J=1,N
            Y(J,I) = WORK(IM17)
            IM17 = IM17+1
   20    CONTINUE
   25 CONTINUE
      IER = 0
      IF (JERROR.GT.0) IER = JERROR+128
      IF (IER.GT.0) CALL UERTST(IER,6HDVCPR )
      RETURN
      END
