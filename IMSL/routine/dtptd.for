C   IMSL ROUTINE NAME   - DTPTD
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1980
C
C   PURPOSE             - NUCLEUS CALLED ONLY BY IMSL SUBROUTINE DTPTB
C
C   REQD. IMSL ROUTINES - DTPTE,UERSET,UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE DTPTD  (FCNI,FCNJ,STEMP,N,FCN,X,Y,XEND,TOL,IND,C,
     1                   NW,W,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N,IND,NW,IER
      REAL               X,Y(N),XEND,TOL,C(1),W(NW,9),STEMP(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            K,LEVEL,LEVOLD
      REAL               ZERO,ONE,TWO,THREE,FOUR,FIVE,SEVEN,TEN,HALF,P9
      REAL               C4D15,C2D3,C5D6,C1D6,C1D15,C2D96,TEMP
      REAL               RK(39),REPS,RTOL
      EXTERNAL           FCNI,FCNJ
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
      LEVEL = 0
      CALL UERSET(LEVEL,LEVOLD)
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
      CALL FCN (FCNI,FCNJ,STEMP,N,X,Y,W(1,1))
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
      CALL FCN (FCNI,FCNJ,STEMP,N,X+C(18)*C1D6,W(1,9),W(1,2))
C
      DO 165 K=1,N
         W(K,9) = Y(K)+C(18)*(W(K,1)*RK(2)+W(K,2)*RK(3))
  165 CONTINUE
      CALL FCN (FCNI,FCNJ,STEMP,N,X+C(18)*C4D15,W(1,9),W(1,3))
C
      DO 170 K=1,N
         W(K,9) = Y(K)+C(18)*(W(K,1)*RK(4)-W(K,2)*RK(5)+W(K,3)*RK(6))
  170 CONTINUE
      CALL FCN (FCNI,FCNJ,STEMP,N,X+C(18)*C2D3,W(1,9),W(1,4))
C
      DO 175 K=1,N
         W(K,9) = Y(K)+C(18)*(-W(K,1)*RK(7)+W(K,2)*RK(8)-W(K,3)*RK(9)
     1   +W(K,4)*RK(10))
  175 CONTINUE
      CALL FCN (FCNI,FCNJ,STEMP,N,X+C(18)*C5D6,W(1,9),W(1,5))
C
      DO 180 K=1,N
         W(K,9) = Y(K)+C(18)*(W(K,1)*RK(11)-W(K,2)*RK(12)+W(K,3)*RK(13)
     1   -W(K,4)*RK(14)+W(K,5)*RK(15))
  180 CONTINUE
      CALL FCN (FCNI,FCNJ,STEMP,N,X+C(18),W(1,9),W(1,6))
C
      DO 185 K=1,N
         W(K,9) = Y(K)+C(18)*(-W(K,1)*RK(16)+W(K,2)*RK(17)-W(K,3)
     1   *RK(18)-W(K,4)*RK(19)+W(K,5)*RK(20))
  185 CONTINUE
      CALL FCN (FCNI,FCNJ,STEMP,N,X+C(18)*C1D15,W(1,9),W(1,7))
C
      DO 190 K=1,N
         W(K,9) = Y(K)+C(18)*(W(K,1)*RK(21)-W(K,2)*RK(22)+W(K,3)*RK(23)
     1   -W(K,4)*RK(24)+W(K,5)*RK(25)+W(K,7)*RK(26))
  190 CONTINUE
      CALL FCN (FCNI,FCNJ,STEMP,N,X+C(18),W(1,9),W(1,8))
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
      CALL UERTST (IER,6HDTPTD )
 9005 CONTINUE
      CALL UERSET(LEVOLD,LEVEL)
      RETURN
      END
