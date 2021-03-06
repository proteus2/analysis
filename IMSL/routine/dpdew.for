C   IMSL ROUTINE NAME   - DPDEW
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1982
C
C   PURPOSE             - NUCLEUS CALLED ONLY BY IMSL SUBROUTINE DPDES
C
C   REQD. IMSL ROUTINES - DGRCS,DPDEX,DPDET,DPDEU,LEQT1B,VMULBF,
C                           UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE DPDEW (FCN,BNDRY,Y,YMAX,ERROR,SAVE1,SAVE2,PW,EQUIL,
     *                   IPIV,N0,XX,WORK)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            IPIV(1),N0
      REAL               Y(N0,1),YMAX(1),ERROR(1),SAVE1(1),SAVE2(1),
     *                   PW(1),EQUIL(1),XX(1),WORK(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            N,MF,KFLAG,JSTART,NQUSED,NSTEP,NFE,NJE,NSQ,I,
     *                   METH,MITER,NQ,L,IDOUB,MFOLD,NOLD,IRET,MEO,MIO,
     *                   IWEVAL,MAXDER,LMAX,IREDO,J,NSTEPJ,J1,J2,M,IER,
     *                   NEWQ,NPW,NERROR,NSAVE1,NSAVE2,NEQUIL,NY,MITER1,
     *                   JER,NWK,IDUMMY(2),NJOB(4)
      REAL               TQ(4)
      REAL               T,H,HMIN,HMAX,EPS,UROUND,HUSED,EL(13),OLDL0,
     *                   TOLD,RMAX,RC,CRATE,EPSOLD,HOLD,FN,EDN,E,EUP,
     *                   BND,RH,R1,CON,R,HL0,R0,D,PHL0,PR3,D1,ENQ3,ENQ2,
     *                   PR2,PR1,ENQ1,EPSJ,DUMMY
      EXTERNAL           FCN,BNDRY
      COMMON /DBAND/     NLC,NUC
      INTEGER            NLC,NUC
      COMMON /GEAR/      T,H,HMIN,HMAX,EPS,UROUND,EPSJ,HUSED,EL,OLDL0,
     *                   TOLD,RMAX,RC,CRATE,EPSOLD,HOLD,FN,EDN,E,EUP,
     *                   BND,RH,R1,R,HL0,R0,D,PHL0,PR3,D1,ENQ3,ENQ2,PR2,
     *                   PR1,ENQ1,DUMMY,TQ,N,MF,KFLAG,JSTART,NSQ,NQUSED,
     *                   NSTEP,NFE,NJE,NPW,NERROR,NSAVE1,NSAVE2,NEQUIL,
     *                   NY,I,METH,MITER,NQ,L,IDOUB,MFOLD,NOLD,IRET,MEO,
     *                   MIO,IWEVAL,MAXDER,LMAX,IREDO,J,NSTEPJ,J1,J2,M,
     *                   NEWQ,IDUMMY
C                                  FIRST EXECUTABLE STATEMENT
      KFLAG = 0
      TOLD = T
C                                  THIS ROUTINE PERFORMS ONE STEP OF
C                                    THE INTEGRATION OF AN INITIAL
C                                    VALUE PROBLEM FOR A SYSTEM OF
C                                    ORDINARY DIFFERENTIAL EQUATIONS.
      IF (JSTART.GT.0) GO TO 55
      IF (JSTART.NE.0) GO TO 10
C                                  ON THE FIRST CALL, THE ORDER IS SET
C                                    TO 1 AND THE INITIAL YDOT IS
C                                    CALCULATED. RMAX IS THE MAXIMUM
C                                    RATIO BY WHICH H CAN BE INCREASED
C                                    IN A SINGLE STEP. IT IS INITIALLY
C                                    1.E4 TO COMPENSATE FOR THE SMALL
C                                    INITIAL H, BUT THEN IS NORMALLY
C                                    EQUAL TO 10. IF A FAILURE OCCURS
C                                    (IN CORRECTOR CONVERGENCE OR ERROR
C                                    TEST), RMAX IS SET AT 2 FOR THE
C                                    NEXT INCREASE.
      CALL DPDET(N,T,Y,SAVE1,XX,WORK,FCN,BNDRY,-1)
      DO 5 I=1,N
    5 Y(I,2) = H*SAVE1(I)
      METH = MF/10
      MITER = MF-10*METH
      NQ = 1
      L = 2
      IDOUB = 3
      RMAX = 1.E4
      RC = 0.
      CRATE = 1.
      HOLD = H
      MFOLD = MF
      NSTEP = 0
      NSTEPJ = 0
      NFE = 1
      NJE = 0
      IRET = 3
      GO TO 15
C                                  IF THE CALLER HAS CHANGED METH,
C                                    DGRCS IS CALLED TO SET THE
C                                    COEFFICIENTS OF THE METHOD. IF THE
C                                    CALLER HAS CHANGED N, EPS, OR
C                                    METH, THE CONSTANTS E, EDN, EUP,
C                                    AND BND MUST BE RESET. E IS A
C                                    COMPARISON FOR ERRORS OF THE
C                                    CURRENT ORDER NQ. EUP IS TO TEST
C                                    FOR INCREASING THE ORDER, EDN FOR
C                                    DECREASING THE ORDER. BND IS USED
C                                    TO TEST FOR CONVERGENCE OF THE
C                                    CORRECTOR ITERATES. IF THE CALLER
C                                    HAS CHANGED H, Y MUST BE RESCALED.
C                                    IF H OR METH HAS BEEN CHANGED,
C                                    IDOUB IS RESET TO L + 1 TO PREVENT
C                                    FURTHER CHANGES IN H FOR THAT MANY
C                                    STEPS.
   10 IF (MF.EQ.MFOLD) GO TO 25
      MEO = METH
      MIO = MITER
      METH = MF/10
      MITER = MF-10*METH
      MFOLD = MF
      IF (MITER.NE.MIO) IWEVAL = MITER
      IF (METH.EQ.MEO) GO TO 25
      IDOUB = L+1
      IRET = 1
   15 CALL DGRCS(METH,NQ,EL,TQ,MAXDER)
      LMAX = MAXDER+1
      RC = RC*EL(1)/OLDL0
      OLDL0 = EL(1)
   20 FN = N
      EDN = FN*(TQ(1)*EPS)**2
      E = FN*(TQ(2)*EPS)**2
      EUP = FN*(TQ(3)*EPS)**2
      BND = FN*(TQ(4)*EPS)**2
      EPSOLD = EPS
      NOLD = N
      GO TO (30, 35, 55), IRET
   25 IF ((EPS.EQ.EPSOLD) .AND. (N.EQ.NOLD)) GO TO 30
      IF (N.EQ.NOLD) IWEVAL = MITER
      IRET = 1
      GO TO 20
   30 IF (H.EQ.HOLD) GO TO 55
      RH = H/HOLD
      H = HOLD
      IREDO = 3
      GO TO 40
   35 RH = AMAX1(RH,HMIN/ABS(H))
   40 RH = AMIN1(RH,HMAX/ABS(H),RMAX)
      R1 = 1.
      DO 50 J=2,L
         R1 = R1*RH
         DO 45 I=1,N
   45    Y(I,J) = Y(I,J)*R1
   50 CONTINUE
      H = H*RH
      RC = RC*RH
      IDOUB = L+1
      IF (IREDO.EQ.0) GO TO 310
C                                  THIS SECTION COMPUTES THE PREDICTED
C                                    VALUES BY EFFECTIVELY MULTIPLYING
C                                    THE Y ARRAY BY THE PASCAL TRIANGLE
C                                    MATRIX. RC IS THE RATIO OF NEW TO
C                                    OLD VALUES OF THE COEFFICIENT
C                                    H*EL(1). WHEN RC DIFFERS FROM 1 BY
C                                    MORE THAN 30 PERCENT, OR THE
C                                    CALLER HAS CHANGED MITER, IWEVAL
C                                    IS SET TO MITER TO FORCE THE
C                                    PARTIALS TO BE UPDATED, IF
C                                    PARTIALS ARE USED. IN ANY CASE,
C                                    THE PARTIALS ARE UPDATED AT LEAST
C                                    EVERY 20-TH STEP.
   55 IF (ABS(RC-1.).GT.0.3) IWEVAL = MITER
      IF (NSTEP.GE.NSTEPJ+20) IWEVAL = MITER
      T = T+H
      DO 70 J1=1,NQ
         DO 65 J2=J1,NQ
            J = (NQ+J1)-J2
            DO 60 I=1,N
   60       Y(I,J) = Y(I,J)+Y(I,J+1)
   65    CONTINUE
   70 CONTINUE
C                                  UP TO 3 CORRECTOR ITERATIONS ARE
C                                    TAKEN. A CONVERGENCE TEST IS MADE
C                                    ON THE R.M.S. NORM OF EACH
C                                    CORRECTION, USING BND, WHICH IS
C                                    DEPENDENT ON EPS. THE SUM OF THE
C                                    CORRECTIONS IS ACCUMULATED IN THE
C                                    VECTOR ERROR(I). THE Y ARRAY IS
C                                    NOT ALTERED IN THE CORRECTOR LOOP.
C                                    THE UPDATED Y VECTOR IS STORED
C                                    TEMPORARILY IN SAVE1.
   75 DO 80 I=1,N
   80 ERROR(I) = 0.
      M = 0
      CALL DPDET(N,T,Y,SAVE2,XX,WORK,FCN,BNDRY,-1)
      NFE = NFE+1
      IF (IWEVAL.LE.0) GO TO 100
C                                  IF INDICATED, THE MATRIX P = I -
C                                    H*EL(1)*J IS REEVALUATED BEFORE
C                                    STARTING THE CORRECTOR ITERATION.
C                                    IWEVAL IS SET TO 0 AS AN INDICATOR
C                                    THAT THIS HAS BEEN DONE. IF MITER
C                                    = 1 OR 2, P IS COMPUTED AND
C                                    PROCESSED IN PSET. IF MITER = 3,
C                                    THE MATRIX USED IS P = I -
C                                    H*EL(1)*D, WHERE D IS A DIAGONAL
C                                    MATRIX.
      IWEVAL = 0
      RC = 1.
      NJE = NJE+1
      NSTEPJ = NSTEP
      GO TO (90, 85, 95), MITER
   85 NFE = NFE+N
   90 CON = -H*EL(1)
      MITER1 = MITER
      CALL DPDEX(FCN,BNDRY,Y,N0,CON,MITER1,YMAX,SAVE1,SAVE2,PW,EQUIL,
     *IPIV,XX,WORK,IER)
      IF (IER.NE.0) GO TO 155
      GO TO 115
   95 CONTINUE
      GO TO 135
  100 IF (MITER.NE.0) GO TO (115, 115, 110), MITER
C
C                                  IN THE CASE OF FUNCTIONAL ITERATION,
C                                    UPDATE Y DIRECTLY FROM THE RESULT
C                                    OF THE LAST FCN CALL.
      D = 0.
      DO 105 I=1,N
         R = H*SAVE2(I)-Y(I,2)
         D = D+((R-ERROR(I))/YMAX(I))**2
         SAVE1(I) = Y(I,1)+EL(1)*R
         ERROR(I) = R
  105 CONTINUE
      GO TO 145
C                                  IN THE CASE OF THE CHORD METHOD,
C                                    COMPUTE THE CORRECTOR ERROR, F SUB
C                                    (M), AND SOLVE THE LINEAR SYSTEM
C                                    WITH THAT AS RIGHT-HAND SIDE AND P
C                                    AS COEFFICIENT MATRIX, USING THE
C                                    LU DECOMPOSITION IF MITER = 1 OR
C                                    2. IF MITER = 3, THE COEFFICIENT
C                                    H*EL(1) IN P IS UPDATED.
  110 CONTINUE
      GO TO 135
  115 DO 120 I=1,N
  120 SAVE1(I) = H*SAVE2(I)-(Y(I,2)+ERROR(I))
      IF (NLC.EQ.-1) GO TO 130
      NWK = (NLC+NUC+1)*N0+1
      DO 125 I=1,N
  125 EQUIL(I) = SAVE1(I)
      NJOB(1) = N
      NJOB(2) = NLC
      NJOB(3) = NUC
      NJOB(4) = 1
      CALL VMULBF(WORK,N0,EQUIL,N0,NJOB,SAVE1,N0)
      CALL LEQT1B(PW,N,NLC,NUC,N0,SAVE1,1,N0,2,PW(NWK),JER)
      GO TO 135
  130 CONTINUE
  135 D = 0.
      DO 140 I=1,N
         ERROR(I) = ERROR(I)+SAVE1(I)
         D = D+(SAVE1(I)/YMAX(I))**2
         SAVE1(I) = Y(I,1)+EL(1)*ERROR(I)
  140 CONTINUE
C                                  TEST FOR CONVERGENCE. IF M.GT.0, THE
C                                    SQUARE OF THE CONVERGENCE RATE
C                                    CONSTANT IS ESTIMATED AS CRATE,
C                                    AND THIS IS USED IN THE TEST.
  145 IF (M.NE.0) CRATE = AMAX1(.9*CRATE,D/D1)
      IF ((D*AMIN1(1.,2.*CRATE)).LE.BND) GO TO 180
      D1 = D
      M = M+1
      IF (M.EQ.3) GO TO 150
      CALL DPDET(N,T,SAVE1,SAVE2,XX,WORK,FCN,BNDRY,-1)
      GO TO 100
C                                  THE CORRECTOR ITERATION FAILED TO
C                                    CONVERGE IN 3 TRIES. IF PARTIALS
C                                    ARE INVOLVED BUT ARE NOT UP TO
C                                    DATE, THEY ARE REEVALUATED FOR THE
C                                    NEXT TRY. OTHERWISE THE Y ARRAY IS
C                                    RETRACTED TO ITS VALUES BEFORE
C                                    PREDICTION, AND H IS REDUCED, IF
C                                    POSSIBLE. IF NOT, A NO-CONVERGENCE
C                                    EXIT IS TAKEN.
  150 NFE = NFE+2
      IF (IWEVAL.EQ.-1) GO TO 175
  155 T = TOLD
      RMAX = 2.
      DO 170 J1=1,NQ
         DO 165 J2=J1,NQ
            J = (NQ+J1)-J2
            DO 160 I=1,N
  160       Y(I,J) = Y(I,J)-Y(I,J+1)
  165    CONTINUE
  170 CONTINUE
      IF (ABS(H).LE.HMIN*1.00001) GO TO 305
      RH = .25
      IREDO = 1
      GO TO 35
  175 IWEVAL = MITER
      GO TO 75
C                                  THE CORRECTOR HAS CONVERGED. IWEVAL
C                                    IS SET TO -1 IF PARTIAL
C                                    DERIVATIVES WERE USED, TO SIGNAL
C                                    THAT THEY MAY NEED UPDATING ON
C                                    SUBSEQUENT STEPS. THE ERROR TEST
C                                    IS MADE AND CONTROL PASSES TO
C                                    STATEMENT 190 IF IT FAILS.
  180 IF (MITER.NE.0) IWEVAL = -1
      NFE = NFE+M
      D = 0.
      DO 185 I=1,N
  185 D = D+(ERROR(I)/YMAX(I))**2
      IF (D.GT.E) GO TO 205
C                                  AFTER A SUCCESSFUL STEP, UPDATE THE
C                                    Y ARRAY. CONSIDER CHANGING H IF
C                                    IDOUB = 1. OTHERWISE DECREASE
C                                    IDOUB BY 1. IF IDOUB IS THEN 1 AND
C                                    NQ .LT. MAXDER, THEN ERROR IS
C                                    SAVED FOR USE IN A POSSIBLE ORDER
C                                    INCREASE ON THE NEXT STEP. IF A
C                                    CHANGE IN H IS CONSIDERED, AN
C                                    INCREASE OR DECREASE IN ORDER BY
C                                    ONE IS CONSIDERED ALSO. A CHANGE
C                                    IN H IS MADE ONLY IF IT IS BY A
C                                    FACTOR OF AT LEAST 1.1. IF NOT,
C                                    IDOUB IS SET TO 10 TO PREVENT
C                                    TESTING FOR THAT MANY STEPS.
      KFLAG = 0
      IREDO = 0
      NSTEP = NSTEP+1
      HUSED = H
      NQUSED = NQ
      DO 195 J=1,L
         DO 190 I=1,N
  190    Y(I,J) = Y(I,J)+EL(J)*ERROR(I)
  195 CONTINUE
      IF (IDOUB.EQ.1) GO TO 225
      IDOUB = IDOUB-1
      IF (IDOUB.GT.1) GO TO 315
      IF (L.EQ.LMAX) GO TO 315
      DO 200 I=1,N
  200 Y(I,LMAX) = ERROR(I)
      GO TO 315
C                                  THE ERROR TEST FAILED. KFLAG KEEPS
C                                    TRACK OF MULTIPLE FAILURES.
C                                    RESTORE T AND THE Y ARRAY TO THEIR
C                                    PREVIOUS VALUES, AND PREPARE TO
C                                    TRY THE STEP AGAIN. COMPUTE THE
C                                    OPTIMUM STEP SIZE FOR THIS OR ONE
C                                    LOWER ORDER.
  205 KFLAG = KFLAG-1
      T = TOLD
      DO 220 J1=1,NQ
         DO 215 J2=J1,NQ
            J = (NQ+J1)-J2
            DO 210 I=1,N
  210       Y(I,J) = Y(I,J)-Y(I,J+1)
  215    CONTINUE
  220 CONTINUE
      RMAX = 2.
      IF (ABS(H).LE.HMIN*1.00001) GO TO 295
      IF (KFLAG.LE.-3) GO TO 285
      IREDO = 2
      PR3 = 1.E+20
      GO TO 235
C                                  REGARDLESS OF THE SUCCESS OR FAILURE
C                                    OF THE STEP, FACTORS PR1, PR2, AND
C                                    PR3 ARE COMPUTED, BY WHICH H COULD
C                                    BE DIVIDED AT ORDER NQ - 1, ORDER
C                                    NQ, OR ORDER NQ + 1, RESPECTIVELY.
C                                    IN THE CASE OF FAILURE, PR3 =
C                                    1.E20 TO AVOID AN ORDER INCREASE.
C                                    THE SMALLEST OF THESE IS
C                                    DETERMINED AND THE NEW ORDER
C                                    CHOSEN ACCORDINGLY. IF THE ORDER
C                                    IS TO BE INCREASED, WE COMPUTE ONE
C                                    ADDITIONAL SCALED DERIVATIVE.
  225 PR3 = 1.E+20
      IF (L.EQ.LMAX) GO TO 235
      D1 = 0.
      DO 230 I=1,N
  230 D1 = D1+((ERROR(I)-Y(I,LMAX))/YMAX(I))**2
      ENQ3 = .5/(L+1)
      PR3 = ((D1/EUP)**ENQ3)*1.4+1.4E-6
  235 ENQ2 = .5/L
      PR2 = ((D/E)**ENQ2)*1.2+1.2E-6
      PR1 = 1.E+20
      IF (NQ.EQ.1) GO TO 245
      D = 0.
      DO 240 I=1,N
  240 D = D+(Y(I,L)/YMAX(I))**2
      ENQ1 = .5/NQ
      PR1 = ((D/EDN)**ENQ1)*1.3+1.3E-6
  245 IF (PR2.LE.PR3) GO TO 250
      IF (PR3.LT.PR1) GO TO 260
      GO TO 255
  250 IF (PR2.GT.PR1) GO TO 255
      NEWQ = NQ
      RH = 1./PR2
      GO TO 275
  255 NEWQ = NQ-1
      RH = 1./PR1
      IF (KFLAG.NE.0 .AND. RH.GT.1.) RH = 1.
      GO TO 275
  260 NEWQ = L
      RH = 1./PR3
      IF (RH.LT.1.1) GO TO 270
      DO 265 I=1,N
  265 Y(I,NEWQ+1) = ERROR(I)*EL(L)/L
      GO TO 280
  270 IDOUB = 10
      GO TO 315
  275 IF ((KFLAG.EQ.0) .AND. (RH.LT.1.1)) GO TO 270
C
C                                  IF THERE IS A CHANGE OF ORDER, RESET
C                                    NQ, L, AND THE COEFFICIENTS. IN
C                                    ANY CASE H IS RESET ACCORDING TO
C                                    RH AND THE Y ARRAY IS RESCALED.
C                                    THEN EXIT FROM 285 IF THE STEP WAS
C                                    OK, OR REDO THE STEP OTHERWISE.
      IF (NEWQ.EQ.NQ) GO TO 35
  280 NQ = NEWQ
      L = NQ+1
      IRET = 2
      GO TO 15
C                                  CONTROL REACHES THIS SECTION IF 3 OR
C                                    MORE FAILURES HAVE OCCURED. IT IS
C                                    ASSUMED THAT THE DERIVATIVES THAT
C                                    HAVE ACCUMULATED IN THE Y ARRAY
C                                    HAVE ERRORS OF THE WRONG ORDER.
C                                    HENCE THE FIRST DERIVATIVE IS
C                                    RECOMPUTED, AND THE ORDER IS SET
C                                    TO 1. THEN H IS REDUCED BY A
C                                    FACTOR OF 10, AND THE STEP IS
C                                    RETRIED. AFTER A TOTAL OF 7
C                                    FAILURES, AN EXIT IS TAKEN WITH
C                                    KFLAG = -2.
  285 IF (KFLAG.EQ.-7) GO TO 300
      RH = .1
      RH = AMAX1(HMIN/ABS(H),RH)
      H = H*RH
      CALL DPDET(N,T,Y,SAVE1,XX,WORK,FCN,BNDRY,-1)
      NFE = NFE+1
      DO 290 I=1,N
  290 Y(I,2) = H*SAVE1(I)
      IWEVAL = MITER
      IDOUB = 10
      IF (NQ.EQ.1) GO TO 55
      NQ = 1
      L = 2
      IRET = 3
      GO TO 15
C                                  ALL RETURNS ARE MADE THROUGH THIS
C                                    SECTION. H IS SAVED IN HOLD TO
C                                    ALLOW THE CALLER TO CHANGE H ON
C                                    THE NEXT STEP.
  295 KFLAG = -1
      GO TO 315
  300 KFLAG = -2
      GO TO 315
  305 KFLAG = -3
      GO TO 315
  310 RMAX = 10.
  315 HOLD = H
      JSTART = NQ
      RETURN
      END
