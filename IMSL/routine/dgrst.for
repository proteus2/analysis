C   IMSL ROUTINE NAME   - DGRST
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1982
C
C   PURPOSE             - NUCLEUS CALLED ONLY BY IMSL SUBROUTINE DGEAR
C
C   REQD. IMSL ROUTINES - DGRCS,DGRPS,LUDATF,LUELMF,LEQT1B,UERTST,
C                           UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE DGRST  (FCN,FCNJ,Y,YMAX,ERROR,SAVE1,SAVE2,PW,EQUIL,
     1                   IPIV,N0)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            IPIV(1),N0
      REAL               Y(N0,1),YMAX(1),ERROR(1),SAVE1(1),SAVE2(1),
     1                   PW(1),EQUIL(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            N,MF,KFLAG,JSTART,NQUSED,NSTEP,NFE,NJE,NSQ,
     1                   I,METH,MITER,NQ,L,IDOUB,MFOLD,NOLD,IRET,MEO,
     2                   MIO,IWEVAL,MAXDER,LMAX,IREDO,J,NSTEPJ,J1,J2,
     3                   M,IER,NEWQ,NPW,NERROR,NSAVE1,NSAVE2,NEQUIL,NY,
     4                   MITER1,IDUMMY(2),NLC,NUC,NWK,JER
      REAL               TQ(4)
      REAL               T,H,HMIN,HMAX,EPS,UROUND,HUSED,EL(13),OLDL0,
     1                   TOLD,RMAX,RC,CRATE,EPSOLD,HOLD,FN,EDN,E,EUP,
     2                   BND,RH,R1,CON,R,HL0,R0,D,PHL0,PR3,D1,ENQ3,ENQ2,
     3                   PR2,PR1,ENQ1,EPSJ,DUMMY
      EXTERNAL           FCN,FCNJ
      COMMON /DBAND/     NLC,NUC
      COMMON /GEAR/      T,H,HMIN,HMAX,EPS,UROUND,EPSJ,HUSED,
     1                   EL,OLDL0,TOLD,RMAX,RC,CRATE,EPSOLD,HOLD,FN,
     2                   EDN,E,EUP,BND,RH,R1,R,HL0,R0,D,PHL0,PR3,D1,
     3                   ENQ3,ENQ2,PR2,PR1,ENQ1,DUMMY,TQ,
     4                   N,MF,KFLAG,JSTART,NSQ,NQUSED,NSTEP,NFE,NJE,
     5                   NPW,NERROR,NSAVE1,NSAVE2,NEQUIL,NY,
     6                   I,METH,MITER,NQ,L,IDOUB,MFOLD,NOLD,IRET,MEO,
     7                   MIO,IWEVAL,MAXDER,LMAX,IREDO,J,NSTEPJ,J1,J2,
     8                   M,NEWQ,IDUMMY
C                                  FIRST EXECUTABLE STATEMENT
      KFLAG = 0
      TOLD = T
C                                  THIS ROUTINE PERFORMS ONE STEP OF
C                                    THE INTEGRATION OF AN INITIAL
C                                    VALUE PROBLEM FOR A SYSTEM OF
C                                    ORDINARY DIFFERENTIAL EQUATIONS.
      IF (JSTART.GT.0) GO TO 50
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
      CALL FCN (N,T,Y,SAVE1)
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
   15 CALL DGRCS (METH,NQ,EL,TQ,MAXDER)
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
      GO TO (30,35,50), IRET
   25 IF ((EPS.EQ.EPSOLD).AND.(N.EQ.NOLD)) GO TO 30
      IF (N.EQ.NOLD) IWEVAL = MITER
      IRET = 1
      GO TO 20
   30 IF (H.EQ.HOLD) GO TO 50
      RH = H/HOLD
      H = HOLD
      IREDO = 3
      GO TO 40
   35 RH = AMAX1(RH,HMIN/ABS(H))
   40 RH = AMIN1(RH,HMAX/ABS(H),RMAX)
      R1 = 1.
      DO 45 J=2,L
         R1 = R1*RH
      DO 45 I=1,N
   45 Y(I,J) = Y(I,J)*R1
      H = H*RH
      RC = RC*RH
      IDOUB = L+1
      IF (IREDO.EQ.0) GO TO 285
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
   50 IF (ABS(RC-1.).GT.0.3) IWEVAL = MITER
      IF (NSTEP.GE.NSTEPJ+20) IWEVAL = MITER
      T = T+H
      DO 55 J1=1,NQ
      DO 55 J2=J1,NQ
         J = (NQ+J1)-J2
      DO 55 I=1,N
   55 Y(I,J) = Y(I,J)+Y(I,J+1)
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
   60 DO 65 I=1,N
   65 ERROR(I) = 0.
      M = 0
      CALL FCN (N,T,Y,SAVE2)
      NFE = NFE+1
      IF (IWEVAL.LE.0) GO TO 95
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
      GO TO (75,70,80), MITER
   70 NFE = NFE+N
   75 CON = -H*EL(1)
      MITER1 = MITER
      CALL DGRPS (FCN,FCNJ,Y,N0,CON,MITER1,YMAX,SAVE1,SAVE2,PW,EQUIL,
     1 IPIV,IER)
      IF (IER.NE.0) GO TO 155
      GO TO 125
   80 R = EL(1)*.1
      DO 85 I=1,N
   85 PW(I) = Y(I,1)+R*(H*SAVE2(I)-Y(I,2))
      CALL FCN (N,T,PW,SAVE1)
      NFE = NFE+1
      HL0 = H*EL(1)
      DO 90 I=1,N
         R0 = H*SAVE2(I)-Y(I,2)
         PW(I) = 1.
         D = .1*R0-H*(SAVE1(I)-SAVE2(I))
         SAVE1(I) = 0.
         IF (ABS(R0).LT.UROUND*YMAX(I)) GO TO 90
         IF (ABS(D).EQ.0.) GO TO 155
         PW(I) = .1*R0/D
         SAVE1(I) = PW(I)*R0
   90 CONTINUE
      GO TO 135
   95 IF (MITER.NE.0) GO TO (125,125,105), MITER
C
C                                  IN THE CASE OF FUNCTIONAL ITERATION,
C                                    UPDATE Y DIRECTLY FROM THE RESULT
C                                    OF THE LAST FCN CALL.
      D = 0.
      DO 100 I=1,N
         R = H*SAVE2(I)-Y(I,2)
         D = D+((R-ERROR(I))/YMAX(I))**2
         SAVE1(I) = Y(I,1)+EL(1)*R
  100 ERROR(I) = R
      GO TO 145
C                                  IN THE CASE OF THE CHORD METHOD,
C                                    COMPUTE THE CORRECTOR ERROR, F SUB
C                                    (M), AND SOLVE THE LINEAR SYSTEM
C                                    WITH THAT AS RIGHT-HAND SIDE AND P
C                                    AS COEFFICIENT MATRIX, USING THE
C                                    LU DECOMPOSITION IF MITER = 1 OR
C                                    2. IF MITER = 3, THE COEFFICIENT
C                                    H*EL(1) IN P IS UPDATED.
  105 PHL0 = HL0
      HL0 = H*EL(1)
      IF (HL0.EQ.PHL0) GO TO 115
      R = HL0/PHL0
      DO 110 I=1,N
         D = 1.-R*(1.-1./PW(I))
         IF (ABS(D).EQ.0.) GO TO 165
  110 PW(I) = 1./D
  115 DO 120 I=1,N
  120 SAVE1(I) = PW(I)*(H*SAVE2(I)-(Y(I,2)+ERROR(I)))
      GO TO 135
  125 DO 130 I=1,N
  130 SAVE1(I) = H*SAVE2(I)-(Y(I,2)+ERROR(I))
      IF (NLC .EQ. -1) GO TO 131
      NWK = (NLC+NUC+1)*N0+1
      CALL LEQT1B(PW,N,NLC,NUC,N0,SAVE1,1,N0,2,PW(NWK),JER)
      GO TO 135
  131 CALL LUELMF (PW,SAVE1,IPIV,N,N0,SAVE1)
  135 D = 0.
      DO 140 I=1,N
         ERROR(I) = ERROR(I)+SAVE1(I)
         D = D+(SAVE1(I)/YMAX(I))**2
  140 SAVE1(I) = Y(I,1)+EL(1)*ERROR(I)
C                                  TEST FOR CONVERGENCE. IF M.GT.0, THE
C                                    SQUARE OF THE CONVERGENCE RATE
C                                    CONSTANT IS ESTIMATED AS CRATE,
C                                    AND THIS IS USED IN THE TEST.
  145 IF (M.NE.0) CRATE = AMAX1(.9*CRATE,D/D1)
      IF ((D*AMIN1(1.,2.*CRATE)).LE.BND) GO TO 170
      D1 = D
      M = M+1
      IF (M.EQ.3) GO TO 150
      CALL FCN (N,T,SAVE1,SAVE2)
      GO TO 95
C                                  THE CORRECTOR ITERATION FAILED TO
C                                    CONVERGE IN 3 TRIES. IF PARTIUARE OF THE CONVERGENCE RATE
C                                    CONSTANT IS ESTIMATED AS CRATE,
C                                    AND THIS IS USED IN THE TEST.
  145 IF (M.NE.0) CRATE = AMAX1(.9*CRATE,D/D1)
      IF ((D*AMIN1(1.,2.*CRATE)).LE.BND) GO TO 170
      D1 = D
      M = M+1
      IF (M.EQ.3) GO TO 150
      CALL FCN (N,T,SAVE1,SAVE2)
      GO TO 95
C                                  THE CORRECTOR ITERATION FAILED TO
C                                    CONVERGE IN 3 TRIES. IF PARTIF (IWEVAL.EQ.-1) GO TO 165
  155 T = TOLD
      RMAX = 2.
      DO 160 J1=1,NQ
      DO 160 J2=J1,NQ
         J = (NQ+J1)-J2
      DO 160 I=1,N
  160 Y(I,J) = Y(I,J)-Y(I,J+1)
      IF (ABS(H).LE.HMIN*1.00001) GO TO 280
      RH = .25
      IREDO = 1
      GO TO 35
  165 IWEVAL = MITER
      GO TO 60
C                                  THE CORRECTOR HAS CONVERGED. IWEVAL
C                                    IS SET TO -1 IF PARTIAL
C                                    DERIVATIVES WERE USED, TO SIGNAL
C                                    THAT THEY MAY NEED UPDATING ON
C                                    SUBSEQUENT STEPS. THE ERROR TEST
C                                    IS MADE AND CONTROL PASSES TO
C                                    STATEMENT 190 IF IT FAILS.
  170 IF (MITER.NE.0) IWEVAL = -1
      NFE = NFE+M
      D = 0.
      DO 175 I=1,N
  175 D = D+(ERROR(I)/YMAX(I))**2
      IF (D.GT.E) GO TO 190
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
      DO 180 J=1,L
      DO 180 I=1,N
  180 Y(I,J) = Y(I,J)+EL(J)*ERROR(I)
      IF (IDOUB.EQ.1) GO TO 200
      IDOUB = IDOUB-1
      IF (IDOUB.GT.1) GO TO 290
      IF (L.EQ.LMAX) GO TO 290
      DO 185 I=1,N
  185 Y(I,LMAX) = ERROR(I)
      GO TO 290
C                                  THE ERROR TEST FAILED. KFLAG KEEPS
C                                    TRACK OF MULTIPLE FAILURES.
C                                    RESTORE T AND THE Y ARRAY TO THEIR
C                                    PREVIOUS VALUES, AND PREPARE TO
C                                    TRY THE STEP AGAIN. COMPUTE THE
C                                    OPTIMUM STEP SIZE FOR THIS OR ONE
C                                    LOWER ORDER.
  190 KFLAG = KFLAG-1
      T = TOLD
      DO 195 J1=1,NQ
      DO 195 J2=J1,NQ
         J = (NQ+J1)-J2
      DO 195 I=1,N
  195 Y(I,J) = Y(I,J)-Y(I,J+1)
      RMAX = 2.
      IF (ABS(H).LE.HMIN*1.00001) GO TO 270
      IF (KFLAG.LE.-3) GO TO 260
      IREDO = 2
      PR3 = 1.E+20
      GO TO 210
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
  200 PR3 = 1.E+20
      IF (L.EQ.LMAX) GO TO 210
      D1 = 0.
      DO 205 I=1,N
  205 D1 = D1+((ERROR(I)-Y(I,LMAX))/YMAX(I))**2
      ENQ3 = .5/(L+1)
      PR3 = ((D1/EUP)**ENQ3)*1.4+1.4E-6
  210 ENQ2 = .5/L
      PR2 = ((D/E)**ENQ2)*1.2+1.2E-6
      PR1 = 1.E+20
      IF (NQ.EQ.1) GO TO 220
      D = 0.
      DO 215 I=1,N
  215 D = D+(Y(I,L)/YMAX(I))**2
      ENQ1 = .5/NQ
      PR1 = ((D/EDN)**ENQ1)*1.3+1.3E-6
  220 IF (PR2.LE.PR3) GO TO 225
      IF (PR3.LT.PR1) GO TO 235
      GO TO 230
  225 IF (PR2.GT.PR1) GO TO 230
      NEWQ = NQ
      RH = 1./PR2
      GO TO 250
  230 NEWQ = NQ-1
      RH = 1./PR1
      IF (KFLAG.NE.0.AND.RH.GT.1.) RH = 1.
      GO TO 250
  235 NEWQ = L
      RH = 1./PR3
      IF (RH.LT.1.1) GO TO 245
      DO 240 I=1,N
  240 Y(I,NEWQ+1) = ERROR(I)*EL(L)/L
      GO TO 255
  245 IDOUB = 10
      GO TO 290
  250 IF ((KFLAG.EQ.0).AND.(RH.LT.1.1)) GO TO 245
C
C                                  IF THERE IS A CHANGE OF ORDER, RESET
C                                    NQ, L, AND THE COEFFICIENTS. IN
C                                    ANY CASE H IS RESET ACCORDING TO
C                                    RH AND THE Y ARRAY IS RESCALED.
C                                    THEN EXIT FROM 285 IF THE STEP WAS
C                                    OK, OR REDO THE STEP OTHERWISE.
      IF (NEWQ.EQ.NQ) GO TO 35
  255 NQ = NEWQ
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
  260 IF (KFLAG.EQ.-7) GO TO 275
      RH = .1
      RH = AMAX1(HMIN/ABS(H),RH)
      H = H*RH
      CALL FCN (N,T,Y,SAVE1)
      NFE = NFE+1
      DO 265 I=1,N
  265 Y(I,2) = H*SAVE1(I)
      IWEVAL = MITER
      IDOUB = 10
      IF (NQ.EQ.1) GO TO 50
      NQ = 1
      L = 2
      IRET = 3
      GO TO 15
C                                  ALL RETURNS ARE MADE THROUGH THIS
C                                    SECTION. H IS SAVED IN HOLD TO
C                                    ALLOW THE CALLER TO CHANGE H ON
C                                    THE NEXT STEP.
  270 KFLAG = -1
      GO TO 290
  275 KFLAG = -2
      GO TO 290
  280 KFLAG = -3
      GO TO 290
  285 RMAX = 10.
  290 HOLD = H
      JSTART = NQ
      RETURN
      END
