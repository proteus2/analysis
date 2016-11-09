C   IMSL ROUTINE NAME   - ZANLYT
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - NOVEMBER 1, 1984
C
C   PURPOSE             - ZEROS OF AN ANALYTIC COMPLEX FUNCTION
C                           USING THE MULLER METHOD WITH DEFLATION
C
C   USAGE               - CALL ZANLYT (F,EPS,NSIG,KN,NGUESS,N,X,ITMAX,
C                           INFER,IER)
C
C   ARGUMENTS    F      - A COMPLEX FUNCTION SUBPROGRAM, F(Z), WRITTEN
C                           BY THE USER SPECIFYING THE EQUATION WHOSE
C                           ROOTS ARE TO BE FOUND.  F MUST APPEAR IN
C                           AN EXTERNAL STATEMENT IN THE CALLING PRO-
C                           GRAM.
C                EPS    - 1ST STOPPING CRITERION.  LET FP(Z)=F(Z)/P
C                           WHERE P = (Z-Z(1))*(Z-Z(2))*,,,*(Z-Z(K-1))
C                           AND Z(1),...,Z(K-1) ARE PREVIOUSLY FOUND
C                           ROOTS.  IF ((CABS(F(Z)).LE.EPS) .AND.
C                           (CABS(FP(Z)).LE.EPS)), THEN Z IS ACCEPTED
C                           AS A ROOT. (INPUT)
C                NSIG   - 2ND STOPPING CRITERION.  A ROOT IS ACCEPTED
C                           IF TWO SUCCESSIVE APPROXIMATIONS TO A GIVEN
C                           ROOT AGREE IN THE FIRST NSIG DIGITS. (INPUT)
C                             NOTE. IF EITHER OR BOTH OF THE STOPPING
C                             CRITERIA ARE FULFILLED, THE ROOT IS
C                             ACCEPTED.
C                KN     - THE NUMBER OF KNOWN ROOTS WHICH MUST BE STORED
C                           IN X(1),...,X(KN), PRIOR TO ENTRY TO ZANLYT
C                NGUESS - THE NUMBER OF INITIAL GUESSES PROVIDED. THESE
C                           GUESSES MUST BE STORED IN X(KN+1),...,
C                           X(KN+NGUESS).  NGUESS MUST BE SET EQUAL
C                           TO ZERO IF NO GUESSES ARE PROVIDED. (INPUT)
C                N      - THE NUMBER OF NEW ROOTS TO BE FOUND BY
C                           ZANLYT (INPUT)
C                X      - A COMPLEX VECTOR OF LENGTH KN+N.  X(1),...,
C                           X(KN) ON INPUT MUST CONTAIN ANY KNOWN
C                           ROOTS.  X(KN+1),..., X(KN+N) ON INPUT MAY,
C                           ON USER OPTION, CONTAIN INITIAL GUESSES FOR
C                           THE N NEW ROOTS WHICH ARE TO BE COMPUTED.
C                           IF THE USER DOES NOT PROVIDE AN INITIAL
C                           GUESS, ZERO IS USED.
C                           ON OUTPUT, X(KN+1),...,X(KN+N) CONTAIN THE
C                           APPROXIMATE ROOTS FOUND BY ZANLYT.
C                ITMAX  - THE MAXIMUM ALLOWABLE NUMBER OF ITERATIONS
C                           PER ROOT (INPUT)
C                INFER  - AN INTEGER VECTOR OF LENGTH KN+N.  ON
C                           OUTPUT INFER(J) CONTAINS THE NUMBER OF
C                           ITERATIONS USED IN FINDING THE J-TH ROOT
C                           WHEN CONVERGENCE WAS ACHIEVED.  IF
C                           CONVERGENCE WAS NOT OBTAINED IN ITMAX
C                           ITERATIONS, INFER(J) WILL BE GREATER THAN
C                           ITMAX (OUTPUT).
C                IER    - ERROR PARAMETER (OUTPUT)
C                         WARNING ERROR
C                           IER = 33 INDICATES FAILURE TO CONVERGE WITH-
C                             IN ITMAX ITERATIONS FOR AT LEAST ONE OF
C                             THE (N) NEW ROOTS.
C
C   REQD. IMSL ROUTINES - UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS      ZANLYT ALWAYS RETURNS THE LAST APPROXIMATION FOR ROOT J
C                IN X(J). IF THE CONVERGENCE CRITERION IS SATISFIED,
C                THEN INFER(J) IS LESS THAN OR EQUAL TO ITMAX. IF THE
C                CONVERGENCE CRITERION IS NOT SATISIFIED, THEN INFER(J)
C                IS SET TO EITHER ITMAX+1 OR ITMAX+K, WITH K GREATER
C                THAN 1. INFER(J) = ITMAX+1 INDICATES THAT ZANLYT DID
C                NOT OBTAIN CONVERGENCE IN THE ALLOWED NUMBER OF ITER-
C                ATIONS. IN THIS CASE, THE USER MAY WISH TO SET ITMAX
C                TO A LARGER VALUE. INFER(J) = ITMAX+K MEANS THAT CON-
C                VERGENCE WAS OBTAINED (ON ITERATION K) FOR THE DEFLA-
C                TED FUNCTION
C                              FP(Z) = F(Z)/((Z-Z(1)...(Z-Z(J-1)))
C
C                BUT FAILED FOR F(Z). IN THIS CASE, BETTER INITIAL
C                GUESSES MIGHT HELP OR, IT MIGHT BE NECESSARY TO RELAX
C                THE CONVERGENCE CRITERION.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE ZANLYT (F,EPS,NSIG,KN,NGUESS,N,X,ITMAX,INFER,IER)
C
      DIMENSION           X(1),INFER(1)
      REAL                RZERO,RTEN,RHUN,RP01,AX,EPS1,QZ,EPS,TPQ
      COMPLEX             X,D,DD,DEN,FPRT,FRT,H,RT,T1,T2,T3,
     1                    TEM,Z0,Z1,Z2,BI,F,XX,XL,Y0,Y1,Y2,X0,
     2                    ZERO,P1,ONE,FOUR,P5
      DATA                ZERO/(0.0,0.0)/,P1/(0.1,0.0)/,
     1                    ONE/(1.0,0.0)/,FOUR/(4.0,0.0)/,
     2                    P5/(0.5,0.0)/,
     3                    RZERO/0.0/,RTEN/10.0/,RHUN/100.0/,
     4                    AX/0.1/,ICKMAX/3/,RP01/0.01/
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
      IF (N .LT. 1) GO TO 9005
      EPS1 = RTEN**(-NSIG)
      EPS1 = AMIN1(EPS1,RP01)
C                                  SET NUMBER OF ITERATIONS
      KNP1 = KN+1
      KNPN = KN+N
      KNPNG = KN+NGUESS
      DO 5 I=1,KNPN
         INFER(I) = 0
         IF (I .GT. KNPNG) X(I) = ZERO
    5 CONTINUE
      L= KNP1
   10 JK = 0
      ICK = 0
      XL = X(L)
   15 IC = 0
      H = AX
      H = P1*H
      IF (CABS(XL) .GT. AX) H = P1*XL
C                                  FIRST THREE POINTS ARE
C                                    XL+H,  XL-H,  XL
      RT = XL+H
      NN = 20
      GO TO 50
   20 Z0 = FPRT
      Y0 = FRT
      X0 = RT
      RT = XL-H
      NN = 25
      GO TO 50
   25 Z1 = FPRT
      Y1 = FRT
      H = XL-RT
      D = H/(RT-X0)
      RT = XL
      NN = 30
      GO TO 50
   30 Z2 = FPRT
      Y2 = FRT
C                                  BEGIN MAIN ALGORITHM
   35 DD = ONE + D
      T1 = Z0*D*D
      T2 = Z1*DD*DD
      XX = Z2*DD
      T3 = Z2*D
      BI = T1-T2+XX+T3
      DEN = BI*BI-FOUR*(XX*T1-T3*(T2-XX))
C                                  USE DENOMINATOR OF MAXIMUM AMPLITUDE
      T1 = CSQRT(DEN)
      QZ = RHUN*AMAX1(CABS(BI),CABS(T1))
      T2 = BI + T1
      TPQ = CABS(T2)+QZ
      IF (TPQ .EQ. QZ) T2 = ZERO
      T3 = BI - T1
      TPQ = CABS(T3) + QZ
      IF (TPQ .EQ. QZ) T3 = ZERO
      DEN = T2
      QZ = CABS(T3)-CABS(T2)
      IF (QZ .GT. RZERO) DEN = T3
C                                  TEST FOR ZERO DENOMINATOR
      NN = 30
      IF (CABS(DEN) .LE. RZERO) GO TO 65
      D = -XX/DEN
      D = D+D
      H = D*H
      RT = RT + H
C                                  CHECK CONVERGENCE OF THE FIRST KIND
C
      IF (CABS(H) .LE. EPS1*AMAX1(CABS(RT),AX)) GO TO 70
      IF (IC .NE. 0) GO TO 15
      NN = 40
      GO TO 50
   40 QZ = CABS(FPRT)-CABS(Z2)*RTEN
      IF (QZ .GE. RZERO) GO TO 45
      Z0 = Z1
      Z1 = Z2
      Z2 = FPRT
      Y0 = Y1
      Y1 = Y2
      Y2 = FRT
      GO TO 35
C                                  TAKE REMEDIAL ACTION TO INDUCE
C                                    CONVERGENCE
   45 CONTINUE
      D = D*P5
      H = H*P5
      RT = RT-H
   50 JK = JK+1
      IF (JK .GT. ITMAX) GO TO 75
      FRT = F(RT)
      FPRT = FRT
C                                  TEST TO SEE IF FIRST ROOT IS BEING
C                                     DETERMINED
      IF (L .EQ. 1) GO TO 60
C                                  COMPUTE DEFLATED FUNCTION
      LM1 = L-1
      DO 55 I=1,LM1
         TEM = RT - X(I)
         IF (CABS(TEM) .EQ. RZERO) GO TO 65
   55 FPRT = FPRT/TEM
   60 CONTINUE
C                                  CHECK CONVERGENCE OF THE SECOND KIND
C
      IF (CABS(FPRT) .LE. EPS .AND. CABS(FRT) .LE. EPS) GO TO 80
      IF (NN .EQ. 20) GO TO 20
      IF (NN .EQ. 25) GO TO 25
      IF (NN .EQ. 30) GO TO 30
      IF (NN .EQ. 40) GO TO 40
   65 CONTINUE
      IF (IC .NE. 0) GO TO 15
      TEM = RTEN*EPS1
      IF (CABS(RT) .GT. AX) TEM = TEM*RT
      RT = RT+TEM
      D = (H+TEM)*D/H
      H = H+TEM
      GO TO 50
C                                  CHECK SOLUTION
   70 CONTINUE
      IF (IC .NE. 0) GO TO 80
      IC = 1
      Z0 = Y1
      Z1 = Y2
      Z2 = F(RT)
      XL = RT
      ICK = ICK+1
      IF (ICK .LE. ICKMAX) GO TO 35
C                                  WARNING ERROR, ITMAX = MAXIMUM
      JK = ITMAX + JK
   75 IER = 33
C                                  A ROOT HAS BEEN FOUND
   80 X(L) = RT
      INFER(L) = JK
      L = L+1
      IF (L .LE. KNPN) GO TO 10
      IF (IER .NE. 0) CALL UERTST (IER,6HZANLYT)
 9005 RETURN
      END
