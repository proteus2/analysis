C   IMSL ROUTINE NAME   - RLLMV
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1980
C
C   PURPOSE             - PERFORM LINEAR REGRESSION USING THE
C                           MINIMAX CRITERION
C
C   USAGE               - CALL RLLMV (XY,IXY,N,M,IOPT,BETA,
C                           REMAX,ITER,IRANK,WK,IER)
C
C   ARGUMENTS    XY     - INPUT/OUTPUT.  ON INPUT XY IS AN N BY (M+5)
C                           MATRIX CONTAINING IN THE FIRST M COLUMNS
C                           THE VALUES OF THE INDEPENDENT VARIABLES
C                           AND IN THE (M+1)TH COLUMN THE VALUES OF
C                           THE DEPENDENT VARIABLE.  THE LAST FOUR
C                           COLUMNS ARE WORKSPACE.  ON OUTPUT, THE
C                           RESIDUALS ARE STORED IN THE (M+5)TH COLUMN
C                           OF XY.  THESE RESIDUALS MAY NOT BE COMPUTED
C                           ACCURATELY, HOWEVER.  SEE REMARKS.  THE
C                           REST OF XY IS DESTROYED ON OUTPUT.
C                IXY    - INPUT ROW DIMENSION OF XY EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT
C                           IN THE CALLING PROGRAM.
C                N      - INPUT NUMBER OF OBSERVATIONS.
C                M      - INPUT NUMBER OF INDEPENDENT VARIABLES.
C                IOPT   - INPUT OPTION PARAMETER SPECIFYING, IF
C                           IOPT=0, THAT AN INTERCEPT TERM IS TO
C                           BE INCLUDED IN THE MODEL. OTHERWISE
C                           NO INTERCEPT TERM IS INCLUDED.
C                BETA   - OUTPUT VECTOR CONTAINING THE ESTIMATED
C                           COEFFICIENTS.  IF IOPT=0, BETA IS OF
C                           LENGTH M+1 WITH THE INTERCEPT ESTIMATE
C                           IN THE (M+1)TH POSITION. OTHERWISE
C                           BETA IS OF LENGTH M.
C                REMAX  - OUTPUT.  THE MAGNITUDE OF THE LARGEST
C                           RESIDUAL.
C                ITER   - OUTPUT.  THE NUMBER OF ITERATIONS PERFORMED.
C                IRANK  - OUTPUT.  RANK OF THE MATRIX OF INDEPENDENT
C                           VARIABLES (INCLUDING THE CONSTANT COLUMN
C                           FOR THE INTERCEPT IF IOPT=0).
C                WK     - WORK VECTOR OF LENGTH M+2.
C                IER    - ERROR PARAMETER.  (OUTPUT)
C                           WARNING ERROR
C                             IER = 33 INDICATES THE SOLUTION
C                               OBTAINED IS PROBABLY NONUNIQUE.
C                           TERMINAL ERROR
C                             IER = 129 INDICATES CALCULATIONS
C                               TERMINATED PREMATURELY DUE TO
C                               ROUNDING.
C
C   REQD. IMSL ROUTINES - UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS      DUE TO THE METHOD THAT THE ROUTINE USES TO COMPUTE
C                  THE RESIDUALS, THE VALUES OUTPUT IN THE (M+5)TH
C                  COLUMN OF XY MAY NOT BE ACCURATE.  THE USER
C                  MAY NEED TO COMPUTE THE RESIDUALS USING THE
C                  ORIGINAL DATA AND THE COEFFICIENTS IN BETA.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE RLLMV  (XY,IXY,N,M,IOPT,BETA,REMAX,ITER,IRANK,WK,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            IXY,N,M,IOPT,ITER,IRANK,IER
      REAL               XY(IXY,1),BETA(1),REMAX,WK(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            IPCOL,IPROW,IRANK1,IROW,I,J,KP1,K,LEV,MM,MODE,
     1                   MP1MK,MP1MR,MP1,MP2,MP3,MP5,NM1,NP1MR
      REAL               BIG,DD,D,FMM,ONE,PIVOT,RELERR,RELTMP,TOL,
     1                   TPIVOT,TWO,VAL,ZERO
      DATA               TOL/Z3C100000/
      DATA               BIG/Z7FFFFFFF/
      DATA               ZERO/0.0/,ONE/1.0/,TWO/2.0/,RELERR/0.0/
C                                  FIRST EXECUTABLE STATEMENT
      MM = M
      MP1 = M+1
      MP2 = M+2
      MP3 = M+3
      MP5 = M+5
      MP1MR = 1
      IRANK = MM
      RELTMP = RELERR
      DO 5 I=1,N
         XY(I,MP5) = XY(I,MP1)
    5 CONTINUE
      IF (IOPT.NE.0) GO TO 15
      DO 10 I=1,N
         XY(I,MP1) = ONE
   10 CONTINUE
      MM = M+1
      MP1 = MM+1
      MP2 = MM+2
      MP3 = MM+3
      IRANK = MM
   15 FMM = MM
      DO 20 J=1,N
         XY(J,MP1) = ONE
         XY(J,MP2) = -XY(J,MP5)
         XY(J,MP3) = MM+J
   20 CONTINUE
      WK(MP1) = ZERO
      ITER = 0
      IER = 0
      DO 25 I=1,MM
         BETA(I) = ZERO
         WK(I) = I
   25 CONTINUE
C                                  LEVEL 1
      LEV = 1
      K = 0
   30 K = K+1
      KP1 = K+1
      MP1MK = MP1-K
      MODE = 0
      DO 35 J=K,N
         XY(J,MP5) = ONE
   35 CONTINUE
C                                  DETERMINE THE VECTOR TO ENTER THE
C                                    BASIS
   40 D = -BIG
      DO 45 J=K,N
         IF (XY(J,MP5).EQ.ZERO) GO TO 45
         DD = ABS(XY(J,MP2))
         IF (DD.LE.D) GO TO 45
         IPCOL = J
         D = DD
   45 CONTINUE
      IF (K.GT.1) GO TO 50
C                                  TEST FOR ZERO RIGHT-HAND SIDE
      IF (D.GT.TOL) GO TO 50
      REMAX = ZERO
      MODE = 2
      GO TO 205
C                                  DETERMINE THE VECTOR TO LEAVE THE
C                                    BASIS
   50 D = TOL
      DO 55 I=1,MP1MK
         DD = ABS(XY(IPCOL,I))
         IF (DD.LE.D) GO TO 55
         IPROW = I
         D = DD
   55 CONTINUE
      IF (D.GT.TOL) GO TO 180
C                                  CHECK FOR LINEAR DEPENDENCE IN
C                                    LEVEL 1
      XY(IPCOL,MP5) = ZERO
      IF (MODE.EQ.1) GO TO 40
      DO 65 J=K,N
         IF (XY(J,MP5).EQ.ZERO) GO TO 65
         DO 60 I=1,MP1MK
            IF (ABS(XY(J,I)).LE.TOL) GO TO 60
            MODE = 1
            GO TO 40
   60    CONTINUE
   65 CONTINUE
      IRANK = K-1
      MP1MR = MP1-IRANK
      IER = 33
      GO TO 95
   70 IF (IPCOL.EQ.K) GO TO 80
C                                  INTERCHANGE COLUMNS IN LEVEL 1
      DO 75 I=1,MP3
         D = XY(IPCOL,I)
         XY(IPCOL,I) = XY(K,I)
         XY(K,I) = D
   75 CONTINUE
   80 IF (IPROW.EQ.MP1MK) GO TO 90
C                                  INTERCHANGE ROWS IN LEVEL 1
      DO 85 J=1,N
         D = XY(J,IPROW)
C
         XY(J,IPROW) = XY(J,MP1MK)
         XY(J,MP1MK) = D
   85 CONTINUE
      D = WK(IPROW)
      WK(IPROW) = WK(MP1MK)
      WK(MP1MK) = D
   90 IF (K.LT.MM) GO TO 30
   95 IF (IRANK.EQ.N) GO TO 205
      IRANK1 = IRANK+1
C                                  LEVEL 2
      LEV = 2
C                                  DETERMINE THE VECTOR TO ENTER THE
C                                    BASIS
      D = TOL
      DO 100 J=IRANK1,N
         DD = ABS(XY(J,MP2))
         IF (DD.LE.D) GO TO 100
         IPCOL = J
         D = DD
  100 CONTINUE
C                                  COMPARE CHEBYSHEV ERROR WITH TOL
      IF (D.GT.TOL) GO TO 105
      REMAX = ZERO
      MODE = 3
      GO TO 205
  105 IF (XY(IPCOL,MP2).LT.-TOL) GO TO 115
      XY(IPCOL,MP1) = TWO-XY(IPCOL,MP1)
      DO 110 I=MP1MR,MP3
         IF (I.EQ.MP1) GO TO 110
         XY(IPCOL,I) = -XY(IPCOL,I)
  110 CONTINUE
C                                  ARRANGE FOR ALL ENTRIES IN PIVOT
C                                    COLUMN (EXCEPT PIVOT) TO BE
C                                    NEGATIVE
  115 DO 125 I=MP1MR,MM
         IF (XY(IPCOL,I).LT.TOL) GO TO 125
         DO 120 J=1,N
            XY(J,MP1) = XY(J,MP1)+TWO*XY(J,I)
            XY(J,I) = -XY(J,I)
  120    CONTINUE
         WK(I) = -WK(I)
  125 CONTINUE
      IPROW = MP1
      GO TO 180
  130 IF (IRANK1.EQ.N) GO TO 205
      IF (IPCOL.EQ.N) GO TO 140
C                                  INTERCHANGE COLUMNS IN LEVEL 2
      DO 135 I=MP1MR,MP3
         D = XY(IPCOL,I)
         XY(IPCOL,I) = XY(N,I)
         XY(N,I) = D
  135 CONTINUE
  140 NM1 = N-1
C                                  LEVEL 3
      LEV = 3
C                                  DETERMINE THE VECTOR TO ENTER THE
C                                    BASIS
  145 D = -TOL
      VAL = TWO*XY(N,MP2)
      DO 155 J=IRANK1,NM1
         IF (XY(J,MP2).GE.D) GO TO 150
         IPCOL = J
         D = XY(J,MP2)
         MODE = 0
         GO TO 155
  150    DD = VAL-XY(J,MP2)
         IF (DD.GE.D) GO TO 155
         MODE = 1
         IPCOL = J
         D = DD
  155 CONTINUE
      IF (D.GE.-TOL) GO TO 205
      DD = -D/XY(N,MP2)
      IF (DD.GE.RELTMP) GO TO 160
      RELERR = DD
      MODE = 4
      GO TO 205
  160 IF (MODE.EQ.0) GO TO 170
      DO 165 I=MP1MR,MP1
         XY(IPCOL,I) = TWO*XY(N,I)-XY(IPCOL,I)
  165 CONTINUE
      XY(IPCOL,MP2) = D
      XY(IPCOL,MP3) = -XY(IPCOL,MP3)
C                                  DETERMINE THE VECTOR TO LEAVE THE
C                                    BASIS
  170 D = BIG
      DO 175 I=MP1MR,MP1
         IF (XY(IPCOL,I).LE.TOL) GO TO 175
         DD = XY(N,I)/XY(IPCOL,I)
         IF (DD.GE.D) GO TO 175
         IPROW = I
         D = DD
  175 CONTINUE
      IF (D.LT.BIG) GO TO 180
      IER = 129
      GO TO 9000
C                                  PIVOT ON XY(IPCOL,IPROW)
  180 PIVOT = XY(IPCOL,IPROW)
      DO 185 J=1,N
         XY(J,IPROW) = XY(J,IPROW)/PIVOT
  185 CONTINUE
      DO 195 J=1,N
         IF (J.EQ.IPCOL) GO TO 195
         D = XY(J,IPROW)
         DO 190 I=MP1MR,MP2
            IF (I.EQ.IPROW) GO TO 190
            XY(J,I) = XY(J,I)-D*XY(IPCOL,I)
  190    CONTINUE
  195 CONTINUE
      TPIVOT = -PIVOT
      DO 200 I=MP1MR,MP2
         XY(IPCOL,I) = XY(IPCOL,I)/TPIVOT
  200 CONTINUE
      XY(IPCOL,IPROW) = ONE/PIVOT
      D = WK(IPROW)
      WK(IPROW) = XY(IPCOL,MP3)
      XY(IPCOL,MP3) = D
      ITER = ITER+1
      GO TO (70,130,145), LEV
C                                  PREPARE OUTPUT
  205 DO 210 J=1,N
         XY(J,MP5) = ZERO
  210 CONTINUE
      IF (MODE.EQ.2) GO TO 240
      DO 215 J=1,IRANK
         K = XY(J,MP3)
         BETA(K) = XY(J,MP2)
  215 CONTINUE
      IF (MODE.EQ.3.OR.IRANK.EQ.N) GO TO 240
      DO 220 I=MP1MR,MP1
         K = ABS(WK(I))-FMM
         XY(K,MP5) = XY(N,MP2)*SIGN(ONE,WK(I))
  220 CONTINUE
      IF (IRANK1.EQ.N) GO TO 230
      DO 225 J=IRANK1,NM1
         K = ABS(XY(J,MP3))-FMM
         XY(K,MP5) = (XY(N,MP2)-XY(J,MP2))*SIGN(ONE,XY(J,MP3))
  225 CONTINUE
C                                  TEST FOR NON-UNIQUE SOLUTION
  230 DO 235 I=MP1MR,MP1
         IF (ABS(XY(N,I)).GT.TOL) GO TO 235
         IER = 33
         GO TO 240
  235 CONTINUE
  240 IF (MODE.NE.2.AND.MODE.NE.3) REMAX = XY(N,MP2)
      IF (IRANK.EQ.N) REMAX = ZERO
      IF (MODE.EQ.4) REMAX = REMAX-D
      IF (IER.EQ.0) GO TO 9005
 9000 CONTINUE
      CALL UERTST (IER,6HRLLMV )
 9005 RETURN
      END
