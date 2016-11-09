C   IMSL ROUTINE NAME   - BEMDP
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1980
C
C   PURPOSE             - MEDIAN POLISH OF A TWO-WAY TABLE
C
C   USAGE               - CALL BEMDP (TAB,IR,IC,IRTAB,MAXIT,WK,IER)
C
C   ARGUMENTS    TAB    - INPUT/OUTPUT MATRIX OF SIZE (IR+1) BY (IC+1).
C                         ON INPUT, TAB CONTAINS THE ENTRIES OF THE
C                           TWO-WAY TABLE IN THE FIRST IR ROWS AND IC
C                           COLUMNS.
C                         ON OUTPUT, TAB CONTAINS THE CELL RESIDUALS
C                           OF THE MEDIAN POLISH IN THE FIRST IR ROWS
C                           AND IC COLUMNS.  THE (IR+1)TH ROW AND THE
C                           (IC+1)TH COLUMN CONTAIN THE MARGINAL
C                           RESIDUALS.
C                IR     - THE NUMBER OF ROWS IN THE TWO-WAY TABLE.
C                           (INPUT)
C                IC     - THE NUMBER OF COLUMNS IN THE TWO-WAY TABLE.
C                           (INPUT)
C                IRTAB  - ROW DIMENSION OF MATRIX TAB EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM. (INPUT)
C                MAXIT  - ON INPUT, MAXIT IS THE MAXIMUM NUMBER OF
C                           ITERATIONS TO BE PERFORMED.
C                         ON OUTPUT, MAXIT IS THE NUMBER OF ITERATIONS
C                           PERFORMED.
C                WK     - WORK VECTOR OF LENGTH MAX(IR,IC).
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         WARNING ERROR
C                           IER=37 INDICATES CONVERGENCE WAS NOT
C                             ACHIEVED WITH MAXIT ITERATIONS.
C                           IER=38 INDICATES IR OR IC IS LESS THAN 2.
C
C   REQD. IMSL ROUTINES - UERTST,UGETIO,VSRTA
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS      THE INPUT TABLE TAB MAY CONTAIN MISSING VALUES.
C                BEMDP CONSIDERS A MISSING VALUE TO BE A VALUE OF
C                -999999.. BEMDP ASSUMES THAT THE INPUT DATA HAS BEEN
C                SCANNED FOR MISSING VALUES AND THAT ALL MISSING
C                VALUES HAVE BEEN REPLACED BY -999999..
C
C-----------------------------------------------------------------------
C
      SUBROUTINE BEMDP  (TAB,IR,IC,IRTAB,MAXIT,WK,IER)
C
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            IR,IC,IRTAB,MAXIT,IER
      REAL               TAB(IRTAB,1),WK(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            ICC,IRR,I,J,NUMIT
      REAL               BIGC,C,EPS,MED,XMISS,ZERO
      DATA               XMISS/-999999./,ZERO/0.0/
      DATA               EPS/Z3C100000/
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
      NUMIT = 0
      IF (IR.LE.1.OR.IC.LE.1) IER = 38
      ICP1 = IC+1
      IRP1 = IR+1
      DO 5 I=1,IR
    5 TAB(I,ICP1) = ZERO
      DO 10 I=1,ICP1
   10 TAB(IRP1,I) = ZERO
   15 NUMIT = NUMIT+1
      BIGC = 0.
C                                  ITERATION FOR ROWS
      DO 30 I=1,IR
         ICC = IC
         DO 20 J=1,IC
            WK(J) = TAB(I,J)
   20    IF (WK(J).EQ.XMISS) ICC = ICC+1
         IF (ICC.GE.(2*IC)) GO TO 30
         CALL VSRTA (WK,IC)
         I1 = (ICC+1)/2
         I2 = I1+ICC-2*I1+1
         MED = (WK(I1)+WK(I2))*.5
         TAB(I,ICP1) = TAB(I,ICP1)+MED
         DO 25 J=1,IC
            C = ABS(MED)
            BIGC = AMAX1(BIGC,C)
            IF (TAB(I,J).EQ.XMISS) GO TO 25
            TAB(I,J) = TAB(I,J)-MED
   25    CONTINUE
   30 CONTINUE
C                                  ITERATION FOR COLUMNS
      DO 45 J=1,ICP1
         IRR = IR
         DO 35 I=1,IR
            WK(I) = TAB(I,J)
   35    IF (WK(I).EQ.XMISS) IRR = IRR+1
         IF (IRR.GE.(2*IR)) GO TO 45
         CALL VSRTA (WK,IR)
         I1 = (IRR+1)/2
         I2 = I1+IRR-2*I1+1
         MED = (WK(I1)+WK(I2))*.5
         TAB(IRP1,J) = TAB(IRP1,J)+MED
         DO 40 I=1,IR
            C = ABS(MED)
            BIGC = AMAX1(BIGC,C)
            IF (TAB(I,J).EQ.XMISS) GO TO 40
            TAB(I,J) = TAB(I,J)-MED
   40    CONTINUE
   45 CONTINUE
C                                  CHECK IMPROVEMENT BETWEEN ITERATIONS
      IF (BIGC.LE.EPS) GO TO 50
C                                  CHECK IF MAXIMUM ITERATION REACHED
      IF (NUMIT.LT.MAXIT) GO TO 15
      IER = 37
   50 MAXIT = NUMIT
 9000 IF (IER.EQ.0) GO TO 9005
      CALL UERTST (IER,'BEMDP ')
 9005 RETURN
      END
