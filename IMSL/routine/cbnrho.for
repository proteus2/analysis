C   IMSL ROUTINE NAME   - CBNRHO
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - ESTIMATION OF THE BIVARIATE NORMAL
C                           CORRELATION COEFFICIENT USING A
C                           CONTINGENCY TABLE
C
C   USAGE               - CALL CBNRHO (AP,IA,IB,IRCV,EPS,RHO,VAR,IER)
C
C   ARGUMENTS    AP     - INPUT/OUTPUT ARRAY OF DIMENSION (IR+1,IC+1,3),
C                           WHERE IR = IRCV(1), IC = IRCV(2), AND
C                           CONTAINS IN
C                           AP(I,J,1)    - ON INPUT THE CONTINGENCY
C                             TABLE OF OBSERVED FREQUENCIES ARE AUG-
C                             MENTED BY THE ROW, COLUMN, AND GRAND
C                             TOTALS ON OUTPUT. I=1,2,...,IR AND
C                             J=1,2,...,IC.
C                           AP(I,IC+1,2) - ON OUTPUT THE POINTS OF POLY-
C                             TOMY, IN MONOTONE INCREASING ORDER, FOR
C                             THE VARIABLE CORRESPONDING TO THE ROWS
C                             OF AP. I=1,2,...,IR-1.
C                           AP(IR+1,J,2) - ON OUTPUT THE POINTS OF POLY-
C                             TOMY, IN MONOTONE INCREASING ORDER, FOR
C                             THE VARIABLE CORRESPONDING TO THE COLUMNS
C                             OF AP. J=1,2,...,IC-1.
C                           AP(I,J,2)    - ON OUTPUT THE BIVARIATE
C                             NORMAL PROBABILITIES, CORRESPONDING TO
C                             RHO. I=1,2,...,IR AND J=1,2,...,IC.
C                           AP(I,J,3)    - ON OUTPUT, THE PARTIAL DERIV-
C                             ATIVES OF THE BIVARIATE NORMAL PROBAB-
C                             ILITIES WITH RESPECT TO RHO. I=1,2,...,IR
C                             AND J=1,2,...,IC.
C                IA     - ROW DIMENSION OF MATRIX AP EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM. (INPUT)
C                IB     - COLUMN DIMENSION OF MATRIX AP EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM. (INPUT)
C                IRCV   - INPUT VECTOR OF LENGTH 2. IRCV(I) CONTAINS
C                           WHEN,
C                           I=1,NUMBER OF ROWS IN THE INPUT AP TABLE.
C                           I=2,NUMBER OF COLUMNS IN THE INPUT AP TABLE.
C                EPS    - INPUT SCALAR CONTAINING THE CONVERGENCE
C                           CRITERION FOR RHO. SEE REMARKS.
C                RHO    - OUTPUT SCALAR CONTAINING THE MAXIMUM
C                           LIKELIHOOD ESTIMATE OF THE CORRELATION
C                           COEFFICIENT.
C                VAR    - OUTPUT SCALAR CONTAINING THE ESTIMATE OF
C                           ASYMPTOTIC VARIANCE OF RHO.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         WARNING ERROR (WITH FIX)
C                           IER=65 INDICATES AT LEAST ONE BIVARIATE
C                             NORMAL PROBABILITY IS LESS THAN THE
C                             SMALLEST POSITIVE REAL NUMBER SUCH THAT
C                             1.0+E .GT. 1.0 AND 1.0 .GT. 1.0-E.
C                             THAT PROBABILITY IS SET TO THIS MINIMUM
C                             VALUE.
C                           IER=66 INDICATES THAT THE INTERVAL HAS BEEN
C                             REDUCED AS FAR AS NUMERICALLY POSSIBLE.
C                             THIS IS DUE TO EPS BEING TOO SMALL. WHEN
C                             THIS OCCURS, RHO GIVES THE LOCATION OF
C                             THE MINIMUM AS ACCURATELY AS POSSIBLE.
C                         TERMINAL ERROR
C                           IER=130 INDICATES AT LEAST ONE ROW OR COLUMN
C                             MARGINAL TOTAL OF THE CONTINGENCY TABLE
C                             WAS LESS THAN OR EQUAL TO ZERO.
C                           IER=131 INDICATES THE NUMBER OF ROWS OR
C                             COLUMNS OF THE TABLE WAS SPECIFIED LESS
C                             THAN 2.
C                           IER=132 INDICATES THE INPUT PARAMETER EPS IS
C                             GREATER THAN OR EQUAL TO 2.
C
C   REQD. IMSL ROUTINES - SINGLE/CTBNLL,MDBNOR,MDNOR,MDNRIS,MDTNF,
C                           MERFI,MERRC=ERFC,UERSET,UERTST,UGETIO,ZXGSP
C                       - DOUBLE/CTBNLL,MDBNOR,MDNOR,MDNORD,MDNRIS,
C                           MDTNF,MERFI,MERRC=ERFC,MERRCD=DERFC,
C                           UERSET,UERTST,UGETIO,ZXGSP
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS  1.  THE COMPUTED MAXIMUM LIKELIHOOD ESTIMATE, RHO, WILL BE
C                WITHIN EPS OF THE TRUE MAXIMUM LIKELIHOOD ESTIMATE.
C                FOR EXAMPLE, IF THE TRUE ESTIMATE IS 0.5 AND EPS IS
C                0.01, THEN RHO WILL BE IN THE INTERVAL (0.49,0.51).
C                EPS MUST BE POSITIVE AND LESS THAN 2. EXECUTION TIME
C                INCREASES AS EPS DECREASES.
C            2.  THE INPUT CONTINGENCY TABLE MUST BE ARRANGED SO THAT
C                THE VARIABLE CORRESPONDING TO THE ROWS HAS ITS MAXIMUM
C                VALUE FOR ROW ONE AND DECREASES MONOTONICALLY. THE
C                VARIABLE CORRESPONDING TO TNE COLUMNS MUST HAVE ITS
C                MINIMUM VALUE IN COLUMN ONE AND INCREASE MONOTONICALLY.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE CBNRHO (AP,IA,IB,IRCV,EPS,RHO,VAR,IER)
C
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            IA,IB,IRCV(2),IER
      REAL               AP(IA,IB,1),EPS,RHO,VAR
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,IC,ICM1,ICP1,IR,IRM1,IRP1,I1,J,J1,JER,KER
      REAL               XMIN,ZER(1),APSQD,B,ONE,REPS,S1D2PI,TOTAL
      REAL               T0,T1,T2,T3,T3MAP,T4,X,A,ZERO
      REAL               XTEMP,APTEMP
      DOUBLE PRECISION   SUM
      EXTERNAL           CTBNLL
      DATA               REPS/Z3C100000/
      DATA               S1D2PI/.1591549/,
     2                   ZERO/0.0/,ONE/1.0/
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
      JER = 0
      KER = 0
C                                  SET THE INTERVAL
      A = REPS - ONE
      B = -A
      IF (B - A .GE. EPS) GO TO 5
C                                  TERMINAL - INVALID EPS VALUE
      IER = 132
      GO TO 9000
   5  CONTINUE
      IR = IRCV(1)
      IC = IRCV(2)
      IF (IR .GE. 2 .AND. IC .GE. 2) GO TO 10
C                                  TERMINAL - INVALID NUMBER OF ROWS
C                                             AND/OR COLUMNS
      IER = 131
      GO TO 9000
  10  IRP1 = IR + 1
      ICP1 = IC + 1
      IRM1 = IR - 1
      ICM1 = IC - 1
      TOTAL = ZERO
C                                  INITIALIZE MARGINALS
      DO 15 I = 1,IR
         AP(I,ICP1,1) = ZERO
  15  CONTINUE
      DO 20 I = 1,IC
         AP(IRP1,I,1) = ZERO
  20  CONTINUE
C                                  COMPUTE THE MARGINALS
      AP(IRP1,ICP1,1) = ZERO
      DO 30 I = 1,IR
         DO 25 J = 1,IC
            AP(I,ICP1,1) = AP(I,ICP1,1) + AP(I,J,1)
            AP(IRP1,J,1) = AP(IRP1,J,1) + AP(I,J,1)
  25     CONTINUE
         AP(IRP1,ICP1,1) = AP(IRP1,ICP1,1) + AP(I,ICP1,1)
  30  CONTINUE
      DO 35 I = 1,IR
         IF (AP(I,ICP1,1) .GT. ZERO) GO TO 35
C                                  TERMINAL - A ROW MARGINAL IS ZERO
C                                             OR NEGATIVE
         IER = 130
         GO TO 9000
  35  CONTINUE
      DO 40 J = 1,IC
         IF (AP(IRP1,J,1) .GT. ZERO) GO TO 40
C                                  TERMINAL - A COLUMN MARGINAL IS ZERO
C                                             OR NEGATIVE
         IER = 130
         GO TO 9000
  40  CONTINUE
      TOTAL = ONE/AP(IRP1,ICP1,1)
      X = ZERO
C                                  COMPUTE THE POINTS OF POLYTOMY FOR
C                                  THE ROWS
      DO 45 I = 1,IRM1
         I1 = IRP1 - I
         X = (X + AP(I1,ICP1,1))
         XTEMP = X*TOTAL
         CALL MDNRIS(XTEMP,APTEMP,IER)
         AP(I,ICP1,2) = APTEMP
  45  CONTINUE
      X = ZERO
C                                  COMPUTE THE POINTS OF POLYTOMY FOR
C                                  THE COLUMNS
      DO 50 J = 1,ICM1
         X = (X + AP(IRP1,J,1))
         XTEMP = X*TOTAL
         CALL MDNRIS(XTEMP,APTEMP,IER)
         AP(IRP1,J,2) = APTEMP
  50  CONTINUE
      LEVEL=0
      CALL UERSET(LEVEL,LEVOLD)
C                                  GET RHO
      CALL ZXGSP(CTBNLL,AP,ZER,IRCV,IA,IB,A,B,EPS,XMIN,KER)
      CALL UERSET(LEVOLD,LEVOLD)
C                                  WARNING - CONVERGENCE FAILURE IN
C                                            IMSL ROUTINE CTBNLL
      IF (ZER(1) .NE. 0.0) JER = 65
      IF (KER.EQ.132) IER=66
C                                  NOW GET THE RIGHT RHO
      RHO = XMIN
      T0 = ONE - RHO*RHO
      T1 = S1D2PI/SQRT(T0)
      T2  = -0.5/T0
      T3  = -RHO - RHO
C                                  COMPUTE PHI
C                                             I,J
      DO 70 I = 1,IRM1
         I1 = IRP1 - I
         AP(I1,IC,3) = ZERO
         APSQD = AP(I,ICP1,2)*AP(I,ICP1,2)
         T3MAP = T3*AP(I,ICP1,2)
         DO 65 J = 1,ICM1
            T4 = T2*(APSQD+T3MAP*AP(IRP1,J,2)+AP(IRP1,J,2)*AP(IRP1,J,2))
            AP(I1,J,3) = T1*EXP(T4)
  65     CONTINUE
  70  CONTINUE
      DO 75 J = 1,IC
         AP(1,J,3) = ZERO
  75  CONTINUE
C                                  COMPUTE PARTIAL(PHI   )
C                                                     I,J
      DO 85 J = 1,ICM1
         J1 = ICP1 - J
         DO 80 I = 1,IRM1
            AP(I,J1,3) = AP(I,J1,3) + AP(I+1,J1-1,3) - AP(I,J1-1,3)
     *                  -AP(I+1,J1,3)
  80     CONTINUE
         AP(IR,J1,3) = AP(IR,J1,3) - AP(IR,J1-1,3)
  85  CONTINUE
      DO 90 I = 1,IRM1
         AP(I,1,3) = AP(I,1,3) - AP(I+1,1,3)
  90  CONTINUE
      SUM = 0.D0
C                                  COMPUTE THE VARIANCE
      DO 100 I = 1,IR
         DO 95 J = 1,IC
            SUM = SUM + AP(I,J,3)*AP(I,J,3)/AP(I,J,2)
  95     CONTINUE
 100  CONTINUE
      VAR = TOTAL/SUM
 9000 CONTINUE
      IF (JER .NE. 0) CALL UERTST(JER,'CBNRHO')
      IF (IER .NE. 0) CALL UERTST(IER,'CBNRHO')
      IER = MAX0(IER,JER)
 9005 RETURN
      END
