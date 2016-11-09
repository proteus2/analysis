C   IMSL ROUTINE NAME   - RSMITZ
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JULY 1, 1983
C
C   PURPOSE             - LEAST SQUARES FIT OF THE NON-LINEAR
C                           REGRESSION MODEL
C                           Y(I) = ALPHA+BETA*GAMMA**X(I)+E(I)
C
C   USAGE               - CALL RSMITZ (XY,IXY,IV,PAR,STAT,VCV,ITER,IER)
C
C   ARGUMENTS    XY     - INPUT/OUTPUT IV(1) BY 3 MATRIX.
C                         ON INPUT XY(I,1) CONTAINS THE INDEPENDENT
C                           VARIABLE SETTINGS FOR I=1,2,...,IV(1).
C                           ALL VALUES OF XY(I,1) EITHER MUST BE
C                           NON-NEGATIVE OR MUST BE NON-POSITIVE.
C                         ON INPUT XY(I,2) CONTAINS THE CORRESPONDING
C                           RESPONSE VALUES FOR I=1,2,...,IV(1).
C                         ON OUTPUT XY(I,3) CONTAINS STAT(3)**XY(I,1)
C                           FOR I=1,2,...,IV(1).
C                IXY    - INPUT ROW DIMENSION OF THE MATRIX XY EXACTLY
C                           AS SPECIFIED IN THE DIMENSION STATEMENT IN
C                           THE CALLING PROGRAM.
C                IV     - INPUT VECTOR OF LENGTH 3.
C                         IV(1) CONTAINS THE NUMBER OF DATA POINTS IN
C                           XY. WHEN THE MODEL INCLUDES ALPHA (IV(2) IS
C                           EQUAL TO ZERO), IV(1) MUST BE
C                           GREATER THAN OR EQUAL TO 4.
C                           WHEN THE MODEL DOES NOT INCLUDE ALPHA (IV(2)
C                           IS NON-ZERO), IV(1) MUST BE
C                           GREATER THAN OR EQUAL TO 3.
C                         IV(2) CONTAINS THE MODEL TYPE OPTION.
C                           IV(2) EQUAL TO ZERO IMPLIES THAT THE
C                             MODEL INCLUDES ALPHA.
C                           IV(2) NON-ZERO IMPLIES THAT THE
C                             MODEL DOES NOT INCLUDE ALPHA.
C                         IV(3) CONTAINS THE INTERVAL OPTION.
C                           IV(3) EQUAL TO ZERO IMPLIES THE
C                             PERMISSIBLE RANGE FOR THE ESTIMATE OF
C                             GAMMA WILL BE THE EXCLUSIVE RANGE, (0,1).
C                           IV(3) NON-ZERO IMPLIES A SUBSET OF
C                             (0,1) MUST BE SPECIFIED. SEE
C                             DESCRIPTION OF PARAMETER, PAR, BELOW.
C                PAR    - INPUT VECTOR OF LENGTH 3. IF IV(3) IS EQUAL TO
C                           ZERO, ONLY PAR(3) NEED BE INITIALIZED.
C                         PAR(1) CONTAINS THE LOWER LIMIT OF THE
C                           PERMISSIBLE INTERVAL FOR THE GAMMA ESTIMATE
C                           WHEN IV(3) IS NON-ZERO.
C                         PAR(2) CONTAINS THE UPPER LIMIT OF THE
C                           PERMISSIBLE INTERVAL WHEN IV(3) IS NON-ZERO.
C                         PAR(3) CONTAINS THE MAXIMUM DEVIATION OF
C                           STAT(3) FROM THE LEAST SQUARES ESTIMATE OF
C                           GAMMA. PAR(3) MUST BE POSITIVE.
C                STAT   - OUTPUT VECTOR OF LENGTH 6.
C                         STAT(1) CONTAINS THE ALPHA ESTIMATE
C                         STAT(2) CONTAINS THE BETA ESTIMATE
C                         STAT(3) CONTAINS THE GAMMA ESTIMATE
C                         STAT(4) CONTAINS THE AVERAGE OF XY(I,2),
C                           I=1,2,...,IV(1).
C                         STAT(5) CONTAINS THE ERROR SUM OF SQUARES
C                         STAT(6) CONTAINS THE AVERAGE OF XY(I,3),
C                           I=1,2,...,IV(1).
C                VCV    - OUTPUT VECTOR CONTAINING THE ESTIMATE OF THE
C                           LARGE SAMPLE VARIANCE-COVARIANCE MATRIX OF
C                           THE PARAMETER ESTIMATES STORED IN SYMMETRIC
C                           STORAGE MODE.
C                           FOR IV(2) EQUAL TO ZERO, VCV IS OF LENGTH 6
C                           (A 3 BY 3 SYMMETRIC MATRIX).
C                           FOR IV(2) NON-ZERO, VCV IS OF LENGTH 3
C                           (2 BY 2 SYMMETRIC MATIRX).
C                ITER   - OUTPUT VECTOR OF LENGTH TWO.
C                         ITER(1) CONTAINS THE NUMBER OF ITERATIONS
C                           REQUIRED TO FIND STAT(3) TO THE ACCURACY
C                           SPECIFIED BY PAR(3).
C                         ITER(2) IS USED AS WORK STORAGE AND ON
C                           OUTPUT CONTAINS THE SAME VALUE AS ITER(1).
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         WARNING ERROR
C                           IER=33 INDICATES THAT IMSL ROUTINE ZXGSP
C                             DETERMINED THAT THE FUNCTION INTERNAL
C                             TO RSMITZ WHICH COMPUTES THE ERROR SUM
C                             OF SQUARES FOR THE MODEL WAS NOT UNIMODAL
C                             DUE TO ROUNDING ERRORS. THE RESULTS
C                             RETURNED ARE CORRECT.
C                         TERMINAL ERROR
C                           IER=130 INDICATES THE NUMBER OF DATA POINTS,
C                             IV(1), IS LESS THAN 4 WHEN THE MODEL
C                             INCLUDES ALPHA OR LESS THAN 3 WHEN THE
C                             MODEL DOES NOT INCLUDE ALPHA.
C                           IER=131 INDICATES AT LEAST TWO ELEMENTS OF
C                             COLUMN 1 OF THE MATRIX XY WERE OF
C                             DIFFERENT SIGNS.
C                           IER=132 INDICATES THAT, FOR IV(3) NON-ZERO,
C                             EITHER THE INTERVAL SPECIFIED IN PAR(1)
C                             AND PAR(2) WAS NOT A SUBSET OF (0,1),
C                             OR PAR(2) - PAR(1) WAS NOT POSITIVE.
C                           IER=133 INDICATES THAT PAR(3) IS NOT
C                             POSITIVE OR THAT PAR(3) EXCEEDS THE
C                             PERMISSIBLE RANGE FOR THE GAMMA ESTIMATE.
C                           IER=134 INDICATES THAT THE LEAST SQUARES
C                             PARAMETER ESTIMATES OR THE DATA WERE OF
C                             SUCH A MAGNITUDE (FOR EXAMPLE, STAT(3)
C                             VERY CLOSE TO ONE, OR XY(I,1) VERY
C                             LARGE) THAT COMPUTATIONS FOR VCV COULD
C                             NOT BE CARRIED OUT. IN THIS CASE IT IS
C                             LIKELY THAT THE MODEL IS INCORRECT FOR
C                             THE DATA.
C
C   REQD. IMSL ROUTINES - SINGLE/RSMSSE,UERTST,UERSET,UGETIO,ZXGSP
C                       - DOUBLE/RSMSSE,UERTST,UERSET,UGETIO,VXADD,
C                           VXMUL,VXSTO,ZXGSP
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE RSMITZ (XY,IXY,IV,PAR,STAT,VCV,ITER,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            IXY,IER,IV(3),ITER(2)
      REAL               XY(IXY,3),PAR(3),STAT(6),VCV(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I1,I2,I3,I4,I5,I6,INTVL,I,JER,KER,LEVEL,LEVOLD,
     1                   MODEL,N
      REAL               AA,BB,BSQ,BSTR,GEST,ONE,Q,RDELP,RN,SGSQ,SPR,
     1                   SPSPR,T2,T3,TEN,THREE,TOL,TWO,XMIN,XNF,ZERO
      DOUBLE PRECISION   A,B,C,D,F,T1
      DOUBLE PRECISION   DZERO
      EXTERNAL           RSMSSE
      DATA               I1 /1/,I2 /2/,I3 /3/,I4 /4/,I5 /5/,I6 /6/
      DATA               ZERO /0.0/,ONE /1.0/,TWO /2.0/,
     *                   THREE /3.0/,TEN /10.0/
      DATA               DZERO /0.0D0/
      DATA               RDELP/Z3C100000/
      DATA               SPR/Z00100000/
      DATA               XNF/Z7FFFFFFF/
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
      N = IV(1)
      MODEL = IV(2)
      INTVL = IV(3)
      RN = N
      IF (N.GT.3 .OR. (N.EQ.3 .AND. MODEL.NE.0)) GO TO 5
C                                  TERMINAL - INSUFFICIENT NUMBER OF
C                                             DATA POINTS
      IER = 130
      GO TO 9000
    5 IF (XY(1,1).GE.ZERO) GO TO 15
      DO 10 I=2,N
         IF (XY(I,1).LT.ZERO) GO TO 10
C                                  TERMINAL - AT LEAST 2 ELEMENTS OF
C                                    XY(I,1) WERE OF DIFFERENT SIGNS
         IER = 131
         GO TO 9000
   10 CONTINUE
      GO TO 25
   15 DO 20 I=2,N
         IF (XY(I,1).GE.ZERO) GO TO 20
C                                  TERMINAL - AT LEAST 2 ELEMENTS OF
C                                    XY(I,1) WERE OF DIFFERENT SIGNS
         IER = 131
         GO TO 9000
   20 CONTINUE
   25 IF (INTVL.EQ.0) GO TO 40
      IF (PAR(2)-PAR(1)) 35, 35, 30
   30 IF (PAR(1).GT.ZERO .AND. PAR(2).LT.ONE) GO TO 45
C                                  TERMINAL - INVALID INTERVAL
   35 IER = 132
      GO TO 9000
   40 PAR(1) = RDELP
      PAR(2) = ONE-RDELP
   45 SPSPR = SQRT(SPR*TEN)
C                                  MINIMIZE THE FUNCTION
      IF (PAR(3).GT.ZERO .AND. PAR(3).LE.PAR(2)-PAR(1)) GO TO 50
C                                  TERMINAL - INVALID RANGE FOR GAMMA
C                                             ESTIMATE
      IER = 133
      GO TO 9000
   50 AA = PAR(1)
      BB = PAR(2)
      TOL = PAR(3)
      LEVEL = 0
      CALL UERSET(LEVEL,LEVOLD)
      CALL ZXGSP(RSMSSE,XY,STAT,IV,IXY,JER,AA,BB,TOL,XMIN,KER)
      CALL UERSET(LEVOLD,LEVEL)
      ITER(1) = (ALOG(PAR(3)/(PAR(2)-PAR(1))))/ALOG((SQRT(5.0)-1.0)
     */2.0)+1.0
      ITER(2) = ITER(1)
      IF (JER.LE.128) GO TO 55
      IER = JER
      GO TO 9000
C                                  WARNING - FAILURE TO CONVERGE IN
C                                    ZXGSP DUE TO ROUNDING ERROR
   55 IF (KER.GE.131) IER = 33
      GEST = STAT(3)
      SGSQ = STAT(5)
      BSTR = -STAT(2)/GEST
      BSQ = BSTR**2
      A = RN*STAT(6)
      B = DZERO
      C = DZERO
      D = DZERO
      F = DZERO
      DO 80 I=1,N
         T1 = DBLE(XY(I,3))*DBLE(XY(I,3))
         B = B+T1
         C = C+DBLE(XY(I,1))*DBLE(XY(I,3))
         D = D+XY(I,1)*T1
         F = F+DBLE(XY(I,1))**2*T1
   80 CONTINUE
      IF (B.GT.SPSPR .AND. F.GT.SPSPR .AND. D.GT.SPSPR) GO TO 85
      IER = 134
      GO TO 9000
   85 IF (MODEL.NE.0) GO TO 100
      IF (A.GT.SPSPR .AND. C.GT.SPSPR) GO TO 90
      IER = 134
      GO TO 9000
   90 Q = RN*B*F-A**2*F-B*C**2+TWO*A*C*D-RN*D**2
      IF (Q.LT.ONE) GO TO 95
      IF (BSQ.LT.XNF/Q) GO TO 95
      IER = 134
      GO TO 9000
   95 Q = BSQ*Q
      T2 = SGSQ/(Q*(RN-THREE))
      T3 = T2*BSQ
C                                  VARIANCE-COVARIANCE MATRIX ESTIMATE
C                                    FOR ALPHA MODEL
      VCV(I1) = T3*(B*F-D**2)
      VCV(I2) = T3*(A*F-C*D)
      VCV(I3) = T3*(RN*F-C**2)
      VCV(I4) = T2*BSTR*(B*C-A*D)
      VCV(I5) = T2*BSTR*(A*C-RN*D)
      VCV(I6) = T2*(RN*B-A**2)
      GO TO 9000
  100 Q = BSQ*(D**2-B*F)
      T2 = SGSQ/(Q*(RN-TWO))
C                                  VARIANCE-COVARIANCE MATRIX ESTIMATE
C                                    FOR MODEL WITHOUT ALPHA
      VCV(I1) = -BSQ*T2*F
      VCV(I2) = T2*BSTR*D
      VCV(I3) = -T2*B
 9000 CONTINUE
      IF (IER.NE.0) CALL UERTST(IER,6HRSMITZ)
      RETURN
      END
