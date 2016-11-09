C   IMSL ROUTINE NAME   - RSMSSE
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1982
C
C   PURPOSE             - NUCLEUS USED TO PROVIDE FUNCTIONAL
C                           COMMUNICATION BETWEEN IMSL SUBROUTINES
C                           RSMITZ AND ZXGSP
C
C   REQD. IMSL ROUTINES - SINGLE/NONE REQUIRED
C                       - DOUBLE/VXADD,VXMUL,VXSTO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      REAL FUNCTION RSMSSE (G,XY,STAT,IOP,IXY,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            IXY,IER,IOP(1)
      REAL               G,XY(IXY,1),STAT(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,MODEL,N
      REAL               ABD,ABN,ONE,RNLGZ,RN,SPR,SRSPR,T1,TEN,
     1                   XETA,XMAX,XNF,XXX,ZAVG1,ZERO
      DOUBLE PRECISION   ESSQ,ZAVG,YAVG,BN,BD,BEST,AEST
      DATA               XNF/Z7FFFFFFF/
      DATA               XETA/Z00100000/
      DATA               SPR/Z00100000/
      DATA               ONE/1.0/,TEN/10.0/,ZERO/0.0/
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
      SRSPR = SQRT(SPR*TEN)
      RNLGZ = ALOG(SRSPR)
      N     = IOP(1)
      RN    = N
      AEST = ZERO
      MODEL = IOP(2)
      IF (N - 2) 10,5,15
   5  IF (MODEL .NE. 0) GO TO 15
C                                  TERMINAL - TOO FEW POINTS
  10  IER = 130
      GO TO 9000
  15  XMAX = ZERO
      DO 25 I = 1,N
         IF (XY(I,1)*ALOG(G) .GT. RNLGZ) GO TO 20
C                                  WARNING(WITH FIX) - PREVENTING
C                                    UNDERFLOW
         IER = 65
         XY(I,3) = ZERO
         GO TO 25
  20     XY(I,3) = G**XY(I,1)
C        IF (XY(I,3).LE.SRSPR) XY(I,3) = ZERO
         IF (XY(I,3).GT.XMAX) XMAX = XY(I,3)
  25  CONTINUE
      YAVG = 0.0D0
      ZAVG = 0.0D0
C                                  SUM Z
C                                  SUM Y
      DO 35 I = 1,N
         YAVG = YAVG + XY(I,2)
         ZAVG = ZAVG + XY(I,3)
  35  CONTINUE
  40  YAVG = YAVG/RN
      ZAVG = ZAVG/RN
      ZAVG1 = ZAVG
      IF (MODEL.NE.0) ZAVG1 = ZERO
      BN = 0.0D0
      BD = 0.0D0
      ESSQ = 0.0D0
C                                  SUM THE NUMERATOR AND
C                                  DENOMINATOR FOR BETA
      DO 50 I = 1,N
         T1 = XY(I,3) - ZAVG1
         IF (ABS(T1).LE.SRSPR) GO TO 50
         BN = BN + DBLE(T1)*DBLE(XY(I,2))
         BD = BD + DBLE(T1)*DBLE(T1)
  50  CONTINUE
      ABN = DABS(BN)
      ABD = DABS(BD)
      IF (ABD.GT.ONE) GO TO 55
      XXX = XNF*ABD
      IF (ABN.GE.XXX) GO TO 60
      IF (XMAX.GE.ONE) ABN = ABN*XMAX*XMAX
      IF (ABN.GE.XXX) GO TO 60
      GO TO 65
  55  IF (ABN.GT.XETA*ABD) GO TO 65
  60  BEST = ZERO
      IER = 134
      GO TO 70
C                                  BETA ESTIMATE
  65  BEST = BN/BD
C                                  ALPHA ESTIMATE
      IF (MODEL.NE.0) GO TO 75
  70  AEST = YAVG - BEST*ZAVG
C                                  SUM OF SQUARES
  75  DO 80 I = 1,N
         ESSQ = ESSQ + (XY(I,2) - AEST - BEST*XY(I,3))**2
  80  CONTINUE
      RSMSSE = ESSQ
      STAT(1) = AEST
      STAT(2) = BEST
      STAT(3) = G
      STAT(4) = YAVG
      STAT(5) = ESSQ
      STAT(6) = ZAVG
 9000 CONTINUE
 9005 RETURN
      END
