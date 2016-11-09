C   IMSL ROUTINE NAME   - FTGEN
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1982
C
C   PURPOSE             - GENERATION OF A TIME SERIES FROM A GIVEN ARIMA
C                           (BOX-JENKINS) STOCHASTIC MODEL
C
C   USAGE               - CALL FTGEN (ARPS,PMAS,PMAC,START,WNV,DSEED,
C                           IP,IQ,LW,W,WA)
C
C   ARGUMENTS    ARPS   - INPUT VECTOR OF LENGTH IP CONTAINING THE
C                           AUTOREGRESSIVE PARAMETERS OF THE MODEL.
C                PMAS   - INPUT VECTOR OF LENGTH IQ CONTAINING
C                           THE MOVING AVERAGE PARAMETERS OF THE MODEL.
C                PMAC   - INPUT OVERALL MOVING AVERAGE PARAMETER.
C                START  - INPUT VECTOR OF LENGTH IP CONTAINING STARTING
C                           VALUES WITH WHICH TO GENERATE THE TIME
C                           SERIES.
C                WNV    - INPUT WHITE NOISE VARIANCE.
C                DSEED  - INPUT/OUTPUT.  DOUBLE PRECISION VARIABLE
C                           ASSIGNED AN INTEGER VALUE IN THE
C                           EXCLUSIVE RANGE (1.D0, 2147483647.D0).
C                           DSEED IS REPLACED BY A NEW VALUE TO BE
C                           USED IN A SUBSEQUENT CALL. DSEED IS USED
C                           IN SIMULATING THE TIME SERIES, SO IT SHOULD
C                           BE CHOSEN RANDOMLY.
C                IP     - INPUT NUMBER OF AUTOREGRESSIVE PARAMETERS
C                           IN THE MODEL.
C                IQ     - INPUT NUMBER OF MOVING AVERAGE PARAMETERS IN
C                           THE MODEL.
C                LW     - INPUT LENGTH OF THE TIME SERIES TO BE
C                           GENERATED.
C                W      - OUTPUT VECTOR OF LENGTH LW CONTAINING THE
C                           GENERATED TIME SERIES.
C                WA     - WORK AREA VECTOR OF LENGTH
C                           LW+MAXIMUM(IP,IQ).
C
C   REQD. IMSL ROUTINES - GGNML,GGUBS,MDNRIS,MERFI,UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS  1.  EITHER OF THE INPUT VALUES IP OR IQ MAY BE ZERO.
C            2.  IF WNV IS EQUAL TO ZERO, THE MODEL REDUCES TO
C
C                  W(I) = THE SUM FROM J = 1 TO IP OF
C                         ARPS(J)*W(I-J)+PMAC
C                  FOR I=1,...,LW AND THE VALUES START(K),K=1,...,IP
C                  CORRESPONDING TO W(-IP+1),...,W(0).
C
C-----------------------------------------------------------------------
C
      SUBROUTINE FTGEN  (ARPS,PMAS,PMAC,START,WNV,DSEED,IP,IQ,LW,W,WA)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            IP,IQ,LW
      REAL               ARPS(1),PMAS(1),PMAC,START(1),WNV,W(LW),
     1                   WA(1)
      DOUBLE PRECISION   DSEED
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,IIP,IIQ,J,LT,LZ
      REAL               SWNV
C                                  FIRST EXECUTABLE STATEMENT
      LZ=LW+IQ
C                                  GENERATE WHITE NOISE SERIES
      LT=LZ
      CALL GGNML (DSEED,LT,WA)
      SWNV=SQRT(WNV)
      DO 5  I=1,LZ
         WA(I)=WA(I)*SWNV
    5 CONTINUE
      DO 15 I=1,LW
         IIQ=I+IQ
         W(I)=PMAC+WA(IIQ)
         IF (IQ.EQ.0) GO TO 15
C                                  COMPUTE MOVING AVERAGE CONTRIBUTIONS
         DO 10 J=1,IQ
            W(I)=W(I)-PMAS(J)*WA(IIQ-J)
   10    CONTINUE
   15 CONTINUE
      IF (IP.EQ.0) GO TO 35
      DO 20 I=1,IP
         WA(I)=START(I)
   20 CONTINUE
C                                  COMPUTE AUTOREGRESSIVE CONTRIBUTIONS
      DO 25 I=1,LW
         IIP=I+IP
         WA(IIP)=W(I)
         DO 25 J=1,IP
            WA(IIP)=WA(IIP)+ARPS(J)*WA(IIP-J)
   25 CONTINUE
      DO 30 I=1,LW
         W(I)=WA(I+IP)
   30 CONTINUE
   35 RETURN
      END
