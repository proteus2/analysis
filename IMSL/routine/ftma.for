C   IMSL ROUTINE NAME   - FTMA
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1982
C
C   PURPOSE             - PRELIMINARY ESTIMATION OF THE MOVING AVERAGE
C                           PARAMETERS IN AN ARIMA STOCHASTIC MODEL
C
C   USAGE               - CALL FTMA (ACV,ARPS,IP,IQ,PMAS,WNV,WA,IER)
C
C   ARGUMENTS    ACV    - INPUT VECTOR OF LENGTH IP+IQ+1 CONTAINING THE
C                           AUTOCOVARIANCES OF THE TIME SERIES
C                           BEING MODELED. ACV(I) IS THE AUTOCOVARIANCE
C                           CORRESPONDING TO A TIME LAG OF I-1 UNITS OF
C                           TIME.
C                ARPS   - INPUT VECTOR OF LENGTH IP.  THIS VECTOR
C                           CONTAINS THE PRELIMINARY ESTIMATES OF THE
C                           AUTOREGRESSIVE PARAMETERS OF THE MODEL.
C                           THESE ESTIMATES CAN BE COMPUTED BY CALLING
C                           SUBROUTINE FTARPS PRIOR TO CALLING FTMA.
C                IP     - INPUT NUMBER OF AUTOREGRESSIVE PARAMETERS IN
C                           THE MODEL. IP MUST BE GREATER THAN OR EQUAL
C                           TO ZERO.
C                IQ     - INPUT NUMBER OF MOVING AVERAGE PARAMETERS IN
C                           THE MODEL. IQ MUST BE GREATER THAN OR EQUAL
C                           TO 1.
C                PMAS   - OUTPUT VECTOR OF LENGTH IQ.  PMAS CONTAINS
C                           THE MOVING AVERAGE PARAMETERS.
C                WNV    - OUTPUT WHITE NOISE VARIANCE.
C                WA     - WORK AREA VECTOR OF LENGTH
C                           ((IQ+1)*(3*IQ+22)/2).
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER = 130 INDICATES FAILURE OF THE ALGORITHM
C                             TO CONVERGE TO ACCEPTABLE ESTIMATES
C                             FOR THE MOVING AVERAGE PARAMETERS.
C                           IER = 131 INDICATES ACV(1) IS LESS THAN OR
C                             EQUAL TO ZERO.
C                           IER = 132 INDICATES A COMPUTED WHITE NOISE
C                             VARIANCE IS LESS THAN OR EQUAL TO ZERO.
C                             THIS IMPLIES THAT IQ IS LESS THAN OR
C                             EQUAL TO ZERO.
C
C   REQD. IMSL ROUTINES - FTMA1,UERTST,UGETIO,ZSPOW,ZSPWA,ZSPWB,ZSPWC,
C                           ZSPWD,ZSPWE,ZSPWF,ZSPWG
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE FTMA  (ACV,ARPS,IP,IQ,PMAS,WNV,WA,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            IP,IQ,IER
      REAL               ACV(1),ARPS(1),PMAS(IQ),WA(1),WNV
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,IM,IPP1,IQ1,IQ2,IQP1,IQP2,ITMA1,ITMAX,
     *                   J,K,L,NSIG
      REAL               AJ,AK,FNORM,ONEM,ZERO
      EXTERNAL           FTMA1
      DATA               ZERO/0.0/,ONEM/-1.0/,NSIG/3/,ITMAX/200/
C                                  FIRST EXECUTABLE STATEMENT
      ITMA1 = ITMAX
      IER = 0
      IQP1 = IQ+1
      IF (IP .GT. 0) GO TO 10
      DO 5  I=1,IQP1
         WA(I) = ACV(I)
    5 CONTINUE
      GO TO 30
C                                  COMPUTE MODIFIED COVARIANCES
   10 IPP1 = IP+1
      DO 25 I=1,IQP1
         IM = I
         WA(I) = ZERO
         AJ = ONEM
         DO 20 J=1,IPP1
            IF (J .GT. 1) AJ = ARPS(J-1)
            AK = ONEM
            DO 15  K=1,IPP1
               L = IABS(IM-K)+1
               IF (K .GT. 1) AK = ARPS(K-1)
               WA(I) = WA(I) + (AJ*AK*ACV(L))
   15       CONTINUE
            IM = IM+1
   20    CONTINUE
   25 CONTINUE
C                                  SET INITIAL GUESSES FOR ZSPOW
   30 IQP2 = IQP1+1
      IF(WA(1) .LE. ZERO) GO TO 60
      WA(IQP2) = SQRT(WA(1))
      IQ1 = IQP2+1
      IQ2 = IQ+IQ+2
      DO 35  I=IQ1,IQ2
         WA(I) = ZERO
   35 CONTINUE
      CALL ZSPOW (FTMA1,NSIG,IQP1,ITMA1,WA,WA(IQP2),FNORM,WA(IQ2+1),IER)
      IF (IER.EQ.0) GO TO 40
      IER=130
      GO TO 9000
   40 WNV = WA(IQP2)
      IF (WNV) 45,45,50
C                                  ERROR, WHITE NOISE VARIANCE .LE. 0
   45 IER = 132
      GO TO 9000
C                                  COMPUTE MOVING AVERAGE PARAMETERS
   50 DO 55  I=1,IQ
         PMAS(I) = -WA(IQP2+I)/WNV
   55 CONTINUE
      WNV = WNV*WNV
      GO TO 9005
C                                  DATA ERROR, FIRST COVARIANCE .LE. 0
   60 IER = 131
 9000 CONTINUE
      CALL UERTST(IER,'FTMA  ')
 9005 RETURN
      END
