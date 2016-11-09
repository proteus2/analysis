C   IMSL ROUTINE NAME   - FTARPS
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1982
C
C   PURPOSE             - PRELIMINARY ESTIMATION OF THE AUTOREGRESSIVE
C                           PARAMETERS IN AN ARIMA STOCHASTIC MODEL.
C
C   USAGE               - CALL FTARPS (ACV,WBAR,IP,IQ,ARPS,PMAC,WA,IER)
C
C   ARGUMENTS    ACV    - INPUT VECTOR OF LENGTH IP+IQ+1 CONTAINING THE
C                           AUTOCOVARIANCES OF THE TIME SERIES
C                           BEING MODELED. ACV(I) IS THE AUTOCOVARIANCE
C                           CORRESPONDING TO A TIME LAG OF I-1 UNITS OF
C                           TIME.
C                WBAR   - INPUT MEAN OF THE TIME SERIES.
C                IP     - INPUT NUMBER OF AUTOREGRESSIVE PARAMETERS IN
C                           THE MODEL. IP MUST BE GREATER THAN OR
C                           EQUAL TO 1.
C                IQ     - INPUT NUMBER OF MOVING AVERAGE PARAMETERS IN
C                           THE MODEL. IQ MUST BE GREATER THAN OR EQUAL
C                           TO ZERO.
C                ARPS   - OUTPUT VECTOR OF LENGTH IP.  ARPS CONTAINS
C                           THE PRELIMINARY ESTIMATES OF THE
C                           AUTOREGRESSIVE PARAMETERS OF THE MODEL.
C                PMAC   - OUTPUT OVERALL MOVING AVERAGE CONSTANT.
C                WA     - WORK AREA VECTOR OF LENGTH IP**2+5*IP.
C                IER    - ERROR PARAMETER (OUTPUT).
C                         TERMINAL ERROR
C                           IER = 129 INDICATES THE INPUT DATA
C                             GENERATED A SINGULAR SYSTEM. THIS USUALLY
C                             INDICATES INVALID VALUES IN THE ACV
C                             VECTOR OR A TIME SERIES WHICH IS
C                             NONSTATIONARY.
C
C   REQD. IMSL ROUTINES - LEQT1F,LUDATN,LUELMN,UERTST,UGETIO,VABMXF
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE FTARPS (ACV,WBAR,IP,IQ,ARPS,PMAC,WA,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            IP,IQ,IER
      REAL               ACV(1),ARPS(IP),WA(1),WBAR,PMAC
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            IDGT,I,KK,J,N,M,K,L,KTY
      REAL               ZERO,DET,VMAX1,VMAX2,ONE,HUND
      DATA               ZERO,ONE/0.,1./,IDGT/4/,HUND/100.0/
C                                  FIRST EXECUTABLE STATEMENT
      K = IQ + 1
      IER = 0
      IF (IP .GT. 2) GO TO 10
      IF (IP .EQ. 2) GO TO 5
      IF (ACV(K) .EQ. ZERO) GO TO 45
C                                  1 BY 1 SYSTEM
      ARPS(1) = ACV(K+1)/ACV(K)
      GO TO 35
C                                  2 BY 2 SYSTEM
    5 L = IABS(IQ-1) + 1
      DET = ACV(K)**2 - ACV(L)*ACV(K+1)
      IF (DET .EQ. ZERO) GO TO 45
      ARPS(1) = (ACV(K+1)*ACV(K) - ACV(L)*ACV(K+2))/DET
      ARPS(2) = (ACV(K)*ACV(K+2) - ACV(K+1)**2)/DET
      GO TO 35
C                                  SET UP MATRIX OF AUTO-
C                                  COVARIANCES (3 BY 3 OR MORE SYSTEM)
   10 KK = IP**2
      DO 15 I = 1,IP
         WA(I) = ACV(IQ+I)
         WA(KK+I) = ACV(K+I)
   15 CONTINUE
      M = 0
      N = IP + 1
      DO 25 J = 2,IP
         DO 20 I = 2,IP
            M = M+1
            N = N+1
            WA(N) = WA(M)
   20    CONTINUE
         M = M + 1
         N = N + 1
         KTY = IABS(K-J)+1
         WA(M+1) = ACV(KTY)
   25 CONTINUE
      N = KK + IP + 1
      CALL VABMXF(WA(KK+1),IP,1,J,VMAX1)
C                                  COMPUTE THE SOLUTION
      CALL LEQT1F(WA,1,IP,IP,WA(KK+1),IDGT,WA(N),IER)
      IF (IER .NE. 0) GO TO 45
      CALL VABMXF(WA(KK+1),IP,1,J,VMAX2)
      IF(VMAX2.GT.HUND*VMAX1) GO TO 45
C                                  MOVE THE SOLUTION TO ARPS
      DO 30  I = 1,IP
         ARPS(I) = WA(KK+I)
   30 CONTINUE
C                                  COMPUTE PMAC FORM ARPS
   35 PMAC = ONE
      DO 40 I=1,IP
         PMAC = PMAC - ARPS(I)
   40 CONTINUE
      PMAC = PMAC*WBAR
      GO TO 9005
   45 IER = 129
 9000 CONTINUE
      CALL UERTST(IER,'FTARPS')
 9005 RETURN
      END
