C   IMSL ROUTINE NAME   - ZCPQLG
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - NUCLEUS CALLED ONLY BY IMSL SUBROUTINE
C                           ZCPOLY
C
C   REQD. IMSL ROUTINES - NONE REQUIRED
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE ZCPQLG (NN,SR,SI,PR,PI,QR,QI,PVR,PVI)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            NN
      REAL               PR(NN),PI(NN),QR(NN),QI(NN),
     1                   SR,SI,PVR,PVI
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I
      REAL               T
C                                  FIRST EXECUTABLE STATEMENT
      QR(1) = PR(1)
      QI(1) = PI(1)
      PVR = QR(1)
      PVI = QI(1)
C                                  EVALUATE A POLYNOMIAL P AT S BY THE
C                                    HORNER RECURRENCE PLACING THE
C                                    PARTIAL SUMS IN Q AND THE COMPUTED
C                                    VALUE IN PV.
      DO 5 I=2,NN
         T = PVR*SR-PVI*SI+PR(I)
         PVI = PVR*SI+PVI*SR+PI(I)
         PVR = T
         QR(I) = PVR
         QI(I) = PVI
    5 CONTINUE
      RETURN
      END
