C   IMSL ROUTINE NAME   - ZCPQLH
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - NUCLEUS CALLED ONLY BY IMSL SUBROUTINE
C                           ZCPOLY
C
C   REQD. IMSL ROUTINES - ZCPQLL
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      REAL FUNCTION ZCPQLH (NN,QR,QI,RMS,RMP,ARE,RMRE)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            NN
      REAL               QR(NN),QI(NN),RMS,RMP,ARE,RMRE
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I
      REAL               E
      REAL               ZCPQLL
C                                  FIRST EXECUTABLE STATEMENT
      E = ZCPQLL(QR(1),QI(1))*RMRE/(ARE+RMRE)
C                                  BOUNDS THE ERROR IN EVALUATING THE
C                                    POLYNOMIAL BY THE HORNER
C                                    RECURRENCE
C                                  QR,QI - THE PARTIAL SUMS
C                                  RMS - MODULUS OF THE POINT
C                                  RMP - MODULUS OF POLYNOMIAL VALUE
C                                  ARE,RMRE - ERROR BOUNDS ON COMPLEX
C                                    ADDITION AND MULTIPLICATION
      DO 5 I=1,NN
         E = E*RMS+ZCPQLL(QR(I),QI(I))
    5 CONTINUE
      ZCPQLH = E*(ARE+RMRE)-RMP*RMRE
      RETURN
      END
