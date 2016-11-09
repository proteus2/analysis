C   IMSL ROUTINE NAME   - RLEAP2
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - NUCLEUS CALLED ONLY BY IMSL SUBROUTINE
C                           RLEAP
C
C   REQD. IMSL ROUTINES - NONE REQUIRED
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE RLEAP2 (XI,KP,N)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            KP,N
      REAL               XI(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,J,IN,NN,ICR,II,JN,IJ,ICJ
      REAL               B,ONE,XINN
      DATA               ONE/1.0/
C                                  FIRST EXECUTABLE STATEMENT
      IN = (N*(N-1))/2+1
      NN = IN+N-1
      XINN = -ONE/XI(NN)
      ICR = 1
      II = 0
      DO 15 I = 1,KP
         II = II+I
         IF (I .EQ. N) GO TO 10
         B = XI(IN)*XINN
         JN = IN
         IJ = II
         ICJ = 1
         DO 5 J = I,KP
            IF (J .NE. N) XI(IJ) = XI(IJ)+B*XI(JN)
            IJ = IJ+J
            IF (J .GE. N) ICJ = J
            JN = JN+ICJ
    5    CONTINUE
         XI(IN) = B
   10    IF (I .GE. N) ICR = I
         IN = IN+ICR
   15 CONTINUE
      XI(NN) = XINN
      RETURN
      END
