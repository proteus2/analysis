C   IMSL ROUTINE NAME   - VNDIST
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1982
C
C   PURPOSE             - NUCLEUS CALLED BY IMSL ROUTINES WHICH SORT
C                           A MATRIX WITH RESPECT TO VECTORS.
C
C   REQD. IMSL ROUTINES - VNINI
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE VNDIST (IR,LA,IRD,LD)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            LA,IRD,LD,IR(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,M,MP1,NZ
C                                  FIRST EXECUTABLE STATEMENT
      IF (LA.LE.0) RETURN
C                                  CLEAN-UP LOOP SO REMAINING VECTOR
C                                  LENGTH IS A MULTIPLE OF 7.
      LD = 0
      M = LA-(LA/7)*7
      IF (M.EQ.0) GO TO 10
      DO 5 I=1,M
         IF (IR(I).GE.0) GO TO 5
         IR(I) = -IR(I)
         LD = LD+1
         IR(IRD-1+LD) = I
    5 CONTINUE
      IF (LA.LT.7) GO TO 50
   10 MP1 = M+1
      DO 45 I=MP1,LA,7
         IF (IR(I).GE.0) GO TO 15
         IR(I) = -IR(I)
         LD = LD+1
         IR(IRD-1+LD) = I
   15    IF (IR(I+1).GE.0) GO TO 20
         IR(I+1) = -IR(I+1)
         LD = LD+1
         IR(IRD-1+LD) = I+1
   20    IF (IR(I+2).GE.0) GO TO 25
         IR(I+2) = -IR(I+2)
         LD = LD+1
         IR(IRD-1+LD) = I+2
   25    IF (IR(I+3).GE.0) GO TO 30
         IR(I+3) = -IR(I+3)
         LD = LD+1
         IR(IRD-1+LD) = I+3
   30    IF (IR(I+4).GE.0) GO TO 35
         IR(I+4) = -IR(I+4)
         LD = LD+1
         IR(IRD-1+LD) = I+4
   35    IF (IR(I+5).GE.0) GO TO 40
         IR(I+5) = -IR(I+5)
         LD = LD+1
         IR(IRD-1+LD) = I+5
   40    IF (IR(I+6).GE.0) GO TO 45
         IR(I+6) = -IR(I+6)
         LD = LD+1
         IR(IRD-1+LD) = I+6
   45 CONTINUE
   50 CONTINUE
      IR(IRD+LA) = LD
      NZ = LA-LD
      CALL VNINI(NZ,IR(IRD+LD),1,0,0)
      RETURN
      END
