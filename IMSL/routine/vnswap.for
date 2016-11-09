C   IMSL ROUTINE NAME   - VNSWAP
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1982
C
C   PURPOSE             - NUCLEUS CALLED BY IMSL ROUTINES WHICH SORT
C                           A MATRIX WITH RESPECT TO VECTORS.
C
C   REQD. IMSL ROUTINES - NONE REQUIRED
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE VNSWAP (N,KX,INCX,KY,INCY)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N,INCX,INCY,KX(1),KY(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,IX,IY,KTEMP1,KTEMP2,KTEMP3,M,MP1,NS
C                                  FIRST EXECUTABLE STATEMENT
      IF (N.LE.0) RETURN
      IF (INCX.EQ.INCY) IF (INCX-1) 5, 15, 35
    5 CONTINUE
C                                  CODE FOR UNEQUAL OR NONPOSITIVE
C                                    INCREMENTS.
      IX = 1
      IY = 1
      IF (INCX.LT.0) IX = (-N+1)*INCX+1
      IF (INCY.LT.0) IY = (-N+1)*INCY+1
      DO 10 I=1,N
         KTEMP1 = KX(IX)
         KX(IX) = KY(IY)
         KY(IY) = KTEMP1
         IX = IX+INCX
         IY = IY+INCY
   10 CONTINUE
      RETURN
C                                  CODE FOR BOTH INCREMENTS EQUAL TO 1
C                                    CLEAN-UP LOOP SO REMAINING VECTOR
C                                    LENGTH IS A MULTIPLE OF 3.
   15 M = N-(N/3)*3
      IF (M.EQ.0) GO TO 25
      DO 20 I=1,M
         KTEMP1 = KX(I)
         KX(I) = KY(I)
         KY(I) = KTEMP1
   20 CONTINUE
      IF (N.LT.3) RETURN
   25 MP1 = M+1
      DO 30 I=MP1,N,3
         KTEMP1 = KX(I)
         KTEMP2 = KX(I+1)
         KTEMP3 = KX(I+2)
         KX(I) = KY(I)
         KX(I+1) = KY(I+1)
         KX(I+2) = KY(I+2)
         KY(I) = KTEMP1
         KY(I+1) = KTEMP2
         KY(I+2) = KTEMP3
   30 CONTINUE
      RETURN
   35 CONTINUE
C                                  CODE FOR EQUAL, POSITIVE, NONUNIT
C                                    INCREMENTS.
      NS = N*INCX
      DO 40 I=1,NS,INCX
         KTEMP1 = KX(I)
         KX(I) = KY(I)
         KY(I) = KTEMP1
   40 CONTINUE
      RETURN
      END
