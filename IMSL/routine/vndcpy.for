C   IMSL ROUTINE NAME   - VNDCPY
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
      SUBROUTINE VNDCPY (N,KX,INCX,SY,INCY)
C
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N,INCX,INCY,KX(1)
      DOUBLE PRECISION   SY(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,IX,IY,M,MP1,NS
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
         SY(IY) = KX(IX)
         IX = IX+INCX
         IY = IY+INCY
   10 CONTINUE
      RETURN
C                                  CODE FOR BOTH INCREMENTS EQUAL TO 1
C                                    CLEAN-UP LOOP SO REMAINING VECTOR
C                                    LENGTH IS A MULTIPLE OF 7.
   15 M = N-(N/7)*7
      IF (M.EQ.0) GO TO 25
      DO 20 I=1,M
         SY(I) = KX(I)
   20 CONTINUE
      IF (N.LT.7) RETURN
   25 MP1 = M+1
      DO 30 I=MP1,N,7
         SY(I) = KX(I)
         SY(I+1) = KX(I+1)
         SY(I+2) = KX(I+2)
         SY(I+3) = KX(I+3)
         SY(I+4) = KX(I+4)
         SY(I+5) = KX(I+5)
         SY(I+6) = KX(I+6)
   30 CONTINUE
      RETURN
C                                  CODE FOR EQUAL, POSITIVE, NONUNIT
C                                    INCREMENTS.
   35 CONTINUE
      NS = N*INCX
      DO 40 I=1,NS,INCX
         SY(I) = KX(I)
   40 CONTINUE
      RETURN
      END
