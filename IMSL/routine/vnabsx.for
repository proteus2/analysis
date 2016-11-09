C   IMSL ROUTINE NAME   - VNABSX
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
      SUBROUTINE VNABSX (N,KX,INCX)
C
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N,INCX,KX(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,IX,IY,M,MP1,NS
C                                  FIRST EXECUTABLE STATEMENT
      IF (N.LE.0) RETURN
      IF (INCX-1) 5, 15, 35
C                                  CODE FOR UNEQUAL OR NONPOSITIVE
C                                    INCREMENTS.
    5 IX = 1
      IF (INCX.LT.0) IX = (-N+1)*INCX+1
      DO 10 I=1,N
         IF (KX(IX).LT.0) KX(IX) = -KX(IX)
         IX = IX+INCX
   10 CONTINUE
      RETURN
C                                  CODE FOR BOTH INCREMENTS EQUAL TO 1
C                                    CLEAN-UP LOOP SO REMAINING VECTOR
C                                    LENGTH IS A MULTIPLE OF 7.
   15 M = N-(N/7)*7
      IF (M.EQ.0) GO TO 25
      DO 20 I=1,M
         IF (KX(I).LT.0) KX(I) = -KX(I)
   20 CONTINUE
      IF (N.LT.7) RETURN
   25 MP1 = M+1
      DO 30 I=MP1,N,7
         IF (KX(I).LT.0) KX(I) = -KX(I)
         IF (KX(I+1).LT.0) KX(I+1) = -KX(I+1)
         IF (KX(I+2).LT.0) KX(I+2) = -KX(I+2)
         IF (KX(I+3).LT.0) KX(I+3) = -KX(I+3)
         IF (KX(I+4).LT.0) KX(I+4) = -KX(I+4)
         IF (KX(I+5).LT.0) KX(I+5) = -KX(I+5)
         IF (KX(I+6).LT.0) KX(I+6) = -KX(I+6)
   30 CONTINUE
      RETURN
C                                  CODE FOR EQUAL, POSITIVE, NONUNIT
C                                    INCREMENTS.
   35 CONTINUE
      NS = N*INCX
      DO 40 I=1,NS,INCX
         IF (KX(I).LT.0) KX(I) = -KX(I)
   40 CONTINUE
      RETURN
      END
