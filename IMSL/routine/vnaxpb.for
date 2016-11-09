C   IMSL ROUTINE NAME   - VNAXPB
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
      SUBROUTINE VNAXPB (N,KA,KX,INCX,KB,INCB)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N,KA,INCX,KB,INCB,KX(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,IB,IX,M,MP1,NS
C                                  FIRST EXECUTABLE STATEMENT
      IF (N.LE.0) RETURN
      IF (INCX.EQ.INCB) IF (INCX-1) 5, 15, 35
    5 CONTINUE
C                                  CODE FOR NONEQUAL OR NONPOSITIVE
C                                    INCREMENTS.
      IX = 1
      IB = 0
      IF (INCX.LT.0) IX = (-N+1)*INCX+1
      DO 10 I=1,N
         KX(IX) = KB+IB+KA*KX(IX)
         IX = IX+INCX
         IB = IB+INCB
   10 CONTINUE
      RETURN
C                                  CODE FOR BOTH INCREMENTS EQUAL TO 1
C                                    CLEAN-UP LOOP SO REMAINING VECTOR
C                                    LENGTH IS A MULTIPLE OF 4.
   15 M = N-(N/4)*4
      IF (M.EQ.0) GO TO 25
      DO 20 I=1,M
         KX(I) = KB+I-1+KA*KX(I)
   20 CONTINUE
      IF (N.LT.4) RETURN
   25 MP1 = M+1
      DO 30 I=MP1,N,4
         KX(I) = KB+I-1+KA*KX(I)
         KX(I+1) = KB+I+KA*KX(I+1)
         KX(I+2) = KB+I+1+KA*KX(I+2)
         KX(I+3) = KB+I+2+KA*KX(I+3)
   30 CONTINUE
      RETURN
C                                  CODE FOR EQUAL, POSITIVE, NONUNIT
C                                    INCREMENTS.
   35 CONTINUE
      NS = N*INCX
      DO 40 I=1,NS,INCX
         KX(I) = KA*KX(I)+KB+I-1
   40 CONTINUE
      RETURN
      END
