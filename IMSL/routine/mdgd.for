C   IMSL ROUTINE NAME   - MDGD
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1980
C
C   PURPOSE             - NUCLEUS CALLED ONLY BY IMSL SUBROUTINES
C                           MDGC AND MDGCI
C
C   REQD. IMSL ROUTINES - DCSQDU,IQHSCU,UERSET,UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE MDGD   (F,M,IOPT,B,C,A,IER)
C
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            IER,IOPT,M
      REAL               A,B(M),C(1),F(M)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,IM1,K,LEVEL,LEVOLD,MM1
      REAL               B1,B2,H,DELTA,FX,X,HALF
      DATA               HALF/.5E0/
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
      A = 0.0
      GO TO (5,25,45,55), IOPT
C                                  EQUALLY SPACED - LINEAR
    5 B1 = B(1)
      B2 = B(2)
      H = (B2-B1)/(M-1)
      X = C(1)
      K = (X-B1)/H
      IF (K.GE.2) GO TO 10
      K = 0
      GO TO 20
   10 DO 15 I=2,K
         A = A+F(I)
   15 CONTINUE
      A = H*(A+HALF*F(1))
   20 DELTA = X-(K*H+B(1))
      A = A+(DELTA+HALF*H)*F(K+1)
      IF (X.EQ.B2) GO TO 60
      A = A+(HALF*DELTA*DELTA*(F(K+2)-F(K+1))/H)
      GO TO 60
C                                  UNEQUALLY SPACED - LINEAR
   25 IM1 = 1
      X = C(1)
      DO 35 I=2,M
         IF (B(IM1).LT.B(I)) GO TO 30
         IER = 131
         GO TO 60
   30    IF (X.LT.B(I)) GO TO 40
         A = A+(F(IM1)+F(I))*(B(I)-B(IM1))
         IM1 = I
   35 CONTINUE
      IM1 = IM1-1
   40 A = HALF*A
      I = IM1+1
      IF (X.EQ.B(I)) GO TO 60
      FX = F(IM1)+(X-B(IM1))*(F(I)-F(IM1))/(B(I)-B(IM1))
      A = A+HALF*(F(IM1)+FX)*(X-B(IM1))
      GO TO 60
C                                  EQUALLY SPACED - FITTED CURVE
   45 MM1 = M-1
      H = (B(2)-B(1))/MM1
      B(M) = B(2)
      IF (MM1.LT.3) GO TO 55
      IM1 = 1
      DO 50 I=2,M
         B(I) = B(1)+H*IM1
         IM1 = I
   50 CONTINUE
C                                  UNEQUALLY SPACED - FITTED CURVE
   55 X = C(1)
      CALL IQHSCU (B,F,M,C,M,IER)
      IF (IER.NE.0) GO TO 60
      LEVEL = 0
      CALL UERSET (LEVEL,LEVOLD)
      CALL DCSQDU (B,F,M,C,M,B(1),X,A,IER)
      IER = 0
      CALL UERSET (LEVOLD,LEVEL)
   60 RETURN
      END
