C   IMSL ROUTINE NAME   - ZRPQLH
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - NUCLEUS CALLED ONLY BY IMSL SUBROUTINE
C                           ZRPOLY
C
C   REQD. IMSL ROUTINES - NONE REQUIRED
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE ZRPQLH (NN,U,V,P,Q,RA,RB)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            NN
      REAL               P(NN),Q(NN),U,V,RA,RB
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I
      REAL               C
C                                  DIVIDES P BY THE QUADRATIC 1,U,V
C                                    PLACING THE QUOTIENT IN Q AND THE
C                                    REMAINDER IN A,B
C                                  FIRST EXECUTABLE STATEMENT
      RB = P(1)
      Q(1) = RB
      RA = P(2)-U*RB
      Q(2) = RA
      DO 5 I=3,NN
         C = P(I)-U*RA-V*RB
         Q(I) = C
         RB = RA
         RA = C
    5 CONTINUE
      RETURN
      END
