C   IMSL ROUTINE NAME   - FTMA1
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1982
C
C   PURPOSE             - NUCLEUS USED ONLY BY IMSL ROUTINE FTMA
C
C   REQD. IMSL ROUTINES - NONE REQUIRED
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE FTMA1 (X,F,N,PAR)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N
      REAL               X(N),F(N),PAR(N)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            IM1,J,K
C                                  FIRST EXECUTABLE STATEMENT
      DO 15 I=1,N
         IM1 = I-1
         F(I) = -PAR(I)
         J = N-I+1
    5    DO 10 K=1,J
            F(I) = F(I)+X(K)*X(IM1+K)
   10    CONTINUE
   15 CONTINUE
      RETURN
      END
