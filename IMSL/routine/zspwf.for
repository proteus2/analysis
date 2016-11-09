C   IMSL ROUTINE NAME   - ZSPWF
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1982
C
C   PURPOSE             - NUCLEUS CALLED ONLY BY IMSL SUBROUTINE ZSPOW
C
C   REQD. IMSL ROUTINES - NONE REQUIRED
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE ZSPWF (M,N,Q,LDQ,WA)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            M,N,LDQ
      REAL               Q(LDQ,M),WA(M)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,JM1,J,K,L,MINMN,NP1
      REAL               ONE,SUM,TEMP,ZERO
      DATA               ONE,ZERO /1.0E0,0.0E0/
C                                  ZERO OUT UPPER TRIANGLE OF Q IN THE
C                                  FIRST MIN(M,N) COLUMNS.
C                                  FIRST EXECUTABLE STATEMENT
      MINMN = MIN0(M,N)
      IF (MINMN.LT.2) GO TO 15
      DO 10 J=2,MINMN
         JM1 = J-1
         DO 5 I=1,JM1
            Q(I,J) = ZERO
    5    CONTINUE
   10 CONTINUE
   15 CONTINUE
C                                  INITIALIZE REMAINING COLUMNS TO THOSE
C                                  OF THE IDENTITY MATRIX.
      NP1 = N+1
      IF (M.LT.NP1) GO TO 30
      DO 25 J=NP1,M
         DO 20 I=1,M
            Q(I,J) = ZERO
   20    CONTINUE
         Q(J,J) = ONE
   25 CONTINUE
   30 CONTINUE
C                                  ACCUMULATE Q FROM ITS FACTORED FORM.
      DO 60 L=1,MINMN
         K = MINMN-L+1
         DO 35 I=K,M
            WA(I) = Q(I,K)
            Q(I,K) = ZERO
   35    CONTINUE
         Q(K,K) = ONE
         IF (WA(K).EQ.ZERO) GO TO 55
         DO 50 J=K,M
            SUM = ZERO
            DO 40 I=K,M
               SUM = SUM+Q(I,J)*WA(I)
   40       CONTINUE
            TEMP = SUM/WA(K)
            DO 45 I=K,M
               Q(I,J) = Q(I,J)-TEMP*WA(I)
   45       CONTINUE
   50    CONTINUE
   55    CONTINUE
   60 CONTINUE
      RETURN
      END
