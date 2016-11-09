C   IMSL ROUTINE NAME   - GGPER
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1980
C
C   PURPOSE             - GENERATE A RANDOM PERMUTATION OF THE
C                          INTEGERS 1 TO K.
C
C   USAGE               - CALL GGPER(DSEED,K,IPER)
C
C   ARGUMENTS    DSEED  - INPUT/OUTPUT DOUBLE PRECISION VARIABLE
C                           ASSIGNED AN INTEGER VALUE IN THE
C                           EXCLUSIVE RANGE (1.D0,2147483647.D0).
C                           DSEED IS REPLACED BY A NEW VALUE TO BE
C                           USED IN A SUBSEQUENT CALL.
C                K      - INPUT NUMBER OF INTEGERS TO BE PERMUTED.
C                IPER   - OUTPUT VECTOR OF LENGTH K CONTAINING THE
C                           RANDOM PERMUTATION OF THE INTEGERS 1 TO K.
C
C   REQD. IMSL ROUTINES - GGUBFS
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE GGPER  (DSEED,K,IPER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            K,IPER(K)
      DOUBLE PRECISION   DSEED
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,ITEMP,KM1,M
C                                  FIRST EXECUTABLE STATEMENT
      M = K
      KM1 = K-1
C                                  FILL IPER WITH INTEGERS 1 TO K.
      DO 5 I=1,K
         IPER(I) = I
    5 CONTINUE
C                                  NOW RANDOMLY PERMUTE IPER.
      DO 10 I=1,KM1
         J = 1+GGUBFS(DSEED)*M
         ITEMP = IPER(M)
         IPER(M) = IPER(J)
         IPER(J) = ITEMP
         M = M-1
   10 CONTINUE
      RETURN
      END
