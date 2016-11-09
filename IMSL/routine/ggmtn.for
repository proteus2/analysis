C   IMSL ROUTINE NAME   - GGMTN
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1980
C
C   PURPOSE             - MULTINOMIAL RANDOM DEVIATE GENERATOR
C
C   USAGE               - CALL GGMTN (DSEED,NR,NIND,K,P,IIR,IR)
C
C   ARGUMENTS    DSEED  - INPUT/OUTPUT DOUBLE PRECISION VARIABLE
C                           ASSIGNED AN INTEGER VALUE IN THE
C                           EXCLUSIVE RANGE (1.D0, 2147483647.D0).
C                           DSEED IS REPLACED BY A NEW VALUE TO BE
C                           USED IN A SUBSEQUENT CALL.
C                NR     - INPUT NUMBER OF K-DEVIATE VECTORS TO BE
C                           GENERATED.
C                NIND   - INPUT MULTINOMIAL PARAMETER INDICATING THE
C                           NUMBER OF INDEPENDENT TRIALS.
C                K      - INPUT PARAMETER INDICATING THE NUMBER OF
C                           MUTUALLY EXCLUSIVE OUTCOMES ON ANY TRIAL.
C                P      - INPUT VECTOR OF LENGTH K CONTAININDG, IN P(I),
C                           THE PROBABILITY THAT EVENT I WILL OCCUR
C                           WHERE I = 1,2,...,K.   THE ELEMENTS OF
C                           VECTOR P MUST BE POSITIVE AND MUST SUM TO
C                           ONE.
C                IIR    - INPUT ROW DIMENSION OF MATRIX IR EXACTLY AS
C                           SPECIFIED IN THE CALLING PROGRAM.
C                IR     - OUTPUT NR BY K MATRIX OF MULTINOMIAL
C                           DEVIATES.
C
C   REQD. IMSL ROUTINES - GGBN,GGBTR,GGUBFS,GGUBS,UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE GGMTN   (DSEED,NR,NIND,K,P,IIR,IR)
C
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            NIND,NR,IIR,K,IR(IIR,K)
      REAL               P(K)
      DOUBLE PRECISION   DSEED
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
C
      INTEGER            I,IB(1),IBNRV,J,NLEFT,NUSED
      REAL               PR,PTOTL,ONE
C                                  DATA INITIALIZATION
      DATA               ONE/1.0/
C                                  FIRST EXECUTABLE STATEMENT
      KM1 = K-1
      DO 10 L=1,NR
         NUSED = 0
         PTOTL = ONE
         DO 5 I=1,KM1
            PR = P(I)/PTOTL
            PTOTL = PTOTL-P(I)
            NLEFT = NIND-NUSED
            CALL GGBN(DSEED,1,NLEFT,PR,IB)
            IBNRV = IB(1)
            NUSED = NUSED+IBNRV
            IR(L,I) = IBNRV
    5    CONTINUE
         IR(L,K) = NIND-NUSED
   10 CONTINUE
      RETURN
      END
