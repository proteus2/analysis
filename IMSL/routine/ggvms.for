C   IMSL ROUTINE NAME   - GGVMS
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1982
C
C   PURPOSE             - VON MISES RANDOM DEVIATE GENERATOR
C
C   USAGE               - CALL GGVMS (DSEED,C,NR,R)
C
C   ARGUMENTS    DSEED  - INPUT/OUTPUT.  DOUBLE PRECISION VARIABLE
C                           ASSIGNED AN INTEGER VALUE IN THE
C                           EXCLUSIVE RANGE (1.D0, 2147483647.D0).
C                           DSEED IS REPLACED BY A NEW VALUE TO BE
C                           USED IN A SUBSEQUENT CALL.
C                C      - INPUT.  PARAMETER FOR THE VON MISES DISTRI-
C                           BUTION.  C MUST BE POSITIVE.
C                NR     - INPUT.  NUMBER OF DEVIATES TO BE GENERATED.
C                R      - OUTPUT.  VECTOR OF LENGTH NR CONTAINING THE
C                           VON MISES DEVIATES.
C
C   REQD. IMSL ROUTINES - GGUBFS
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE GGVMS (DSEED,C,NR,R)
C
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            NR
      REAL               C,R(NR)
      DOUBLE PRECISION   DSEED
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I
      REAL               CALG,CC,CMUD,D,E,F,FF,FOUR,ONE,PORN,PT25,PT5,
     *                   RR,S,T,TWO,U,V,Z,ZERO
C                                  SPECIFICATIONS FOR FUNCTIONS
      REAL               GGUBFS
      DATA               ZERO /0.0/,PT25 /0.25/,PT5 /0.5/
      DATA               ONE /1.0/,TWO /2.0/,FOUR /4.0/
C                                  FIRST EXECUTABLE STATEMENT
      RR = ONE+SQRT(ONE+FOUR*C*C)
      RR = (RR-SQRT(TWO*RR))/(TWO*C)
      RR = (ONE+RR*RR)/(TWO*RR)
      DO 15 I=1,NR
    5    V = GGUBFS(DSEED)-PT5
         D = V**2
         E = (GGUBFS(DSEED)-PT5)**2
         S = D+E
         IF (S.GT.PT25) GO TO 5
         T = D/E
         Z = (ONE-T)/(ONE+T)
         F = (ONE+RR*Z)/(RR+Z)
         CC = C*(RR-F)
         U = GGUBFS(DSEED)
         CMUD = CC*(TWO-CC)-U
         IF (CMUD.GT.ZERO) GO TO 10
         CALG = ALOG(CC/U)+ONE
         IF (CALG.LT.CC) GO TO 5
   10    FF = SQRT(ONE-F*F)
         R(I) = ATAN2(FF,F)
         IF (V.LT.ZERO) R(I) = -R(I)
   15 CONTINUE
      RETURN
      END
