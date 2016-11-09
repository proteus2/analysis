C   IMSL ROUTINE NAME   - GGNLG
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1980
C
C   PURPOSE             - LOG-NORMAL RANDOM DEVIATE GENERATOR
C
C   USAGE               - CALL GGNLG (DSEED,NR,XM,S,R)
C
C   ARGUMENTS    DSEED  - INPUT/OUTPUT DOUBLE PRECISION VARIABLE
C                           ASSIGNED AN INTEGER VALUE IN THE
C                           EXCLUSIVE RANGE (1.D0, 2147483647.D0).
C                           DSEED IS REPLACED BY A NEW VALUE TO BE
C                           USED IN A SUBSEQUENT CALL.
C                NR     - INPUT NUMBER OF DEVIATES TO BE GENERATED.
C                XM     - INPUT PARAMETER. SEE REMARKS.
C                S      - INPUT PARAMETER. SEE REMARKS.
C                R      - OUTPUT VECTOR OF LENGTH NR+1 CONTAINING THE
C                           LOG NORMAL DEVIATES IN THE FIRST NR
C                           LOCATIONS. R(NR+1) IS WORK STORAGE.
C
C   REQD. IMSL ROUTINES - GGNPM,GGUBS
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS      NORMAL (0,1) PSEUDO RANDOM DEVIATES U ARE TRANSFORMED
C                TO LOG-NORMAL DEVIATES X, AS X = EXP(XM + S*U).
C
C-----------------------------------------------------------------------
C
      SUBROUTINE GGNLG  (DSEED,NR,XM,S,R)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            NR
      REAL               R(1),XM,S
      DOUBLE PRECISION   DSEED
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,NN,M
C                                  FIRST EXECUTABLE STATEMENT
      NN = NR
      M = MOD(NR,2)
      IF(M.NE.0) NN = NN+1
      CALL GGNPM(DSEED,NN,R)
      DO 5 I = 1,NR
         R(I) = EXP(XM+S*R(I))
    5 CONTINUE
      RETURN
      END
