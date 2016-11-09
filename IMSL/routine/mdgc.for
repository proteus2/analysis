C   IMSL ROUTINE NAME   - MDGC
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1980
C
C   PURPOSE             - GENERAL CUMULATIVE PROBABILITY DISTRIBUTION
C                           FUNCTION, GIVEN ORDINATES OF THE DENSITY
C
C   USAGE               - CALL MDGC (X,F,M,IOPT,B,C,P,IER)
C
C   ARGUMENTS    X      - INPUT VALUE AT WHICH THE FUNCTION IS TO BE
C                           EVALUATED.
C                F      - INPUT ARRAY OF LENGTH M CONTAINING THE
C                           PROBABILITY DENSITY FUNCTION ORDINATES.
C                M      - INPUT NUMBER OF ORDINATES SUPPLIED. M MUST BE
C                           GREATER THAN 1 IF LINEAR INTERPOLATION IS
C                           DESIRED AND GREATER THAN 3 IF A CURVE
C                           FITTED THROUGH THE ORDINATES IS DESIRED.
C                IOPT   - INPUT OPTION PARAMETER.
C                           IOPT = 1 IMPLIES LINEAR INTERPOLATION IS
C                             DESIRED AND THE DATA ARE EQUALLY SPACED.
C                           IOPT = 2 IMPLIES LINEAR INTERPOLATION IS
C                             DESIRED AND THE DATA ARE NOT NECESSARILY
C                             EQUALLY SPACED.
C                           IOPT = 3 IMPLIES A CURVE FITTED THROUGH THE
C                             DATA IS DESIRED AND THE DATA ARE EQUALLY
C                             SPACED.
C                           IOPT = 4 IMPLIES A CURVE FITTED THROUGH THE
C                             DATA IS DESIRED AND THE DATA ARE NOT
C                             NECESSARILY EQUALLY SPACED.
C                B      - INPUT ARRAY OF LENGTH M IF IOPT = 2, 3, OR 4.
C                           OTHERWISE, B IS OF LENGTH 2. IF IOPT = 2 OR
C                           4, B(I) IS THE ABSCISSA CORRESPONDING TO
C                           F(I). THE ABSCISSAE MUST BE SPECIFIED IN
C                           INCREASING ORDER. IF IOPT = 1 OR 3, B(1)
C                           IS THE LOWER ENDPOINT OF THE SUPPORT OF THE
C                           DISTRIBUTION (CORRESPONDING TO F(1)) AND
C                           B(2) IS THE UPPER ENDPOINT (CORRESPONDING
C                           TO F(M)).
C                C      - WORK AREA. IF IOPT = 3 OR 4, THEN C MUST BE
C                           OF LENGTH AT LEAST 3*M. OTHERWISE, ONLY
C                           C(1) IS USED.
C                P      - OUTPUT PROBABILITY THAT THE RANDOM VARIABLE
C                           IS LESS THAN OR EQUAL TO X.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER = 129 INDICATES AN ERROR OCCURRED IN
C                             IMSL SUBROUTINE IQHSCU.
C                           IER = 130 INDICATES AN ELEMENT OF F WAS
C                             NEGATIVE.
C                           IER = 131 INDICATES THE ELEMENTS OF B ARE
C                             NOT IN INCREASING ORDER.
C                           IER = 132 INDICATES M IS OUT OF RANGE.
C                           IER = 133 INDICATES IOPT IS OUT OF RANGE.
C
C   REQD. IMSL ROUTINES - DCSQDU,IQHSCU,MDGD,UERSET,UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE MDGC   (X,F,M,IOPT,B,C,P,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            IER,IOPT,M
      REAL               B(M),C(1),F(M),P,X
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,LEVEL,LEVOLD,L
      REAL               P1
C                                  FIRST EXECUTABLE STATEMENT
      P = 1.0
      IER = 0
      IF (M.GT.1) GO TO 5
      IER = 132
      GO TO 9000
    5 IF (IOPT.GT.0.AND.IOPT.LT.5) GO TO 10
      IER = 133
      GO TO 9000
   10 L = 2
      IF (IOPT.EQ.2.OR.IOPT.EQ.4) L = M
      IF (B(L).GE.B(1)) GO TO 15
      IER = 131
      GO TO 9000
   15 IF (X.GE.B(L)) GO TO 9005
      P = 0.0
      IF (X.LE.B(1)) GO TO 9005
      DO 20 I=1,M
         IF (F(I).GE.0.0) GO TO 20
         IER = 130
         GO TO 9000
   20 CONTINUE
      GO TO (25,25,35,35), IOPT
C                                  LINEAR CASE
   25 IF (M.GT.1) GO TO 30
      IER = 132
      GO TO 9000
   30 C(1) = X
      CALL MDGD (F,M,IOPT,B,C,P,IER)
      IF (IER.NE.0) GO TO 9000
      C(1) = B(L)
      CALL MDGD (F,M,IOPT,B,C,P1,IER)
      IF (IER.NE.0) GO TO 9000
      P = P/P1
      GO TO 45
C                                  FITTED CURVE CASE
   35 IF (M.GT.3) GO TO 40
      IER = 132
      GO TO 9000
   40 C(1) = X
      CALL MDGD (F,M,IOPT,B,C,P,IER)
      IF (IER.NE.0) GO TO 9000
      LEVEL = 0
      CALL UERSET (LEVEL,LEVOLD)
      CALL DCSQDU (B,F,M,C,M,X,B(M),P1,IER)
      IER = 0
      CALL UERSET (LEVOLD,LEVEL)
      P1 = P+P1
      P = P/P1
      IF (IOPT.EQ.3) B(2) = B(M)
   45 IF (P.GT.1.0E0) P = 1.0E0
      GO TO 9005
 9000 CONTINUE
      CALL UERTST (IER,6HMDGC  )
 9005 RETURN
      END
