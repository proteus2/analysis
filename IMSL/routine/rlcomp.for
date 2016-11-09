C   IMSL ROUTINE NAME   - RLCOMP
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - GENERATION OF AN ORTHOGONAL CENTRAL
C                           COMPOSITE DESIGN
C
C   USAGE               - CALL RLCOMP (N,XMNX,IX,IOP,M,DSN,IS,IER)
C
C   ARGUMENTS    N      - INPUT VECTOR OF LENGTH TWO.
C                         N(1) CONTAINS THE NUMBER OF VARIABLES.
C                           N(1) MUST BE GREATER THAN OR EQUAL TO TWO
C                           AND LESS THAN OR EQUAL TO EIGHT.
C                         N(2) CONTAINS THE NUMBER OF CENTER POINTS.
C                           N(2) MUST BE GREATER THAN ZERO.
C                XMNX   - INPUT N(1) BY 2 MATRIX CONTAINING THE MINIMUM
C                           AND MAXIMUM VALUE FOR EACH OF THE N(1)
C                           VARIABLES. THE MATRIX XMNX MUST CONTAIN
C                           THE MINIMUM FOR VARIABLE I IN ROW I, COLUMN
C                           1, AND THE MAXIMUM FOR VARIABLE I IN ROW I,
C                           COLUMN 2, FOR I=1,2,...,N(1).
C                IX     - INPUT ROW DIMENSION OF THE MATRIX XNMN EXACTLY
C                           AS SPECIFIED IN THE DIMENSION STATEMENT IN
C                           THE CALLING PROGRAM.
C                IOP    - INPUT OPTION, REQUIRED ONLY WHEN N(1) = 5,6,7,
C                           OR 8, INDICATING WHICH FRACTIONAL REPLICA-
C                           TION OF THE 2**N(1) PLAN IS DESIRED. SEE
C                           THE SECOND TABLE IN THE ALGORITHM SECTION
C                           IN THE MANUAL DOCUMENT FOR FURTHER DETAILS.
C                M      - OUTPUT TOTAL NUMBER OF DESIGN POINTS.
C                DSN    - OUTPUT M BY N(1) MATRIX CONTAINING THE
C                           ORTHOGONAL CENTRAL COMPOSITE DESIGN.
C                           THE DESIGN SETTINGS FOR VARIABLE I ARE
C                           CONTAINED IN COLUMN I OF THE MATRIX DSN,
C                           FOR I=1,2,...,N(1).
C                IS     - INPUT ROW DIMENSION OF THE MATRIX DSN EXACTLY
C                           AS SPECIFIED IN THE DIMENSION STATEMENT IN
C                           THE CALLING PROGRAM.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER=129 MEANS THAT N(1) OR N(2) WAS OUT OF
C                             RANGE
C                           IER=130 MEANS THAT XMNX(I,1) IS GREATER
C                             THAN OR EQUAL TO XMNX(I,2) FOR SOME I,
C                             I=1,2,...,N(1).
C
C   REQD. IMSL ROUTINES - UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE RLCOMP (N,XMNX,IX,IOP,M,DSN,IS,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N(1),IX,IOP,M,IS,IER
      REAL               XMNX(IX,1),DSN(IS,1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            IND(8),ISUM(8),ICY(8),I,IC,ICT,ID,IFF,II,
     1                   IJ1,IJ2,IN,IN1,IN2,J,JT,J1,J2,L,N1,N1I,N2
      REAL               F,XM,Q,ALFA,V1,V2,S,ONE,TWO,FOUR,FOURTH
      DATA               ONE,TWO,FOUR,FOURTH/1.0,2.0,4.0,0.25/
      DATA               ICY/1,2,4,8,16,32,64,128/
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
      N1 = N(1)
      N2 = N(2)
C                                  TERMINAL ERROR CHECK FOR N(2)
      IF (N2 .LT. 1) GO TO 130
C                                  TERMINAL ERROR CHECK FOR N(1)
      IF (N1 .LT. 2 .OR. N1 .GT. 8) GO TO 130
C                                  TERMINAL ERROR CHECK FOR XMNX
      DO 5 I = 1,N1
         IF (XMNX(I,1) .LT. XMNX(I,2)) GO TO 5
         IER = 130
         GO TO 9000
    5 CONTINUE
C                                  COMPUTE NUMBER OF FACTORIAL POINTS
      IFF = 2**N1
      IF (N1 .LE. 4) GO TO 10
      IFF = IFF/2
      IF (N1 .EQ. 8) IFF = IFF/2
C                                  COMPUTE TOTAL NUMBER OF DESIGN POINTS
   10 M = 2*N1 + N2 + IFF
      XM = M
      F = IFF
      Q = (SQRT(XM) - SQRT(F))**TWO
      ALFA = ((Q*F)/FOUR)**FOURTH
C                                  SET INDICATORS FOR FRACTIONAL
C                                  REPLICATE DESIGNS
      IF (N1 .LE. 4) GO TO 35
      IF (N1 .EQ. 8) GO TO 15
      ID = 0
      IF (IOP .NE. 0) ID = 1
      GO TO 35
   15 IJ1 = 0
      IF (IOP .NE. 0) GO TO 20
      IJ2 = 0
      GO TO 35
   20 IF (IOP .NE. 1) GO TO 25
      IJ2 = 1
      GO TO 35
   25 IJ1 = 1
      IF (IOP .NE. 2) GO TO 30
      IJ2 = 0
      GO TO 35
   30 IJ2 = 1
C                                  COMPUTE SETTINGS FOR EACH VARIABLE
   35 DO 125 I = 1,N1
         V1 = (XMNX(I,2) - XMNX(I,1))/(TWO*ALFA)
         V2 = (XMNX(I,2) + XMNX(I,1))/TWO
         S = -ONE
         N1I = N1 - I + 1
         L = ICY(N1I)
         J = 0
         ICT = 0
         IF (N1 .LE. 4) GO TO 45
         DO 40 II = 1,N1
            ISUM(II) = 1
            IND(II) = 1
   40    CONTINUE
   45    J = J + 1
         IF (J - ICT .GT. IFF) GO TO 110
         IC = 0
         S = -S
         GO TO 55
   50    J = J + 1
   55    JT = J
         IF (N1 .LE. 4) GO TO 100
         IF (N1 .EQ. 8) GO TO 65
         IN = 0
         DO 60 II = 1,N1
            IN = IN + IND(II)
   60    CONTINUE
         IN = MOD(IN,2)
         GO TO 80
   65    IN1 = 0
         IN2 = 0
         DO 70 II = 1,5
            IN1 = IN1 + IND(II)
   70    CONTINUE
         DO 75 II = 4,8
            IN2 = IN2 + IND(II)
   75    CONTINUE
         IN1 = MOD(IN1,2)
         IN2 = MOD(IN2,2)
   80    DO 85 II = 1,N1
            ISUM(II) = ISUM(II) + 1
            N1I = N1 - II + 1
            IF (ISUM(II) .LE. ICY(N1I)) GO TO 85
            IND(II) = IND(II) + 1
            IND(II) = MOD(IND(II),2)
            ISUM(II) = 1
   85    CONTINUE
         IF (N1 .EQ. 8) GO TO 90
         IF (IN .EQ. ID) GO TO 95
         ICT = ICT + 1
         GO TO 105
   90    IF (IN1 .EQ. IJ1 .AND. IN2 .EQ. IJ2) GO TO 95
         ICT = ICT + 1
         GO TO 105
   95    J = J - ICT
  100    DSN(J,I) = V2 + SIGN(V1,S)
         J = JT
  105    IC = IC + 1
         IF (IC .EQ. L) GO TO 45
         GO TO 50
  110    J1 = IFF + 1
         J2 = IFF + N2
         DO 115 J = J1,J2
            DSN(J,I) = V2
  115    CONTINUE
         J1 = J2 + 1
         J2 = J2 + 2*N1
         DO 120 J = J1,J2
            DSN(J,I) = V2
  120    CONTINUE
         J = J1 + 2*I - 2
         DSN(J,I) = XMNX(I,1)
         J = J + 1
         DSN(J,I) = XMNX(I,2)
  125 CONTINUE
      GO TO 9005
  130 IER = 129
 9000 CALL UERTST(IER,6HRLCOMP)
 9005 RETURN
      END
