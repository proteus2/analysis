C   IMSL ROUTINE NAME   - AGLMOD
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1982
C
C   PURPOSE             - GENERAL LINEAR MODEL ANALYSIS
C
C   USAGE               - CALL AGLMOD (X,IX,NV,Y,XL,IL,BETA,SS,VARB,WK,
C                           IPR,IER)
C
C   ARGUMENTS    X      - INPUT NV(1) BY NV(2) MATRIX CONTAINING THE
C                           DESIGN (INDEPENDENT VARIABLE SETTINGS).
C                IX     - INPUT ROW DIMENSION OF THE MATRIX X EXACTLY
C                           AS SPECIFIED IN THE DIMENSION STATEMENT IN
C                           THE CALLING PROGRAM.
C                NV     - INPUT VECTOR OF LENGTH 3.
C                           NV(1) CONTAINS THE NUMBER OF DATA POINTS
C                             (ROWS OF X). NV(1) MUST BE GREATER THAN
C                             OR EQUAL TO ONE.
C                           NV(2) CONTAINS THE NUMBER OF PARAMETERS
C                             (COLUMNS OF X). NV(2) MUST BE GREATER
C                             THAN OR EQUAL TO ONE AND LESS THAN OR
C                             EQUAL TO NV(1).
C                           NV(3) CONTAINS THE NUMBER OF INDEPENDENT
C                             LINEAR RESTRICTIONS ON THE PARAMETERS.
C                             NV(3) MUST BE GREATER THAN OR EQUAL TO
C                             ZERO AND LESS THAN OR EQUAL TO NV(2).
C                Y      - INPUT VECTOR OF LENGTH NV(1) CONTAINING THE
C                           RESPONSES.
C                XL     - INPUT NV(3) BY NV(2) MATRIX CONTAINING THE
C                           COEFFICIENTS FOR THE LINEAR RESTRICTIONS
C                           ON THE PARAMETERS.
C                IL     - INPUT ROW DIMENSION OF THE MATRIX XL EXACTLY
C                           AS SPECIFIED IN THE DIMENSION STATEMENT IN
C                           THE CALLING PROGRAM.
C                BETA   - OUTPUT VECTOR OF LENGTH NV(2) CONTAINING THE
C                           PARAMETER ESTIMATES.
C                SS     - OUTPUT SCALAR CONTAINING THE SUM OF SQUARES
C                           ATTRIBUTABLE TO THE PARAMETERS IN THE
C                           MODEL.
C                VARB   - OUTPUT VECTOR OF LENGTH NV(2)*(NV(2)+1)/2
C                           CONTAINING THE NV(2) BY NV(2) VARIANCE-
C                           COVARIANCE MATRIX ESTIMATE DIVIDED BY THE
C                           ERROR MEAN SQUARE FOR THE PARAMETER
C                           ESTIMATES. VARB IS STORED IN SYMMETRIC
C                           STORAGE MODE.
C                WK     - WORK AREA OF LENGTH 3*(IPR**2) (SEE DES-
C                           CRIPTION OF PARAMETER IPR, BELOW).
C                IPR    - INPUT SCALAR CONTAINING THE VALUE OF
C                           NV(2)+NV(3). IPR IS USED TO PROVIDE PROPER
C                           DIMENSIONING OF THE WORK AREA, WK.
C                           IPR MUST BE GREATER THAN 2.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER=129 INDICATES THAT EITHER THE NUMBER
C                             OF DATA POINTS, NV(1), OR THE NUMBER OF
C                             PARAMETERS, NV(2), WAS LESS THAN ONE OR
C                             THAT THE NUMBER OF LINEAR RESTRICTIONS,
C                             NV(3), WAS LESS THAN ZERO.
C                           IER=130 INDICATES THAT EITHER THE NUMBER OF
C                             PARAMETERS, NV(2), WAS GREATER THAN THE
C                             NUMBER OF DATA POINTS, NV(1), OR THAT
C                             THE NUMBER OF LINEAR RESTRICTIONS,
C                             NV(3), WAS GREATER THAN THE NUMBER OF
C                             PARAMETERS.
C                           IER=131 INDICATES MATRIX SINGULARITY
C                             DISCOVERED BY IMSL ROUTINE LINV3F.
C
C   REQD. IMSL ROUTINES - SINGLE/LINV3F,LUDATN,LUELMN,UERTST,UGETIO,
C                           VCVTFS,VCVTSF,VIPRFF,VMULFF,VMULFS,VTPROF
C                       - DOUBLE/LINV3F,LUDATN,LUELMN,UERTST,UGETIO,
C                           VCVTFS,VCVTSF,VIPRFF,VMULFF,VMULFS,VTPROF,
C                           VXADD,VXMUL,VXSTO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE AGLMOD (X,IX,NV,Y,XL,IL,BETA,SS,VARB,WK,IPR,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            IX,NV(3),IL,IPR,IER
      REAL               X(IX,1),Y(1),XL(IL,1),BETA(1),SS,VARB(1),
     1                   WK(IPR,1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,IP,IPR2,IR,I1,J,N
      REAL               EPS,D1,D2,DUMY(1)
      DOUBLE PRECISION   SUM
      DATA               EPS/1.E-20/
      DATA               D1/-1.E0/
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
      N  = NV(1)
      IP = NV(2)
      IR = NV(3)
      IPR2 = IPR + IPR
      IF (N .GE. 1 .AND. IP .GE. 1 .AND. IR .GE. 0) GO TO 5
C                                  TERMINAL - INVALID NUMBER OF INPUT
C                                             DATA POINTS SPECIFIED
      IER = 129
      GO TO 9000
   5  IF (IP .LE. N .AND. IR .LE. IP) GO TO 10
C                                  TERMINAL - INVALID RELATIONSHIP
C                                             AMONG THE NUMBER OF
C                                             DATA POINTS SPECIFIED IN
C                                             NV(1),NV(2),NV(3)
      IER = 130
      GO TO 9000
C                                  GET XTX AND STORE IN Z
  10  CALL VTPROF(X,N,IP,IX,WK(1,IPR2+1))
      CALL VCVTSF(WK(1,IPR2+1),IP,WK,IPR)
C                                  MOVE XL AND XL TRANSPOSE INTO Z
      IF (IR .LE. 0) GO TO 35
      DO 20 I = 1,IR
         I1 = IP + I
         DO 15 J = 1,IP
            WK(I1,J) = XL(I,J)
            WK(J,I1) = XL(I,J)
  15     CONTINUE
  20  CONTINUE
C                                  FILL IN THE REST OF Z WITH ZEROS
      DO 30 I = 1,IR
         DO 25 J = 1,IR
            WK(IP+I,IP+J) = 0.0
  25     CONTINUE
  30  CONTINUE
C                                           -1
C                                  NOW GET Z
  35  CALL LINV3F(WK,DUMY,1,IPR,IPR,D1,D2,WK(1,IPR+1),IER)
C                                  TERMINAL - MATRIX SINGULAR
      IF(IER .LT. 128) GO TO 37
      IER = 131
      GO TO 9000
C                                  CHECK FOR POSSIBLE UNDERFLOW
  37  DO 45 I = 1,IP
         DO 40 J = 1,I
            IF ( ABS(WK(I,J)) .GT. EPS) GO TO 40
            WK(I,J) = 0.0
            WK(J,I) = 0.0
  40     CONTINUE
  45  CONTINUE
C                                          -1          -1
C                                  VARB = Z   * XTX * Z
C                                          11          11
C
C                                        -1
C                                  GET  Z   * XTX  FIRST
C                                        11
      CALL VMULFS(WK,WK(1,IPR2+1),IP,IP,IPR,WK(1,IPR+1),IPR)
C                                  NOW GET VARB IN FULL STORAGE MODE
      CALL VMULFF(WK(1,IPR+1),WK,IP,IP,IP,IPR,IPR,WK(1,IPR2+1),IPR,IER)
C                                  GET VARB IN SYMMETRIC MODE
      CALL VCVTFS(WK(1,IPR2+1),IP,IPR,VARB)
C                                  COMPUTE XTY
      DO 55 I = 1,IP
         SUM = 0.D0
         DO 50 J = 1,N
            SUM = SUM + DBLE(X(J,I))*DBLE(Y(J))
  50     CONTINUE
         WK(I,IPR+1) = SUM
  55  CONTINUE
C                                          -1
C                                  BETA = Z   * XTY
C                                          11
      CALL VMULFF(WK,WK(1,IPR+1),IP,IP,1,IPR,IPR,BETA,IP,IER)
C
C                                  SS = BETA * XTY
      SUM = 0.0D0
      DO 65 I = 1,IP
         SUM = SUM + DBLE(BETA(I))*DBLE(WK(I,IPR+1))
  65  CONTINUE
      SS = SUM
 9000 CONTINUE
      IF (IER .NE. 0) CALL UERTST(IER,'AGLMOD')
 9005 RETURN
      END
