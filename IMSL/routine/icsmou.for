C   IMSL ROUTINE NAME   - ICSMOU
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - ONE-DIMENSIONAL DATA SMOOTHING BY ERROR
C                           DETECTION
C
C   USAGE               - CALL ICSMOU (X,Y,NX,DIS,SC,MAXIT,WK,IER)
C
C   ARGUMENTS    X      - VECTOR OF LENGTH NX CONTAINING THE ABSCISSAE
C                           OF THE NX DATA POINTS (X(I),Y(I)),I=1,...,
C                           NX. (INPUT) X MUST BE ORDERED SO THAT
C                           X(I) .LT. X(I+1).
C                Y      - VECTOR OF LENGTH NX CONTAINING THE ORDINATES
C                           (OR FUNCTION VALUES) OF THE NX DATA POINTS.
C                           (INPUT) ON OUTPUT Y CONTAINS THE SMOOTHED
C                           ORDINATES. NOTE - FOR PURPOSES OF
C                           NUMERICAL STABILITY, THE ROUTINE ASSUMES
C                           THAT Y(I),I=1,2,3,NX-2,NX-1,NX ARE NOT IN
C                           ERROR. THE ROUTINE NEVER MOVES THESE POINTS.
C                NX     - NUMBER OF ELEMENTS IN X AND Y. (INPUT)
C                           NX MUST BE .GE. 7.
C                DIS    - PROPORTION OF THE DISTANCE THE ORDINATE IN
C                           ERROR IS MOVED TO ITS INTERPOLATING CURVE.
C                           (INPUT) 0.0 .LE. DIS .LE. 1.0. A SUGGESTED
C                           VALUE IS DIS=1.0. SEE PROGRAMMING NOTES FOR
C                           FURTHER DETAILS.
C                SC     - STOPPING CRITERION. (INPUT)
C                           IF THE ABSOLUTE VALUE OF THE DISTANCE FROM
C                           Y(I) TO ITS INTERPOLATING CURVE IS LESS
C                           THAN OR EQUAL TO SC*((X(I+3)-X(I-3))/6.0)
C                           THEN Y(I) IS NOT MOVED, FOR I=4,...,NX-3.
C                           A SUGGESTED VALUE IS SC=0.0. SEE
C                           PROGRAMMING NOTES FOR FURTHER DETAILS.
C                MAXIT  - MAXIMUM NUMBER OF ITERATIONS ALLOWED. (INPUT)
C                           ON OUTPUT MAXIT CONTAINS THE ACTUAL NUMBER
C                           OF ITERATIONS TAKEN. A SUGGESTED VALUE IS
C                           MAXIT=ESTIMATED NUMBER OF POINTS IN ERROR.
C                           SEE PROGRAMMING NOTES FOR FURTHER DETAILS.
C                WK     - WORK AREA. WK IS A VECTOR OF DIMENSION
C                           2*NX-12.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER = 129, NX IS LESS THAN 7.
C                           IER = 130, INPUT ABSCISSAE ARE NOT ORDERED
C                             SO THAT X(1) .LT. X(2) ... .LT. X(NX).
C
C   REQD. IMSL ROUTINES - ICSEVU,IQHSCU,UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS      THIS ROUTINE WAS DESIGNED TO SMOOTH A DATA SET WHICH
C                IS MILDLY CONTAMINATED WITH ISOLATED ERRORS. IN
C                GENERAL, THE ROUTINE WILL NOT WORK WELL IF MORE THAN
C                25 PERCENT OF THE DATA POINTS ARE IN ERROR.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE ICSMOU (X,Y,NX,DIS,SC,MAXIT,WK,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            NX,MAXIT,IER
      REAL               X(NX),Y(NX),DIS,SC,WK(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,NXM3,NXM6,NXM9,IT,ISTART,ISTOP,II,J,IIER
      REAL               C(6,3),XTEMP(6),YTEMP(6),T(1),
     1                   D,DD,ZERO,SIX
      DATA               ZERO/0.0/,SIX/6.0/
C                                  FIRST EXECUTABLE STATEMENT
      IF (NX .GT. 6) GO TO 5
C                                  TERMINAL - TOO FEW DATA POINTS
      IER = 129
      GO TO 9000
    5 DO 10 I=2,NX
         IF (X(I) .GT. X(I-1)) GO TO 10
C                                  TERMINAL - X NOT MONOTONE INCREASING
         IER = 130
         GO TO 9000
   10 CONTINUE
      IER = 0
      NXM3 = NX-3
      NXM6 = NX-6
      NXM9 = NX-9
      IT = MAXIT
C                                  CALCULATE STOPPING CRITERION.
      D = SC/SIX
      DO 15 I=1,NXM6
         WK(I) = D*(X(I+6)-X(I))
   15 CONTINUE
      ISTART = 4
      ISTOP = NXM3
      MAXIT = 0
C                                  EACH ITERATION BEGINS HERE.
C                                  CALCULATE POINT ENERGY FOR AFFECTED
C                                    POINTS.
   20 DO 30 I=ISTART,ISTOP
         II = -3
C                                  SET UP TEMPORARY VECTORS FOR INPUT
C                                    TO IQHSCU.
         DO 25 J=1,6
            XTEMP(J) = X(I+II)
            YTEMP(J) = Y(I+II)
            II = II+1
            IF (II .EQ. 0) II = 1
   25    CONTINUE
C                                  CALL IQHSCU TO GET THE COEFFICIENTS.
         CALL IQHSCU(XTEMP,YTEMP,6,C,6,IIER)
C                                  CALL ICSEVU TO EVALUATE AT X(I).
         CALL ICSEVU(XTEMP,YTEMP,6,C,6,X(I),T,1,IIER)
C                                  CALCULATE THE POINT ENERGY OF Y(I).
         WK(NXM9+I) = Y(I)-T(1)
   30 CONTINUE
C                                  CHECK TO SEE IF THE MAXIMUM NUMBER
C                                    OF ITERATIONS HAS BEEN EXCEEDED.
C                                    IF SO RETURN.
      IF (MAXIT .GE. IT) GO TO 9005
C                                  FIND THE LARGEST DEVIATION.
      II = 0
      D = ZERO
      DO 35 I=1,NXM6
         DD = ABS(WK(NXM6+I))
         IF (DD .LE. WK(I)) GO TO 35
         IF (DD .LE. D) GO TO 35
         D = DD
         II = I+3
   35 CONTINUE
C                                  CHECK TO SEE IF STOPPING CRITERION
C                                    HAS BEEN MET. IF SO RETURN.
      IF (II .LE. 0) GO TO 9005
C                                  STOPPING CRITERION HAS NOT BEEN
C                                    SATISFIED. MOVE Y(II)
C                                    DIS*POINT ENERGY.
      MAXIT = MAXIT+1
      Y(II) = Y(II)-(WK(NXM9+II)*DIS)
C                                  RECOMPUTE DEVIATIONS AFFECTED BY
C                                    MOVEMENT OF Y(II).
      ISTART = MAX0(II-3,4)
      ISTOP = MIN0(NXM3,II+3)
      GO TO 20
C                                  RETURN SECTION
 9000 CONTINUE
      CALL UERTST(IER,'ICSMOU')
 9005 RETURN
      END
