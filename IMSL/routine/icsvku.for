C   IMSL ROUTINE NAME   - ICSVKU
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - LEAST SQUARES APPROXIMATION BY CUBIC
C                           SPLINES - VARIABLE KNOTS
C
C   USAGE               - CALL ICSVKU(X,F,NX,XK,NXK,Y,C,IC,ERROR,WK,IER)
C
C   ARGUMENTS    X      - VECTOR OF LENGTH NX CONTAINING THE ABSCISSA OF
C                           THE NX POINTS (X(I),F(I)) I=1,...,NX.
C                           X MUST BE ORDERED SO THAT X(1) .LE. X(2) ...
C                           .LE. X(NX). (INPUT)
C                F      - VECTOR OF LENGTH NX CONTAINING THE ORDINATES
C                           (OR FUNCTION VALUES) OF THE NX POINTS.
C                           (INPUT)
C                NX     - NUMBER OF ELEMENTS IN X AND F. (INPUT)
C                XK     - ON INPUT, VECTOR OF LENGTH NXK CONTAINING
C                           INITIAL GUESSES AT THE KNOT LOCATIONS. XK
C                           MUST BE ORDERED SO THAT XK(1) .LT. XK(2) ...
C                           .LT. XK(NXK). ALSO, XK(1) .LE. X(1) AND
C                           XK(NXK) .GE. X(NX) MUST HOLD. ON OUTPUT,
C                           XK CONTAINS THE KNOT LOCATIONS AS DETERMINED
C                           BY ICSVKU.
C                NXK    - NUMBER OF KNOTS (ELEMENTS OF XK). (INPUT)
C                Y,C    - SPLINE COEFFICIENTS. (OUTPUT) Y IS A VECTOR OF
C                           LENGTH NXK-1 AND C IS AN NXK-1 BY 3 MATRIX.
C                           THE VALUE OF THE SPLINE APPROXIMATION AT T
C                           IS S(T) = ((C(I,3)*D+C(I,2))*D+C(I,1))*D
C                           +Y(I) WHERE XK(I) .LE. T .LT. XK(I+1) AND
C                           D = T-XK(I).
C                IC     - ROW DIMENSION OF MATRIX C EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT
C                           IN THE CALLING PROGRAM. (INPUT)
C                ERROR  - LEAST SQUARES ERROR OF CUBIC SPLINE
C                           APPROXIMATION. (OUTPUT)
C                WK     - WORK AREA OF DIMENSION NX*(NXK+6).
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         WARNING ERROR
C                           IER = 33, ITERATION FAILED TO CONVERGE.
C                             KNOT LOCATIONS RETURNED IN XK ARE THE
C                             BEST THAT ICSVKU COULD DETERMINE, BUT THE
C                             CONVERGENCE CRITERIA WAS NOT SATISFIED.
C                         TERMINAL ERROR
C                           IER = 130, NXK IS GREATER THAN 28.
C                           IER = 131, INPUT ABSCISSAE ARE NOT ORDERED
C                             SO THAT X(1) .LE. X(2) ... .LE. X(NX).
C                           IER = 132, INPUT ABSCISSAE DO NOT LIE
C                             BETWEEN THE BOUNDARY KNOTS.
C                           IER = 133, KNOTS ARE NOT ORDERED SO THAT
C                             XK(1) .LT. XK(2) ... .LT. XK(NXK).
C
C   REQD. IMSL ROUTINES - ICSFKU,ICSFKV,UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS  1.  THE ERROR WHICH ICSFKU MINIMIZES IS DEFINED AS
C                ERROR = SQRT(R(1)**2*W(1)+...+R(NX)**2*W(NX)) WHERE
C                R(I) = F(I)-S(X(I)), I=1,...,NX
C                W(1) = (X(2)-X(1))/(X(NX)-X(1))
C                W(I) = (X(I+1)-X(I-1))/(X(NX)-X(1)), I=2,...,NX-1
C                W(NX) = (X(NX)-X(NX-1))/(X(NX)-X(1))
C                (X(I),F(I)), I=1,...,NX IS THE GIVEN SET OF POINTS
C                AND S IS THE LEAST SQUARES CUBIC SPLINE APPROXIMATION
C                TO THAT SET OF POINTS.
C            2.  ICSVKU IS PRIMARILY INTENDED FOR OBTAINING LOW TO
C                MEDIUM ACCURACY (2 TO 5 SIGNIFICANT DIGITS)
C                APPROXIMATIONS OF FUNCTIONS OF A MORE OR LESS ARBITRARY
C                NATURE. IT SHOULD NOT BE USED FOR HIGH ACCURACY
C                APPROXIMATIONS TO MATHEMATICALLY DEFINED FUNCTIONS.
C            3.  ICSVKU ATTEMPTS TO LOCATE THE KNOTS TO GIVE THE LEAST
C                SQUARES ERROR TO ONE SIGNIFICANT DIGIT.
C            4.  THE INITIAL KNOT LOCATIONS INPUT TO ICSVKU ARE
C                IMPORTANT TO DETERMINING THE LEAST SQUARES CUBIC
C                SPLINE APPROXIMATION. IT IS POSSIBLE FOR ICSVKU TO
C                REACH A LOCAL MINIMUM FAR FROM THE BEST APPROXIMATION.
C                IN SOME CASES, IT MIGHT BE WORTHWHILE TO EXPERIMENT
C                WITH SEVERAL INITIAL KNOT PLACEMENTS.
C            5.  COMMON BLOCK /ICSFK1/ APPEARS IN ICSVKU SO THAT IT
C                WILL ALWAYS BE PRESERVED BETWEEN CALLS TO ICSFKU.
C            6.  Y(NXK) CAN BE DEFINED BY THE FOLLOWING FORMULA
C                  Y(NXK) = ((C(NXK-1),3)*D+C(NXK-1,2))*D+C(NXK-1,1))*
C                           D+Y(NXK-1)
C                WHERE D = X(NXK)-X(NXK-1).
C                Y(NXK) CAN ALSO BE DEFINED (FOR MODE=0) BY CALLING
C                IMSL SUBROUTINE ICSEVU IN THE FOLLOWING MANNER
C                  CALL ICSEVU(XK,Y,NXK,C,IC,XK(NXK),Y(NXK),1,IER)
C
C-----------------------------------------------------------------------
C
      SUBROUTINE ICSVKU (X,F,NX,XK,NXK,Y,C,IC,ERROR,WK,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            NX,NXK,IC,IER
      REAL               X(NX),F(NX),XK(NXK),Y(IC),C(IC,3),ERROR,
     1                   WK(NX,1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,INDLP,ITER,IXX,J,KVARY,LPCNT,NXK1
      REAL               A,AA,ABEST,ACC,AEPS,AHIGH,ALEFT,ALOW,ARIGHT,
     1                   BD,CHANGE,DEL,DIFF,DIST,DXLEFT,DXRGHT,E,
     2                   EBEST,ELEFT,EPRED,EPSERA,EPSERR,ERIGHT,
     3                   ETRY,FMAX,H,HALF,PREVER,P4,Q,RNXK2,
     4                   XL,XR,ZERO,DXX,DYLEFT,DYRGHT,
     5                   VALEFT(1),VARIGH(1),VAA(1),VABEST(1)
      LOGICAL            MODE3
      COMMON /ICSFK1/    DXX(3809),IXX(69),MODE3
      EQUIVALENCE        (ALEFT,VALEFT(1)),(ARIGHT,VARIGH(1)),
     2                   (AA,VAA(1)),(ABEST,VABEST(1))
      DATA               INDLP/10/
      DATA               ZERO/0.0/,P4/0.4/,HALF/0.5/,
     1                   ACC/0.1/,AEPS/1.0E-3/,DIST/0.0625/
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
C                                  CHECK KNOT ORDERING
      DO 5 I=2,NXK
         IF (XK(I).GT.XK(I-1)) GO TO 5
         IER = 133
         GO TO 9000
    5 CONTINUE
C                                  INITIALIZE NUMERICAL CONTROL
C                                    VARIABLES
      RNXK2 = NXK-2
      ITER = 15
      EPSERR = P4*ACC
      CHANGE = P4*RNXK2
      FMAX = ZERO
      DO 10 I=1,NX
   10 FMAX = AMAX1(FMAX,ABS(F(I)))
      EPSERA = FMAX*AEPS
C                                  COMPUTE B.A. WITH INITIAL KNOTS
C
      CALL ICSFKU (X,F,NX,0,XK,NXK,Y,C,IC,ERROR,WK,IER)
      IF (IER.NE.0) GO TO 9000
      IF (NXK.LE.2) GO TO 9005
   15 CONTINUE
C                                  INITIALIZE FOR A NEW CYCLE
      KVARY = NXK-1
      XR = XK(KVARY+1)
      A = XK(KVARY)
      XL = XK(KVARY-1)
      PREVER = ERROR
      BD = ACC*ERROR/RNXK2
      Q = CHANGE/RNXK2
      CHANGE = ZERO
C                                  OPTIMIZE WRTO KNOT KVARY
   20 CONTINUE
C                                  FIND THREE POINTS ALEFT, A, AND
C                                    ARIGHT SO THAT
C                                    E(ALEFT) .GE. E(A) AND
C                                    E(ARIGHT) .GE. E(A)
C                                    WHERE E(Z) DENOTES THE L2 ERROR
C                                    WITH XK(KVARY) = Z
      H = XR-XL
      ALOW = XL+DIST*H
      AHIGH = XR-DIST*H
      LPCNT = 0
      E = ERROR
      ALEFT = A+Q*(XL-A)
      CALL ICSFKU (X,F,NX,3,VALEFT,1,Y,C,IC,ELEFT,WK,IER)
      IF (ELEFT.GE.E) GO TO 30
      GO TO 40
C                                  SEARCH FOR A NEW KNOT TO THE RIGHT
   25 ALEFT = A
      ELEFT = E
      A = ARIGHT
      E = ERIGHT
   30 ARIGHT = A+Q*(XR-A)
C                                  BUFFER TO PREVENT COALESCING OF KNOTS
      IF (AHIGH.GE.ARIGHT) GO TO 35
      AA = AHIGH
      GO TO 90
   35 CONTINUE
      CALL ICSFKU (X,F,NX,3,VARIGH,1,Y,C,IC,ERIGHT,WK,IER)
      IF (E.LE.ERIGHT) GO TO 50
C                                  CHECK TERMINATION CRITERION
      IF (E-ERIGHT.LE.BD) GO TO 85
      LPCNT = LPCNT+1
      IF (LPCNT.GE.INDLP) GO TO 85
      IF (ELEFT.GT.E) GO TO 25
C                                  SEARCH FOR A NEW KNOT TO THE LEFT
   40 ARIGHT = A
      ERIGHT = E
      A = ALEFT
      E = ELEFT
      ALEFT = A+Q*(XL-A)
C                                  BUFFER TO PREVENT COALESCING OF KNOTS
      IF (ALEFT.GE.ALOW) GO TO 45
      AA = ALOW
      GO TO 90
   45 CONTINUE
      CALL ICSFKU (X,F,NX,3,VALEFT,1,Y,C,IC,ELEFT,WK,IER)
      IF (E.LE.ELEFT) GO TO 50
C                                  CHECK TERMINATION CRITERION
      IF (E-ELEFT.LE.BD) GO TO 80
      LPCNT = LPCNT+1
      IF (LPCNT.GE.INDLP) GO TO 80
      GO TO 40
C                                  THE NEW GUESS ,ABEST, IS DETERMINED
C                                    TO MINIMIZE THE QUADRATIC P(T)
C                                    THAT INTERPOLATES (ALEFT,ELEFT),
C                                    (A,E), AND (ARIGHT,ERIGHT)
   50 DXLEFT = ALEFT-A
      DXRGHT = ARIGHT-A
      DYLEFT = (ELEFT-E)/DXLEFT
      DYRGHT = (ERIGHT-E)/DXRGHT
      DIFF = DYLEFT-DYRGHT
      IF (DIFF.EQ.ZERO) GO TO 95
      DEL = HALF*(DXRGHT*DYLEFT-DXLEFT*DYRGHT)/DIFF
      IF (DEL.EQ.ZERO) GO TO 95
      EPRED = E+DEL*(DYRGHT+(DXRGHT-DEL)/(ARIGHT-ALEFT)*DIFF)
      ABEST = A+DEL
      CALL ICSFKU (X,F,NX,3,VABEST,1,Y,C,IC,EBEST,WK,IER)
C
C                                  SWITCH THREE POINTS AND (ABEST,EBEST)
      IF (EBEST.LE.E) GO TO 60
      IF (DEL.GT.ZERO) GO TO 55
      ALEFT = ABEST
      ELEFT = EBEST
      GO TO 75
   55 ARIGHT = ABEST
      ERIGHT = EBEST
      GO TO 75
   60 CONTINUE
      IF (DEL.GT.ZERO) GO TO 65
      ARIGHT = A
      ERIGHT = E
      GO TO 70
   65 ALEFT = A
      ELEFT = E
   70 A = ABEST
      E = EBEST
C                                  CHECK TERMINATION CRITERION
   75 CONTINUE
      IF (ABS(EPRED-EBEST).LT.HALF*BD) GO TO 95
      LPCNT = LPCNT+1
      IF (LPCNT.GT.INDLP) GO TO 95
C                                  TERMINATE ITERATION
      GO TO 50
   80 A = ALEFT
      E = ELEFT
      GO TO 95
   85 A = ARIGHT
      E = ERIGHT
      GO TO 95
   90 CONTINUE
      CALL ICSFKU (X,F,NX,3,VAA,1,Y,C,IC,ETRY,WK,IER)
      IF (E.LE.ETRY) GO TO 95
      A = AA
      E = ETRY
   95 CONTINUE
      CHANGE = CHANGE+ABS(A-XK(NXK-1))/H
      ERROR = E
      KVARY = KVARY-1
C                                  REORDER KNOTS
      XK(NXK-1) = XK(KVARY)
      XK(KVARY) = A
      IF (KVARY.EQ.1) GO TO 100
      XR = A
      A = XK(NXK-1)
      XL = XK(KVARY-1)
C                                  REMOVE KNOTS
C
      CALL ICSFKU (X,F,NX,1,XK,KVARY-NXK,Y,C,IC,ERROR,WK,IER)
C
C                                  INSERT KNOTS IN NEW ORDER
C                                    KNOT KVARY LAST
C
      CALL ICSFKU (X,F,NX,1,XK(KVARY),NXK-KVARY,Y,C,IC,ERROR,WK,IER)
      GO TO 20
C                                  ORDER OF KNOTS IN XK IS
C                                    2,3,...,NXK-1,1,NXK
C                                    REORDER TO 1,2,...,NXK
  100 CONTINUE
      XR = XK(NXK-1)
      NXK1 = NXK+1
      DO 105 I=3,NXK
         J = NXK1-I
         XK(J+1) = XK(J)
  105 CONTINUE
      XK(1) = XR
      CALL ICSFKU (X,F,NX,1,XK,2-NXK,Y,C,IC,ERROR,WK,IER)
      CALL ICSFKU (X,F,NX,1,XK(2),NXK-2,Y,C,IC,ERROR,WK,IER)
C
C                                  CHECK FOR RELATIVE CONVERGENCE
C
      IF (PREVER-ERROR.LE.EPSERR*PREVER) GO TO 9005
C
C                                  CHECK FOR ABSOLUTE CONVERGENCE
C
      IF (ERROR.LE.EPSERA) GO TO 9005
      ITER = ITER-1
      IF (ITER.GT.0) GO TO 15
      IER = 33
 9000 CONTINUE
      CALL UERTST (IER,6HICSVKU)
 9005 RETURN
      END
