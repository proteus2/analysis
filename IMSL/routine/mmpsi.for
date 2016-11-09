C   IMSL ROUTINE NAME   - MMPSI
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1982
C
C   PURPOSE             - LOGARITHMIC DERIVATIVE OF THE GAMMA FUNCTION
C
C   USAGE               - FUNCTION MMPSI (ARG,IER)
C
C   ARGUMENTS    MMPSI  - OUTPUT VALUE OF THE FUNCTION AT ARG. MMPSI
C                           MUST BE TYPED APPROPRIATELY IN THE CALLING
C                           PROGRAM. (SEE THE PRECISION/HARDWARE
C                           SECTION.)
C                ARG    - INPUT ARGUMENT. THE ABSOLUTE VALUE OF ARG
C                           MUST BE GREATER THAN XMIN, -ARG MUST BE
C                           LESS THAN XMAX, AND ARG MUST NOT BE A
C                           NEGATIVE INTEGER. XMIN IS OF THE ORDER OF
C                           10**(-38) AND XMAX IS OF THE ORDER OF
C                           10**10. THE EXACT VALUES OF XMIN AND XMAX
C                           MAY ALLOW LARGER RANGES FOR ARG ON SOME
C                           COMPUTERS. SEE THE PROGRAMMING NOTES IN
C                           THE MANUAL FOR THE EXACT VALUES.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER = 129 INDICATES ARG IS LESS THAN OR
C                             EQUAL TO -XMAX. MMPSI IS SET TO ZERO.
C                           IER = 130 INDICATES ARG IS A NEGATIVE
C                             INTEGER (OR ZERO) OR -XMIN IS LESS THAN
C                             ARG WHICH IS LESS THAN ZERO. MMPSI IS SET
C                             TO MACHINE INFINITY.
C                           IER = 131 INDICATES ZERO IS LESS THAN ARG
C                             WHICH IS LESS THAN XMIN. MMPSI IS SET TO
C                             NEGATIVE MACHINE INFINITY.
C
C   REQD. IMSL ROUTINES - UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      DOUBLE PRECISION FUNCTION MMPSI (ARG,IER)
C
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            IER
      DOUBLE PRECISION   ARG
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,M,N,NQ,IEND,IPEND,IQEND,JEND,JENDP1
      DOUBLE PRECISION   P1(8),P2(5),Q1(7),Q2(5)
      DOUBLE PRECISION   AUG,DEN,FUDGE,PIOV4,SCALE,SGN,XSMALL,UPPER
      DOUBLE PRECISION   W,X01,X02,XINF,XMINS,Z,XLARGE,XLGER,XMAX1,
     *                   X,XMX
      DOUBLE PRECISION   XDINT
C                                  MACHINE DEPENDENT CONSTANTS
C                                    PIOV4 = PI / 4
C                                    SCALE = SCALING FACTOR
C                                    XINF = LARGEST POSITIVE MACHINE
C                                    NUMBER
C                                    XMAX1 = SMALLEST POSITIVE
C                                    FLOATING-POINT CONSTANT WITH
C                                    ENTIRELY INTEGER REPRESENTATION.
C                                    ALSO USED AS NEGATIVE OF LOWER
C                                    BOUND ON ACCEPTABLE NEGATIVE
C                                    ARGUMENTS
C                                    XLGER = THE POSITIVE
C                                    ARGUMENT BEYOND WHICH PSI MAY BE
C                                    REPRESENTED AS DLOG(ARG)
C                                    XLARGE = ONE FOURTH THE VALUE
C                                    OF XLGER.  USED TO ALLOW THE
C                                    VALUE OF XLGER TO BE ON THE ORDER
C                                    OF 10**17.
C                                    XMINS= SMALLEST ACCEPTABLE ABSOLUTE
C                                    ARGUMENT SCALED BY SCALE
C                                    XSMALL = ABSOLUTE ARGUMENT BELOW
C                                    WHICH PI*COTAN(PI*ARG) MAY BE
C                                    REPRESENTED BY 1/ARG
C                                    X01 + X02 = ZERO OF PSI TO
C                                    EXTENDED PRECISION
      DATA               PIOV4/.7853981633974483D0/
      DATA               SCALE/1.0D0/
      DATA               XINF/.723700557733226D+76/
      DATA               XMAX1/.4503599627370496D16/
      DATA               XMINS/Z0210000000000001/
      DATA               XSMALL/Z3940000000000000/
      DATA               X01/Z411762D860000000/
      DATA               X02/.6219074969459455D-08/
      DATA               FUDGE/Z4100000000000001/
      DATA               XLARGE/.9007199254740992D16/
      DATA               XDINT/Z4E10000000000000/
C                                  COEFFICIENTS FOR RATIONAL
C                                  APPROXIMATION OF PSI(ARG) / (ARG -
C                                    X0), 0.5 .LE. ARG .LE. 3.0
      DATA               P1/-.4937716493081015D00,-.2855677555817539D02,
     1                   -.3816935645399708D03,-.1211830788802065D04,
     2                   .3833868258616269D04,.2450758805205975D05,
     3                   .3540131456802392D05,.1352499966737940D05/
      DATA               Q1/.3371456475818930D02,.6223738892835428D03,
     1                   .4540983303742759D04,.1469514355996634D05,
     2                   .2062758041767692D05,.9884287131523368D04,
     3                   .3469455587688172D-06/
C                                  COEFFICIENTS FOR RATIONAL
C                                    APPROXIMATION OF PSI(ARG) - LN(ARG)
C                                    + 1 / (2*ARG), ARG .GT. 3.0
      DATA               P2/-.2431393158434656D01,
     1                   -.1077240563464793D02,-.1042268336388353D02,
     2                   -.3050247680803867D01,-.2461513967345629D0/
      DATA               Q2/.3868046608354867D02,.1405216313263703D03,
     1                   .1286213778152642D03,.3689835384569604D02,
     2                   .2953816760814839D01/
C                                  FIRST EXECUTABLE STATEMENT
      X = ARG
      IEND = 6
      IPEND = 8
      IQEND = 7
      JEND = 4
      JENDP1 = 5
      IER = 0
      AUG = 0.0D0
      XLGER = XLARGE*4.0D0
      IF (X.GE.0.5D0) GO TO 25
C                                  ARG .LT. 0.5, USE REFLECTION
C                                    FORMULA PSI(1-ARG) = PSI(ARG) + PI
C                                    * COTAN(PI*ARG)
      IF (DABS(X).GT.XSMALL) GO TO 5
      IF (DABS(X)*SCALE.LT.XMINS) GO TO 55
C                                  XMIN .LT. DABS(ARG) .LE. XSMALL. USE
C                                    1/ARG AS A SUBSTITUTE FOR
C                                    PI*COTAN(PI*ARG)
      AUG = -1.0D0/X
      GO TO 20
C                                  REDUCTION OF ARGUMENT FOR COTAN
    5 W = -X
      SGN = PIOV4
      IF (W.GT.0.0D0) GO TO 10
      W = -W
      SGN = -SGN
C                                  MAKE AN ERROR EXIT IF ARG .LE.
C                                    -XMAX1
   10 XMX = XMAX1
      IF (W.GE.XMX) GO TO 50
      Z = W+XDINT
      W = W-(Z-XDINT)
      NQ = IDINT(W*4.0D0)
      W = 4.0D0 * (W - (DBLE(FLOAT(NQ))*.25D0))
C                                  W IS NOW RELATED TO THE FRACTIONAL
C                                    PART OF 4.0 * ARG. ADJUST ARGUMENT
C                                    TO CORRESPOND TO VALUES IN FIRST
C                                    QUADRANT AND DETERMINE SIGN
      N = NQ/2
      IF ((N+N).NE.NQ) W = 1.0D0-W
      Z = PIOV4*W
      M = N/2
      IF ((M+M).NE.N) SGN = -SGN
C                                  DETERMINE FINAL VALUE FOR
C                                    -PI*COTAN(PI*ARG)
      N = (NQ+1)/2
      M = N/2
      M = M+M
      IF (M.NE.N) GO TO 15
C                                  CHECK FOR SINGULARITY
      IF (Z.EQ.0.0D0) GO TO 55
C                                  USE COS/SIN AS A SUBSTITUTE FOR
C                                    COTAN, AND SIN/COS AS A SUBSTITUTE
C                                    FOR TAN
      AUG = SGN*((DCOS(Z)/DSIN(Z))*4.0D0)
      GO TO 20
   15 AUG = SGN*((DSIN(Z)/DCOS(Z))*4.0D0)
   20 X = 1.0D0-X
   25 IF (X.GT.3.0D0) GO TO 35
C                                  0.5 .LE. ARG .LE. 3.0
      DEN = X*0.5D0
      UPPER = P1(1)*X
      DO 30 I=1,IEND
         DEN = (DEN+Q1(I))*X
         UPPER = (UPPER+P1(I+1))*X
   30 CONTINUE
      DEN = (UPPER+P1(IPEND))/(DEN+Q1(IQEND))
      X = (X - X01)*.5D0 - X02
      MMPSI = DEN * X + X + AUG
      GO TO 9005
C                                  IF ARG .GE. XMAX1, PSI = LN(ARG)
   35 IF (AUG.EQ.0.0D0) AUG = FUDGE
      IF (X.GE.XLGER) GO TO 45
C                                  3.0 .LT. ARG .LT. XMAX1
      W = 1.0D0/(X*X)
      DEN = W
      UPPER = P2(1)*W
      DO 40 I=1,JEND
         DEN = (DEN+Q2(I))*W
         UPPER = (UPPER+P2(I+1))*W
   40 CONTINUE
      AUG = UPPER/(DEN+Q2(JENDP1))-0.5D0/X+AUG
   45 MMPSI = AUG+DLOG(X)
      GO TO 9005
C                                  ERROR RETURN FOR ARG .LE. -XMAX1
   50 MMPSI = 0.0D0
      IER = 129
      GO TO 9000
C                                  ERROR RETURN FOR -ARG AN INTEGER OR
C                                    DABS(ARG) .LT. XMIN
   55 MMPSI = XINF
      IER = 130
      IF (X .LE. 0.0D0) GO TO 9000
      MMPSI = -XINF
      IER = 131
C                                  UPDATE ERROR COUNTS, ETC.
 9000 CONTINUE
      CALL UERTST (IER,6HMMPSI )
 9005 RETURN
      END
