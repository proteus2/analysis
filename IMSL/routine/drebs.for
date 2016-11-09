C   IMSL ROUTINE NAME   - DREBS
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - NOVEMBER 1, 1979
C
C   PURPOSE             - DIFFERENTIAL EQUATION SOLVER -
C                           EXTRAPOLATION METHOD
C
C   USAGE               - CALL DREBS (FCN,Y,X,N,JM,IND,JSTART,H,HMIN,
C                           TOL,R,S,WK,IER)
C
C   ARGUMENTS    FCN    - NAME OF SUBROUTINE FOR EVALUATING FUNCTIONS.
C                           (INPUT)
C                           THE SUBROUTINE ITSELF MUST ALSO BE PROVIDED
C                             BY THE USER AND IT SHOULD BE OF THE
C                             FOLLOWING FORM
C                               SUBROUTINE FCN(N,X,Y,YPRIME)
C                               REAL Y(N),YPRIME(N)
C                                    .
C                                    .
C                                    .
C                           FCN SHOULD EVALUATE YPRIME(1),...,YPRIME(N)
C                             GIVEN N,X, AND Y(1),...,Y(N).  YPRIME(I)
C                             IS THE FIRST DERIVATIVE OF Y(I) WITH
C                             RESPECT TO X.
C                           FCN MUST APPEAR IN AN EXTERNAL STATEMENT IN
C                             THE CALLING PROGRAM AND N,X,Y(1),...,Y(N)
C                             MUST NOT BE ALTERED BY FCN.
C                Y      - DEPENDENT VARIABLES, VECTOR OF LENGTH N.
C                           (INPUT AND OUTPUT)
C                         ON INPUT, Y(1),...,Y(N) SUPPLY INITIAL
C                           VALUES.
C                         ON OUTPUT, Y(1),...,Y(N) ARE REPLACED WITH
C                           AN APPROXIMATE SOLUTION AT X (AS SET ON
C                           OUTPUT).
C                X      - INDEPENDENT VARIABLE. (INPUT AND OUTPUT)
C                           ON INPUT, X SUPPLIES THE INITIAL VALUE.
C                           ON OUTPUT, X IS REPLACED WITH THE UPDATED
C                             VALUE OF THE INDEPENDENT VARIABLE.
C                N      - THE NUMBER OF EQUATIONS. (INPUT)
C                JM     - THE MAXIMUM ORDER OF THE RATIONAL APPROX-
C                           IMATION. (INPUT) JM MUST BE LESS THAN 7.
C                           A SUGGESTED VALUE IS JM=6. SEE REMARKS.
C                IND    - CONVERGENCE TYPE INDICATOR. (INPUT)
C                           IND = 1 SPECIFIES THE STANDARD ERROR TEST
C                           IND = 2 SPECIFIES THE RELATIVE ERROR TEST
C                           IND = 3 SPECIFIES THE ABSOLUTE ERROR TEST
C                           SEE REMARKS FOR FURTHER DETAILS.
C                JSTART - INDICATOR. (INPUT)
C                           THE USER MUST SET JSTART TO 0 OR -1
C                           JSTART = 0 IMPLIES PERFORM A STEP.
C                             THE FIRST STEP MUST BE DONE WITH THIS
C                             VALUE OF JSTART SO THAT THE SUBROUTINE
C                             CAN INITIALIZE ITSELF.
C                           JSTART = -1 IMPLIES REPEAT THE LAST STEP
C                             WITH A NEW VALUE OF H OR JM.
C                             THE INITIAL VALUES OF Y, S,
C                             AND X ARE SET TO THE INITIAL VALUES OF Y,
C                             S, AND X FROM THE MOST RECENT CALL TO
C                             DREBS WITH JSTART = 0.
C                H      - STEP SIZE. (INPUT AND OUTPUT)
C                         ON INPUT, H IS AN INITIAL GUESS FOR THE STEP
C                           SIZE.
C                         ON OUTPUT, H IS REPLACED BY A SUGGESTED STEP
C                           SIZE FOR THE NEXT STEP. THE SUGGESTED VALUE
C                           MAY BE LARGER OR SMALLER THAN THE ORIGINAL
C                           STEP SIZE.
C                HMIN   - THE SMALLEST PERMISSIBLE STEP SIZE. (INPUT)
C                           DREBS WILL DECREASE THE STEP SIZE
C                           UNTIL CONVERGENCE CAN BE OBTAINED.
C                TOL    - TOLERANCE FOR ERROR CONTROL. (INPUT)
C                R      - VECTOR OF LENGTH N. (OUTPUT)
C                         ON OUTPUT, R CONTAINS THE ABSOLUTE
C                           ERRORS IN EACH COMPONENT FOR
C                           THE CURRENT STEP.
C                S      - VECTOR OF LENGTH N. (INPUT AND OUTPUT)
C                         IF IND = 1,
C                           BEFORE THE FIRST CALL TO THE ROUTINE,
C                           S(I) SHOULD BE SET TO Y(I), I=1,...,N.
C                           ON OUTPUT, S CONTAINS THE LARGEST VALUE
C                           OF EACH Y COMPUTED SINCE THE START OF THE
C                           INTEGRATION.
C                         IF IND = 2,
C                           BEFORE THE FIRST CALL TO THE ROUTINE,
C                           S(I) SHOULD BE SET TO Y(I), I=1,...,N.
C                           ON OUTPUT, S CONTAINS THE LARGEST VALUE
C                           OF EACH Y COMPUTED DURING THE CURRENT STEP.
C                         IF IND = 3,
C                           BEFORE THE FIRST CALL TO THE ROUTINE,
C                           S(I) SHOULD BE SET TO 1, I=1,...,N.
C                           ON OUTPUT, S IS UNCHANGED.
C                WK     - WORK VECTOR OF LENGTH 29*N.
C                           WK MUST REMAIN UNCHANGED BETWEEN SUCCESSIVE
C                           CALLS DURING INTEGRATION.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER = 129 INDICATES CONVERGENCE COULD NOT BE
C                             OBTAINED WITH CURRENT VALUES OF H AND
C                             HMIN. Y,X, AND H HAVE BEEN UPDATED.
C                         WARNING ERROR (WITH FIX)
C                           IER = 66 INDICATES JM IS LESS THAN 1 OR
C                             GREATER THAN 6. JM IS RESET TO 6.
C
C   REQD. IMSL ROUTINES - UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS  1.  THE SOLUTION Y, THE INDEPENDENT VARIABLE X, AND THE
C                SUGGESTED STEP SIZE H ARE ALWAYS UPDATED EVEN IF
C                CONVERGENCE IS NOT OBTAINED.
C            2.  IN GENERAL, HMIN SHOULD BE MUCH SMALLER THAN H TO
C                ALLOW THE PROGRAM TO ADJUST FOR RAPIDLY CHANGING
C                SOLUTIONS (WITH RESPECT TO H).
C            3.  AT EACH STEP OF THE INTEGRATION, THE EXTRAPOLATION
C                PROCESS IS CONSIDERED TO HAVE CONVERGED WHEN EACH
C                Y(I), I=1,...,N, HAS SATISFIED A CONVERGENCE
C                CRITERION SPECIFIED BY THE USER. THE USER MAY CHOOSE
C                ONE OF THREE CONVERGENCE CRITERIA. IN TESTING FOR
C                CONVERGENCE, TWO SUCCESSIVE EXTRAPOLATED VALUES
C                (FOR EACH COMPONENT) AT THE POINT IN QUESTION ARE
C                COMPARED. LET THE DIFFERENCE BETWEEN THE TWO FOR
C                THE J-TH COMPONENT BE CALLED D(J). THE THREE
C                CONVERGENCE CRITERIA CAN BE STATED IN THE FOLLOWING
C                MANNER;
C
C                A. STANDARD ERROR
C                   LET YMAX(J) BE THE LARGEST ABSOLUTE VALUE ATTAINED
C                   SO FAR IN THE INTEGRATION BY THE DEPENDENT VARIABLE
C                   Y(J). THE CONVERGENCE REQUIREMENT IS
C                     ABS(D(J)/YMAX(J)) .LE. TOL, FOR J=1,...,N
C                   IF YMAX(J) IS LESS THAN TOL, IT IS REPLACED IN THE
C                   TEST BY TOL.
C
C                B. RELATIVE ERROR
C                   LET Y(J) BE THE CURRENT APPROXIMATION TO THE
C                   RESPECTIVE DEPENDENT VARIABLE. THE CONVERGENCE
C                   REQUIREMENT IS
C                     ABS(D(J)/Y(J)) .LE. TOL, FOR J=1,...,N
C                   IF ABS(Y(J)) IS LESS THAN TOL, IT IS REPLACED IN THE
C                   TEST BY TOL.
C
C                C. ABSOLUTE ERROR
C                   THE CONVERGENCE REQUIREMENT IS
C                     ABS(D(J)) .LE. TOL, J=1,...,N
C            4.  JM, THE ORDER OF THE RATIONAL APPROXIMATION, DOES NOT
C                HAVE TO EQUAL 1, THE ORDER OF THE DIFFERENTIAL
C                EQUATIONS. AT EACH INTEGRATION STEP, AS MANY AS JM
C                APPLICATIONS OF THE MIDPOINT RULE ARE COMPUTED FOR
C                SUCCESSIVELY SMALLER VALUES OF H AND EXTRAPOLATED TO
C                H=0 IN ATTEMPTING TO ACHIEVE CONVERGENCE.
C                TYPICAL USAGE OF DREBS WOULD BE WITH JM SET TO 6.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE DREBS  (FCN,Y,X,N,JM,IND,JSTART,H,HMIN,TOL,R,S,WK,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N,JM,IND,JSTART,IER
      REAL               Y(N),X,H,HMIN,TOL,R(N),S(N),WK(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,IJK4,J,JMAX,L,IJK2,IJK1,K,M,IJK7,IJK6,IJK5
      INTEGER            IP5,N2,N3,N4,N5,N6,N7,N8,N2P1,N3P1,N8P1,JR,JS
      INTEGER            JJ,IJ6,IJ7,I6,I7,IJK3,KH,IJK8,IK5
      REAL               D(7),ZOTUP,A,G,B,XU,U,C,TA,B1,V,UST,XOLD,HALF
      REAL               EP,SRH,POS,FC
      LOGICAL            KONVF,KONV,BO,BH
      DATA               EP/Z3C100000/
      DATA               SRH/.7071068/
      DATA               POS/.1666667/
      DATA               HALF/.5/,XOLD/0.0/
C                                  FIRST EXECUTABLE STATEMENT
      IER=0
      N2=N+N
      N3=N2+N
      N4=N3+N
      N5=N4+N
      N6=N5+N5+N2
      N7=N6+N5+N3
      N8=N7+N5+N3
      N2P1=N2+1
      N3P1=N3+1
      N8P1=N8+1
      ZOTUP = 64.0/EP
      IF(JM.GT.0.AND.JM.LE.6) GO TO 5
      IER = 66
      JM=6
C                                  FOR AN EXTRAPOLATION OF ORDER JM,
C                                  JM+1 APPROXIMATIONS ARE REQUIRED.
C                                  THREE MORE ARE ALLOWED IN ATTEMPTING
C                                  TO ACHIEVE CONVERGENCE.
    5 JMAX = JM+4
      IF(JSTART.NE.0) GO TO 135
C                                  INITIALIZATION
C                                  SAVE THE INITIAL VALUES FOR THE DEP-
C                                  ENDENT VARIABLES AND THE ERROR TEST
C                                  VECTOR FOR THE STEP
      XOLD = X
      IJK4=N4
      DO 10 I = 1,N
         WK(I) = Y(I)
         IJK4=IJK4+1
         WK(IJK4) = S(I)
   10 CONTINUE
C                                  USE THE FUNCTION ROUTINE TO OBTAIN
C                                  THE INITIAL SLOPES
C                                  DZ = DY/DX
      CALL FCN(N,X,Y,WK(N3P1))
C                                  THE LOGICAL VARIABLE,BH, DETERMINES
C                                  WHETHER THE STEPSIZE HAS BEEN
C                                  HALVED. INITIALLY FALSE.
C                                  LATER BH IS FALSE IF THE STEPSIZE IS
C                                  CUT BY A FACTOR NOT 2
   15 BH = .FALSE.
C                                  PRESET THE CONVERGENCE SUCCESS FLAG
C                                     TRUE
      KONVF = .TRUE.
C                                  ADVANCE THE INDEPENDENT VARIABLE BY
C                                  THE STEPSIZE, H
   20 A = H+X
C                                  SET THE SWITCH BO FOR THE FIRST SET
C                                  OF COEFFICIENTS, D
      BO = .FALSE.
C                                  INITIALIZE THE H SEQUENCE...
C                                     H/M,H/JR,H/JS
      M = 1
      JR = 2
      JS = 3
C                                  JJ IS THE INDEX FOR THE ARRAY OF
C                                  VALUES SAVED IN CASE THE INTERVAL
C                                  MUST BE HALVED
      JJ = 0
C                                  INTEGRATION STEP
C                                  MIDPOINT + EXTRAPOLATION
      IJ6=N6-N
      IJ7=N7-N
      I6=IJ6
      I7=IJ7
      DO 125 J = 1,JMAX
C                                  SET THE VALUES OF THE EXTRAPOLATION
C                                  COEFFICIENTS TO THEIR CORRECT VALUES
C                                  FOR THIS EXTRAPOLATION STEP
         IJ6=IJ6+N
         IJ7=IJ7+N
         IF (.NOT. BO) GO TO 25
         D(2)=1.777778
         D(4)=7.111111
         D(6)=28.44444
         GO TO 30
   25    D(2)=2.25
         D(4) = 9.0
         D(6) = 36.0
C                                  IF THE ORDER OF THE EXTRAPOLATION
C                                  STEP BEING COMPUTED IS LESS THAN JM/2
C                                  SET KONV FALSE
   30    KONV = .TRUE.
         IF (J .LE. (JM/2)) KONV = .FALSE.
         IF (J .LE. (JM+1)) GO TO 35
C                                  RESTRICT THE ORDER OF THE EXTRAPOL-
C                                  ATION TO JM
C                                  ADJUST THE EXTRAPOLATION COEFFICIENT
         L = JM+1
         D(L) = 4.0*D(L-2)
C                                  DISCOURAGE THE STEP-INCREASING FACTOR
C                                  FC, BY A FACTOR OF SQRT(.5) SINCE
C                                  CONVERGENCE WAS NOT OBTAINED
C                                  IN JM EXTRAPOLATIONS
         FC = SRH*FC
         GO TO 40
C                                  THE NUMBER, J, OF EXTRAPOLATIONS HAS
C                                  NOT EXCEEDED JM
C                                  FIND D(J) = ((H DIVIDED BY H/M)**2
C                                  ADJUST THE FACTOR,FC, USED TO ADJUST
C                                  THE STEPSIZE FOR THE NEXT STEP TO BE
C                                  TAKEN
   35    L = J
         D(L) = M*M
         FC = 1.0+(JM+1-J)*POS
C                                  MODIFIED MIDPOINT RULE USED TO FIND
C                                  FIRST VALUE FOR THIS EXTRAPOLATION
C                                  STEP
   40    M = M+M
         G = H/M
         B = G+G
C                                  IF THE STEPSIZE HAS NOT BEEN HALVED
C                                  OR IF THE ORDER OF THE EXTRAPOLATION
C                                  STEP EXCEEDS THAT FOR WHICH PREVIOUS-
C                                  LY COMPUTED VALUES WERE SAVED, THEY
C                                  MUST BE COMPUTED
         IF ((.NOT. BH) .OR. (J .GE. (JMAX-1))) GO TO 50
C                                  OTHERWISE THE VALUES HAVE BEEN SAVED
C                                  AND CAN BE RESTORED
         IJK1=N
         IJK2=N2
         IJK6=IJ6
         IJK7=IJ7
         DO 45 I = 1,N
            IJK2=IJK2+1
            IJK7=IJK7+1
            WK(IJK2) = WK(IJK7)
            IJK1=IJK1+1
            IJK6=IJK6+1
            WK(IJK1) = WK(IJK6)
   45    CONTINUE
         GO TO 75
C                                  COMPUTE STARTING VALUES FOR THE MODI-
C                                  FIED MIDPOINT RULE
   50    IJK1=N
         IJK2=N2
         IJK3=N3
         DO 55 I = 1,N
            IJK1=IJK1+1
            WK(IJK1) = WK(I)
            IJK2=IJK2+1
            IJK3=IJK3+1
            WK(IJK2) = WK(I)+G*WK(IJK3)
   55    CONTINUE
         KH = M/2
         XU=X
C                                  THE MEMBER OF THE H SEQUENCE
C                                  BEING USED BY THE MIDPOINT INTEGRA-
C                                  TION RULE IS H/M. COMPUTE THE END
C                                  OF THE STEP  FOR EACH DEPENDENT VARI-
C                                  ABLE
         DO 70 K = 2,M
            XU = XU + G
            CALL FCN(N,XU,WK(N2P1),WK(N8P1))
            IJK1=N
            IJK2=N2
            IJK8=N8
            DO  60 I = 1,N
               IJK1=IJK1+1
               IJK8=IJK8+1
               U = WK(IJK1) + B*WK(IJK8)
               IJK2=IJK2+1
               WK(IJK1) = WK(IJK2)
               WK(IJK2) = U
   60       CONTINUE
C                                  IN CASE THE INTERVAL MUST BE HALVED
C                                  NEXT TIME, SAVE THE VALUES AT HALFWAY
C                                  ALONG (KH=M/2) THE STEP UNLESS K=3
            IF ((K .NE. KH) .OR. (K .EQ. 3)) GO TO 70
            JJ = 1+JJ
            I6=I6+N
            I7=I7+N
            IJK1=N
            IJK2=N2
            IJK6=I6
            IJK7=I7
            DO 65 I = 1,N
               IJK2=IJK2+1
               IJK7=IJK7+1
               WK(IJK7) = WK(IJK2)
               IJK6=IJK6+1
               IJK1=IJK1+1
               WK(IJK6) = WK(IJK1)
   65       CONTINUE
   70    CONTINUE
   75    CALL FCN(N,A,WK(N2P1),WK(N8P1))
         IJK1=N
         IJK2=N2
         IJK5=N5
         IK5=N5
         IJK8=N8
         DO 115 I = 1,N
C                                  V IS USED TO SAVE THE VALUE OBTAINED
C                                  BY THE MIDPOINT RULE USING THE PREVI-
C                                  OUS MEMBER OF THE H SEQUENCE
C                                  (THE FIRST TIME THROUGH THIS VALUE IS
C                                  MEANINGLESS, BUT IT IS NOT USED SINCE
C                                  L IS LESS THAN 2)
            IJK1=IJK1+1
            IJK2=IJK2+1
            IJK5=IJK5+1
            IK5=IK5+1
            IF (L .GE. 2) V = WK(IJK5)
C                                  COMPUTE THE FINAL VALUE OBTAINED FOR
C                                  THIS MEMBER OF THE H SEQUENCE BY THE
C                                  MODIFIED MIDPOINT RULE
            IJK8=IJK8+1
            WK(IJK5) = (WK(IJK2) + WK(IJK1) + G * WK(IJK8))*HALF
            C = WK(IJK5)
            TA = C
C                                  AT LEAST TWO VALUES ARE NEEDED TO
C                                  START EXTRAPOLATION
            IF (L .LT. 2) GO TO 90
C                                  IF THE VALUE JUST COMPUTED BY THE
C                                  MIDPOINT RULE SHOWS A LARGE JUMP FROM
C                                  THE PREVIOUS, HALVE THE INTERVAL
            IF ((ABS(V)*ZOTUP .LT. ABS(C)) .AND. (H .NE. HMIN) .AND.
     1           (J.GT.JM/2+1)) GO TO 130
C                                  PERFORM THE L STEPS FOR THE CURRENT
C                                  LTH ORDER EXTRAPOLATION STEP. IF THE
C                                  DENOMINATOR OF THE RATIONAL FUNCTION
C                                  GOES TO ZERO AT ANY STEP, SET DT AT
C                                  THAT STEP TO ITS VALUE JUST BEFORE
            IP5=IK5
            DO 85 K = 2,L
               IP5=IP5+N
               B1 = D(K) * V
               B = B1-C
               U = V
               IF (B .EQ. 0.) GO TO 80
               B = (C-V)/B
               U = C*B
               C = B1*B
   80          IF(K.LT.L) V=WK(IP5)
               WK(IP5) = U
               TA = U + TA
   85       CONTINUE
C                                  USE THE ERROR ROUTINE FOR EACH
C                                  DEPENDENT VARIABLE TO CHECK WHETHER
C                                  CONVERGENCE HAS BEEN ACHIEVED
   90       GO TO (95,100,105),IND
   95       UST=ABS(TA)
            IF(UST.GT.S(I)) S(I)=UST
            GO TO 110
  100       S(I)=ABS(Y(I))
            GO TO 110
  105       S(I)=1.0
  110       R(I)=ABS(Y(I)-TA)
            Y(I)=TA
            IF(S(I).LT.TOL) S(I)=TOL
            IF(R(I).GT.TOL*S(I)) KONV=.FALSE.
  115    CONTINUE
         IF(KONV) GO TO 155
C                                  RESET THE EXTRAPOLATION COEFFICIENTS
         D(3) = 4.
         D(5) = 16.
C                                  FLIP THE BO SWITCH FOR THE NEXT SET
C                                  OF COEFFICIENTS
         BO = (.NOT.BO)
C                                  RESET S
         IJK4=N4
         DO 120 I=1,N
            IJK4=IJK4+1
            S(I)=WK(IJK4)
  120    CONTINUE
C                                  TAKE THE NEXT MEMBER
C                                  OF THE H SEQUENCE
         M = JR
         JR = JS
         JS = M+M
C                                  AND GO BACK FOR THE
C                                  NEXT EXTRAPOLATION
  125 CONTINUE
C                                  IF, AFTER ALL THE EXTRAPOLATIONS
C                                  ALLOWED, CONVERGENCE HAS NOT BEEN
C                                  ACHIEVED, ATTEMPT TO HALVE H SO THAT
C                                  THE SAVED VALUES CAN BE USED (SET BH
C                                  TRUE FOR THIS PURPOSE)
C                                  IF HALVING H MAKES IT LESS THAN HMIN,
C                                  SET H = HMIN
C                                  IN THIS CASE THE SAVED VALUES CANNOT
C                                  BE USED
C                                  IF H HAD ALREADY BEEN AT HMIN,
C                                  CONVERGENCE CANNOT BE ACHIEVED FOR
C                                  THIS HMIN AND THIS TOL CRITERION.
C                                  SET KONV FALSE
      BH = (.NOT. BH)
  130 IF (ABS(H) .LE. ABS(HMIN)) GO TO 150
      H = H * HALF
      IF (ABS(H) .GE. ABS(HMIN)) GO TO 20
      H = SIGN(HMIN,H)
      GO TO 15
  135 DO 140 I = 1,N
         Y(I) = WK(I)
  140 CONTINUE
      IJK4=N4
      DO 145 I=1,N
         IJK4=IJK4+1
         S(I)=WK(IJK4)
  145 CONTINUE
      X = XOLD
      GO TO 15
  150 KONVF = .FALSE.
C                                  WHETHER OR NOT CONVERGENCE HAS BEEN
C                                  ACHIEVED SET A NEW SUGGESTED STEPSIZE
C                                  FOR THE NEXT STEP
C                                  ASSIGN THE END OF STEP VALUE TO THE
C                                  INDEPENDENT VARIABLE
  155 H = FC*H
      X=A
      IF(KONVF) GO TO 160
      IER=129
  160 IF(IER.EQ.0) GO TO 9005
 9000 CONTINUE
      CALL UERTST (IER,6HDREBS )
 9005 RETURN
      END
