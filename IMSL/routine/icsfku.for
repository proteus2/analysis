C   IMSL ROUTINE NAME   - ICSFKU
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - LEAST SQUARES APPROXIMATION BY CUBIC SPLINES -
C                           FIXED KNOTS.
C
C   USAGE               - CALL ICSFKU (X,F,NX,MODE,XK,NXK,Y,C,IC,ERROR,
C                           WK,IER)
C
C   ARGUMENTS    X      - VECTOR OF LENGTH NX CONTAINING THE ABSCISSAE
C                           OF THE NX DATA POINTS (X(I),F(I)) I=1,...,
C                           NX. X MUST BE ORDERED SO THAT X(1) .LE. X(2)
C                           ... .LE. X(NX). (INPUT)
C                F      - VECTOR OF LENGTH NX CONTAINING THE ORDINATE
C                           (OR FUNCTION VALUES) OF THE NX DATA POINTS.
C                           (INPUT)
C                NX     - NUMBER OF ELEMENTS IN X AND F. (INPUT)
C                MODE   - AN INTEGER BETWEEN 0 AND 4 THAT SPECIFIES THE
C                           FUNCTION TO BE PERFORMED. (INPUT)
C                             MODE = 0, IN THIS CASE ICSFKU COMPUTES THE
C                               LEAST SQUARES APPROXIMATION TO THE DATA
C                               (X(I),F(I)), I=1,...,NX BY A CUBIC
C                               SPLINE WITH NXK KNOTS LOCATED AT XK(1),
C                               ...,XK(NXK).  THE KNOTS NEED NOT BE OR-
C                               DERED BUT MUST BE DISTINCT AND MUST ALL
C                               LIE BETWEEN THE BOUNDARY KNOTS, XK(1)
C                               AND XK(NXK).  ALSO, THE ABSCISSAE OF
C                               THE NX DATA POINTS MUST LIE BETWEEN THE
C                               BOUNDARY KNOTS (I.E., XK(1) .LE. X(I)
C                               .LE. XK(NXK) FOR I=1,...,NX.).
C                             MODE = 1, IN THIS CASE ICSFKU CAN COMPUTE
C                               A NEW LEAST SQUARES APPROXIMATION BY
C                               ADDING KNOTS OR BY DELETING KNOTS FROM
C                               THE CURRENT KNOT SET (SEE PROGRAMMING
C                               NOTES).  IF NXK IS NEGATIVE ICSFKU
C                               WILL REMOVE -NXK KNOTS.
C                               KNOTS ARE REMOVED IN THE REVERSE ORDER
C                               OF INSERTION (LAST IN FIRST OUT).  THE
C                               BOUNDARY KNOTS ARE NEVER REMOVED.  IF
C                               NXK IS POSITIVE ICSFKU WILL COMPUTE A
C                               NEW SPLINE APPROXIMATION WITH THE NEW
C                               KNOTS XK(1),...,XK(NXK) ADDED TO THE
C                               CURRENT KNOT SET.  THESE NEW KNOTS NEED
C                               NOT BE ORDERED BUT MUST LIE BETWEEN THE
C                               BOUNDARY KNOTS.
C                             MODE = 2, SAME OPTIONS AS MODE = 1 BUT OLD
C                               BASIS FUNCTIONS ARE USED IF THEY HAVE
C                               BEEN COMPUTED IN PREVIOUS CALLS.
C                             MODE = 3, IN THIS CASE ICSFKU WILL COMPUTE
C                               A NEW LEAST SQUARES APPROXIMATION WITH A
C                               NEW VALUE FOR THE LAST KNOT ENTERED INTO
C                               THE KNOT SET.  XK(1) IS ASSUMED TO CON-
C                               TAIN THE NEW KNOT LOCATION.
C                             MODE = 4, IN THIS CASE ICSFKU WILL RETURN
C                               THE CURRENT KNOT SET (PROPERLY ORDERED)
C                               AND THE NUMBER OF KNOTS IN THAT SET.
C                               ON OUTPUT, NXK CONTAINS THE NUMBER OF
C                               KNOTS IN THE CURRENT KNOT SET AND XK(I)
C                               CONTAINS THE VALUE OF THE I-TH KNOT FOR
C                               I=1,2,...,NXK.  NOTE THAT IN THIS CASE
C                               NXK AND XK ARE USED AS OUTPUT
C                               PARAMETERS, MODE IS USED AS AN INPUT
C                               PARAMETER, AND NO OTHER PARAMETERS ARE
C                               USED.
C                XK     - VECTOR OF LENGTH NXK CONTAINING THE KNOT LOCA-
C                           TIONS. (INPUT)
C                NXK    - NUMBER OF ELEMENTS IN XK. (INPUT)
C                Y,C    - SPLINE COEFFICIENTS (OUTPUT).  Y IS A VECTOR
C                           OF LENGTH NXKMAX-1 AND C IS AN NXKMAX-1 BY
C                           3 MATRIX.
C                           NXKMAX IS THE MAXIMUM NUMBER OF KNOTS IN
C                           THE SPLINE APPROXIMATION.
C                           IF ONLY MODE=0 IS USED, NXKMAX=NXK.
C                             THE VALUE OF THE SPLINE APPROXIMATION AT
C                             T IS
C                               S(T) = ((C(I,3)*D+C(I,2))*D+C(I,1))*D+
C                                      Y(I)
C                             WHERE XK(I) .LE. T .LT. XK(I+1) AND
C                             D = T - XK(I).(THIS FORMULA ASSUMES
C                             THAT KNOTS ARE ORDERED).
C                             THE SPLINE COEFFICIENTS ARE NOT RETURNED
C                             WHEN MODE = 3.
C                IC     - ROW DIMENSION OF MATRIX C EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM. (INPUT)
C                ERROR  - LEAST SQUARES ERROR OF THE CUBIC SPLINE
C                           APPROXIMATION. (OUTPUT)
C                WK     - WORK AREA OF DIMENSION NX*(NXKMAX+6)
C                           NXKMAX IS THE MAXIMUM NUMBER OF KNOTS IN
C                           THE SPLINE APPROXIMATION.
C                           IF ONLY MODE=0 IS USED, NXKMAX=NXK.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER = 129 ATTEMPT TO INSERT A KNOT WHICH IS
C                             NOT WITHIN THE BOUNDARY KNOTS OR IS
C                             COINCIDENT WITH A PREVIOUS KNOT.
C                           IER = 130 ATTEMPT TO USE MORE THAN 28 KNOTS
C                             TO APPROXIMATE DATA.
C                           IER = 131 INPUT ABSCISSAE ARE NOT ORDERED SO
C                             THAT X(1) .LE. X(2) ... .LE. X(NX)
C                           IER = 132 INPUT ABSCISSAE DO NOT LIE
C                             BETWEEN THE BOUNDARY KNOTS
C                           IER = 133 KNOTS ARE NOT DISTINCT OR DO NOT
C                             LIE BETWEEN THE BOUNDARY KNOTS.
C
C   REQD. IMSL ROUTINES - ICSFKV,UERTST,UGETIO
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
C            2.  IF THE KNOTS ARE ORDERED, THEN INTERVAL I HAS
C                ENDPOINTS XK(I) AND XK(I+1). THE CUBIC SPLINE FUNCTION
C                IS GIVEN BY
C                  S(T) = ((C(I,3)*D+C(I,2))*D+C(I,1))*D+Y(I)
C                WHERE T IS IN INTERVAL I AND D = T-XK(I). THE SPLINE
C                COEFFICIENTS ARE ALWAYS COMPUTED AS IF THE KNOTS ARE
C                ORDERED, AND ORDERING THE KNOTS IS ADVISED FOR EASE
C                OF OUTPUT USAGE.
C            3.  ICSFKU IS INTENDED TO BE USED FOR FUNCTIONS THAT CAN
C                BE APPROXIMATED ADEQUATELY WITH RELATIVELY FEW KNOTS.
C                IT HAS A BUILT-IN LIMIT OF 28 KNOTS.
C            4.  THE ARRAYS X, F AND THE WORK AREA WK MUST BE PRESERVED
C                BETWEEN CALLS TO ICSFKU. OTHER INFORMATION STORED IN
C                COMMON/ICSFK1/ MUST ALSO BE PRESERVED.
C            5.  THE MOST COMMON USE OF ICSFKU WILL BE WITH MODE=0.
C                MODE 1, 2, AND 3 ARE PROVIDED PRIMARILY FOR THE USE OF
C                ICSFKU BY THE IMSL VARIABLE KNOT SUBROUTINE ICSVKU.
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
      SUBROUTINE ICSFKU (X,F,NX,MODE,XK,NXK,Y,C,IC,ERROR,WK,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            NX,MODE,NXK,IC,IER
      REAL               X(NX),F(NX),XK(1),Y(IC),C(IC,3),WK(NX,1),
     1                   ERROR
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,ICLAST,ICUBE,IFCT,IFCTL,ILAST,ILM3,IMODE,
     1                   INSERT,INSIRT,INTERV,IORDER,ITRPZ,IUERR,J,JF,
     2                   JL,K,KLOC,KNOT,L,LL,MKN,MKNP2,NXKP
      REAL               BC,COEF,COEFL,ERBUT1,FXDKNT,VORD,VORDL,XI,XIL,
     1                   DS,DUM1,DUM2,DX,HALF,XSCALE,ZERO
      LOGICAL            MODE3
      COMMON /ICSFK1/    COEFL(4,27),XIL(28),VORDL(28,2),BC(30),
     1                   COEF(4,381),XI(381),VORD(30,28,2),FXDKNT,
     2                   ERBUT1,KNOT,INTERV,ILAST,ICLAST,
     3                   IFCTL,IUERR,ITRPZ,ICUBE,IFCT,MKN,
     4                   INSERT,INSIRT(30),IORDER(28),MODE3
      DATA               ZERO/0.0/,HALF/0.5/
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
      IF (MODE .EQ. 4) GO TO 140
      IF (MODE.EQ.3) GO TO 40
      JF = 1
      JL = NXK
      IF (MODE.EQ.1.OR.MODE.EQ.2) GO TO 55
C                                  MODE = 0
      MKN = 28
      ILAST = 0
C                                  WORK AREA ASSIGNMENTS
C                                    WK(I,IFCTL) = FCTL(I)
C                                    WK(I,IUERR) = UERROR(I)
C                                    WK(I,ITRPZ) = TRPZWT(I)
C                                    WK(I,ICUBE) = CUBERR(I)
C                                    WK(I,IFCT+J) = FCT(I,J)
      IFCTL = 1
      IUERR = 2
      ITRPZ = 3
      ICUBE = 4
      IFCT = 4
      MKNP2 = MKN+2
      DO 5 I=1,MKNP2
    5 INSIRT(I) = 0
C                                  INITIAL ERROR
      DO 10 L=1,NX
   10 WK(L,IUERR) = F(L)
C                                  TRAPEZOIDAL WEIGHTS
      XSCALE = HALF/(X(NX)-X(1))
      DO 15 L=3,NX
   15 WK(L-1,ITRPZ) = (X(L)-X(L-2))*XSCALE
      WK(1,ITRPZ) = (X(2)-X(1))*XSCALE
      WK(NX,ITRPZ) = (X(NX)-X(NX-1))*XSCALE
      DO 20 L=2,NX
         IF (X(L) .GE. X(L-1)) GO TO 20
C                                  TERMINAL - X NOT MONOTONE INCREASING
         IER = 131
         GO TO 9000
   20 CONTINUE
C                                  CHECK FOR DISTINCT KNOTS LYING
C                                    BETWEEN THE BOUNDARY KNOTS
      IER = 133
      IF (XK(1) .GE. XK(NXK)) GO TO 9000
      IF (NXK .LT. 3) GO TO 23
      L = NXK-1
      IF (XK(1) .GE. XK(L) .OR. XK(NXK) .LE. XK(L)) GO TO 9000
      IF (NXK .LT. 4) GO TO 23
      LL = L-1
      DO 22 I=2,LL
         IF (XK(1) .GE. XK(I) .OR. XK(NXK) .LE. XK(I)) GO TO 9000
         J = I+1
         DO 21 K=J,L
            IF (XK(I) .EQ. XK(K)) GO TO 9000
   21    CONTINUE
   22 CONTINUE
C                                  CHECK FOR ABSCISSAE BETWEEN BOUNDARY
C                                    KNOTS
   23 IER = 132
      IF (X(1) .LT. XK(1) .OR. X(NX) .GT. XK(NXK)) GO TO 9000
C
      IER = 0
C                                  COMPUTE BASIS FUNCTIONS 1 THROUGH 4
C                                    AND THE B. A. WITH RESPECT TO THEM
C
      DO 30 I=1,4
         CALL ICSFKV (MODE,XK,NXK,X,NX,WK(1,IUERR),WK(1,ITRPZ),
     1   WK(1,IFCTL),WK(1,IFCT+1),IER)
         IF (IER.NE.0) GO TO 9000
         DO 25 L=1,NX
   25    WK(L,IUERR) = WK(L,IUERR)-BC(I)*WK(L,IFCT+I)
   30 CONTINUE
C                                  SAVE ERROR FROM B. A. BY CUBICS
      DO 35 L=1,NX
   35 WK(L,ICUBE) = WK(L,IUERR)
      JL = JL-1
      JF = JF+1
      GO TO 95
   40 CONTINUE
C                                  MODE = 3
C                                  REPLACE THE LAST KNOT INTRODUCED BY
C                                    XK(1) AND RECOMPUTE ERROR
      IF (MODE3) GO TO 50
      MODE3 = .TRUE.
      ERBUT1 = FXDKNT+BC(ILAST)*BC(ILAST)
      DO 45 L=1,NX
   45 WK(L,IUERR) = WK(L,IUERR)+BC(ILAST)*WK(L,IFCT+ILAST)
   50 CONTINUE
         CALL ICSFKV (MODE,XK,NXK,X,NX,WK(1,IUERR),WK(1,ITRPZ),
     1   WK(1,IFCTL),WK(1,IFCT+1),IER)
      IF (IER.NE.0) GO TO 9000
      FXDKNT = ERBUT1-BC(ILAST)*BC(ILAST)
      IF (FXDKNT.LT.ZERO) FXDKNT=ZERO
      ERROR = SQRT(FXDKNT)
      GO TO 9005
   55 CONTINUE
C                                  MODE= 1 OR 2
      IF (NXK.LT.0) GO TO 65
      IF (.NOT.MODE3) GO TO 95
      DO 60 L=1,NX
   60 WK(L,IUERR) = WK(L,IUERR)-BC(ILAST)*WK(L,IFCT+ILAST)
      GO TO 95
   65 CONTINUE
C                                  REMOVE -NXK KNOTS
      DO 70 L=1,NX
   70 WK(L,IUERR) = WK(L,ICUBE)
      KNOT = KNOT+NXK
      INTERV = KNOT-1
      IF (KNOT.LE.2) GO TO 90
      NXKP = IABS(NXK)
      DO 80 I=1,NXKP
         INSERT = INSIRT(ILAST)
         ILM3 = ILAST-3
         DO 75 K=INSERT,ILM3
            IORDER(K) = IORDER(K+1)
            XIL(K) = XIL(K+1)
   75    CONTINUE
         ILAST = ILAST-1
   80 CONTINUE
      ICLAST = (ILAST*(ILAST-7))/2+10
C                                  COMPUTE CURRENT ERROR
      DO 85 I=5,ILAST
      DO 85 L=1,NX
   85 WK(L,IUERR) = WK(L,IUERR)-BC(I)*WK(L,IFCT+I)
      GO TO 110
   90 XIL(2) = XIL(ILAST-2)
      IORDER(2) = 2
      KNOT = 2
      INTERV = 1
      ILAST = 4
      ICLAST = 4
      GO TO 110
   95 CONTINUE
      IF (JF.GT.JL) GO TO 110
C                                  ADD BASIS FUNCTIONS CORRESPONDING TO
C                                    KNOTS XK(JF),XK(JF+1),...,XK(JL)
      IMODE = MAX0(1,MODE)
      DO 105 I=JF,JL
         IF (KNOT.EQ.MKN) GO TO 135
         CALL ICSFKV (IMODE,XK(I),1,X,NX,WK(1,IUERR),WK(1,ITRPZ),
     1   WK(1,IFCTL),WK(1,IFCT+1),IER)
         IF (IER.NE.0) GO TO 9000
C                                  RECOMPUTE ERROR
         DO 100 L=1,NX
  100    WK(L,IUERR) = WK(L,IUERR)-BC(ILAST)*WK(L,IFCT+ILAST)
  105 CONTINUE
  110 CONTINUE
C                                  COMPUTE SCALED L2 ERROR
      FXDKNT = ZERO
      DO 115 L=1,NX
         FXDKNT = FXDKNT+WK(L,IUERR)*WK(L,IUERR)*WK(L,ITRPZ)
  115 CONTINUE
      IF (FXDKNT.LT.ZERO) FXDKNT=ZERO
      ERROR = SQRT(FXDKNT)
C                                  COMPUTE COEFFICIENTS OF B.A.
      DO 125 K=1,KNOT
         KLOC = IORDER(K)
      DO 125 I=1,2
         DS = ZERO
         DO 120 J=1,ILAST
  120    DS = DS+BC(J)*VORD(J,KLOC,I)
         VORDL(K,I) = DS
  125 CONTINUE
      DO 130 I=1,INTERV
         Y(I) = VORDL(I,1)
         C(I,1) = VORDL(I,2)
         DX = XIL(I+1)-XIL(I)
         DUM1 = (VORDL(I+1,1)-VORDL(I,1))/DX
         DUM2 = VORDL(I,2)+VORDL(I+1,2)-DUM1-DUM1
         C(I,2) = (DUM1-DUM2-VORDL(I,2))/DX
         C(I,3) = DUM2/(DX*DX)
  130 CONTINUE
      MODE3 = .FALSE.
      GO TO 9005
  135 CONTINUE
      IER = 130
      GO TO 9000
  140 CONTINUE
      NXK = KNOT
      DO 145 I=1,KNOT
         XK(I) = XIL(I)
  145 CONTINUE
      GO TO 9005
 9000 CONTINUE
      CALL UERTST (IER,'ICSFKU')
 9005 RETURN
      END
