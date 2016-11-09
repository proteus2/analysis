C   IMSL ROUTINE NAME   - BECTR
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1982
C
C   PURPOSE             - TETRACHORIC CORRELATION COEFFICIENT ESTIMATION
C
C   USAGE               - CALL BECTR (N,U,V,HU,HV,IOPT,R,RS,K,IER)
C
C   ARGUMENTS    N      - SAMPLE SIZE. (INPUT)
C                U      - INPUT VECTOR CONTAINING DATA PERTINENT TO THE
C                           FIRST VARIABLE. IF IOPT.EQ.0, U HAS LENGTH
C                           N. IF IOPT.NE.0, U HAS LENGTH 2.
C                V      - INPUT VECTOR CONTAINING DATA PERTINENT TO THE
C                           SECOND VARIABLE. IF IOPT.EQ.0, V HAS LENGTH
C                           N. IF IOPT.NE.0, V HAS LENGTH 2.
C                HU     - INPUT CONSTANT USED TO CATEGORIZE VALUES OF U.
C                           IF ANY VALUE OF U IS EQUAL TO OR GREATER
C                           THAN HU, IT WILL BE CLASSIFIED IN THE
C                           HIGHER CATEGORY.  HU IS USED ONLY IF IOPT
C                           IS EQUAL TO 0.
C                HV     - INPUT CONSTANT USED TO CATEGORIZE VALUES OF V.
C                           SEE HU.
C                IOPT   - OPTION INDICATING WHETHER U AND V ARE TO BE
C                           CATEGORIZED. (INPUT)  IN A TWO BY TWO TABLE
C                           A,B,C,D, LET A AND C CONSTITUTE ROW 1 WHILE
C                           B AND D CONSTITUTE ROW 2.  A WILL CONTAIN
C                           COUNTS OF ALL OBSERVATIONS SUCH THAT U(I)
C                           IS GREATER THAN OR EQUAL TO HU, AND V(I)
C                           IS GREATER THAN OR EQUAL TO HV.  SIMILARLY,
C                              B -- U(I).LT.HU, V(I).GE.HV,
C                              C -- U(I).GE.HU, V(I).LT.HV,
C                              D -- U(I).LT.HU, V(I).LT.HV.
C                           IF IOPT = 0, THE DATA STREAM U,V IS
C                           CATEGORIZED AS NOTED ABOVE.  IF IOPT ON
C                           ENTRY IS NON-ZERO, THEN, ON ENTRY TO BECTR,
C                              U(1) MUST CONTAIN COUNT A,
C                              U(2) MUST CONTAIN COUNT B,
C                              V(1) MUST CONTAIN COUNT C, AND
C                              V(2) MUST CONTAIN COUNT D.
C                R      - VECTOR OF LENGTH 7 CONTAINING K ROOTS (ESTI-
C                           MATES OF RHO) IN THE RANGE (-1.,1.)
C                           (OUTPUT). R(K+1),...,R(7) ARE SET TO MACHINE
C                           INFINITY.
C                RS     - STANDARD ERROR OF THE ESTIMATE OF RHO.
C                           (OUTPUT)
C                K      - NUMBER OF REAL ROOTS (ESTIMATES OF RHO)
C                           FOUND IN (-1.,1.). (OUTPUT)
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER = 129 INDICATES THAT AN ERROR OCCURRED
C                             IN IMSL SUBROUTINE MDNRIS..
C                           IER = 130 INDICATES THAT ONE OR MORE THAN
C                             ONE CELL OF THE TWO BY TWO TABLE IS ZERO.
C                           IER = 131 INDICATES THAT AN ERROR OCCURRED
C                             IN THE IMSL SUBROUTINE ZRPOLY.
C                           IN THE ABOVE ERROR SITUATIONS, R AND RS
C                             ARE SET TO MACHINE INFINITY.
C                         WARNING ERROR
C                           IER = 33 INDICATES THAT MULTIPLE ROOTS OR
C                             NO ROOTS WERE FOUND IN THE INTERVAL
C                             (-1.,1.).  THAT IS, K IS NOT EQUAL TO
C                             ONE.  R(I),I=1,K,...,7 IS SET TO
C                             INFINITY.
C
C   REQD. IMSL ROUTINES - MDNRIS,MERFI,UERTST,UGETIO,ZRPOLY,ZRPQLB,
C                           ZRPQLC,ZRPQLD,ZRPQLE,ZRPQLF,ZRPQLG,ZRPQLH,
C                           ZRPQLI
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE BECTR  (N,U,V,HU,HV,IOPT,R,RS,K,IER)
C
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N,IOPT,K,IER
      REAL               U(1),V(1),HU,HV,R(7),RS
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,IC,IR,J,M
      REAL               QD1,QD2,Q1(2),XD1,XD2,X1(2)
      REAL               A,B,C,D,FN,ONE,REPS,RINFP,SE,TEMP
      REAL               EX1,EX2,EX3,EX4,EX5,SIDSQR,ZERO
      REAL               XD12,XD22,Y1
      REAL               XCOF(8),ZZ(14),COUNT(2,2)
      EQUIVALENCE        (A,COUNT(1,1)),(B,COUNT(2,1)),
     1                   (C,COUNT(1,2)),(D,COUNT(2,2)),(X1(1),XD1),
     2                   (X1(2),XD2),(Q1(1),QD1),(Q1(2),QD2)
      DATA               ZERO/0.0/,ONE/1.0/
      DATA               RINFP/Z7FFFFFFF/
      DATA               SE/2.718282/
      DATA               SIDSQR/.1591549/
      DATA               EX5/.1666667/
      DATA               EX4/4.166667E-2/
      DATA               EX3/8.333333E-3/
      DATA               EX2/1.388889E-3/
      DATA               EX1/1.984127E-4/
      DATA               REPS/Z3C100000/
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
      K = 1
      RS = RINFP
      R(1) = 0.0
      DO 3 I=2,7
    3 R(I) = RINFP
      A = ZERO
      B = ZERO
      C = ZERO
      D = ZERO
      M = 7
      IF(IOPT.NE.0) GO TO 10
C                                  CATEGORIZE DATA (U,V) INTO CELLS
      DO 5 I=1,N
         IR = 2
         IF(U(I).GE.HU) IR = 1
         IC = 2
         IF(V(I).GE.HV) IC = 1
         COUNT(IR,IC) = COUNT(IR,IC) + ONE
    5 CONTINUE
      GO TO 15
C                                  ASSIGN PRE-CATEGORIZED DATA TO CELLS
   10 A = U(1)
      B = U(2)
      C = V(1)
      D = V(2)
   15 TEMP = A*B*C*D
C                                  TEST IF ANY CELL EQUALS ZERO
      IF(TEMP.NE.ZERO) GO TO 20
      IER = 130
      GO TO 60
C                                  FIND STANDARD NORMAL DEVIATES AT QD1
C                                  AND QD2 AND THE ORDINATES AT THOSE
C                                  POINTS
   20 FN = ONE/N
      QD1 =  (B+D)*FN
      QD2 =  (C+D)*FN
      DO 25 I=1,2
         CALL MDNRIS(Q1(I),X1(I),IER)
         IF(IER.NE.0) GO TO 60
   25 CONTINUE
      XD12 = XD1*XD1
      XD22 = XD2*XD2
      Y1 = SIDSQR*SE**(-.5*(XD12+XD22))
      TEMP = XD1*XD2
C                                  COMPUTE THE TETRACHORIC CORRELATION
C                                  COEFFICIENT
      XCOF(8) = (B*C-A*D)*FN*FN/Y1
      XCOF(7) = ONE
      XCOF(6) = TEMP*.5
      XCOF(5) = (XD12-ONE)*(XD22-ONE)*EX5
      XCOF(4) = (TEMP*(XD12-3.0)*(XD22-3.0)*EX4)
      XCOF(3) = (XD12*(XD12-6.0)+3.0)*(XD22*(XD22-6.0)+3.0)*EX3
      XCOF(2) = TEMP*(XD12*(XD12-10.0)+15.0)*(XD22*(XD22-10.0)+15.0)*
     1EX2
      XCOF(1) = (((XD12-15.0)*XD12+45.0)*XD12-15.0)*(((XD22-15.0)*XD22
     1+45.0)*XD22-15.0)*EX1
C                                  FIND THE ROOTS (JENKINS-TRAUB)
      CALL ZRPOLY(XCOF,M,ZZ,IER)
      K = 0
      IF(IER.EQ.0) GO TO 40
      IER = 131
      GO TO 60
C                                  FIND A REAL ROOT
   40 DO 45 I=1,13,2
         J = I+1
         TEMP = ABS(ZZ(I))
         IF(ABS(ZZ(J)).GT.REPS*AMAX1(TEMP,.1)) GO TO 45
C                                  DETERMINE WHICH REAL ROOTS ARE IN
C                                  THE (-1,1) RANGE
         IF(TEMP.GT.ONE) GO TO 45
         R(K+1) = ZZ(I)
         K = K+1
   45 CONTINUE
C                                  COMPUTE THE STANDARD ERROR OF RHO
   50 RS = FN*SQRT((A+C)*(A+B)*FN*QD1*QD2)/Y1
   55 IF(K.NE.1) IER = 33
      IF(IER.EQ.0) GO TO 9005
   60 IF(IER.GE.128) R(1) = RINFP
 9000 CONTINUE
      CALL UERTST(IER,'BECTR ')
 9005 RETURN
      END
