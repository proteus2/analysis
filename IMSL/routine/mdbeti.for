C   IMSL ROUTINE NAME   - MDBETI
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1980
C
C   PURPOSE             - INVERSE BETA PROBABILITY DISTRIBUTION FUNCTION
C
C   USAGE               - CALL MDBETI (P,A,B,X,IER)
C
C   ARGUMENTS    P      - INPUT PROBABILITY IN THE EXCLUSIVE RANGE
C                           (0,1)
C                A      - INPUT FIRST PARAMETER OF THE BETA DISTRIBUTION
C                B      - INPUT SECOND PARAMETER OF THE BETA DISTRI-
C                           BUTION
C                X      - OUTPUT VALUE SUCH THAT P IS THE PROBABILITY
C                           THAT A RANDOM VARIABLE DISTRIBUTED BETA(A,B)
C                           IS LESS THAN OR EQUAL TO X.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER = 129 INDICATES THAT P WAS NOT IN THE
C                             EXCLUSIVE RANGE (0,1).
C                           IER = 130 INDICATES A AND/OR B LESS THAN
C                             OR EQUAL TO 0.
C                           IER = 131 INDICATES THAT THE VALUE X COULD
C                             NOT BE FOUND WITHIN 100 ITERATIONS OF
C                             NEWTONS METHOD.
C
C   REQD. IMSL ROUTINES - H32/MDBETA,MLGAMD=DLGAMA,UERTST,UGETIO
C                       - H36,H48,H60/MDBETA,MLGAMA=ALGAMA,UERTST,
C                           UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C
C-----------------------------------------------------------------------
C
      SUBROUTINE MDBETI (P,A,B,X,IER)
C
      DOUBLE PRECISION   DLGAMA
      DATA               EPS/.0001/,SIG/1.E-5/,ITMAX/30/,ZERO/0.0/
      DATA               SMEXE/-170.0/
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
      XL = AMIN1(A,B)
      IF (XL.LE.0.0) GO TO 70
      IF (XL.LE.1.0) GO TO 20
      XR = AMAX1(A,B)
      IF (10.0*XL.LE.XR) GO TO 20
      IC = 0
      XL = 0.0
      XR = 1.0
      FXL = -P
      FXR = 1.0-P
      IF (FXL*FXR.GT.ZERO) GO TO 65
C                                  BISECTION METHOD
    5 X = (XL+XR)*.5
      CALL MDBETA (X,A,B,P1,IER)
      IF (IER.NE.0) GO TO 20
      FCS = P1-P
      IF (FCS*FXL.GT.ZERO) GO TO 10
      XR = X
      FXR = FCS
      GO TO 15
   10 XL = X
      FXL = FCS
   15 XRMXL = XR-XL
      IF (XRMXL.LE.SIG.AND.ABS(FCS).LE.EPS) GO TO 9005
      IC = IC+1
      IF (IC.LE.ITMAX) GO TO 5
C                                  ERROR RETURNED FROM MDBETA
C                                  USE NEWTONS METHOD FOR SKEWED CASES
   20 IF (P.LE.0.0.OR.P.GE.1.0) GO TO 65
      IF (P.GT..5) GO TO 25
      AA = A
      BB = B
      Q0 = ALOG(P)
      GO TO 30
   25 Q0 = ALOG(1.0-P)
      AA = B
      BB = A
   30 XT = AA/(AA+BB)
      DTEMP = DLGAMA(DBLE(AA+BB))-DLGAMA(DBLE(AA))-DLGAMA(DBLE(BB))
      DTEMP = DTEMP-(AA+BB)*ALOG(AA+BB)+(AA-.5)*ALOG(AA)+(BB-.5)
     1*ALOG(BB)
      DTEMP = DTEMP+.5*ALOG(BB/AA)+AA*ALOG(1.0+BB/AA)+BB*ALOG(1.+AA/BB)
      DO 45 NC=1,100
         TEMP = ALOG(15.0+AA+BB)
         FN = 0.7*TEMP*TEMP+AMAX1(XT*(AA+BB)-AA,0.0)
         TEMP = AA+FN+FN
         AFN = AINT(FN)+1.0
         C = 1.0-(AA+BB)*XT/TEMP
         ZI = 2.0/(C+SQRT(C*C-4.0*FN*(FN-BB)*XT/(TEMP*TEMP)))
   35    AFN = AFN-1.0
         IF (AFN.LT..5) GO TO 40
         TEMP = AA+AFN+AFN
         ZI = (TEMP-2.0)*(TEMP-1.0-AFN*(AFN-BB)*XT*ZI/TEMP)
         TEMP = AA+AFN-1.0
         ZI = 1.0/(1.0-TEMP*(TEMP+BB)*XT/ZI)
         GO TO 35
   40    ZZ = ZI
         TEMP = ALOG(XT)
         IF (TEMP.LE.SMEXE) GO TO 50
         QX = DTEMP+AA*TEMP+BB*ALOG(1.-XT)+ALOG(ZZ)
         XC = (Q0-QX)*(1.-XT)*ZZ/AA
         XC = AMAX1(XC,-.99)
         TEMP = .5/XT-.5
         XC = AMIN1(XC,TEMP)
         XT = XT*(1.+XC)
         IF (ABS(XC).LT.SIG) GO TO 55
   45 CONTINUE
      IER = 131
      GO TO 9000
   50 XT = 0.0
   55 IF (P.GT..5) GO TO 60
      X = XT
      GO TO 9005
   60 X = 1.0-XT
      GO TO 9005
   65 IER = 129
      GO TO 9000
   70 IER = 130
 9000 CONTINUE
      CALL UERTST (IER,6HMDBETI)
 9005 RETURN
      END
