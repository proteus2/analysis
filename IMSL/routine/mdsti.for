C   IMSL ROUTINE NAME   - MDSTI
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1981
C
C   PURPOSE             - INVERSE OF A MODIFICATION OF STUDENTS T
C                           PROBABILITY DISTRIBUTION FUNCTION
C
C   USAGE               - CALL MDSTI (Q,F,X,IER)
C
C   ARGUMENTS    Q      - INPUT PROBABILITY IN THE EXCLUSIVE RANGE
C                           (0,1). (THE SUM OF THE AREAS (EQUAL) IN
C                           BOTH TAILS OF THE T DISTRIBUTION.)
C                F      - INPUT DEGREES OF FREEDOM FOR T DISTRIBUTION.
C                           F MUST BE GREATER THAN ZERO.
C                X      - OUTPUT VALUE SUCH THAT THE PROBABILITY OF THE
C                           ABSOLUTE VALUE OF T BEING GREATER THAN X IS
C                           Q.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         WARNING ERROR
C                           IER = 35 INDICATES OVERFLOW WOULD HAVE
C                             OCCURRED. X IS SET TO MACHINE INFINITY.
C                         TERMINAL ERROR
C                           IER = 129 INDICATES THAT F (DEGREES OF
C                             FREEDOM) IS LESS THAN OR EQUAL TO ZERO.
C                           IER = 130 INDICATES THAT Q IS OUT OF THE
C                             EXCLUSIVE RANGE (0,1).
C                           IER = 131 INDICATES THAT AN ERROR OCCURRED
C                             IN IMSL ROUTINE MDNRIS, THE INVERSE NORMAL
C                             PROBABILITY DISTRIBUTION FUNCTION.
C                           IER = 132 INDICATES THAT CONVERGENCE WAS
C                             NOT ACHIEVED. THIS CAN ONLY OCCUR FOR
C                             F LESS THAN 2.0 .
C   REQD. IMSL ROUTINES - MDNRIS,MERFI,UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS      NOTE THAT MDSTI DOES NOT PROVIDE THE ACTUAL T INVERSE.
C                FOR P EQUAL TO THE PROBABILITY THAT A STUDENTS T
C                RANDOM VARIABLE IS LESS THAN X, THAT INVERSE CAN BE
C                OBTAINED BY THE FOLLOWING RULES.
C                  A.  FOR P IN THE RANGE (0.0,0.5), CALL MDSTI WITH
C                      Q = 2*P AND NEGATE THE RESULT X.
C                  B.  FOR P IN THE RANGE (0.5,1.0), CALL MDSTI WITH
C                      Q = 2*(1-P).
C
C-----------------------------------------------------------------------
C
      SUBROUTINE MDSTI  (Q,F,X,IER)
C
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            IER
      REAL               F,Q,X
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            NC
      REAL               A,AA,AFN,B,BB,C,CON1,CON2,CON3,CON4,CON5,CON6,
     1                   CON7,CON8,D,DTEMP,FN,HPI,P,PP,Q0,QX,SIG,SMEXE,
     2                   SDELP,TEMP,TWOPI,XC,XINF,XT,XX,Y,ZI,ZZ
      DATA               SDELP/.9999999E0/
      DATA               XINF/.7237005E+76/
      DATA               SMEXE/-170.0/
      DATA               SIG/1.E-5/
      DATA               HPI/1.570796/
      DATA               TWOPI/6.283185/
      DATA               CON1/.8333333E-1/
      DATA               CON2/.3333333E-1/
      DATA               CON3/.2523810E0/
      DATA               CON4/.5256065E0/
      DATA               CON5/.1011523E01/
      DATA               CON6/.1517474E01/
      DATA               CON7/.2269489E01/
      DATA               CON8/.1534264E0/
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
      IF (F.GT.0.0) GO TO 5
C                                  TERMINAL ERROR - F OUT OF RANGE
      IER = 129
      GO TO 9000
    5 IF (Q.LE.0.0.OR.Q.GE.1.0) GO TO 90
C                                  EXACT INTEGRAL FOR 2 DEGREES OF
C                                  FREEDOM
      IF (ABS(F-2.0).GT..000001) GO TO 10
      X = SQRT(2.0/(Q*(2.0-Q))-2.0)
      GO TO 9005
C                                  EXACT INTEGRAL FOR 1 DEGREE OF
C                                  FREEDOM
   10 IF (ABS(F-1.0).GT..000001) GO TO 15
      A = Q*HPI
      X = COS(A)/SIN(A)
      GO TO 9005
   15 IF (F.GT.2.0) GO TO 75
C                                  F IS LESS THAN 2.0
      Y = .5*F
      D = Y
      XX = -2.0+(D+.5)*ALOG(1.0+1.0/D)+(D+1.5)*ALOG(1.0+1.0/(D+1.0))
      A = D+2.0
      C = CON7/A
      C = CON4/(A+CON5/(A+CON6/(A+C)))
      C = CON1/(A+CON2/(A+CON3/(A+C)))
      B = C+XX
      D = D+.5
      IF(D.GE.1.0) GO TO 20
      XX = -2.0+(D+.5)*ALOG(1.0+1.0/D)+(D+1.5)*ALOG(1.0+1.0/(D+1.0))
      A = D+2.0
      GO TO 25
   20 XX = -1.0+(D+.5)*ALOG(1.0+1.0/D)
      A = D+1.0
   25 C = CON7/A
      C = CON4/(A+CON5/(A+CON6/(A+C)))
      C = CON1/(A+CON2/(A+CON3/(A+C)))
      C = C+XX
      A = CON8
      P = 1.0 - Q
      PP = P
      AA = .5
      BB = Y
      IF (P.LE..5) GO TO 30
      PP = Q
      AA = Y
      BB = .5
   30 Q0 = ALOG(PP)
      XT = AA/(AA+BB)
      DTEMP = C-A-B-.5*ALOG((AA+BB)*TWOPI)
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
      IER = 132
      GO TO 9000
   50 XT = 0.0
   55 X = XT
      IF (P.GT..5) GO TO 60
      IF (X.GE.SDELP) GO TO 70
      X = F*X/(1.-X)
      GO TO 65
   60 IF (X .EQ. 0.0) GO TO 70
      X = F/X - F
   65 X = SQRT(X)
      GO TO 9005
   70 X = XINF
      IER = 35
      GO TO 9000
C                                  EXPANSION FOR F GREATER THAN 2
   75 A = 1.0/(F-0.5)
      B = 48.0/(A*A)
      C = ((20700.*A/B-98.)*A-16.)*A+96.36
      D = ((94.5/(B+C)-3.0)/B+1.0)*SQRT(A*HPI)*F
      XX = D*Q
      Y = XX**(2.0/F)
      IF (Y.GT.A+.05) GO TO 85
      Y = ((1.0/(((F+6.0)/(F*Y)-0.089*D-0.822)*(F+2.0)*3.0)+0.5/(F
     1+4.0))*Y-1.0)*(F+1.0)/(F+2.0)+1.0/Y
   80 X = SQRT(F*Y)
      GO TO 9005
C                                  ASYMPTOTIC INVERSE EXPANSION ABOUT
C                                  NORMAL
   85 X = .5*Q
      CALL MDNRIS (X,XX,IER)
      IF (IER.NE.0) GO TO 95
      Y = XX*XX
      IF (F.LT.5.) C = C+0.3*(F-4.5)*(XX+0.6)
      C = (((.05*D*XX-5.0)*XX-7.0)*XX-2.0)*XX+B+C
      Y = (((((0.4*Y+6.3)*Y+36.)*Y+94.5)/C-Y-3.0)/B+1.0)*XX
      Y = A*Y*Y
      D = Y
      IF (D.LE..002) Y = .5*Y*Y+Y
      IF (D.GT..002) Y = EXP(Y)-1.0
      GO TO 80
C                                  Q IS OUT OF RANGE
   90 IER = 130
      GO TO 9000
C                                  ERROR OCCURRED IN SUBROUTINE MDNRIS
   95 IER = 131
 9000 CONTINUE
      CALL UERTST (IER,6HMDSTI )
 9005 RETURN
      END
