C   IMSL ROUTINE NAME   - ICSCCU
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1980
C
C   PURPOSE             - CUBIC SPLINE INTERPOLATION
C                           (EASY-TO-USE VERSION)
C
C   USAGE               - CALL ICSCCU (X,Y,NX,C,IC,IER)
C
C   ARGUMENTS    X      - VECTOR OF LENGTH NX CONTAINING THE ABSCISSAE
C                           OF THE NX DATA POINTS (X(I),Y(I)) I=1,...,
C                           NX. (INPUT) X MUST BE ORDERED SO THAT
C                           X(I) .LT. X(I+1).
C                Y      - VECTOR OF LENGTH NX CONTAINING THE ORDINATES
C                           (OR FUNCTION VALUES) OF THE NX DATA POINTS.
C                           (INPUT)
C                NX     - NUMBER OF ELEMENTS IN X AND Y. (INPUT) NX
C                           MUST BE .GE. 2.
C                C      - SPLINE COEFFICIENTS. (OUTPUT) C IS AN NX-1 BY
C                           3 MATRIX. THE VALUE OF THE SPLINE
C                           APPROXIMATION AT T IS
C                           S(T) = ((C(I,3)*D+C(I,2))*D+C(I,1))*D+Y(I)
C                           WHERE X(I) .LE. T .LT. X(I+1) AND
C                           D = T-X(I).
C                IC     - ROW DIMENSION OF MATRIX C EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN
C                           THE CALLING PROGRAM. (INPUT)
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER = 129, IC IS LESS THAN NX-1.
C                           IER = 130, NX IS LESS THAN 2.
C                           IER = 131, INPUT ABSCISSA ARE NOT ORDERED
C                             SO THAT X(1) .LT. X(2) ... .LT. X(NX).
C
C   REQD. IMSL ROUTINES - UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE ICSCCU (X,Y,NX,C,IC,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            NX,IC,IER
      REAL               X(NX),Y(NX),C(IC,3)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            IM1,I,JJ,J,MM1,MP1,M,NM1,NM2
      REAL               DIVDF1,DIVDF3,DTAU,G,CNX(3)
C                                  FIRST EXECUTABLE STATEMENT
      NM1 = NX-1
      IER = 129
      IF (IC .LT. NM1) GO TO 9000
      IER = 130
      IF (NX .LT. 2) GO TO 9000
      IER = 131
      IF (NX .EQ. 2) GO TO 45
C                                  COMPUTE NOT-A-KNOT SPLINE
      DO 5 M = 2,NM1
         MM1=M-1
         C(M,2) = X(M)-X(MM1)
         IF (C(M,2).LE.0.0) GO TO 9000
         C(M,3) = (Y(M)-Y(MM1))/C(M,2)
    5 CONTINUE
      CNX(2) = X(NX)-X(NM1)
      IF (CNX(2).LE.0.0) GO TO 9000
      CNX(3) = (Y(NX)-Y(NM1))/CNX(2)
      IER = 0
      NM2 = NX-2
      IF (NX .GT. 3) GO TO 10
      C(1,3) = CNX(2)
      C(1,2) = C(2,2)+CNX(2)
      C(1,1) = ((C(2,2)+2.*C(1,2))*C(2,3)*CNX(2)+C(2,2)**2*CNX(3))
     1/C(1,2)
      GO TO 20
   10 C(1,3) = C(3,2)
      C(1,2) = C(2,2)+C(3,2)
      C(1,1) = ((C(2,2)+2.*C(1,2))*C(2,3)*C(3,2)+C(2,2)**2*C(3,3))
     1/C(1,2)
      DO 15 M=2,NM2
         MP1=M+1
         MM1=M-1
         G = -C(MP1,2)/C(MM1,3)
         C(M,1) = G*C(MM1,1)+3.*C(M,2)*C(MP1,3)+3.*C(MP1,2)*C(M,3)
         C(M,3) = G*C(MM1,2)+2.*C(M,2)+2.*C(MP1,2)
   15 CONTINUE
   20 G = -CNX(2)/C(NM2,3)
      C(NM1,1) = G*C(NM2,1)+3.*C(NM1,2)*CNX(3)+3.*CNX(2)*C(NM1,3)
      C(NM1,3) = G*C(NM2,2)+2.*C(NM1,2)+2.*CNX(2)
      IF (NX.GT.3) GO TO 25
      CNX(1)=2.*CNX(3)
      CNX(3)=1.
      G=-1./C(NM1,3)
      GO TO 30
   25 G = C(NM1,2)+CNX(2)
      CNX(1) = ((CNX(2)+2.*G)*CNX(3)*C(NM1,2)+CNX(2)**2*
     1(Y(NM1)-Y(NX-2))/C(NM1,2))/G
      G = -G/C(NM1,3)
      CNX(3) = C(NM1,2)
   30 CNX(3) = G*C(NM1,2)+CNX(3)
      CNX(1) = (G*C(NM1,1)+CNX(1))/CNX(3)
      C(NM1,1) = (C(NM1,1)-C(NM1,2)*CNX(1))/C(NM1,3)
      DO 35 JJ=1,NM2
         J = NM1-JJ
         C(J,1) = (C(J,1)-C(J,2)*C(J+1,1))/C(J,3)
   35 CONTINUE
      DO 40 I=2,NM1
         IM1 = I-1
         DTAU = C(I,2)
         DIVDF1 = (Y(I)-Y(IM1))/DTAU
         DIVDF3 = C(IM1,1)+C(I,1)-2.*DIVDF1
         C(IM1,2) = (DIVDF1-C(IM1,1)-DIVDF3)/DTAU
         C(IM1,3) = DIVDF3/DTAU**2
   40 CONTINUE
      DTAU = CNX(2)
      DIVDF1 = (Y(NX)-Y(NM1))/DTAU
      DIVDF3 = C(NM1,1)+CNX(1)-2.*DIVDF1
      C(NM1,2) = (DIVDF1-C(NM1,1)-DIVDF3)/DTAU
      C(NM1,3) = DIVDF3/DTAU**2
      GO TO 9005
   45 IF (X(1) .GE. X(2)) GO TO 9000
      IER = 0
      C(1,1) = (Y(2)-Y(1))/(X(2)-X(1))
      C(1,2) = 0.0
      C(1,3) = 0.0
      GO TO 9005
 9000 CONTINUE
      CALL UERTST(IER,'ICSCCU')
 9005 RETURN
      END
