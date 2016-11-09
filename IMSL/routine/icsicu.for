C   IMSL ROUTINE NAME   - ICSICU
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - INTERPOLATORY APPROXIMATION BY CUBIC SPLINES
C                           WITH ARBITRARY SECOND DERIVATIVE END
C                           CONDITIONS.
C
C   USAGE               - CALL ICSICU (X,Y,NX,BPAR,C,IC,IER)
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
C                BPAR   - VECTOR OF LENGTH 4 CONTAINING THE END
C                           CONDITION PARAMETERS. (INPUT)
C                           2.0*SPP(1)+BPAR(1)*SPP(2) = BPAR(2),
C                           BPAR(3)*SPP(NX-1)+2.0*SPP(NX) = BPAR(4),
C                           WHERE SPP(I) = SECOND DERIVATIVE OF THE
C                           CUBIC SPLINE FUNCTION S EVALUATED AT X(I).
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
      SUBROUTINE ICSICU (X,Y,NX,BPAR,C,IC,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            NX,IC,IER
      REAL               X(NX),Y(NX),BPAR(4),C(IC,3)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,J,NXM1
      REAL               DX,DXJ,DXJP1,DXP,DYJ,DYJP1,HALF,ONE,PJ,
     1                   SIX,SIXI,TWO,YPPA,YPPB,ZERO
      EQUIVALENCE        (DXJ,YPPB),(PJ,SIXI),(DXJP1,YPPA)
      DATA               ZERO/0.0/,HALF/0.5/,ONE/1.0/,
     1                   TWO/2.0/,SIX/6.0/
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
C                                  CHECK ERROR CONDITIONS
      NXM1 = NX-1
      IF (IC .LT. NXM1) GO TO 30
      IF (NX .LT. 2) GO TO 35
      IF (NX .EQ. 2) GO TO 10
C                                  COMPUTE COEFFICIENTS AND RIGHT
C                                  HAND SIDE OF THE TRIDIAGONAL
C                                  SYSTEM DEFINING THE SECOND
C                                  DERIVATIVES OF THE SPLINE
C                                  INTERPOLANT FOR (X,Y)
C                                  C(J,1) = LAMBDA(J)
C                                  C(J,2) = MU(J)
C                                  C(J,3) = D(J)
      DXJ = X(2)-X(1)
      IF (DXJ .LE. ZERO) GO TO 40
      DYJ = Y(2)-Y(1)
      DO 5 J=2,NXM1
         DXJP1 = X(J+1)-X(J)
         IF (DXJP1 .LE. ZERO) GO TO 40
         DYJP1 = Y(J+1)-Y(J)
         DXP = DXJ+DXJP1
         C(J,1) = DXJP1/DXP
         C(J,2) = ONE-C(J,1)
         C(J,3) = SIX*(DYJP1/DXJP1-DYJ/DXJ)/DXP
         DXJ = DXJP1
         DYJ = DYJP1
    5 CONTINUE
C                                  FACTOR THE TRIDIAGONAL MATRIX
C                                  AND SOLVE FOR U
C                                  C(J,2) = U(J)
C                                  C(J,1) = Q(J)
C                                  BPAR(1) = LAMBDA(1)
C                                  BPAR(2) = D(1)
C                                  BPAR(3) = MU(NX)
C                                  BPAR(4) = D(NX)
   10 C(1,1) = -BPAR(1)*HALF
      C(1,2) = BPAR(2)*HALF
      IF (NX .EQ. 2) GO TO 20
      DO 15 J=2,NXM1
         PJ = C(J,2)*C(J-1,1)+TWO
         C(J,1) = -C(J,1)/PJ
         C(J,2) = (C(J,3)-C(J,2)*C(J-1,2))/PJ
   15 CONTINUE
C                                  SOLVE FOR CUBIC COEFFICIENTS
C                                  OF SPLINE INTERPOLANT
C                                  C(J,1), C(J,2), AND C(J,3)
   20 YPPB = (BPAR(4)-BPAR(3)*C(NXM1,2))/(BPAR(3)*C(NXM1,1)+TWO)
      SIXI = ONE/SIX
      DO 25 I=1,NXM1
         J = NX-I
         YPPA = C(J,1)*YPPB+C(J,2)
         DX = X(J+1)-X(J)
         C(J,3) = SIXI*(YPPB-YPPA)/DX
         C(J,2) = HALF*YPPA
         C(J,1) = (Y(J+1)-Y(J))/DX-(C(J,2)+C(J,3)*DX)*DX
         YPPB = YPPA
   25 CONTINUE
      GO TO 9005
   30 IER = 129
      GO TO 9000
   35 IER = 130
      GO TO 9000
   40 IER = 131
 9000 CONTINUE
      CALL UERTST(IER,'ICSICU')
 9005 RETURN
      END
