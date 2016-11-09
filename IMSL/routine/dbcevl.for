C   IMSL ROUTINE NAME   - DBCEVL
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1982
C
C   PURPOSE             - BICUBIC SPLINE MIXED PARTIAL DERIVATIVE
C                           EVALUATOR
C
C   USAGE               - CALL DBCEVL (X,NX,Y,NY,C,IC,XL,YL,PDS,IER)
C
C   ARGUMENTS    X      - VECTOR OF LENGTH NX. (INPUT) X MUST BE
C                           ORDERED SO THAT X(I) .LT. X(I+1) FOR
C                           I=1,...,NX-1.
C                NX     - NUMBER OF ELEMENTS IN X. (INPUT) NX MUST BE
C                           .GE. 2.
C                Y      - VECTOR OF LENGTH NY. (INPUT) Y MUST BE
C                           ORDERED SO THAT Y(J) .LT. Y(J+1) FOR
C                           J=1,...,NY-1.
C                NY     - NUMBER OF ELEMENTS IN Y. (INPUT) NY MUST BE
C                           .GE. 2.
C                         NOTE - THE COORDINATE PAIRS (X(I),Y(J)), FOR
C                           I=1,...,NX AND J=1,...,NY, GIVE THE POINTS
C                           WHERE THE FUNCTION VALUES ARE DEFINED.
C                C      - ARRAY OF SPLINE COEFFICIENTS. (INPUT)
C                           C IS OF DIMENSION 2 BY NX BY 2 BY NY.
C                           THE SPLINE COEFFICIENTS CAN BE COMPUTED BY
C                           IMSL SUBROUTINE IBCCCU.
C                           (NOTE - C IS TREATED INTERNALLY AS A
C                            2 BY NX BY 2*NY ARRAY BECAUSE CERTAIN
C                            ENVIRONMENTS DO NOT PERMIT QUADRUPLY-
C                            DIMENSIONED ARRAYS.  IN THESE
C                            ENVIRONMENTS THE CALLING PROGRAM MAY
C                            DIMENSION C IN THE SAME MANNER.)
C                IC     - SECOND DIMENSION OF ARRAY C EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT
C                           (INPUT).  IC MUST BE .GE. NX.
C                XL,YL  - (XL,YL) IS THE POINT AT WHICH THE MIXED
C                           PARTIAL DERIVATIVES OF THE SPLINE ARE TO BE
C                           EVALUATED. (INPUT)
C                PDS    - VECTOR OF LENGTH 6 CONTAINING THE PARTIAL
C                           DERIVATIVES OF THE BICUBIC SPLINE, S(X,Y),
C                           EVALUATED AT X=XL AND Y=YL. (OUTPUT)
C                             PDS(1) = S(XL,YL)
C                             PDS(2) = DS/DX
C                             PDS(3) = DS/DY
C                             PDS(4) = D(DS/DX)/DY
C                             PDS(5) = D(DS/DX)/DX
C                             PDS(6) = D(DS/DY)/DY.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         WARNING ERROR
C                           IER = 33, XL IS LESS THAN X(1).
C                           IER = 34, YL IS LESS THAN Y(1).
C                           IER = 35, XL IS GREATER THAN X(NX).
C                           IER = 36, YL IS GREATER THAN Y(NY).
C
C   REQD. IMSL ROUTINES - UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE DBCEVL (X,NX,Y,NY,C,IC,XL,YL,PDS,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            NX,NY,IC,IER
      REAL               X(1),Y(1),C(2,IC,1),XL,YL,PDS(6)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,J,K,KM1,KP1,KP2,LXPL,LX,LY,L,LXP1
      REAL               HX,HY,SUX(2),SUY(2),SU(2),SVX(2),SV(2),SXY(2),
     *                   U,V,SPLN0,SPLN1,SPLN2,S0,SH,SP0,SPH,H,D
      SPLN0(S0,SH,SP0,SPH,H,D) = S0+D*(H*SP0+D*(3.*(SH-S0)-
     * (SPH+2.*SP0)*H+D*(2.*(S0-SH)+(SPH+SP0)*H)))
      SPLN1(S0,SH,SP0,SPH,H,D) = SP0+D*(6.*(SH-S0)/H-2.*
     * (SPH+2.*SP0)+3.*D*(2.*(S0-SH)/H+(SPH+SP0)))
      SPLN2(S0,SH,SP0,SPH,H,D) = 6.*(SH-S0)/H**2-2.*
     * (SPH+2.*SP0)/H+D*(2.*(S0-SH)/H**2+(SPH+SP0)/H)*6.
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
      IF (XL.LT.X(1)) IER = 33
      DO 5 I=2,NX
         LX = I-1
         IF (XL.LE.X(I)) GO TO 10
    5 CONTINUE
      IER = 35
   10 IF (YL.LT.Y(1)) IER = 34
      DO 15 J=2,NY
         LY = J-1
         IF (YL.LE.Y(J)) GO TO 20
   15 CONTINUE
      IER = 36
   20 LXP1 = LX+1
      HX = X(LXP1)-X(LX)
      HY = Y(LY+1)-Y(LY)
      U = (XL-X(LX))/HX
      V = (YL-Y(LY))/HY
      K = 2*LY
      KP1 = K+1
      KP2 = K+2
      KM1 = K-1
      DO 25 L=1,2
         LXPL = LX-1+L
         I = 2*(LY-1+L)
         J = I-1
         SUX(L) = SPLN1(C(1,LX,J),C(1,LXP1,J),C(2,LX,J),
     *   C(2,LXP1,J),HX,U)
         SXY(L) = SPLN1(C(1,LX,I),C(1,LXP1,I),C(2,LX,I),
     *   C(2,LXP1,I),HX,U)
         SU(L) = SPLN0(C(1,LX,J),C(1,LXP1,J),C(2,LX,J),
     *   C(2,LXP1,J),HX,U)
         SUY(L) = SPLN0(C(1,LX,I),C(1,LXP1,I),C(2,LX,I),
     *   C(2,LXP1,I),HX,U)
         SV(L) = SPLN0(C(1,LXPL,KM1),C(1,LXPL,KP1),C(1,LXPL,K),
     *   C(1,LXPL,KP2),HY,V)
         SVX(L) = SPLN0(C(2,LXPL,KM1),C(2,LXPL,KP1),C(2,LXPL,K),
     *   C(2,LXPL,KP2),HY,V)
   25 CONTINUE
      PDS(1) = SPLN0(SV(1),SV(2),SVX(1),SVX(2),HX,U)
      PDS(2) = SPLN1(SV(1),SV(2),SVX(1),SVX(2),HX,U)
      PDS(3) = SPLN1(SU(1),SU(2),SUY(1),SUY(2),HY,V)
      PDS(4) = SPLN1(SUX(1),SUX(2),SXY(1),SXY(2),HY,V)
      PDS(5) = SPLN2(SV(1),SV(2),SVX(1),SVX(2),HX,U)
      PDS(6) = SPLN2(SU(1),SU(2),SUY(1),SUY(2),HY,V)
      IF (IER.GT.0) CALL UERTST(IER,'DBCEVL')
      RETURN
      END
