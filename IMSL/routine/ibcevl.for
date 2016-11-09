C   IMSL ROUTINE NAME   - IBCEVL
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1982
C
C   PURPOSE             - EVALUATION OF A BICUBIC SPLINE
C
C   USAGE               - CALL IBCEVL (X,NX,Y,NY,C,IC,XL,YL,FL,IER)
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
C                XL,YL  - (XL,YL) IS THE POINT AT WHICH THE SPLINE IS
C                           TO BE EVALUATED. (INPUT)
C                FL     - THE VALUE OF THE SPLINE APPROXIMATION AT
C                           (XL,YL). (OUTPUT)
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         WARNING ERROR
C                           IER = 33, XL IS LESS THAN X(1).
C                           IER = 34, YL IS LESS THAN Y(1).
C                           IER = 35, XL IS GREATER THAN X(NX).
C                           IER = 36, YL IS GREATER THAN Y(NY).
C
C
C   REQD. IMSL ROUTINES - UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE IBCEVL (X,NX,Y,NY,C,IC,XL,YL,FL,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            NX,NY,IC,IER
      REAL               X(1),Y(1),C(2,IC,1),XL,YL,FL
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,J,LX,LYPL,LY,L
      REAL               HX,HY,SUY(2),SU(2),U,V,SPLN,S0,SH,SP0,SPH,H,D
      SPLN(S0,SH,SP0,SPH,H,D) = S0+D*(H*SP0+D*(3.*(SH-S0)-
     * (SPH+2.*SP0)*H+D*(2.*(S0-SH)+(SPH+SP0)*H)))
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
   20 I = LX+1
      HX = X(I)-X(LX)
      HY = Y(LY+1)-Y(LY)
      U = (XL-X(LX))/HX
      V = (YL-Y(LY))/HY
      DO 25 L=1,2
         LYPL = LY-1+L
         J = 2*LYPL-1
         SU(L) = SPLN(C(1,LX,J),C(1,I,J),C(2,LX,J),C(2,I,J),HX,U)
         J = 2*LYPL
         SUY(L) = SPLN(C(1,LX,J),C(1,I,J),C(2,LX,J),C(2,I,J),HX,U)
   25 CONTINUE
      FL = SPLN(SU(1),SU(2),SUY(1),SUY(2),HY,V)
      IF (IER.GT.0) CALL UERTST(IER,'IBCEVL')
      RETURN
      END
