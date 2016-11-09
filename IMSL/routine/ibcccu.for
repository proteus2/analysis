C   IMSL ROUTINE NAME   - IBCCCU
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1982
C
C   PURPOSE             - BICUBIC SPLINE TWO-DIMENSIONAL COEFFICIENT
C                           CALCULATOR
C
C   USAGE               - CALL IBCCCU (F,X,NX,Y,NY,C,IC,WK,IER)
C
C   ARGUMENTS    F      - NX BY NY MATRIX CONTAINING THE FUNCTION
C                           VALUES. (INPUT) F(I,J) IS THE FUNCTION VALUE
C                           AT THE POINT (X(I),Y(J)) FOR I=1,...,NX AND
C                           J=1,...,NY.
C                X      - VECTOR OF LENGTH NX. (INPUT) X MUST BE
C                           ORDERED SO THAT X(I) .LT. X(I+1) FOR
C                           I=1,...,NX-1.
C                NX     - NUMBER OF ELEMENTS IN X. (INPUT) NX MUST BE
C                           .GE. 4.
C                Y      - VECTOR OF LENGTH NY. (INPUT) Y MUST BE
C                           ORDERED SO THAT Y(J) .LT. Y(J+1) FOR
C                           J=1,...,NY-1.
C                NY     - NUMBER OF ELEMENTS IN Y. (INPUT) NY MUST BE
C                           .GE. 4.
C                         NOTE - THE COORDINATE PAIRS (X(I),Y(J)), FOR
C                           I=1,...,NX AND J=1,...,NY, GIVE THE POINTS
C                           WHERE THE FUNCTION VALUES F(I,J) ARE
C                           DEFINED.
C                C      - ARRAY OF SPLINE COEFFICIENTS. (OUTPUT)
C                           C IS OF DIMENSION 2 BY NX BY 2 BY NY.
C                           AT THE POINT (X(I),Y(J))
C                             C(1,I,1,J) = S
C                             C(2,I,1,J) = DS/DX
C                             C(1,I,2,J) = DS/DY
C                             C(2,I,2,J) = D(DS/DX)/DY
C                           WHERE S(X,Y) IS THE SPLINE APPROXIMATION.
C                           (NOTE - C IS TREATED INTERNALLY AS A
C                             2 BY NX BY 2*NY ARRAY BECAUSE CERTAIN
C                             ENVIRONMENTS DO NOT PERMIT QUADRUPLY-
C                             DIMENSIONED ARRAYS.  IN THESE
C                             ENVIRONMENTS THE CALLING PROGRAM MAY
C                             DIMENSION C IN THE SAME MANNER.)
C                IC     - ROW DIMENSION OF MATRIX F AND SECOND
C                           DIMENSION OF ARRAY C EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT.
C                           (INPUT). IC MUST BE .GE. NX.
C                WK     - WORK VECTOR OF LENGTH
C                           2*NX*NY+2*MAX(NX,NY)
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER = 129, IC IS LESS THAN NX
C                           IER = 130, NX IS LESS THAN 4
C                           IER = 131, NY IS LESS THAN 4
C                           IER = 132, X OR Y ARE NOT ORDERED SO THAT
C                             X(I) .LT. X(I+1) AND
C                             Y(I) .LT. Y(I+1)
C
C
C   REQD. IMSL ROUTINES - IBCDCU,UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE IBCCCU (F,X,NX,Y,NY,C,IC,WK,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            NX,NY,IC,IER
      REAL               F(IC,1),X(1),Y(1),C(2,IC,1),WK(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            IWK
C                                  FIRST EXECUTABLE STATEMENT
      IER = 129
      IF (IC .LT. NX) GO TO 9000
      IER = 130
      IF (NX .LT. 4) GO TO 9000
      IER = 131
      IF (NY .LT. 4) GO TO 9000
      IWK = 2*NY*NX
      CALL IBCDCU(X,F,NX,NY,WK(IWK+1),WK,IC,NY,IER)
      IF (IER .GT. 0) GO TO 9000
      CALL IBCDCU(Y,WK,NY,2*NX,WK(IWK+1),C,NY,2*IC,IER)
      IF (IER .EQ. 0) GO TO 9005
 9000 CONTINUE
      CALL UERTST(IER,'IBCCCU')
 9005 RETURN
      END
