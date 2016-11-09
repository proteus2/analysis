C   IMSL ROUTINE NAME   - ICSPLN
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1980
C
C   PURPOSE             - CUBIC SPLINE INTERPOLATION WITH
C                           PERIODIC END CONDITIONS
C
C   USAGE               - CALL ICSPLN (X,Y,NX,C,IC,WK,IER)
C
C   ARGUMENTS    X      - VECTOR OF LENGTH NX CONTAINING THE
C                           ABSCISSAE OF THE DATA POINTS
C                           (X(I),Y(I)),I=1...NX. (INPUT) OFTEN
C                           THE POLAR COORDINATE THETA IN
C                           APPLICATIONS. X MUST BE ORDERED
C                           SO THAT X(I).LT.X(I+1).
C                Y      - VECTOR OF LENGTH NX CONTAINING THE
C                           ORDINATES (OR FUNCTION VALUES) OF
C                           THE DATA POINTS. (INPUT) OFTEN THE
C                           POLAR COORDINATE R IN APPLICATIONS.
C                           THE ROUTINE ASSUMES THAT Y(1) IS EQUAL TO
C                           Y(NX) AND A WARNING ERROR IS GIVEN IF THIS
C                           CONDITION DOES NOT EXIST.
C                NX     - THE NUMBER OF ELEMENTS IN X AND Y. (INPUT)
C                           NX MUST BE GREATER THAN OR EQUAL TO 4.
C                C      - SPLINE COEFFICIENTS. (OUTPUT)C IS AN
C                           NX-1 BY 3 MATRIX. THE VALUE OF THE SPLINE
C                           APPROXIMATION AT T IS
C                           S(T) = ((C(I,3)*D+C(I,2))*D+C(I,1))*D+Y(I)
C                           WHERE X(I) .LE. T .LT. X(I+1) AND
C                           D = T-X(I).
C                IC     - ROW DIMENSION OF MATRIX C EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN
C                           THE CALLING PROGRAM. (INPUT)
C                WK     - WORK ARRAY OF LENGTH 6*NX.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         WARNING ERROR
C                           IER = 33, Y(1) IS NOT EQUAL TO Y(NX).
C                         TERMINAL ERROR
C                           IER = 129, IC IS LESS THAN NX-1.
C                           IER = 130, NX IS LESS THAN 4.
C                           IER = 131, INPUT ABSCISSAE ARE NOT ORDERED
C                             SO THAT X(1) .LT. X(2) ... .LT. X(NX).
C
C   REQD. IMSL ROUTINES - UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS  1.  THE CUBIC SPLINE CALCULATED PASSES THROUGH EACH DATA
C                POINT AND IS PERIODIC. THAT IS, ITS VALUE AND ITS FIRST
C                AND SECOND DERIVATIVES ARE EQUAL AT THE POINTS X(1)
C                AND X(NX).
C            2.  A WARNING ERROR IS GIVEN IF Y(1) IS NOT EQUAL TO Y(NX).
C                Y(1) IS USED IN THE COMPUTATIONS.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE ICSPLN (X,Y,NX,C,IC,WK,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            NX,IC,IER
      REAL               X(1),Y(1),C(IC,3),WK(6,1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            II,IM1,IMJ,IP1,IP2,I,J,N,NM1
      REAL               A2,A4,B2,B3,B4,C2,C3,C4,D1,D2,D3,D4,DET1,DET2,
     1                   DET3,DET4,HH(6),H(7),R,HT1,HT2,HT3,HT4
C                                  FIRST EXECUTABLE STATEMENT
      IER = 129
      N = NX-1
      IF (IC .LT. N) GO TO 9000
      IER = 130
      IF (NX .LT. 4) GO TO 9000
      NM1 = 0
      IF (Y(1) .NE. Y(NX)) NM1 = 33
C                                  CALCULATE COEFFICIENTS OF MATRIX A
C                                    WHERE A(I,J)=JTH SPLINE BASIS
C                                    FUNCTION EVALUATED AT X(I)
      IER = 131
      DO 10 I=1,N
         DO 5 J=1,4
            IMJ = MOD(N+I+J-4,N)+1
            H(J) = X(IMJ+1)-X(IMJ)
            IF (H(J).LE.0.0) GO TO 9000
    5    CONTINUE
         HT1 = H(3)+H(4)
         HT2 = H(1)+H(2)
         HT3 = H(2)+H(3)
         DET1 = (HT1+H(3))*H(2)*HT2+(HT2+H(2))*H(3)
     1   *HT1
         WK(1,I) = 0.0
         A2 = H(1)**2*(HT3+H(4))*HT3/HT2/DET1
         WK(2,I) = A2
         WK(3,I) = 1.0
         A4 = H(4)**2*HT3*(H(1)+HT3)/HT1/DET1
         WK(4,I) = A4
         WK(5,I) = 0.0
         WK(6,I) = Y(I)
   10 CONTINUE
      IER = NM1
      NM1 = N-1
      WK(1,1) = WK(4,N)
      WK(5,1) = WK(2,1)
      WK(1,NM1) = WK(2,N)
      WK(5,NM1) = WK(4,NM1)
      WK(1,N) = 1.0
C                                  SOLVE LINEAR SYSTEM FOR COEFFICIENTS
C                                    OF SPLINE BASIS FUNCTIONS
      DO 15 I=1,NM1
         R = WK(5,I)/WK(3,I)
         WK(1,N) = WK(1,N)-R*WK(1,I)
         WK(6,N) = WK(6,N)-R*WK(6,I)
         IF (I.EQ.NM1) GO TO 15
         IP1=I+1
         WK(5,IP1) = WK(5,IP1)-R*WK(2,IP1)
         R = WK(4,I)/WK(3,I)
         WK(3,IP1) = WK(3,IP1)-R*WK(2,IP1)
         WK(1,IP1) = WK(1,IP1)-R*WK(1,I)
         WK(6,IP1) = WK(6,IP1)-R*WK(6,I)
   15 CONTINUE
      WK(6,N) = WK(6,N)/WK(1,N)
      WK(6,NM1) = (WK(6,NM1)-WK(1,NM1)*WK(6,N))/WK(3,NM1)
      DO 20 II=2,NM1
         I = N-II
         IP1=I+1
         WK(6,I) = (WK(6,I)-WK(2,IP1)*WK(6,IP1)-WK(1,I)*WK(6,N))/WK(3,I)
   20 CONTINUE
C                                  IN EACH INTERVAL CALCULATE
C                                    COEFFICIENTS C(I,J),J=1,3 FROM
C                                    BASIS FUNCTION COEFFICIENTS
      DO 35 I=1,N
         IM1 = MOD(N+I-2,N)+1
         IP1 = MOD(N+I,N)+1
         IP2 = MOD(N+I+1,N)+1
         DO 25 J=1,7
            IMJ = MOD(N+I+J-5,N)+1
            H(J) = X(IMJ+1)-X(IMJ)
   25    CONTINUE
         DO 30 J=1,6
            HH(J) = H(J)+H(J+1)
   30    CONTINUE
         HT1 = H(4)+HH(4)
         HT2 = H(3)+HH(2)
         HT3 = H(5)+HH(5)
         HT4 = H(4)+HH(3)
         DET1 = (H(3)+HH(3))*H(2)*HH(1)+(H(2)+HH(1))*H(3)*HH(3)
         DET2 = HT1*H(3)*HH(2)+HT2*H(4)*HH(4)
         DET3 = HT3*H(4)*HH(3)+HT4*H(5)*HH(5)
         DET4 = (H(6)+HH(6))*H(5)*HH(4)+(H(5)+HH(4))*H(6)*HH(6)
         D1 = HH(5)*(HH(5)+H(7))/H(4)/HH(4)/DET4
         B2 = 3.*H(3)*HH(4)*(HH(4)+H(6))/HH(3)/DET3
         C2 = B2/H(3)
         D2 = -(HT3*HT4+3.*HH(3)*H(4)+H(3)*H(3)+H(5)
     1   *HH(5))/HH(3)/H(4)/DET3
         B3 = 3.*(H(4)*HH(4)-H(3)*HH(2))/DET2
         C3 = -3.*(HT1+HT2)/DET2
         D3 = (HT1*HT2+3.*H(4)*HH(4)+H(5)*H(5)+H(3)
     1   *HH(2))/HH(4)/H(4)/DET2
         B4 = -3.*H(4)*HH(2)*(HH(1)+H(3))/HH(3)/DET1
         C4 = -B4/H(4)
         D4 = -C4/H(4)/3.0
         C(I,1) = WK(6,IM1)*B4+WK(6,I)*B3+WK(6,IP1)*B2
         C(I,2) = WK(6,IM1)*C4+WK(6,I)*C3+WK(6,IP1)*C2
         C(I,3) = WK(6,IM1)*D4+WK(6,I)*D3+WK(6,IP1)*D2+WK(6,IP2)*D1
   35 CONTINUE
 9000 CONTINUE
      IF (IER .NE. 0) CALL UERTST(IER,'ICSPLN')
 9005 RETURN
      END
