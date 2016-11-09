C   IMSL ROUTINE NAME   - ICSSCV
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1980
C
C   PURPOSE             - CUBIC SPLINE DATA SMOOTHER
C                           (EASY-TO-USE VERSION)
C
C   USAGE               - CALL ICSSCV (X,F,NX,Y,C,IC,IJOB,WK,IER)
C
C   ARGUMENTS    X      - VECTOR OF LENGTH NX CONTAINING THE
C                           ABSCISSAE OF THE NX DATA POINTS
C                           (X(I),F(I)) I=1..NX. (INPUT) X
C                           MUST BE ORDERED SO THAT
C                           X(I) .LT. X(I+1).
C                F      - VECTOR OF LENGTH NX CONTAINING THE
C                           ORDINATES (OR FUNCTION VALUES)
C                           OF THE NX DATA POINTS (INPUT).
C                NX     - NUMBER OF DATA POINTS (INPUT).
C                           NX MUST BE .GE. 4.
C                Y,C    - SPLINE COEFFICIENTS. (OUTPUT) Y
C                           IS A VECTOR OF LENGTH NX. C IS
C                           AN NX-1 BY 3 MATRIX. THE VALUE
C                           OF THE SPLINE APPROXIMATION AT T IS
C                           S(T)=((C(I,3)*D+C(I,2))*D+C(I,1))*D+Y(I)
C                           WHERE X(I).LE.T.LT.X(I+1) AND
C                           D = T-X(I).
C                IC     - ROW DIMENSION OF MATRIX C EXACTLY
C                           AS SPECIFIED IN THE DIMENSION
C                           STATEMENT IN THE CALLING PROGRAM. (INPUT)
C                IJOB   - JOB SELECTION PARAMETER. (INPUT)
C                         IJOB = 1 SHOULD BE SELECTED WHEN
C                           NX IS SMALL (LESS THAN ABOUT 20)
C                           OR WHEN UNEQUALLY SPACED ABSCISSAE
C                           (X(1),X(2),...) ARE USED.
C                         IJOB = 2 SHOULD BE SELECTED WHEN
C                           NX IS LARGE AND THE ABSCISSAE ARE
C                           EQUALLY SPACED.
C                WK     - WORK VECTOR OF LENGTH
C                           NX*(3*NX+5) IF IJOB=1,
C                           OR
C                           6*NX IF IJOB=2.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER = 129, IC IS LESS THAN NX-1.
C                           IER = 130, NX IS LESS THAN 4.
C                           IER = 131, INPUT ABSCISSAE ARE NOT
C                             ORDERED SO THAT X(I).LT.X(I+1).
C                           IER = 132, FAILURE TO CONVERGE ON
C                             A MINIMUM OF THE CROSS-VALIDATION
C                             FUNCTION.
C                           IER = 133, IJOB = 2 AND THE
C                             ABSCISSAE ARE NOT EQUALLY SPACED.
C
C   REQD. IMSL ROUTINES - EQRT2S,ICSSD,ICSSE,ICSSF,ICSSG,ICSSH,
C                           LEQT1B,LSVDB,LSVDF,LSVG1,LSVG2,VHS12,
C                           UERSET,UERTST,UGETIO,ZXGSP
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS      THE NUMBER OF ARITHMETIC OPERATIONS REQUIRED BY
C                THE ROUTINE IS PROPORTIONAL TO
C                NX**3 WHEN IJOB IS EQUAL TO 1 AND
C                NX    WHEN IJOB IS EQUAL TO 2.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE ICSSCV (X,F,NX,Y,C,IC,IJOB,WK,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            NX,IC,IJOB,IER
      REAL               X(1),F(1),Y(1),C(IC,3),WK(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I1,I2,I3,IWK1,IWK2,IWK3,IWK4,IWK5,IWK6,IWK7,
     1                   IWK8,I,NXL1,NXL2,IM1,IP1
      REAL               CONS,DIF,HMAX,HMIN,H,PI,RO,VP,VSTORE(10)
      DATA               PI/3.141593/
      DATA               VSTORE/5.61873E-3,-2.4733E-4,1.0679E-5,
     1                   -4.62E-7,2.0E-8,-1.0E-9,4*0.0/
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
      IF (IC.GE.NX-1) GO TO 5
      IER = 129
      GO TO 9000
    5 CONTINUE
      IF (NX.GE.4) GO TO 10
      IER = 130
      GO TO 9000
   10 CONTINUE
      NXL1 = NX-1
      IF (IJOB.EQ.2) GO TO 20
C                                  IJOB = 1 SELECTED
      IWK1 = 1
      IM1 = NX*NX
      IWK2 = IWK1+IM1
      IWK3 = IWK2+IM1
      IWK4 = IWK3+IM1
      IWK5 = IWK4+NX
      IWK6 = IWK5+NX
      IWK7 = IWK6+NX
      IWK8 = IWK7+NX
      CALL ICSSD (NX,NX,X,F,Y,WK(IWK4),WK(IWK5),WK(IWK6),WK(IWK7),
     1 WK(IWK8),WK(IWK1),WK(IWK2),WK(IWK3),IER)
      DO 15 I=1,NXL1
         IM1 = I-1
         I1 = IWK4+IM1
         C(I,1) = WK(I1)
         I2 = IWK5+IM1
         C(I,2) = WK(I2)
         I3 = IWK6+IM1
         C(I,3) = WK(I3)
   15 CONTINUE
      GO TO 9000
C                                  IJOB = 2 SELECTED
   20 HMAX = X(2)-X(1)
      HMIN = HMAX
      DO 25 I=2,NXL1
         H = X(I+1)-X(I)
         IF (H.GT.HMAX) HMAX = H
         IF (H.LT.HMIN) HMIN = H
   25 CONTINUE
      IF (ABS(HMIN-HMAX).LT.0.01*HMAX) GO TO 30
      IER = 133
      GO TO 9000
   30 CONTINUE
      H = (HMIN+HMAX)/2.0
      HMIN = NX
      CONS = PI/HMIN/H
      NXL2 = NX-2
      DO 35 I=1,NXL2
         HMIN = I
         VP = HMIN+0.5
         IF (I.LE.10) VP = VP+VSTORE(I)
         VP = CONS*VP*CONS*VP
         WK(NX+I) = VP*VP*H
   35 CONTINUE
      CALL ICSSG (F,Y,WK,WK(NX+1),NX,H,RO,WK(2*NX+1))
      DO 40 I=1,NX
         Y(I) = F(I)
   40 CONTINUE
      CALL ICSSF (NX,RO,H,Y,WK,WK(2*NX+1))
      H = 1./H
      HMIN = H*H
      DO 45 I=1,NXL1
         IP1 = I+1
         DIF = (Y(IP1)-Y(I))*H
         C(I,1) = WK(I)
         C(I,2) = (3.0*DIF-2.0*WK(I)-WK(IP1))*H
         C(I,3) = (-2.0*DIF+WK(I)+WK(IP1))*HMIN
   45 CONTINUE
 9000 CONTINUE
      IF (IER .NE. 0) CALL UERTST (IER,6HICSSCV)
 9005 RETURN
      END
