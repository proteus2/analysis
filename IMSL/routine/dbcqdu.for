C   IMSL ROUTINE NAME   - DBCQDU
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1980
C
C   PURPOSE             - BICUBIC SPLINE QUADRATURE
C
C   USAGE               - CALL DBCQDU (F,IFD,X,NX,Y,NY,A,B,C,D,Q,WK,IER)
C
C   ARGUMENTS    F      - NX BY NY MATRIX CONTAINING THE FUNCTION
C                           VALUES. (INPUT) F(I,J) IS THE FUNCTION VALUE
C                           AT THE POINT (X(I),Y(J)) FOR I=1,...,NX AND
C                           J=1,...,NY.
C                IFD    - ROW DIMENSION OF THE MATRIX F EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT
C                           IN THE CALLING PROGRAM. (INPUT)
C                X      - VECTOR OF LENGTH NX. (INPUT) X MUST BE
C                           ORDERED SO THAT X(I) .LT. X(I+1) FOR
C                           I=1,...,NX-1.
C                NX     - NUMBER OF ELEMENTS IN X. (INPUT) NX MUST BE
C                           .GE. 2.
C                Y      - VECTOR OF LENGTH NY. (INPUT) Y MUST BE
C                           ORDERED SO THAT Y(J) .LT. Y(J+1) FOR
C                           J=1,...,NY-1.
C                NY     - NUMBER OF ELEMENTS IN Y. (INPUT) NY MUST BE
C                           .GE. 2.
C                         NOTE - THE COORDINATE PAIRS (X(I),Y(J)),FOR
C                           I=1,...,NX AND J=1,...,NY, GIVE THE POINTS
C                           WHERE THE FUNCTION VALUES F(I,J) ARE
C                           DEFINED.
C                A,B    - X-DIRECTION LIMITS OF INTEGRATION. (INPUT)
C                C,D    - Y-DIRECTION LIMITS OF INTEGRATION. (INPUT)
C                Q      - DOUBLE INTEGRAL FROM A TO B AND C TO D.
C                           (OUTPUT)
C                WK     - WORK VECTOR OF LENGTH
C                           (NY+5)*NX+NY-1+MAX(5*NX-4,5*NY-2)
C                           (SEE PROGRAMMING NOTES FOR FURTHER DETAILS)
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER = 129, IFD IS LESS THAN NX.
C                           IER = 130, NX IS LESS THAN 2.
C                           IER = 131, NY IS LESS THAN 2.
C                         WARNING ERROR
C                           IER = 36, A AND/OR B IS LESS THAN X(1).
C                           IER = 37, A AND/OR B IS GREATER THAN X(NX).
C                           IER = 38, C AND/OR D IS LESS THAN Y(1).
C                           IER = 39, C AND/OR D IS GREATER THAN Y(NY).
C
C   REQD. IMSL ROUTINES - ICSEVU,ICSCCU,UERSET,UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS      WHEN THE LIMITS OF INTEGRATION ARE OUTSIDE OF THE
C                RECTANGLE (X(1),X(NX)) X (Y(1),Y(NY)), THE
C                BICUBIC SPLINE IS EXTENDED USING THE BOUNDARY PIECES.
C                INTEGRATION IS PERFORMED OVER THE EXTENDED SPLINE.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE DBCQDU (F,IFD,X,NX,Y,NY,A,B,C,D,Q,WK,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            IFD,NX,NY,IER
      REAL               F(IFD,NY),X(NX),Y(NY),A,B,C,D,Q,WK(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,IA,IAP1,IB,IC,ICP1,ID,IDM1,II,IP1,IPT,IV,
     1                   IYL,J,JER,JJER,KER,LER,MER,NCOEF,NER,NFI,NFI1,
     2                   NFL,NFTEMP,NFTMP1,NN,NWTXX,NXM1,NXX,NYL,NYLS,
     3                   NYLSP,NYM1,LEVEL,LEVOLD
      REAL               AA,BB,CC,DD,DIST1,DIST2,FOUR,
     1                   FRDTRD,HALF,SUM,THIRD,V,WTYY,ZERO
      DATA               ZERO/0.0/,HALF/.5/,FOUR/4.0/
      DATA               THIRD/.3333333/,
     1                   FRDTRD/1.333333/
C                                  FIRST EXECUTABLE STATEMENT
      JER = 0
      KER = 0
      LER = 0
      MER = 0
      NER = 129
      IF (IFD.LT.NX) GO TO 140
      NER = 130
      IF (NX.LT.2) GO TO 140
      NER = 131
      IF (NY.LT.2) GO TO 140
      NER = 0
      LEVEL = 2
      CALL UERSET (LEVEL,LEVOLD)
      NXM1 = NX-1
      NYM1 = NY-1
C                                  FIND THE X INTERVALS
      IA = 1
      IPT = 1
      V = AMIN1(A,B)
    5 DIST1 = V-X(IA)
      DO 10 I=IA,NXM1
         IV = I
         DIST2 = V-X(I+1)
         IF (DIST2.LT.ZERO) GO TO 15
         IF (I.LT.NXM1) DIST1 = DIST2
   10 CONTINUE
      IV = NXM1
C                                  IF V .GT. X(NX) - WARNING
      IF (DIST2.GT.ZERO) KER = 37
   15 CONTINUE
C                                  CHECK FOR V .LT. X(1)
      IF (DIST1.LT.ZERO) JER = 36
      IF (IPT.EQ.2) GO TO 20
      IPT = 2
      IA = IV
      V = AMAX1(A,B)
      GO TO 5
   20 IB = IV
C                                  FIND THE Y INTERVALS
      IC = 1
      IPT = 1
      V = AMIN1(C,D)
   25 DIST1 = V-Y(IC)
      DO 30 I=IC,NYM1
         IV = I
         DIST2 = V-Y(I+1)
         IF (DIST2.LT.ZERO) GO TO 35
         IF (I.LT.NYM1) DIST1 = DIST2
   30 CONTINUE
      IV = NYM1
C                                  IF V .GT. Y(NY) - WARNING
      IF (DIST2.GT.ZERO) MER = 39
   35 CONTINUE
C                                  CHECK FOR V .LT. Y(1)
      IF (DIST1.LT.ZERO) LER = 38
      IF (IPT.EQ.2) GO TO 40
      IPT = 2
      IC = IV
      V = AMAX1(C,D)
      GO TO 25
   40 ID = IV
      AA = AMIN1(A,B)
      BB = AMAX1(A,B)
C                                  DEFINE XX(I),I=1,...,NXX
      NXX = 2*(IB-IA)+3
      WK(1) = AA
      WK(NXX) = BB
      IF (NXX.NE.3) GO TO 45
      WK(2) = (AA+BB)*HALF
      GO TO 55
   45 IAP1 = IA+1
      WK(2) = (AA+X(IAP1))*HALF
      WK(3) = X(IAP1)
      WK(NXX-1) = (X(IB)+BB)*HALF
      IF (NXX.EQ.5) GO TO 55
      NWTXX = NXX-3
      II = IAP1
      DO 50 I=4,NWTXX,2
         WK(I) = (X(II)+X(II+1))*HALF
         II = II+1
         WK(I+1) = X(II)
   50 CONTINUE
C                                  COMPUTE WTXX(I),I=1,...,NXX
   55 WK(NXX+1) = (WK(3)-WK(2))*THIRD
      WK(NXX+2) = FOUR*WK(NXX+1)
      WK(NXX+NXX) = (WK(NXX)-WK(NXX-1))*THIRD
      IF (NXX.EQ.3) GO TO 65
      NWTXX = NXX-2
      DO 60 I=3,NWTXX,2
         IP1 = I+1
         WK(NXX+I) = (WK(IP1)-WK(I-1))*THIRD
         WK(NXX+IP1) = FRDTRD*(WK(IP1)-WK(I))
   60 CONTINUE
   65 CC = AMIN1(C,D)
      DD = AMAX1(C,D)
C                                  DEFINE YL(I),I=1,...,NYL
      NYL = ID-IC+3
C                                  NYLS AND NYLSP ARE THE STARTING AND
C                                    ENDING LOCATIONS OF YL
      NYLS = NXX+NXX+1
      NYLSP = NYLS+NYL-1
      WK(NYLS) = CC
      WK(NYLSP) = DD
      IF (NYL.NE.3) GO TO 70
      WK(NYLS+1) = (CC+DD)*HALF
      GO TO 80
   70 ICP1 = IC+1
      WK(NYLS+1) = (CC+Y(ICP1))*HALF
      WK(NYLSP-1) = (Y(ID)+DD)*HALF
      IF (NYL.EQ.4) GO TO 80
      IDM1 = ID-1
      IYL = NYLS+2
      DO 75 I=ICP1,IDM1
         WK(IYL) = (Y(I)+Y(I+1))*HALF
         IYL = IYL+1
   75 CONTINUE
C                                  INTERPOLATE AT YL(I),J=1,...,NYL FOR
C                                    EACH X(I),I=1,...,NX. NFL, NFI,
C                                    NFCOEF, AND NFTEMP ARE THE STARTING
C                                    LOCATIONS OF FL, FI, COEF, AND
C                                    TEMP.
   80 NFL = NYLS+NYL
      NFI = NFL+NX*NYL
      NFI1 = NFI-1
      NCOEF = NFI+NY
      NFTEMP = NCOEF+NYM1*3
      NFTMP1 = NFTEMP+NYL-1
      DO 95 I=1,NX
C                                  MOVE F(I,J),J=1,...,NY TO FI(J)
         DO 85 J=1,NY
            WK(NFI1+J) = F(I,J)
   85    CONTINUE
C                                  INTERPOLATE
         CALL ICSCCU (Y,WK(NFI),NY,WK(NCOEF),NYM1,JJER)
C                                  EVALUATE
         CALL ICSEVU (Y,WK(NFI),NY,WK(NCOEF),NYM1,WK(NYLS),WK(NFTEMP),
     1   NYL,JJER)
C                                  MOVE FTEMP(J),J=1,...,NYL TO FL(I,J)
         NN = NYLSP+I
         DO 90 J=NFTEMP,NFTMP1
            WK(NN) = WK(J)
            NN = NN+NX
   90    CONTINUE
   95 CONTINUE
C                                  INTEGRATE
      Q = ZERO
C                                  RECOMPUTE NCOEF AND NFTEMP TO SAVE
C                                    STORAGE
      NCOEF = NFL+NX*NYL
      NFTEMP = NCOEF+NXM1*3
      NFTMP1 = NFTEMP-1-NXX
C                                  NWTXX AND NN POINT TO THE STARTING
C                                    AND ENDING LOCATIONS OF WTXX
      NWTXX = NXX+1
      NN = NXX+NXX
      NFI = NFL+NX
      NFI1 = NFL+(NYL-1)*NX
C                                  INTERPOLATE AT XX(I),I=1,NXX FOR EACH
C                                    YY WHERE YY IS THE UNION OF THE
C                                    SETS YL AND Y(J), WHERE
C                                    IC .LT. J .LT. ID
C                                    YL(1), YL(2), YL(NYL) ARE SPECIAL
C                                    CASES
C                                  COMPUTE THE WEIGHT OF YL(1)
      WTYY = (WK(NYLS+1)-WK(NYLS))*THIRD
      I = NFL
  100 SUM = ZERO
C                                  INTERPOLATE
      CALL ICSCCU (X,WK(I),NX,WK(NCOEF),NXM1,JJER)
C                                  EVALUATE
      CALL ICSEVU (X,WK(I),NX,WK(NCOEF),NXM1,WK(1),WK(NFTEMP),NXX,JJER)
C                                  COMPUTE THE SUM FOR THE INTERVAL
      DO 105 J=NWTXX,NN
         SUM = SUM+WK(J)*WK(NFTMP1+J)
  105 CONTINUE
C                                  ACCUMULATE THE INTEGRAL
      Q = Q+WTYY*SUM
      IF (I.EQ.NFI1) GO TO 115
      IF (I.EQ.NFI) GO TO 110
C                                  COMPUTE THE WEIGHT FOR YL(2)
      WTYY = FOUR*WTYY
      I = NFI
      GO TO 100
C                                  COMPUTE THE WEIGHT FOR YL(NYL)
  110 WTYY = (WK(NYLSP)-WK(NYLSP-1))*THIRD
      I = NFI1
      GO TO 100
  115 IF (NYL.EQ.3) GO TO 135
      II = ICP1
      NFL = NFI
      NFI = NYLS-1
      NFI1 = NYL-2
      DO 130 I=2,NFI1
         SUM = ZERO
C                                  INTERPOLATE
         CALL ICSCCU (X,F(1,II),NX,WK(NCOEF),NXM1,JJER)
C                                  EVALUATE
         CALL ICSEVU (X,F(1,II),NX,WK(NCOEF),NXM1,WK(1),WK(NFTEMP),NXX,
     1   JJER)
C                                  COMPUTE THE SUM FOR THE INTERVAL
         DO 120 J=NWTXX,NN
            SUM = SUM+WK(J)*WK(NFTMP1+J)
  120    CONTINUE
C                                  COMPUTE THE WEIGHT FOR Y(II)
         WTYY = (WK(NYLS+I)-WK(NFI+I))*THIRD
C                                  ACCUMULATE THE INTEGRAL
         Q = Q+WTYY*SUM
         SUM = ZERO
C                                  INTERPOLATE
         NFL = NFL+NX
         CALL ICSCCU (X,WK(NFL),NX,WK(NCOEF),NXM1,JJER)
C                                  EVALUATE
         CALL ICSEVU (X,WK(NFL),NX,WK(NCOEF),NXM1,WK(1),WK(NFTEMP),NXX,
     1   JJER)
C                                  COMPUTE THE SUM FOR THE INTERVAL
         DO 125 J=NWTXX,NN
            SUM = SUM+WK(J)*WK(NFTMP1+J)
  125    CONTINUE
C                                  COMPUTE THE WEIGHT FOR YL(I+1)
         WTYY = FRDTRD*(WK(NYLS+I)-Y(II))
C                                  ACCUMULATE THE INTEGRAL
         Q = Q+WTYY*SUM
         II = II+1
  130 CONTINUE
C                                  IF THE LIMITS OF INTEGRATION ARE
C                                    REVERSED, CHANGE THE SIGN OF THE
C                                    INTEGRAL
  135 IF (B.LT.A) Q = -Q
      IF (D.LT.C) Q = -Q
      CALL UERSET (LEVOLD,LEVEL)
  140 IER = MAX0(JER,KER,LER,MER,NER)
 9000 CONTINUE
      IF (NER.NE.0) CALL UERTST (NER,'DBCQDU')
      IF (NER.NE.0) RETURN
      IF (JER.NE.0) CALL UERTST (JER,'DBCQDU')
      IF (KER.NE.0) CALL UERTST (KER,'DBCQDU')
      IF (LER.NE.0) CALL UERTST (LER,'DBCQDU')
      IF (MER.NE.0) CALL UERTST (MER,'DBCQDU')
 9005 RETURN
      END
