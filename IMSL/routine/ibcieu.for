C   IMSL ROUTINE NAME   - IBCIEU
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1980
C
C   PURPOSE             - BICUBIC SPLINE TWO-DIMENSIONAL INTERPOLATOR
C
C   USAGE               - CALL IBCIEU (F,IFD,X,NX,Y,NY,XL,NXL,YL,NYL,
C                           FL,IFLD,WK,IER)
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
C                         NOTE - THE COORDINATE PAIRS (X(I),Y(J)), FOR
C                           I=1,...,NX AND J=1,...,NY, GIVE THE POINTS
C                           WHERE THE FUNCTION VALUES F(I,J) ARE
C                           DEFINED.
C                XL     - VECTOR OF LENGTH NXL. (INPUT)
C                NXL    - NUMBER OF ELEMENTS IN XL. (INPUT)
C                YL     - VECTOR OF LENGTH NYL. (INPUT)
C                NYL    - NUMBER OF ELEMENTS IN YL. (INPUT)
C                         NOTE - THE COORDINATE PAIRS (XL(I),YL(J)),
C                           FOR I=1,...,NXL AND J=1,...,NYL, GIVE THE
C                           POINTS AT WHICH THE INTERPOLATORY BICUBIC
C                           SPLINE IS TO BE EVALUATED.
C                FL     - NXL BY NYL MATRIX CONTAINING THE INTERPOLATORY
C                           BICUBIC SPLINE VALUES. (OUTPUT) FL(I,J) IS
C                           SET TO THE VALUE OF THE INTERPOLATORY
C                           BICUBIC SPLINE AT (XL(I),YL(J)) FOR
C                           I=1,...,NXL AND J=1,...,NYL. NOTE THAT THE
C                           NUMBER OF COLUMNS IN FL MUST BE .GE.
C                           MAX(NYL,NY) SINCE FL IS ALSO USED AS
C                           WORKING STORAGE (OF SIZE NXL BY NY) DURING
C                           THE COMPUTATION.
C                IFLD   - ROW DIMENSION OF THE MATRIX FL EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT
C                           IN THE CALLING PROGRAM. (INPUT)
C                WK     - WORK VECTOR OF LENGTH
C                           MAX((NX-1)*3,(NY-1)*3+NY).
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER = 129, IFD IS LESS THAN NX.
C                           IER = 130, IFLD IS LESS THAN NXL.
C                           IER = 131, NX IS LESS THAN 2.
C                           IER = 132, NY IS LESS THAN 2.
C                           IER = 133, VECTOR X IS NOT ORDERED SO THAT
C                             X(1) .LT. X(2) ... .LT. X(NX).
C                           IER = 134, VECTOR Y IS NOT ORDERED SO THAT
C                             Y(1) .LT. Y(2) ... .LT. Y(NY).
C                         WARNING ERROR
C                           IER = 37, XL(I) IS LESS THAN X(1) OR
C                             GREATER THAN X(NX).
C                           IER = 38, YL(I) IS LESS THAN Y(1) OR
C                             GREATER THAN Y(NY).
C
C
C   REQD. IMSL ROUTINES - ICSEVU,ICSCCU,UERSET,UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE IBCIEU (F,IFD,X,NX,Y,NY,XL,NXL,YL,NYL,FL,IFLD,WK,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            IFD,NX,NY,NXL,NYL,IFLD,IER
      REAL               F(IFD,NY),X(NX),Y(NY),XL(NXL),YL(NYL),
     1                   FL(IFLD,1),WK(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            KER,LER,MER,NXM1,NYM1,IY,JER,KYL,KYLP1,IXL,
     1                   IYL,LEVEL,LEVOLD
C                                  FIRST EXECUTABLE STATEMENT
      KER = 0
      LER = 0
C                                  CHECK IFD .GE. NX
      IF (IFD .LT. NX) GO TO 30
C                                  CHECK IFLD .GE. NX.
      IF (IFLD .LT. NXL) GO TO 35
C                                  CHECK NX .GE. 2
      IF (NX .LT. 2) GO TO 36
C                                  CHECK NY .GE. 2
      IF (NY .LT. 2) GO TO 37
C                                  CALL UERSET TO SILENCE WARNING
C                                    MESSAGES FROM ICSEVU
      LEVEL = 2
      CALL UERSET (LEVEL,LEVOLD)
      MER = 0
      NXM1 = NX-1
      NYM1 = NY-1
C                                  INTERPOLATE IN THE X-DIRECTION
      DO 10 IY=1,NY
C                                  CALCULATE THE COEFFICIENTS
         CALL ICSCCU (X,F(1,IY),NX,WK(1),NXM1,JER)
C                                  CHECK FOR ERROR IN ICSCCU
         IF (JER .NE. 0) GO TO 40
C                                  EVALUATE
         CALL ICSEVU (X,F(1,IY),NX,WK(1),NXM1,XL,FL(1,IY),NXL,JER)
C                                  CHECK FOR ERROR IN ICSEVU
         IF (JER .NE. 0) KER = 37
   10 CONTINUE
      KYL = NYM1*3
      KYLP1 = KYL+1
C                                  INTERPOLATE IN THE Y-DIRECTION
      DO 25 IXL=1,NXL
         DO 15 IY=1,NY
            WK(KYL+IY) = FL(IXL,IY)
   15    CONTINUE
C                                  CALCULATE THE COEFFICIENTS
         CALL ICSCCU (Y,WK(KYLP1),NY,WK(1),NYM1,JER)
C                                  CHECK FOR ERROR IN ICSCCU
         IF (JER .NE. 0) GO TO 45
C                                  EVALUATE
         DO 20 IYL=1,NYL
            CALL ICSEVU (Y,WK(KYLP1),NY,WK(1),NYM1,YL(IYL),FL(IXL,IYL),
     1                   1,JER)
C                                  CHECK FOR ERROR IN ICSEVU
            IF (JER .NE. 0) LER = 38
   20    CONTINUE
   25 CONTINUE
      GO TO 46
C                                  HANDLE ERRORS
   30 MER = 129
      GO TO 50
   35 MER = 130
      GO TO 50
   36 MER = 131
      GO TO 50
   37 MER = 132
      GO TO 50
   40 MER = 133
      GO TO 46
   45 MER = 134
   46 CALL UERSET (LEVOLD,LEVEL)
   50 IER = MAX0(MER,KER,LER)
 9000 CONTINUE
      IF (MER .NE. 0) CALL UERTST(MER,'IBCIEU')
      IF (KER .NE. 0) CALL UERTST(KER,'IBCIEU')
      IF (LER .NE. 0) CALL UERTST(LER,'IBCIEU')
 9005 RETURN
      END
 
R; T=0.05/0.50 21:06:02
 
