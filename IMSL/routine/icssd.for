C   IMSL ROUTINE NAME   - ICSSD
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1980
C
C   PURPOSE             - NUCLEUS CALLED ONLY BY IMSL SUBROUTINE ICSSCV
C
C   REQD. IMSL ROUTINES - EQRT2S,ICSSE,LEQT1B,LSVDB,LSVDF,LSVG1,
C                           LSVG2,VHS12,UERSET,UERTST,UGETIO,ZXGSP
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE ICSSD  (N,NMAX,X,Y,A,B,C,D,RES,H,Z,TS,U,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N,NMAX,IER
      REAL               X(1),Y(1),A(1),B(1),C(1),D(1),RES(1),H(1),
     1                   Z(NMAX,1),TS(NMAX,1),U(NMAX,1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            IL1,IP1,IDUM(1),I,JER,J,KTER,K,LEVOLD,LEV,
     1                   MAXIT,NL1,NL2
      REAL               DD1,DD2,F1,F2,F3,HSUM,PA,PB,PG,PT,P,SGN,TOL,VA,
     1                   VB,VPA,VPB,ICSSE,THIRD,TEMP
      EXTERNAL           ICSSE
      COMMON /VCOM/      VPB
      DATA               THIRD/.3333333/
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
      TOL = 0.01
      MAXIT = 30
      PG = 0.0
      NL1 = N-1
      NL2 = N-2
      DO 5 I=1,NL1
         H(I) = X(I+1)-X(I)
         IF (H(I).GT.0.0) GO TO 5
C                                  DATA POINTS NOT ORDERED CORRECTLY
         IER = 131
         GO TO 150
    5 CONTINUE
      DO 10 I=1,NL2
         J = I+1
         A(I) = 2.0*(H(I)+H(J))*THIRD
         B(J) = H(J)*THIRD
   10 CONTINUE
      DO 20 I=1,NL2
         DO 15 J=I,NL2
            Z(I,J) = 0.0
            Z(J,I) = 0.0
   15    CONTINUE
         Z(I,I) = 1.0
   20 CONTINUE
C                                  FIND EIGENVALUES AND EIGENVECTORS
C                                    OF THE TRIDIAGONAL MATRIX (A,B)
      CALL EQRT2S (A,B,NL2,Z,NMAX,JER)
      IF (JER.LE.128) GO TO 25
      IER = 132
      GO TO 150
   25 CONTINUE
      DO 30 I=1,NL2
         A(I) = SQRT(A(I))
   30 CONTINUE
      DO 40 I=1,NL2
      DO 40 J=I,NL2
         TEMP = 0.0
         DO 35 K=1,NL2
            TEMP = TEMP+Z(I,K)*Z(J,K)/A(K)
   35    CONTINUE
         TS(I,J) = TEMP
         TS(J,I) = TEMP
   40 CONTINUE
      DO 80 I=1,N
         IL1 = I-1
         F1 = 0.0
         F2 = 0.0
         F3 = 0.0
         IF (I.GE.3) F1 = 1.0/H(IL1)
         IF (I.GE.2.AND.I.LE.NL1) F2 = -(1.0/H(I)+1.0/H(IL1))
         IF (I.LE.NL2) F3 = 1.0/H(I)
         IF (I.GE.3) GO TO 55
         DO 45 J=1,NL2
            Z(I,J)=F3*TS(I,J)
   45    CONTINUE
         IF (I.EQ.1) GO TO 80
         DO 50 J=1,NL2
            Z(I,J)=Z(I,J)+F2*TS(IL1,J)
   50    CONTINUE
         GO TO 80
   55    CONTINUE
         IF (I.LT.NL1) GO TO 70
         DO 60 J=1,NL2
            Z(I,J)=F1*TS(I-2,J)
   60    CONTINUE
         IF (I.EQ.N) GO TO 80
         DO 65 J=1,NL2
            Z(I,J)=Z(I,J)+F2*TS(IL1,J)
   65    CONTINUE
         GO TO 80
   70    CONTINUE
         DO 75 J=1,NL2
            Z(I,J) = F1*TS(I-2,J)+F2*TS(IL1,J)+F3*TS(I,J)
   75    CONTINUE
   80 CONTINUE
      DO 90 I=1,N
         DO 85 J=1,N
            U(I,J) = 0.0
   85    CONTINUE
         U(I,I) = 1.0
   90 CONTINUE
C                                  SINGULAR VALUE DECOMPOSITION OF Z
      CALL LSVDF (Z,NMAX,N,NL2,U,NMAX,N,D,TS,JER)
      IF (JER.LT.128) GO TO 95
      IER = 132
      GO TO 150
   95 DO 105 J=1,NL2
         TEMP = 0.0
         DO 100 K=1,N
            TEMP = TEMP+U(J,K)*Y(K)
  100    CONTINUE
         C(J) = TEMP
         D(J) = D(J)*D(J)
  105 CONTINUE
      PA = PG
      VA = ICSSE(PA,D,C,IDUM,N,NMAX)
      VPA = VPB
      IF (VPA.EQ.0.0) GO TO 125
      SGN = 3.0*SIGN(1.0,VPA)
      DO 110 KTER=1,MAXIT
         PB = PA-SGN
         VB = ICSSE(PB,D,C,IDUM,N,NMAX)
         IF (VB.GE.VA) GO TO 115
         IF (VPA*VPB.LE.0.0) GO TO 115
         PA = PB
         VA = VB
         VPA = VPB
  110 CONTINUE
C                                  ATTEMPT TO BRACKET MINIMUM OF
C                                    CROSS VALIDATION FUNCTION FAILS
      IER = 132
      GO TO 150
  115 CONTINUE
      IF (PA.LT.PB) GO TO 120
      PT = PA
      PA = PB
      PB = PT
  120 CONTINUE
      LEV = 0
      CALL UERSET (LEV,LEVOLD)
C                                  FIND MINIMUM OF CROSS VALIDATION
C                                    FUNCTION
      CALL ZXGSP (ICSSE,D,C,IDUM,N,NMAX,PA,PB,TOL,PG,JER)
      CALL UERSET (LEVOLD,LEV)
      IF (JER.LT.128) GO TO 125
      IF (JER.EQ.131) GO TO 125
      IER = 132
      GO TO 150
  125 CONTINUE
      P = EXP(PG)
      DO 135 I=1,N
         TEMP = Y(I)
         DO 130 J=1,NL2
            TEMP = TEMP-D(J)/(P+D(J))*C(J)*U(J,I)
  130    CONTINUE
         A(I) = TEMP
         RES(I) = Y(I)-TEMP
  135 CONTINUE
C                                  CALCULATE SPLINE COEFFICIENTS
      C(1) = 0.0
      C(N) = 0.0
      DO 140 I=1,NL2
         IP1 = I+1
         TEMP = H(IP1)
         Z(I,2) = 2.0
         HSUM = 1.0/(TEMP+H(I))
         Z(I,3) = TEMP*HSUM
         Z(I,1) = 1.0-Z(I,3)
         DD2 = (A(I+2)-A(IP1))/TEMP
         DD1 = (A(IP1)-A(I))/H(I)
         C(IP1) = 3.0*(DD2-DD1)*HSUM
  140 CONTINUE
      CALL LEQT1B (Z,NL2,1,1,NMAX,C(2),1,NL2,0,TS,JER)
      DO 145 I=1,NL1
         IP1 = I+1
         TEMP = H(I)
         D(I) = (C(IP1)-C(I))*THIRD/TEMP
         B(I) = (A(IP1)-A(I))/TEMP-TEMP*(C(I)+D(I)*TEMP)
  145 CONTINUE
      B(N) = B(NL1)+H(NL1)*(2.0*C(NL1)+3.0*D(NL1)*H(NL1))
      D(N) = D(NL1)
  150 RETURN
      END
