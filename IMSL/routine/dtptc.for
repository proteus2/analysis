C   IMSL ROUTINE NAME   - DTPTC
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1982
C
C   PURPOSE             - NUCLEUS CALLED ONLY BY IMSL SUBROUTINE DTPTB
C
C   REQD. IMSL ROUTINES - SINGLE/DTPTD,DTPTE,LEQT2F,LUDATN,LUELMN,
C                           LUREFN,UERSET,UERTST,UGETIO
C                       - DOUBLE/DTPTD,DTPTE,LEQT2F,LUDATN,LUELMN,
C                           LUREFN,UERSET,UERTST,UGETIO,VXADD,VXMUL,
C                           VXSTO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE DTPTC  (FCNI,FCNJ,FCNB,XA,XB,X,MAXR,Y,IY,N,NITER,STEMP,
     1                   WORK,TOL,EPS,A,F,ATEMP,BTEMP,C,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            MAXR,IY,N,NITER,IER
      REAL               XA,XB,X(1),Y(IY,1),STEMP(1),WORK(1),
     1                   TOL,EPS,A(N,N,1),F(N,1),ATEMP(N,1),BTEMP(N,1),
     2                   C(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            IDGT,II,INDA,INDB,IND,ITER,ITRY,I,JER,J,K,L,
     1                   MAX,MM1,M,N2P1,N2,NEWTON,NM,NP1,IP1
      INTEGER            NCOM
      REAL               ANR,ASUM,BNR,BSUM,EBIG,ECOND,ERRMX,HAK,HBK,H,
     1                   TOLM,XP,XSTRT,TEMP,TEMP1,EECOND,RECOND
      REAL               ALPHA,ASGN
      EXTERNAL           FCNI,FCNJ,DTPTE
      COMMON /DALPHA/    ALPHA
      COMMON /DTPTZ /    ASGN,NCOM
      DATA               EECOND/Z3C100000/
      DATA               EBIG/0.8507057E38/
C                                  FIRST EXECUTABLE STATEMENT
      ECOND = 1.0/SQRT(EECOND)
C                                  ECOND=1/SQRT(MACHINE TOLERANCE)
C                                  EBIG=SQRT(LARGEST REAL NUMBER)
      RECOND = 1.0/ECOND
      MAX = MAXR
      NCOM = N
      NP1 = N+1
      ITRY = 0
      M = IABS(MAX)
    5 ALPHA = 1.0
      H = (XB-XA)*0.01
      N2P1 = N*NP1
      C(1) = 4.
      C(2) = 0.
      C(3) = (XB-XA)*1.0E-6
      C(4) = H
      C(5) = 0.
      C(6) = (XB-XA)*0.1
      C(7) = 0.
      C(8) = 0.
      C(9) = 0.
      N2 = N*N
      DO 10 J=1,N2
         C(J+30) = EBIG
   10 CONTINUE
      DO 15 J=1,N
         L=N2+J
         C(L+30) = 1.0
   15 CONTINUE
      IF (MAX.LT.0) GO TO 25
      DO 20 J=1,N
      DO 20 I=1,M
         Y(J,I) = 0.0
   20 CONTINUE
   25 CONTINUE
      IF (MAX.LE.0.OR.ITRY.NE.0) GO TO 70
C                                  PROGRAM CALCULATES SHOOTING PTS
      M = 1
      X(1) = XA
      L = 0
   30 DO 40 J=1,N
         DO 35 K=1,N
            ATEMP(J,K) = 0.0
            BTEMP(J,K) = 0.0
   35    CONTINUE
         ATEMP(J,NP1) = Y(J,M)
         BTEMP(J,NP1) = Y(J,M)
         ATEMP(J,J) = 1.0
         BTEMP(J,J) = 1.0
   40 CONTINUE
   45 TEMP = L
      XP = XA+TEMP*H
      INDA = 2
      ASGN = 1.0
      CALL DTPTD (FCNI,FCNJ,STEMP,N2P1,DTPTE,XP,ATEMP,XP+H,TOL,INDA,
     1 C,N2P1,WORK,JER)
      XP = XA+TEMP*H
      INDB = 2
      ASGN = -1.0
      CALL DTPTD (FCNI,FCNJ,STEMP,N2P1,DTPTE,XP,BTEMP,XP+H,TOL,INDB,
     1 C,N2P1,WORK,JER)
      IF (INDA.EQ.3.AND.INDB.EQ.3) GO TO 50
      IER = 129
      RETURN
C                                  MONITOR CONDITION OF JACOBIAN.
C                                    INSERT SHOOTING POINTS WHEN IT
C                                    BECOMES LARGE
   50 ANR = 0.0
      BNR = 0.0
      DO 60 I=1,N
         ASUM = 0.0
         BSUM = 0.0
         DO 55 J=1,N
            ASUM = ASUM+ABS(ATEMP(I,J))
            BSUM = BSUM+ABS(BTEMP(I,J))
   55    CONTINUE
         IF (ASUM.GT.ANR) ANR = ASUM
         IF (BSUM.GT.BNR) BNR = BSUM
   60 CONTINUE
      L = L+1
      IF (ANR*BNR.LT.ECOND.AND.L.LT.100) GO TO 45
      M = M+1
      IF (M.LE.MAX) GO TO 65
      TEMP = MAX
      MAXR = TEMP*(XB-XA)/(X(M-1)-XA)+5
      IER = 130
      RETURN
   65 X(M) = XP
      IF (L.LT.100) GO TO 30
      MAXR = M
C                                  BEGIN NEWTONS ITERATION
   70 NM = N*M
      MM1 = M-1
      IF (NITER.GE.10) TEMP1 = 1.0/AMAX1(FLOAT(NITER-9),1.0)
      DO 170 ITER=1,NITER
         TEMP = ITER-1
         IF (NITER.GE.10) ALPHA = AMIN1(TEMP*TEMP1,1.0)
C                                  FROM XA TO XB SOLVE N**2+N
C                                    SIMULTANEOUS INITIAL VALUE ODES
C                                    TO CALCULATE F VECTOR AND JACOBIAN
C                                    FOR NEWTON METHOD.
         DO 95 I=1,MM1
            IP1 = I+1
            DO 80 J=1,N
               DO 75 K=1,N
                  ATEMP(J,K) = 0.0
   75          CONTINUE
               ATEMP(J,J) = 1.0
               ATEMP(J,NP1) = Y(J,I)
   80       CONTINUE
            XSTRT = X(I)
            ASGN = -1.0
            TOLM = TOL*0.5
            IND = 2
            CALL DTPTD (FCNI,FCNJ,STEMP,N2P1,DTPTE,XSTRT,ATEMP,X(IP1),
     1      TOLM,IND,C,N2P1,WORK,JER)
            IF (IND.EQ.3) GO TO 85
            IER = 129
            GO TO 200
   85       DO 90 J=1,N
               F(J,I) = ATEMP(J,NP1)-Y(J,IP1)
               DO 90 K=1,N
                  A(J,K,I) = ATEMP(J,K)
   90       CONTINUE
   95    CONTINUE
         DO 110 K=1,N
            HAK = AMAX1(ABS(Y(K,1)),0.1)*RECOND
            Y(K,1) = Y(K,1)+HAK
            CALL FCNB (N,Y(1,1),Y(1,M),ATEMP(1,NP1))
            Y(K,1) = Y(K,1)-2.*HAK
            CALL FCNB (N,Y(1,1),Y(1,M),BTEMP(1,NP1))
            Y(K,1) = Y(K,1)+HAK
            DO 100 J=1,N
               ATEMP(J,K) = (ATEMP(J,NP1)-BTEMP(J,NP1))/(2.0*HAK)
  100       CONTINUE
            HBK = AMAX1(ABS(Y(K,M)),0.1)*RECOND
            Y(K,M) = Y(K,M)+HBK
            CALL FCNB (N,Y(1,1),Y(1,M),ATEMP(1,NP1))
            Y(K,M) = Y(K,M)-2.0*HBK
            CALL FCNB (N,Y(1,1),Y(1,M),BTEMP(1,NP1))
            Y(K,M) = Y(K,M)+HBK
            DO 105 J=1,N
               BTEMP(J,K) = (ATEMP(J,NP1)-BTEMP(J,NP1))/(2.0*HBK)
  105       CONTINUE
  110    CONTINUE
         CALL FCNB (N,Y(1,1),Y(1,M),F(1,M))
         ERRMX = 0.0
         DO 115 I=1,M
         DO 115 J=1,N
            IF (ABS(F(J,I)).GT.ERRMX) ERRMX = ABS(F(J,I))
  115    CONTINUE
C                                  SOLVE LINEAR SYSTEM TO DO ONE
C                                    ITERATION OF NEWTONS METHOD.
         DO 135 I=1,MM1
            DO 125 J=1,N
            DO 125 K=1,N
               A(J,K,M) = 0.0
               DO 120 L=1,N
                  A(J,K,M) = A(J,K,M)+ATEMP(J,L)*A(L,K,I)
  120          CONTINUE
  125       CONTINUE
            DO 130 J=1,N
            DO 130 K=1,N
               ATEMP(J,K) = A(J,K,M)
               F(J,M) = F(J,M)-ATEMP(J,K)*F(K,I)
  130       CONTINUE
  135    CONTINUE
         DO 140 J=1,N
         DO 140 K=1,N
            ATEMP(J,K) = ATEMP(J,K)+BTEMP(J,K)
  140    CONTINUE
         IDGT = 0
         CALL LEQT2F (ATEMP,1,N,N,F(1,M),IDGT,WORK,JER)
         IF (JER.NE.129.AND.JER.NE.131) GO TO 145
         IER = 132
         GO TO 200
  145    DO 160 II=1,MM1
            I = M-II
            DO 150 J=1,N
               ATEMP(J,NP1) = F(J,I)+F(J,I+1)
  150       CONTINUE
            DO 155 J=1,N
               F(J,I) = 0.0
               DO 155 K=1,N
                  F(J,I) = F(J,I)+A(J,K,I)*ATEMP(K,NP1)
  155       CONTINUE
  160    CONTINUE
         DO 165 I=1,M
         DO 165 J=1,N
            Y(J,I) = Y(J,I)-F(J,I)
  165    CONTINUE
         IF (ALPHA.LT.1.0) GO TO 170
         IF (ERRMX.GT.0.5*EPS) GO TO 170
C                                  NEWTONS METHOD HAS CONVERGED.
         NEWTON = 1
         GO TO 175
  170 CONTINUE
C                                  NEWTONS METHOD FAILS TO CONVERGE
      NEWTON = 0
C                                  CHECK IF SOLUTION TO INITIAL VALUE
C                                    PROBLEM WITH INITIAL VALUES
C                                    Y(J,1),J=1-N SATISFIES BOUNDARY
C                                    CONDITIONS.
  175 DO 180 I=1,N
         F(I,1) = Y(I,1)
  180 CONTINUE
      ALPHA = 1.0
      XSTRT = XA
      IND = 2
      DO 185 J=1,N
         C(J+30) = 1.0
  185 CONTINUE
      CALL DTPTD (FCNI,FCNJ,STEMP,N,DTPTE,XSTRT,F(1,1),XB,TOL,IND,C,
     1 N,WORK,JER)
      IF (IND.NE.3) GO TO 195
      CALL FCNB (N,Y(1,1),F(1,1),F(1,M))
      DO 190 I=1,N
         IF (ABS(F(I,M)).GT.EPS) GO TO 195
  190 CONTINUE
      IER = 0
      RETURN
  195 IER = 33
      IF(NEWTON.EQ.1) RETURN
      IER = 131
  200 IF(MAX.LT.0.OR.M.EQ.MAX) RETURN
      M = MAX
      MAXR = M
      ITRY=1
      TEMP = 1.0/(MAX-1)
      DO 205 I=1,MAX
         TEMP1 = I-1
         X(I) = XA + TEMP1*TEMP*(XB-XA)
  205 CONTINUE
      GO TO 5
      END
