C   IMSL ROUTINE NAME   - MMDEN
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - NOVEMBER 1, 1984
C
C   PURPOSE             - EXPONENTIAL INTEGRALS OF INTEGER ORDER
C                           FOR REAL ARGUMENT X SCALED BY EXP(X)
C
C   USAGE               - CALL MMDEN (X,N,F,IER)
C
C   ARGUMENTS    X      - INPUT ARGUMENT. X MUST BE GREATER THAN ZERO. X
C                           MUST BE TYPED APPROPRIATELY IN THE CALLING
C                           PROGRAM. (SEE PRECISION/HARDWARE SECTION.)
C                N      - INPUT POSITIVE INTEGER SPECIFYING THE MAXIMUM
C                           ORDER FOR WHICH THE EXPONENTIAL INTEGRAL IS
C                           TO BE CALCULATED.
C                F      - OUTPUT VECTOR OF LENGTH N CONTAINING THE
C                           COMPUTED EXPONENTIAL INTEGRALS SCALED BY
C                           EXP(X). F MUST BE TYPED APPROPRIATELY IN
C                           THE CALLING PROGRAM.
C                           (SEE PRECISION/HARDWARE SECTION.)
C                           F(1) WILL CONTAIN THE COMPUTED VALUE FOR
C                           ORDER 1, F(2) FOR ORDER 2, ETC.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER = 129 INDICATES THAT EITHER X OR N IS
C                             OUT OF RANGE. F(I), (I=1,N) IS SET TO
C                             MACHINE INFINITY.
C                           IER = 130 INDICATES THAT AN ERROR OCCURRED
C                             IN MMDEI. F(I), (I=1,N) IS SET TO
C                             MACHINE INFINITY.
C
C   REQD. IMSL ROUTINES - MMDEI,UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE MMDEN(X,N,F,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N,IER
      DOUBLE PRECISION   X,F(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,IM1,IOPT,IP1,JER,N1,K,K1,M,NM1
      DOUBLE PRECISION   EPS,UE,VE,WE,WE1,U0,V0,W0,W01,WSEPS,W,R,S,
     *                   MMDEI,RI,RN1,RK,WOMWE,RM1,XINF
      DATA               EPS/Z3410000000000000/
      DATA               XINF/.723700557733226D+76/
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
      IF (X.GT.0.0D0 .AND. N.GT.0) GO TO 10
      IER = 129
      F(1) = XINF
      IF (N .LT. 2) GO TO 9000
      DO 5 I = 2, N
        F(I) = XINF
    5 CONTINUE
      GO TO 9000
   10 IF (X .LE. 2.0D0/EPS) GO TO 20
      DO 15 I = 1, N
        F(I) = 1.0D0/X
   15 CONTINUE
      GO TO 9005
   20 IF (X .GT. 1.5D0) GO TO 40
      IOPT = 2
      F(1) = MMDEI(IOPT,X,JER)
      F(1) = F(1)*DEXP(X)
      IF (JER.EQ.0) GO TO 30
      DO 25 I = 1, N
        F(I) = XINF
   25 CONTINUE
      IER = 130
      GO TO 9000
   30 CONTINUE
C                                  IF MAXIMUM ORDER IS 1, NO FURTHER
C                                    VALUES ARE TO BE FIGURED. RETURN.
      IF (N .EQ. 1) GO TO 9005
      NM1 = N - 1
      IP1 = 2
      DO 35 I = 1, NM1
        RI = I
        F(IP1) = (1.0D0-X*F(I))/RI
        IP1 = IP1 + 1
   35 CONTINUE
      GO TO 9005
   40 RN1 = X + 0.5
      N1 = N
      IF (RN1 .LE. DBLE(FLOAT(N))) N1 = X + 0.5
      RN1 = N1
      UE = 1.0D0
      VE = 1.0D0/(X+RN1)
      WE = VE
      WE1 = 0.0D0
      U0 = 1.0D0
      V0 = -RN1/(X*(X+RN1+1.0D0))
      W01 = 1.0D0/X
      W0 = V0 + W01
      W = (WE+W0)*.5D0
      K1 = 1
C                                  BEGIN ITERATION
   45 CONTINUE
      RK = K1
      WE1 = WE
      W01 = W0
      R = RN1 + RK
      S = R + X + RK
      UE = 1.0D0/(1.0D0-RK*(R-1.0D0)*UE/((S-2.0D0)*S))
      U0 = 1.0D0/(1.0D0-RK*R*U0/(S*S-1.0D0))
      VE = VE*(UE-1)
      V0 = V0*(U0-1.0D0)
      WE = WE + VE
      W0 = W0 + V0
      W = (WE+W0)*.5D0
      K1 = K1 + 1
      WOMWE = W0 - WE
      WSEPS = EPS*W
      IF (WOMWE.GT.WSEPS .AND. WE.GT.WE1 .AND. W0.LT.W01) GO TO 45
      F(N1) = W
      M = N1
      DO 50 I = 1, M
        IM1 = M - I
        IF (IM1 .EQ. 0) GO TO 55
        RM1 = IM1
        W = (1.0D0-RM1*W)/X
        F(IM1) = W
   50 CONTINUE
   55 IF (N1.EQ.N) GO TO 9005
      NM1 = N - 1
      DO 60 I = N1, NM1
        RI = I
        IP1 = I + 1
        F(IP1) = (1.0D0-X*F(I))/RI
   60 CONTINUE
      GO TO 9005
 9000 CONTINUE
      CALL UERTST(IER,6HMMDEN )
 9005 RETURN
      END
