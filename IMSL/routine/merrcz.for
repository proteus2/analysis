C   IMSL ROUTINE NAME   - MERRCZ
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1982
C
C   PURPOSE             - EVALUATE A FUNCTION RELATED TO THE
C                           COMPLEMENTED ERROR FUNCTION FOR A
C                           COMPLEX ARGUMENT
C
C   USAGE               - CALL MERRCZ (Z,W,IER)
C
C   ARGUMENTS    Z      - INPUT COMPLEX ARGUMENT.  LET X AND Y
C                           REPRESENT THE REAL AND IMAGINARY PARTS OF Z
C                           RESPECTIVELY.  THEN X**2 + Y**2 MUST BE
C                           LESS THAN MACHINE INFINITY.  ALSO, IF Z IS
C                           LOCATED IN THE SECOND, THIRD OR FOURTH
C                           QUADRANT OF THE COMPLEX PLANE, OR ON THE
C                           BOUNDARY OF THE QUADRANTS, THEN Y**2 - X**2
C                           MUST BE LESS THAN OR EQUAL TO THE LARGEST
C                           ACCEPTABLE ARGUMENT FOR THE FORTRAN
C                           EXPONENTIAL FUNCTION.
C                W      - OUTPUT COMPLEX VALUE
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER = 129 INDICATES THAT THE ABSOLUTE VALUE
C                             OF INPUT ARGUMENT Z IS GREATER THAN THE
C                             SQUARE ROOT OF MACHINE INFINITY.  W IS SET
C                             TO 0.
C                           IER = 130 INDICATES THAT Z = (X,Y) IS NOT IN
C                             THE FIRST QUADRANT AND Y**2 - X**2 IS
C                             GREATER THAN THE LARGEST ACCEPTABLE
C                             ARGUMENT FOR THE FORTRAN EXPONENTIAL
C                             FUNCTION. W IS SET TO MACHINE INFINITY.
C
C   REQD. IMSL ROUTINES - UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS       ON THOSE MACHINES WHICH DO NOT OFFER THE APPROPRIATE
C                 COMPLEX DATA TYPE, INPUT ARGUMENT Z AND OUTPUT
C                 ARGUMENT W ARE TREATED AS VECTORS OF LENGTH 2.  THE
C                 USER SHOULD DIMENSION THEM AS SUCH ON THOSE MACHINES.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE MERRCZ (Z,W,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            IER
      COMPLEX*16         Z,W
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            IFLAG,N,NCAPN,NM1,NN,NP1,NU,NUP1,NUP2
      DOUBLE PRECISION   C,C1,H,H2,R1,R2,RE,RIMAG,RLMBDA,S,S1,S2,SRINF,
     *                   T1,T2,X,XINF,XLARGE,XSMALL,XX,Y,YY,XSMARG,
     *                   SRSMAG
      DOUBLE PRECISION   DIMAG,DREAL
      COMPLEX*16         ZDUM
      COMPLEX*16         Q,WNU,WW,WZ,ZNU,ZSQ
      LOGICAL            B
      DATA               XSMALL/ -180.2182D0/
      DATA               XLARGE/ 174.0D0/
      DATA               XINF/.723700557733226D+76/
      DATA               SRINF/ .8507059173023461D+38/
      DATA               XSMARG/Z0010000000000000/
      DREAL(ZDUM) = ZDUM
      DIMAG(ZDUM) = (0.D0,-1.D0)*ZDUM
C                                  FIRST EXECUTABLE STATEMENT
      C1 = 1.12837916709551D0
      B = .FALSE.
      X = DREAL(Z)
      Y = DIMAG(Z)
      IF(DABS(X).GT.SRINF .OR. DABS(Y).GT.SRINF) GO TO 5
      SRSMAG = DSQRT(XSMARG) * 5.5D0
      IF(DABS(X).LT.SRSMAG) X = 0.D0
      IF(DABS(Y).LT.SRSMAG) Y = 0.D0
      XX = X*X
      YY = XINF-(Y*Y)
      IF(XX .LE. YY) GO TO 10
    5 IER = 129
      W = DCMPLX(0.D0,0.D0)
      GO TO 9000
   10 XX = Y*Y-X*X
      IF(XX.LE.XLARGE) GO TO 15
      IER = 130
      W = DCMPLX(XINF,XINF)
      GO TO 9000
   15 XX = X
      YY = Y
      IFLAG = 1
      IF(X.GE.0.D0) GO TO 20
      XX = -X
      IFLAG = 2
C                                  X IS NEGATIVE
   20 IF(Y.GE.0.D0) GO TO 25
      YY = -Y
      IF(IFLAG.EQ.1) IFLAG = 3
      IFLAG = IFLAG+1
   25 IF(YY.GE.4.29D0 .OR. XX.GE.5.33D0) GO TO 30
      S = (1.D0-YY/4.29D0)*DSQRT(1.D0-XX*XX/28.41D0)
      H = 1.6D0*S
      H2 = H+H
      NCAPN = 6.D0+23.D0*S
      RLMBDA = H2**NCAPN
      NU = 9.D0+21.D0*S
      GO TO 35
   30 H = 0.D0
      NCAPN = 0
      NU = 8
   35 IF(H.EQ.0.D0.OR.RLMBDA.EQ.0.D0) B = .TRUE.
      R1 = 0.D0
      R2 = 0.D0
      S1 = 0.D0
      S2 = 0.D0
      NUP1 = NU+1
      NUP2 = NU+2
      DO 40  NN=1,NUP1
         N = NUP2-NN
         NM1 = N-1
         T1 = YY+H+N*R1
         T2 = XX-N*R2
         C = .5D0/(T1*T1+T2*T2)
         R1 = C*T1
         R2 = C*T2
         IF(H.LE.0.D0.OR.NM1.GT.NCAPN) GO TO 40
         T1 = RLMBDA+S1
         S1 = R1*T1-R2*S2
         S2 = R2*T1+R1*S2
         RLMBDA = RLMBDA/H2
   40 CONTINUE
      IF(B) GO TO 45
      RE = C1*S1
      RIMAG = C1*S2
      GO TO 50
   45 RE = C1*R1
      RIMAG = C1*R2
   50 IF(YY.EQ.0.D0) RE = 0.D0
      IF(YY.EQ.0.D0 .AND. -XX*XX.GT.XSMALL) RE = DEXP(-XX*XX)
      IF(IFLAG.EQ.1) GO TO 60
      ZNU = DCMPLX(XX,YY)
      WZ = DCMPLX(RE,RIMAG)
      ZSQ = ZNU*ZNU
      WNU = CDEXP(-ZSQ)
      WW = WNU+WNU
      Q = WW-WZ
      IF(IFLAG.EQ.3) GO TO 55
      RE = DREAL(Q)
      RIMAG = -DIMAG(Q)
      IF(IFLAG.EQ.4) GO TO 60
      ZNU = DCMPLX(XX,-YY)
      ZSQ = ZNU*ZNU
      WNU = CDEXP(-ZSQ)
      WW = WNU+WNU
      Q = WW-DCMPLX(RE,RIMAG)
   55 RE = DREAL(Q)
      RIMAG = DIMAG(Q)
   60 W = DCMPLX(RE,RIMAG)
      GO TO 9005
 9000 CALL UERTST(IER,6HMERRCZ)
 9005 RETURN
      END
