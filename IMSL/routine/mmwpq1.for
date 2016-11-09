C   IMSL ROUTINE NAME   - MMWPQ1
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1982
C
C   PURPOSE             - FIRST DERIVATIVE OF THE WEIERSTRASS
C                           P-FUNCTION IN THE EQUIANHARMONIC CASE FOR
C                           COMPLEX ARGUMENT WITH UNIT PERIOD
C                           PARALLELOGRAM
C
C   USAGE               - CALL MMWPQ1 (Z,PEQ1,IER)
C
C   ARGUMENTS    Z      - INPUT COMPLEX ARGUMENT OF THE WEIERSTRASS
C                           P-FUNCTION.
C                PEQ1   - OUTPUT COMPLEX DERIVATIVE OF THE WEIERSTRASS
C                           P-FUNCTION IN THE EQUIANHARMONIC CASE.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER = 129 INDICATES THAT INPUT ARGUMENT
C                             Z CORRESPONDS TO A LATTICE POINT AT WHICH
C                             THE WEIERSTRASS P-FUNCTION HAS A POLE.
C                             PEQ1 IS SET TO MACHINE INFINITY.
C
C   REQD. IMSL ROUTINES - UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS      ON THOSE MACHINES WHICH DO NOT OFFER THE APPROPRIATE
C                COMPLEX DATA TYPE, INPUT ARGUMENT Z AND OUTPUT ARGUMENT
C                PEQ1 ARE TREATED AS VECTORS OF LENGTH 2. THE USER
C                SHOULD DIMENSION THEM AS SUCH ON THOSE MACHINES.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE MMWPQ1 (Z,PEQ1,IER)
C
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            IER
      COMPLEX*16         Z,PEQ1
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      DOUBLE PRECISION   ZR,ZI,C1,C2,C3,C4,C5,C6,C7,C8,C9,C10,C11,DM,DN
      DOUBLE PRECISION   C12,C13,C0,C,TZERO
      DOUBLE PRECISION   DINT,DUMMY
      DOUBLE PRECISION   DREAL,DIMAG
      COMPLEX*16         ZDUM
      COMPLEX*16         Z2,Z3,Z4,CXINF,PTMP,Z6
      DATA               C/.8660254037844387D0/
      DATA               C0/1.154700538379252D0/
      DATA               C1/-2.95539175D-9/
      DATA               C2/2.6764693031D-7/
      DATA               C3/2.402192743346D-5/
      DATA               C4/1.9656661451391D-4/
      DATA               C5/1.760135529461036D-2/
      DATA               C6/8.102624349882264D-1/
      DATA               C7/2.739366131491968D0/
      DATA               C8/4.6397763D-10/
      DATA               C9/5.413482233D-8/
      DATA               C10/1.56293298374D-6/
      DATA               C11/1.0393701076352D-4/
      DATA               C12/9.5553182532237D-4/
      DATA               C13/9.131106969640212D-2/
      DINT(DUMMY) = DBLE(FLOAT(IDINT(DUMMY)))
      DREAL(ZDUM) = ZDUM
      DIMAG(ZDUM) = (0.D0,-1.D0)*ZDUM
C                                  FIRST EXECUTABLE STATEMENT
      CXINF = (.723700557733226D+76,0.D0)
      IER = 0
C                                  REDUCTION TO FUNDAMENTAL
C                                    PARALLELOGRAM
      ZI = C0 * DIMAG(Z) + .5D0
      DM = DINT(ZI)
      IF(ZI.LT.0.D0) DM = DM - 1.D0
      ZR = DREAL(Z) - .5D0*DM + .5D0
      DN = DINT(ZR)
      IF (ZR.LT.0.D0) DN = DN - 1.D0
      Z3 = Z - DN - DCMPLX(.5D0,C) * DM
C                                  IF Z3=0 THEN Z COINCIDES WITH A
C                                    LATTICE POINT. SINCE P HAS POLES
C                                    AT THE LATTICE POINTS, A DIVISION
C                                    ERROR WILL OCCUR.
      TZERO = CDABS(Z3)
      IF(TZERO.NE.0.D0)GO TO 5
      IER = 129
      PEQ1 = CXINF
      GO TO 9000
    5 Z3 = Z3*Z3*Z3
      Z6 = Z3*Z3
      PTMP = ((((C1*Z6-C2)*Z6 + C3)*Z6 + C5)*Z6 + C6)*Z6 - C7
      PEQ1 = (((((C8*Z6+C9)*Z6-C10)*Z6-C11)*Z6+C12)*Z6+C13)*Z6+1.D0
      PEQ1 = PTMP/PEQ1 * Z3
      PTMP = 1.D0-Z6
      PTMP = PTMP*PTMP*PTMP
      PEQ1 = PEQ1 + (((14.D0*Z6+294.D0)*Z6+126.D0)*Z6-2.D0)/(Z3*PTMP)
      GO TO 9005
 9000 CONTINUE
      CALL UERTST(IER,6HMMWPQ1 )
 9005 RETURN
      END
