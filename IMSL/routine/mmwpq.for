C   IMSL ROUTINE NAME   - MMWPQ
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1982
C
C   PURPOSE             - WEIERSTRASS P-FUNCTION IN THE EQUIANHARMONIC
C                           CASE FOR COMPLEX ARGUMENT WITH UNIT PERIOD
C                           PARALLELOGRAM
C
C   USAGE               - CALL MMWPQ (Z,PEQ,IER)
C
C   ARGUMENTS    Z      - INPUT COMPLEX ARGUMENT OF THE WEIERSTRASS
C                           P-FUNCTION.
C                PEQ    - OUTPUT COMPLEX WEIERSTRASS P-FUNCTION IN THE
C                           EQUIANHARMONIC CASE.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER = 129 INDICATES THAT INPUT ARGUMENT
C                             Z CORRESPONDS TO A LATTICE POINT AT WHICH
C                             THE WEIERSTRASS P-FUNCTION HAS A POLE.
C                             PEQ IS SET TO MACHINE INFINITY.
C
C   REQD. IMSL ROUTINES - UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS      ON THOSE MACHINES WHICH DO NOT OFFER THE APPROPRIATE
C                COMPLEX DATA TYPE, INPUT ARGUMENT Z AND OUTPUT ARGUMENT
C                PEQ ARE TREATED AS VECTORS OF LENGTH 2. THE USER
C                SHOULD DIMENSION THEM AS SUCH ON THOSE MACHINES.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE MMWPQ (Z,PEQ,IER)
C
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            IER
      COMPLEX*16         Z,PEQ
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      DOUBLE PRECISION   ZR,ZI,C1,C2,C3,C4,C5,C6,C7,C8,C9,C10,C11,DM,DN
      DOUBLE PRECISION   C,C0,TZERO
      DOUBLE PRECISION   DINT,DUMMY
      DOUBLE PRECISION   DREAL,DIMAG
      COMPLEX*16         ZDUM
      COMPLEX*16         Z2,Z4,CXINF,PTMP,Z6
      DATA               C/.8660254037844387D0/
      DATA               C0/1.154700538379252D0/
      DATA               C1/-2.6427662D-10/
      DATA               C2/1.610954818D-8/
      DATA               C3/7.38610752879D-6/
      DATA               C4/4.3991444671178D-4/
      DATA               C5/7.477288220490697D-2/
      DATA               C6/6.848415328729920D-1/
      DATA               C7/6.2252191D-10/
      DATA               C8/2.553314573D-7/
      DATA               C9/2.619832920421D-5/
      DATA               C10/5.6444801847646D-4/
      DATA               C11/4.565553484820106D-2/
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
      ZR = DREAL(Z) - 0.5D0*DM + 0.5D0
      DN = DINT(ZR)
      IF (ZR.LT.0.D0) DN = DN - 1.D0
      Z2 = Z - DN - DCMPLX(.5D0,C) * DM
C                                  IF Z2=0 THEN Z COINCIDES WITH A
C                                    LATTICE POINT. SINCE P HAS POLES
C                                    AT THE LATTICE POINTS, A DIVISION
C                                    ERROR WILL OCCUR.
      TZERO = CDABS(Z2)
      IF (TZERO.NE.0.D0) GO TO 5
      IER = 129
      PEQ = CXINF
      GO TO 9000
    5 Z2 = Z2*Z2
      Z4 = Z2*Z2
      Z6 = Z4*Z2
      PEQ = ((((C1*Z6+C2)*Z6 + C3)*Z6 + C4)*Z6 + C5)*Z6 - C6
      PTMP = ((((C7*Z6+C8)*Z6 - C9)*Z6 - C10)*Z6 + C11)*Z6 + 1.D0
      PTMP = PEQ/PTMP * Z4
      PEQ = 1.D0/Z2 + 6.D0*Z4*(5.D0+Z6)/(1.D0-Z6)**2 + PTMP
      GO TO 9005
 9000 CONTINUE
      CALL UERTST(IER,6HMMWPQ )
 9005 RETURN
      END
