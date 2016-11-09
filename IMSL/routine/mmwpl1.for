C   IMSL ROUTINE NAME   - MMWPL1
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1982
C
C   PURPOSE             - FIRST DERIVATIVE OF THE WEIERSTRASS
C                           P-FUNCTION IN THE LEMNISCATIC CASE FOR
C                           COMPLEX ARGUMENT WITH UNIT PERIOD
C                           PARALLELOGRAM
C
C   USAGE               - CALL MMWPL1 (Z,PLEM1,IER)
C
C   ARGUMENTS    Z      - INPUT COMPLEX ARGUMENT OF THE WEIERSTRASS
C                           P-FUNCTION.
C                PLEM1  - OUTPUT COMPLEX DERIVATIVE OF THE WEIERSTRASS
C                           P-FUNCTION IN THE LEMNISCATE CASE.
C
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER = 129 INDICATES THAT INPUT ARGUMENT
C                             Z CORRESPONDS TO A LATTICE POINT AT WHICH
C                             THE WEIERSTRASS P-FUNCTION HAS A POLE.
C                             PLEM1 IS SET TO MACHINE INFINITY.
C
C   REQD. IMSL ROUTINES - UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS      ON THOSE MACHINES WHICH DO NOT OFFER THE APPROPRIATE
C                COMPLEX DATA TYPE, INPUT ARGUMENT Z AND OUTPUT ARGUMENT
C                PLEM1 ARE TREATED AS VECTORS OF LENGTH 2. THE USER
C                SHOULD DIMENSION THEM AS SUCH ON THOSE MACHINES.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE MMWPL1 (Z,PLEM1,IER)
C
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            IER
      COMPLEX*16         Z,PLEM1
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      DOUBLE PRECISION   ZR,ZI,C1,C2,C3,C4,C5,C6,C7,C8,C9,C10,C11,DM,DN
      DOUBLE PRECISION   C12,C13,C14,C15,C16,C17,C18,C19,C20,C21,TZERO
      DOUBLE PRECISION   DINT,DUMMY
      DOUBLE PRECISION   DREAL,DIMAG
      COMPLEX*16         ZDUM
      COMPLEX*16         Z1,Z3,Z4,CXINF
      DATA               C1/-3.9046302D-9/
      DATA               C2/1.001487137D-8/
      DATA               C3/5.9573043092D-7/
      DATA               C4/2.482518130524D-5/
      DATA               C5/1.4557266595395D-4/
      DATA               C6/4.56633655643206D-3/
      DATA               C7/6.224782572111135D-2/
      DATA               C8/1.038527937794269D-2/
      DATA               C9/1.198046208026379D0/
      DATA               C10/6.427914396838117D0/
      DATA               C11/5.092727987076615D0/
      DATA               C12/4.726888D-11/
      DATA               C13/3.0667983D-9/
      DATA               C14/1.0087596089D-7/
      DATA               C15/8.060683451D-8/
      DATA               C16/1.184299251664D-5/
      DATA               C17/2.3096723361547D-4/
      DATA               C18/2.90730903142055D-3/
      DATA               C19/1.338392411135511D-2/
      DATA               C20/2.309863932002143D-1/
      DATA               C21/8.471988096455415D-1/
      DINT(DUMMY) = DBLE(FLOAT(IDINT(DUMMY)))
      DREAL(ZDUM) = ZDUM
      DIMAG(ZDUM) = (0.D0,-1.D0)*ZDUM
C                                  FIRST EXECUTABLE STATEMENT
      CXINF = (.723700557733226D+76,0.D0)
      IER = 0
C                                  REDUCTION TO FUNDAMENTAL
C                                    PARALLELOGRAM
      ZR = DREAL(Z) + .5D0
      ZI = DIMAG(Z) + .5D0
      DM = DINT(ZR)
      DN = DINT(ZI)
      IF (ZR.LT.0.D0) DM = DM - 1.D0
      IF (ZI.LT.0.D0) DN = DN - 1.D0
      Z1 = Z - DM - (0.D0,1.D0)*DN
C                                  IF Z1=0 THEN Z COINCIDES WITH A
C                                    LATTICE POINT. SINCE P HAS POLES
C                                    AT THE LATTICE POINTS, A DIVISION
C                                    ERROR WILL OCCUR
      TZERO = CDABS(Z1)
      IF(TZERO.NE.0.D0) GO TO 5
      PLEM1 = CXINF
      IER = 129
      GO TO 9000
    5 Z3 = Z1*Z1*Z1
      Z4 = Z3*Z1
      Z3 = (((1.D1*Z4+9.D1)*Z4+3.D1)*Z4-2.D0)/(Z1*(1.D0-Z4))**3
      PLEM1 = (((((((((C1*Z4-C2)*Z4 + C3)*Z4 - C4)*Z4 + C5)*Z4 + C6)*Z4
     *         + C7)*Z4 + C8)*Z4 + C9)*Z4 + C10)*Z4 - C11
      PLEM1 = Z1*PLEM1
      Z1 = (((((((((C12*Z4-C13)*Z4 + C14)*Z4 - C15)*Z4 + C16)*Z4 -
     *      C17)*Z4 - C18)*Z4 + C19)*Z4 + C20)*Z4 + C21)*Z4 + 1.D0
      PLEM1 = Z3 + PLEM1/Z1
      GO TO 9005
 9000 CONTINUE
      CALL UERTST(IER,6HMMWPL1)
 9005 RETURN
      END
