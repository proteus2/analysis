C   IMSL ROUTINE NAME   - MMWPL
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1982
C
C   PURPOSE             - WEIERSTRASS P-FUNCTION IN THE LEMNISCATIC
C                           CASE FOR COMPLEX ARGUMENT WITH UNIT PERIOD
C                           PARALLELOGRAM
C
C   USAGE               - CALL MMWPL (Z,PLEM,IER)
C
C   ARGUMENTS    Z      - INPUT COMPLEX ARGUMENT OF THE WEIERSTRASS
C                           P-FUNCTION.
C                PLEM   - OUTPUT COMPLEX WEIERSTRASS P-FUNCTION IN THE
C                           LEMNISCATE CASE.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER = 129 INDICATES THAT INPUT ARGUMENT
C                             Z CORRESPONDS TO A LATTICE POINT AT WHICH
C                             THE WEIERSTRASS P-FUNCTION HAS A POLE.
C                             PLEM IS SET TO MACHINE INFINITY.
C
C   REQD. IMSL ROUTINES - UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS      ON THOSE MACHINES WHICH DO NOT OFFER THE APPROPRIATE
C                COMPLEX DATA TYPE, INPUT ARGUMENT Z AND OUTPUT ARGUMENT
C                PLEM ARE TREATED AS VECTORS OF LENGTH 2. THE USER
C                SHOULD DIMENSION THEM AS SUCH ON THOSE MACHINES.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE MMWPL (Z,PLEM,IER)
C
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            IER
      COMPLEX*16         Z,PLEM
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      DOUBLE PRECISION   ZR,ZI,C1,C2,C3,C4,C5,C6,C7,C8,C9,C10,C11,DM,DN
      DOUBLE PRECISION   C12,C13,C14,C15,C16,C17,TZERO
      DOUBLE PRECISION   DINT,DUMMY
      DOUBLE PRECISION   DREAL,DIMAG
      COMPLEX*16         ZDUM
      COMPLEX*16         Z2,Z4,CXINF,PTMP
      DATA               C1/7.233108D-11/
      DATA               C2/1.714197273D-8/
      DATA               C3/2.5369036492D-7/
      DATA               C4/7.98710206868D-6/
      DATA               C5/6.4850606909737D-4/
      DATA               C6/7.39624629362938D-3/
      DATA               C7/2.012382768497244D-2/
      DATA               C8/7.117729754313660D-1/
      DATA               C9/2.546363993538307D0/
      DATA               C10/5.1161516D-10/
      DATA               C11/6.61289408D-9/
      DATA               C12/4.4618987048D-7/
      DATA               C13/8.42694918892D-6/
      DATA               C14/4.42886829095D-6/
      DATA               C15/4.22629935217101D-3/
      DATA               C16/2.577496871700433D-2/
      DATA               C17/4.235994048227707D-1/
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
      Z2 = Z - DM - (0.D0,1.D0)*DN
C                                  IF Z2=0 THEN Z COINCIDES WITH A
C                                    LATTICE POINT. SINCE P HAS POLES
C                                    AT THE LATTICE POINTS, A DIVISION
C                                    ERROR WILL OCCUR
      TZERO = CDABS(Z2)
      IF(TZERO.NE.0.D0) GO TO 5
      PLEM  = CXINF
      IER = 129
      GO TO 9000
    5 Z2 = Z2*Z2
      Z4 = Z2*Z2
      PLEM = (((((((-C1*Z4+C2)*Z4 - C3)*Z4 - C4)*Z4 + C5)*Z4 + C6)*Z4
     *        + C7)*Z4 + C8)*Z4 - C9
      PLEM = Z2 * PLEM
      PTMP = (((((((C10*Z4 + C11)*Z4 + C12)*Z4 - C13)*Z4 + C14)*Z4
     *        - C15)*Z4 + C16)*Z4 + C17)*Z4 + 1.D0
      PLEM = PLEM/PTMP
      PLEM = 1.D0/Z2 + 4.D0*Z2*(3.D0+Z4)/((1.D0-Z4)*(1.D0-Z4))+PLEM
      GO TO 9005
 9000 CONTINUE
      CALL UERTST(IER,6HMMWPL )
 9005 RETURN
      END
