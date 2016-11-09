C   IMSL ROUTINE NAME   - ZQADC
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - ZEROS OF A QUADRATIC WITH COMPLEX COEFFICIENTS
C
C   USAGE               - CALL ZQADC (A,B,C,ZSM,ZLG,IER)
C
C   ARGUMENTS    A      - COEFFICIENT OF THE QUADRATIC EQUATION. (INPUT)
C                B      - COEFFICIENT OF THE QUADRATIC EQUATION. (INPUT)
C                C      - COEFFICIENT OF THE QUADRATIC EQUATION. (INPUT)
C                         A, B, AND C MUST BE DECLARED TYPE COMPLEX.
C                         THE QUADRATIC EQUATION IS OF THE FORM
C                           A*Z**2+B*Z+C = 0.0
C                ZSM    - ROOT OF THE QUADRATIC EQUATION. (OUTPUT)
C                ZLG    - ROOT OF THE QUADRATIC EQUATION. (OUTPUT)
C                         ZSM AND ZLG MUST BE DECLARED TYPE COMPLEX.
C                           FOR THE ROOTS ZSM AND ZLG THE FOLLOWING
C                             CONDITION HOLDS
C                             CABS(ZSM) .LE. CABS(ZLG)
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         WARNING (WITH FIX)
C                           IER = 65, IMPLIES A=B=0.0
C                             IN THIS CASE, THE LARGE ROOT,
C                             ZLG = CMPLX(FINITY,0.0), AND
C                             THE SMALL ROOT, ZSM = -ZLG, WHERE
C                             FINITY = LARGEST NUMBER WHICH CAN BE
C                             REPRESENTED IN THE MACHINE.
C                           IER = 66, IMPLIES A=0.0
C                             IN THIS CASE, THE LARGE ROOT,
C                             ZLG = CMPLX(FINITY,0.0), WHERE
C                             FINITY = LARGEST NUMBER WHICH CAN BE
C                             REPRESENTED IN THE MACHINE.
C
C   REQD. IMSL ROUTINES - SINGLE/UERTST,UGETIO
C                       - DOUBLE/UERTST,UGETIO,VXADD,VXMUL,VXSTO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE ZQADC  (A,B,C,ZSM,ZLG,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            IER
      COMPLEX            A,B,C,ZSM,ZLG
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            IS
      REAL               RADIX,RNLGRX,SCALE
      REAL               AA(2),BB(2),CC(2),ZSA(2),ZLA(2),ZSR,ZLR,ZSI
      REAL               ZLI,HALF,ZERO,FINITY,AR,BR,CR,AI,BI,CI,TEMP
      DOUBLE PRECISION   DR,DI,D,D1
      COMPLEX            A0,B0,C0,ZS,ZL
      EQUIVALENCE        (AA(1),A0),(BB(1),B0),(CC(1),C0),(ZSA(1),ZS),
     1                   (ZLA(1),ZL),(AA(1),AR),(BB(1),BR),(CC(1),CR),
     2                   (ZSA(1),ZSR),(ZLA(1),ZLR),(AA(2),AI),
     3                   (BB(2),BI),(CC(2),CI),(ZSA(2),ZSI),(ZLA(2),ZLI)
      DATA               FINITY/Z7FFFFFFF/
      DATA               RADIX/16.0/
      DATA               RNLGRX/2.772589/
      DATA               ZERO/0.0/,HALF/0.5/
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
C                                  PUT THE COEFFICIENTS IN TEMPORARY TO
C                                    SAVE EXECUTION TIME.
      A0 = A
      B0 = -B
      C0 = C
C                                  CHECK FOR A=ZERO OR C=ZERO.
      IF(AR .NE. ZERO .OR. AI .NE. ZERO) GO TO 5
      IER = 65
      ZL = CMPLX(FINITY,ZERO)
      ZS = -ZL
      IF(BR .EQ. ZERO .AND. BI .EQ. ZERO) GO TO 35
      IER = 66
      ZS = C0/B0
      GO TO 35
    5 IF(CR .NE. ZERO .OR. CI .NE. ZERO) GO TO 10
      ZS = CMPLX(ZERO,ZERO)
      GO TO 30
C                                  SCALING TO AVOID OVERFLOW OR
C                                    UNDERFLOW. SCALE THE COEFFICIENTS
C                                    SO THAT A*C IS APPROXIMATELY ONE.
C                                    THE SCALE FACTOR CSQRT(A*C) FITS
C                                    THIS REQUIREMENT BUT MAY CAUSE
C                                    OVERFLOW OR UNDERFLOW IN THE
C                                    SCALING PROCEDURE.
C                                    LET AMAX1(ABS(AR),ABS(AI)) BE
C                                    REPRESENTED BY RADIX**IA AND LET
C                                    AMAX1(ABS(CR),ABS(CI)) BE
C                                    REPRESENTED BY RADIX**IC.
C                                    THE SCALE FACTOR, SCALE, IS DEFINED
C                                    BY THE FOLLOWING FORMULA,
C                                    SCALE=RADIX**IS, WHERE
C                                    IS=ENTIER((IA+IC+1)/2) AND
C                                    ENTIER IS THE MATHEMATICAL GREATEST
C                                    INTEGER FUNCTION.
   10 IS = (ALOG(AMAX1(ABS(AR),ABS(AI)))+ALOG(AMAX1(ABS(CR),ABS(CI)))+
     1      RNLGRX)/(RNLGRX+RNLGRX)
      SCALE = RADIX**IS
C                                  IF THE SCALE FACTOR .LE.
C                                    DEPS*MAX(ABS(BR),ABS(BI))
C                                    DO NOT SCALE THE COEFFICIENTS.
      TEMP = AMAX1(ABS(BR),ABS(BI))
      D1 = DBLE(TEMP)
      D = D1+SCALE
      D = D-D1
      IF (SNGL(D) .EQ. ZERO) GO TO 25
C                                  IF MAX(ABS(BR),ABS(BI)) .GE.
C                                    DEPS*SCALE FACTOR THEN SCALE
C                                    B0. OTHERWISE SET B0 = ZERO.
      D = D1+SCALE
      D = D-SCALE
      IF (SNGL(D) .NE. ZERO) GO TO 15
      BR = ZERO
      BI = ZERO
      GO TO 20
   15 BR = (BR/SCALE)*HALF
      BI = (BI/SCALE)*HALF
   20 AR = AR/SCALE
      AI = AI/SCALE
      CR = CR/SCALE
      CI = CI/SCALE
C                                  SOLVE A0*Z**2-2.0*B0*Z+C0=ZERO
      DR = DBLE(BR)**2
      DI = DBLE(BI)*(2.0D0*DBLE(BR))
      ZS = CMPLX(SNGL(((DR-DBLE(BI)**2)-DBLE(AR)*DBLE(CR))+DBLE(AI)*
     1        DBLE(CI)),SNGL((DI-DBLE(AI)*DBLE(CR))-DBLE(AR)*DBLE(CI)))
      ZS = CSQRT(ZS)
C                                  CHOOSE THE SIGN OF ZS SUCH THAT
C                                  CABS(B)=AMAX1(CABS(B+ZS),CABS(B-ZS)).
      IF(DBLE(ZSR)*DBLE(BR)+DBLE(ZSI)*DBLE(BI) .LE. ZERO) ZS = -ZS
      B0 = B0+ZS
C                                  PERFORM THE FINAL COMPLEX OPERATION
C                                    FOR THE ZEROS.
   25 ZS = C0/B0
   30 ZL = B0/A0
   35 ZSM = ZS
      ZLG = ZL
 9000 CONTINUE
      IF(IER .NE. 0) CALL UERTST(IER,6HZQADC )
 9005 RETURN
      END
