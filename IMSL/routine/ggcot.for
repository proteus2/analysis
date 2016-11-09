C   IMSL ROUTINE NAME   - GGCOT
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1982
C
C   PURPOSE             - NUCLEUS CALLED ONLY BY IMSL ROUTINE GGCOR
C
C   REQD. IMSL ROUTINES - NONE REQUIRED
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE GGCOT  (A,B,C,ZSM,ZLG,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            IER
      REAL               A,B,C
      COMPLEX            ZSM,ZLG
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            IS
      REAL               RADIX,RNLGRX,SCALE
      REAL               HALF,S,ZERO,FINITY,A0,B0,C0,B1,DD
      DOUBLE PRECISION   D,D1
      COMPLEX            ZS,ZL
      DATA               FINITY/Z7FFFFFFF/
      DATA               RADIX/16.0/
      DATA               RNLGRX/2.772589/
      DATA               ZERO/0.0/,HALF/0.5/
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
C                                  PUT THE COEFFICIENTS IN TEMPORARY TO
C                                    SAVE EXECUTION TIME.
      A0 = A
      B1 = -B
      C0 = C
C                                  CHECK FOR A=ZERO OR C=ZERO.
      IF(A0 .NE. ZERO) GO TO 5
      IER = 65
      ZL = CMPLX(FINITY,ZERO)
      ZS = -ZL
      IF(B1 .EQ. ZERO) GO TO 30
      IER = 66
      ZL = CMPLX(SIGN(FINITY,B1),ZERO)
      ZS = CMPLX(C0/B1,ZERO)
      GO TO 30
    5 IF(C0 .NE. ZERO) GO TO 10
      ZS = CMPLX(ZERO,ZERO)
      GO TO 25
C                                  SCALING TO AVOID OVERFLOW OR
C                                    UNDERFLOW. SCALE THE COEFFICIENTS
C                                    SO THAT A*C IS APPROXIMATELY ONE.
C                                    THE SCALE FACTOR SQRT(A*C) FITS
C                                    THIS REQUIREMENT BUT MAY CAUSE
C                                    OVERFLOW OR UNDERFLOW IN THE
C                                    SCALING PROCEDURE.
C                                    LET A=RADIX**IA AND C=RADIX**IC.
C                                    THE SCALE FACTOR, SCALE, IS DEFINED
C                                    BY THE FOLLOWING FORMULA,
C                                    SCALE=RADIX**IS, WHERE
C                                    IS=ENTIER((IA+IC+1)/2) AND
C                                    ENTIER IS THE MATHEMATICAL GREATEST
C                                    INTEGER FUNCTION.
   10 IS = (ALOG(ABS(A0))+ALOG(ABS(C0))+RNLGRX)/(RNLGRX+RNLGRX)
      SCALE = RADIX**IS
C                                  IF THE SCALE FACTOR .LE.
C                                    DEPS*ABS(B1) DO NOT SCALE
C                                    THE COEFFICIENTS.
      D1 = DBLE(ABS(B1))
      D = D1+SCALE
      D = D-D1
      IF (SNGL(D) .EQ. ZERO) GO TO 20
C                                  IF ABS(B1) .GE. DEPS*SCALE FACTOR
C                                    THEN SCALE B0. OTHERWISE SET
C                                    B0 = ZERO.
      B0 = ZERO
      D = D1+SCALE
      D = D-SCALE
      IF (SNGL(D) .NE. ZERO) B0 = (B1/SCALE)*HALF
      A0 = A0/SCALE
      C0 = C0/SCALE
C                                  SOLVE A0*Z**2-2.0*B0*Z+C0=ZERO
      DD = DBLE(B0)**2-DBLE(A0)*DBLE(C0)
      S = SQRT(ABS(DD))
      IF(DD .GT. ZERO) GO TO 15
C                                  COINCIDENT OR COMPLEX ROOTS
C                                    (D .LE. ZERO).
      ZL = CMPLX(B0/A0,ABS(S/A0))
      ZS = CONJG(ZL)
      GO TO 30
C                                  DISTINCT REAL ROOTS (D .GT. ZERO).
   15 B1 = SIGN(S,B0)+B0
   20 ZS = CMPLX(C0/B1,ZERO)
   25 ZL = CMPLX(B1/A0,ZERO)
      IF(ABS(REAL(ZL)) .LT. ABS(REAL(ZS))) ZS = -ZL
   30 ZSM = ZS
      ZLG = ZL
      RETURN
      END
