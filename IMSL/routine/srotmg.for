C   IMSL ROUTINE NAME   - VBLA=SROTMG
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1980
C
C   PURPOSE             - CONSTRUCT A MODIFIED GIVENS PLANE ROTATION
C                           (SINGLE PRECISION)
C
C   USAGE               - CALL SROTMG (SD1,SD2,SX1,SY1,SPARAM)
C
C   ARGUMENTS    SD1    - ON INPUT, FIRST SCALE FACTOR.
C                         ON OUTPUT, SD1 IS REPLACED WITH THE UPDATE
C                           SCALE FACTOR.
C                SD2    - ON INPUT, SECOND SCALE FACTOR.
C                         ON OUTPUT, SD2 IS REPLACED WITH THE UPDATE
C                           SCALE FACTOR.
C                SX1    - ON INPUT, FIRST COMPONENT OF VECTOR TO BE
C                           TRANSFORMED.
C                         ON OUTPUT, SX1 IS REPLACED WITH ITS
C                           TRANSFORMED VALUE.
C                SY1    - ON INPUT, SECOND COMPONENT OF VECTOR TO BE
C                           TRANSFORMED. SINCE THIS COMPONENT IS
C                           ZEROED BY THE TRANSFORMATION, IT IS LEFT
C                           UNCHANGED IN STORAGE.
C                SPARAM - REAL VECTOR OF LENGTH 5 WHICH DEFINES THE
C                           TRANSFORMATION MATRIX H. SEE REMARKS.
C                           (OUTPUT)
C
C   REQD. IMSL ROUTINES - NONE REQUIRED
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS      SROTMG CONSTRUCTS THE MODIFIED GIVENS TRANSFORMATION H
C                AND UPDATED SCALE FACTORS (SD1, AND SD2) WHICH ZERO
C                SY1. THE TRANSFORMED VALUE OF SD1 REPLACES SD1 IN
C                STORAGE. THAT IS,
C                ON INPUT,  SW1 = SQRT(SD1)*SX1
C                           SZ1 = SQRT(SD2)*SY1.
C                ON OUTPUT, (C  S)(SW1)=(C*SW1+S*SZ1)=(SQRT(SD1)*SX1)
C                           (-S C)(SZ1) (     0     ) (      0      )
C                WHERE C AND S DEFINE A GIVENS ROTATION.
C
C                IF SPARAM(1)=-2.0
C                  SPARAM(2) = UNCHANGED  SPARAM(4) = UNCHANGED
C                  SPARAM(3) = UNCHANGED  SPARAM(5) = UNCHANGED
C                  H IS THE IDENTITY MATRIX, BUT IT IS NOT STORED.
C                IF SPARAM(1)=-1.0
C                  SPARAM(2) = H11        SPARAM(4) = H12
C                  SPARAM(3) = H21        SPARAM(5) = H22
C                IF SPARAM(1)=0.0
C                  SPARAM(2) = UNCHANGED  SPARAM(4) = H12
C                  SPARAM(3) = H21        SPARAM(5) = UNCHANGED
C                  H11=1.0 AND H22=1.0, BUT THEY ARE NOT STORED.
C                IF SPARAM(1)=1.0
C                  SPARAM(2) = H11        SPARAM(4) = UNCHANGED
C                  SPARAM(3) = UNCHANGED  SPARAM(5) = H22
C                  H12=1.0 AND H21=-1.0, BUT THEY ARE NOT STORED.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE SROTMG (SD1,SD2,SX1,SY1,SPARAM)
C
C                                  SPECIFICATIONS FOR ARGUMENTS
      REAL               SPARAM(5),SD1,SD2,SX1,SY1
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            IGO
      REAL               GAM,GAMSQ,ONE,RGAMSQ,SFLAG,SH11,SH12,
     *                   SH21,SH22,SP1,SP2,SQ1,SQ2,STEMP,SU,TWO,ZERO
      DATA ZERO,ONE,TWO /0.E0,1.E0,2.E0/
      DATA GAM,GAMSQ,RGAMSQ/4096.E0,1.678E7,5.960E-8/
C                                  FIRST EXECUTABLE STATEMENT
      IF (.NOT.SD1.LT.ZERO) GO TO 5
C                                  GO ZERO-H-D-AND-SX1..
      GO TO 30
    5 CONTINUE
C                                  CASE-SD1-NONNEGATIVE
      SP2 = SD2*SY1
      IF (.NOT.SP2.EQ.ZERO) GO TO 10
      SFLAG = -TWO
      GO TO 135
C                                  REGULAR-CASE..
   10 CONTINUE
      SP1 = SD1*SX1
      SQ2 = SP2*SY1
      SQ1 = SP1*SX1
C
      IF (.NOT.ABS(SQ1).GT.ABS(SQ2)) GO TO 20
      SH21 = -SY1/SX1
      SH12 = SP2/SP1
C
      SU = ONE-SH12*SH21
C
      IF (.NOT.SU.LE.ZERO) GO TO 15
C                                  GO ZERO-H-D-AND-SX1..
      GO TO 30
   15 CONTINUE
      SFLAG = ZERO
      SD1 = SD1/SU
      SD2 = SD2/SU
      SX1 = SX1*SU
C                                  GO SCALE-CHECK..
      GO TO 50
   20 CONTINUE
      IF (.NOT.SQ2.LT.ZERO) GO TO 25
C                                  GO ZERO-H-D-AND-SX1..
      GO TO 30
   25 CONTINUE
      SFLAG = ONE
      SH11 = SP1/SP2
      SH22 = SX1/SY1
      SU = ONE+SH11*SH22
      STEMP = SD2/SU
      SD2 = SD1/SU
      SD1 = STEMP
      SX1 = SY1*SU
C                                  GO SCALE-CHECK
      GO TO 50
C                                  PROCEDURE..ZERO-H-D-AND-SX1..
   30 CONTINUE
      SFLAG = -ONE
      SH11 = ZERO
      SH12 = ZERO
      SH21 = ZERO
      SH22 = ZERO
C
      SD1 = ZERO
      SD2 = ZERO
      SX1 = ZERO
C                                  RETURN..
      GO TO 115
C                                  PROCEDURE..FIX-H..
   35 CONTINUE
      IF (.NOT.SFLAG.GE.ZERO) GO TO 45
C
      IF (.NOT.SFLAG.EQ.ZERO) GO TO 40
      SH11 = ONE
      SH22 = ONE
      SFLAG = -ONE
      GO TO 45
   40 CONTINUE
      SH21 = -ONE
      SH12 = ONE
      SFLAG = -ONE
   45 CONTINUE
      GO TO IGO, (65,80,95,110)
C                                  PROCEDURE..SCALE-CHECK
   50 CONTINUE
   55 CONTINUE
      IF (.NOT.SD1.LE.RGAMSQ) GO TO 70
   60 CONTINUE
      IF (SD1.EQ.ZERO) GO TO 85
      ASSIGN 65 TO IGO
C                                  FIX-H..
      GO TO 35
   65 CONTINUE
      SD1 = SD1*GAM**2
      SX1 = SX1/GAM
      SH11 = SH11/GAM
      SH12 = SH12/GAM
      GO TO 55
   70 CONTINUE
   75 CONTINUE
      IF (.NOT.SD1.GE.GAMSQ) GO TO 85
      ASSIGN 80 TO IGO
C                                  FIX-H..
      GO TO 35
   80 CONTINUE
      SD1 = SD1/GAM**2
      SX1 = SX1*GAM
      SH11 = SH11*GAM
      SH12 = SH12*GAM
      GO TO 75
   85 CONTINUE
   90 CONTINUE
      IF (.NOT.ABS(SD2).LE.RGAMSQ) GO TO 100
      IF (SD2.EQ.ZERO) GO TO 115
      ASSIGN 95 TO IGO
C                                  FIX-H..
      GO TO 35
   95 CONTINUE
      SD2 = SD2*GAM**2
      SH21 = SH21/GAM
      SH22 = SH22/GAM
      GO TO 90
  100 CONTINUE
  105 CONTINUE
      IF (.NOT.ABS(SD2).GE.GAMSQ) GO TO 115
      ASSIGN 110 TO IGO
C                                  FIX-H..
      GO TO 35
  110 CONTINUE
      SD2 = SD2/GAM**2
      SH21 = SH21*GAM
      SH22 = SH22*GAM
      GO TO 105
  115 CONTINUE
      IF (SFLAG) 130,120,125
  120 CONTINUE
      SPARAM(3) = SH21
      SPARAM(4) = SH12
      GO TO 135
  125 CONTINUE
      SPARAM(2) = SH11
      SPARAM(5) = SH22
      GO TO 135
  130 CONTINUE
      SPARAM(2) = SH11
      SPARAM(3) = SH21
      SPARAM(4) = SH12
      SPARAM(5) = SH22
  135 CONTINUE
      SPARAM(1) = SFLAG
      RETURN
      END
