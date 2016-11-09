C   IMSL ROUTINE NAME   - VBLA=DROTMG
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1980
C
C   PURPOSE             - CONSTRUCT A MODIFIED GIVENS PLANE ROTATION
C                           (DOUBLE PRECISION)
C
C   USAGE               - CALL DROTMG (DD1,DD2,DX1,DY1,DPARAM)
C
C   ARGUMENTS    DD1    - ON INPUT, FIRST SCALE FACTOR.
C                         ON OUTPUT, DD1 IS REPLACED WITH THE UPDATE
C                           SCALE FACTOR.
C                DD2    - ON INPUT, SECOND SCALE FACTOR.
C                         ON OUTPUT, DD2 IS REPLACED WITH THE UPDATE
C                           SCALE FACTOR.
C                DX1    - ON INPUT, FIRST COMPONENT OF VECTOR TO BE
C                           TRANSFORMED.
C                         ON OUTPUT, DX1 IS REPLACED WITH ITS
C                           TRANSFORMED VALUE.
C                DY1    - ON INPUT, SECOND COMPONENT OF VECTOR TO BE
C                           TRANSFORMED. SINCE THIS COMPONENT IS
C                           ZEROED BY THE TRANSFORMATION, IT IS LEFT
C                           UNCHANGED IN STORAGE.
C                DPARAM - DOUBLE PRECISION VECTOR OF LENGTH 5 WHICH
C                           DEFINES THE TRANSFORMATION MATRIX H. SEE
C                           REMARKS. (OUTPUT)
C
C   REQD. IMSL ROUTINES - NONE REQUIRED
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS      DROTMG CONSTRUCTS THE MODIFIED GIVENS TRANSFORMATION H
C                AND UPDATED SCALE FACTORS (DD1, AND DD2) WHICH ZERO
C                DY1. THE TRANSFORMED VALUE OF DD1 REPLACES DD1 IN
C                STORAGE. THAT IS,
C                ON INPUT,  DW1 = DSQRT(DD1)*DX1
C                           DZ1 = DSQRT(DD2)*DY1.
C                ON OUTPUT, (C  S)(DW1)=(C*DW1+S*DZ1)=(DSQRT(DD1)*DX1)
C                           (-S C)(DZ1) (     0     ) (       0      )
C                WHERE C AND S DEFINE A GIVENS ROTATION.
C                H TAKES THE FORM,
C
C                DPARAM(1)=-2.0
C                  DPARAM(2) = UNCHANGED  DPARAM(4) = UNCHANGED
C                  DPARAM(3) = UNCHANGED  DPARAM(5) = UNCHANGED
C                DPARAM(1)=-1.0
C                  DPARAM(2) = H11        DPARAM(4) = H12
C                  DPARAM(3) = H21        DPARAM(5) = H22
C                DPARAM(1)=0.0
C                  DPARAM(2) = UNCHANGED  DPARAM(4) = H12
C                  DPARAM(3) = H21        DPARAM(5) = UNCHANGED
C                DPARAM(1)=1.0
C                  DPARAM(2) = H11        DPARAM(4) = UNCHANGED
C                  DPARAM(3) = UNCHANGED  DPARAM(5) = H22
C
C-----------------------------------------------------------------------
C
      SUBROUTINE DROTMG (DD1,DD2,DX1,DY1,DPARAM)
C
C                                  SPECIFICATIONS FOR ARGUMENTS
      DOUBLE PRECISION   DD1,DD2,DX1,DY1,DPARAM(5)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      DOUBLE PRECISION   GAM,ONE,RGAMSQ,DH11,DH21,DP2,DQ2,DU,ZERO,
     1                   GAMSQ,DFLAG,DH12,DH22,DP1,DQ1,DTEMP,TWO
      DATA ZERO,ONE,TWO /0.D0,1.D0,2.D0/
      DATA GAM,GAMSQ,RGAMSQ/4096.D0,1.678D7,5.960D-8/
C                                  FIRST EXECUTABLE STATEMENT
      IF (.NOT.DD1.LT.ZERO) GO TO 5
C                                  GO ZERO-H-D-AND-DX1..
      GO TO 30
    5 CONTINUE
C                                  CASE-DD1-NONNEGATIVE
      DP2 = DD2*DY1
      IF (.NOT.DP2.EQ.ZERO) GO TO 10
      DFLAG = -TWO
      GO TO 135
C                                  REGULAR-CASE..
   10 CONTINUE
      DP1 = DD1*DX1
      DQ2 = DP2*DY1
      DQ1 = DP1*DX1
C
      IF (.NOT.DABS(DQ1).GT.DABS(DQ2)) GO TO 20
      DH21 = -DY1/DX1
      DH12 = DP2/DP1
C
      DU = ONE-DH12*DH21
C
      IF (.NOT.DU.LE.ZERO) GO TO 15
C                                  GO ZERO-H-D-AND-DX1..
      GO TO 30
   15 CONTINUE
      DFLAG = ZERO
      DD1 = DD1/DU
      DD2 = DD2/DU
      DX1 = DX1*DU
C                                  GO SCALE-CHECK..
      GO TO 50
   20 CONTINUE
      IF (.NOT.DQ2.LT.ZERO) GO TO 25
C                                  GO ZERO-H-D-AND-DX1..
      GO TO 30
   25 CONTINUE
      DFLAG = ONE
      DH11 = DP1/DP2
      DH22 = DX1/DY1
      DU = ONE+DH11*DH22
      DTEMP = DD2/DU
      DD2 = DD1/DU
      DD1 = DTEMP
      DX1 = DY1*DU
C                                  GO SCALE-CHECK
      GO TO 50
C                                  PROCEDURE..ZERO-H-D-AND-DX1..
   30 CONTINUE
      DFLAG = -ONE
      DH11 = ZERO
      DH12 = ZERO
      DH21 = ZERO
      DH22 = ZERO
C
      DD1 = ZERO
      DD2 = ZERO
      DX1 = ZERO
C                                  RETURN..
      GO TO 115
C                                  PROCEDURE..FIX-H..
   35 CONTINUE
      IF (.NOT.DFLAG.GE.ZERO) GO TO 45
C
      IF (.NOT.DFLAG.EQ.ZERO) GO TO 40
      DH11 = ONE
      DH22 = ONE
      DFLAG = -ONE
      GO TO 45
   40 CONTINUE
      DH21 = -ONE
      DH12 = ONE
      DFLAG = -ONE
   45 CONTINUE
      GO TO IGO, (65,80,95,110)
C                                  PROCEDURE..SCALE-CHECK
   50 CONTINUE
   55 CONTINUE
      IF (.NOT.DD1.LE.RGAMSQ) GO TO 70
   60 CONTINUE
      IF (DD1.EQ.ZERO) GO TO 85
      ASSIGN 65 TO IGO
C                                  FIX-H..
      GO TO 35
   65 CONTINUE
      DD1 = DD1*GAM**2
      DX1 = DX1/GAM
      DH11 = DH11/GAM
      DH12 = DH12/GAM
      GO TO 55
   70 CONTINUE
   75 CONTINUE
      IF (.NOT.DD1.GE.GAMSQ) GO TO 85
      ASSIGN 80 TO IGO
C                                  FIX-H..
      GO TO 35
   80 CONTINUE
      DD1 = DD1/GAM**2
      DX1 = DX1*GAM
      DH11 = DH11*GAM
      DH12 = DH12*GAM
      GO TO 75
   85 CONTINUE
   90 CONTINUE
      IF (.NOT.DABS(DD2).LE.RGAMSQ) GO TO 100
      IF (DD2.EQ.ZERO) GO TO 115
      ASSIGN 95 TO IGO
C                                  FIX-H..
      GO TO 35
   95 CONTINUE
      DD2 = DD2*GAM**2
      DH21 = DH21/GAM
      DH22 = DH22/GAM
      GO TO 90
  100 CONTINUE
  105 CONTINUE
      IF (.NOT.DABS(DD2).GE.GAMSQ) GO TO 115
      ASSIGN 110 TO IGO
C                                  FIX-H..
      GO TO 35
  110 CONTINUE
      DD2 = DD2/GAM**2
      DH21 = DH21*GAM
      DH22 = DH22*GAM
      GO TO 105
  115 CONTINUE
      IF (DFLAG) 130,120,125
  120 CONTINUE
      DPARAM(3) = DH21
      DPARAM(4) = DH12
      GO TO 135
  125 CONTINUE
      DPARAM(2) = DH11
      DPARAM(5) = DH22
      GO TO 135
  130 CONTINUE
      DPARAM(2) = DH11
      DPARAM(3) = DH21
      DPARAM(4) = DH12
      DPARAM(5) = DH22
  135 CONTINUE
      DPARAM(1) = DFLAG
      RETURN
      END
