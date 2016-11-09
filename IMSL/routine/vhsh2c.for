C   IMSL ROUTINE NAME   - VHSH2C
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - COMPLEX HOUSEHOLDER TRANSFORMATION TO ZERO A
C                           SINGLE ELEMENT OF A MATRIX
C
C   USAGE               - CALL VHSH2C (AJR,AJI,AJP1R,AJP1I,C,SR,SI)
C
C   ARGUMENTS    AJR    - THE REAL AND IMAGINARY PARTS OF THE J-TH
C                AJI        ENTRY OF THE COLUMN OF A CONTAINING THE
C                           ELEMENT TO BE ZEROED. (INPUT)
C                AJP1R  - THE REAL AND IMAGINARY PARTS OF THE (J+1)-TH
C                AJP1I      ENTRY OF THE COLUMN OF A CONTAINING THE
C                           ELEMENT TO BE ZEROED. (INPUT)
C                C      - THE REAL AND IMAGINARY PARTS OF THE HOUSE-
C                SR         HOLDER TRANSFORMATION WHICH WILL ZERO THE
C                SI         (J+1)-TH ENTRY OF A. (OUTPUT)
C
C   REQD. IMSL ROUTINES - NONE REQUIRED
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE VHSH2C (AJR,AJI,AJP1R,AJP1I,C,SR,SI)
C
      REAL               AJR,AJI,AJP1R,AJP1I,C,SR,SI,R,ZERO,ONE
      DATA               ZERO,ONE/0.,1./
C                                  FIRST EXECUTABLE STATEMENT
      IF(AJP1R.EQ.ZERO.AND.AJP1I.EQ.ZERO) GO TO 5
      IF(AJR.EQ.ZERO.AND.AJI.EQ.ZERO) GO TO 10
      R=SQRT(AJR*AJR+AJI*AJI)
      C=R
      SR=(AJR*AJP1R+AJI*AJP1I)/R
      SI=(AJR*AJP1I-AJI*AJP1R)/R
      R=ONE/SQRT(C*C+SR*SR+SI*SI)
      C=C*R
      SR=SR*R
      SI=SI*R
      GO TO 9005
    5 C=ONE
      SR=ZERO
      SI=ZERO
      GO TO 9005
   10 C=ZERO
      SR=ONE
      SI=ZERO
 9005 RETURN
      END
