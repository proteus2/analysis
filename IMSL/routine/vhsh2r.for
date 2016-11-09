C   IMSL ROUTINE NAME   - VHSH2R
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - REAL HOUSEHOLDER TRANSFORMATION TO ZERO A
C                           SINGLE ELEMENT OF A MATRIX
C
C   USAGE               - CALL VHSH2R (AJ,AJP1,UJ,UJP1,VJ,VJP1)
C
C   ARGUMENTS    AJ     - THE J-TH AND (J+1)-TH ENTRIES OF THE COLUMN
C                AJP1       OF A CONTAINING THE ELEMENT TO BE ZEROED.
C                           AJP1 IS THE ELEMENT TO BE ZEROED. (INPUT)
C                UJ     - UJ,UJP1,VJ, AND VJP1 SPECIFY THE HOUSEHOLDER
C                UJP1       TRANSFORMATION P = I + (UJ,UJP1)*(VJ,VJP1)
C                VJ         **T WHICH WILL ZERO AJP1. USING THIS
C                VJP1       OUTPUT, THE USER IS RESPONSIBLE FOR FORMING
C                           P AND TRANSFORMING A. (OUTPUT)
C
C   REQD. IMSL ROUTINES - NONE REQUIRED
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE VHSH2R (AJ,AJP1,UJ,UJP1,VJ,VJP1)
C
      REAL               AJ,AJP1,UJ,UJP1,VJ,VJP1,S,R,ZERO,ONE
      DATA               ZERO,ONE/0.,1./
C                                  FIRST EXECUTABLE STATEMENT
      UJ=ZERO
      UJP1=ZERO
      VJ=ZERO
      VJP1=ZERO
      IF(AJP1.EQ.ZERO) GO TO 9005
      S=ABS(AJ)+ABS(AJP1)
      UJ=AJ/S
      UJP1=AJP1/S
      R=SQRT(UJ*UJ+UJP1*UJP1)
      IF(UJ.LT.ZERO) R=-R
      VJ=-(UJ+R)/R
      VJP1=-UJP1/R
      UJ=ONE
      UJP1=VJP1/VJ
 9005 RETURN
      END
