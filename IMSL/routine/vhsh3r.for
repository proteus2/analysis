C   IMSL ROUTINE NAME   - VHSH3R
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - REAL HOUSEHOLDER TRANSFORMATION TO ZERO TWO
C                           ELEMENTS OF A MATRIX
C
C   USAGE               - CALL VHSH3R (AJ,AJP1,AJP2,UJ,UJP1,UJP2,VJ,
C                           VJP1,VJP2)
C
C   ARGUMENTS    AJ     - THE J-TH, (J+1)-TH, AND (J+2)-TH ENTRIES OF
C                AJP1       THE COLUMN OF A CONTAINING THE ELEMENTS TO
C                AJP2       BE ZEROED. AJP1 AND AJP2 ARE THE ENTRIES TO
C                           BE ZEROED. (INPUT)
C                UJ     - UJ,UJP1,UJP2,VJ,VJP1, AND VJP2 SPECIFY THE
C                UJP1       HOUSEHOLDER TRANSFORMATION
C                UJP2       P = I+(UJ,UJP1,UJP2)*(VJ,VJP1,VJP2)**T
C                VJ         WHICH WILL ZERO AJP1 AND AJP2. USING THIS
C                VJP1       OUTPUT, THE USER IS RESPONSIBLE FOR FORMING
C                VJP2       P AND TRANSFORMING A. (OUTPUT)
C
C
C   REQD. IMSL ROUTINES - NONE REQUIRED
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE VHSH3R (AJ,AJP1,AJP2,UJ,UJP1,UJP2,VJ,VJP1,VJP2)
C
      REAL               AJ,AJP1,AJP2,UJ,UJP1,UJP2,VJ,VJP1,VJP2,ZERO
      REAL               ONE,S,R,RD
      DATA               ZERO,ONE/0.,1./
C                                  FIRST EXECUTABLE STATEMENT
      UJ=ZERO
      UJP1=ZERO
      UJP2=ZERO
      VJ=ZERO
      VJP1=ZERO
      VJP2=ZERO
      IF(AJP1.EQ.ZERO.AND.AJP2.EQ.ZERO) GO TO 9005
      S=ABS(AJ)+ABS(AJP1)+ABS(AJP2)
      RD=ONE/S
      UJ=AJ*RD
      UJP1=AJP1*RD
      UJP2=AJP2*RD
      R=SQRT(UJ*UJ+UJP1*UJP1+UJP2*UJP2)
      IF(UJ.LT.ZERO) R=-R
      RD=ONE/R
      VJ=-(UJ+R)*RD
      VJP1=-UJP1*RD
      VJP2=-UJP2*RD
      UJ=ONE
      UJP1=VJP1/VJ
      UJP2=VJP2/VJ
 9005 RETURN
      END
