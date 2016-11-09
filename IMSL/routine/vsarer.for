C   IMSL ROUTINE NAME   - VSARER
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - NOVEMBER 1, 1984
C
C   PURPOSE             - NUCLEUS CALLED BY IMSL ROUTINES WHICH SORT
C                           A MATRIX WITH RESPECT TO VECTORS.
C
C   REQD. IMSL ROUTINES - VNSWAP,VNINI,VSENA
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE VSARER (IA,NR,NC,IOP,KPOS,NK,LR,LC,LK,IR,LIR,LJR,LH,IK,
     *                   ID,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            IA,NR,NC,NK,LR,LC,LK,LIR,LJR,LH,IK,ID,IER,
     *                   IOP(1),KPOS(1),IR(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,IU,JIR,KP,LKH,LKM,LKP,NK1
C                                  FIRST EXECUTABLE STATEMENT
C                                  CHECK VALIDTY OF IOP
      IER = 0
      LR = 1
      LC = 1
      LK = 1
      LIR = 1
      LH = 1
      IK = 1
      ID = 1
      IF (IOP(1).LT.0 .OR. IOP(1).GT.1) IER = 129
      IF (IOP(2).LT.0 .OR. IOP(2).GT.4) IER = 129
      IF (IOP(3).LT.0 .OR. IOP(3).GT.1) IER = 129
      IF (IOP(4).LT.0 .OR. IOP(4).GT.5) IER = 129
C                                  CHECK VALIDITY OF IA,NR,NC,NK
      IF (IA.LE.0 .OR. NR.LE.0 .OR. NC.LE.0 .OR. NK.LE.0) IER = 129
      IF (IER.EQ.129) GO TO 9000
      LC = NC
      LR = NR
      IF (NR.GT.IA) LR = IA
      IF (NR.GT.IA) IER = 65
      JIR = (1-IOP(1))*LR+IOP(1)*LC
      IF (IOP(4).NE.0) LIR = (1-IOP(1))*LC+IOP(1)*LR
      IF (IOP(4).NE.0) LH = (LIR+1)/2
      LKP = (1-IOP(1))*LR+IOP(1)*LC
      LJR = LKP
C                                  START OF INDEX OF DISTINCT RECORDS
      ID = 1
      IF (IOP(4).GT.3) ID = LIR+1
C                                  START OF INDEX OF KEYS PERMUTATION
      IK = 1
      IF (IOP(4).GT.0) IK = LIR+1
C                                  INITIALIZE PERMUTATION VECTOR IR
      IF (IOP(4).NE.0) CALL VNINI(LIR,IR,1,1,1)
      CALL VNINI(LKP,IR(IK),1,1,1)
C                                  CHECK FOR NEGATIVE KEYS IN KPOS
C                                  AND BUILD KEYS PERMUTATION VECTOR
      NK1 = 0
      DO 5 I=1,NK
         KP = KPOS(I)
         IF (KP.LE.0 .OR. KP.GT.JIR) NK1 = NK1+1
         IF (KP.GT.0 .AND. KP.LE.JIR) IR(IK+KP-1) = -I
    5 CONTINUE
      IF (NK1.EQ.NK) IER = 129
      IF (NK1.LT.NK .AND. NK1.GT.0) IER = 65
      IF (IER.EQ.129) GO TO 9000
      LK = NK-NK1
      CALL VSENA(IR(IK),LKP)
      LKH = LK/2
      LKM = (LK+1)/2+1
      CALL VNSWAP(LKH,IR(IK),1,IR(LKM+IK-1),-1)
      DO 10 I=1,LK
         IU = -IR(I+IK-1)
         IR(I+IK-1) = KPOS(IU)
   10 CONTINUE
 9000 RETURN
      END
