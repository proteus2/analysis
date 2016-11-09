C   IMSL ROUTINE NAME   - EHOUSH
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1982
C
C   PURPOSE             - NUCLEUS CALLED ONLY BY IMSL ROUTINE EIGCH
C
C   REQD. IMSL ROUTINES - NONE REQUIRED
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE EHOUSH (AR,AI,N,D,E,TAU)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N
      REAL               AR(1),AI(1),D(1),E(1),TAU(2,1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            NM1,NN,I,NR,NRM1,L,INDX,J,JJ,INX1,INX2,JP1,KK,
     *                   IX,IM1
      REAL               RHO,TOLER,ZERO,ONE,T1,T2,TESTBB,VR,ROOT,DELTA,
     *                   RATIO,RDELP,Q1,Q2,X1,X2,TT1,TT2,BB
      DATA               ZERO/0.0/,ONE/1.0/
      DATA               RDELP/Z3C100000/
C                                  FIRST EXECUTABLE STATEMENT
      NM1=N-1
      TOLER=ZERO
      NN=(N*(N+1))/2
      DO 5 I=1,NN
         T1=ABS(AR(I))
         T2=ABS(AI(I))
         IF(T2.GT.T1) T1=T2
         IF (T1.GT.TOLER) TOLER=T1
    5 CONTINUE
      TESTBB=RDELP*TOLER
      IF (N.LE.2) GO TO 65
C                                  PERFORM N - 2 SIMILARITY
C                                    TRANSFORMATIONS
      DO 60 NR=2,NM1
         NRM1=NR-1
         VR=ZERO
         TAU(1,NR)=ZERO
         TAU(2,NR)=ZERO
         TAU(2,1)=ZERO
         DO 10 L=NR,N
            INDX=(L*(L-1))/2+NRM1
            VR=AR(INDX)**2+AI(INDX)**2+VR
   10    CONTINUE
         INDX=(NR*NRM1)/2+NRM1
         IF ((TESTBB)**2 .GE. VR) GO TO 60
         ROOT = CABS(CMPLX(AR(INDX),AI(INDX)))*SQRT(VR)
         IF(ROOT.NE.ZERO) GO TO 15
         AR(INDX)=SQRT(VR)
         DELTA=VR
         TAU(1,1)=-AR(INDX)
         GO TO 20
   15    DELTA=VR+ROOT
         RATIO=VR/ROOT
         TAU(1,1)=-RATIO*AR(INDX)
         TAU(2,1)= RATIO*AI(INDX)
         AR(INDX)=(RATIO+ONE)*AR(INDX)
         AI(INDX)=(RATIO+ONE)*AI(INDX)
C                                  THE MATRIX TO BE USED IN THE
C                                    SIMILARITY TRANSFORMATION HAS
C                                    BEEN DETERMINED. THE TRANSFOR-
C                                    MATION FOLLOWS
   20    DO 35 J=NR,N
            JJ=(J*(J-1))/2
            INDX=JJ+NRM1
            TAU(1,J)=AR(INDX)/DELTA
            TAU(2,J)=AI(INDX)/DELTA
            D(J)=ZERO
            E(J)=ZERO
            DO 25 L=NR,J
               INX1=(L*(L-1))/2+NRM1
               INX2=JJ+L
               D(J)= D(J)+AR(INX2)*AR(INX1)-AI(INX2)*AI(INX1)
               E(J)= E(J)+AR(INX2)*AI(INX1)+AI(INX2)*AR(INX1)
   25       CONTINUE
            JP1=J+1
            IF (JP1 .GT. N) GO TO 40
            DO 30 L=JP1,N
               KK=(L*(L-1))/2
               INX1=KK+NRM1
               INX2=KK+J
               D(J)=D(J)+AR(INX2)*AR(INX1)+AI(INX2)*AI(INX1)
               E(J)=E(J)+AR(INX2)*AI(INX1)-AI(INX2)*AR(INX1)
   30       CONTINUE
   35    CONTINUE
   40    RHO=ZERO
         DO 45 L=NR,N
            RHO=RHO+D(L)*TAU(1,L)+E(L)*TAU(2,L)
   45    CONTINUE
         IX=(NRM1*(NR-2))/2
         DO 55 I=NR,N
            IX=IX+I-1
            INX2=IX+NRM1
            DO 50 J=NR,I
               INX1=IX+J
               X1=TAU(1,I)*D(J)+TAU(2,I)*E(J)
               X2=TAU(2,I)*D(J)-TAU(1,I)*E(J)
               Q1=D(I)-RHO*AR(INX2)
               Q2=E(I)-RHO*AI(INX2)
               T1=Q1*TAU(1,J)+Q2*TAU(2,J)
               T2=Q2*TAU(1,J)-Q1*TAU(2,J)
               AR(INX1)=AR(INX1)-X1-T1
               AI(INX1)=AI(INX1)-X2-T2
   50       CONTINUE
   55    CONTINUE
         TAU(1,NR)=TAU(1,1)
         TAU(2,NR)=TAU(2,1)
   60 CONTINUE
C                                  THE MATRIX HAS BEEN REDUCED TO TRI-
C                                    DIAGONAL HERMITIAN FORM. THE SUB-
C                                    DIAGONAL HAS BEEN TEMPORARILY
C                                    STORED IN VECTOR TAU. STORE THE
C                                    DIAGONAL OF THE REDUCED MATRIX IN D
   65 INDX=0
      DO 70 I=1,N
         INDX=INDX+I
         D(I)=AR(INDX)
   70 CONTINUE
C                                  PERFORM THE DIAGONAL UNITARY SIMILA-
C                                    RITY TRANSFORMATION
      TAU(1,1)=ONE
      TAU(2,1)=ZERO
      E(1)=ZERO
      IF (N .EQ. 1) GO TO 85
      INDX=(N*NM1)/2+NM1
      TAU(1,N)=AR(INDX)
      TAU(2,N)=-AI(INDX)
C                                  CALCULATE SUBDIAGONAL E OF THE REAL
C                                    SYMMETRIC TRIDIAGONAL MATRIX. CAL-
C                                    CULATE TAU, THE DIAGONAL OF THE
C                                    DIAGONAL UNITARY MATRIX
      INDX=1
      DO 80 I=2,N
         INDX=INDX+I
         IM1=I-1
         BB= SQRT(TAU(1,I)**2+TAU(2,I)**2)
         E(I)=BB
         AI(INDX)=BB
         IF (TESTBB .LT. BB) GO TO 75
         TAU(1,I)=ONE
         TAU(2,I)=ZERO
         BB=ONE
   75    TT1=TAU(1,I)*TAU(1,IM1)-TAU(2,I)*TAU(2,IM1)
         TT2=TAU(1,I)*TAU(2,IM1)+TAU(2,I)*TAU(1,IM1)
         TAU(1,I)=TT1/BB
         TAU(2,I)=TT2/BB
   80 CONTINUE
   85 RETURN
      END
