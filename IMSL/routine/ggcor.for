C   IMSL ROUTINE NAME   - GGCOR
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - NOVEMBER 1, 1984
C
C   PURPOSE             - GENERATE A RANDOM ORTHOGONAL MATRIX AND A
C                           RANDOM CORRELATION MATRIX
C
C   USAGE               - CALL GGCOR (DSEED,N,E,A,IA,COR,IWK,WK,IER)
C
C   ARGUMENTS    DSEED  - INPUT/OUTPUT DOUBLE PRECISION VARIABLE
C                           ASSIGNED AN INTEGER VALUE IN THE
C                           EXCLUSIVE RANGE (1.D0, 2147483647.D0).
C                           DSEED IS REPLACED BY A NEW VALUE TO BE
C                           USED IN A SUBSEQUENT CALL.
C                N      - INPUT ORDER OF MATRICES TO BE GENERATED.
C                E      - INPUT VECTOR OF LENGTH N CONTAINING THE
C                           EIGENVALUES OF THE CORRELATION MATRIX
C                           TO BE GENERATED.  THE ELEMENTS OF E MUST
C                           BE POSITIVE AND MUST SUM TO N.
C                A      - OUTPUT N BY N ORTHOGONAL MATRIX IN FULL
C                           STORAGE MODE.
C                IA     - INPUT ROW DIMENSION OF MATRIX A EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN
C                           THE CALLING PROGRAM.
C                COR    - OUTPUT N BY N CORRELATION MATRIX IN SYMMETRIC
C                           STORAGE MODE.
C                IWK    - WORK VECTOR OF LENGTH N.
C                WK     - WORK VECTOR OF LENGTH 2*N+2.
C                IER    - ERROR PARAMETER.  (OUTPUT)
C                           TERMINAL ERROR
C                             IER = 129 INDICATES THAT THE EIGENVALUES
C                               IN E ARE INCORRECT (NONPOSITIVE OR DO
C                               NOT SUM TO N).
C                             IER = 130 INDICATES THAT N IS LESS THAN 2.
C                           WARNING WITH FIX ERROR
C                             IER = 67 INDICATES THAT AN ERROR OCCURRED
C                               IN THE ROTATIONS USED TO FORM THE
C                               CORRELATION MATRIX.  THE EIGENVALUES
C                               OF COR MAY DIFFER FROM THE VALUES
C                               SPECIFIED IN E BY MORE THAN JUST
C                               ROUNDING DIFFERENCES.
C
C   REQD. IMSL ROUTINES - GGCOT,GGNPM,GGUBFS,UERTST,UGETIO,VBLA=SNRM2,
C                           VSRTR
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE GGCOR (DSEED,N,E,A,IA,COR,IWK,WK,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N,IA,IWK(1),IER
      REAL               E(1),A(IA,1),COR(1),WK(1)
      DOUBLE PRECISION   DSEED
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,I1,I1COL,I1M1,I1ROW,I2,I2COL,I2P1,I2ROW,IAL,
     *                   IBET,IBETA,IERQ,IGAM,ILARGE,IP1,ISMALL,J,K,L,
     *                   LP1,LPN,NM1,NMI,NP1,NP1L,NP2L,NTRI
      REAL               AA,ADJ,ALNEW,ALPHA,BB,BBS,BETA,C,CI1,CI2,
     *                   COEFF0,COEFF1,COEFF2,CSQ,CSQ1,CSQ2,DD,ESUM,GAM,
     *                   HV,HW,S,S1,SSQ,TOL,UL,WIL,WL,WWIL,XEPS,XN
      REAL               GGUBFS
      REAL               SNRM2
      DOUBLE PRECISION   TEMP
      COMPLEX            ZROOT1,ZROOT2
      EQUIVALENCE        (ZROOT1,CSQ1),(ZROOT2,CSQ2)
      DATA               XEPS /Z3C100000/
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
      IF (N.GE.2) GO TO 5
      IER = 130
      GO TO 9000
    5 ESUM = E(1)
      DO 10 I=2,N
         ESUM = ESUM+E(I)
   10 CONTINUE
      TOL = 5.0*N*N*XEPS
      XN = N
      IF (ESUM.LT.XN+TOL .AND. ESUM.GT.XN-TOL) GO TO 15
      IER = 129
      GO TO 9000
   15 NP1 = N+1
C                                  GENERATE A RANDOM ORTHOGONAL MATRIX
C                                  FROM A HAAR INVARIANT DISTRIBUTION.
      DO 25 I=1,N
         DO 20 J=1,N
            A(I,J) = 0.0
   20    CONTINUE
         A(I,I) = 1.0
   25 CONTINUE
      DO 65 L=1,N
         NP1L = NP1-L
         LPN = L+N
         NP2L = NP1L+1
   30    CALL GGNPM(DSEED,NP2L,WK(LPN))
         UL = SNRM2(NP1L,WK(LPN),1)
         IF (UL.LE.0.0) GO TO 30
         IF (L.EQ.N) GO TO 55
C                                  THIS ALGORITHM IS NOT CODED EXACTLY
C                                  AS HEIBERGER (JRSSC 27/2 (199-206)).
C                                  SO THE FOLLOWING FIX SUGGESTED BY
C                                  TANNER AND THISTED (JRSSC 31/2
C                                  (190-192)) IS NOT IMPLEMENTED.
C                                  CHANGE THE NEXT LINE AS FOLLOWS
C                                  FROM   WL = SIGN(UL,-WK(LPN))
C                                  TO     WL = SIGN(UL, WK(LPN)).
         WL = SIGN(UL,-WK(LPN))
         DO 35 J=L,N
            WK(J) = WK(J+N)
   35    CONTINUE
         LP1 = L+1
         WIL = 1.0/WL
         WWIL = 1.0/(WL-WK(L))
         DO 50 I=1,N
            TEMP = 0.0D0
            DO 40 J=L,N
               TEMP = TEMP+DBLE(A(I,J))*DBLE(WK(J))
   40       CONTINUE
            HV = TEMP*WIL
            HW = (HV-A(I,L))*WWIL
            A(I,L) = HV
            DO 45 J=LP1,N
               A(I,J) = A(I,J)-HW*WK(J)
   45       CONTINUE
   50    CONTINUE
   55    IF (WK(N+L).GT.0.0) GO TO 65
         DO 60 I=1,N
            A(I,L) = -A(I,L)
   60    CONTINUE
   65 CONTINUE
C                                  FORM A(TRANSPOSE)*E*A
   70 K = 0
      DO 85 I=1,N
         DO 80 J=1,I
            K = K+1
            TEMP = 0.0D0
            DO 75 L=1,N
               TEMP = DBLE(A(L,I))*DBLE(E(L))*DBLE(A(L,J))+TEMP
   75       CONTINUE
            COR(K) = TEMP
   80    CONTINUE
   85 CONTINUE
C                                  STORE DIAGONAL ELEMENTS AND SORT
      L = 0
      DO 90 I=1,N
         L = L+I
         IWK(I) = I
         WK(I) = COR(L)
   90 CONTINUE
      CALL VSRTR(WK,N,IWK)
C                                  PERFORM PLANAR ROTATIONS BEGINNING
C                                  WITH SMALLEST AND LARGEST DIAGONALS,
C                                  SO THAT SMALLER DIAGONAL BECOMES 1.0.
C                                  THEN REPEAT, USING THE SMALLEST AND
C                                  LARGEST OF THE DIAGONALS THAT REMAIN.
      NM1 = N-1
      NTRI = (N*(N+1))/2
      DO 135 I=1,NM1
         ISMALL = IWK(I)
         IAL = (ISMALL*(ISMALL+1))/2
         ALPHA = COR(IAL)
         ILARGE = IWK(N)
         IGAM = (ILARGE*(ILARGE+1))/2
         GAM = COR(IGAM)
         IF (ALPHA.GE.1.0) GO TO 140
         IF (GAM.LE.1.0) GO TO 140
C                                  I1 IS THE FIRST ROW TO BE
C                                  TRANSFORMED AND I2 IS THE
C                                  SECOND.
         I1 = MIN0(ISMALL,ILARGE)
         I2 = MAX0(ISMALL,ILARGE)
         IBETA = (I2*(I2-1))/2+I1
         BETA = COR(IBETA)
C                                  DETERMINE SINE AND COSINE TO
C                                  ROTATE I1(TH) COLUMN INTO A
C                                  VECTOR WITH A ONE ON THE
C                                  DIAGONAL.
         AA = 1.0-GAM
         BB = ALPHA-GAM
         BBS = BB*BB
         DD = 2.0*BETA
         COEFF0 = AA*AA
         COEFF1 = -(2.0*AA*BB+DD*DD)
         COEFF2 = BB*BB+DD*DD
         CALL GGCOT(COEFF2,COEFF1,COEFF0,ZROOT1,ZROOT2,IERQ)
         IF (CSQ1.GE.0.0 .AND. CSQ1.LE.1.0) GO TO 95
         IF (CSQ1.LT.0.0 .AND. CSQ2-1.0.LT.-CSQ1) CSQ = CSQ2
         IF (CSQ.GE.0.0 .AND. CSQ.LE.1.0) GO TO 100
         IER = 67
         IF (CSQ.LT.0.0) CSQ = 0.0
         IF (CSQ.GT.1.0) CSQ = 1.0
         GO TO 100
   95    CSQ = CSQ1
  100    SSQ = 1.0-CSQ
         C = SQRT(CSQ)
         S = SQRT(SSQ)
         IF (GGUBFS(DSEED).GT.0.5) C = -C
         ALNEW = CSQ*ALPHA+SSQ*GAM
         ADJ = 2.*S*C*BETA
         IF (ABS(ALNEW-ADJ-1.0).LT.ABS(ALNEW+ADJ-1.0)) S = -S
         S1 = S
         IF (I1.EQ.ILARGE) S1 = -S1
         IF (N.EQ.2) GO TO 130
C                                  FIX ELEMENTS OF COR IN I1(TH)
C                                  AND I2(TH) ROWS AND COLUMNS.
         IF (I1.EQ.1) GO TO 110
         I1ROW = (I1*(I1-1))/2+1
         I2ROW = (I2*(I2-1))/2+1
         I1M1 = I1-1
         DO 105 J=1,I1M1
            CI1 = COR(I1ROW)
            CI2 = COR(I2ROW)
            COR(I1ROW) = C*CI1+S1*CI2
            COR(I2ROW) = -S1*CI1+C*CI2
            I1ROW = I1ROW+1
            I2ROW = I2ROW+1
  105    CONTINUE
  110    IF (I2.EQ.N) GO TO 120
         I2COL = ((I2+1)*(I2+2))/2-1
         I1COL = I2COL-I2+I1
         I2P1 = I2+1
         DO 115 J=I2P1,N
            CI1 = COR(I1COL)
            CI2 = COR(I2COL)
            COR(I1COL) = C*CI1+S1*CI2
            COR(I2COL) = -S1*CI1+C*CI2
            I1COL = I1COL+J
            I2COL = I2COL+J
  115    CONTINUE
  120    IBET = I2-I1-1
         IF (IBET.EQ.0) GO TO 130
         I1COL = ((I1+1)*(I1+2))/2-1
         I2ROW = (I2*(I2+1))/2-I2+I1+1
         DO 125 J=1,IBET
            CI1 = COR(I1COL)
            CI2 = COR(I2ROW)
            COR(I1COL) = C*CI1+S1*CI2
            COR(I2ROW) = -S1*CI1+C*CI2
            I1COL = I1COL+J+I1
            I2ROW = I2ROW+1
  125    CONTINUE
  130    COR(IAL) = CSQ*ALPHA+SSQ*GAM+2.*S*C*BETA
         COR(IGAM) = CSQ*GAM+SSQ*ALPHA-2.*S*C*BETA
         COR(IBETA) = C*S*(GAM-ALPHA)+BETA*(CSQ-SSQ)
         IP1 = I+1
         NMI = N-I
         WK(ILARGE) = COR(IGAM)
         CALL VSRTR(WK(IP1),NMI,IWK(IP1))
  135 CONTINUE
C                                  FIX DIAGONALS.
  140 L = 0
      DO 145 I=1,N
         L = L+I
         COR(L) = 1.0
  145 CONTINUE
      IF (IER.EQ.0) GO TO 9005
 9000 CONTINUE
      CALL UERTST(IER,'GGCOR ')
 9005 RETURN
      END
