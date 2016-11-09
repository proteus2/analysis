C   IMSL ROUTINE NAME   - LSVDB
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - NOVEMBER 1, 1984
C
C   PURPOSE             - SINGULAR VALUE DECOMPOSITION OF A BIDIAGONAL
C                           MATRIX.
C
C   USAGE               - CALL LSVDB (D,E,N,V,IV,NRV,C,IC,NCC,IER)
C
C   ARGUMENTS    D      - VECTOR OF LENGTH N. (INPUT/OUTPUT)
C                         ON INPUT, D CONTAINS THE DIAGONAL ELEMENTS
C                           OF THE BIDIAGONAL MATRIX B. D(I)=B(I,I),
C                           I=1,...,N.
C                         ON OUTPUT, D CONTAINS THE N (NONNEGATIVE)
C                           SINGULAR VALUES OF B IN NONINCREASING
C                           ORDER.
C                E      - VECTOR OF LENGTH N. (INPUT/OUTPUT)
C                         ON INPUT, E CONTAINS THE SUPERDIAGONAL
C                           ELEMENTS OF B. E(1) IS ARBITRARY,
C                           E(I)=B(I-1,I), I=2,...,N.
C                         ON OUTPUT, THE CONTENTS OF E ARE MODIFIED
C                           BY THE SUBROUTINE.
C                N      - ORDER OF THE MATRIX B. (INPUT)
C                V      - NRV BY N MATRIX. (INPUT/OUTPUT)
C                           IF NRV.LE.0, V IS NOT USED. OTHERWISE,
C                           V IS REPLACED BY THE NRV BY N PRODUCT
C                           MATRIX V*VB. SEE REMARKS.
C                IV     - ROW DIMENSION OF MATRIX V EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM. (INPUT)
C                NRV    - NUMBER OF ROWS OF V. (INPUT)
C                C      - N BY NCC MATRIX. (INPUT/OUTPUT)
C                           IF NCC.LE.0 C IS NOT USED. OTHERWISE, C
C                           IS REPLACED BY THE N BY NCC PRODUCT
C                           MATRIX UB**(T) * C. SEE REMARKS.
C                IC     - ROW DIMENSION OF MATRIX C EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM. (INPUT)
C                NCC    - NUMBER OF COLUMNS IN C. (INPUT)
C                IER    - ERROR PARAMETER. (INPUT)
C                         WARNING ERROR
C                           IER=33 INDICATES THAT MATRIX B IS NOT FULL
C                             RANK OR VERY ILL-CONDITIONED. SMALL
C                             SINGULAR VALUES MAY NOT BE VERY ACCURATE.
C                         TERMINAL ERROR
C                           IER=129 INDICATES THAT CONVERGENCE WAS
C                             NOT ATTAINED AFTER 10*N QR SWEEPS.
C                             (CONVERGENCE USUALLY OCCURS IN ABOUT
C                             2*N SWEEPS).
C
C   REQD. IMSL ROUTINES - SINGLE/LSVG2,VHS12,UERTST,UGETIO,VBLA=SROTG
C                       - DOUBLE/LSVG2,VHS12,UERTST,UGETIO,VBLA=DROTG
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS      LSVDB COMPUTES THE SINGULAR VALUE DECOMPOSITION OF
C                AN N BY N BIDIAGONAL MATRIX
C                     B = UB * S * VB**(T)    WHERE
C                UB AND VB ARE N BY N ORTHOGONAL MATRICES AND
C                S IS DIAGONAL.
C                IF ARGUMENTS V AND C ARE N BY N IDENTITY MATRICES,
C                ON EXIT THEY ARE REPLACED BY VB AND UB**T,
C                RESPECTIVELY.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE LSVDB (D,E,N,V,IV,NRV,C,IC,NCC,IER)
C
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N,IV,NRV,IC,NCC,IER
      REAL               D(N),E(N),V(IV,1),C(IC,1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,II,J,K,KK,L,LL,LP1,NQRS,N10
      LOGICAL            WNTV,HAVERS,FAIL
      REAL               DNORM,ZERO,ONE,TWO,CS,F,SQINF,FTEMP,G,H,HTEMP,
     *                   SN,T,X,Y,Z
      DATA               SQINF/0.8507057E38/
      DATA               ZERO/0.0/,ONE/1.0/,TWO/2.0/
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
      IF (N.LE.0) GO TO 9005
      N10 = 10*N
      WNTV = NRV.GT.0
      HAVERS = NCC.GT.0
      FAIL = .FALSE.
      NQRS = 0
      E(1) = ZERO
      DNORM = ZERO
      DO 5 J=1,N
    5 DNORM = AMAX1(ABS(D(J))+ABS(E(J)),DNORM)
      DO 100 KK=1,N
         K = N+1-KK
C                                  TEST FOR SPLITTING OR RANK
C                                    DEFICIENCIES FIRST MAKE TEST FOR
C                                    LAST DIAGONAL TERM, D(K), BEING
C                                    SMALL.
   10    IF (K.EQ.1) GO TO 25
         T = DNORM+D(K)
         IF (T.NE.DNORM) GO TO 25
C
C                                  SINCE D(K) IS SMALL WE WILL MAKE A
C                                    SPECIAL PASS TO TRANSFORM E(K) TO
C                                    ZERO.
         CS = ZERO
         SN = -ONE
         DO 20 II=2,K
            I = K+1-II
            F = -SN*E(I+1)
            E(I+1) = CS*E(I+1)
            T = D(I)
            FTEMP = F
            CALL SROTG (D(I),FTEMP,CS,SN)
C                                  TRANSFORMATION CONSTRUCTED TO ZERO
C                                    POSITION (I,K).
            IF (.NOT.WNTV) GO TO 20
            DO 15 J=1,NRV
   15       CALL LSVG2 (CS,SN,V(J,I),V(J,K))
C
C                                  ACCUMULATE RT. TRANSFORMATIONS IN V.
   20    CONTINUE
C                                  THE MATRIX IS NOW BIDIAGONAL, AND OF
C                                    LOWER ORDER SINCE E(K) .EQ. ZERO
   25    DO 30 LL=1,K
            L = K+1-LL
            T = DNORM+E(L)
            IF (T.EQ.DNORM) GO TO 50
            T = DNORM+D(L-1)
            IF (T.EQ.DNORM) GO TO 35
   30    CONTINUE
C                                  THIS LOOP CANT COMPLETE SINCE E(1) =
C                                    ZERO.
         GO TO 50
C                                  CANCELLATION OF E(L), L.GT.1.
   35    CS = ZERO
         SN = -ONE
         DO 45 I=L,K
            F = -SN*E(I)
            E(I) = CS*E(I)
            T = DNORM+F
            IF (T.EQ.DNORM) GO TO 50
            T = D(I)
            FTEMP = F
            CALL SROTG (D(I),FTEMP,CS,SN)
            IF (.NOT.HAVERS) GO TO 45
            DO 40 J=1,NCC
   40       CALL LSVG2 (CS,SN,C(I,J),C(L-1,J))
   45    CONTINUE
C                                  TEST FOR CONVERGENCE
   50    Z = D(K)
         IF (L.EQ.K) GO TO 85
C                                  SHIFT FROM BOTTOM 2 BY 2 MINOR OF
C                                    B**(T)*B.
         X = D(L)
         Y = D(K-1)
         G = E(K-1)
         H = E(K)
         F = ((Y-Z)*(Y+Z)+(G-H)*(G+H))/(TWO*H*Y)
         G = ABS(F)
         IF (ABS(F) .LT. SQINF) G = SQRT(ONE+F**2)
         IF (F.LT.ZERO) GO TO 55
         T = F+G
         GO TO 60
   55    T = F-G
   60    F = ((X-Z)*(X+Z)+H*(Y/T-H))/X
C                                  NEXT QR SWEEP
         CS = ONE
         SN = ONE
         LP1 = L+1
         DO 80 I=LP1,K
            G = E(I)
            Y = D(I)
            H = SN*G
            G = CS*G
            HTEMP = H
            CALL SROTG (F,HTEMP,CS,SN)
            E(I-1) = F
            F = X*CS+G*SN
            G = -X*SN+G*CS
            H = Y*SN
            Y = Y*CS
            IF (.NOT.WNTV) GO TO 70
C                                  ACCUMULATE ROTATIONS (FROM THE
C                                    RIGHT) IN V
            DO 65 J=1,NRV
   65       CALL LSVG2 (CS,SN,V(J,I-1),V(J,I))
            HTEMP = H
   70       CALL SROTG (F,HTEMP,CS,SN)
            D(I-1) = F
            F = CS*G+SN*Y
            X = -SN*G+CS*Y
            IF (.NOT.HAVERS) GO TO 80
            DO 75 J=1,NCC
   75       CALL LSVG2 (CS,SN,C(I-1,J),C(I,J))
C
C                                  APPLY ROTATIONS FROM THE LEFT TO
C                                    RIGHT HAND SIDES IN C
   80    CONTINUE
         E(L) = ZERO
         E(K) = F
         D(K) = X
         NQRS = NQRS+1
         IF (NQRS.LE.N10) GO TO 10
C                                  RETURN TO TEST FOR SPLITTING.
         FAIL = .TRUE.
C                                  CUTOFF FOR CONVERGENCE FAILURE. NQRS
C                                    WILL BE 2*N USUALLY.
   85    IF (Z.GE.ZERO) GO TO 95
         D(K) = -Z
         IF (.NOT.WNTV) GO TO 95
         DO 90 J=1,NRV
   90    V(J,K) = -V(J,K)
   95    CONTINUE
C                                  CONVERGENCE. D(K) IS MADE
C                                    NONNEGATIVE
  100 CONTINUE
      IF (N.EQ.1) GO TO 140
      DO 105 I=2,N
         IF (D(I).GT.D(I-1)) GO TO 110
  105 CONTINUE
      GO TO 140
C                                  EVERY SINGULAR VALUE IS IN ORDER
  110 DO 135 I=2,N
         T = D(I-1)
         K = I-1
         DO 115 J=I,N
            IF (T.GE.D(J)) GO TO 115
            T = D(J)
            K = J
  115    CONTINUE
         IF (K.EQ.I-1) GO TO 135
         D(K) = D(I-1)
         D(I-1) = T
         IF (.NOT.HAVERS) GO TO 125
         DO 120 J=1,NCC
            T = C(I-1,J)
            C(I-1,J) = C(K,J)
  120    C(K,J) = T
  125    IF (.NOT.WNTV) GO TO 135
         DO 130 J=1,NRV
            T = V(J,I-1)
            V(J,I-1) = V(J,K)
  130    V(J,K) = T
  135 CONTINUE
C                                  END OF ORDERING ALGORITHM.
  140 IER = 129
      IF (FAIL) GO TO 9000
C                                  CHECK FOR POSSIBLE RANK DEFICIENCY
      IER = 33
      T = 0.0
      IF (D(1).NE.ZERO) T=D(N)/D(1)
      F=100.0+T
      IF (F.EQ.100.0) GO TO 9000
      IER = 0
      GO TO 9005
 9000 CONTINUE
      CALL UERTST (IER,6HLSVDB )
 9005 RETURN
      END
