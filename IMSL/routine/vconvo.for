C   IMSL ROUTINE NAME   - VCONVO
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - VECTOR CONVOLUTION
C
C   USAGE               - CALL VCONVO (A,B,LA,LB,IWK)
C
C   ARGUMENTS    A      - ON INPUT, A MUST CONTAIN (IN ITS FIRST LA
C                           ELEMENTS) THE FIRST DATA VECTOR.
C                         ON OUTPUT, A WILL CONTAIN (IN THE FIRST
C                           LA+LB-1 ELEMENTS)  THE RESULT OF THE
C                           CONVOLUTION. THE LENGTH OF VECTOR A MUST BE
C                           AT LEAST 2**(M+1), WHERE M IS THE SMALLEST
C                           INTEGER SUCH THAT LA+LB-1 IS LESS THAN OR
C                           EQUAL TO 2**M.
C                         NOTE - THE ROUTINE TREATS A AND B AS COMPLEX
C                           VECTORS OF LENGTH 2**M AND 2**(M-1)
C                           RESPECTIVELY. AN APPROPRIATE EQUIVALENCE
C                           STATEMENT MAY BE REQUIRED. SEE DOCUMENT
C                           EXAMPLE.
C                B      - B CONTAINS (IN ITS FIRST LB ELEMENTS) THE
C                           DATA (FILTER) VECTOR.
C                         ON OUTPUT, B IS DESTROYED. THE LENGTH OF
C                           VECTOR B MUST BE AT LEAST 2**M, WHERE M IS
C                           DEFINED ABOVE.
C                LA     - LENGTH OF THE FIRST DATA SEQUENCE. (INPUT)
C                           (SEE DESCRIPTION OF VECTOR A, ABOVE.)
C                LB     - LENGTH OF THE SECOND (FILTER) DATA SEQUENCE.
C                           (INPUT) (SEE DESCRIPTION OF VECTOR B.)
C                IWK    - WORK AREA OF LENGTH M + 1 WHERE M IS
C                           DEFINED ABOVE.
C
C   REQD. IMSL ROUTINES - FFT2C
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE VCONVO (A,B,LA,LB,IWK)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            LA,LB,IWK(1)
      COMPLEX            A(1),B(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            IMAX,NN,NREM,MTWO,M,I,MM,LAP1,LBP1,ND2,
     *                   LAP1S2,LBP1S2,ND4,NDT,NP2,NMK,N2,J
      REAL               PIE,ZERO,ONE,THETA,TP,TEMP,DN,A1(2),GA(2),
     *                   GB(2),HALF
      COMPLEX            XIMAG,ALPH,BETA,GAMA,GAMB,CZERO,C1,S1
      EQUIVALENCE        (GA(1),GAMA),(GB(1),GAMB),(ALPH,A1(1)),
     1                   (ZERO,CZERO)
      DATA               CZERO/(0.,0.)/,IMAX/24/,ONE/1.0/,HALF/0.5/
      DATA               PIE/3.141593/
C                                  FIRST EXECUTABLE STATEMENT
      NN = LA+LB-1
      NREM = (NN+1)/2
C                                  DETERMINE THE SMALLEST INTEGER M
C                                  SUCH THAT LA + LB - 1 IS LESS THAN
C                                  OR EQUAL TO 2**M.
      MTWO = 2
      M = 1
      DO 5 I=1,IMAX
         IF (NN.LE.MTWO) GO TO 10
         MTWO = MTWO+MTWO
         M = M+1
    5 CONTINUE
   10 NN = 2**M
      MM = M-1
C                                  ZERO THE TRAILING ELEMENTS OF BOTH A
C                                  AND B.
      LAP1 = LA+1
      LBP1 = LB+1
      DN = ONE/NN
      ND2 = NN/2
      LAP1S2 = LAP1/2
      LBP1S2 = LBP1/2
      IF (LAP1S2*2.NE.LAP1) GO TO 15
      TEMP = A(LAP1S2)
      A(LAP1S2) = TEMP
   15 LAP1S2 = LAP1S2+1
      DO 20 I=LAP1S2,ND2
         A(I) = CZERO
   20 CONTINUE
      IF (LBP1S2*2.NE.LBP1) GO TO 25
      TEMP = B(LBP1S2)
      B(LBP1S2) = TEMP
   25 LBP1S2 = LBP1S2+1
      DO 30 I=LBP1S2,ND2
         B(I) = CZERO
   30 CONTINUE
      ND4 = ND2/2
      NDT = ND4+1
      THETA = PIE/ND2
      XIMAG = CMPLX(ZERO,ONE)
      NP2 = ND2+2
C                                  COMPUTE THE CENTER ELEMENTS FOR
C                                  BOTH THE A AND B VECTORS.
      GAMA = CZERO
      GAMB = CZERO
      DO 35 I=1,ND2
         GAMA = GAMA+A(I)
         GAMB = GAMB+B(I)
   35 CONTINUE
      TEMP = GA(1)-GA(2)
      GAMA = TEMP
      TEMP = GB(1)-GB(2)
      GAMB = TEMP
C                                  COMPUTE THE FAST FOURIER TRANSFORM
C                                  OF BOTH A AND B AS COMPLEX VECTORS.
      CALL FFT2C (A,MM,IWK)
      CALL FFT2C (B,MM,IWK)
      ALPH = A(1)
      A(1) = A1(1)+A1(2)
      ALPH = B(1)
      B(1) = A1(1)+A1(2)
      TP = THETA
C                                  COMPUTE THE FAST FOURIER TRANSFORMS
C                                  OF BOTH A AND B AS REAL VECTORS
C                                  USING THE COMPLEX RESULT ABOVE.
      DO 40 K=2,ND4
         NMK = NP2-K
         S1 = CONJG(A(NMK))
         ALPH = A(K)+S1
         BETA = XIMAG*(S1-A(K))
         S1 = CMPLX(COS(THETA),SIN(THETA))
         A(K) = (ALPH+BETA*S1)*HALF
         A(NMK) = CONJG(ALPH-BETA*S1)*HALF
         C1 = CONJG(B(NMK))
         ALPH = B(K)+C1
         BETA = XIMAG*(C1-B(K))
         B(K) = (ALPH+BETA*S1)*HALF
         B(NMK) = CONJG(ALPH-BETA*S1)*HALF
         THETA = THETA+TP
   40 CONTINUE
      N2 = NN+2
C                                  MULTIPLY TRANSFORMS OF A AND B AND
C                                  TAKE THE INVERSE TRANSFORM.
      DO 45 I=1,ND2
         A(I) = CONJG(A(I)*B(I))
   45 CONTINUE
      DO 50 I=2,ND2
         J = N2-I
         A(J) = CONJG(A(I))
   50 CONTINUE
      A(ND2+1) = GAMA*GAMB
      CALL FFT2C (A,M,IWK)
C                                  REDUCE THE COMPLEX RESULT TO ITS
C                                  REAL COMPONENTS ONLY.
      J = -1
      DO 55 I=1,NREM
         J = J+2
         A1(1) = A(J)
         A1(2) = A(J+1)
         A(I) = ALPH*DN
   55 CONTINUE
      RETURN
      END
