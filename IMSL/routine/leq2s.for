C   IMSL ROUTINE NAME   - LEQ2S
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - LINEAR EQUATION SOLUTION - INDEFINITE MATRIX
C                           - SYMMETRIC STORAGE MODE - HIGH ACCURACY
C                           SOLUTION
C
C   USAGE               - CALL LEQ2S (A,N,B,M,IB,IJOB,ICHNG,DET,IER)
C
C   ARGUMENTS    A      - THE COEFFICIENT MATRIX OF THE EQUATION
C                           AX = B, WHERE A IS ASSUMED TO BE AN N BY N
C                           SYMMETRIC MATRIX. A IS STORED IN SYMMETRIC
C                           STORAGE MODE AND THEREFORE HAS DIMENSION
C                           N*(N+1)/2. (INPUT)
C                N      - ORDER OF A AND THE NUMBER OF ROWS IN B.
C                           (INPUT)
C                B      - INPUT/OUTPUT MATRIX OF DIMENSION N BY M.
C                         ON INPUT, B CONTAINS THE M RIGHT HAND SIDES
C                           OF THE EQUATION AX = B.
C                         ON OUTPUT, THE SOLUTION MATRIX X REPLACES B.
C                           IF IJOB = 1, B IS NOT USED.
C                M      - NUMBER OF RIGHT HAND SIDES (COLUMNS IN B).
C                           (INPUT)
C                IB     - ROW DIMENSION OF MATRIX B EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM. (INPUT)
C                IJOB   - INPUT OPTION PARAMETER. IJOB = I IMPLIES:
C                           I = 0, FACTOR THE MATRIX A AND SOLVE THE
C                             EQUATION AX = B.
C                           I = 1, FACTOR THE MATRIX A. THE FACTORIZED
C                             FORM OF A IS STORED IN THE FIRST
C                             N*(N+1)/2 LOCATIONS OF DET.
C                           I = 2, SOLVE THE EQUATION AX = B. THIS
C                             OPTION IMPLIES THAT MATRIX A HAS ALREADY
C                             BEEN FACTORED BY LEQ2S USING
C                             IJOB = 0 OR 1. IN THIS CASE, THE
C                             INFORMATION CONTAINED IN DET AND ICHNG
C                             MUST HAVE BEEN SAVED FOR REUSE IN THE
C                             CALL TO LEQ2S.
C                ICHNG  - WORK AREA OF LENGTH 2N.
C                DET    - WORK AREA OF LENGTH N*(N+1)/2+3N.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER = 129 INDICATES THAT MATRIX A IS
C                             ALGORITHMICALLY SINGULAR. (SEE THE
C                             CHAPTER L PRELUDE)
C                           IER = 130 INDICATES THAT ITERATIVE
C                             IMPROVEMENT FAILED TO CONVERGE. THE
C                             MATRIX IS TOO ILL-CONDITIONED.
C
C   REQD. IMSL ROUTINES - SINGLE/LEQ1S,UERTST,UGETIO
C                       - DOUBLE/LEQ1S,UERTST,UGETIO,VXADD,VXMUL,
C                           VXSTO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE LEQ2S  (A,N,B,M,IB,IJOB,ICHNG,DET,IER)
C
      INTEGER            ICHNG(1)
      REAL               A(1),DET(1),B(IB,M),ZERO,XNORM,DXNORM
      DOUBLE PRECISION   SUM
      DATA               ZERO/0.0/,ITMAX/50/
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
      NL = (N*(N+1))/2
      NLP1 = NL+1
      NLPN = NL+N
      NLP2N = NLPN+N
      IF (IJOB .EQ. 2) GO TO 10
C                                  SAVE MATRIX A
      DO 5 I = 1,NL
         DET(I) = A(I)
    5 CONTINUE
C                                  FACTOR MATRIX A
      CALL LEQ1S (DET,N,B,M,IB,1,ICHNG,DET(NLP1),IER)
      IF (IER .NE. 0) GO TO 9000
      IF (IJOB .EQ. 1) GO TO 9005
C                                  SAVE THE RIGHT HAND SIDES
   10 DO 55 J = 1,M
         DO 15 I = 1,N
            DET(NLPN+I) = B(I,J)
   15    CONTINUE
C                                  OBTAIN A SOLUTION
         CALL LEQ1S (DET,N,DET(NLPN+1),1,IB,2,ICHNG,DET(NLP1),IER)
C                                  COMPUTE THE NORM OF THE SOLUTION
         XNORM = ZERO
         DO 20 I = 1,N
            XNORM = AMAX1(XNORM, ABS(DET(NLPN+I)))
   20    CONTINUE
         IF (XNORM .EQ. ZERO) GO TO 55
C                                  COMPUTE THE RESIDUALS
         DO 40 ITER = 1,ITMAX
            L = 1
            DO 30 I = 1,N
               K = L
               SUM = B(I,J)
               INC = 1
               DO 25 JJ = 1,N
                  SUM = SUM-DBLE(A(K))*DBLE(DET(NLPN+JJ))
                  IF (JJ .GE. I) INC = JJ
                  K = K+INC
   25          CONTINUE
               DET(NLP2N+I) = SUM
               L = L+I
   30       CONTINUE
            CALL LEQ1S (DET,N,DET(NLP2N+1),1,IB,2,ICHNG,DET(NLP1),IER)
            DXNORM = ZERO
            DO 35 I = 1,N
               DET(NLPN+I) = DET(NLPN+I)+DET(NLP2N+I)
               DXNORM = AMAX1(DXNORM, ABS(DET(NLP2N+I)))
   35       CONTINUE
            DXNORM = DXNORM+XNORM
            IF (DXNORM .EQ. XNORM) GO TO 45
   40    CONTINUE
         IER = 130
C                                  STORE THE SOLUTION
   45    DO 50 I = 1,N
            B(I,J) = DET(NLPN+I)
   50    CONTINUE
         IF (IER .NE. 0) GO TO 9000
   55 CONTINUE
      GO TO 9005
 9000 CONTINUE
      CALL UERTST(IER,6HLEQ2S )
 9005 RETURN
      END
