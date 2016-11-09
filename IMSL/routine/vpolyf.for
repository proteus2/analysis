C   IMSL ROUTINE NAME   - VPOLYF
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - MATRIX POLYNOMIAL (FULL STORAGE MODE)
C
C   USAGE               - CALL VPOLYF (A,N,IA,COEF,M,B,IB,WKAREA)
C
C   ARGUMENTS    A      - INPUT N BY N MATRIX FOR WHICH THE
C                           POLYNOMIAL IS TO BE EVALUATED.
C                N      - INPUT.  ORDER OF MATRIX A.
C                IA     - INPUT.  ROW DIMENSION OF MATRIX A EXACTLY
C                           AS SPECIFIED IN THE DIMENSION STATEMENT IN
C                           THE CALLING PROGRAM.
C                COEF   - INPUT VECTOR OF LENGTH M CONTAINING THE
C                           POLYNOMIAL COEFFICIENTS.  COEF(M+1-K)
C                           MULTIPLIES A TO THE POWER K-1, FOR
C                           K=1,2,...,M.
C                M      - INPUT.  NUMBER OF COEFFICIENTS.
C                B      - OUTPUT N BY N MATRIX CONTAINING THE VALUE OF
C                           THE POLYNOMIAL.
C                IB     - INPUT.  ROW DIMENSION OF MATRIX B EXACTLY
C                           AS SPECIFIED IN THE DIMENSION STATEMENT IN
C                           THE CALLING PROGRAM.
C                WKAREA - WORK AREA OF LENGTH AT LEAST N.
C
C   REQD. IMSL ROUTINES - SINGLE/NONE REQUIRED
C                       - DOUBLE/VXADD,VXMUL,VXSTO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE VPOLYF (A,N,IA,COEF,M,B,IB,WKAREA)
C
      DIMENSION          A(IA,N),COEF(M),B(IB,N),WKAREA(N)
      DOUBLE PRECISION   TEMP
      DATA               ZERO/0.0/
C                                  FIRST EXECUTABLE STATEMENT
      DO 5 I=1,N
         DO 5 J=1,N
            B(I,J) = ZERO
    5 CONTINUE
      DO 10 I=1,N
   10 B(I,I)=COEF(1)
      IF(M .LT. 2) GO TO 35
      DO 30 K=2,M
         DO 25 K1=1,N
            DO 20 K2=1,N
               TEMP = 0.0D0
               DO 15 K3 = 1,N
                  TEMP = DBLE(B(K1,K3))*DBLE(A(K3,K2))+TEMP
   15          CONTINUE
               WKAREA(K2) = TEMP
   20       CONTINUE
            DO 25 K0=1,N
   25    B(K1,K0)=WKAREA(K0)
         DO 30 I=1,N
            B(I,I)=B(I,I)+COEF(K)
   30 CONTINUE
   35 CONTINUE
      RETURN
      END
