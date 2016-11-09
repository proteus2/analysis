C   IMSL ROUTINE NAME   - FFT3D
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1980
C
C   PURPOSE             - COMPUTE THE FAST FOURIER TRANSFORM OF
C                           A COMPLEX VALUED 1,2 OR 3 DIMENSIONAL
C                           ARRAY
C
C   USAGE               - CALL FFT3D (A,IA1,IA2,N1,N2,N3,IJOB,IWK,RWK,
C                           CWK)
C
C   ARGUMENTS    A      - COMPLEX ARRAY. A MAY BE A THREE
C                           DIMENSIONAL ARRAY OF DIMENSION N1 BY N2
C                           BY N3, A TWO DIMENSIONAL ARRAY OF
C                           DIMENSION N1 BY N2, OR A VECTOR OF
C                           LENGTH N1. ON INPUT A CONTAINS THE
C                           ARRAY TO BE TRANSFORMED. ON OUTPUT
C                           A IS REPLACED BY THE FOURIER OR
C                           INVERSE FOURIER TRANSFORM (DEPENDING ON
C                           THE VALUE OF INPUT PARAMETER IJOB).
C                IA1    - FIRST DIMENSION OF THE ARRAY A EXACTLY
C                           AS SPECIFIED IN THE DIMENSION STATEMENT
C                           IN THE CALLING PROGRAM. (INPUT)
C                IA2    - SECOND DIMENSION OF THE ARRAY A EXACTLY
C                           AS SPECIFIED IN THE DIMENSION STATEMENT
C                           IN THE CALLING PROGRAM. (INPUT)
C                N1     - LIMITS ON THE FIRST, SECOND, AND THIRD
C                N2         SUBSCRIPTS OF THE ARRAY A, RESPECTIVELY.
C                N3         (INPUT)
C                IJOB   - INPUT OPTION PARAMETER.
C                           IF IJOB IS POSITIVE, THE FAST FOURIER
C                             TRANSFORM OF A IS TO BE CALCULATED.
C                           IF IJOB IS NEGATIVE, THE INVERSE
C                             FAST FOURIER TRANSFORM OF A IS TO BE
C                             CALCULATED.
C                IWK    - INTEGER WORK VECTOR OF LENGTH
C                           6*MAX(N1,N2,N3)+150.
C                RWK    - REAL WORK VECTOR OF LENGTH
C                           6*MAX(N1,N2,N3)+150.
C                CWK    - COMPLEX WORK VECTOR OF LENGTH
C                           MAX(N2,N3).
C
C   REQD. IMSL ROUTINES - FFTCC
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS  1.  IF IJOB IS POSITIVE, FFT3D CALCULATES THE FOURIER
C                TRANSFORM, X, ACCORDING TO THE FOLLOWING FORMULA
C
C                  X(I+1,J+1,K+1)=TRIPLE SUM OF A(L+1,M+1,N+1)*
C                  EXP(2*PI*SQRT(-1)*(I*L/N1+J*M/N2+K*N/N3))
C                  WITH L=0...N1-1, M=0...N2-1, N=0...N3-1
C                  AND PI=3.1415...
C
C                IF IJOB IS NEGATIVE, FFT3D CALCULATES THE INVERSE
C                FOURIER TRANSFORM, X, ACCORDING TO THE FOLLOWING
C                FORMULA
C
C                  X(I+1,J+1,K+1)=1/(N1*N2*N3)*TRIPLE SUM OF
C                  A(L+1,M+1,N+1)*
C                  EXP(-2*PI*SQRT(-1)*(I*L/N1+J*M/N2+K*N/N3))
C                  WITH L=0...N1-1, M=0...N2-1, N=0...N3-1
C                  AND PI=3.1415...
C
C                NOTE THAT X OVERWRITES A ON OUTPUT.
C            2.  IF A IS A TWO DIMENSIONAL ARRAY, SET N3 = 1.
C                IF A IS A ONE DIMENSIONAL ARRAY (VECTOR),
C                SET IA2 = N2 = N3 = 1.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE FFT3D  (A,IA1,IA2,N1,N2,N3,IJOB,IWK,RWK,CWK)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            IA1,IA2,N1,N2,N3,IJOB,IWK(1)
      REAL               RWK(1)
      COMPLEX            A(IA1,IA2,N3),CWK(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,J,K,L,M,N
      REAL               R123
      COMPLEX            C123
C                                  FIRST EXECUTABLE STATEMENT
      IF (IJOB.GT.0) GO TO 10
C                                  INVERSE TRANSFORM
      DO 5 I=1,N1
      DO 5 J=1,N2
      DO 5 K=1,N3
         A(I,J,K) = CONJG(A(I,J,K))
    5 CONTINUE
C                                  TRANSFORM THIRD SUBSCRIPT
   10 DO 25 L=1,N1
      DO 25 M=1,N2
         DO 15 N=1,N3
            CWK(N) = A(L,M,N)
   15    CONTINUE
         CALL FFTCC (CWK,N3,IWK,RWK)
         DO 20 K=1,N3
            A(L,M,K) = CWK(K)
   20    CONTINUE
   25 CONTINUE
C                                  TRANSFORM SECOND SUBSCRIPT
      DO 40 L=1,N1
      DO 40 K=1,N3
         DO 30 M=1,N2
            CWK(M) = A(L,M,K)
   30    CONTINUE
         CALL FFTCC (CWK,N2,IWK,RWK)
         DO 35 J=1,N2
            A(L,J,K) = CWK(J)
   35    CONTINUE
   40 CONTINUE
C                                  TRANSFORM FIRST SUBSCRIPT
      DO 45 J=1,N2
      DO 45 K=1,N3
         CALL FFTCC (A(1,J,K),N1,IWK,RWK)
   45 CONTINUE
      IF (IJOB.GT.0) GO TO 55
      R123 = N1*N2*N3
      C123 = CMPLX(R123,0.0)
      DO 50 I=1,N1
      DO 50 J=1,N2
      DO 50 K=1,N3
         A(I,J,K) = CONJG(A(I,J,K))/C123
   50 CONTINUE
   55 RETURN
      END
