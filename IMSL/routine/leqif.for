C   IMSL ROUTINE NAME   - LEQIF
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1982
C
C   PURPOSE             - LINEAR EQUATION SOLUTION - FULL MATRICES
C                           (VIRTUAL MEMORY VERSION)
C
C   USAGE               - CALL LEQIF(A,IA,N,MA,B,IB,M,IJOB,WK,IER)
C
C   ARGUMENTS    A      - INPUT N BY N MATRIX CONTAINING THE
C                           COEFFICIENT MATRIX OF THE EQUATION AX=B.
C                         ON OUTPUT, A IS REPLACED BY THE LU
C                           DECOMPOSITION OF A ROWWISE PERMUTATION
C                           OF A.   (INPUT/OUTPUT)
C                IA     - ROW DIMENSION OF A EXACTLY AS SPECIFIED
C                           IN THE DIMENSION STATEMENT OF THE CALLING
C                           PROGRAM. (INPUT)
C                N      - ORDER OF MATRIX A. (INPUT)
C                MA     - NUMBER OF COLUMNS PER BLOCK (INPUT).  THE
C                           CHOICE OF MA WILL AFFECT THE SOLUTION
C                           SPEED AS FOLLOWS.  AS MA IS INCREASED,
C                           THE ALGORITHM WILL RUN FASTER, UNTIL A
C                           POINT IS REACHED BEYOND WHICH 2*MA*IA
C                           WORKING PRECISION WORDS CANNOT BE HELD
C                           IN MAIN MEMORY, WITHOUT PAGING.  MA
C                           MUST BE LESS THAN OR EQUAL TO N.
C                B      - INPUT N BY M MATRIX CONTAINING THE M
C                           RIGHT HAND SIDES OF THE EQUATION AX = B.
C                         ON OUTPUT, THE SOLUTION MATRIX X REPLACES B.
C                           (INPUT/OUTPUT)
C                IB     - ROW DIMENSION OF B EXACTLY AS SPECIFIED
C                           IN THE DIMENSION STATEMENT OF THE CALLING
C                           PROGRAM. (INPUT)
C                M      - NUMBER OF RIGHT HAND SIDES (COLUMNS IN B).
C                           (INPUT)
C                IJOB   - OPTION PARAMETER.  (INPUT) IJOB=I IMPLIES WHEN
C                           I=0, FACTOR THE MATRIX AND SOLVE THE
C                             EQUATION AX=B.
C                           I=1, SOLVE THE EQUATION AX=B.  THIS
C                             OPTION IMPLIES THAT LEQIF HAS ALREADY
C                             BEEN CALLED USING IJOB=0 SO THAT
C                             THE MATRIX A HAS ALREADY BEEN
C                             FACTORED, AND THAT WK HAS NOT
C                             BEEN ALTERED SINCE THAT CALL.
C                WK     - REAL WORK AREA OF LENGTH 3*N.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER=129 INDICATES THAT MATRIX A IS
C                             ALGORITHMICALLY SINGULAR.  (SEE THE
C                             CHAPTER L PRELUDE.)
C
C   REQD. IMSL ROUTINES - SINGLE/VBLA=SAXPY,UERTST,UGETIO
C                         DOUBLE/VBLA=DAXPY,UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE LEQIF (A,IA,N,MA,B,IB,M,IJOB,WK,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            IA,IB,N,MA,M,IJOB,IER
      REAL               A(IA,N),B(IB,M),WK(N,3)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,IBIG,IIK,IP1,IZ,J,JJ,JJJ,JZ,KB,KBP1,
     *                   KKB,LB,LIMK,LIMK0,LIMK1,LIMK2,LIML,LIML0,
     *                   LIML1,LIML2,NBLOCK
      REAL               ABIG,AM,REPS,TEMP
      DATA               REPS/Z3C100000/
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
      DO 5 I=1,N
    5 WK(I,1) = I
      NBLOCK = (N-1)/MA+1
C                                  BEGIN GAUSSIAN ELIMINATION
      DO 80 KB=1,NBLOCK
         LIMK0 = (KB-1)*MA
         LIMK1 = LIMK0+1
         LIMK2 = MIN0(LIMK0+MA,N)
         IF (KB.GT.1 .OR. IJOB.EQ.1) GO TO 20
         DO 15 J=LIMK1,LIMK2
            TEMP = 0.0
            DO 10 I=1,N
               TEMP = AMAX1(TEMP,ABS(A(I,J)))
   10       CONTINUE
            WK(J,2) = TEMP
   15    CONTINUE
   20    CONTINUE
C                                  FACTOR FIRST BUFFER
         DO 40 I=LIMK1,LIMK2
            ABIG = 0.0
            DO 25 J=I,N
               JZ = WK(J,1)
               IF (ABS(A(JZ,I)).LE.ABIG) GO TO 25
               ABIG = ABS(A(JZ,I))
               IBIG = J
   25       CONTINUE
            IF (ABIG.EQ.0.0) GO TO 9000
            IF (IJOB.EQ.1) GO TO 30
            IF (ABIG.LE.10.0*REPS*WK(I,2)) GO TO 9000
   30       CONTINUE
            TEMP = WK(IBIG,1)
            WK(IBIG,1) = WK(I,1)
            WK(I,1) = TEMP
            IF (I.GE.N) GO TO 40
            IZ = WK(I,1)
            IP1 = I+1
            DO 35 J=IP1,N
               JZ = WK(J,1)
               AM = A(JZ,I)/A(IZ,I)
               IF (AM.EQ.0.0) GO TO 35
               CALL SAXPY(M,-AM,B(IZ,1),IA,B(JZ,1),IA)
               IF (IJOB.EQ.1) GO TO 35
               IF (I.GE.LIMK2) GO TO 35
               CALL SAXPY(LIMK2-I,-AM,A(IZ,I+1),IA,A(JZ,I+1),IA)
   35       CONTINUE
   40    CONTINUE
         IF (IJOB.EQ.1) GO TO 80
         IF (KB.GE.NBLOCK) GO TO 75
         KBP1 = KB+1
         DO 70 LB=KBP1,NBLOCK
            LIML0 = (LB-1)*MA
            LIML1 = LIML0+1
            LIML2 = MIN0(LIML0+MA,N)
            LIML = LIML2-LIML0
            IF (KB.GT.1) GO TO 55
            DO 50 J=LIML1,LIML2
               TEMP = 0.0
               DO 45 I=1,N
                  TEMP = AMAX1(TEMP,ABS(A(I,J)))
   45          CONTINUE
               WK(J,2) = TEMP
   50       CONTINUE
   55       CONTINUE
C                                  DO ELIMINATION ON SECOND BLOCK
C                                    USING FACTORS SAVED IN FIRST
            DO 65 I=LIMK1,LIMK2
               IF (I.GE.N) GO TO 65
               IZ = WK(I,1)
               IP1 = I+1
               DO 60 J=IP1,N
                  JZ = WK(J,1)
                  AM = A(JZ,I)/A(IZ,I)
                  IF (AM.EQ.0.0) GO TO 60
                  CALL SAXPY(LIML,-AM,A(IZ,LIML1),IA,A(JZ,LIML1),IA)
   60          CONTINUE
   65       CONTINUE
   70    CONTINUE
   75    CONTINUE
   80 CONTINUE
C                                  BACK SUBSTITUTION
      DO 105 KKB=1,NBLOCK
         KB = NBLOCK+1-KKB
         LIMK0 = (KB-1)*MA
         LIMK2 = MIN0(LIMK0+MA,N)
         LIMK = LIMK2-LIMK0
         DO 100 IIK=1,LIMK
            I = LIMK2+1-IIK
            IZ = WK(I,1)
            TEMP = A(IZ,I)
            DO 85 JJJ=1,M
               B(IZ,JJJ) = B(IZ,JJJ)/TEMP
   85       CONTINUE
            IF (I.EQ.1) GO TO 100
            DO 95 JJ=2,I
               J = I+1-JJ
               JZ = WK(J,1)
               TEMP = A(JZ,I)
               DO 90 JJJ=1,M
                  B(JZ,JJJ) = B(JZ,JJJ)-TEMP*B(IZ,JJJ)
   90          CONTINUE
   95       CONTINUE
  100    CONTINUE
  105 CONTINUE
C                                  SORT SOLUTION VECTOR
      DO 120 JJJ=1,M
         DO 110 I=1,N
            IZ = WK(I,1)
            WK(I,3) = B(IZ,JJJ)
  110    CONTINUE
         DO 115 I=1,N
            B(I,JJJ) = WK(I,3)
  115    CONTINUE
  120 CONTINUE
      GO TO 9005
 9000 IER = 129
      CALL UERTST(IER,6HLEQIF )
 9005 RETURN
      END
