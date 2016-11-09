C   IMSL ROUTINE NAME   - LEQOF
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1982
C
C   PURPOSE             - LINEAR EQUATION SOLUTION - FULL MATRICES
C                           (OUT-OF-CORE VERSION)
C
C   USAGE               - CALL LEQOF(IUNITA,IUNITW,N,MA,B,IB,M,IJOB,WK,
C                                    IER)
C
C   ARGUMENTS    IUNITA - UNIT IDENTIFIER OF THE FILE WHERE MATRIX
C                           A IS STORED (INPUT).  A IS ASSUMED TO BE
C                           STORED IN UNFORMATTED FORM, I.E. AS OUTPUT
C                           BY UNFORMATTED WRITE STATEMENTS.  THE FIRST
C                           MA COLUMNS OF A MUST BE STORED ON THE FIRST
C                           RECORD, THE SECOND MA COLUMNS ON THE SECOND
C                           RECORD, ETC.  THE RECORD LENGTH MUST BE
C                           N*MA WORKING PRECISION WORDS.  THE
C                           LAST RECORD NEED NOT BE FULL, I.E. N NEED
C                           NOT BE DIVISIBLE BY MA.
C                         ON OUTPUT, THE MATRIX ON UNIT IUNITA WILL BE
C                           OVERWRITTEN BY THE L-U DECOMPOSITION OF
C                           A ROWWISE PERMUTATION OF THE ORIGINAL MATRIX
C                IUNITW - UNIT IDENTIFIER OF ANOTHER FILE WITH THE
C                           SAME CHARACTERISTICS (RECORD LENGTH, ETC)
C                           AS IUNITA (INPUT).  WHEN IJOB=0, THIS FILE
C                           IS USED FOR WORKING SPACE.  WHEN IJOB=1,
C                           THIS FILE IS NOT USED.
C                N      - ORDER OF MATRIX A. (INPUT)
C                MA     - NUMBER OF COLUMNS STORED PER RECORD (INPUT).
C                B      - INPUT N BY M MATRIX CONTAINING THE M
C                           RIGHT HAND SIDES OF THE EQUATION AX = B.
C                         ON OUTPUT, THE SOLUTION MATRIX X REPLACES B.
C                           (INPUT/OUTPUT)
C                IB     - ROW DIMENSION OF MATRIX B EXACTLY AS SPECIFIED
C                           IN THE DIMENSION STATEMENT OF THE CALLING
C                           PROGRAM (INPUT).
C                M      - NUMBER OF RIGHT HAND SIDES (COLUMNS IN B).
C                           (INPUT)
C                IJOB   - OPTION PARAMETER. (INPUT)  IJOB=I IMPLIES WHEN
C                           I=0, FACTOR THE MATRIX AND SOLVE THE
C                             EQUATION AX=B.
C                           I=1, SOLVE THE EQUATION AX=B.  THIS
C                             OPTION IMPLIES THAT LEQOF HAS ALREADY
C                             BEEN CALLED USING IJOB=0 SO THAT
C                             THE MATRIX STORED ON UNIT IUNITA HAS
C                             ALREADY BEEN FACTORED.
C                WK     - REAL WORK AREA OF LENGTH 2*N*(MA+1).
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
      SUBROUTINE LEQOF (IUNITA,IUNITW,N,MA,B,IB,M,IJOB,WK,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            IUNITA,IUNITW,N,MA,M,IB,IJOB,IER
      REAL               B(IB,M),WK(N,2,1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            IA,IBIG,IFIN,IFOUT,IIK,IKP1,IK,IP1,I,JA,JJJ,JJ,
     *                   J,KBP1,KB,KKB,KK,LB,LDBLE,LIMK,LIML,NBLOCK,N2,
     *                   MAP1
      REAL               ABIG,AM,REPS,TEMP
      DATA               REPS/Z3C100000/
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
      LDBLE = 1
      N2 = 2*N
      MAP1 = MA+1
      REWIND IUNITA
      REWIND IUNITW
      DO 5 I=1,N
    5 WK(I,1,MAP1) = I
      NBLOCK = (N-1)/MA+1
C                                  BEGIN GAUSSIAN ELIMINATION
      DO 90 KB=1,NBLOCK
         IFIN = IUNITA
         IF (IJOB.EQ.1) GO TO 10
         IFOUT = IUNITW
         IF (MOD(KB,2).EQ.0) IFIN = IUNITW
         IF (MOD(KB,2).EQ.0) IFOUT = IUNITA
   10    LIMK = MIN0(MA,N-(KB-1)*MA)
C                                  READ A BLOCK OF COLUMNS INTO FIRST
C                                    BUFFER
         READ (IFIN) ((WK(I,1,J),I=1,N),J=1,LIMK)
         IF (KB.GT.1 .OR. IJOB.EQ.1) GO TO 25
         DO 20 J=1,LIMK
            TEMP = 0.0
            DO 15 I=1,N
               TEMP = AMAX1(TEMP,ABS(WK(I,1,J)))
   15       CONTINUE
            WK(J,2,MAP1) = TEMP
   20    CONTINUE
   25    CONTINUE
C                                  FACTOR FIRST BUFFER
         DO 45 IK=1,LIMK
            I = (KB-1)*MA+IK
            ABIG = 0.0
            DO 30 J=I,N
               JA = WK(J,1,MAP1)
               IF (ABS(WK(JA,1,IK)).LE.ABIG) GO TO 30
               ABIG = ABS(WK(JA,1,IK))
               IBIG = J
   30       CONTINUE
            IF (ABIG.EQ.0.0) GO TO 9000
            IF (IJOB.EQ.1) GO TO 35
            IF (ABIG.LE.10.0*REPS*WK(I,2,MAP1)) GO TO 9000
   35       CONTINUE
            TEMP = WK(IBIG,1,MAP1)
            WK(IBIG,1,MAP1) = WK(I,1,MAP1)
            WK(I,1,MAP1) = TEMP
            IF (I.GE.N) GO TO 45
            IA = WK(I,1,MAP1)
            IP1 = I+1
            DO 40 J=IP1,N
               JA = WK(J,1,MAP1)
               AM = WK(JA,1,IK)/WK(IA,1,IK)
               IF (AM.EQ.0.0) GO TO 40
               CALL SAXPY(M,-AM,B(IA,1),IB,B(JA,1),IB)
               IF (IJOB.EQ.1) GO TO 40
               IF (IK.GE.LIMK) GO TO 40
               IKP1 = IK+1
               CALL SAXPY(LIMK-IK,-AM,WK(IA,1,IKP1),N2,WK(JA,1,IKP1),N2)
   40       CONTINUE
   45    CONTINUE
         IF (IJOB.EQ.1) GO TO 90
         WRITE (IFOUT) ((WK(I,1,J),I=1,N),J=1,LIMK)
         IF (KB.GE.NBLOCK) GO TO 80
         KBP1 = KB+1
         DO 75 LB=KBP1,NBLOCK
            LIML = MIN0(MA,N-(LB-1)*MA)
C                                  READ A BLOCK OF COLUMNS INTO
C                                    SECOND BUFFER
            READ (IFIN) ((WK(I,2,J),I=1,N),J=1,LIML)
            IF (KB.GT.1) GO TO 60
            DO 55 KK=1,LIML
               J = (LB-1)*MA+KK
               TEMP = 0.0
               DO 50 I=1,N
                  TEMP = AMAX1(TEMP,ABS(WK(I,2,KK)))
   50          CONTINUE
               WK(J,2,MAP1) = TEMP
   55       CONTINUE
   60       CONTINUE
C                                  DO ELIMINATION ON SECOND BUFFER
C                                    USING FACTORS SAVED IN FIRST
            DO 70 IK=1,LIMK
               I = (KB-1)*MA+IK
               IF (I.GE.N) GO TO 70
               IA = WK(I,1,MAP1)
               IP1 = I+1
               DO 65 J=IP1,N
                  JA = WK(J,1,MAP1)
                  AM = WK(JA,1,IK)/WK(IA,1,IK)
                  IF (AM.EQ.0.0) GO TO 65
                  CALL SAXPY(LIML,-AM,WK(IA,2,1),N2,WK(JA,2,1),N2)
   65          CONTINUE
   70       CONTINUE
            WRITE (IFOUT) ((WK(I,2,J),I=1,N),J=1,LIML)
   75    CONTINUE
   80    CONTINUE
         DO 85 LB=KB,NBLOCK
            BACKSPACE IUNITA
            BACKSPACE IUNITW
   85    CONTINUE
         READ (IFOUT) ((WK(I,1,J),I=1,N),J=1,LIMK)
         WRITE (IFIN) ((WK(I,1,J),I=1,N),J=1,LIMK)
   90 CONTINUE
C                                  BACK SUBSTITUTION
      DO 115 KKB=1,NBLOCK
         KB = NBLOCK+1-KKB
         LIMK = MIN0(MA,N-(KB-1)*MA)
         BACKSPACE IUNITA
         READ (IUNITA) ((WK(I,1,J),I=1,N),J=1,LIMK)
         BACKSPACE IUNITA
         DO 110 IIK=1,LIMK
            IK = LIMK+1-IIK
            I = (KB-1)*MA+IK
            IA = WK(I,1,MAP1)
            TEMP = WK(IA,1,IK)
            DO 95 JJJ=1,M
               B(IA,JJJ) = B(IA,JJJ)/TEMP
   95       CONTINUE
            IF (I.EQ.1) GO TO 110
            DO 105 JJ=2,I
               J = I+1-JJ
               JA = WK(J,1,MAP1)
               TEMP = WK(JA,1,IK)
               DO 100 JJJ=1,M
                  B(JA,JJJ) = B(JA,JJJ)-TEMP*B(IA,JJJ)
  100          CONTINUE
  105       CONTINUE
  110    CONTINUE
  115 CONTINUE
C                                  SORT SOLUTION VECTOR
      DO 130 JJJ=1,M
         DO 120 I=1,N
            IA = WK(I,1,MAP1)
            WK(I,1,1) = B(IA,JJJ)
  120    CONTINUE
         DO 125 I=1,N
            B(I,JJJ) = WK(I,1,1)
  125    CONTINUE
  130 CONTINUE
      GO TO 9005
 9000 IER = 129
      CALL UERTST(IER,6HLEQOF )
 9005 RETURN
      END
