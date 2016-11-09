C   IMSL ROUTINE NAME   - LEQT2F
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1982
C
C   PURPOSE             - LINEAR EQUATION SOLUTION - FULL STORAGE
C                           MODE - HIGH ACCURACY SOLUTION
C
C   USAGE               - CALL LEQT2F (A,M,N,IA,B,IDGT,WKAREA,IER)
C
C   ARGUMENTS    A      - INPUT MATRIX OF DIMENSION N BY N CONTAINING
C                           THE COEFFICIENT MATRIX OF THE EQUATION
C                           AX = B.
C                M      - NUMBER OF RIGHT-HAND SIDES. (INPUT)
C                N      - ORDER OF A AND NUMBER OF ROWS IN B. (INPUT)
C                IA     - ROW DIMENSION OF A AND B EXACTLY AS SPECIFIED
C                           IN THE DIMENSION STATEMENT IN THE CALLING
C                           PROGRAM. (INPUT)
C                B      - INPUT MATRIX OF DIMENSION N BY M CONTAINING
C                           THE RIGHT-HAND SIDES OF THE EQUATION AX = B.
C                           ON OUTPUT, THE N BY M MATRIX OF SOLUTIONS
C                           REPLACES B.
C                IDGT   - INPUT OPTION.
C                         IF IDGT IS GREATER THAN 0, THE ELEMENTS OF
C                           A AND B ARE ASSUMED TO BE CORRECT TO IDGT
C                           DECIMAL DIGITS AND THE ROUTINE PERFORMS
C                           AN ACCURACY TEST.
C                         IF IDGT EQUALS 0, THE ACCURACY TEST IS
C                           BYPASSED.
C                         ON OUTPUT, IDGT CONTAINS THE APPROXIMATE
C                           NUMBER OF DIGITS IN THE ANSWER WHICH
C                           WERE UNCHANGED AFTER IMPROVEMENT.
C                WKAREA - WORK AREA OF DIMENSION GREATER THAN OR EQUAL
C                           TO N**2+3N.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         WARNING ERROR
C                           IER = 34 INDICATES THAT THE ACCURACY TEST
C                             FAILED. THE COMPUTED SOLUTION MAY BE IN
C                             ERROR BY MORE THAN CAN BE ACCOUNTED FOR
C                             BY THE UNCERTAINTY OF THE DATA. THIS
C                             WARNING CAN BE PRODUCED ONLY IF IDGT IS
C                             GREATER THAN 0 ON INPUT. (SEE THE
C                             CHAPTER L PRELUDE FOR FURTHER DISCUSSION.)
C                         TERMINAL ERROR
C                           IER = 129 INDICATES THAT THE MATRIX IS
C                             ALGORITHMICALLY SINGULAR. (SEE THE
C                             CHAPTER L PRELUDE).
C                           IER = 131 INDICATES THAT THE MATRIX IS TOO
C                             ILL-CONDITIONED FOR ITERATIVE IMPROVEMENT
C                             TO BE EFFECTIVE.
C
C   REQD. IMSL ROUTINES - SINGLE/LUDATN,LUELMN,LUREFN,UERTST,UGETIO
C                       - DOUBLE/LUDATN,LUELMN,LUREFN,UERTST,UGETIO,
C                           VXADD,VXMUL,VXSTO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE LEQT2F (A,M,N,IA,B,IDGT,WKAREA,IER)
C
      DIMENSION          A(IA,1),B(IA,1),WKAREA(1)
C                                  FIRST EXECUTABLE STATEMENT
C                                  INITIALIZE IER
      IER=0
      JER=0
      J = N*N+1
      K = J+N
      MM = K+N
      KK = 0
      MM1 = MM-1
      JJ=1
      DO 5 L=1,N
         DO 5 I=1,N
            WKAREA(JJ)=A(I,L)
            JJ=JJ+1
    5 CONTINUE
C                                  DECOMPOSE A
      CALL LUDATN (WKAREA,N,N,A,IA,IDGT,D1,D2,WKAREA(J),WKAREA(K),
     *             WA,IER)
      IF (IER.GT.128) GO TO 25
      IF (IDGT .EQ. 0 .OR. IER .NE. 0) KK = 1
      DO 15 I = 1,M
C                                  PERFORMS THE ELIMINATION PART OF
C                                  AX = B
         CALL LUELMN (A,IA,N,B(1,I),WKAREA(J),WKAREA(MM))
C                                  REFINEMENT OF SOLUTION TO AX = B
         IF (KK .NE. 0)
     *   CALL LUREFN (WKAREA,N,N,A,IA,B(1,I),IDGT,WKAREA(J),WKAREA(MM),
     *                WKAREA(K),WKAREA(K),JER)
         DO 10 II=1,N
            B(II,I) = WKAREA(MM1+II)
   10    CONTINUE
         IF (JER.NE.0) GO TO 20
   15 CONTINUE
      GO TO 25
   20 IER = 131
   25 JJ=1
      DO 30 J = 1,N
         DO 30 I = 1,N
            A(I,J)=WKAREA(JJ)
            JJ=JJ+1
   30 CONTINUE
      IF (IER .EQ. 0) GO TO 9005
 9000 CONTINUE
      CALL UERTST (IER,6HLEQT2F)
 9005 RETURN
      END
 
R; T=0.03/0.28 22:00:29
