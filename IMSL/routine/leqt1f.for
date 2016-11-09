C   IMSL ROUTINE NAME   - LEQT1F
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1982
C
C   PURPOSE             - LINEAR EQUATION SOLUTION - FULL STORAGE
C                           MODE - SPACE ECONOMIZER SOLUTION.
C
C   USAGE               - CALL LEQT1F (A,M,N,IA,B,IDGT,WKAREA,IER)
C
C   ARGUMENTS    A      - INPUT MATRIX OF DIMENSION N BY N CONTAINING
C                           THE COEFFICIENT MATRIX OF THE EQUATION
C                           AX = B.
C                         ON OUTPUT, A IS REPLACED BY THE LU
C                           DECOMPOSITION OF A ROWWISE PERMUTATION OF
C                           A.
C                M      - NUMBER OF RIGHT-HAND SIDES. (INPUT)
C                N      - ORDER OF A AND NUMBER OF ROWS IN B. (INPUT)
C                IA     - ROW DIMENSION OF A AND B EXACTLY AS SPECIFIED
C                           IN THE DIMENSION STATEMENT OF THE CALLING
C                           PROGRAM. (INPUT)
C                B      - INPUT MATRIX OF DIMENSION N BY M CONTAINING
C                           RIGHT-HAND SIDES OF THE EQUATION AX = B.
C                         ON OUTPUT, THE N BY M SOLUTION X REPLACES B.
C                IDGT   - INPUT OPTION.
C                         IF IDGT IS GREATER THAN 0, THE ELEMENTS OF
C                           A AND B ARE ASSUMED TO BE CORRECT TO IDGT
C                           DECIMAL DIGITS AND THE ROUTINE PERFORMS
C                           AN ACCURACY TEST.
C                         IF IDGT EQUALS ZERO, THE ACCURACY TEST IS
C                           BYPASSED.
C                WKAREA - WORK AREA OF DIMENSION GREATER THAN OR EQUAL
C                           TO N.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER = 129 INDICATES THAT MATRIX A IS
C                             ALGORITHMICALLY SINGULAR. (SEE THE
C                             CHAPTER L PRELUDE).
C                         WARNING ERROR
C                           IER = 34 INDICATES THAT THE ACCURACY TEST
C                             FAILED.  THE COMPUTED SOLUTION MAY BE IN
C                             ERROR BY MORE THAN CAN BE ACCOUNTED FOR
C                             BY THE UNCERTAINTY OF THE DATA.  THIS
C                             WARNING CAN BE PRODUCED ONLY IF IDGT IS
C                             GREATER THAN 0 ON INPUT.  (SEE CHAPTER L
C                             PRELUDE FOR FURTHER DISCUSSION).
C
C   REQD. IMSL ROUTINES - LUDATN,LUELMN,UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE LEQT1F (A,M,N,IA,B,IDGT,WKAREA,IER)
C
      DIMENSION          A(IA,1),B(IA,1),WKAREA(1)
C                                  INITIALIZE IER
C                                  FIRST EXECUTABLE STATEMENT
      IER=0
C                                  DECOMPOSE A
      CALL LUDATN (A,IA,N,A,IA,IDGT,D1,D2,WKAREA,WKAREA,WA,IER)
      IF (IER .GT. 128) GO TO 9000
C                                  CALL ROUTINE LUELMN (FORWARD AND
C                                  BACKWARD SUBSTITUTIONS)
      DO 10 J=1,M
         CALL LUELMN (A,IA,N,B(1,J),WKAREA,B(1,J))
   10 CONTINUE
      IF (IER .EQ. 0) GO TO 9005
 9000 CONTINUE
      CALL UERTST (IER,6HLEQT1F)
 9005 RETURN
      END
