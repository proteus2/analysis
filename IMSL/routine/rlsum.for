C   IMSL ROUTINE NAME   - RLSUM
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - REORDERING OF THE ROWS AND CORRESPONDING
C                           COLUMNS OF A SYMMETRIC MATRIX STORED IN
C                           SYMMETRIC STORAGE MODE
C
C   USAGE               - CALL RLSUM (AA,MM,IH,M,A,IER)
C
C   ARGUMENTS    AA     - INPUT VECTOR OF LENGTH MM*(MM+1)/2 CONTAINING
C                           THE MM BY MM SYMMETRIC MATRIX STORED IN
C                           SYMMETRIC STORAGE MODE.
C                MM     - INPUT ORDER OF MATRIX AA.
C                IH     - INPUT AND OUTPUT VECTOR OF LENGTH 2*M+MM.
C                         ON INPUT, THE FIRST M LOCATIONS ARE USED TO
C                           INDICATE HOW THE ROWS AND CORRESPONDING
C                           COLUMNS OF AA ARE TO BE REORDERED AND
C                           RETURNED IN A. IH(J) = K, FOR ANY J IN
C                           1,...,M IMPLIES THAT ROW-COLUMN K OF MATRIX
C                           AA WILL BE ROW-COLUMN J OF MATRIX A.
C                           THESE VALUES MUST CONSIST OF A SUBSET OF THE
C                           INTEGERS 1,2,...,MM.
C                         ON OUTPUT THE LAST MM LOCATIONS SPECIFY HOW
C                           THE ENTIRE AA MATRIX HAS BEEN REORDERED AND
C                           RETURNED IN A, USING THE SAME CONVENTION AS
C                           IS USED FOR THE FIRST M ELEMENTS OF IH.
C                         THE REMAINING LOCATIONS ARE WORK STORAGE.
C                M      - INPUT NUMBER, NOT EXCEEDING MM, OF ROW-
C                           COLUMNS OF AA THAT ARE TO BE REORDERED AND
C                           RETURNED AS THE INITIAL SET OF ROW-COLUMNS
C                           OF A.
C                A      - OUTPUT VECTOR OF LENGTH MM*(MM+1)/2 CONTAINING
C                           THE MM BY MM REORDERED SYMMETRIC MATRIX
C                           STORED IN SYMMETRIC STORAGE MODE.
C                           THE SUBSET OF ROW-COLUMNS SPECIFIED BY THE
C                           FIRST M ELEMENTS OF IH APPEARS AS A
C                           SYMMETRIC SUBMATRIX STORED IN THE FIRST
C                           M*(M+1)/2 LOCATIONS OF A.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER=129 INDICATES THAT THE FIRST M ELEMENTS
C                             OF VECTOR IH DID NOT CONSTITUTE A SUBSET
C                             OF THE INTEGERS (1,2,...,MM) OR THAT M
C                             WAS NOT POSITIVE.
C
C   REQD. IMSL ROUTINES - UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS  1.  RLSUM IS INTENDED PRIMARILY FOR USE IN ASSOCIATION
C                WITH IMSL ROUTINES BECOVM AND RLMUL. IN SOME CASES
C                USAGE IN ASSOCIATION WITH BECOVM AND RLSTP MAY BE
C                APPROPRIATE.
C            2.  THE ELEMENTS OF IH CONTROL THE ROW-COLUMN SELECTING
C                AND REORDERING. BOTH THE POSITION AND VALUE OF THE
C                ELEMENTS OF IH ARE IMPORTANT. VALUE SPECIFIES THE
C                THE ROW-COLUMN NUMBER IN INPUT AA. POSITION SPECIFIES
C                THE CORRESPONDING ROW-COLUMN NUMBER IN OUTPUT A.
C            3.  THE REORDERING MAY BE PERFORMED IN PLACE BY USING THE
C                SAME ACTUAL ARGUMENTS FOR DUMMY ARGUMENTS AA AND A.
C            4.  GENERALIZING THE REGRESSION APPLICATION DISCUSSED IN
C                THE ALGORITHM SECTION IN THE MANUAL DOCUMENT AND
C                REQUIRING THAT THE M-TH VARIABLE IN OUTPUT A BE THE
C                DEPENDENT (RESPONSE) VARIABLE, THE FOLLOWING RESULTS
C                ARE AVAILABLE.
C                (A)  THE SUMS OF CROSSPRODUCTS VECTOR REQUIRED FOR THE
C                     RIGHT HAND SIDE OF THE NORMAL EQUATIONS HAS ITS
C                     INITIAL ELEMENT STORED IN A(K), WHERE
C                     K = (M-1)*M/2+1.
C                (B)  THE I-TH DIAGONAL ELEMENT (CONTAINING A CORRECTED
C                     SUM OF SQUARES FOR AN INDEPENDENT VARIABLE) IS
C                     STORED IN A(L), WHERE L = I*(I+1)/2 AND I IS IN
C                     THE SET (1,2,...,M-1).
C
C-----------------------------------------------------------------------
C
      SUBROUTINE RLSUM (AA,MM,IH,M,A,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            MM,IH(1),M,IER
      REAL               AA(1),A(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,II1,IM1,INX,INXJ,ITEMP,J,JI,JJ,JK,K,KK,
     1                   KNXJ,KX,K0,K1,M2
      REAL               TEMP
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
C                                  ERROR - M IS LESS THAN 1
      IF (M .LT. 1) GO TO 65
      M2 = M+M
      DO 5 I = 1,MM
         IH(M2+I) = I
    5 CONTINUE
      DO 20 I = 1,M
         DO 10 J = I,MM
            JJ = J
C                                  ERROR - IH VALUE OUT OF RANGE
            IF (IH(M2+JJ) .EQ. IH(I)) GO TO 15
   10    CONTINUE
         GO TO 65
   15    IH(M+I) = JJ
         ITEMP = IH(M2+I)
         IH(M2+I) = IH(M2+JJ)
         IH(M2+JJ) = ITEMP
   20 CONTINUE
      K = (MM*(MM+1))/2
      DO 25 I = 1,K
         A(I) = AA(I)
   25 CONTINUE
      INX = 0
      DO 60 I = 1,M
         KK = IH(I+M)
         IF (KK .EQ. I) GO TO 55
         K1 = KK+1
         JK = (K1*KK)/2
         KX = JK-KK
         IF (K1 .GT. MM) GO TO 35
C                                  INTERCHANGE COLUMNS
         JI = JK+I
         JK = JK+KK
         DO 30 J = K1,MM
            TEMP = A(JK)
            A(JK) = A(JI)
            A(JI) = TEMP
            JK = JK+J
            JI = JI+J
   30    CONTINUE
   35    K0 = KK-1
         II1 = I+1
         IF (II1 .GT. K0) GO TO 45
         JK = KX+II1
         JI = INX+I+I
         DO 40 J = II1,K0
            TEMP = A(JI)
            A(JI) = A(JK)
            A(JK) = TEMP
            JK = JK+1
            JI = JI+J
   40    CONTINUE
C                                  INTERCHANGE DIAGONAL ELEMENTS
   45    INXJ = INX+I
         KNXJ = KX+KK
         TEMP = A(INXJ)
         A(INXJ) = A(KNXJ)
         A(KNXJ) = TEMP
         IF (I .EQ. 1) GO TO 55
C                                  INTERCHANGE ROWS
         IM1 = I-1
         DO 50 J = 1,IM1
            INXJ = INX+J
            KNXJ = KX+J
            TEMP = A(INXJ)
            A(INXJ) = A(KNXJ)
            A(KNXJ) = TEMP
   50    CONTINUE
   55    INX = INX+I
   60 CONTINUE
      GO TO 9005
   65 IER = 129
 9000 CONTINUE
      CALL UERTST(IER,6HRLSUM )
 9005 RETURN
      END
