C   IMSL ROUTINE NAME   - RLSUBM
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - RETRIEVAL OF A SYMMETRIC SUBMATRIX FROM A
C                           MATRIX STORED IN SYMMETRIC STORAGE MODE
C                           BY RLSTP
C
C   USAGE               - CALL RLSUBM (A,M,IH,S,N)
C
C   ARGUMENTS    A      - M BY M SYMMETRIC MATRIX STORED IN SYMMETRIC
C                           MODE. (INPUT) A IS A VECTOR OF LENGTH
C                           M*(M+1)/2.
C                M      - ORDER OF THE MATRIX A. (INPUT)
C                IH     - VECTOR OF LENGTH M. (INPUT)
C                         IF IH(I)=IH(J)=1 WHERE J AND I = 1,2,...,M,
C                           THE (I,J)-TH ELEMENT OF A WILL BE INCLUDED
C                           IN THE SUBMATRIX S.
C                         OTHERWISE, THE (I,J)-TH ELEMENT OF A WILL NOT
C                           BE IN THE SUBMATRIX S ON OUTPUT.
C                S      - SYMMETRIC SUBMATRIX OF MATRIX A. (OUTPUT)
C                           S IS A VECTOR OR LENGTH N*(N+1)/2.
C                           A AND S MAY SHARE THE SAME STORAGE IF IT IS
C                           NOT NECESSARY TO RETAIN THE ORIGINAL
C                           MATRIX. SEE REMARKS FOR RLSTP.
C                N      - ORDER OF THE SUBMATRIX S. (OUTPUT)
C
C   REQD. IMSL ROUTINES - NONE REQUIRED
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE RLSUBM (A,M,IH,S,N)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            M,IH(M),N
      REAL               A(1),S(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,II,J,JG,JJ,JJJ,J1,K,L,LL
C                                  FIRST EXECUTABLE STATEMENT
      N = 0
      DO 5 J=1,M
         I = J
         IF (IH(J).EQ.1) GO TO 10
    5 CONTINUE
      GO TO 35
   10 J1 = I+1
      II = I
      L = (I*J1)/2
      S(1) = A(L)
      K = 2
      N = 1
      L = 2
      JG = J1
      IF (JG .GT. M) GO TO 35
      DO 30 JJ=JG,M
         IF (IH(JJ).NE.1) GO TO 20
         J1 = (JJ*I)/2
         N = N+1
         LL = 0
         DO 15 JJJ=II,JJ
            LL = LL+1
            IF (IH(JJJ).NE.1) GO TO 15
            S(K) = A(J1+JJJ)
            K = K+1
            IF (LL.EQ.L) GO TO 20
   15    CONTINUE
   20    L = L+1
   25    I = JJ
   30 CONTINUE
   35 RETURN
      END
