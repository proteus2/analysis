C   IMSL ROUTINE NAME   - VUASQ
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1981
C
C   PURPOSE             - MATRIX ADDITION (SYMMETRIC + BAND SYMMETRIC
C                           MATRICES)
C
C   USAGE               - CALL VUASQ (A,N,B,NC,IB,C)
C
C   ARGUMENTS    A      - INPUT VECTOR OF LENGTH N(N+1)/2 CONTAINING
C                           AN N BY N MATRIX IN SYMMETRIC STORAGE MODE.
C                N      - ORDER OF THE MATRICES A, B, AND C. (INPUT)
C                B      - N BY N SYMMETRIC BAND MATRIX STORED IN BAND
C                           SYMMETRIC STORAGE MODE.  B IS STORED AS A
C                           MATRIX WITH DIMENSION N BY (NC + 1). (INPUT)
C                NC     - NUMBER OF UPPER OR LOWER CO-DIAGONALS IN B.
C                           (INPUT)
C                IB     - ROW DIMENSION OF B EXACTLY AS SPECIFIED IN
C                           THE DIMENSION STATEMENT IN THE CALLING
C                           PROGRAM. (INPUT)
C                C      - VECTOR OF LENGTH N(N+1)/2 CONTAINING AN
C                           N BY N MATRIX IN SYMMETRIC STORAGE MODE.
C                           C = A + B. (OUTPUT)
C
C   REQD. IMSL ROUTINES - NONE REQUIRED
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS      MATRIX SUBTRACTION MAY BE DONE VIA THIS ROUTINE IF
C                (PRIOR TO ENTRY) THE USER MANIPULATES THE SIGNS OF THE
C                MATRICES TO GIVE THE DESIRED RESULT.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE VUASQ  (A,N,B,NC,IB,C)
C
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N,NC,IB
      REAL               A(1),B(IB,1),C(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            II,JSAV,J,K,LL,L,M,NCP1
C                                  FIRST EXECUTABLE STATEMENT
      K = (N*(N + 1))/2
      NCP1 = NC + 1
      II = 1
      J = NCP1
      JSAV = J
      L = 1
      M = 1
C                                  ADD THE ELEMENTS WHICH FALL
C                                  WITHIN THE BAND WIDTH OF B.
    5 C(L) = A(L) + B(II,J)
      IF (L .EQ. K) GO TO 20
      L = L + 1
      IF (J .LT. NCP1) GO TO 15
      II = II + 1
      J = MAX0(1,JSAV-1)
      JSAV = J
      IF (II .LE. NCP1) GO TO 5
      DO 10 LL = 1,M
C                                  MOVE ELEMENTS OUTSIDE THE BAND
C                                  WIDTH OF B, FROM A TO C.
         C(L) = A(L)
         L = L + 1
   10 CONTINUE
      M = M + 1
      GO TO 5
   15 J = J + 1
      GO TO 5
   20 RETURN
      END
