C   IMSL ROUTINE NAME   - VUAFQ
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1981
C
C   PURPOSE             - MATRIX ADDITION (FULL + BAND SYMMETRIC
C                           MATRICES)
C
C   USAGE               - CALL VUAFQ (A,N,IA,B,NC,IB,C,IC)
C
C   ARGUMENTS    A      - N BY N MATRIX IN FULL STORAGE MODE. (INPUT)
C                N      - ORDER OF MATRICES A,B, AND C. (INPUT)
C                IA     - ROW DIMENSION OF MATRIX A EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM. (INPUT)
C                B      - N BY N BAND SYMMETRIC MATRIX.  B IS STORED IN
C                           BAND SYMMETRIC STORAGE MODE AND THERE-
C                           FORE HAS DIMENSION N BY (NC + 1). (INPUT)
C                NC     - NUMBER OF UPPER OR LOWER CO-DIAGONALS IN
C                           MATRIX B. (INPUT)
C                IB     - ROW DIMENSION OF MATRIX B EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM. (INPUT)
C                C      - N BY N MATRIX CONTAINING THE SUM C = A+B.
C                           MATRIX C IS STORED IN FULL STORAGE MODE.
C                           (OUTPUT)
C                IC     - ROW DIMENSION OF MATRIX C EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM. (INPUT)
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
      SUBROUTINE VUAFQ  (A,N,IA,B,NC,IB,C,IC)
C
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N,IA,NC,IB,IC
      REAL               A(IA,1),B(IB,1),C(IC,1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,J,K,L,NC1,NC2
C                                  FIRST EXECUTABLE STATEMENT
      NC1 = NC+1
      NC2 = NC+2
      IF (NC2 .GT. N) GO TO 15
C                                  SAVE THE ELEMENTS OF A WHICH FALL
C                                  OUTSIDE THE BAND WIDTH OF B
      L = 1
      DO 10 J = NC2,N
         DO 5 I = 1,L
            C(I,J) = A(I,J)
            C(J,I) = A(J,I)
    5    CONTINUE
         L = L+1
   10 CONTINUE
C                                  ADD THE ELEMENTS OF A AND B
C                                  WITHIN THE BAND WIDTH OF B
   15 DO 25 I = 1,N
         K = MAX0(1,I-NC)
         L = MAX0(1,NC2-I)
         DO 20 J = L,NC1
            C(I,K) = A(I,K)+B(I,J)
            IF (J .NE. NC1) C(K,I) = A(K,I)+B(I,J)
            K = K+1
   20    CONTINUE
   25 CONTINUE
      RETURN
      END
