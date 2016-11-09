C   IMSL ROUTINE NAME   - VSRTU
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - INTERCHANGE THE ROWS OR COLUMNS OF A MATRIX
C                           USING A PERMUTATION VECTOR SUCH AS THE ONE
C                           OBTAINED FROM IMSL ROUTINES VSRTP OR
C                           VSRTR
C
C   USAGE               - CALL VSRTU (Z,IZ,N,M,IND,IR,WK)
C
C   ARGUMENTS    Z      - INPUT MATRIX OF DIMENSION N BY M TO BE
C                           INTERCHANGED.  ON OUTPUT, Z CONTAINS
C                           THE INTERCHANGED MATRIX.
C                IZ     - ROW DIMENSION OF MATRIX Z EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM. (INPUT)
C                N      - NUMBER OF ROWS IN Z.  (INPUT)
C                M      - NUMBER OF COLUMNS IN Z.  (INPUT)
C                IND    - INPUT OPTION PARAMETER.
C                           IF IND IS GREATER THAN ZERO, THE ROWS OF Z
C                             WILL BE INTERCHANGED ACCORDING TO THE
C                             INFORMATION IN VECTOR IR.
C                           OTHERWISE, THE COLUMNS OF Z WILL BE
C                             INTERCHANGED ACCORDING TO THE
C                             INFORMATION IN VECTOR IR.
C                IR     - INPUT INTEGER PERMUTATION VECTOR OF LENGTH
C                           N, IF IND IS POSITIVE, AND OF LENGTH M
C                           OTHERWISE.  IR CONTAINS THE FIRST N OR M
C                           POSITIVE INTEGERS.  SEE PROGRAMMING NOTES.
C                           IR IS DESTROYED ON OUTPUT.
C                WK     - WORK VECTOR OF LENGTH M IF IND IS POSITIVE
C                           AND OF LENGTH N OTHERWISE.
C
C   REQD. IMSL ROUTINES - NONE REQUIRED
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE VSRTU (Z,IZ,N,M,IND,IR,WK)
C
      DIMENSION          Z(IZ,1),WK(1),IR(1)
C                                  FIRST EXECUTABLE STATEMENT
      IPTR = 1
      IF (IND .GT. 0) GO TO 45
C                                  SORT Z BY COLUMNS
C                                  CHECK IF ALL COLUMNS ARE SORTED
    5 IF (IPTR .GE. M) GO TO 85
C                                  CHECK IF COLUMN IPTR HAS BEEN SORTED
      IF (IR(IPTR) .GT. 0) GO TO 15
   10 IPTR = IPTR + 1
      GO TO 5
C                                  CHECK IF COLUMN IPTR NEED BE MOVED
   15 IF (IR(IPTR) .EQ. IPTR) GO TO 10
      K = IPTR
C                                  STORE COLUMN IPTR IN TEMPORARY VECTOR
      DO 20 I = 1,N
         WK(I) = Z(I,K)
   20 CONTINUE
   25 L = IR(K)
C                                  CHECK IF TEMPORARY VECTOR NEEDED HERE
      IF (L .EQ. IPTR) GO TO 35
C                                  INSERT COLUMN L INTO COLUMN K
      DO 30 I = 1,N
         Z(I,K) = Z(I,L)
   30 CONTINUE
C                                  MARK COLUMN K AS ALREADY SORTED
      IR(K) = 0
      K = L
      GO TO 25
C                                  INSERT TEMPORARY VECTOR IN COLUMN K
   35 DO 40 I = 1,N
         Z(I,K) = WK(I)
   40 CONTINUE
      IR(K) = 0
      GO TO 10
C                                  SORT Z BY ROWS
   45 IF (IPTR .GE. N) GO TO 85
      IF (IR(IPTR) .GT. 0) GO TO 55
   50 IPTR = IPTR + 1
      GO TO 45
   55 IF (IR(IPTR) .EQ. IPTR) GO TO 50
      K = IPTR
      DO 60 I = 1,M
         WK(I) = Z(K,I)
   60 CONTINUE
   65 L = IR(K)
      IF (L .EQ. IPTR) GO TO 75
      DO 70 I = 1,M
         Z(K,I) = Z(L,I)
   70 CONTINUE
      IR(K) = 0
      K = L
      GO TO 65
   75 DO 80 I = 1,M
         Z(K,I) = WK(I)
   80 CONTINUE
      IR(K) = 0
      GO TO 50
   85 RETURN
      END
