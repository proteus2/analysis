C   IMSL ROUTINE NAME   - NMRANK
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - NUMERICAL RANKING
C
C   USAGE               - CALL NMRANK (X,N,EPS,IR,R,RANK,S,T)
C
C   ARGUMENTS    X      - INPUT VECTOR OF LENGTH N CONTAINING THE VALUES
C                           TO BE RANKED
C                N      - INPUT NUMBER OF VALUES TO BE RANKED
C                EPS    - INPUT VALUE USED FOR DETECTING TIES.
C                           X VALUES ARE CONSIDERED EQUAL IF THE
C                           DIFFERENCE BETWEEN THEM IS LESS THAN OR
C                           EQUAL TO EPS.
C                IR     - WORK VECTOR OF LENGTH N
C                R      - WORK VECTOR OF LENGTH N
C                RANK   - OUTPUT VECTOR OF LENGTH N CONTAINING THE RANKS
C                           FOR THE VALUES IN THE X VECTOR
C                S      - OUTPUT STATISTIC. S EQUALS THE SUMMATION OF
C                           TI**2-TI WHERE TI IS THE NUMBER OF ELEMENTS
C                           IN THE I-TH GROUP OF TIED VALUES.
C                T      - OUTPUT STATISTIC. T EQUALS THE SUMMATION OF
C                           TI**3-TI WHERE TI IS THE NUMBER OF ELEMENTS
C                           IN THE I-TH GROUP OF TIED VALUES.
C
C   REQD. IMSL ROUTINES - SINGLE/VSRTR
C                       - DOUBLE/VSRTRD
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS  1.  IF IT IS NOT NECESSARY TO RETAIN THE INPUT VECTOR X,
C                X MAY SHARE THE SAME STORAGE LOCATIONS AS EITHER R OR
C                RANK. R AND RANK CAN NEVER SHARE THE SAME STORAGE.
C            2.  IF AN EPS VALUE OF ZERO IS ENTERED, THEN A TIE IN
C                RANKINGS IMPLIES THAT THE TIED OBSERVATIONS OF X ARE
C                EQUAL TO MACHINE PRECISION.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE NMRANK (X,N,EPS,IR,R,RANK,S,T)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N,IR(1)
      REAL               X(1),EPS,R(1),RANK(1),S,T
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,J,JDON,JJ,JJJ,J2,K,K1,L,N1
      REAL               Y
C                                  FIRST EXECUTABLE STATEMENT
      DO 5 I = 1,N
         R(I) = X(I)
    5 IR(I) = I
C                                  SORT ELEMENTS OF VECTOR R INTO
C                                  ASCENDING SEQUENCE SAVING
C                                  PERMUTATIONS
      CALL VSRTR (R,N,IR)
      S = 0.0
      T = 0.0
      N1 = N-1
      L = 1
   10 DO 30 J = L,N1
         JJ = J
         Y = R(J)
         IF (ABS(Y-R(J+1)) .GT. EPS) GO TO 28
C                                  COUNT THE NUMBER OF TIES
         K = 1
         J2 = J+2
         IF (J2 .GT. N) GO TO 20
         DO 15 I = J2,N
            IF (ABS(Y-R(I)) .GT. EPS) GO TO 20
   15       K = K+1
   20    Y = J+.5*K
      K1 = K+1
      DO 25 I = 1,K1
         JJ = J+I-1
         JJJ = IR(JJ)
   25    RANK(JJJ) = Y
         I = K*(K+1)
         S = S+I
         T = T+I*(K+2)
         GO TO 35
   28 JDON = IR(J)
   30 RANK(JDON) = J
   35 L = JJ+1
      IF (L .LE. N1) GO TO 10
      IF (L .NE. N) GO TO 40
      JDON = IR(N)
      RANK(JDON) = N
   40 CONTINUE
      RETURN
      END
