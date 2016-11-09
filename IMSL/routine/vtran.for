C   IMSL ROUTINE NAME   - VTRAN
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - TRANSPOSE A RECTANGULAR MATRIX
C
C   USAGE               - CALL VTRAN (A,N,M)
C
C   ARGUMENTS    A      - REAL MATRIX OF DIMENSION N BY M. ON INPUT A
C                           CONTAINS THE MATRIX TO BE TRANSPOSED.
C                           ON OUTPUT A CONTAINS THE TRANSPOSED MATRIX.
C                N      - INTEGER NUMBER OF ROWS OF A. (INPUT)
C                M      - INTEGER NUMBER OF COLUMNS OF A. (INPUT)
C
C   REQD. IMSL ROUTINES - VDCPS
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE VTRAN  (A,N,M)
C
      DIMENSION          A(1),IPF(13),IEXP(13),IPWR(13),ITEMP(13)
      REAL               A,TEMP,TEMP1
C                                  FIRST EXECUTABLE STATEMENT
C                                  IF A IS ONE COLUMN OR ONE ROW RETURN
      IF (N .LE. 1 .OR. M .LE. 1) GO TO 9005
C                                  CHECK TO SEE IF A IS SQUARE
      NM = N*M-1
      IF (N .NE. M) GO TO 15
C                                  SQUARE MATRICES ARE HANDLED
C                                    SEPARATELY
      K = 2
      NP1 = N+1
      NM1 = N-1
      DO 10 I=N,NM,N
         L = K+NM1
         DO 5 J=K,I
            TEMP = A(J)
            A(J) = A(L)
            A(L) = TEMP
            L = L+N
    5    CONTINUE
         K = K+NP1
   10 CONTINUE
      GO TO 9005
C                                  A IS NOT SQUARE. TRANSPOSITION OF
C                                    AN N BY M MATRIX AMOUNTS TO
C                                    REPLACING THE ELEMENT AT VECTOR
C                                    POSITION I WITH THE ELEMENT AT
C                                    POSITION N*I (MOD N*M-1). EACH
C                                    SUBCYCLE OF THIS PERMUTATION IS
C                                    COMPLETED IN ORDER.
C                                  DECOMPOSE NM INTO ITS PRIME FACTORS.
   15 CALL VDCPS(NM,NPF,IPF,IEXP,IPWR)
      DO 20 I=1,NPF
         ITEMP(I) = 0
   20 CONTINUE
C                                  GENERATE DIVISORS OF NM LESS THAN
C                                    NM/2
      ID = 1
   25 IF (ID .GE. NM/2) GO TO 9005
      J = NM/ID
      DO 30 I=1,NPF
         IF (ITEMP(I) .NE. IEXP(I)) J = (J/IPF(I))*(IPF(I)-1)
   30 CONTINUE
C                                  THE STARTING POINT OF A SUBCYCLE IS
C                                    DIVISIBLE ONLY BY ID AND MUST NOT
C                                    APPEAR IN ANY OTHER SUBCYCLE
      ISTART = ID
   35 NMIS = NM-ISTART
      IF (ISTART .EQ. ID) GO TO 50
      ISD = ISTART/ID
      DO 40 I=1,NPF
         IF (ITEMP(I) .EQ. IEXP(I)) GO TO 40
         IF (MOD(ISD,IPF(I)) .EQ. 0) GO TO 70
   40 CONTINUE
      IT = ISTART
   45 IT = MOD(N*IT,NM)
      IF (IT .LT. ISTART .OR. IT .GT. NMIS) GO TO 70
      IF (IT .GT. ISTART .AND. IT .LT. NMIS) GO TO 45
   50 TEMP = A(ISTART+1)
      TEMP1 = A(NMIS+1)
      K = ISTART
   55 KK = MOD(N*K,NM)
      L = NM-K
      LL = NM-KK
      J = J-2
C                                  MOVE TWO ELEMENTS. THE SECOND FROM
C                                    THE NEGATIVE SUBCYCLE. CHECK TO
C                                    SEE IF THE SUBCYCLE IS COMPLETE.
      IF (KK .EQ. ISTART) GO TO 60
      IF (LL .EQ. ISTART) GO TO 65
C                                  SUBCYCLE NOT COMPLETE. MOVE ELEMENTS
C                                    AND UPDATE POINTERS
      A(K+1) = A(KK+1)
      A(L+1) = A(LL+1)
      K = KK
      GO TO 55
C                                  SUBCYCLE COMPLETE. RECOMPUTE ID
   60 A(K+1) = TEMP
      A(L+1) = TEMP1
      GO TO 70
   65 A(K+1) = TEMP1
      A(L+1) = TEMP
   70 ISTART = ISTART+ID
      IF (J .GT. 0) GO TO 35
      DO 80 I=1,NPF
         IF (ITEMP(I) .EQ. IEXP(I)) GO TO 75
         ITEMP(I) = ITEMP(I)+1
         ID = ID*IPF(I)
         GO TO 25
   75    ITEMP(I) = 0
         ID = ID/IPWR(I)
   80 CONTINUE
 9005 RETURN
      END
