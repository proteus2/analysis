C   IMSL ROUTINE NAME   - ZSPWE
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1982
C
C   PURPOSE             - NUCLEUS CALLED ONLY BY IMSL SUBROUTINE ZSPOW
C
C   REQD. IMSL ROUTINES - NONE REQUIRED
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE ZSPWE (M,N,S,LS,U,V,W,SING)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            M,N,LS
      REAL               S(LS),U(M),V(N),W(M)
      LOGICAL            SING
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,JJ,J,L,NM1,NMJ
      REAL               TEMP1,TEMP2,GIANT,ONE,P25,P5,TEMP3,SPMPAR,
     *                   TEMP4,TAU,TEMP,ZERO
      DATA               GIANT /Z7FFFFFFF/
      DATA               ONE,P5,P25,ZERO /1.0E0,5.0E-1,2.5E-1,0.0E0/
C                                  INITIALIZE THE DIAGONAL ELEMENT
C                                  POINTER.
C                                  FIRST EXECUTABLE STATEMENT
      JJ = (N*(2*M-N+1))/2-(M-N)
C                                  MOVE THE NONTRIVIAL PART OF THE LAST
C                                  COLUMN OF S INTO W.
      L = JJ
      DO 5 I=N,M
         W(I) = S(L)
         L = L+1
    5 CONTINUE
C                                  ROTATE THE VECTOR V INTO A MULTIPLE
C                                  OF THE N-TH UNIT VECTOR IN SUCH A WAY
C                                  THAT A SPIKE IS INTRODUCED INTO W.
      NM1 = N-1
      IF (NM1.LT.1) GO TO 35
      DO 30 NMJ=1,NM1
         J = N-NMJ
         JJ = JJ-(M-J+1)
         W(J) = ZERO
         IF (V(J).EQ.ZERO) GO TO 25
C                                  DETERMINE A GIVENS ROTATION WHICH
C                                  ELIMINATES THE J-TH ELEMENT OF V.
         IF (ABS(V(N)).GE.ABS(V(J))) GO TO 10
         TEMP2 = V(N)/V(J)
         TEMP3 = P5/SQRT(P25+P25*TEMP2**2)
         TEMP1 = TEMP3*TEMP2
         TAU = ONE
         IF (ABS(TEMP1)*GIANT.GT.ONE) TAU = ONE/TEMP1
         GO TO 15
   10    CONTINUE
         TEMP4 = V(J)/V(N)
         TEMP1 = P5/SQRT(P25+P25*TEMP4**2)
         TEMP3 = TEMP1*TEMP4
         TAU = TEMP3
   15    CONTINUE
C                                  APPLY THE TRANSFORMATION TO V AND
C                                  STORE THE INFORMATION NECESSARY TO
C                                  RECOVER THE GIVENS ROTATION.
         V(N) = TEMP3*V(J)+TEMP1*V(N)
         V(J) = TAU
C                                  APPLY THE TRANSFORMATION TO S AND
C                                  EXTEND THE SPIKE IN W.
         L = JJ
         DO 20 I=J,M
            TEMP = TEMP1*S(L)-TEMP3*W(I)
            W(I) = TEMP3*S(L)+TEMP1*W(I)
            S(L) = TEMP
            L = L+1
   20    CONTINUE
   25    CONTINUE
   30 CONTINUE
   35 CONTINUE
C                                  ADD THE SPIKE FROM THE RANK 1 UPDATE
C                                  TO W.
      DO 40 I=1,M
         W(I) = W(I)+V(N)*U(I)
   40 CONTINUE
C                                  ELIMINATE THE SPIKE.
      SING = .FALSE.
      IF (NM1.LT.1) GO TO 70
      DO 65 J=1,NM1
         IF (W(J).EQ.ZERO) GO TO 60
C                                  DETERMINE A GIVENS ROTATION WHICH
C                                  ELIMINATES THE J-TH ELEMENT OF THE
C                                  SPIKE.
         IF (ABS(S(JJ)).GE.ABS(W(J))) GO TO 45
         TEMP2 = S(JJ)/W(J)
         TEMP3 = P5/SQRT(P25+P25*TEMP2**2)
         TEMP1 = TEMP3*TEMP2
         TAU = ONE
         IF (ABS(TEMP1)*GIANT.GT.ONE) TAU = ONE/TEMP1
         GO TO 50
   45    CONTINUE
         TEMP4 = W(J)/S(JJ)
         TEMP1 = P5/SQRT(P25+P25*TEMP4**2)
         TEMP3 = TEMP1*TEMP4
         TAU = TEMP3
   50    CONTINUE
C                                  APPLY THE TRANSFORMATION TO S AND
C                                  REDUCE THE SPIKE IN W.
         L = JJ
         DO 55 I=J,M
            TEMP = TEMP1*S(L)+TEMP3*W(I)
            W(I) = -TEMP3*S(L)+TEMP1*W(I)
            S(L) = TEMP
            L = L+1
   55    CONTINUE
C                                  STORE THE INFORMATION NECESSARY TO
C                                  RECOVER THE GIVENS ROTATION.
         W(J) = TAU
   60    CONTINUE
C                                  TEST FOR ZERO DIAGONAL ELEMENTS IN
C                                  THE OUTPUT S.
         IF (S(JJ).EQ.ZERO) SING = .TRUE.
         JJ = JJ+(M-J+1)
   65 CONTINUE
   70 CONTINUE
C                                  MOVE W BACK INTO THE LAST COLUMN OF
C                                  THE OUTPUT S.
      L = JJ
      DO 75 I=N,M
         S(L) = W(I)
         L = L+1
   75 CONTINUE
      IF (S(JJ).EQ.ZERO) SING = .TRUE.
      RETURN
      END
