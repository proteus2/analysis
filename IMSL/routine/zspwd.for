C   IMSL ROUTINE NAME   - ZSPWD
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
      SUBROUTINE ZSPWD (M,N,A,LDA,V,W)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            M,N,LDA
      REAL               A(LDA,N),V(N),W(N)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,J,NM1,NMJ
      REAL               TEMP1,ONE,TEMP2,TEMP
      DATA               ONE /1.0E0/
C                                  APPLY THE FIRST SET OF GIVENS
C                                  ROTATIONS TO A.
C                                  FIRST EXECUTABLE STATEMENT
      NM1 = N-1
      IF (NM1.LT.1) GO TO 25
      DO 10 NMJ=1,NM1
         J = N-NMJ
         IF (ABS(V(J)).GT.ONE) TEMP1 = ONE/V(J)
         IF (ABS(V(J)).GT.ONE) TEMP2 = SQRT(ONE-TEMP1**2)
         IF (ABS(V(J)).LE.ONE) TEMP2 = V(J)
         IF (ABS(V(J)).LE.ONE) TEMP1 = SQRT(ONE-TEMP2**2)
         DO 5 I=1,M
            TEMP = TEMP1*A(I,J)-TEMP2*A(I,N)
            A(I,N) = TEMP2*A(I,J)+TEMP1*A(I,N)
            A(I,J) = TEMP
    5    CONTINUE
   10 CONTINUE
C                                  APPLY THE SECOND SET OF GIVENS
C                                  ROTATIONS TO A.
      DO 20 J=1,NM1
         IF (ABS(W(J)).GT.ONE) TEMP1 = ONE/W(J)
         IF (ABS(W(J)).GT.ONE) TEMP2 = SQRT(ONE-TEMP1**2)
         IF (ABS(W(J)).LE.ONE) TEMP2 = W(J)
         IF (ABS(W(J)).LE.ONE) TEMP1 = SQRT(ONE-TEMP2**2)
         DO 15 I=1,M
            TEMP = TEMP1*A(I,J)+TEMP2*A(I,N)
            A(I,N) = -TEMP2*A(I,J)+TEMP1*A(I,N)
            A(I,J) = TEMP
   15    CONTINUE
   20 CONTINUE
   25 CONTINUE
      RETURN
      END
