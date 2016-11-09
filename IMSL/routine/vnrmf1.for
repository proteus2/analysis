C   IMSL ROUTINE NAME   - VNRMF1
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - 1-NORM OF MATRICES (FULL STORAGE MODE)
C
C   USAGE               - CALL VNRMF1 (A,N,IA,XNRMA)
C
C   ARGUMENTS    A      - MATRIX (IN FULL STORAGE MODE) FOR WHICH THE
C                           NORM IS TO BE COMPUTED. (INPUT)
C                N      - ORDER OF MATRIX A. (INPUT)
C                IA     - ROW DIMENSION OF MATRIX A EXACTLY AS
C                           SPECIFIED IN  DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM. (INPUT)
C                XNRMA  - VALUE OF THE 1-NORM OF A. (OUTPUT)
C                           XNRMA = MAX OVER J OF (SUM FROM I=1 TO N OF
C                           ABS(A(I,J)).
C
C   REQD. IMSL ROUTINES - VABSMF
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE VNRMF1 (A,N,IA,XNRMA)
C
      DIMENSION          A(IA,1)
C                                  FIRST EXECUTABLE STATEMENT
      XNRMA=0.0
         DO 5 J=1,N
         CALL VABSMF(A(1,J),N,1,VSUM)
         IF(VSUM.GE.XNRMA) XNRMA=VSUM
    5    CONTINUE
      RETURN
      END
