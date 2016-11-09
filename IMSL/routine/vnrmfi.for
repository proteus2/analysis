C   IMSL ROUTINE NAME   - VNRMFI
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - INFINITY-NORM MATRICES (FULL STORAGE MODE)
C
C   USAGE               - CALL VNRMFI (A,N,IA,XNRMA)
C
C   ARGUMENTS    A      - MATRIX FOR WHICH THE INFINITY-NORM IS TO BE
C                           COMPUTED. (INPUT)
C                N      - ORDER OF MATRIX A. (INPUT)
C                IA     - ROW DIMENSION OF MATRIX A EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM. (INPUT)
C                XNRMA  - VALUE OF THE INFINITY-NORM. (OUTPUT)
C
C   REQD. IMSL ROUTINES - VABSMF
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE VNRMFI (A,N,IA,XNRMA)
      DIMENSION          A(IA,1)
C                                  FIRST EXECUTABLE STATEMENT
      XNRMA=0.0
         DO 5 I=1,N
         CALL VABSMF(A(I,1),N,IA,VSUM)
         IF(VSUM.GE.XNRMA) XNRMA=VSUM
    5    CONTINUE
      RETURN
      END
