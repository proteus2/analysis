C   IMSL ROUTINE NAME   - VNRMF2
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - EUCLIDEAN-NORM OF MATRICES (FULL STORAGE
C                           MODE)
C
C   USAGE               - CALL VNRMF2 (A,N,IA,XNRMA)
C
C   ARGUMENTS    A      - MATRIX FOR WHICH THE NORM IS TO BE COMPUTED.
C                           (INPUT)
C                N      - ORDER OF MATRIX A. (INPUT)
C                IA     - ROW DIMENSION OF MATRIX A EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM. (INPUT)
C                XNRMA  - EUCLIDEAN-NORM OF A. (OUTPUT)
C
C   REQD. IMSL ROUTINES - SINGLE/VBLA=SNRM2
C                       - DOUBLE/VBLA=DNRM2
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE VNRMF2 (A,N,IA,XNRMA)
C
C                                  SPECIFICATIONS FOR ARGUEMENTS
      INTEGER            IA,N
      REAL               A(IA,N),XNRMA
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            J
      REAL               COL(2)
      REAL               SNRM2
C                                  FIRST EXECUTABLE STATEMENT
      COL(1) = SNRM2(N,A(1,1),1)
      IF (N.LE.1) GO TO 10
      DO 5 J=2,N
         COL(2) = SNRM2(N,A(1,J),1)
         COL(1) = SNRM2(2,COL,1)
    5 CONTINUE
   10 XNRMA = COL(1)
      RETURN
      END
