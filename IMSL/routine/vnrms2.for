C   IMSL ROUTINE NAME   - VNRMS2
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - EUCLIDEAN-NORM OF MATRICES (SYMMETRIC
C                           STORAGE MODE)
C
C   USAGE               - CALL VNRMS2 (A,N,XNRMA)
C
C   ARGUMENTS    A      - SYMMETRIC MATRIX FOR WHICH THE NORM IS TO BE
C                           COMPUTED. A IS STORED IN SYMMETRIC STORAGE
C                           MODE. (INPUT)
C                N      - ORDER OF THE MATRIX A. (INPUT)
C                XNRMA  - EUCLIDEAN-NORM OF MATRIX A. (OUTPUT)
C
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
      SUBROUTINE VNRMS2 (A,N,XNRMA)
C
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N
      REAL               A(1),XNRMA
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            J,K,L,M
      REAL               COL(2)
      REAL               SNRM2
C                                  FIRST EXECUTABLE STATEMENT
      M = (N*(N+1))/2
      COL(1) = SNRM2(M,A,1)
      IF (N.LE.1) GO TO 10
      L = 1
      K = 2
      DO 5 J=2,N
         COL(2) = SNRM2(L,A(K),1)
         COL(1) = SNRM2(2,COL,1)
         K = K+J
         L = L+1
    5 CONTINUE
   10 XNRMA = COL(1)
      RETURN
      END
