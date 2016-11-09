C   IMSL ROUTINE NAME   - VTPROS
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - TRANSPOSE PRODUCT OF A MATRIX (SYMMETRIC
C                           STORAGE MODE)
C
C   USAGE               - CALL VTPROS (A,N,ATA)
C
C   ARGUMENTS    A      - N BY N MATRIX STORED IN SYMMETRIC MODE AS
C                           A VECTOR. (INPUT)
C                N      - ORDER OF MATRIX A. (INPUT)
C                ATA    - N BY N MATRIX STORED IN SYMMETRIC MODE AS
C                           A VECTOR. (OUTPUT)
C
C   REQD. IMSL ROUTINES - SINGLE/VIPRSS
C                       - DOUBLE/VIPRSS,VXADD,VXMUL,VXSTO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE VTPROS (A,N,ATA)
C
      DIMENSION          A(1),ATA(1)
C                                  FIRST EXECUTABLE STATEMENT
      DO 5 I=1,N
         K=(I*(I-1))/2
         DO 5 J=1,I
            K=K+1
    5 CALL VIPRSS(A,A,N,I,J,ATA(K))
      RETURN
      END
