C   IMSL ROUTINE NAME   - VTPROF
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - TRANSPOSE PRODUCT OF MATRIX (FULL STORAGE
C                           MODE)
C
C   USAGE               - CALL VTPROF (A,L,M,IA,ATA)
C
C   ARGUMENTS    A      - L BY M MATRIX IN FULL STORAGE MODE. (INPUT)
C                L      - MAXIMUM VALUE OF FIRST SUBSCRIPT OF A.
C                           (INPUT)
C                M      - MAXIMUM VALUE OF SECOND SUBSCRIPT OF A.
C                           (INPUT)
C                IA     - ROW DIMENSION OF A EXACTLY AS SPECIFIED IN
C                           THE DIMENSION STATEMENT IN THE CALLING
C                           PROGRAM. (INPUT)
C                ATA    - M BY M MATRIX STORED IN SYMMETRIC MODE AS A
C                           VECTOR. (OUTPUT)
C
C   REQD. IMSL ROUTINES - SINGLE/VIPRFF
C                       - DOUBLE/VIPRFF,VXADD,VXMUL,VXSTO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE VTPROF (A,L,M,IA,ATA)
C
      DIMENSION          A(IA,1),ATA(1)
C                                  FIRST EXECUTABLE STATEMENT
      DO 5 I=1,M
         K=I*(I-1)/2
            DO 5 J=1,I
            K=K+1
    5 CALL VIPRFF(A(1,I),A(1,J),L,1,1,ATA(K))
      RETURN
      END
