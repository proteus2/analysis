C   IMSL ROUTINE NAME   - VMULFM
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1980
C
C   PURPOSE             - MATRIX MULTIPLICATION OF THE TRANSPOSE OF
C                           MATRIX A BY MATRIX B (FULL STORAGE MODE)
C
C   USAGE               - CALL VMULFM (A,B,L,M,N,IA,IB,C,IC,IER)
C
C   ARGUMENTS    A      - L BY M MATRIX STORED IN FULL STORAGE MODE.
C                           (INPUT)
C                B      - L BY N MATRIX STORED IN FULL STORAGE MODE.
C                           (INPUT)
C                L      - NUMBER OF ROWS IN A AND B. (INPUT)
C                M      - NUMBER OF COLUMNS IN MATRIX A AND NUMBER
C                           OF ROWS IN MATRIX C. (INPUT)
C                N      - NUMBER OF COLUMNS IN B AND C. (INPUT)
C                IA     - ROW DIMENSION OF A EXACTLY AS SPECIFIED IN
C                           IN THE DIMENSION STATEMENT IN THE CALLING
C                           PROGRAM. (INPUT)
C                IB     - ROW DIMENSION OF B EXACTLY AS SPECIFIED IN
C                           THE DIMENSION STATEMENT OF THE CALLING
C                           PROGRAM. (INPUT)
C                C      - M BY N MATRIX CONTAINING THE PRODUCT
C                           C = (A-TRANSPOSE) * B. (OUTPUT)
C                IC     - ROW DIMENSION OF C EXACTLY AS SPECIFIED IN
C                           THE DIMENSION STATEMENT IN THE CALLING
C                           PROGRAM. (INPUT)
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         IER=0 IMPLIES NO ERROR.
C                         TERMINAL ERROR
C                           IER=129 INDICATES A,B, OR C WAS DIMENSIONED
C                             INCORRECTLY.
C
C   REQD. IMSL ROUTINES - UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE VMULFM (A,B,L,M,N,IA,IB,C,IC,IER)
C
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            L,M,N,IA,IB,IC,IER
      REAL               A(IA,1),B(IB,1),C(IC,1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      DOUBLE PRECISION   TEMP
C                                  FIRST EXECUTABLE STATEMENT
      IF (IA.GE.L .AND. IB.GE.L .AND. IC.GE.M) GO TO 5
C                                  TERMINAL ERROR
      IER = 129
      GO TO 9000
C                                  ROW INDICATOR
    5 IER = 0
      DO 20 I = 1,M
C                                  COLUMN INDICATOR
         DO 15 J = 1,N
            TEMP = 0.0D0
C                                  VECTOR DOT PRODUCT
            DO 10 K = 1,L
               TEMP = TEMP + DBLE(A(K,I))*DBLE(B(K,J))
   10       CONTINUE
            C(I,J) = TEMP
   15    CONTINUE
   20 CONTINUE
      GO TO 9005
 9000 CONTINUE
      CALL UERTST (IER,6HVMULFM)
 9005 RETURN
      END
