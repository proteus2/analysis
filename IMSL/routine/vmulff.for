C   IMSL ROUTINE NAME   - VMULFF
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - NOVEMBER 1, 1979
C
C   PURPOSE             - MATRIX MULTIPLICATION (FULL STORAGE MODE)
C
C   USAGE               - CALL VMULFF (A,B,L,M,N,IA,IB,C,IC,IER)
C
C   ARGUMENTS    A      - L BY M MATRIX STORED IN FULL STORAGE MODE.
C                           (INPUT)
C                B      - M BY N MATRIX STORED IN FULL STORAGE MODE.
C                           (INPUT)
C                L      - NUMBER OF ROWS IN A. (INPUT)
C                M      - NUMBER OF COLUMNS IN A (SAME AS NUMBER OF
C                           ROWS IN B). (INPUT)
C                N      - NUMBER OF COLUMNS IN B. (INPUT)
C                IA     - ROW DIMENSION OF MATRIX A EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM. (INPUT)
C                IB     - ROW DIMENSION OF MATRIX B EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM. (INPUT)
C                C      - L BY N MATRIX CONTAINING THE PRODUCT
C                           C = A*B. (OUTPUT)
C                IC     - ROW DIMENSION OF MATRIX C EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM. (INPUT)
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER=129 INDICATES A,B,OR C WAS DIMENSIONED
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
      SUBROUTINE VMULFF (A,B,L,M,N,IA,IB,C,IC,IER)
C
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            L,M,N,IA,IB,IC,IER
      REAL               A(IA,M),B(IB,N),C(IC,N)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      DOUBLE PRECISION   TEMP
C                                  FIRST EXECUTABLE STATEMENT
      IF (IA .GE. L .AND. IB .GE. M .AND. IC .GE. L) GO TO 5
C                                  TERMINAL ERROR
      IER=129
      GO TO 9000
C                                  ROW INDICATOR
    5 IER = 0
      DO 15 I=1,L
C                                  COLUMN INDICATOR
         DO 15 J=1,N
            TEMP=0.0D0
C                                  VECTOR DOT PRODUCT
            DO 10 K=1,M
               TEMP=DBLE(A(I,K))*DBLE(B(K,J))+TEMP
   10       CONTINUE
            C(I,J)=TEMP
   15 CONTINUE
      GO TO 9005
 9000 CONTINUE
      CALL UERTST (IER,6HVMULFF)
 9005 RETURN
      END
