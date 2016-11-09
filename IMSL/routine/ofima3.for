C   IMSL ROUTINE NAME   - OFIMA3
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - LEAST SQUARES SOLUTION TO THE MATRIX
C                           EQUATION AT = B.
C
C   USAGE               - CALL OFIMA3 (A,IA,B,IB,NV,NS,NF,T,IT,WK,IER)
C
C   ARGUMENTS    A      - INPUT MATRIX OF DIMENSION NV BY NF OF
C                           COLUMN RANK NF.
C                IA     - INPUT ROW DIMENSION OF MATRIX A EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM.
C                B      - INPUT NV BY NS MATRIX.
C                IB     - INPUT ROW DIMENSION OF MATRIX B EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM.
C                NV     - INPUT NUMBER OF ROWS IN MATRICES A AND B.
C                NS     - INPUT NUMBER OF COLUMNS IN MATRICES T AND B.
C                NF     - INPUT NUMBER OF ROWS IN MATRIX T AND
C                           NUMBER OF COLUMNS IN MATRIX A.
C                T      - OUTPUT MATRIX OF DIMENSION NF BY NS CONTAINING
C                           THE LEAST SQUARES EQUATION SOLUTION.
C                IT     - INPUT ROW DIMENSION OF MATRIX T EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM.
C                WK     - WORK VECTOR OF LENGTH (NF+1)*NF/2.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER = 129 INDICATES THAT THE RANK OF A IS
C                             LESS THAN NF NUMERICALLY.
C                           IER = 130 INDICATES AT LEAST ONE OF IA, IB,
C                             OR IT WAS SPECIFIED INCORRECTLY.
C
C
C   REQD. IMSL ROUTINES - SINGLE/LEQT1P,LUDECP,LUELMP,UERTST,UGETIO,
C                           VIPRFF,VMULFM,VTPROF
C                       - DOUBLE/LEQT1P,LUDECP,LUELMP,UERTST,UGETIO,
C                           VIPRFF,VMULFM,VTPROF,VXADD,VXMUL,VXSTO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE OFIMA3 (A,IA,B,IB,NV,NS,NF,T,IT,WK,IER)
C
      REAL               D1,D2
      REAL               A,B,T,WK
      DIMENSION          A(IA,NF),B(IB,NS),T(IT,NS),WK(1)
C                                  ERROR CHECKS
C                                  FIRST EXECUTABLE STATEMENT
      IF (NV.LE.IA.AND.NV.LE.IB.AND.NF.LE.IT.AND.NF.LE.NV) GO TO 5
      IER = 130
      GO TO 9000
C                                  CALCULATE (A-TRANSPOSE) * A
    5 CALL VTPROF (A,NV,NF,IA,WK)
C                                  CALCULATE (A-TRANSPOSE) * B
      CALL VMULFM (A,B,NV,NF,NS,IA,IB,T,IT,IER)
C                                  SOLVE FOR T MATRIX
      CALL LEQT1P (WK,NS,NF,T,IT,0,D1,D2,IER)
      IF (IER.EQ.0) GO TO 9005
 9000 CONTINUE
      CALL UERTST (IER,6HOFIMA3)
 9005 RETURN
      END
