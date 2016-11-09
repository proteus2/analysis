C   IMSL ROUTINE NAME   - OFRESI
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - COMMUNALITIES AND NORMALIZED FACTOR
C                           RESIDUAL CORRELATION MATRIX CALCULATION
C
C   USAGE               - CALL OFRESI (R,NV,NF,A,IA,Y,S,WK)
C
C   ARGUMENTS    R      - INPUT VECTOR OF LENGTH (NV+1)*NV/2 CONTAINING
C                           THE NV BY NV CORRELATION MATRIX IN
C                           SYMMETRIC STORAGE MODE.
C                NV     - INPUT NUMBER OF VARIABLES.
C                NF     - INPUT NUMBER OF FACTORS
C                A      - INPUT NV BY NF FACTOR LOADING MATRIX.
C                IA     - INPUT ROW DIMENSION OF MATRIX A EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM.
C                Y      - OUTPUT VECTOR OF LENGTH NV CONTAINING THE
C                           COMMUNALITIES OF THE VARIABLES.
C                S      - OUTPUT VECTOR OF LENGTH (NV+1)*NV/2 CONTAINING
C                           THE NV BY NV NORMALIZED RESIDUAL CORRELATION
C                           MATRIX IN SYMMETRIC STORAGE MODE.
C                WK     - WORK VECTOR OF LENGTH NV.
C
C   REQD. IMSL ROUTINES - NONE REQUIRED
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
      SUBROUTINE OFRESI  (R,NV,NF,A,IA,Y,S,WK)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            NV,NF,IA
      REAL               A(IA,NF),R(1),S(1),Y(NV),WK(NV)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            IR,I,J,K,NVV
      REAL               DD,ONE,ZERO
      DOUBLE PRECISION   TEMP
      DATA               ZERO,ONE/0.0E0,1.0E0/
C                                  CALCULATE REPRODUCED CORRELATION
C                                  MATRIX, COMMUNALITIES ON DIAGONAL
C                                  FIRST EXECUTABLE STATEMENT
      IR = 0
      DO 15 I=1,NV
         DO 10 J=1,I
            IR = IR+1
            TEMP = 0.0D0
            DO 5 K=1,NF
               TEMP = TEMP+DBLE(A(I,K))*A(J,K)
    5       CONTINUE
            S(IR) = R(IR)-TEMP
   10    CONTINUE
         Y(I) = TEMP
   15 CONTINUE
      NVV = ((NV+1)*NV)/2
C                                  SCALE RESIDUAL CORRELATIONS
      IR = 0
      DO 25 I=1,NV
         WK(I) = S(IR+I)
         IF (WK(I).LE.ZERO) WK(I) = ONE
         WK(I) = SQRT(WK(I))
         DD = WK(I)
         DO 20 J=1,I
            IR = IR+1
            S(IR) = S(IR)/(DD*WK(J))
   20    CONTINUE
   25 CONTINUE
      RETURN
      END
