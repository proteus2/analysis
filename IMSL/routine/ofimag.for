C   IMSL ROUTINE NAME   - OFIMAG
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - COMPUTE AN UNROTATED FACTOR LOADING MATRIX
C                           ACCORDING TO AN IMAGE MODEL
C
C   USAGE               - CALL OFIMAG (R,NV,NF,A,IA,CI,IC,Y,S,G,IS,
C                           WK,IER)
C
C   ARGUMENTS    R      - INPUT VECTOR OF LENGTH (NV+1)*NV/2 CONTAINING
C                           THE NV BY NV CORRELATION MATRIX IN
C                           SYMMETRIC STORAGE MODE. (SEE REMARKS).
C                NV     - INPUT NUMBER OF VARIABLES (AFTER DATA
C                           TRANSFORMATION, IF ANY WAS PERFORMED).
C                NF     - ON INPUT, NF POSITIVE IMPLIES THAT NF FACTORS
C                           ARE TO BE RETAINED.  OTHERWISE, ON OUTPUT,
C                           NF WILL BE DETERMINED BY THE KAISER-GUTTMAN
C                           CRITERION OF GREATER THAN UNITY EIGENVALUES.
C                           NF MUST BE LESS THAN OR EQUAL TO NV.
C                A      - OUTPUT NV BY NV MATRIX CONTAINING THE
C                           IMAGE FACTOR PATTERN IN THE FIRST NF
C                           COLUMNS. THE REMAINING LOCATIONS ARE WORK
C                           STORAGE.
C                IA     - INPUT ROW DIMENSION OF MATRIX A EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM.
C                CI     - OUTPUT NV BY NF MATRIX CONTAINING THE
C                           IMAGE SCORE COEFFICIENT MATRIX.
C                IC     - INPUT ROW DIMENSION OF MATRIX CI EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM.
C                Y      - OUTPUT VECTOR OF LENGTH NV CONTAINING THE
C                           EIGENVALUES OF S IN ASCENDING ORDER.
C                S      - OUTPUT VECTOR OF LENGTH (NV+1)*NV/2 CONTAINING
C                           THE (D INVERSE)*R*(D INVERSE) MATRIX IN
C                           SYMMETRIC STORAGE MODE.
C                           (SEE ALGORITHM SECTION).
C                G      - OUTPUT SCALING VECTOR OF LENGTH NV, CONTAINING
C                           THE ELEMENTS OF THE DIAGONAL MATRIX D.
C                           (SEE ALGORITHM SECTION).
C                IS     - INTEGER WORKING VECTOR OF LENGTH NV.
C                WK     - WORK VECTOR OF LENGTH NV+((NV+1)*NV/2).
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER = 129 INDICATES THAT AT LEAST ONE OF NV,
C                             NF, IA, OR IC WAS SPECIFIED INCORRECTLY.
C                           IER = 130 INDICATES THAT ALL THE EIGENVALUES
C                             OF R COULD NOT BE CALCULATED.
C                           IER = 131 INDICATES THAT R WAS NOT POSITIVE
C                             DEFINITE.
C
C   REQD. IMSL ROUTINES - SINGLE/EHOBKS,EHOUSS,EIGRS,EQRT2S,LINV1P,
C                           LUDECP,LUELMP,UERTST,UGETIO,VSRTU
C                       - DOUBLE/EHOBKS,EHOUSS,EIGRS,EQRT2S,LINV1P,
C                           LUDECP,LUELMP,UERTST,UGETIO,VSRTUD
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS      R IS ASSUMED TO BE A POSITIVE DEFINITE MATRIX. THAT IS,
C                ALL NV EIGENVALUES OF R ARE POSITIVE. OFTEN THE DIA-
C                GONAL ELEMENTS OF R ARE REPLACED BY THE SQUARED MULTI-
C                PLE CORRELATIONS OF THE VARIABLES PRIOR TO CALLING
C                OFIMAG.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE OFIMAG  (R,NV,NF,A,IA,CI,IC,Y,S,G,IS,WK,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            NV,NF,IA,IC,IS(2),IER
      REAL               R(1),G(NV),CI(IC,1),A(IA,1),Y(NV),S(1),WK(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            NVP1,NVV,IR,I,J
      REAL               DD,ONE,VVV
      DATA               ONE/1.0E0/
C                                  ERROR CHECKS
C                                  FIRST EXECUTABLE STATEMENT
      IF (NF.LE.NV.AND.NV.LE.IA.AND.NV.LE.IC) GO TO 5
      IER = 129
      GO TO 9000
    5 NVP1 = NV+1
      NVV = ((NV+1)*NV)/2
      DO 10 I=1,NVV
         WK(I) = R(I)
   10 CONTINUE
C                                  INVERT CORRELATION MATRIX
      CALL LINV1P (WK,NV,S,0,DD,VVV,IER)
      IF (IER.EQ.0) GO TO 15
      IER = 131
      GO TO 9000
C                                  FORM (S-INVERSE)*R*(S-INVERSE) MATRIX
   15 IR = 0
      DO 25 I=1,NV
         IS(I) = NVP1-I
         VVV = SQRT(S(IR+I))
         WK(NVV+I) = VVV
         G(I) = ONE/VVV
         DO 20 J=1,I
            IR = IR+1
            S(IR) = R(IR)*VVV*WK(NVV+J)
            WK(IR) = S(IR)
   20    CONTINUE
   25 CONTINUE
C                                  CALCULATE EIGENPAIR OF S
      CALL EIGRS (WK,NV,1,Y,A,IA,WK(NVV+1),IER)
      IF (IER.EQ.0) GO TO 30
      IER = 130
      GO TO 9000
   30 CALL VSRTU  (A,IA,NV,NV,0,IS,WK)
C                                  DETERMINE NF
      IF (NF.GT.0) GO TO 40
      NF = 0
      DO 35 I=1,NV
         IF (Y(NVP1-I).LE.ONE) GO TO 40
         NF = NF+1
   35 CONTINUE
C                                  FORM UNROTATED IMAGE SCORE
C                                  COEFFICIENT MATRIX CI, AND SCALE
C                                  UNROTATED FACTOR LOADING MATRIX
   40 DO 50 J=1,NF
         VVV = ONE/SQRT(Y(NVP1-J))
         DD = ABS(ONE-Y(NVP1-J))*VVV
         DO 45 I=1,NV
            CI(I,J) = A(I,J)*VVV/G(I)
            A(I,J) = A(I,J)*G(I)*DD
   45    CONTINUE
   50 CONTINUE
      GO TO 9005
 9000 CONTINUE
      CALL UERTST (IER,6HOFIMAG)
 9005 RETURN
      END
