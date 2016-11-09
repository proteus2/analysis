C   IMSL ROUTINE NAME   - OFSCOR
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - COMPUTE A SET OF FACTOR SCORES GIVEN THE
C                           FACTOR SCORE COEFFICIENT MATRIX
C
C   USAGE               - CALL OFSCOR (C,IC,NV,NF,NT,Z,IZ,ZBAR,STD,
C                           FMEAN,SS,WK,IER)
C
C   ARGUMENTS    C      - INPUT NV BY NF FACTOR SCORE COEFFICIENT
C                           MATRIX.
C                IC     - INPUT ROW DIMENSION OF MATRIX C EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM.
C                NV     - INPUT NUMBER OF VARIABLES.
C                NF     - INPUT NUMBER OF FACTORS.
C                NT     - INPUT NUMBER OF CASES OR SUBJECTS IN THE
C                           INPUT DATA MATRIX Z.
C                Z      - ON INPUT, THE NT BY NV MATRIX CONTAINING SOME
C                           OR ALL THE DATA (AFTER ANY TRANSFORMATIONS,
C                           EXCEPTING STANDARDIZATIONS INDICATED IN
C                           INPUT VECTORS ZBAR AND STD).  ON OUTPUT,
C                           THE NT BY NF MATRIX OF FACTOR SCORES.
C                           Z(I,J) CONTAINS THE SCORE OF THE J-TH
C                           FACTOR ON SUBJECT NUMBER I.
C                IZ     - INPUT ROW DIMENSION OF MATRIX Z EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM.
C                ZBAR   - INPUT VECTOR OF LENGTH NV.  ZBAR(I) CONTAINS
C                           THE MEAN OF THE I-TH VARIABLE.
C                STD    - INPUT VECTOR OF LENGTH NV.  STD(I) CONTAINS
C                           THE STANDARD DEVIATION OF THE I-TH VARIABLE
C                           FROM ITS MEAN ZBAR(I).
C                FMEAN  - INPUT CONSTANT.  EACH CALCULATED FACTOR SCORE
C                           IS NORMALIZED TO HAVE FMEAN AS ITS MEAN AND
C                           SS AS ITS STANDARD DEVIATION.  TYPICAL
C                           VALUES FOR FMEAN AND SS ARE 50.0 AND 10.0.
C                SS     - INPUT CONSTANT.  SEE DESCRIPTION OF FMEAN.
C                WK     - WORK VECTOR OF LENGTH NV+NF.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER = 129 INDICATES THAT AT LEAST ONE OF NV,
C                             NF, NT, IZ, OR IC WAS SPECIFIED INCOR-
C                             RECTLY.
C
C   REQD. IMSL ROUTINES - UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE OFSCOR  (C,IC,NV,NF,NT,Z,IZ,ZBAR,STD,FMEAN,SS,WK,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            IC,NV,NF,NT,IZ,IER
      REAL               Z(IZ,NV),ZBAR(NV),STD(NV),WK(1),C(IC,NF)
      REAL               SS,FMEAN
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,J,K
      REAL               SAVE,ZERO
      DOUBLE PRECISION   TEMP
      DATA               ZERO/0.0E0/
C                                  ERROR CHECKS
C                                  FIRST EXECUTABLE STATEMENT
      IF (IZ.GE.NT.AND.IC.GE.NV.AND.NV.GE.NF) GO TO 5
      IER = 129
      GO TO 9000
    5 IER = 0
      DO 10 I=1,NF
         WK(NV+I) = ZERO
   10 CONTINUE
C                                  REPLACE DATA MATRIX WITH SCORES
      DO 30 K=1,NT
         DO 15 I=1,NV
            WK(I) = (Z(K,I)-ZBAR(I))/STD(I)
   15    CONTINUE
         DO 25 J=1,NF
            TEMP = 0.0D0
            DO 20 I=1,NV
               TEMP = TEMP+DBLE(WK(I))*C(I,J)
   20       CONTINUE
            Z(K,J) = TEMP
            WK(NV+J) = WK(NV+J)+TEMP*TEMP
   25    CONTINUE
   30 CONTINUE
      SAVE = DSQRT(1.0D0/NT)
      DO 35 I=1,NF
         WK(I) = SQRT(WK(NV+I))*SAVE
   35 CONTINUE
C                                  SCALE FACTOR SCORES
      DO 45 K=1,NT
         DO 40 J=1,NF
            Z(K,J) = FMEAN+Z(K,J)*SS/WK(J)
   40    CONTINUE
   45 CONTINUE
      GO TO 9005
 9000 CONTINUE
      CALL UERTST (IER,6HOFSCOR)
 9005 RETURN
      END
