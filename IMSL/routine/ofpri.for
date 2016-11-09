C   IMSL ROUTINE NAME   - OFPRI
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JULY 1, 1983
C
C   PURPOSE             - COMPUTE AN UNROTATED FACTOR LOADING MATRIX
C                           ACCORDING TO A PRINCIPAL COMPONENT MODEL
C
C   USAGE               - CALL OFPRI (R,NV,NF,CRIT,A,IA,E,Y,S,G,IS,
C                           WK,IER)
C
C   ARGUMENTS    R      - INPUT VECTOR OF LENGTH (NV+1)*NV/2 CONTAINING
C                           THE NV BY NV CORRELATION MATRIX IN
C                           SYMMETRIC STORAGE MODE.
C                NV     - INPUT NUMBER OF VARIABLES (AFTER DATA
C                           TRANSFORMATION, IF ANY WAS PERFORMED).
C                NF     - INPUT/OUTPUT VALUE SPECIFYING THE NUMBER OF
C                           FACTORS TO BE RETAINED. NF POSITIVE IMPLIES
C                           THAT NF FACTORS ARE TO BE RETAINED. OTHER-
C                           WISE, ON OUTPUT, NF WILL BE DETERMINED BY
C                           THE NUMBER OF EIGENVALUES GREATER THAN
C                           INPUT CRIT. NF MUST BE LESS THAN OR EQUAL
C                           TO NV.
C                CRIT   - INPUT FACTOR RETENTION CRITERION. IF NF IS
C                           NONPOSITIVE, THE NUMBER OF FACTORS TO BE
C                           RETAINED IS DETERMINED BY THE NUMBER OF
C                           EIGENVALUES GREATER THAN CRIT. IN THIS
C                           CASE, CRIT MUST BE GREATER THAN OR EQUAL
C                           TO ZERO. IF THE KAISER-GUTTMAN CRITERION
C                           IS DESIRED, SPECIFY CRIT=1.0. IF NF IS
C                           POSITIVE, CRIT IS NOT USED.
C                A      - OUTPUT NV BY NV MATRIX CONTAINING THE
C                           UNROTATED FACTOR LOADING MATRIX IN THE
C                           FIRST NF COLUMNS. THE REMAINING LOCATIONS
C                           ARE WORK STORAGE.
C                IA     - INPUT ROW DIMENSION OF MATRIX A EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM.
C                E      - OUTPUT VECTOR OF LENGTH NV CONTAINING THE
C                           EIGENVALUES OF R IN DESCENDING ORDER.
C                Y      - OUTPUT VECTOR OF LENGTH NV CONTAINING THE
C                           COMMUNALITIES OF THE VARIABLES.
C                S      - OUTPUT VECTOR OF LENGTH (NV+1)*NV/2 CONTAINING
C                           THE NV BY NV NORMALIZED RESIDUAL CORRELATION
C                           MATRIX IN SYMMETRIC STORAGE MODE.
C                           (SEE OFRESI FOR NORMALIZATION DESCRIPTION).
C                G      - OUTPUT SCALING VECTOR OF LENGTH NV NEEDED
C                           FOR LATER USE IN THE ROTATION PROGRAMS.
C                IS     - INTEGER WORK VECTOR OF LENGTH NV.
C                WK     - WORK VECTOR OF LENGTH NV.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER = 129 INDICATES THAT AT LEAST ONE OF NF,
C                             NV, CRIT, OR IA WAS SPECIFIED INCORRECTLY.
C                           IER = 130 INDICATES THAT ALL THE EIGENVALUES
C                             OF R COULD NOT BE CALCULATED.
C                           IER = 131 INDICATES THAT WHEN NF IS DETER-
C                             MINED BY CRIT, THE LARGEST EIGENVALUE IS
C                             LESS THAN CRIT.
C                         WARNING ERROR (WITH FIX)
C                           IER = 68 INDICATES THAT A NEGATIVE EIGEN-
C                             VALUE WAS ENCOUNTERED BY OFPRI WHILE
C                             USING THE NF SPECIFIED BY THE USER. NF IS
C                             SET TO THE NUMBER OF POSITIVE EIGENVALUES
C                             OBTAINED AND CALCULATIONS CONTINUE.
C
C   REQD. IMSL ROUTINES - SINGLE/EHOBKS,EHOUSS,EIGRS,EQRT2S,OFRESI,
C                           UERTST,UGETIO,VSRTU
C                       - DOUBLE/EHOBKS,EHOUSS,EIGRS,EQRT2S,OFRESI,
C                           UERTST,UGETIO,VSRTUD
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE OFPRI (R,NV,NF,CRIT,A,IA,E,Y,S,G,IS,WK,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            NV,NF,IA,IS(NV),IER
      REAL               R(1),E(NV),A(IA,1),S(1),G(NV),Y(NV),WK(NV),CRIT
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            NVV,NVP1,I,M,K,J,KK
      REAL               DD,EXCH,ONE
      DATA               ONE /1.0E0/
C                                  ERROR CHECKS
C                                  FIRST EXECUTABLE STATEMENT
      IER = 129
      IF (NF.GT.0) GO TO 5
      IF (CRIT.LT.0.0) GO TO 9000
    5 IF (NF.GT.NV .OR. IA.LT.NV) GO TO 9000
      NVV = ((NV+1)*NV)/2
      NVP1 = NV+1
      IER = 0
      KK = 0
      DO 10 I=1,NVV
         S(I) = R(I)
   10 CONTINUE
C                                  CALCULATE EIGENVALUES OF R
      CALL EIGRS(S,NV,1,E,A,IA,WK,IER)
      IF (IER.EQ.0) GO TO 15
      IER = 130
      GO TO 9000
   15 M = NV/2
C                                  PUT EIGENVALUES IN DESCENDING ORDER
      DO 20 I=1,M
         K = NVP1-I
         EXCH = E(K)
         E(K) = E(I)
         E(I) = EXCH
         IS(I) = K
         IS(K) = I
   20 CONTINUE
      IF (M+M.NE.NV) IS(M+1) = M+1
      CALL VSRTU(A,IA,NV,NV,0,IS,WK)
C                                  CALCULATE THE NUMBER OF FACTORS
      IF (NF.GT.0) GO TO 35
      NF = 0
      DO 25 I=1,NV
         IF (E(I).LE.CRIT) GO TO 30
         NF = NF+1
   25 CONTINUE
   30 IF (NF.GT.0) GO TO 35
      IER = 131
      GO TO 9000
C                                  CALCULATE FACTOR LOADING MATRIX
   35 DO 50 J=1,NF
         IF (E(J).GE.0.0) GO TO 40
         IER = 68
         GO TO 55
   40    DD = SQRT(E(J))
         DO 45 I=1,NV
            A(I,J) = A(I,J)*DD
   45    CONTINUE
         KK = J
   50 CONTINUE
C                                  COMMUNALITIES, RESIDUAL CORRELATIONS
   55 NF = KK
      CALL OFRESI(R,NV,NF,A,IA,Y,S,WK)
C                                  SCALING VECTOR IS UNITY
      DO 60 I=1,NV
         G(I) = ONE
   60 CONTINUE
      IF (IER.EQ.0) GO TO 9005
 9000 CONTINUE
      CALL UERTST(IER,6HOFPRI )
 9005 RETURN
      END
