C   IMSL ROUTINE NAME   - NRWMD
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1980
C
C   PURPOSE             - WILCOXON SIGNED RANK TEST
C
C   USAGE               - CALL NRWMD (IOPT,X,Y,N,EPS,DSEED,IR,STAT,IER)
C
C   ARGUMENTS    IOPT   - INPUT OPTION PARAMETER.
C                         IF IOPT IS EQUAL TO ZERO, THE ROUTINE
C                           COMPUTES THE DIFFERENCES X(I)-Y(I), FOR
C                           I=1,...,N.
C                         IF IOPT IS NOT EQUAL TO ZERO, THE ROUTINE
C                           ASSUMES THAT THE DIFFERENCES X(I)-Y(I),
C                           FOR I=1,...,N, HAVE ALREADY BEEN COMPUTED
C                           AND ARE STORED IN X.
C                X      - INPUT/OUTPUT VECTOR OF LENGTH N.
C                         IF IOPT IS EQUAL TO ZERO, X CONTAINS SAMPLES
C                           FROM THE X POPULATION. ON OUTPUT
C                           X IS REPLACED BY THE ORDERED ABSOLUTE VALUES
C                           OF THE DIFFERENCES X(I)-Y(I), FOR I=1,...,N.
C                         IF IOPT IS NOT EQUAL TO ZERO, X MUST CONTAIN
C                           THE DIFFERENCES AND X IS REPLACED BY THE
C                           ORDERED ABSOLUTE DIFFERENCES.
C                Y      - INPUT/OUTPUT VECTOR OF LENGTH N.
C                         IF IOPT IS EQUAL TO ZERO, Y CONTAINS SAMPLES
C                           FROM THE Y POPULATION.
C                         IF IOPT IS NOT EQUAL TO ZERO, Y IS NOT
C                           REQUIRED.
C                         ON OUTPUT, Y IS REPLACED BY THE SIGNED RANKS
C                           OF THE ABSOLUTE DIFFERENCES.
C                N      - INPUT NUMBER OF SAMPLES FOR EACH OF THE
C                           POPULATIONS. N MUST BE GREATER THAN OR
C                           EQUAL TO 3.
C                EPS    - INPUT POSITIVE VALUE TO BE USED TO
C                           DETECT A ZERO DIFFERENCE OF X(I)-Y(I),
C                           FOR I=1,...,N, AND TIES BETWEEN THE
C                           DIFFERENCES. IF ANY DIFFERENCE IS LESS THAN
C                           EPS, IT IS COUNTED AS A ZERO DIFFERENCE.
C                           FURTHERMORE, IF ANY DIFFERENCE IS WITHIN
C                           EPS TIMES ITSELF OF ANOTHER DIFFERENCE,
C                           A TIE IS COUNTED.  IF EPS IS LESS THAN
C                           10**-5, A VALUE OF 10**-5 IS USED IN PLACE
C                           OF EPS.
C                DSEED  - INPUT/OUTPUT DOUBLE PRECISION VARIABLE
C                           ASSIGNED AN INTEGER VALUE IN THE
C                           EXCLUSIVE RANGE (1.D0, 2147483647.D0).
C                           DSEED IS REPLACED BY A NEW VALUE TO BE
C                           USED IN A SUBSEQUENT CALL.
C                IR     - WORK VECTOR OF LENGTH N. ON OUTPUT IR CONTAINS
C                           THE PERMUTATIONS REQUIRED IN SORTING THE
C                           ABSOLUTE DIFFERENCES (ATTACHED TO THE SIGNS
C                           OF THE DIFFERENCES).  SEE REMARKS.
C                STAT   - OUTPUT VECTOR OF LENGTH 5.
C                           STAT(1) CONTAINS W+, THE POSITIVE RANK SUM.
C                           STAT(2) CONTAINS ABS(W-), THE ABSOLUTE VALUE
C                             OF THE NEGATIVE RANK SUM.
C                           STAT(3) CONTAINS THE STANDARDIZED MINIMUM OF
C                             (W+,W-)
C                           STAT(4) CONTAINS THE PROBABILITY OF NOT
C                             EXCEEDING STAT(3) IF THE HYPOTHESIS OF
C                             IDENTICAL POPULATIONS IS TRUE.
C                           STAT(5) CONTAINS THE NUMBER OF ZERO
C                             DIFFERENCES.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER=129 MEANS THAT ALL DIFFERENCES WERE ZERO
C                           IER=130 MEANS THAT N IS LESS THAN 3.
C                         WARNING ERROR
C                           IER=35 MEANS THAT N IS LESS THAN 50, AND
C                             EXACT TABLES SHOULD BE REFERENCED FOR
C                             PROBABILITIES.
C                           IER=36 MEANS THAT TIES EXIST.
C
C   REQD. IMSL ROUTINES - GGUBS,MERRC=ERFC,UERTST,UGETIO,VSRTP
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS      THE VECTOR IR CONTAINS A RECORD OF THE PERMUTATIONS
C                MADE IN SORTING THE X VECTOR OF DIFFERENCES. IR(I)=J
C                IMPLIES THAT ELEMENT I OF THE OUTPUT X VECTOR OF
C                DIFFERENCES WAS IN POSITION J ON INPUT FOR I AND
C                J=1,...,N. THE SIGN OF IR(I) INDICATES THE SIGN
C                OF THE DIFFERENCE ON INPUT.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE NRWMD  (IOPT,X,Y,N,EPS,DSEED,IR,STAT,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            IOPT,N,IR(N),IER
      REAL               X(N),Y(N),EPS,STAT(5)
      DOUBLE PRECISION   DSEED
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,ICT,IER1,I1,J,K,L
      REAL               RNARAY(1),EP1,RN,R1DS24
C                                  FIRST EXECUTABLE STATEMENT
      R1DS24 = 1.0/SQRT(24.0)
      IF (IOPT .NE. 0) GO TO 10
C                                  SET X = X-Y
      DO 5 I=1,N
    5 X(I) = X(I)-Y(I)
C                                  TEST FOR VALID N VALUE
   10 IF(N.LT.3) GO TO 80
      ICT = 0
      IER = 0
      IER1 = 0
      IF(N.LT.50) IER = 35
      EP1 = EPS
      IF(EPS .LT. 1.0E-5) EP1 = 1.0E-5
C                                  INITIALIZE IR FOR SORT ROUTINE
      DO 15 I=1,N
   15 IR(I) = I
      DO 20 I=1,N
         IF(X(I) .LT. 0.0) IR(I) = -IR(I)
         IF(ABS(X(I)) .GE. EP1) GO TO 20
         X(I) = 0.0
         ICT = ICT+1
   20 CONTINUE
      CALL UGETIO (1,NIN,NOUT)
      IF(ICT .EQ. N) GO TO 75
C                                  GO SORT THE X VECTOR
      CALL VSRTP (X,N,IR)
      ISTAR = 1
C                                  ASSIGN PLUS OR MINUS TO RANKS OF
C                                  ZEROS IN X VECTOR
      DO 25 I=1,N
         IF(X(I) .NE. 0.0) GO TO 30
         ISTAR =ISTAR+1
         CALL GGUBS (DSEED,1,RNARAY)
         RN = RNARAY(1)
         IF(RN .LE. 0.5) IR(I) = -IR(I)
   25 CONTINUE
   30 IFQ = 1
      K = N-1
      DO 50 I=ISTAR,K
C                                   DETERMINE TIE SUBSETS
         IF(X(I+1) .GT. X(I)*(1.0+EP1)) GO TO 35
         IFQ = IFQ+1
         GO TO 50
   35    IF (IFQ .EQ. 1) GO TO 50
         I1 = I-IFQ+1
         M = IFQ
C                                   RAMDOMLY PERMUTE ARRAY IR
C                                   FOR THE TIE SUBSETS
   40 CALL GGUBS (DSEED,1,RNARAY)
        J = INT(FLOAT(M)*RNARAY(1))+I1
        TEMP = IR(I1+M-1)
        IR(I1+M-1) = IR(J)
        IR(J) = TEMP
        M = M-1
        IF (M .LE. 1) GO TO 45
      GO TO 40
   45 IFQ = 1
   50 CONTINUE
C                                   NOW GET THE SIGNED RANKS
C                                   FOR X
      DO 60 J=1,N
         K = J
         L = IR(J)
         IF (L .GT. 0) GO TO 55
         K = -K
         L = -L
   55    Y(L) = K
   60 CONTINUE
C
   65 STAT(1) = 0.0
      DO 70 I=1,N
         IF(Y(I) .GT. 0.0) STAT(1) = STAT(1)+Y(I)
   70 CONTINUE
C                                  CALCULATE ABS VALUE OF NEG RANK SUM
      TEMP1 = N
      TEMP2 = TEMP1*(TEMP1+1.0)
      STAT(2) = ABS(TEMP2*.5-STAT(1))
      STAT(3) = AMIN1(STAT(1),STAT(2))
      TEMP3 = TEMP2*(TEMP1+TEMP1+1.0)
      STAT(3) = (STAT(3)-.25*TEMP2)/(SQRT(TEMP3)*R1DS24)
      STAT(4) = 0.5*ERFC(-.7071068*STAT(3))
      STAT(5) = ICT
      GO TO 9000
   75 IER = 129
      STAT(5) = ICT
      GO TO 9000
   80 IER = 130
      IER1 = 0
 9000 CONTINUE
      IF(IER.GT.0) CALL UERTST(IER,6HNRWMD )
      IF(IER1.GT.0) CALL UERTST(IER1,6HNRWMD )
      IER = MAX0(IER,IER1)
 9005 RETURN
      END
