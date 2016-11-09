C   IMSL ROUTINE NAME   - ARCBAN
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - ANALYSIS OF TWO-WAY CLASSIFICATION DESIGN
C                           DATA
C
C   USAGE               - CALL ARCBAN (Y,NR,NB,NT,EM,GM,S,NDF,IER)
C
C   ARGUMENTS    Y      - INPUT VECTOR OF LENGTH NR*NB*NT CONTAINING
C                           THE RESPONSES FOR EACH TREATMENT-BLOCK
C                           COMBINATION. THE RESPONSES ARE ARRANGED
C                           IN A SPECIFIC MANNER. THE FIRST NR
C                           COMPONENTS OF Y CONTAIN THE RESPONSES
C                           FOR TREATMENT ONE, BLOCK ONE. THE SECOND
C                           NR COMPONENTS CONTAIN THE RESPONSES FOR
C                           TREATMENT TWO, BLOCK ONE, AND SO FORTH.
C                NR     - INPUT NUMBER OF REPLICATE RESPONSES FOR EACH
C                           TREATMENT-BLOCK COMBINATION. NR MUST BE
C                           GREATER THAN OR EQUAL TO ONE.
C                NB     - INPUT NUMBER OF BLOCKS. NB MUST BE GREATER
C                           THAN OR EQUAL TO TWO.
C                NT     - INPUT NUMBER OF TREATMENTS. NT MUST BE GREATER
C                           THAN OR EQUAL TO TWO.
C                EM     - OUTPUT VECTOR OF LENGTH NB+NT+NB*NT CONTAINING
C                           THE BLOCK, TREATMENT, AND CELL MEANS,
C                           RESPECTIVELY. THE MEANS ARE ARRANGED IN A
C                           SPECIFIC MANNER. THE FIRST NB COMPONENTS
C                           OF EM CONTAIN THE BLOCK MEANS. THE NEXT NT
C                           COMPONENTS CONTAIN THE TREATMENT MEANS.
C                           IF NR IS GREATER THAN ONE, THE LAST NB*NT
C                           COMPONENTS CONTAIN THE CELL MEANS (SEE
C                           REMARKS). THE ORDERING OF THE MEANS
C                           WITHIN THE THREE GROUPS OF MEANS
C                           CORRESPONDS TO THE ORDERING IN VECTOR Y.
C                GM     - OUTPUT GRAND MEAN OF THE RESPONSES.
C                S      - OUTPUT VECTOR OF LENGTH 5 CONTAINING SUMS OF
C                           SQUARES.
C                           S(1) CONTAINS THE BLOCK SUM OF SQUARES
C                           S(2) CONTAINS THE TREATMENT SUMS OF SQUARES
C                           S(3) CONTAINS THE EXPERIMENTAL ERROR
C                             (TREATMENT BY BLOCK INTERACTION) SUM OF
C                             SQUARES
C                           S(4) CONTAINS THE SAMPLING ERROR (WITHIN
C                             CELL) SUM OF SQUARES
C                           S(5) CONTAINS THE CORRECTED TOTAL SUM OF
C                             SQUARES
C                NDF    - OUTPUT VECTOR OF LENGTH 5 CONTAINING DEGREES
C                           OF FREEDOM CORRESPONDING TO COMPONENTS OF
C                           THE S VECTOR.
C                           NDF(1) CONTAINS THE BLOCK DEGREES OF FREEDOM
C                           NDF(2) CONTAINS THE TREATMENT DEGREES OF
C                             FREEDOM
C                           NDF(3) CONTAINS THE EXPERIMENTAL ERROR
C                             DEGREES OF FREEDOM
C                           NDF(4) CONTAINS THE SAMPLING ERROR DEGREES
C                             OF FREEDOM
C                           NDF(5) CONTAINS THE CORRECTED TOTAL DEGREES
C                             OF FREEDOM
C                IER    - ERROR PARAMETER. (OUTPUT)
C                           TERMINAL ERROR
C                             IER=129 INDICATES THAT NB WAS SPECIFIED
C                               LESS THAN TWO.
C                             IER=130 INDICATES THAT NT WAS SPECIFIED
C                               LESS THAN TWO.
C                             IER=131 INDICATES THAT NR WAS SPECIFIED
C                               LESS THAN ONE.
C
C   REQD. IMSL ROUTINES - UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS      IF NR IS EQUAL TO ONE, THEN S(4) AND NDF(4) ARE
C                EQUAL TO ZERO AND THE CELL MEANS ARE NOT DEFINED. IN
C                THIS CASE THE EM VECTOR WILL CONTAIN ONLY BLOCK AND
C                TREATMENT MEANS AND THUS MAY BE OF LENGTH NB+NT.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE ARCBAN (Y,NR,NB,NT,EM,GM,S,NDF,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            NR,NB,NT,NDF(5),IER
      REAL               Y(1),EM(1),GM,S(5)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,J,K,KK,L,N,NTR,NTRB,NTRB1,NT2,NT3
      DOUBLE PRECISION   Z,SUM,XNR
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
      IF (NB .GE. 2) GO TO 5
C                                  TERMINAL ERROR - NB SPECIFIED
C                                  INCORRECTLY
      IER = 129
      GO TO 9000
    5 IF (NT .GE. 2) GO TO 10
C                                  TERMINAL ERROR - NT SPECIFIED
C                                  INCORRECTLY
      IER = 130
      GO TO 9000
   10 IF (NR .GE. 1) GO TO 15
C                                  TERMINAL ERROR - NR SPECIFIED
C                                  INCORRECTLY
      IER = 131
      GO TO 9000
C                                  INITIALIZE VECTOR OF MEANS
   15 NTRB = NB+NT
      IF (NR .GT. 1) NTRB = NTRB+NB*NT
      DO 20 I = 1,NTRB
   20 EM(I) = 0.
      NT2 = NB+1
      NT3 = NT+NT2
      NTR = NT3-1
      GM = 0.
      K = 1
      KK = NT3
      XNR = 1.D0/NR
      DO 30 I = 1,NB
         DO 30 J = 1,NT
C                                  CALCULATE CELL AVERAGES
            Z = 0.D0
            DO 25 L = 1,NR
               Z = Z+Y(K)
   25       K = K+1
            IF (NR .GT. 1) EM(KK) = Z*XNR
C                                  CALCULATE TREATMENT TOTAL
            N = NB+J
            EM(N) = Z+EM(N)
C                                  CALCULATE BLOCK TOTAL
            EM(I) = EM(I)+Z
C                                  CALCULATE GRAND TOTAL
           GM = GM+Z
   30    KK = KK+1
      NTRB1 = NT*NB*NR
      GM = GM/NTRB1
C                                  CALCULATE SUMS OF SQUARES
C                                  AND DEGREES OF FREEDOM
      SUM = 0.D0
      NTRB = NT*NR
      XNR = 1.D0/NTRB
      DO 35 I = 1,NB
         EM(I) = EM(I)*XNR
         Z = EM(I)-GM
   35 SUM = SUM+Z*Z
      S(1) = SUM*NTRB
      NDF(1) = NB-1
      NTRB = NR*NB
      XNR = 1.D0/NTRB
      SUM = 0.D0
      DO 40 I = NT2,NTR
         EM(I) = EM(I)*XNR
         Z = EM(I)-GM
   40 SUM = SUM+Z*Z
      S(2) = SUM*NTRB
      NDF(2) = NT-1
      SUM = 0.D0
      DO 45 I = 1,NTRB1
         Z = Y(I)-GM
   45 SUM = SUM+Z*Z
      S(5) = SUM
      NDF(5) = NTRB1-1
      IF (NR .EQ. 1) GO TO 60
      K = 1
      SUM = 0.D0
      J = 0
      NTR = NT*NB
      DO 55 I = 1,NTR
         DO 50 L = 1,NR
            Z = Y(K)-EM(NT3+J)
            SUM = SUM+Z*Z
   50    K = K+1
   55 J = J+1
      S(4) = SUM
      NDF(4) = (NR-1)*NB*NT
      GO TO 65
   60 S(4) = 0.
      NDF(4) = 0
   65 S(3) = S(5)-S(1)-S(2)-S(4)
      NDF(3) = NDF(1)*NDF(2)
      GO TO 9005
 9000 CONTINUE
      CALL UERTST (IER,'ARCBAN')
 9005 RETURN
      END
