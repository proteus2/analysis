C   IMSL ROUTINE NAME   - ALSQAN
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - ANALYSIS OF LATIN SQUARE DESIGN DATA
C
C   USAGE               - CALL ALSQAN (Y,NR,NT,IND,EM,GM,S,NDF,IER)
C
C   ARGUMENTS    Y      - INPUT VECTOR OF LENGTH NR*NT*NT CONTAINING
C                           THE RESPONSES FOR THE TREATMENT IN EACH
C                           ROW-COLUMN POSITION OF THE LATIN SQUARE.
C                           THE RESPONSES ARE ARRANGED IN A SPECIFIC
C                           MANNER. THE FIRST NR COMPONENTS OF Y CONTAIN
C                           THE RESPONSES FOR ROW ONE, COLUMN ONE, THE
C                           SECOND NR COMPONENTS CONTAIN THE RESPONSES
C                           FOR ROW ONE, COLUMN TWO, AND SO FORTH.
C                NR     - INPUT NUMBER OF REPLICATE RESPONSES IN EACH
C                           ROW-COLUMN POSITION. NR MUST BE GREATER
C                           THAN OR EQUAL TO ONE.
C                NT     - INPUT NUMBER OF ROWS, COLUMNS, AND TREATMENTS.
C                           NT MUST BE GREATER THAN OR EQUAL TO THREE.
C                IND    - INPUT VECTOR OF LENGTH NT*NT CONTAINING THE
C                           TREATMENT NUMBERS FOR THE NT*NT ROW-COLUMN
C                           POSITIONS. THE TREATMENT NUMBERS ARE ASSUMED
C                           TO BE FROM THE SET 1,2,...,NT. THE ORDERING
C                           OF THE TREATMENT NUMBERS IN IND CORRESPONDS
C                           TO THE ORDERING IN VECTOR Y.
C                EM     - OUTPUT VECTOR OF LENGTH NT*(NT+3) CONTAINING
C                           THE ROW, COLUMN, TREATMENT, AND CELL MEANS.
C                           THE MEANS ARE ARRANGED IN A SPECIFIC
C                           MANNER. THE FIRST NT COMPONENTS OF EM
C                           CONTAIN THE ROW MEANS. THE NEXT NT CONTAIN
C                           THE COLUMN MEANS. THE NEXT NT CONTAIN
C                           THE TREATMENT MEANS. IF NR IS GREATER THAN
C                           ONE, THE LAST NT*NT COMPONENTS CONTAIN
C                           THE CELL MEANS (SEE REMARKS). THE ORDERING
C                           OF THE MEANS WITHIN THE FOUR GROUPS OF
C                           MEANS CORRESPONDS TO THE ORDERING IN
C                           VECTOR Y.
C                GM     - OUTPUT GRAND MEAN OF THE RESPONSES.
C                S      - OUTPUT VECTOR OF LENGTH 6 CONTAINING THE SUMS
C                           OF SQUARES.
C                           S(1) CONTAINS THE ROW SUM OF SQUARES
C                           S(2) CONTAINS THE COLUMN SUM OF SQUARES
C                           S(3) CONTAINS THE TREATMENT SUM OF SQUARES
C                           S(4) CONTAINS THE EXPERIMENTAL ERROR SUM OF
C                             SQUARES
C                           S(5) CONTAINS THE SAMPLING ERROR (WITHIN
C                             CELL) SUM OF SQUARES
C                           S(6) CONTAINS THE CORRECTED TOTAL SUM OF
C                             SQUARES
C                NDF    - OUTPUT VECTOR OF LENGTH 4 CONTAINING THE
C                           DEGREES OF FREEDOM CORRESPONDING TO
C                           COMPONENTS OF THE S VECTOR.
C                           NDF(1) CONTAINS THE ROW, COLUMN, AND
C                             TREATMENT DEGREES OF FREEDOM
C                           NDF(2) CONTAINS THE EXPERIMENTAL ERROR
C                             DEGREES OF FREEDOM
C                           NDF(3) CONTAINS THE SAMPLING ERROR DEGREES
C                             OF FREEDOM
C                           NDF(4) CONTAINS THE CORRECTED TOTAL DEGREES
C                             OF FREEDOM
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER=129 INDICATES THAT NT WAS SPECIFIED
C                             LESS THAN 3 OR THAT NR WAS SPECIFIED LESS
C                             THAN 1.
C                           IER=130 INDICATES THAT THE INPUT
C                             EXPERIMENTAL DESIGN WAS NOT A LATIN
C                             SQUARE.
C
C   REQD. IMSL ROUTINES - UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS      IF NR IS EQUAL TO ONE, THEN S(5) AND NDF(3) ARE
C                EQUAL TO ZERO AND THE CELL MEANS ARE NOT DEFINED. IN
C                THIS CASE THE EM VECTOR WILL CONTAIN ONLY ROW, COLUMN,
C                AND TREATMENT MEANS AND THUS MAY BE OF LENGTH 3*NT.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE ALSQAN (Y,NR,NT,IND,EM,GM,S,NDF,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            NR,NT,IND(1),NDF(4),IER
      REAL               Y(1),EM(1),GM,S(6)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,J,K,KK,KKK,L,N,NS,NSUM,NTR,NT1,NT2,NT3
      REAL               TEMP
      DOUBLE PRECISION   Z,R,C,T,TOT
C                                  FIRST EXECUTABLE STATEMENT
      IER=0
      IF (NT .GE. 3 .AND. NR .GE. 1) GO TO 5
C                                  TERMINAL ERROR - NT WAS LESS THAN 3
C                                  OR NR WAS LESS THAN 1
      IER = 129
      GO TO 9000
C                                  INITIALIZE VECTOR OF MEANS
    5 N = NT*3
      IF (NR .GT. 1) N = N+NT*NT
      DO 10 I = 1,N
         EM(I) = 0.0
   10 CONTINUE
      NT1 = NT+1
      NSUM = (NT*NT1)/2
      DO 25 I = 1,NT
         NS = 0
         K = I
         DO 15 J = 1,NT
            NS = NS+IND(K)
            K = K+NT
   15    CONTINUE
         IF (NS .NE. NSUM) GO TO 80
C                                  TERMINAL ERROR - INPUT DESIGN IS
C                                  NOT A LATIN SQUARE
   25 CONTINUE
      NT2 = NT+NT
      NTR = NT*NR
      NT3 = NT2+NT
      GM = 0.0
      K = 1
      KKK = 1
      KK = NT2+NT1
      TEMP = 1.0/NR
      DO 40 I = 1,NT
         NS = 0
         DO 35 J = 1,NT
            NS = NS+IND(KKK)
            Z = 0.0D0
            DO 30 L = 1,NR
               Z = Z+Y(K)
               K = K+1
   30       CONTINUE
C                                  CALCULATE CELL AVERAGES
            IF (NR .GT. 1) EM(KK) = Z*TEMP
            N = NT+J
C                                  CALCULATE COLUMN TOTALS
            EM(N) = EM(N)+Z
            N = IND(KKK)+NT2
C                                  CALCULATE TREATMENT TOTALS
            EM(N) = EM(N)+Z
C                                  CALCULATE ROW TOTALS
            EM(I) = EM(I)+Z
            KKK = KKK+1
            KK = KK+1
   35    CONTINUE
         IF (NS .NE. NSUM) GO TO 80
         GM = GM+EM(I)
   40 CONTINUE
C                                  CALCULATE AVERAGES
      TEMP = 1.0/NTR
      DO 45 I = 1,NT3
         EM(I) = EM(I)*TEMP
   45 CONTINUE
      GM = GM/(NTR*NT)
      L = 0
      R = 0.0D0
      C = 0.0D0
      T = 0.0D0
      TOT = 0.0D0
C                                  CALCULATE SUMS OF SQUARES
      DO 55 I = 1,NT
C                                  FOR ROWS
         Z = EM(I)-GM
         R = R+Z*Z
C                                  FOR COLUMNS
         Z = EM(NT+I)-GM
         C = C+Z*Z
C                                  FOR TREATMENTS
         Z = EM(NT2+I)-GM
         T = T+Z*Z
         DO 50 J = 1,NTR
            Z = Y(J+L)-GM
            TOT = TOT+Z*Z
   50    CONTINUE
         L = L+NTR
   55 CONTINUE
      S(1) = NTR*R
      S(2) = NTR*C
      S(3) = NTR*T
      NT2 = NT*NT
      IF (NR .EQ. 1) GO TO 70
      K = 1
      R = 0.0D0
      J = 1
      DO 65 I = 1,NT2
         DO 60 L = 1,NR
            Z = Y(K)-EM(NT3+J)
            R = R+Z*Z
            K = K+1
   60    CONTINUE
         J = J+1
   65 CONTINUE
      S(5) = R
      NDF(3) = NT2*(NR-1)
      GO TO 75
   70 S(5) = 0.0
      NDF(3) = 0
   75 S(4) = TOT-S(1)-S(2)-S(3)-S(5)
      S(6) = TOT
C                                  CALCULATE DEGREES OF FREEDOM
      NDF(1) = NT-1
      NDF(2) = NDF(1)*(NT-2)
      NDF(4) = NT*NT*NR-1
      GO TO 9005
   80 IER = 130
 9000 CONTINUE
      CALL UERTST (IER,'ALSQAN')
 9005 RETURN
      END
