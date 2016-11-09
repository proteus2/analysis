C   IMSL ROUTINE NAME   - CTLLF
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - NOVEMBER 1, 1984
C
C   PURPOSE             - LOG-LINEAR FIT OF CONTINGENCY TABLE
C
C   USAGE               - CALL CTLLF (A,NVAR,NVAL,NMARG,MARG,IMARG,
C                           EPS,MAXIT,AFIT,AMAR,INDEX,DEV,WK,IWK,IER)
C
C   ARGUMENTS    A      - INPUT VECTOR OF LENGTH
C                           NVAL(1)*NVAL(2)*...*NVAL(NVAR) CONTAINING
C                           THE TABLE TO BE FIT.
C                NVAR   - INPUT NUMBER OF CLASSIFICATION VARIABLES.
C                NVAL   - INPUT VECTOR OF LENGTH NVAR CONTAINING
C                           THE NUMBER OF CATEGORIES OF EACH
C                           CLASSIFICATION VARIABLE.
C                NMARG  - INPUT NUMBER OF MARGINAL TABLES TO BE FIT.
C                MARG   - INPUT NVAR BY NMARG MATRIX CONTAINING
C                           IN ITS COLUMNS THE SUFFICIENT
C                           CONFIGURATIONS SPECIFYING THE MODEL.
C                           MARG INDICATES THE MARGINAL TABLES
C                           TO BE FIT.
C                IMARG  - INPUT ROW DIMENSION OF MARG EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT
C                           IN THE CALLING PROGRAM.
C                EPS    - INPUT MAXIMUM PERMISSIBLE DEVIATION
C                           BETWEEN AN OBSERVED AND A FITTED
C                           MARGINAL TOTAL.
C                MAXIT  - INPUT/OUTPUT VARIABLE WHICH ON INPUT IS
C                           THE MAXIMUM PERMISSIBLE NUMBER OF
C                           ITERATIONS AND WHICH ON OUTPUT IS THE
C                           NUMBER OF ITERATIONS ACTUALLY PERFORMED.
C                AFIT   - INPUT/OUTPUT VECTOR OF LENGTH
C                           NVAL(1)*NVAL(2)...*NVAL(NVAR).  ON INPUT
C                           AFIT IS USED TO INDICATE STRUCTURAL
C                           ZEROS AND TO PROVIDE INITIAL ESTIMATES.
C                           SEE REMARKS.  ON OUTPUT AFIT IS THE
C                           FITTED TABLE.
C                AMAR   - OUTPUT VECTOR CONTAINING THE MARGINAL
C                           TABLES TO BE FIT.  THE DIMENSION OF
C                           AMAR IS THE SUM OVER J OF THE PRODUCTS
C                           OVER I OF NVAL(MARG(I,J)).  THIS
C                           DIMENSION WILL NEVER BE GREATER
C                           THAN NMARG TIMES THE DIMENSION OF A.
C                INDEX  - OUTPUT VECTOR OF LENGTH NMARG
C                           CONTAINING THE BEGINNING INDICES OF
C                           AMAR FOR THE MARGINAL TABLES.
C                DEV    - OUTPUT VECTOR OF LENGTH MAXIT
C                           CONTAINING THE MAXIMUM DIFFERENCE
C                           ON EACH ITERATION BETWEEN THE
C                           OBSERVED AND FITTED MARGINAL TOTALS.
C                WK     - WORK VECTOR OF LENGTH EQUAL TO THE
C                           LARGEST MARGINAL TABLE SPECIFIED,
C                           THAT IS, THE MAXIMUM OVER J OF THE PRODUCT
C                           OVER I OF NVAL(MARG(I,J)), FOR MARG(I,J)
C                           GREATER THAN ZERO.
C                IWK    - WORK VECTOR OF LENGTH 4*NVAR+1.
C                IER    - ERROR PARAMETER.  (OUTPUT)
C                         WARNING ERROR
C                           IER = 33, INDICATES THE ALGORITHM DID
C                             NOT CONVERGE TO THE DESIRED ACCURACY
C                             WITHIN MAXIT ITERATIONS.
C                         TERMINAL ERROR
C                           IER = 129, INDICATES THAT NVAR WAS
C                             SPECIFIED LESS THAN OR EQUAL TO ONE.
C                           IER = 130, INDICATES THAT NVAL WAS
C                             SPECIFIED LESS THAN OR EQUAL TO ZERO.
C                           IER = 131, INDICATES THAT A WAS SPECIFIED
C                             LESS THAN ZERO.
C                           IER = 132, INDICATES THAT AFIT WAS
C                             SPECIFIED LESS THAN ZERO.
C                           IER = 133, INDICATES THAT EVERY ELEMENT
C                             IN AFIT WAS SPECIFIED EQUAL TO ZERO.
C                           IER = 134, INDICATES THAT MARG CONTAINS
C                             AN ERROR.
C
C   REQD. IMSL ROUTINES - CTLL1,CTLL2,UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS  1.  MARG IS USED TO DESCRIBE THE MARGINAL
C                TABLES TO BE FIT.  FOR EXAMPLE WITH NVAR=3
C                IF THE MARGINAL TABLES FOR VARIABLES 1 AND 2
C                AND VARIABLES 1 AND 3 ARE TO BE FIT, NMARG=2 AND
C                    MARG(1,1) = 1             MARG(1,2) = 1
C                    MARG(2,1) = 2             MARG(2,2) = 3
C                    MARG(3,1) = 0             MARG(3,2) = 0.
C            2.  ON INPUT, IF THE USER DOES NOT WISH TO SPECIFY
C                INITIAL ESTIMATES, AFIT(I) SHOULD BE SET TO 1.0,
C                FOR ANY CELL WHICH MAY HAVE A POSITIVE FITTED
C                VALUE.  A STRUCTURAL ZERO (INCOMPLETE TABLE) IS
C                SPECIFIED BY THE CORRESPONDING VALUE OF AFIT
C                SET TO 0.0 ON INPUT.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE CTLLF  (A,NVAR,NVAL,NMARG,MARG,IMARG,EPS,MAXIT,
     *                   AFIT,AMAR,INDEX,DEV,WK,IWK,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            NVAR,NVAL(NVAR),NMARG,IMARG,MARG(IMARG,1),
     *                   MAXIT,INDEX(NMARG),IWK(1),IER
      REAL               A(1),EPS,AFIT(1),AMAR(1),DEV(1),WK(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            II,IPOINT,ISIZE,I,J,KK,K,NV2P1,N
      REAL               XMAX,X,Y
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
      IF (NVAR.LE.1) GO TO 85
      IF (MAXIT.GT.0) GO TO 5
      MAXIT = 1
    5 NV2P1 = 2*NVAR+1
C                                  LOOK AT TABLE AND FIT CONSTANTS
      ISIZE = 1
      DO 10 J=1,NVAR
         IF (NVAL(J).LE.0) GO TO 90
         ISIZE = ISIZE*NVAL(J)
   10 CONTINUE
      X = 0.0
      Y = 0.0
      DO 15 I=1,ISIZE
         IF (A(I).LT.0.0) GO TO 95
         IF (AFIT(I).LT.0.0) GO TO 100
         X = X+A(I)
         Y = Y+AFIT(I)
   15 CONTINUE
C                                  MAKE A PRELIMINARY ADJUSTMENT TO
C                                    OBTAIN THE FIT TO AN EMPTY
C                                    CONFIGURATION LIST
      IF (Y.EQ.0.0) GO TO 105
      X = X/Y
      DO 20 I=1,ISIZE
         AFIT(I) = X*AFIT(I)
   20 CONTINUE
C                                  ALLOCATE MARGINAL TABLES
      IPOINT = 1
      IF (NMARG.EQ.0) GO TO 9005
      DO 45 I=1,NMARG
         II = I
         IF (MARG(1,I).EQ.0) GO TO 50
C                                  GET MARGINAL TABLE SIZE. AND SEE
C                                    IF THE CONFIGURATION LIST
C                                    CONTAINS DUPLICATIONS OR
C                                    ELEMENTS OUT OF RANGE.
         ISIZE = 1
         DO 25 J=1,NVAR
            IWK(J+NVAR) = 0
   25    CONTINUE
         DO 35 J=1,NVAR
            K = MARG(J,I)
C                                  A ZERO INDICATES THE END OF THE
C                                    STRING
            IF (K.EQ.0) GO TO 40
C                                  SEE IF ELEMENT VALID
            IF (K.GT.0.AND.K.LE.NVAR) GO TO 30
            GO TO 110
C                                  CHECK FOR DUPLICATION
   30       IF (IWK(K+NVAR).EQ.1) GO TO 110
            IWK(K+NVAR) = 1
            ISIZE = ISIZE*NVAL(K)
   35    CONTINUE
C                                  INDEX POINTS TO MARGINAL TABLES TO
C                                    BE PLACED IN MARG
   40    INDEX(I) = IPOINT
         IPOINT = IPOINT+ISIZE
   45 CONTINUE
C                                  GET N, NUMBER OF VALID
C                                    CONFIGURATIONS
      II = NMARG+1
   50 N = II-1
      IF (N.EQ.0) GO TO 110
C                                  OBTAIN MARGINAL TABLES
      DO 60 I=1,N
         DO 55 J=1,NVAR
            IWK(J) = MARG(J,I)
   55    CONTINUE
         CALL CTLL1 (NVAR,A,AMAR,1,INDEX(I),NVAL,IWK(1),IWK(NV2P1))
   60 CONTINUE
C                                  PERFORM ITERATIONS
      DO 75 K=1,MAXIT
C                                  XMAX IS MAXIMUM DEVIATION OBSERVED
C                                    BETWEEN FITTED AND TRUE MARGINAL
C                                    DURING A CYCLE
         KK = K
         XMAX = 0.0
         DO 70 I=1,N
            DO 65 J=1,NVAR
               IWK(J) = MARG(J,I)
   65       CONTINUE
            CALL CTLL1 (NVAR,AFIT,WK,1,1,NVAL,IWK(1),IWK(NV2P1))
            CALL CTLL2 (NVAR,AFIT,WK,AMAR,1,1,INDEX(I),NVAL,IWK(1),XMAX,
     *        IWK(NV2P1))
   70    CONTINUE
C                                  TEST CONVERGENCE
         DEV(K) = XMAX
         IF (XMAX.LT.EPS) GO TO 80
   75 CONTINUE
C                                  NO CONVERGENCE
      IER = 33
      CALL UERTST (IER,6HCTLLF )
   80 MAXIT = KK
      GO TO 9005
   85 IER = 129
      GO TO 9000
   90 IER = 130
      GO TO 9000
   95 IER = 131
      GO TO 9000
  100 IER = 132
      GO TO 9000
  105 IER = 133
      GO TO 9000
  110 IER = 134
 9000 CONTINUE
      CALL UERTST (IER,'CTLLF ')
 9005 RETURN
      END
