C   IMSL ROUTINE NAME   - GGTAB
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1982
C
C   PURPOSE             - GENERATE A RANDOM CONTINGENCY TABLE
C                          WITH GIVEN ROW AND COLUMN TOTALS
C
C   USAGE               - CALL GGTAB(DSEED,NROW,NCOL,NRTOT,NCTOT,IND,
C                                    IIT,ITAB,IWK,IER)
C
C   ARGUMENTS    DSEED  - INPUT/OUTPUT DOUBLE PRECISION VARIABLE
C                           ASSIGNED AN INTEGER VALUE IN THE
C                           EXCLUSIVE RANGE (1.D0,2147483647.D0).
C                           DSEED IS REPLACED BY A NEW VALUE TO BE
C                           USED IN A SUBSEQUENT CALL.
C                NROW   - INPUT NUMBER OF ROWS IN TABLE.
C                NCOL   - INPUT NUMBER OF COLUMNS IN TABLE.
C                NRTOT  - INPUT VECTOR OF LENGTH NROW CONTAINING
C                          THE ROW TOTALS.
C                NCTOT  - INPUT VECTOR OF LENGTH NCOL CONTAINING
C                          THE COLUMN TOTALS.
C                IND    - INPUT INDICATOR.  IND = 0 INDICATES THE
C                          FIRST CALL TO GGTAB FOR A GIVEN PROBLEM.
C                          FOR SUBSEQUENT CALLS IND MAY BE SET TO A
C                          NONZERO VALUE TO AVOID CERTAIN INITIALI-
C                          ZATIONS AND VALIDITY CHECKS.
C                IIT    - INPUT ROW DIMENSION OF ITAB EXACTLY AS
C                          SPECIFIED IN THE DIMENSION STATEMENT
C                          IN THE CALLING PROGRAM.
C                ITAB   - OUTPUT NROW BY NCOL RANDOM MATRIX SATIS-
C                          FYING THE GIVEN ROW AND COLUMN TOTALS.
C                IWK    - WORK VECTOR OF LENGTH EQUAL TO THE SUM OF
C                          THE ELEMENTS IN NRTOT.  IWK SHOULD NOT
C                          BE CHANGED BETWEEN CALLS TO GGTAB.
C                IER    - OUTPUT.  ERROR PARAMETER.
C                         WARNING ERROR
C                           IER = 66 INDICATES THAT THE INPUT VALUES
C                             OF NRTOT AND/OR NCTOT ARE SUCH THAT
C                             THE PROBABILITY DISTRIBUTION OF TABLES
C                             IS DEGENERATE, THAT IS, ONLY ONE SUCH
C                             TABLE IS POSSIBLE.
C                         TERMINAL ERROR
C                           IER = 129  INDICATES AN ERROR IN THE
C                             SPECIFICATION OF NRTOT AND/OR NCTOT.
C                             THE ELEMENTS OF NRTOT AND NCTOT MUST
C                             BE NONNEGATIVE AND MUST SUM TO THE
C                             SAME QUANTITY.
C
C   REQD. IMSL ROUTINES - GGUBFS,UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE GGTAB (DSEED,NROW,NCOL,NRTOT,NCTOT,IND,IIT,ITAB,IWK,
     *                   IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            NROW,NCOL,IND,IIT,IER,NRTOT(NROW),NCTOT(NCOL),
     *                   ITAB(IIT,NCOL),IWK(1)
      DOUBLE PRECISION   DSEED
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,ICZ,IRZ,ITEMPR,IWKK,J,K,L,NC1,NCT,NR1,NRT,
     *                   NRT1
C                                  SPECIFICATIONS FOR FUNCTIONS
      REAL               GGUBFS
      DATA               NRT/0/
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
      IF (IND.NE.0) GO TO 25
C                                  CHECK VALIDITY OF NRTOT AND NCTOT.
      IRZ = 0
      NRT = 0
      DO 5 I=1,NROW
         NRT = NRT+NRTOT(I)
         IF (NRTOT(I).GT.0) GO TO 5
         IRZ = IRZ+1
         IF (NRTOT(I).LT.0) IER = 129
    5 CONTINUE
      ICZ = 0
      NCT = 0
      DO 10 I=1,NCOL
         NCT = NCT+NCTOT(I)
         IF (NCTOT(I).GT.0) GO TO 10
         ICZ = ICZ+1
         IF (NCTOT(I).LT.0) IER = 129
   10 CONTINUE
      IF (NRT.NE.NCT) IER = 129
      NR1 = NROW-1
      NC1 = NCOL-1
      IF (IER.NE.0) GO TO 9000
      IF (IRZ.GE.NR1 .OR. ICZ.GE.NC1) IER = 66
C                                  FILL IWK VECTOR.
      K = 0
      DO 20 I=1,NCOL
         L = NCTOT(I)
         IF (L.EQ.0) GO TO 20
         DO 15 J=1,L
            K = K+1
            IWK(K) = I
   15    CONTINUE
   20 CONTINUE
C                                  RANDOMLY PERMUTE IWK.
   25 K = NRT
      NRT1 = NRT-1
      IF (NRT1.LT.1) GO TO 35
      DO 30 I=1,NRT1
         J = 1+GGUBFS(DSEED)*K
         ITEMPR = IWK(K)
         IWK(K) = IWK(J)
         IWK(J) = ITEMPR
         K = K-1
   30 CONTINUE
C                                  FILL ITAB MATRIX.
   35 K = 0
      DO 50 I=1,NROW
         DO 40 J=1,NCOL
            ITAB(I,J) = 0
   40    CONTINUE
         L = NRTOT(I)
         IF (L.EQ.0) GO TO 50
         DO 45 J=1,L
            K = K+1
            IWKK = IWK(K)
            ITAB(I,IWKK) = ITAB(I,IWKK)+1
   45    CONTINUE
   50 CONTINUE
      IF (IER.EQ.0) GO TO 9005
 9000 CONTINUE
      CALL UERTST(IER,'GGTAB ')
 9005 RETURN
      END
