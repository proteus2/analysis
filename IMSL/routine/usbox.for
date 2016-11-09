C   IMSL ROUTINE NAME   - USBOX
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1981
C
C   PURPOSE             - PRINT A BOXPLOT (K SAMPLES)
C
C   USAGE               - CALL USBOX (X,K,NI,MAXL,IER)
C
C   ARGUMENTS    X      - INPUT VECTOR OF LENGTH NI(1)+NI(2)+...+NI(K)
C                           CONTAINING IN THE FIRST NI(1) POSITIONS
C                           THE OBSERVATIONS FOR THE FIRST SAMPLE, IN
C                           THE NEXT NI(2) POSITIONS, THE OBSERVATIONS
C                           FOR THE SECOND SAMPLE, AND SO ON.
C                         ON OUTPUT, EACH SAMPLE WILL BE SORTED.
C                K      - THE NUMBER OF SAMPLES. (INPUT)
C                NI     - VECTOR OF LENGTH K. (INPUT) NI(I) IS THE
C                         NUMBER OF OBSERVATIONS IN THE I-TH SAMPLE.
C                MAXL   - MAXIMUM DISPLAY WIDTH. (INPUT)
C                           MAXL MUST BE 80 OR 129.
C                IER    - ERROR PARAMETER (OUTPUT)
C                         WARNING ERROR WITH FIX
C                           IER=65 INDICATES MAXL DID NOT EQUAL 80 OR
C                             129. MAXL IS SET TO 80.
C                         TERMINAL ERROR
C                           IER=129 INDICATES NI(I) IS LESS THAN 1.
C                           IER=130 INDICATES K IS LESS THAN 1.
C                           IER=131 INDICATES ALL VALUES IN X ARE
C                             CONSTANT.
C
C   REQD. IMSL ROUTINES - BDLTV,UERTST,UGETIO,USBOX1,VSRTA
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS      IF XMAT IS A USER DATA MATRIX IN WHICH COLUMNS
C                REPRESENT SAMPLES AND IF EACH SAMPLE IS COMPLETE (HAS
C                THE SAME NUMBER OF OBSERVATIONS), THE USER CAN
C                EQUIVALENCE X AND XMAT AND ENTER USBOX WITHOUT
C                EXPLICITLY MOVING XMAT INTO X. HOWEVER, THE ROWS OF
C                XMAT WILL NOT BE PRESERVED.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE USBOX (X,K,NI,MAXL,IER)
C
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            K,NI(1),MAXL,IER
      REAL               X(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,ILL,ISTRT,JER,N,NSTRT,NUM
      REAL               FIVNO(5),SC,XMAX,XMIN
      DATA               NUM /5/
C                                  FIRST EXECUTABLE STATEMENT
      IER = 130
      IF (K.LT.1) GO TO 9000
      IER = 129
      N = 0
      DO 5 I=1,K
         IF (NI(I).LE.0) GO TO 9000
         N = N+NI(I)
    5 CONTINUE
      IER = 0
      IF (MAXL.NE.129 .AND. MAXL.NE.80) IER = 65
      IF (MAXL.NE.129) MAXL = 80
C                                  DETERMINE OVERALL MAX, MIN AND
C                                  SCALE FACTOR.
      XMIN = X(1)
      XMAX = X(1)
      DO 10 I=2,N
         IF (X(I).LT.XMIN) XMIN = X(I)
         IF (X(I).GT.XMAX) XMAX = X(I)
   10 CONTINUE
      IF (XMAX.GT.XMIN) GO TO 15
      IER = 131
      GO TO 9000
   15 CALL UGETIO(1,NIN,NOUT)
      IF (MAXL.EQ.80) WRITE (NOUT,30) XMIN, XMAX
      IF (MAXL.EQ.129) WRITE (NOUT,35) XMIN, XMAX
      SC = (XMAX-XMIN)/(MAXL-3)
C                                  MAIN LOOP TO SORT SAMPLES, GET FIVE-
C                                  NUMBER SUMMARY, AND PRINT BOXPLOTS
      NSTRT = 1
      DO 20 I=1,K
         CALL VSRTA(X(NSTRT),NI(I))
         CALL BDLTV(X(NSTRT),NI(I),NUM,FIVNO,JER)
         ISTRT = (X(NSTRT)-XMIN)*.9999/SC+1.
         ILL = NSTRT+NI(I)-1
         CALL USBOX1(X(NSTRT),NI(I),FIVNO,SC,ISTRT,MAXL)
         WRITE (NOUT,25)
         NSTRT = NSTRT+NI(I)
   20 CONTINUE
   25 FORMAT (///)
   30 FORMAT (//, 1X, E11.4, 56X, E11.4, //)
   35 FORMAT (//, 1X, E11.4, 105X, E11.4, //)
      IF (IER.EQ.0) GO TO 9005
 9000 CONTINUE
      CALL UERTST(IER,6HUSBOX )
 9005 RETURN
      END
