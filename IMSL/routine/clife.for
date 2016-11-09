C   IMSL ROUTINE NAME   - CLIFE
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1982
C
C   PURPOSE             - LIFE TABLE ANALYSIS.
C
C   USAGE               - CALL CLIFE (ITYPE,N,IRAD,AGE,A,IPOP,IDTH,L,
C                           IDL,DR,Q,SQ,POI,SPOI,E,SE,TI,IER)
C
C   ARGUMENTS    ITYPE  - INPUT, TYPE OF LIFE TABLE. ITYPE=0 INDICATES
C                           A CURRENT TABLE. ITYPE=1 INDICATES A COHORT
C                           TABLE.
C                N      - INPUT, NUMBER OF AGE CLASSES.
C                IRAD   - INPUT, REQUIRED ONLY FOR A CURRENT LIFE TABLE.
C                           IRAD IS THE POPULATION SIZE FOR THE FIRST
C                           AGE INTERVAL IN A CURRENT LIFE TABLE, AND
C                           ALL THE ELEMENTS IN L ARE BASED ON THIS
C                           NUMBER.  THE VALUE IS SOMEWHAT ARBITRARY
C                           AND 10000 MAY BE A REASONABLE NUMBER.
C                AGE   - INPUT OR INPUT/OUTPUT VECTOR OF LENGTH N
C                           CONTAINING THE LOWEST AGE IN EACH AGE
C                           INTERVAL.  ORDINARILY AGE(1)=0.0.
C                           FOR COMPLETE TABLES, IN WHICH THE AGE
C                           INTERVALS ARE ALL OF EQUAL LENGTH,
C                           AGE(1) CAN BE SPECIFIED AS A
C                           NEGATIVE QUANTITY EQUAL IN ABSOLUTE VALUE
C                           TO THE CONSTANT LENGTH OF THE INTERVALS,
C                           AND THE OTHER VALUES OF AGE DO NOT NEED
C                           TO BE SPECIFIED.  IN THIS CASE, ON OUTPUT,
C                           AGE WILL CONTAIN THE LOWEST AGE IN EACH
C                           AGE INTERVAL, BEGINNING WITH 0.0.
C                A      - INPUT VECTOR OF LENGTH N CONTAINING THE
C                           FRACTIONS OF LAST AGE INTERVALS OF LIFE
C                           FOR THE AGES IN VECTOR AGE.  SINCE 0.5 IS A
C                           COMMON CHOICE FOR ALL OF THESE, IF A(1) IS
C                           SET TO A NEGATIVE VALUE, CLIFE WILL IGNORE
C                           THE OTHER ELEMENTS OF A AND WILL ASSUME A
C                           CONSTANT VALUE OF 0.5.  IF A(1) IS POSITIVE,
C                           THEN ALL VALUES OF A MUST BE BETWEEN
C                           0 AND 1.
C                IPOP   - INPUT VECTOR OF LENGTH N CONTAINING THE SIZES
C                           OF THE POPULATION.  THE ELEMENTS OF IPOP
C                           MUST BE POSITIVE.  IF ITYPE=0 (CURRENT LIFE
C                           TABLE), THE SIZE IS AT THE MID-INTERVAL.
C                           IF ITYPE=1, THE SIZE IS THE COHORT
C                           POPULATION AT THE BEGINNING OF THE AGE
C                           INTERVAL AND IN THIS CASE THE ELEMENTS OF
C                           IPOP MUST BE DECREASING IN SIZE.
C                IDTH   - INPUT VECTOR OF LENGTH N CONTAINING THE
C                           NUMBER OF DEATHS IN EACH AGE INTERVAL
C                           IF ITYPE=0.  IF ITYPE=1, IDTH IS NOT USED.
C                L      - OUTPUT VECTOR OF LENGTH N CONTAINING THE
C                           NUMBERS LIVING AT THE AGES IN VECTOR AGE.
C                           (IF ITYPE=1, L CONTAINS THE SAME VALUES AS
C                           IPOP EXCEPT POSSIBLY FOR ROUNDING.)
C                IDL    - OUTPUT VECTOR OF LENGTH N CONTAINING THE
C                           NUMBERS OF DEATHS REFERENCED TO THE
C                           NUMBERS LIVING IN L.
C                DR     - OUTPUT VECTOR OF LENGTH N CONTAINING THE
C                           DEATH RATE IN EACH AGE INTERVAL. (DEFINED
C                           ONLY IF ITYPE=0.)
C                Q      - OUTPUT VECTOR OF LENGTH N CONTAINING THE
C                           PROPORTION DYING IN EACH AGE INTERVAL.
C                SQ     - OUTPUT VECTOR OF LENGTH N CONTAINING THE
C                           SAMPLE STANDARD ERRORS OF Q.
C                POI    - OUTPUT VECTOR OF LENGTH N.  THE I-TH ENTRY
C                           CONTAINS THE PROPORTION OF SURVIVORS UP
C                           TO AGE(I).
C                SPOI   - OUTPUT VECTOR OF LENGTH N CONTAINING THE
C                           SAMPLE STANDARD ERRORS OF POI.
C                E      - OUTPUT VECTOR OF LENGTH N CONTAINING THE
C                           ESTIMATED EXPECTED LIFE AT EACH AGE IN
C                           VECTOR AGE.
C                SE     - OUTPUT VECTOR OF LENGTH N CONTAINING THE
C                           ESTIMATED STANDARD ERRORS OF E.
C                TI     - OUTPUT VECTOR OF LENGTH N CONTAINING THE
C                           TOTAL NUMBER OF TIME UNITS LIVED BY ALL
C                           OF THE POPULATION IN EACH INTERVAL.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR.
C                           IER=129 IRAD OR FIRST ELEMENT OF IPOP IS
C                             NEGATIVE.
C                           IER=130 AGE IS NOT IN INCREASING ORDER.
C                           IER=131 SOME ELEMENT OF IPOP IS NEGATIVE,
C                             OR THE ELEMENTS OF IPOP ARE NOT
C                             DECREASING.
C                           IER=132 IDTH IS SPECIFIED INCORRECTLY.  (IF
C                              ITYPE=0, THE ELEMENTS OF IDTH MUST BE
C                              POSITIVE AND LESS THAN THE ELEMENTS OF
C                              IPOP.)
C                           IER=133 SOME ELEMENT OF A (OTHER THAN A(1))
C                             IS NOT IN THE EXCLUSIVE INTERVAL 0.0 TO
C                             1.0.
C                           IER=134 N IS LESS THAN 2.
C
C   REQD. IMSL ROUTINES - UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE CLIFE  (ITYPE,N,IRAD,AGE,A,IPOP,IDTH,L,IDL,DR,Q,SQ,POI,
     *                   SPOI,E,SE,TI,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            ITYPE,N,IRAD,IER,IPOP(N),IDTH(N),L(N),IDL(N)
      REAL               AGE(N),A(N),DR(N),Q(N),SQ(N),POI(N),SPOI(N),
     *                   E(N),SE(N),TI(N)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,IRAD1,J,N1
      REAL               POP,SVAR,TYEARS,XLJ,XLJ2,XN
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
      IF (N.GE.2) GO TO 5
      IER = 134
      GO TO 9000
    5 N1 = N-1
      IF (ITYPE.EQ.0) IRAD1 = IRAD
      IF (ITYPE.EQ.1) IRAD1 = IPOP(1)
      IF (IRAD1.GT.0) GO TO 10
      IER = 129
      GO TO 9000
   10 IF (AGE(1).GE.0.0) GO TO 20
      AGE(2) = -AGE(1)
      AGE(1) = 0.0
      IF (N.EQ.2) GO TO 20
      DO 15 I=3,N
   15 AGE(I) = AGE(I-1)+AGE(2)
   20 IF (A(1).GT.0.0) GO TO 30
      DO 25 I=1,N
   25 A(I) = 0.5
      GO TO 40
   30 DO 35 I=1,N1
         IF (A(I).GT.0.0 .AND. A(I).LT.1.0) GO TO 35
         IER = 133
         GO TO 9000
   35 CONTINUE
   40 DO 75 I=2,N
         J = I-1
         XN = AGE(I)-AGE(J)
         IF (XN.GT.0.0) GO TO 45
         IER = 130
         GO TO 9000
   45    IF (IPOP(J).GT.0) GO TO 50
         IER = 131
         GO TO 9000
   50    IF (ITYPE.EQ.1) GO TO 60
         IF (IDTH(J).GT.0 .AND. IDTH(J).LE.IPOP(J)) GO TO 55
         IER = 132
         GO TO 9000
   55    DR(J) = FLOAT(IDTH(J))/FLOAT(IPOP(J))
C                                  COMPUTE POPULATION AT START OF
C                                  INTERVAL
         POP = FLOAT(IPOP(J))/XN+FLOAT(IDTH(J))*(1.0-A(J))
         Q(J) = FLOAT(IDTH(J))/POP
         SQ(J) = Q(J)*(1.0-Q(J))/POP
         L(J) = IRAD1
         XRAD1 = IRAD1*Q(J)
         IDL(J) = XRAD1+0.5
         GO TO 70
   60    IDL(J) = IPOP(J)-IPOP(I)
         IF (IDL(J).GT.0) GO TO 65
         IER = 131
         GO TO 9000
   65    POP = IPOP(J)
         Q(J) = FLOAT(IDL(J))/POP
         SQ(J) = Q(J)*(1.0-Q(J))/POP
         L(J) = IRAD1
         XRAD1 = IRAD1*Q(J)
   70    IRAD1 = IRAD1-IDL(J)
   75 CONTINUE
      Q(N) = 1.0
      IF (ITYPE.EQ.0) DR(N) = FLOAT(IDTH(N))/FLOAT(IPOP(N))
      L(N) = IRAD1
      IDL(N) = IRAD1
      IF (ITYPE.EQ.1) GO TO 80
      TYEARS = L(N)/DR(N)
      GO TO 85
   80 TYEARS = FLOAT(IPOP(N))*A(N)*(AGE(N)-AGE(N-1))
   85 IF (L(N).LE.0) L(N) = 1
      E(N) = TYEARS/L(N)
      SVAR = 0.0
      TI(N) = TYEARS
      DO 90 I=1,N1
         J = N-I
         XN = AGE(J+1)-AGE(J)
         XLJ = L(J)
         XLJ2 = XLJ*XLJ
         TI(J) = XN*(XLJ+IDL(J)*(A(J)-1.0))
         TYEARS = TYEARS+TI(J)
         E(J) = TYEARS/XLJ
         SVAR = SVAR+XLJ2*SQ(J)*(E(J+1)+XN*(1.0-A(J)))**2
         SE(J) = SVAR/XLJ2
   90 CONTINUE
      SVAR = 0.0
      XN = 1.0/L(1)
      DO 95 I=1,N1
         POI(I) = L(I)*XN
         SPOI(I) = SQRT(POI(I)*POI(I)*SVAR)
         SVAR = SVAR+SQ(I)/(1.0-Q(I))**2
         SQ(I) = SQRT(SQ(I))
         SE(I) = SQRT(SE(I))
   95 CONTINUE
      SQ(N) = 0.0
      SE(N) = 0.0
      POI(N) = L(N)*XN
      SPOI(N) = SQRT(POI(N)*POI(N)*SVAR)
      GO TO 9005
 9000 CONTINUE
      CALL UERTST(IER,'CLIFE ')
 9005 RETURN
      END
