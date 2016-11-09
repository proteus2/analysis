C   IMSL ROUTINE NAME   - MDCHN
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JULY 1, 1983
C
C   PURPOSE             - NON-CENTRAL CHI-SQUARED PROBABILITY DISTRI-
C                           BUTION FUNCTION
C
C   USAGE               - CALL MDCHN (CS,DF,PNONC,P,IER)
C
C   ARGUMENTS    CS     - INPUT VALUE FOR WHICH THE PROBABILITY IS
C                           COMPUTED.  CS MUST BE GREATER THAN OR EQUAL
C                           TO ZERO.
C                DF     - INPUT VALUE CONTAINING NUMBER OF DEGREES OF
C                           FREEDOM OF THE NON-CENTRAL CHI-SQUARED
C                           DISTRIBUTION.  DF MUST BE GREATER THAN
C                           OR EQUAL TO .5 AND LESS THAN OR EQUAL
C                           TO 200,000. SEE REMARKS.
C                PNONC  - INPUT VALUE OF THE NON-CENTRALITY PARAMETER.
C                           PNONC MUST BE GREATER THAN OR EQUAL TO ZERO.
C                           SEE REMARKS.
C                P      - OUTPUT VALUE CONTAINING THE PROBABILITY THAT
C                           A RANDOM VARIABLE WHICH FOLLOWS THE NON-
C                           CENTRAL CHI-SQUARED DISTRIBUTION WITH NON-
C                           CENTRALITY PARAMETER PNONC AND CONTINUOUS
C                           DEGREES OF FREEDOM DF, IS LESS THAN OR
C                           EQUAL TO CS.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER = 129 INDICATES THAT CS OR DF OR PNONC
C                             WAS SPECIFIED INCORRECTLY. SEE REMARKS.
C                           IER = 130 INDICATES THAT CONVERGENCE WAS NOT
C                             OBTAINED. SEE REMARKS.
C                         WARNING ERROR
C                           IER = 34 INDICATES THAT UNDERFLOW WOULD HAVE
C                             OCCURRED IN IMSL ROUTINE MDCH. P IS SET TO
C                             ZERO.
C
C   REQD. IMSL ROUTINES - H32/MDCH,MDNOR,MERRC=ERFC,MGAMAD=DGAMMA,
C                           MLGAMA=ALGAMA,UERTST,UGETIO
C                       - H36,H48,H60/MDCH,MDNOR,MERRC=ERFC,MGAMA=GAMMA,
C                           MLGAMA=ALGAMA,UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS  1.  THE SUM, DF + PNONC, MUST BE LESS THAN OR EQUAL TO
C                200000. OTHERWISE, A TERMINAL ERROR WILL RESULT.
C            2.  THIS SUBROUTINE SUMS TERMS OF AN INFINITE SERIES.
C                SUMMING TERMINATES WHEN EITHER THE CURRENT TERM IS
C                LESS THAN EPS TIMES THE SUM OR NTIRED TERMS HAVE
C                BEEN ADDED. IN THE LATTER CASE A TERMINAL ERROR IS
C                ISSUED. NTIRED AND EPS ARE SET IN A DATA STATEMENT
C                TO 1000 AND 1.0E-5.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE MDCHN  (CS,DF,PNONC,P,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
C     INTEGER            IER
      REAL               CS,DF,PNONC,P
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            ICENT,ITERB,I,ITERF,NTIRED,JER
      REAL               CENTWT,CENTAJ,FACT,SUM,EPS,XNONC,XCENT,
     *                   CHID2,DFD2,DG,PCENT,SUMADJ,ADJ,WT,ACENTW,
     *                   ACENTA,AFACT,PTERM,TERM
      DATA               NTIRED/1000/,EPS/1.0E-5/
      DG(I) = DF+2.0*FLOAT(I)
C                                  ERROR CHECK CS, DF, AND PNONC
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
      IF ((CS.LT.0.0) .OR. (DF.LT.0.5.OR.DF.GT.200000.0) .OR.
     *    (PNONC.LT.0.0)) GO TO 45
C                                  HANDLE CASE WHEN NON-CENTRALITY
C                                  PARAMETER IS (ESSENTIALLY) ZERO
      IF (PNONC.GT.1.0E-10) GO TO 5
      CALL MDCH (CS,DF,P,JER)
      IF (JER.EQ.34) GO TO 50
      RETURN
C                                  THE FOLLOWING CODE CALCULATES THE
C                                  WEIGHT, CHI-SQUARE, AND ADJUSTMENT
C                                  TERM FOR THE CENTRAL TERM IN THE
C                                  INFINITE SERIES. THE CENTRAL TERM IS
C                                  THE ONE IN WHICH THE POISSON WEIGHT
C                                  IS GREATEST. THE ADJUSTMENT TERM IS
C                                  THE AMOUNT THAT MUST BE SUBTRACTED
C                                  FROM THE CHI-SQUARE TO MOVE UP TWO
C                                  DF.
    5 IF (CS.GT.0.0) GO TO 10
      P = 0.0
      GO TO 9005
   10 XNONC = PNONC/2.0
      ICENT = XNONC
      IF (ICENT.EQ.0) ICENT = 1
      XCENT = ICENT
      CHID2 = CS/2.0
C                                  CALCULATE CENTRAL WEIGHT TERM
      AFACT = ALGAMA(FLOAT(ICENT+1))
      ACENTW = -XNONC+ICENT*ALOG(XNONC)-AFACT
      CENTWT = EXP(ACENTW)
C                                  CALCULATE CENTRAL CHI-SQUARE
      CALL MDCH (CS,DG(ICENT),PCENT,JER)
      IF (JER.GT.128) GO TO 45
C                                  CALCULATE CENTRAL ADJUSTMENT TERM
      DFD2 = DG(ICENT)/2.0
      AFACT = ALGAMA((1.0+DFD2))
      ACENTA = DFD2*ALOG(CHID2)-CHID2-AFACT
      CENTAJ = EXP(ACENTA)
      SUM = CENTWT*PCENT
C                                  WE NOW SUM BACKWARDS FROM THE CENTRAL
C                                  TERM TOWARDS ZERO. WE QUIT WHENEVER
C                                  EITHER (1) WE REACH THE ZERO TERM OR
C                                  (2) THE TERM GETS SMALL RELATIVE TO
C                                  THE SUM OR (3) WE SUM MORE THAN
C                                  NTIRED TERMS.
      ITERB = 0
      SUMADJ = 0.0
      ADJ = CENTAJ
      WT = CENTWT
      I = ICENT
      GO TO 20
   15 IF (ITERB.GT.NTIRED .OR. SUM.LT.1.0E-20 .OR. TERM.LT.EPS*SUM .OR.
     *    I.EQ.0) GO TO 25
   20 DFD2 = DG(I)/2.0
C                                  ADJUST CHI-SQUARE FOR TWO FEWER
C                                  DEGREES OF FREEDOM. THE ADJUSTED
C                                  VALUE ENDS UP IN PTERM.
      ADJ = ADJ*DFD2/CHID2
      SUMADJ = SUMADJ+ADJ
      PTERM = PCENT+SUMADJ
C                                  NOW ADJUST POISSON WEIGHT FOR J
C                                  DECREASED BY ONE
      WT = WT*(I/XNONC)
      TERM = WT*PTERM
      SUM = SUM+TERM
      I = I-1
      ITERB = ITERB+1
      GO TO 15
C                                  WE NOW SUM FORWARD FROM THE CENTRAL
C                                  TERM TOWARDS INFINITY. WE QUIT WHEN
C                                  EITHER (1) THE TERM GETS SMALL RELA-
C                                  TIVE TO THE SUM OR (2) WE SUM MORE
C                                  THAN NTIRED TERMS.
   25 ITERF = 0
      SUMADJ = CENTAJ
      ADJ = CENTAJ
      WT = CENTWT
      I = ICENT
      GO TO 35
   30 IF (ITERF.GT.NTIRED .OR. SUM.LT.1.0E-20 .OR. TERM.LT.EPS*SUM)
     *    GO TO 40
C                                  UPDATE WEIGHTS FOR NEXT HIGHER J
   35 WT = WT*(XNONC/(I+1))
C                                  NOW CALCULATE PTERM AND ADD TERM TO
C                                  SUM
      PTERM = PCENT-SUMADJ
      TERM = WT*PTERM
      SUM = SUM+TERM
C                                  NOW UPDATE ADJUSTMENT TERM FOR DF FOR
C                                  NEXT ITERATION
      I = I+1
      DFD2 = DG(I)/2.0
      ADJ = ADJ*CHID2/DFD2
      SUMADJ = SUMADJ+ADJ
      ITERF = ITERF+1
      GO TO 30
   40 P = SUM
      IF (ITERB.LE.NTIRED .AND. ITERF.LE.NTIRED) GO TO 9005
      IER = 130
      GO TO 9000
   45 IER = 129
      GO TO 9000
   50 IER = 34
      P = 0.0
 9000 CONTINUE
      CALL UERTST (IER,6HMDCHN )
 9005 RETURN
      END
