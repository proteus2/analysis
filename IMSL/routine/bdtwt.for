C   IMSL ROUTINE NAME   - BDTWT
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1982
C
C   PURPOSE             - COMPUTATIONS OF A TWO-WAY FREQUENCY TABLE
C
C   USAGE               - CALL BDTWT(N1,N2,K,ITAB,VECVAL,IVEC,VARVAL,
C                                    IVAR,MATFRQ,IM,IRTOT,ICTOT,IALTOT,
C                                    CHISQ,P,IER)
C
C   ARGUMENTS    N1     - INPUT NUMBER OF DISTINCT VALUES OF THE
C                           FIRST VARIABLE.
C                N2     - INPUT NUMBER OF DISTINCT VALUES OF THE
C                           SECOND VARIABLE.
C                K      - INPUT NUMBER OF DISTINCT COMBINATIONS OF THE
C                           TWO VARIABLES (LESS THAN OR EQUAL TO N1*N2)
C                ITAB   - INPUT VECTOR OF LENGTH K CONTAINING THE
C                           FREQUENCIES FOR THE CELLS.
C                VECVAL - INPUT 2 BY K ARRAY CONTAINING THE UNIQUE
C                           VALUES OF THE DATA VECTORS SORTED BY ROWS.
C                IVEC   - INPUT  ROW DIMENSION OF VECVAL EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM.
C                VARVAL - INPUT NMAX BY 2 ARRAY CONTAINING THE DISTINCT
C                           VALUES OF THE INDIVIDUAL VARIABLES. EACH
C                           COLUMN OF VARVAL MUST BE SORTED IN
C                           ASCENDING ORDER.
C                           NMAX IS THE MAXIMUM OF N1 AND N2.
C                IVAR   - INPUT. ROW DIMENSION OF VARVAL EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM.
C                MATFRQ - OUTPUT N1 BY N2 MATRIX OF FREQUENCIES FOR THE
C                           TWO VARIABLES.
C                IM     - INPUT ROW DIMENSION OF MATFRQ EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM.
C                IRTOT  - OUTPUT VECTOR OF LENGTH N1 CONTAINING THE
C                           TOTALS OF THE ROWS OF MATFRQ.
C                ICTOT  - OUTPUT VECTOR OF LENGTH N2 CONTAINING THE
C                           TOTALS OF THE COLUMNS OF MATFRQ.
C                IALTOT - OUTPUT. THE TOTAL NUMBER OF OBSERVATIONS.
C                CHISQ  - INPUT/OUTPUT.
C                           IF CHISQ IS NONNEGATIVE ON INPUT, PEARSONS
C                              CHI-SQUARE IS COMPUTED AND RETURNED
C                              IN CHISQ.
C                           IF CHISQ IS NEGATIVE ON INPUT, THE CHI-
C                              SQUARE COMPUTATIONS ARE NOT PERFORMED
C                              AND CHISQ IS UNCHANGED.
C                P      - OUTPUT APPROXIMATE PROBABILITY OF OBTAINING A
C                           LARGER CHI-SQUARE STATISTIC, IF THERE IS NO
C                           ASSOCIATION BETWEEN THE ROWS AND THE
C                           COLUMNS.
C                IER    - ERROR PARAMETER (OUTPUT)
C                          TERMINAL ERROR
C                             IER=129 INDICATES THAT N1, N2, OR K IS
C                               LESS THAN ONE OR K IS GREATER THAN
C                               N1*N2.
C                             IER=130 INDICATES THAT SOME ROW OR COLUMN
C                               TOTAL IS ZERO. THE TABLE IN MATFRQ IS
C                               CORRECT BUT CHISQ IS NOT COMPUTED.
C                             IER=131 INDICATES THAT N1 OR N2 IS LESS
C                               THAN TWO. THE TABLE IN MATFRQ IS
C                               CORRECT BUT CHISQ IS NOT COMPUTED.
C                             IER=132 INDICATES THAT THE TOTAL OF THE
C                               NUMBERS IN ITAB IS LARGER THAN THE
C                               LARGEST REPRESENTABLE INTEGER.
C                           WARNING ERROR
C                             IER=34 INDICATES THAT 20 PERCENT OF
C                               EXPECTED VALUES ARE LESS THAN 5.
C                             IER=35 INDICATES THAT DEGREES OF FREEDOM
C                               ARE GREATER THAN 30. THE USER SHOULD
C                               CONSIDER USING EXACT MEAN, STANDARD
C                               DEVIATION, AND NORMAL DISTRIBUTION
C                               FUNCTION (ALSO INDICATES THAT IER=34
C                               CONDITION IS MET).
C                             IER=36 INDICATES THAT SOME EXPECTED VALUE
C                               IS LESS THAN 2 AND STATISTICS MAY NOT
C                               BE VALID (ALSO INDICATES THAT IER=34
C                               CONDITION IS MET).
C                             IER=37 INDICATES THAT SOME EXPECTED VALUE
C                               IS LESS THAN 1 AND STATISTICS MAY NOT
C                               BE VALID.
C                             IER=38 INDICATES THAT DEGREES OF FREEDOM
C                               FOR THE CHI-SQUARE STATISTIC ARE
C                               GREATER THAN 200,000. P IS SET TO -1.
C
C   REQD. IMSL ROUTINES - H32/MDCH,MDNOR,MERRC=ERFC,MGAMAD=DGAMMA,
C                           UERTST,UGETIO
C                       - H36,H48,H60/MDCH,MDNOR,MERRC=ERFC,MGAMA=GAMMA,
C                           UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE BDTWT (N1,N2,K,ITAB,VECVAL,IVEC,VARVAL,IVAR,MATFRQ,IM,
     *                   IRTOT,ICTOT,IALTOT,CHISQ,P,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N1,N2,K,IVEC,IVAR,IM,IALTOT,IER,ITAB(1),
     *                   MATFRQ(IM,1),IRTOT(1),ICTOT(1)
      REAL               VECVAL(IVEC,1),VARVAL(IVAR,1),CHISQ,P
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,IC,ICNT,ITWO,J,JER,LARGE
      REAL               D,DF,E,RT,TI,XALTOT,XLARGE
      REAL               SCC,SDF,SP
      DOUBLE PRECISION   CC
      DATA               LARGE/Z7FFFFFFF/
C                                  FIRST EXECUTABLE STATEMENT
      XLARGE = LARGE
      IER = 0
C                                  CHECK INPUT SPECIFICATIONS
      IF (N1.GT.0 .AND. N2.GT.0 .AND. K.GT.0 .AND. K.LE.N1*N2) GO TO 5
      IER = 129
      GO TO 9000
C                                  CHECK FOR OVERFLOW OF IALTOT
    5 XALTOT = 0.0
      DO 10 I=1,K
         XALTOT = XALTOT+ITAB(I)
   10 CONTINUE
      IF (XALTOT.LT.XLARGE) GO TO 15
      IER = 132
      GO TO 9000
   15 DO 25 I=1,N1
         DO 20 J=1,N2
            MATFRQ(I,J) = 0
   20    CONTINUE
   25 CONTINUE
      IC = 1
      DO 35 I=1,N1
         DO 30 J=1,N2
            IF (VECVAL(1,IC).NE.VARVAL(I,1)) GO TO 35
            IF (VECVAL(2,IC).NE.VARVAL(J,2)) GO TO 30
            MATFRQ(I,J) = ITAB(IC)
            IC = IC+1
            IF (IC.GT.K) GO TO 40
   30    CONTINUE
   35 CONTINUE
C                                  COMPUTE MARGINAL TOTALS
   40 DO 45 J=1,N2
   45 ICTOT(J) = 0
      IALTOT = 0
      DO 55 I=1,N1
         IRTOT(I) = 0
         DO 50 J=1,N2
            IRTOT(I) = IRTOT(I)+MATFRQ(I,J)
            ICTOT(J) = ICTOT(J)+MATFRQ(I,J)
            IALTOT = IALTOT+MATFRQ(I,J)
   50    CONTINUE
   55 CONTINUE
      IF (CHISQ.LT.0.0) GO TO 9005
      IF (N1.GT.1 .AND. N2.GT.1) GO TO 60
      IER = 131
      CHISQ = 0.0
      GO TO 9000
C                                  COMPUTE CHI-SQUARE
   60 DO 65 I=1,N1
         IF (IRTOT(I).LE.0) GO TO 75
   65 CONTINUE
      DO 70 I=1,N2
         IF (ICTOT(I).LE.0) GO TO 75
   70 CONTINUE
      GO TO 80
   75 IER = 130
      GO TO 9000
   80 TI = 1.0/XALTOT
      CC = 0.D0
      ICNT = 0
      ITWO = 0
      DO 90 I=1,N1
         RT = IRTOT(I)
         DO 85 J=1,N2
            E = RT*ICTOT(J)*TI
            IF (E.LT.1.0) IER = 37
            IF (E.LT.5.0) ICNT = ICNT+1
            IF (E.LT.2.0) ITWO = 1
            D = MATFRQ(I,J)-E
            CC = CC+DBLE(D)**2/DBLE(E)
   85    CONTINUE
   90 CONTINUE
      IF (ICNT.GT.IALTOT/5 .AND. IER.EQ.0) IER = 34
      IF ((N1-1)*(N2-1).GT.30 .AND. IER.EQ.34) IER = 35
      IF (ITWO.GT.0 .AND. IER.EQ.34) IER = 36
      CHISQ = CC
      DF = (N1-1)*(N2-1)
      IF (DF.GT..5 .AND. DF.LT.2.E5) GO TO 95
      IER=38
      P = -1.
      GO TO 9000
   95 SCC = CC
      SDF = DF
      CALL MDCH(SCC,SDF,SP,JER)
      P = 1.-SP
      IF (IER.LE.0) GO TO 9005
 9000 CONTINUE
      CALL UERTST(IER,'BDTWT ')
 9005 RETURN
      END
