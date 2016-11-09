C   IMSL ROUTINE NAME   - MDTN
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - NON-CENTRAL T PROBABILITY DISTRIBUTION
C                           FUNCTION
C
C   USAGE               - CALL MDTN (TVAL,IDF,D,P,IER)
C
C   ARGUMENTS    TVAL   - INPUT VALUE TO WHICH INTEGRATION IS PERFORMED.
C                IDF    - INPUT. INTEGER DEGREES OF FREEDOM.
C                D      - INPUT. NON-CENTRALITY PARAMETER.
C                P      - OUTPUT.  PROBABILITY THAT A RANDOM VARIABLE
C                           DISTRIBUTED AS T WITH NON-CENTRALITY
C                           PARAMETER D IS LESS THAN OR EQUAL TO TVAL.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER = 129, INDICATES DEGREES OF FREEDOM OUT
C                             OF RANGE. IDF MUST BE GREATER THAN ZERO.
C
C   REQD. IMSL ROUTINES - MDNOR,MDTNF,MERRC=ERFC,UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE MDTN   (TVAL,IDF,D,P,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            IDF,IER
      REAL               TVAL,D,P
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I1,IDFM2,L
      REAL               DF,A,B,SB,DA,DSB,DASB,P1,F2,F1,SUM,A1,A2,
     *                   A3,EPS,AZ,FZ,FKM1
      DATA               A1/.3989423/,A2/.1591549/
      DATA               A3/2.506628/,EPS/.000001/
C                                  FIRST EXECUTABLE STATEMENT
      IF (IDF .LE. 0) GO TO 60
      IER = 0
      DF = IDF
      I1 = IDF-2*(IDF/2)
      A = TVAL/SQRT(DF)
      B = DF/(DF+TVAL*TVAL)
      SB = SQRT(B)
      DA = D*A
      DSB = D*SB
      DASB = A*DSB
      CALL MDNOR(DASB,P1)
      F2 = A*SB*EXP(-.5*DSB*DSB) * P1*A1
      F1 = B*(DA*F2+A*A2*EXP(-.5*D*D))
      SUM = 0.0
      IF (IDF .EQ. 1) GO TO 50
      IF (I1) 5,5,10
    5 SUM = F2
      GO TO 15
   10 SUM = F1
   15 IF (IDF .LT. 4) GO TO 45
   20 IDFM2 = IDF-2
      AZ = 1.0
      FZ = 2.0
      DO 40 L=2,IDFM2,2
         FKM1 = FZ-1.0
         F2 = B*(DA*AZ*F1+F2)*FKM1/FZ
         AZ = 1.0/(AZ*FKM1)
         F1 = B*(DA*AZ*F2+F1)*FZ/(FZ+1.0)
         IF (I1) 25,25,30
   25    SUM = SUM+F2
         GO TO 35
   30    SUM = SUM+F1
   35    AZ = 1.0/(AZ*FZ)
   40 FZ = FZ+2.0
   45 IF (I1) 55,55,50
   50 CALL MDNOR(-DSB,P1)
      CALL MDTNF(DSB,A,EPS,P)
      P = P1+2.0*(P+SUM)
      GO TO 9005
   55 CALL MDNOR(-D,P1)
      P = P1+SUM*A3
      GO TO 9005
   60 IER = 129
 9000 CONTINUE
      CALL UERTST (IER,6HMDTN  )
 9005 CONTINUE
      RETURN
      END
