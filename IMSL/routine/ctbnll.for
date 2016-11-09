C   IMSL ROUTINE NAME   - CTBNLL
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1982
C
C   PURPOSE             - NUCLEUS CALLED ONLY BY IMSL SUBROUTINE CBNRHO
C
C   REQD. IMSL ROUTINES - SINGLE/MDBNOR,MDNOR,MDTNF,MERRC=ERFC,UERTST,
C                           UGETIO
C                       - DOUBLE/MDBNOR,MDNOR,MDNORD,MDTNF,MERRC=ERFC,
C                           MERRCD=DERFC,UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      REAL FUNCTION CTBNLL (RHO,AP,ZER,IRCV,IA,IB)
C
C                                 SPECIFICATIONS FOR ARGUMENTS
      INTEGER            IRCV(1),IA,IB
      REAL               RHO,AP(IA,IB,2),ZER(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,IC,ICM1,ICP1,IER,IR,IRM1,IRP1,I1,J,J1
      REAL               RDELP,X
      DOUBLE PRECISION   SUM
      DATA               RDELP/Z3C100000/
C                                  FIRST EXECUTABLE STATEMENT
      IR = IRCV(1)
      IC = IRCV(2)
      IRP1 = IR + 1
      ICP1 = IC + 1
      IRM1 = IR - 1
      ICM1 = IC - 1
      IER = 0
C                                       *
C                                  THE P 'S ARE CUMULATIVE BIVARIATE
C                                  NORMAL PROBABILITIES
C                                   *
C                                  P (INF,ALFA ) USING MDNOR/MDNORD
C                                             I
C
      DO 45 I = 1,IRM1
         I1 = IR - I + 1
         X = AP(I,ICP1,2)
         IF (X .GT. 0.0) GO TO 30
      CALL MDNOR(X,AP(I1,IC,2))
         GO TO 35
   30 CALL MDNOR(-X,AP(I1,IC,2))
      AP(I1,IC,2) = 1. - AP(I1,IC,2)
C
C                                   *
C                                  P (BETA ,ALFA ) USING MDBNOR
C                                         J     I
C
  35     DO 40 J = 1,ICM1
            CALL MDBNOR(X,AP(IRP1,J,2),RHO,AP(I1,J,2),IER)
      CALL UGETIO(1,NIN,NOUT)
  40     CONTINUE
  45  CONTINUE
C
C                                   *
C                                  P (BETA ,INF) USING MDNOR/MDNORD
C                                         J
C
      DO 55 J = 1,ICM1
         X = AP(IRP1,J,2)
         IF (X .GT. 0.0) GO TO 50
      CALL MDNOR(X,AP(1,J,2))
         GO TO 55
   50 CALL MDNOR(-X,AP(1,J,2))
      AP(1,J,2) = 1. - AP(1,J,2)
   55 CONTINUE
C                                  P(BETA ,ALFA )
C                                        J     I
      AP(1,IC,2) = 1.0
      DO 65 J = 1,ICM1
         J1 = IC + 1 - J
         DO 60 I = 1,IRM1
            AP(I,J1,2) = AP(I,J1,2) + AP(I+1,J1-1,2) - AP(I,J1-1,2)
     *                  -AP(I+1,J1,2)
  60     CONTINUE
         AP(IR,J1,2) = AP(IR,J1,2) - AP(IR,J1-1,2)
  65  CONTINUE
      DO 70 I = 1,IRM1
         AP(I,1,2) = AP(I,1,2) - AP(I+1,1,2)
  70  CONTINUE
      SUM = 0.D0
      DO 85 I = 1,IR
         DO 80 J = 1,IC
            IF (AP(I,J,2) .GE. RDELP) GO TO 75
            AP(I,J,2) = RDELP
C                                  WARNING(WITH FIX) - A BIVARIATE
C                                    NORMAL PROBABILITY TOO SMALL
            IER = 65
  75        SUM = SUM + AP(I,J,1)*ALOG(AP(I,J,2))
  80     CONTINUE
  85  CONTINUE
      CTBNLL =-SUM
 9000 CONTINUE
      IF (IER .EQ. 65) ZER(1)=1
 9005 RETURN
      END
