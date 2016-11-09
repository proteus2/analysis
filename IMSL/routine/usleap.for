C   IMSL ROUTINE NAME   - USLEAP
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1982
C
C   PURPOSE             - PRINT RESULTS OF THE BEST-REGRESSIONS
C                           ANALYSIS PERFORMED BY IMSL ROUTINE RLEAP.
C
C   USAGE               - CALL USLEAP (IJOB,KZ,IXS,STAT,IXV,NVAR,IXB,
C                           BEST,IB)
C
C   ARGUMENTS    IJOB   - INPUT (RLEAP OUTPUT) OPTION AND CONTROL
C                           PARAMETER VECTOR OF LENGTH 4.
C                KZ     - NUMBER OF VARIABLES. (INPUT)
C                IXS    - INPUT (RLEAP OUTPUT) VECTOR CONTAINING THE
C                           LOCATION OF THE FIRST ELEMENT FOR EACH
C                           SUBSET SIZE IN STAT.
C                STAT   - INPUT (RLEAP OUTPUT) VECTOR CONTAINING THE
C                           CRITERION VALUES FOR EACH SUBSET CONSIDERED,
C                           IN INCREASING SUBSET SIZE ORDER.
C                IXV    - INPUT (RLEAP OUTPUT) VECTOR CONTAINING THE
C                           LOCATION OF THE FIRST ELEMENT FOR EACH
C                           SUBSET SIZE IN NVAR.
C                NVAR   - INPUT (RLEAP OUTPUT) VECTOR CONTAINING THE
C                           VARIABLE NUMBERS FOR EACH SUBSET CONSIDERED,
C                           ORDERED CORRESPONDINGLY TO STAT.
C                IXB    - INPUT (RLEAP OUTPUT) VECTOR CONTAINING THE ROW
C                           NUMBER OF THE FIRST ROW FOR EACH SUBSET SIZE
C                           IN BEST.
C                BEST   - INPUT (RLEAP OUTPUT) MATRIX CONTAINING THE
C                           RESULTS FOR THE BEST REGRESSIONS BEGINNING
C                           WITH ONE VARIABLE REGRESSIONS AND INCREASING
C                           THE SUBSET SIZE. COLUMNS 1,2,3, AND 4
C                           CONTAIN VARIABLE NUMBER, REGRESSION
C                           COEFFICIENT, F VALUE, AND TAIL AREA OF THE
C                           F DISTRIBUTION, RESPECTIVELY.
C                IB     - ROW DIMENSION OF MATRIX BEST EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM. (INPUT)
C
C   REQD. IMSL ROUTINES - NONE REQUIRED
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS      OUTPUT IS WRITTEN TO THE UNIT SPECIFIED BY IMSL
C                ROUTINE UGETIO. SEE THE UGETIO DOCUMENT FOR DETAILS.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE USLEAP (IJOB,KZ,IXS,STAT,IXV,NVAR,IXB,BEST,IB)
C
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            KZ,IB,IJOB(4),IXS(1),IXV(1),NVAR(1),IXB(1)
      REAL               STAT(1),BEST(IB,4)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,IBEG,IBIT,IC,ICNT,ICP,IEND,J,JBEG,JEND,JJ,K,
     *                   MBST,MMM,NIN,NOUT
C                                  FIRST EXECUTABLE STATEMENT
      CALL UGETIO(1,NIN,NOUT)
      IBIT = IJOB(2)
      MBST = IJOB(3)
      MMM = KZ-1
      IF (IBIT.GT.0) GO TO 5
      MMM = -IBIT
      IBIT = 1
    5 CONTINUE
      WRITE (NOUT,40)
      ICNT = 0
      DO 20 I=1,MMM
         IF (I.EQ.MMM) IXS(I+1) = IXS(I)+1
         IBEG = IXS(I)
         IEND = IXS(I+1)-1
         JBEG = IXV(I)
         JEND = JBEG+I-1
         ICNT = ICNT+6+IEND-IBEG
         IF (ICNT.LE.55) GO TO 10
         ICNT = IEND-IBEG+6
         WRITE (NOUT,40)
   10    IF (IBIT.EQ.1) WRITE (NOUT,45) I
         IF (IBIT.EQ.2) WRITE (NOUT,50) I
         IF (IBIT.EQ.3) WRITE (NOUT,55) I
         WRITE (NOUT,60)
         DO 15 J=IBEG,IEND
            WRITE (NOUT,65) STAT(J), (NVAR(K),K=JBEG,JEND)
            JBEG = JEND+1
            JEND = JBEG+I-1
   15    CONTINUE
   20 CONTINUE
      WRITE (NOUT,40)
      ICNT = 0
      ICP = 0
      DO 35 I=1,MBST
         IBEG = IXB(I)
         IEND = IXB(I+1)-1
         IC = IEND-IBEG+1
         ICNT = ICNT+3+IEND-IBEG
         IF (ICNT.LE.55) GO TO 25
         ICNT = IEND-IBEG+3
         WRITE (NOUT,40)
   25    IF (IBIT.EQ.1 .AND. IC.NE.ICP) WRITE (NOUT,70) IC
         IF (IBIT.EQ.2 .AND. IC.NE.ICP) WRITE (NOUT,75) IC
         IF (IBIT.EQ.3 .AND. IC.NE.ICP) WRITE (NOUT,80) IC
         IF (IC.NE.ICP) ICNT = ICNT+3
         WRITE (NOUT,85)
         DO 30 J=IBEG,IEND
            K = BEST(J,1)
            WRITE (NOUT,90) K, (BEST(J,JJ),JJ=2,4)
   30    CONTINUE
         ICP = IC
   35 CONTINUE
      RETURN
   40 FORMAT (1H1)
   45 FORMAT (//, 17H REGRESSIONS WITH, I4, 12H VARIABLE(S), 8H (R-SQUA,
     *4HRED), /)
   50 FORMAT (//, 17H REGRESSIONS WITH, I4, 12H VARIABLE(S), 8H (ADJUST,
     *13HED R-SQUARED), /)
   55 FORMAT (//, 17H REGRESSIONS WITH, I4, 12H VARIABLE(S), 8H (MALLOW,
     *6HS  CP), /)
   60 FORMAT (8X, 27HCRITERION         VARIABLES)
   65 FORMAT (E20.6, 5X, 25I3, (/25X, 25I3))
   70 FORMAT (//, 22H BEST REGRESSIONS WITH, I4, 12H VARIABLE(S),
     *12H (R-SQUARED))
   75 FORMAT (//, 22H BEST REGRESSIONS WITH, I4, 12H VARIABLE(S),
     *21H (ADJUSTED R-SQUARED))
   80 FORMAT (//, 22H BEST REGRESSIONS WITH, I4, 12H VARIABLE(S),
     *14H (MALLOWS  CP))
   85 FORMAT (1H0, 42H       VARIABLE  COEFFICIENT          PART,
     *21HIAL F           ALPHA)
   90 FORMAT (7X, I3, 2X, 3E20.6)
      END
