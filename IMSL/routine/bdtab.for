C   IMSL ROUTINE NAME   - BDTAB
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1982
C
C   PURPOSE             - COMPUTATIONS OF FREQUENCIES OF MULTIVARIATE
C                           DATA.
C
C   USAGE               - CALL BDTAB(X,M,KMAX,NOPT,ICNT,K,ITAB,
C                           VECVAL,IVEC,VARVAL,IVAR,IDIST,WK,IER)
C
C   ARGUMENTS    X      - INPUT DATA VECTOR OF LENGTH M
C                M      - INPUT NUMBER OF VARIABLES
C                KMAX   - INPUT MAXIMUM NUMBER OF CELLS.  THIS QUANTITY
C                           DOES NOT HAVE TO BE EXACT, BUT MUST BE AT
C                           LEAST AS LARGE AS THE TRUE VALUE, K.
C                NOPT   - INPUT. OPTION INDICATOR.
C                           IF NOPT IS EQUAL TO ZERO, THE CURRENT X IS
C                             TO BE PROCESSED, BUT NO SORTING IS TO
C                             BE DONE.  (USUALLY, THIS OPTION IS CHOSEN
C                             BECAUSE THERE ARE MORE OBSERVATIONS TO
C                             BE PASSED TO BDTAB.)
C                           IF NOPT IS POSITIVE, THE ELEMENTS
C                             OF VECVAL AND VARVAL ARE TO BE SORTED
C                             AFTER PROCESSING THE CURRENT X.
C                           IF NOPT IS NEGATIVE, THE LAST OBSERVATION IS
C                             ASSUMED TO HAVE BEEN PASSED TO BDTAB ON A
C                             PREVIOUS CALL. X IS IGNORED AND THE
C                             CURRENT CALL ONLY SORTS THE ELEMENTS
C                             IN VECVAL AND VARVAL.
C                ICNT   - INPUT/OUTPUT. ON INITIAL ENTRY ICNT IS SET TO
C                           ZERO IN THE CALLING PROGRAM. ICNT IS
C                           INCREMENTED BY BDTAB EACH TIME THE ROUTINE
C                           IS CALLED AND SHOULD NOT BE CHANGED BETWEEN
C                           CALLS TO BDTAB.
C                K      - OUTPUT NUMBER OF CELLS.  K SHOULD NOT BE
C                           CHANGED BETWEEN CALLS TO BDTAB.
C                ITAB   - OUTPUT VECTOR OF LENGTH K CONTAINING THE
C                           FREQUENCIES FOR THE CELLS.  SINCE IN
C                           GENERAL K WILL NOT BE KNOWN IN ADVANCE,
C                           KMAX (THE UPPER BOUND FOR K) SHOULD BE
C                           USED FOR PURPOSES OF DIMENSIONING ITAB AND
C                           THE FOLLOWING ARRAYS WHOSE DIMENSIONS
C                           INVOLVE K.
C                VECVAL - OUTPUT M BY K ARRAY CONTAINING THE UNIQUE
C                           VALUES OF THE DATA VECTORS.
C                IVEC   - INPUT. ROW DIMENSION OF VECVAL EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM.
C                VARVAL - OUTPUT NMAX BY M ARRAY CONTAINING THE DISTINCT
C                           VALUES OF THE INDIVIDUAL VARIABLES.
C                           NMAX IS THE MAXIMUM ELEMENT IN IDIST (HENCE
C                           IT MAY NOT BE KNOWN IN ADVANCE).
C                IVAR   - INPUT. ROW DIMENSION OF VARVAL EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM.
C                IDIST  - OUTPUT VECTOR OF LENGTH M CONTAINING THE
C                           NUMBER OF DISTINCT VALUES OF EACH VARIABLE.
C                WK     - WORK VECTOR OF LENGTH 2*M.
C                IER    - ERROR PARAMETER (OUTPUT)
C                          TERMINAL ERROR
C                          IER=129 INDICATES NUMBER OF DISTINCT
C                            COMBINATIONS OF VALUES EXCEEDED KMAX.
C                          IER=130 INDICATES NUMBER OF DISTINCT VALUES
C                            OF SOME VARIABLE EXCEEDED IVAR.
C
C   REQD. IMSL ROUTINES - SINGLE/UERTST,UGETIO,VBLA=SCOPY,VSCMP,VSORAP,
C                           VSRTA
C                       - DOUBLE/UERTST,UGETIO,VBLA=DCOPY,VDCMP,VSODAP,
C                           VSRTAD
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE BDTAB (X,M,KMAX,NOPT,ICNT,K,ITAB,VECVAL,IVEC,VARVAL,
     *                   IVAR,IDIST,WK,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            M,KMAX,NOPT,ICNT,K,IVEC,IVAR,IER,ITAB(1),
     *                   IDIST(M)
      REAL               X(M),VECVAL(IVEC,1),VARVAL(IVAR,1),WK(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,ILAST,J,JER
C                                  FIRST EXECUTABLE STATEMENT
      IF (KMAX.GT.0) GO TO 5
      IER = 129
      GO TO 9000
    5 IF (ICNT.GT.0) GO TO 15
      K = 1
      IER = 0
      DO 10 J=1,M
         VECVAL(J,1) = X(J)
         VARVAL(1,J) = X(J)
         IDIST(J) = 1
   10 CONTINUE
      ITAB(1) = 1
      ICNT = 1
      GO TO 9005
   15 IF (NOPT.LT.0) GO TO 60
      ICNT=ICNT+1
      DO 25 I=1,K
         DO 20 J=1,M
            IF (X(J).NE.VECVAL(J,I)) GO TO 25
   20    CONTINUE
         ITAB(I) = ITAB(I)+1
         GO TO 40
   25 CONTINUE
      K = K+1
      IF (K.LE.KMAX) GO TO 30
      IER = 129
      GO TO 9000
   30 DO 35 J=1,M
   35 VECVAL(J,K) = X(J)
      ITAB(K) = 1
   40 DO 55 J=1,M
         ILAST = IDIST(J)
         DO 45 I=1,ILAST
            IF (X(J).EQ.VARVAL(I,J)) GO TO 55
   45    CONTINUE
         IF (ILAST.LT.IVAR) GO TO 50
         IER = 130
         GO TO 9000
   50    IDIST(J) = IDIST(J)+1
         VARVAL(ILAST+1,J) = X(J)
   55 CONTINUE
      IF (NOPT.EQ.0) GO TO 9005
   60 CONTINUE
      CALL VSORAP(VECVAL,IVEC,M,K,M,ITAB,WK,JER)
      DO 65 I=1,M
         CALL VSRTA(VARVAL(1,I),IDIST(I))
   65 CONTINUE
      IF (IER.EQ.0) GO TO 9005
 9000 CONTINUE
      CALL UERTST(IER,'BDTAB ')
 9005 RETURN
      END
