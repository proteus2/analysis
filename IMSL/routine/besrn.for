C   IMSL ROUTINE NAME   - BESRN
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - BISERIAL CORRELATION COEFFICIENT FOR A
C                           QUALITATIVELY DICHOTOMIZED VARIABLE AND
C                           A NUMERICALLY OR QUALITATIVELY CLASSIFIED
C                           VARIABLE
C
C   USAGE               - CALL BESRN (N,A,IA,STAT,IER)
C
C   ARGUMENTS    N      - INPUT NUMBER OF CLASSIFICATIONS.
C                A      - INPUT 2 BY N MATRIX.  A(1,*) CONTAINS THE
C                           FREQUENCIES PER CLASSIFICATION ACROSS ONE
C                           OF THE TWO DICHOTOMIZATIONS.  A(2,*)
C                           PROVIDES THE OTHER SET OF FREQUENCIES PER
C                           CLASSIFICATION.
C                IA     - ROW DIMENSION OF MATRIX A EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM. (INPUT)
C                STAT   - OUTPUT VECTOR OF LENGTH 5.  THE I-TH ELEMENT
C                           OF STAT CONTAINS, WHEN
C                           I = 1, THE SUM OF A(1,J), J=1,...,N.
C                           I = 2, THE SUM OF A(2,J), J=1,...,N.
C                           I = 3, THE SUM OF STAT(1) AND STAT(2).
C                           I = 4, THE ABSOLUTE VALUE OF THE RESULTANT
C                             CORRELATION COEFFICIENT ESTIMATE (RN).
C                           I = 5, RN**2.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER=129 INDICATES THAT N IS LESS THAN 1.
C                           IER=130 INDICATES THAT A(1,I) OR A(2,I) IS
C                             LESS THAN OR EQUAL TO ZERO FOR SOME I.
C
C   REQD. IMSL ROUTINES - MDNRIS,MERFI,UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE BESRN  (N,A,IA,STAT,IER)
C
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N,IA,IER
      REAL               A(IA,N),STAT(5)
      DATA               ZERO/0.0/,ONE/1.0/
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
      SUMX = ZERO
      IF(N.LT.1) GO TO 25
    5 DO 10 I = 1,2
         STAT(I) = ZERO
         DO 10 J=1,N
            IF(A(I,J).LE.ZERO) GO TO 30
            STAT(I) = STAT(I) + A(I,J)
   10 CONTINUE
   15 DO 20 J=1,N
         T = A(1,J) + A(2,J)
         P = A(1,J)/T
         CALL MDNRIS(P,X,IER)
         SUMX = SUMX + (T*X*X)
   20 CONTINUE
      STAT(3) = STAT(1) + STAT(2)
      P = STAT(1)/STAT(3)
      CALL MDNRIS(P,X,IER)
      SUMX = SUMX / STAT(3)
      STAT(5) = (SUMX-X*X)/(SUMX+ONE)
      STAT(4) = SQRT(ABS(STAT(5)))
      GO TO 9005
   25 IER = 129
      GO TO 9000
   30 IER = 130
 9000 CONTINUE
      CALL UERTST(IER,' BESRN')
 9005 RETURN
      END
