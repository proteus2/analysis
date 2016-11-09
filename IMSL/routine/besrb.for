C   IMSL ROUTINE NAME   - BESRB
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1981
C
C   PURPOSE             - BISERIAL AND POINT-BISERIAL CORRELATION
C                           COEFFICIENTS FOR A QUALITATIVELY
C                           DICHOTOMIZED VARIABLE AND A NUMERICALLY
C                           MEASURABLE AND CLASSIFIED VARIABLE
C
C   USAGE               - CALL BESRB (N,A,IA,STAT,IER)
C
C   ARGUMENTS    N      - INPUT NUMBER OF CLASSIFICATIONS.
C                A      - INPUT 3 BY N MATRIX.  A(1,*) CONTAINS THE
C                           FREQUENCIES PER CLASSIFICATION ACROSS ONE OF
C                           THE TWO DICHOTOMIZATIONS.  A(2,*) PROVIDES
C                           THE OTHER SET OF FREQUENCIES PER CLASSIFICA-
C                           TION.  A(3,*) CONTAINS THE VALUES (CLASS
C                           MARKS) TAKEN ON PER CLASSIFICATION BY THE
C                           MEASURED VARIABLE.
C                IA     - ROW DIMENSION OF MATRIX A EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM. (INPUT)
C                STAT   - OUTPUT VECTOR OF LENGTH 11.  THE I-TH ELEMENT
C                           OF STAT CONTAINS, WHEN
C                           I = 1, THE SUM OF A(1,J), J=1,...,N.
C                           I = 2, THE SUM OF A(2,J), J=1,...,N.
C                           I = 3, THE SUM OF STAT(1) AND STAT(2).
C                           I = 4, THE MEAN OF THE MEASURED VARIABLE.
C                           I = 5, THE MEAN OF THE MEASURED VARIABLE
C                             IN CLASS 1 OF THE DICHOTOMY.
C                           I = 6, THE MEAN OF THE MEASURED VARIABLE
C                             IN CLASS 2 OF THE DICHOTOMY.
C                           I = 7, THE STANDARD DEVIATION OF THE
C                             MEASURED VARIABLE.
C                           I = 8, THE BISERIAL CORRELATION COEFFICIENT
C                             ESTIMATE.
C                           I = 9, THE STANDARD DEVIATION ESTIMATE FOR
C                             THE RESULTANT CORRELATION COEFFICIENT
C                             ESTIMATE.
C                           I = 10, THE ASYMPTOTIC PROBABILITY OF
C                             ERRING WHEN REJECTING THE NULL HYPOTHESIS
C                             OF ZERO CORRELATION.
C                           I = 11, THE POINT BISERIAL CORRELATION
C                             COEFFICIENT ESTIMATE.
C                IER    - ERROR PARAMETERS. (OUTPUT)
C                         TERMINAL ERROR
C                           IER=129 INDICATES THAT N IS LESS THAN 2.
C                           IER=130 INDICATES THAT A(1,I) OR A(2,I) IS
C                             LESS THAN ZERO FOR SOME I.
C                           IER=131 INDICATES THAT STAT(1), STAT(2), OR
C                             STAT(7) IS LESS THAN OR EQUAL TO ZERO.
C
C   REQD. IMSL ROUTINES - MDNOR,MDNRIS,MERFI,MERRC=ERFC,UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE BESRB  (N,A,IA,STAT,IER)
C
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N,IA,IER
      REAL               A(IA,N),STAT(11)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,J
      REAL               ZERO,ONE,TWO,PFIV,SE,SIDSQR,P,PQ,SUM5,SUM7,T,X
      DATA               ZERO/0.0/,ONE/1.0/
      DATA               TWO/2.0/,PFIV/0.5/
      DATA               SE/2.718282/
      DATA               SIDSQR/.3989423/
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
      IF(N.LT.2) GO TO 20
      DO 5 I = 1,2
         STAT(I) = ZERO
         DO 5 J=1,N
            IF(A(I,J).LT.ZERO) GO TO 25
            STAT(I) = STAT(I) + A(I,J)
    5 CONTINUE
      IF (STAT(1).LE.ZERO.OR.STAT(2).LE.ZERO) GO TO 30
      STAT(3) = STAT(1) + STAT(2)
      P = ZERO
      PQ = ZERO
      SUM7 = ZERO
      DO 10   J =1,N
         T = A(1,J) + A(2,J)
         P = P + T * A(3,J)
         PQ = PQ + A(1,J) * A(3,J)
         SUM7 = SUM7 + A(2,J) * A(3,J)
   10 CONTINUE
      STAT(4) = P/STAT(3)
      STAT(5) = PQ/STAT(1)
      STAT(6) = SUM7/STAT(2)
      SUM7 = ZERO
      DO 15 J = 1,N
         T = A(1,J) + A(2,J)
         SUM7 = SUM7 + ((A(3,J)-STAT(4))**2) * T
   15 CONTINUE
      STAT(7) = SQRT(SUM7/STAT(3))
      IF (STAT(7).LE.ZERO) GO TO 30
      P = STAT(1)/STAT(3)
      PQ = P * (ONE-P)
      SUM5 = SQRT(PQ)
      CALL MDNRIS(P,X,IER)
      X = SIDSQR*SE**(-PFIV*X*X)
      STAT(8) = (STAT(5)-STAT(6))*PQ/(STAT(7)*X)
      STAT(9) = SQRT((SUM5/X-STAT(8)**2)/STAT(3))
      P = ABS(STAT(8)) / STAT(9)
      CALL MDNOR(P,SUM7)
      STAT(10) = TWO * (ONE-SUM7)
      STAT(11) = X*STAT(8)/SUM5
      GO TO 9005
   20 IER = 129
      GO TO 9000
   25 IER = 130
      GO TO 9000
   30 IER = 131
 9000 CONTINUE
      CALL UERTST(IER,' BESRB')
 9005 RETURN
      END
