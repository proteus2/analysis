C   IMSL ROUTINE NAME   - OPRINC
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - PRINCIPAL COMPONENTS OF A MULTIVARIATE SAMPLE
C                           OF OBSERVATIONS
C
C   USAGE               - CALL OPRINC (S,M,IA,EVAL,EVEC,COMP,VAR,CL,CU,
C                           IER)
C
C   ARGUMENTS    S      - INPUT VARIANCE-COVARIANCE (OR CORRELATION)
C                           MATRIX OF ORDER M.  S IS IN SYMMETRIC
C                           STORAGE MODE (A VECTOR OF LENGTH M(M+1)/2).
C                           S IS DESTROYED ON OUTPUT.
C                M      - INPUT ORDER OF MATRIX S.
C                IA     - INPUT ROW DIMENSION OF MATRICES EVEC AND
C                           COMP EXACTLY AS SPECIFIED IN THE DIMENSION
C                           STATEMENT IN THE CALLING PROGRAM.
C                EVAL   - OUTPUT VECTOR OF LENGTH M, OF EIGENVALUES
C                           OF S, ORDERED MONOTONIC INCREASING.
C                           NEAR ZERO VALUES MADE NEGATIVE BY ROUND
C                           OFF ERRORS ARE SET TO ZERO.
C                EVEC   - OUTPUT MATRIX OF ORTHONORMAL EIGENVECTORS
C                           OF S.  THE ITH COLUMN OF EVEC CORRESPONDS
C                           TO THE ITH EIGENVALUE IN EVAL.  EVEC IS
C                           M BY M.
C                COMP   - OUTPUT M BY M COMPONENT CORRELATION MATRIX.
C                           (CAN BE IN EVEC LOCATIONS, SEE REMARKS).
C                VAR    - OUTPUT VECTOR, LENGTH M, OF PERCENTAGES OF
C                           TOTAL VARIANCE ASSOCIATED WITH THE
C                           COMPONENTS, IN THE ORDER OF THE EIGENVALUES
C                           EVAL.
C                CL     - ON INPUT, CL(1) CONTAINS THE NUMBER OF
C                           OBSERVATIONS N ON WHICH S WAS BASED.
C                           ON OUTPUT, CL CONTAINS THE LEFT 95 PERCENT
C                           CONFIDENCE BOUNDS ON THE EIGENVALUES
C                           EVAL(I),I=1,2,...,M.
C                CU     - RIGHT 95 PERCENT CONFIDENCE BOUNDS ON THE
C                           EIGENVALUES EVAL(I),I=1,2,...,M.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           N = 129 MEANS AN ERROR WAS RETURNED FROM
C                             SUBROUTINE EQRT2S
C                         WARNING ERROR
C                           IER = 34 INDICATES THAT CU(M) CANNOT BE
C                             CALCULATED. ON INPUT, CL(1) WAS LESS THAN
C                             9. CALCULATION CONTINUES AS IF CL(1) WERE
C                             9. SEE DOCUMENTATION FOR FURTHER DETAILS.
C
C   REQD. IMSL ROUTINES - EHOBKS,EHOUSS,EQRT2S,UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS  1.  PRIOR TO ENTRY, INPUT S MAY BE COMPUTED USING THE
C                BASIC STATISTICS CHAPTER SUBROUTINES. IF S IS A
C                CORRELATION MATRIX, THE CONFIDENCE LIMITS ARE NOT
C                VALID. SAMPLING PROPERTIES AS USED HERE ARE BASED
C                ON THE ASSUMPTION OF NORMALITY AND ON THE VARIANCE-
C                COVARIANCE MATRIX.
C            2.  IF EIGENVECTORS ARE NOT REQUIRED, COMP CAN BE IN THE
C                SAME LOCATIONS AS GIVEN TO EVEC. IN THIS CASE, HOWEVER,
C                ESTIMATES OF COMMUNALITIES (SUMS OF SQUARES OF THE ROWS
C                OF EVEC) WILL NOT BE AVAILABLE TO THE USER ON EXIT
C                FROM OPRINC.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE OPRINC  (S,M,IA,EVAL,EVEC,COMP,VAR,CL,CU,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            M,IA,IER
      REAL               S(1),EVAL(1),VAR(1),CL(1),CU(1),EVEC(IA,1),
     1                   COMP(IA,1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,J,K
      REAL               ZERO,ONE,RNINE,HUND,SCALE,AN,SUM,SUMR,ANP,ANM
      DATA               ZERO/0./,ONE/1./,RNINE/9./,HUND/100./
      DATA               SCALE/2.771859/
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
      AN=CL(1)
C                                  INITIALIZE EVEC AS THE N X N
C                                  IDENTITY MATRIX FOR THE SUBROUTINE
C                                  EQRT2S
      DO 10 I = 1,M
         DO 5 J=1,M
            EVEC(I,J) = ZERO
    5    CONTINUE
         EVEC(I,I) = ONE
   10 CONTINUE
      K = 0
      DO 15 I=1,M
         K=K+I
         CL(I)= SQRT(S(K))
   15 CONTINUE
C                                  COMPUTE EIGENVALUES AND EIGEN-
C                                  VECTORS OF S
      CALL EHOUSS (S,M,EVAL,VAR,CU)
      CALL EQRT2S (EVAL,VAR,M,EVEC,IA,IER)
      IF (IER.NE.0) GO TO 9000
      CALL EHOBKS (S,M,1,M,EVEC,IA)
C                                  COMPUTE THE VAR VECTOR
      SUM = ZERO
      DO  20  I=1,M
         IF (EVAL(I) .LT. ZERO) EVAL(I) = ZERO
         SUM = SUM + EVAL(I)
   20 CONTINUE
      SUMR = HUND/SUM
      DO 25  I=1,M
         VAR(I) = EVAL(I)*SUMR
   25 CONTINUE
C                                  COMPUTE COMP (CORRELATION) MATRIX
      DO 30  I=1,M
         SUMR = ONE/CL(I)
         DO 30  J=1,M
         COMP(I,J) = SQRT(EVAL(J))*EVEC(I,J)*SUMR
   30 CONTINUE
      IF (AN .GE. RNINE) GO TO 35
      AN = RNINE
      IER = 34
C                                  COMPUTE VECTOR CL
   35 AN = SQRT(AN-ONE)
C                                  SCALE=1.96*SQRT(2.)
      ANP = AN/(AN+SCALE)
      ANM = AN/(AN-SCALE)
      DO 40  I=1,M
         CL(I) = EVAL(I)*ANP
         CU(I) = EVAL(I)*ANM
   40 CONTINUE
      IF(IER.EQ.0) GO TO 9005
 9000 CONTINUE
      CALL UERTST(IER,6HOPRINC)
 9005 RETURN
      END
