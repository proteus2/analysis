C   IMSL ROUTINE NAME   - ZXMWD
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - NOVEMBER 1, 1984
C
C   PURPOSE             - GLOBAL MINIMUM (WITH CONSTRAINTS) OF A
C                           FUNCTION OF N VARIABLES
C
C   USAGE               - CALL ZXMWD(FCN,N,NSIG,A,B,NSRCH,X,F,WORK,
C                           IWORK,IER)
C
C   ARGUMENTS    FCN    - A USER SUPPLIED SUBROUTINE WHICH CALCULATES
C                           F GIVEN X(1),X(2),...,X(N).
C                           FCN IS REFERENCED AS FOLLOWS,
C                             CALL FCN(N,X,F)
C                             WHERE X IS A VECTOR OF LENGTH N
C                             FCN MUST APPEAR IN AN EXTERNAL STATEMENT
C                               IN THE CALLING PROGRAM.
C                             FCN MUST NOT ALTER THE VALUES OF
C                               X(I),I=1,...,N,  OR N.
C                N      - THE NUMBER OF UNKNOWN PARAMETERS. (INPUT)
C                NSIG   - CONVERGENCE CRITERION. (INPUT)  NSIG IS THE
C                           NUMBER OF DIGITS OF ACCURACY REQUIRED IN
C                           THE PARAMETER ESTIMATES.
C                A,B    - CONSTRAINT VECTORS OF LENGTH N. (INPUT)
C                           X(I) IS REQUIRED TO SATISFY -
C                                A(I) .LE. X(I) .LE. B(I)
C                NSRCH  - NUMBER OF STARTING POINTS TO BE GENERATED.
C                           (INPUT) SUGGESTED VALUE = MIN(2**N+5,100)
C                X      - VECTOR OF LENGTH N CONTAINING THE FINAL
C                           PARAMETER ESTIMATES. (OUTPUT)
C                F      - VALUE OF THE FUNCTION AT THE FINAL
C                           PARAMETER ESTIMATES. (OUTPUT)
C                WORK   - REAL WORK VECTOR OF LENGTH N*(N+1)/2+11*N
C                           AN ESTIMATE OF THE NUMBER OF SIGNIFICANT
C                           DIGITS IN THE FINAL PARAMETER ESTIMATES IS
C                           RETURNED IN WORK(1).
C                IWORK  - INTEGER WORK VECTOR OF LENGTH N
C                IER    - ERROR PARAMETER (OUTPUT)
C                         TERMINAL ERROR
C                           IER = 129 INDICATES THAT THE ALGORITHM
C                             HAS CONVERGED TO A POINT WHICH MAY
C                             ONLY BE A SADDLE POINT.
C                           IER = 130 INDICATES THAT IT WAS NOT
C                             POSSIBLE TO CALCULATE THE SOLUTION
C                             TO NSIG DIGITS.  SEE REMARKS.
C                           IER = 131 INDICATES THAT THE ITERATION
C                             WAS TERMINATED AFTER 200*(N+1)
C                             FUNCTION EVALUATIONS.  SEE REMARKS.
C                           IER = 132 INDICATES THAT A(I).GE.B(I)
C                             FOR SOME I=1,...,N. NO ATTEMPT IS MADE
C                             TO FIND THE MINIMUM IN THIS CASE.
C
C   REQD. IMSL ROUTINES - UERSET,UERTST,UGETIO,ZSRCH,ZXMJN,ZXMWE
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS      WHEN IER IS RETURNED AS 130 OR 131, THE PARAMETER
C                ESTIMATES IN X MAY NOT BE RELIABLE.  FURTHER CHECKING
C                SHOULD BE PERFORMED.  USE OF A LARGER NSRCH VALUE MAY
C                PRODUCE MORE RELIABLE PARAMETER ESTIMATES.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE ZXMWD (FCN,N,NSIG,A,B,NSRCH,X,F,WORK,IWORK,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N,NSIG,NSRCH,IER,IWORK(1)
      REAL               FCN,A(N),B(N),X(N),F,WORK(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I1,I2,I5,IOPT,IP,I,J1,J2,JER,JERL,J,KK,K,
     *                   LEVOLD,NG,NH,NSRCH5,NS,NW,NX,IW(9)
      REAL               ASIGSV,BIG,FAST(5),FQ,PI,TEMP
      EXTERNAL           FCN
      DATA               PI /3.141593/
      DATA               BIG /Z7FFFFFFF/
C                                  MAXIMUM NO. OF SEARCH POINTS
      DATA               J2 /5/
C                                  FIRST EXECUTABLE STATEMENT
      CALL UERSET(0,LEVOLD)
      ASIGSV = 0.0
      DO 5 I=1,5
         FAST(I) = BIG
    5 CONTINUE
      NH = N+1
      NG = NH+N*(N+1)/2
      NX = NG+N
      NS = NX+N
      NW = NS+5*N
      DO 10 I=NS,NW
         WORK(I) = 0.0
   10 CONTINUE
      IP = 0
      NSRCH5 = MAX0(NSRCH,5)
      DO 35 I=1,NSRCH5
C                                  GENERATE STARTING POINTS
         CALL ZSRCH(A,B,N,NSRCH5,IP,WORK,IWORK,IW,JER)
         IER = 132
         DO 15 J=1,N
            IF (A(J).GE.B(J)) GO TO 9000
            WORK(J) = PI/2.*(WORK(J)-A(J))/(B(J)-A(J))
   15    CONTINUE
         IER = 0
C                                  DO 4 ITERATIONS WITH EACH
         CALL ZXMWE(FCN,N,NSIG,4*(N+1)+2*N+1,2,WORK,WORK(NH),WORK(NG),
     *   FQ,WORK(NW),A,B,WORK(NX),JER)
         IF (FQ.GE.FAST(5)) GO TO 35
         FAST(5) = FQ
         DO 20 K=1,N
            I5 = NS+4*N+K-1
            WORK(I5) = AMOD(WORK(K),PI)
   20    CONTINUE
         DO 30 J1=1,4
            IF (FAST(J2).GE.FAST(J1)) GO TO 30
            TEMP = FAST(J1)
            FAST(J1) = FAST(J2)
            FAST(J2) = TEMP
            DO 25 K=1,N
               I1 = NS+(J1-1)*N+K-1
               I2 = NS+(J2-1)*N+K-1
               TEMP = WORK(I1)
               WORK(I1) = WORK(I2)
               WORK(I2) = TEMP
   25       CONTINUE
   30    CONTINUE
   35 CONTINUE
      F = BIG
      DO 60 I=1,5
         DO 40 K=1,N
            KK = NS+(I-1)*N+K-1
            WORK(K) = WORK(KK)
   40    CONTINUE
C                                  DO 200 ITERATIONS WITH THE 5
C                                    STARTING POINTS WHICH GAVE
C                                    THE SMALLEST SUM OF SQUARES
C                                    AFTER 4 ITERATIONS
         IOPT = 0
         CALL ZXMWE(FCN,N,NSIG,200*(N+1),IOPT,WORK,WORK(NH),WORK(NG),FQ,
     *   WORK(NW),A,B,WORK(NX),JER)
         IF (FQ.GE.F) GO TO 60
         IOPT = 3
         CALL ZXMWE(FCN,N,NSIG,200*(N+1),IOPT,WORK,WORK(NH),WORK(NG),FQ,
     *   WORK(NW),A,B,WORK(NX),JER)
         IER = JER
         F = FQ
         DO 55 J=1,N
            X(J) = A(J)+(B(J)-A(J))*SIN(WORK(J))**2
   55    CONTINUE
C                                  SAVE NSIG ESTIMATE FOR BEST X
         ASIGSV = WORK(NW+2)
   60 CONTINUE
C                                  RETURN NSIG ESTIMATE IN WORK(1)
 9000 CALL UERSET(LEVOLD,LEVOLD)
      WORK(1) = ASIGSV
      IF (IER.GT.0) CALL UERTST(IER,6HZXMWD )
      RETURN
      END
