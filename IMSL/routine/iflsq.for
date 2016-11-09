C   IMSL ROUTINE NAME   - IFLSQ
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1980
C
C   PURPOSE             - LEAST SQUARES APPROXIMATION WITH USER
C                           SUPPLIED FUNCTIONS
C
C   USAGE               - CALL IFLSQ (F,X,Y,M,A,N,WK,IER)
C
C   ARGUMENTS    F      - NAME OF THE REAL FUNCTION SUBPROGRAM FOR
C                           EVALUATING THE BASIS FUNCTIONS. (INPUT)
C                           THE FUNCTION ITSELF MUST BE SUPPLIED BY THE
C                           USER AND IT SHOULD BE OF THE FOLLOWING FORM
C                             REAL FUNCTION F(K,X)
C                             INTEGER K
C                             REAL X
C                                .
C                                .
C                                .
C                           IT IS ASSUMED THAT THE USER WANTS TO
C                           APPROXIMATE THE DATA BY A SERIES OF THE
C                           FORM
C                             A(1)*F(1,X)+A(2)*F(2,X)+...+A(N)*F(N,X)
C                           F MUST APPEAR IN AN EXTERNAL STATEMENT IN
C                           THE CALLING PROGRAM AND K AND X MUST NOT BE
C                           ALTERED BY F.
C                X      - VECTOR OF LENGTH M CONTAINING THE
C                           ABSCISSAE OF THE DATA POINTS
C                           (X(I),Y(I)),I=1,...,M. (INPUT)
C                Y      - VECTOR OF LENGTH M CONTAINING THE
C                           ORDINATES (OR FUNCTION VALUES) OF
C                           THE DATA POINTS. (INPUT)
C                M      - NUMBER OF DATA POINTS. (INPUT)
C                A      - VECTOR OF LENGTH N CONTAINING THE
C                           COEFFICIENTS OF THE BASIS FUNCTIONS.
C                           (OUTPUT)
C                N      - NUMBER OF BASIS FUNCTIONS. (INPUT)
C                           N MUST BE LESS THAN OR EQUAL TO M.
C                WK     - WORK ARRAY OF LENGTH N*(N+3).
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER = 129, THE PROBLEM DOES NOT HAVE
C                             A UNIQUE SOLUTION. THIS WILL ALWAYS
C                             OCCUR IF N IS GREATER THAN M.
C                           IER = 130, THE PROBLEM IS ILL-CONDITIONED.
C                             THAT IS, THE BASIS FUNCTIONS (RESTRICTED
C                             TO THE DATA SET) ARE ALMOST LINEARLY
C                             DEPENDENT.
C
C   REQD. IMSL ROUTINES - SINGLE/LEQT2P,LUDECP,LUELMP,LUREFP,UERTST,
C                           UGETIO
C                       - DOUBLE/LEQT2P,LUDECP,LUELMP,LUREFP,UERTST,
C                           UGETIO,VXADD,VXMUL,VXSTO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE IFLSQ  (F,X,Y,M,A,N,WK,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            M,N,IER
      REAL               F,X(1),Y(1),A(1),WK(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,ICNT,IDGT,J,K,L,NWK
      REAL               D1,D2
C                                  FIRST EXECUTABLE STATEMENT
      IER = 129
      IF (N .GT. M) GO TO 9000
      IER = 0
      ICNT = 0
      DO 20 I=1,N
         DO 10 J=1,I
            ICNT=ICNT+1
            WK(ICNT)=0.0
            DO 5 K=1,M
               WK(ICNT)=WK(ICNT)+F(I,X(K))*F(J,X(K))
    5       CONTINUE
   10    CONTINUE
         A(I)=0.0
         DO 15 L=1,M
            A(I)=A(I)+F(I,X(L))*Y(L)
   15    CONTINUE
   20 CONTINUE
      NWK=N*(N+1)/2+1
      CALL LEQT2P(WK,1,N,N,A,IDGT,D1,D2,WK(NWK),IER)
 9000 CONTINUE
      IF (IER .NE. 0) CALL UERTST(IER,6HIFLSQ )
 9005 RETURN
      END
