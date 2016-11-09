C   IMSL ROUTINE NAME   - ZSRCH
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - GENERATE POINTS IN AN N DIMENSIONAL SPACE
C
C   USAGE               - CALL ZSRCH (A,B,N,K,IP,S,M,IW,IER)
C
C   ARGUMENTS    A,B,N  - PARAMETERS WHICH DEFINE THE RECTANGULAR REGION
C                           IN N DIMENSIONAL SPACE. (INPUT)
C                           A AND B ARE VECTORS OF LENGTH N.
C                           GENERATED POINTS SATISFY
C                           A(I) .LT. S(I) .LT. B(I) FOR I=1,2,...,N.
C                           NOTE THAT IF B(I) .LT. A(I), THEN
C                           B(I) .LT. S(I) .LT. A(I).
C                K      - NUMBER OF POINTS TO BE GENERATED. (INPUT)
C                IP     - INITIALIZATION PARAMETER. (INPUT)
C                           IP MUST BE SET TO 0 FOR THE FIRST CALL.
C                           ZSRCH RESETS IP TO 1 AND RETURNS THE
C                           FIRST GENERATED POINT IN S.
C                           SUBSEQUENT CALLS SHOULD BE MADE WITH
C                           IP = 1.
C                S      - VECTOR OF LENGTH N CONTAINING THE
C                           GENERATED POINT. (OUTPUT)
C                           EACH CALL RESULTS IN THE NEXT GENERATED
C                           POINT BEING STORED IN S
C                           (THAT IS S(1),S(2),...,S(N)).
C                M      - WORK VECTOR OF LENGTH N.
C                           M MUST BE PRESERVED BETWEEN CALLS TO ZSRCH.
C                IW     - WORK VECTOR OF LENGTH 9.
C                           IW MUST BE PRESERVED BETWEEN CALLS TO ZSRCH.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER=129, ATTEMPT TO GENERATE MORE THAN
C                             K POINTS.
C
C   REQD. IMSL ROUTINES - UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS      ZSRCH MAY BE USED WITH ANY NONLINEAR OPTIMIZATION
C                ROUTINE THAT REQUIRES STARTING POINTS. THE RECTANGLE
C                TO BE SEARCHED (DEFINED BY PARAMETERS A, B, AND N)
C                MUST BE DETERMINED AND THE NUMBER OF STARTING POINTS,
C                K, MUST BE CHOSEN. ONE POSSIBLE USE FOR ZSRCH WOULD
C                BE TO CALL ZSRCH TO GENERATE A POINT IN THE CHOSEN
C                RECTANGLE. THEN CALL THE NONLINEAR OPTIMIZATION
C                ROUTINE USING THIS POINT AS AN INITIAL GUESS FOR THE
C                SOLUTION. REPEAT THIS PROCESS K TIMES. THE NUMBER OF
C                ITERATIONS THAT THE OPTIMIZATION ROUTINE IS ALLOWED
C                TO PERFORM SHOULD BE QUITE SMALL (5 TO 10) DURING
C                THIS SEARCH PROCESS. THE BEST (OR BEST SEVERAL)
C                POINT(S) FOUND DURING THE SEARCH MAY BE USED AS AN
C                INITIAL GUESS TO ALLOW THE OPTIMIZATION ROUTINE TO
C                DETERMINE THE OPTIMUM MORE ACCURATELY. IN THIS
C                MANNER, AN N DIMENSIONAL RECTANGLE MAY BE
C                EFFECTIVELY SEARCHED FOR A GLOBAL OPTIMUM OF A
C                NONLINEAR FUNCTION. THE CHOICE OF K DEPENDS UPON
C                THE NONLINEARITY OF THE FUNCTION BEING OPTIMIZED.
C                ONE WITH MANY LOCAL OPTIMA REQUIRES A LARGER VALUE
C                THAN ONE WITH ONLY A FEW LOCAL OPTIMA.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE ZSRCH  (A,B,N,K,IP,S,M,IW,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N,K,IP,M(N),IW(9),IER
      REAL               A(N),B(N),S(N)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            L,I,NS,J,JP1,L1,NML1,L2,I1,I2,NJP1,MS,JG,JM,
     1                   NCHALF
      REAL               E
      REAL               C,RJP1,RMI,ZERO,HALF,ONE,RK
      DATA               ZERO,HALF,ONE/0.0,0.5,1.0/
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
C                                  IP .EQ. 0 SIGNALS INITIAL ENTRY TO
C                                    INITIALIZE AND GENERATE FIRST
C                                    POINT
      IF (IP.GT.0) GO TO 15
      E = ONE/N
C                                  COMPUTE J, L1 AND L2 SO THAT
C                                    K = J**(N-L1) * (J+1)**L1 + L2
      RK = K
      J = RK**E
      JP1 = J+1
      IF (JP1**N.LE.K) J = JP1
      JP1 = J+1
      DO 5 L=1,N
         L1 = L
         IF (J**(N-L)*(JP1)**L.GT.K) GO TO 10
    5 CONTINUE
   10 L1 = L1-1
      NML1 = N-L1
      L2 = K-J**NML1*JP1**L1
C                                  I1 IS THE NUMBER OF POINTS THAT CAN
C                                    BE GENERATED WITHOUT ESTABLISHING
C                                    A NEW BASE FOR M
      I1 = J**(NML1-1)*JP1**L1
      IF (L2.EQ.0.AND.L1.GT.0) I1 = (J*I1)/JP1
C
C                                  I2 COUNTS(MODULO I1) THE POINTS AS
C                                    THEY ARE GENERATED
      I2 = 0
      NJP1 = NML1
      IF (L2.EQ.0) NJP1 = NML1+1
      MS = 0
      NS = MIN0(N,NJP1)
C                                  M IS MOVED TO THE NEXT DIAGONAL
C                                    AFTER EACH GROUP OF JG POINTS
      JG = JP1
      IF (NS.EQ.N) JG = J
C                                  PREVENT REJECTION OF POINTS WHEN
C                                    L2 .EQ. 0
      IF (L2.EQ.0) L2 = I1
      IP = 1
      GO TO 20
   15 CONTINUE
C                                  ENTRY TO GENERATE NEXT POINT
C                                    RESTORE LOCAL VARIABLES
      J = IW(1)
      L1 = IW(2)
      L2 = IW(3)
      I1 = IW(4)
      I2 = IW(5)
      NJP1 = IW(6)
      MS = IW(7)
      NS = IW(8)
      JG = IW(9)
      JP1 = J+1
      NML1 = N-L1
   20 CONTINUE
C                                  INCREMENT I2
      I2 = MOD(I2,I1)+1
      IF (I2.GT.1) GO TO 30
C                                  ESTABLISH BASE VALUE FOR M
      DO 25 I=1,N
   25 M(I) = 1
      MS = MS+1
      JM = J
      IF (L1.GT.0.OR.L2.GT.0) JM = JP1
C                                  IF MS .GT. JM ALL K POINTS HAVE BEEN
C                                    GENERATED
      IF (MS.GT.JM) GO TO 9000
      M(NS) = MS
      GO TO 60
C                                  MOVE M ALONG CURRENT DIAGONAL
   30 JM = J
      DO 35 I=1,N
         IF (I.EQ.NJP1) JM = JP1
         M(I) = MOD(M(I),JM)+1
   35 CONTINUE
C                                  CHECK FOR MOVE TO NEW DIAGONAL
      IF (MOD(I2,JG).NE.1) GO TO 60
      IF (NS.EQ.N) GO TO 45
C                                  MOVE M BACK TO LAST DIAGONAL
      DO 40 I=1,NS
         M(I) = M(I)-1
         IF (M(I).EQ.0) M(I) = J
   40 CONTINUE
      M(NS) = MS
C                                  MOVE M TO NEXT DIAGONAL
   45 JM = J
      DO 55 I=1,N
         IF (I.NE.NS) GO TO 50
         JM = JP1
         GO TO 55
   50    M(I) = MOD(M(I),JM)+1
         IF (M(I).GT.1) GO TO 60
   55 CONTINUE
C                                  IF THIS LOOP TERMINATES ALL K
C                                    POINTS HAVE BEEN GENERATED
      GO TO 9000
C                                  REJECT M(NS) .EQ. J+1 WHEN
C                                    I2 .GT. L2
C
   60 IF (I2.GT.L2.AND.M(NS).EQ.JP1) GO TO 20
C
C                                  MAP M INTO S
      C = ZERO
      NCHALF = NML1+1
      IF (I2.LE.L2) NCHALF = NJP1
      RJP1 = ONE/JP1
      DO 65 I=1,N
         IF (I.EQ.NCHALF) C = HALF
         RMI = M(I)
         S(I) = A(I)+(RMI-C)*(B(I)-A(I))*RJP1
   65 CONTINUE
C                                  SAVE LOCAL VARIABLES
      IW(1) = J
      IW(2) = L1
      IW(3) = L2
      IW(4) = I1
      IW(5) = I2
      IW(6) = NJP1
      IW(7) = MS
      IW(8) = NS
      IW(9) = JG
      GO TO 9005
 9000 CONTINUE
      IP = 0
      IER = 129
      CALL UERTST (IER,6HZSRCH )
 9005 RETURN
      END
