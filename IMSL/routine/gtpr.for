C   IMSL ROUTINE NAME   - GTPR
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - TALLY OF COORDINATES OF PAIRS (OR LAGGED
C                           PAIRS) OF RANDOM NUMBERS
C
C   USAGE               - CALL GTPR (R,N,K,L,A,IA,IER)
C
C   ARGUMENTS    R      - INPUT VECTOR OF LENGTH N CONTAINING A
C                           SEQUENCE OF FLOATING POINT PSEUDO-
C                           UNIFORM (0,1) RANDOM NUMBERS. THE PAIRS,
C                           (R(I),R(I+L)), WILL BE CONSIDERED BY GTPR.
C                N      - INPUT LENGTH OF VECTOR R. IF A PAIRS TALLY IS
C                           BEING USED (L = 1), THEN N MUST BE EVEN.
C                K      - INPUT NUMBER OF SUBDIVISIONS OF THE UNIT LINE.
C                           K SHOULD BE GREATER THAN OR EQUAL TO 2 AND
C                           LESS THAN OR EQUAL TO 448. K MAY BE OBTAINED
C                           FROM IMSL ROUTINE GTCN AS THE SQUARE ROOT OF
C                           OUTPUT SCALAR K.
C                L      - INPUT SCALAR VALUE. COORDINATES (R(I),R(I+L))
C                           WILL BE CONSIDERED FOR I=1,2,3,...,N-L IN
C                           THE TALLY. THEREFORE, L = 1 WOULD IMPLY A
C                           PAIRS TALLY AND L GREATER THAN 1 WOULD
C                           IMPLY A LAGGED PAIRS TALLY.
C                A      - INPUT/OUTPUT K BY K TALLY MATRIX. ON ANY CALL
C                           TO GTPR, A MAY CONTAIN PARTIAL TALLIES IF
C                           THE USER SO DESIRES. IF, ON THE FIRST
C                           CALL TO GTPR, A DOES NOT CONTAIN PARTIAL
C                           TALLIES, THEN THE ELEMENTS OF A SHOULD BE
C                           SET TO ZERO.
C                IA     - INPUT ROW DIMENSION OF MATRIX A EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN
C                           THE CALLING PROGRAM.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR.
C                           IER = 129 INDICATES THAT N IS ODD IF L IS
C                             EQUAL TO ONE OR THAT N IS LESS THAN OR
C                             EQUAL TO L IF L IS GREATER THAN ONE.
C
C   REQD. IMSL ROUTINES - UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE GTPR  (R,N,K,L,A,IA,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N,K,IA,IER,L
      REAL               R(N),A(IA,1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,J,M,NEND,NINC
      REAL               C
C                                  FIRST EXECUTABLE STATEMENT
      IF (L.EQ.1) GO TO 5
C                                  LAGGED PAIRS TALLY
      IF (L.GE.N) GO TO 20
      NINC = 1
      NEND = N-L
      GO TO 10
C                                  PAIRS TALLY
    5 IF ((N/2)*2.NE.N) GO TO 20
      NINC = 2
      NEND = N
   10 IER = 0
      C = K
C                                  UPDATE THE TALLY MATRIX A
      DO 15 I=1,NEND,NINC
         J = R(I)*C+1
         M = I+L
         M = R(M)*C+1
   15 A(J,M) = A(J,M)+1
      GO TO 9005
   20 IER = 129
      CALL UERTST (IER,'GTPR  ')
 9005 RETURN
      END
