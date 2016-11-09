C   IMSL ROUTINE NAME   - GTTRT
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - TALLY FOR TRIPLETS TEST
C
C   USAGE               - CALL GTTRT (R,N,K,A,IA1,IA2,IER)
C
C   ARGUMENTS    R      - INPUT VECTOR OF LENGTH N CONTAINING THE
C                           SEQUENCE OF FLOATING POINT PSEUD0-UNIFORM
C                           (0,1) RANDOM NUMBERS
C                N      - INPUT LENGTH OF VECTOR R AND MUST BE EVENLY
C                           DIVISIBLE BY 3.
C                K      - INPUT NUMBER OF SUBDIVISIONS OF THE UNIT LINE.
C                           K MAY BE THE CUBE ROOT OF THE K OBTAINED
C                           FROM IMSL ROUTINE GTCN. K SHOULD BE GREATER
C                           THAN OR EQUAL TO 2 AND LESS THAN OR EQUAL
C                           TO 58
C                A      - OUTPUT THREE DIMENSIONAL TALLY MATRIX. A IS
C                           DIMENSIONED K BY K BY K. A MAY CONTAIN
C                           PARTIAL SUMS ON ENTRY. ON THE FIRST CALL
C                           TO GTTRT, THE ENTRIES OF A MUST BE SET TO
C                           ZERO.
C                IA1    - INPUT FIRST DIMENSION OF MATRIX A EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM.
C                IA2    - INPUT SECOND DIMENSION OF MATRIX A EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER = 129 INDICATES THAT N IS NOT DIVISIBLE
C                             BY 3.
C
C   REQD. IMSL ROUTINES - UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE GTTRT  (R,N,K,A,IA1,IA2,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N,K,IA1,IA2,IER
      REAL               A(IA1,IA2,1),R(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,J,L,M
      REAL               C
C                                  FIRST EXECUTABLE STATEMENT
      IF ((N/3)*3 .EQ. N) GO TO 5
C                                  TERMINAL ERROR
      IER=129
      GO TO 9000
    5 IER=0
      C=K
C                                  UPDATE TALLY MATRIX
         DO 10 I=1,N,3
         J=R(I)*C+1
         L=R(I+1)*C+1
         M=R(I+2)*C+1
   10    A(J,L,M)=A(J,L,M)+1
      GO TO 9005
 9000 CONTINUE
      CALL UERTST(IER,'GTTRT')
 9005 RETURN
      END
