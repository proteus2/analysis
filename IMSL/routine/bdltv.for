C   IMSL ROUTINE NAME   - BDLTV
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1980
C
C   PURPOSE             - PRODUCE LETTER VALUE SUMMARY
C
C   USAGE               - CALL BDLTV (X,N,NUM,SUMRY,IER)
C
C   ARGUMENTS    X      - INPUT VECTOR OF LENGTH N CONTAINING THE DATA
C                           SORTED IN ASCENDING VALUE.
C                N      - INPUT NUMBER OF OBSERVATIONS.
C                NUM    - INPUT NUMBER OF SUMMARY NUMBERS DESIRED. NUM
C                           MUST BE AN ODD INTEGER GREATER THAN OR
C                           EQUAL TO 3. A COMMON VALUE FOR NUM IS 5.
C                SUMRY  - OUTPUT VECTOR OF LENGTH NUM CONTAINING THE
C                           LETTER VALUES. IF NUM IS 5, FOR EXAMPLE,
C                           SUMRY WILL CONTAIN THE MINIMUM, THE LOWER
C                           HINGE (QUARTILE), THE MEDIAN, THE UPPER
C                           HINGE (QUARTILE), AND THE MAXIMUM.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER = 129 INDICATES THAT NUM IS NOT AN
C                             ODD INTEGER GREATER THAN OR EQUAL TO 3.
C
C   REQD. IMSL ROUTINES - UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE BDLTV  (X,N,NUM,SUMRY,IER)
C
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N,NUM,IER
      REAL               X(N),SUMRY(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            ITEMP1,ITEMP2,L,I1,I2,L1,L2,N1,N2
      REAL               XI1,HALF
      DATA               HALF/0.5/
C                                  FIRST EXECUTABLE STATEMENT
      IER=129
      ITEMP1=NUM/2
      ITEMP2=ITEMP1*2
      IF (NUM.EQ.ITEMP2) GO TO 9000
      IF (NUM.LT.3) GO TO 9000
      IER=0
      L=NUM/2 + 1
      XI1 = .9+N/2.
      I1 = XI1
      I2 = N/2+1
      SUMRY(L) = HALF*(X(I1)+X(I2))
      IF (L.LE.2) GO TO 10
      DO 5 J=3,L
         XI1 = .9+I1/2.
         I2 = I1/2+1
         I1 = XI1
         L1=L+2-J
         L2=L-2+J
         N1=N-I1+1
         N2=N-I2+1
         SUMRY(L1) = HALF*(X(I1)+X(I2))
         SUMRY(L2) = HALF*(X(N1)+X(N2))
    5 CONTINUE
   10 SUMRY(1) = X(1)
      SUMRY(NUM) = X(N)
      IF (IER.EQ.0) GO TO 9005
 9000 CONTINUE
      CALL UERTST(IER,'BDLTV ')
 9005 RETURN
      END
