C   IMSL ROUTINE NAME   - NMKN
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1982
C
C   PURPOSE             - KENDALLS TEST FOR CORRELATION
C                          (RANK CORRELATION COEFFICIENT)
C
C   USAGE               - CALL NMKN (X,Y,N,STAT,XSTAT,IWK,WK,IER)
C
C   ARGUMENTS    X      - INPUT VECTOR OF LENGTH N CONTAINING
C                           SAMPLE 1.
C                           ON OUTPUT, X CONTAINS THE SORTED RANKS
C                           OF INPUT X.
C                Y      - INPUT VECTOR OF LENGTH N CONTAINING
C                           SAMPLE 2.
C                           ON OUTPUT, Y CONTAINS THE SORTED RANKS
C                           OF INPUT Y.
C                N      - INPUT SIZE OF SAMPLES 1 AND 2. N MUST BE
C                           GREATER THAN OR EQUAL TO 3.
C                STAT   - INPUT/OUTPUT VECTOR OF LENGTH 7.
C                         ON INPUT, STAT(1) MUST CONTAIN THE VALUE OF
C                           EPSILON, WHICH IS USED TO DETECT TIES IN X
C                           OR Y. IF INPUT ELEMENTS DIFFER BY NO MORE
C                           THAN EPSILON, A TIE IS COUNTED.
C                         ON OUTPUT
C                         STAT(1) CONTAINS THE KENDALL TAU(A)
C                         STAT(2) CONTAINS THE KENDALL TAU(B)
C                         STAT(3) CONTAINS THE TIES STATISTIC FROM
C                           SAMPLE X (0.0 IF NO TIES DETECTED)
C                         STAT(4) CONTAINS THE TIES STATISTIC FROM
C                           SAMPLE Y (0.0 IF NO TIES DETECTED)
C                         STAT(5) CONTAINS THE EXACT PROBABILITY OF
C                           ACHIEVING A SCORE AT LEAST AS LARGE AS S
C                           (CALCULATED WHEN N IS NOT TOO LARGE FOR IMSL
C                           ROUTINE NMKSF AND THERE ARE NO TIES. IN
C                           EITHER OF THOSE CASES, IT IS SET TO -1.0.)
C                         STAT(6) CONTAINS THE SAME PROBABILITY,
C                           ASSUMING NORMALITY (CALCULATED WHEN N IS
C                           GREATER THAN 7, SET TO -1.0 OTHERWISE)
C                         STAT(7) CONTAINS THE SAME PROBABILITY,
C                           ASSUMING NORMALITY, WITH CONTINUITY
C                           CORRECTION (CALCULATED WHEN N IS GREATER
C                           THAN 7, SET TO -1.0 OTHERWISE)
C                XSTAT  - OUTPUT VECTOR OF LENGTH N*(N-1)/2+3
C                         XSTAT(1) CONTAINS S, THE KENDALL SCORE FROM
C                           THE TEST. S IS IN THE RANGE
C                           (-N*(N-1)/2,N*(N-1)/2), INCLUSIVELY.
C                         XSTAT(2) + 2 IS THE INDEX I WHICH
C                           CORRESPONDS TO THE ELEMENT XSTAT(I)
C                           THAT CONTAINS THE FREQUENCY CORRESPONDING
C                           TO THE S OBSERVED, FOR I=3,...,N*(N-1)/2+3.
C                           SEE PROGRAMMING NOTES IN THE MANUAL
C                           DOCUMENT FOR FURTHER DETAILS.
C                         XSTAT(3),XSTAT(4),...,XSTAT(N*(N-1)/2+3)
C                           CONTAIN THE POSSIBLE FREQUENCIES OF
C                           OCCURRENCE OF S, IN AN N-SIZED SAMPLE
C                         XSTAT(2) THROUGH XSTAT(N*(N-1)/2+3) ONLY
C                           CONTAIN VALID INFORMATION WHEN EXACT
C                           PROBABILITIES ARE CALCULATED (WHEN WARNING
C                           ERROR IER=34 DOES NOT OCCUR. SEE DESCRIPTION
C                           OF IER BELOW.)
C                IWK    - WORK VECTOR OF LENGTH N.
C                WK     - WORK VECTOR OF LENGTH (N-1)*(N-2)/2+1.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER=129 INDICATES SAMPLE SIZE WAS LESS
C                             THAN 3.
C                         WARNING ERROR
C                           IER=34 INDICATES THAT TIES EXIST OR THAT
C                             N WAS TOO LARGE TO PERFORM ALL OF THE
C                             COMPUTATIONS IN NMKSF. STAT(5) IS SET
C                             TO -1.0.  SEE ALGORITHM DESCRIPTION AND
C                             PROGRAMMING NOTES.
C                           IER=35 INDICATES N IS LESS THAN OR EQUAL
C                             TO 7 (STAT(6)=STAT(7)=-1.0)
C
C   REQD. IMSL ROUTINES - MERRC=ERFC,NMKSF,NMTIE,NNUC,UERTST,
C                           UGETIO,VSRTR
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE NMKN  (X,Y,N,STAT,XSTAT,IWK,WK,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N,IWK(1),IER
      REAL               X(1),Y(1),STAT(1),XSTAT(1),WK(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,IP,IS,K,M,NMAX
      REAL               PZF,S,T,XINF,XN,XNM,XXX
      DATA               PZF/.5555556E-1/
      DATA               XINF/Z7FFFFFFF/
C                                  FIRST EXECUTABLE STATEMENT
      IF (N .GE. 3) GO TO 5
      IER = 129
      GO TO 9000
    5 IER = 0
      NMAX = 33
      XXX = XINF/1.E11
      IF (XXX.GT.1.E38) NMAX = 55
C                                  ACCUMULATE K=S BY COMPARING DIRECTION
C                                  OF DIFFERENCES WITHIN CORRESPONDING
C                                  X AND Y PAIRS
      CALL NNUC (X,N,STAT(1),IWK,WK,X,S,S)
      CALL NNUC (Y,N,STAT(1),IWK,WK,Y,S,S)
      K = 0
      IS = N - 1
      DO 35 I=1,IS
         IP = I + 1
         DO 30 M=IP,N
            IF (X(I)-X(M)) 10,30,15
   10       IF (Y(I)-Y(M)) 20,30,25
   15       IF (Y(I)-Y(M)) 25,30,20
   20       K = K + 1
            GO TO 30
   25       K = K - 1
   30    CONTINUE
   35 CONTINUE
      S = K
      XN = N
      XNM = XN*(XN-1.0)
C                                  SORT X AND Y VECTORS
      CALL VSRTR (X,N,IWK)
      CALL VSRTR (Y,N,IWK)
C                                  DETERMINE TIE STATISTICS FOR X RANKS
      CALL NMTIE(X,N,STAT(1),STAT(3))
      T = 0.5 * XNM
      CALL NMTIE(Y,N,STAT(1),XSTAT)
      STAT(4) = XSTAT(1)
C                                  COMPUTE TAU(A)
      STAT(1) = (S+S)/XNM
      STAT(7) = PZF*(XNM*(XN+XN+5.0)-STAT(5)-XSTAT(3))+XSTAT(4)*
     1   STAT(6)/(9.0*XNM*(XN-2.0))+XSTAT(1)*STAT(3)/T
C                                  COMPUTE TAU(B)
      STAT(2) = S/SQRT((T-STAT(3))*(T-STAT(4)))
      IF(N.LE.NMAX .AND. STAT(3).EQ.0. .AND. STAT(4).EQ.0.) GO TO 40
      IER = 34
      STAT(5) = -1.0
      GO TO 45
C                                  INVOKE NMKSF FOR PROBABILITY
C                                  OF S OR GREATER AND FREQUENCIES
   40 CALL NMKSF(K,N,XSTAT(3),WK,STAT(5))
      XSTAT(2) = WK(1)
   45 IF (N .GT. 7) GO TO 50
      IER = 35
      STAT(6) = -1.0
      STAT(7) = -1.0
      GO TO 55
   50 STAT(7) = SQRT(STAT(7))
      T = S/STAT(7)
C                                  PROB(S OR GREATER) ASSUMING NORMALITY
      STAT(6) = 0.5 * ERFC(T*.7071068)
C                                  SAME PROB W/CONTINUITY CORRECTION
      M = IABS(K) - 1
      T = ISIGN(M,K)
      T = T/STAT(7)
      STAT(7) = 0.5 * ERFC(T*.7071068)
C                                  PUT S AND ITS INDEX INTO XSTAT
   55 XSTAT(1) = S
      IF (IER .EQ. 0) GO TO 9005
 9000 CONTINUE
      CALL UERTST (IER,6HNMKN  )
 9005 RETURN
      END
