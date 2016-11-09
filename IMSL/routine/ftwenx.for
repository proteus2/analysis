C   IMSL ROUTINE NAME   - FTWENX D
C
C----------------------------------------------------------------------
C
C   LATEST REVISION     - NOVEMBER 1, 1984
C
C   PURPOSE             - MAXIMUM LIKELIHOOD PARAMETER ESTIMATES FOR A
C                           MULTICHANNEL, SINGLE OUTPUT TIME SERIES
C                           MODEL
C
C   USAGE               - CALL FTWENX (X,IX,NS,LS,IP,LAG,ID,R,IR,T,PMAC,
C                           WK,IER)
C
C   ARGUMENTS    X      - INPUT MATRIX OF DIMENSION LS BY NS CONTAINING
C                           NS TIME SERIES, OF WHICH X(*,1) IS THE BASE
C                           TIME SERIES OR OUTPUT CHANNEL.
C                           X IS DESTROYED ON OUTPUT.
C                IX     - INPUT ROW DIMENSION OF THE MATRIX X EXACTLY
C                           AS SPECIFIED IN THE DIMENSION STATEMENT IN
C                           THE CALLING PROGRAM.
C                NS     - INPUT NUMBER OF TIME SERIES. NS MUST BE
C                           GREATER THAN OR EQUAL TO ONE.
C                LS     - INPUT LENGTH OF EACH TIME SERIES.
C                           LS MUST BE GREATER THAN OR EQUAL TO ONE.
C                IP     - INPUT VECTOR OF LENGTH NS. IP(I) CONTAINS THE
C                           NUMBER OF REGRESSIVE PARAMETERS (IN THE
C                           DIFFERENCED FORM OF THE MODEL) ASSOCIATED
C                           WITH TIME SERIES I, I=1,2,...,NS. IP(I)
C                           MUST BE GREATER THAN OR EQUAL TO ZERO FOR
C                           I=1,2,...,NS.
C                LAG    - INPUT VECTOR OF LENGTH NS. LAG(I) CONTAINS THE
C                           LAG OF TIME SERIES I WITH RESPECT TO THE
C                           BASE TIME SERIES. LAG(I) MUST BE GREATER
C                           THAN ZERO FOR I=1,2,...,NS.
C                         NOTE THAT LS MUST BE GREATER THAN OR EQUAL TO
C                           THE MAXIMUM OF (LAG(I)+IP(I))
C                           FOR I=1,2,...,NS.
C                ID     - INPUT VECTOR OF LENGTH NS. IF ID(I) = 0,
C                           THE TIME SERIES I IS TO BE TRANSFORMED BY
C                           REMOVING THE MEAN. IF ID(I) IS POSITIVE,
C                           TIME SERIES I IS DIFFERENCED ID(I) TIMES.
C                           NO TRANSFORMATION IS APPLIED IF ID(I)
C                           IS NEGATIVE. I=1,2,...,NS.
C                R      - WORK AREA MATRIX OF DIMENSION L BY L WHERE L
C                           EQUALS THE SUM OF IP(I), I=1,...,NS.
C                IR     - INPUT ROW DIMENSION OF THE MATRIX R EXACTLY
C                           AS SPECIFIED IN THE DIMENSION STATEMENT IN
C                           IN THE CALLING PROGRAM.
C                T      - OUTPUT VECTOR OF LENGTH N CONTAINING
C                           PARAMETER ESTIMATES OF THE UNDIFFERENCED
C                           MODEL FORM, WHERE N EQUALS THE SUM OF IADD
C                           PLUS THE NS ELEMENTS (IP(I) + MAX(0,ID(I))),
C                           AND IADD IS THE NUMBER OF ELEMENTS ADDED AT
C                           THE BEGINNING BECAUSE OF DIFFERENCING
C                           THE PREDICTED SERIES. (IADD = MIN
C                           (LAG(1)-1,ID(1)). THE PARAMETER
C                           ESTIMATES FOR THE FIRST TIME SERIES ARE
C                           CONTAINED IN T(1),T(2),...,T(IADD+IP(1)+
C                           MAX(0,ID(1))). THE PARAMETER ESTIMATES
C                           FOR THE SECOND TIME SERIES ARE CONTAINED
C                           IN T(I+1),T(I+2),...,T(I+IP(2)+MAX(0,ID(2)),
C                           WHERE I=IADD+IP(1)+MAX(0,ID(1)). THE
C                           PARAMETER ESTIMATES FOR TIME SERIES I ARE
C                           CONTAINED IN T(J+1),T(J+2),...,T(J+IP(I)+
C                           MAX(0,ID(I))), WHERE J=IADD+SUM OVER K=1,
C                           ...,I-1 OF (IP(K)+MAX(0,ID(K)). (I.GT.1).
C                           (SEE REMARK 3.)
C                PMAC   - OUTPUT MOVING AVERAGE CONSTANT IN THE
C                           UNDIFFERENCED MODEL.
C                WK     - WORK AREA VECTOR OF LENGTH MAX(IADD,NS+L),
C                           WHERE L EQUALS THE SUM OF IP(I),I=1,...,NS,
C                           AND IADD IS DEFINED ABOVE WITH VARIABLE T.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER=129 INDICATES AN ERROR OCCURRED IN
C                             IMSL ROUTINE LEQT1F.
C                           IER=130 INDICATES THAT ONE OF IP(I),LAG(I),
C                             IX,NS, OR LS WAS OUT OF RANGE.
C
C   REQD. IMSL ROUTINES - FTCRXY,LEQT1F,LUDATN,LUELMN,UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS  1.  IN PICKING SIGNIFICANT LAGS, CROSS-VARIANCE STUDIES
C                AVAILABLE FROM IMSL ROUTINES FTCROS AND FTFREQ MIGHT
C                BE USED.
C            2.  NOTE THAT THE STRUCTURE OF FTWENX APPARENTLY ALLOWS
C                ONLY PARAMETER ESTIMATION AT LAGS BY EACH OTHER. THIS
C                MAY BE EXTENDED AS IN THE FOLLOWING EXAMPLE;
C                SUPPOSE THAT TWO PARAMETERS ARE DESIRED WITH SERIES
C                Z AT LAGS 5 AND 10. TIME SERIES Z SHOULD BE INCLUDED
C                TWICE IN MATRIX X. FIRST IN COLUMN ONE OF MATRIX X
C                WITH IP(1)=1 AND LAG(1)=5 AND AGAIN IN COLUMN 2
C                WITH IP(2)=1 AND LAG(2)=10.
C            3.  IF BOTH IP(1) AND ID(1) ARE GREATER THAN ZERO, AND
C                LAG(1) IS GREATER THAN ID(1), THEN THE SERIES ONE
C                PARAMETERS MAY HAVE A JUMP IN LAG INDICES.  T(1)
C                WILL CORRESPOND TO SERIES 1 AT LAG 1, T(2) WILL
C                CORRESPOND TO SERIES ONE AT LAG 2, ETC., UP TO
C                T(ID(1)).  T(ID(1)+1) WILL CORRESPOND TO SERIES 1
C                AT LAG=LAG(1), T(ID(1)+2) WILL CORRESPOND TO SERIES 1
C                AT LAG=LAG(1)+1, ETC..
C
C-----------------------------------------------------------------------
C
      SUBROUTINE FTWENX (X,IX,NS,LS,IP,LAG,ID,R,IR,T,PMAC,WK,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            IX,NS,LS,IP(NS),LAG(NS),ID(NS),IR,IER
      REAL               X(IX,1),R(IR,1),T(1),PMAC,WK(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,ICHECK,ICOL,IADD,ID1,IDD,IDIF,IDMAX,IID,IKM,
     1                   ILINE,IPP,IRS,ISAVE,ISUM,ITEMP,IUP,I0,I1,J,K,L,
     2                   LMI,LM1,M,MCOL,MLAG,MLINE,N
      REAL               XBAR,SAVE,Z0,SUM
      DOUBLE PRECISION   TEMP
      DATA               I0/0/,I1/1/,Z0/0.0/
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
      ISUM = 0
      ICHECK = 0
      IDMAX = 0
      IRS = 0
      DO 5 I = 1,NS
         ISUM = ISUM + IP(I)
         ITEMP = LAG(I) + IP(I)
         ICHECK = MAX0(ICHECK,ITEMP)
         IDMAX = MAX0(IDMAX,ID(I))
         IRS = MAX0(IRS,ITEMP)
         IF (IP(I).LT.0 .OR. LAG(I).LT.1) IER = 130
    5 CONTINUE
      IF (IX.LT.LS .OR. NS.LT.1 .OR. LS.LT.ICHECK) IER = 130
      IF (IER .NE. 0) GO TO 9000
C                                  TRANSFORM SERIES
      DO 35 M = 1,NS
         WK(M) = 0.0
         IF (ID(M) .LT. 0) GO TO 35
         N = LS
         IF (ID(M) .EQ. 0) GO TO 20
         IID = ID(M)
         DO 15 K = 1,IID
            N = N - 1
            DO 10 I = 1,N
               X(I,M) = X(I+1,M) - X(I,M)
   10       CONTINUE
   15    CONTINUE
         GO TO 35
   20    TEMP = 0.0D0
         DO 25 I = 1,N
            TEMP = TEMP + DBLE(X(I,M))
   25    CONTINUE
         XBAR = TEMP / N
         WK(M) = XBAR
         DO 30 I = 1,N
            X(I,M) = X(I,M) - XBAR
   30    CONTINUE
   35 CONTINUE
C                                  SUM BOUNDS
      N = LS - IDMAX
C                                  CALCULATE MATRIX R, VECTOR T
C                                  OF MULTI-CHANNEL COVARIANCES
      MLINE = 0
      DO 75 M = 1,NS
         IF (IP(M) .EQ. 0) GO TO 75
         MCOL = 0
         DO 65 K = 1,M
            IKM = LAG(K) - LAG(M)
            IF (IP(K) .EQ. 0) GO TO 65
            IUP = IP(M)
            LM1 = 0
            DO 45 L = 1,IUP
               MLAG = IKM - LM1
               CALL FTCRXY (X(1,K),X(1,M),N,Z0,Z0,MLAG,IRS,SAVE,IER)
               ITEMP = LM1 + MLINE
               J = MIN0(IP(K),IP(M)-LM1)
               DO 40 I = 1,J
                  ILINE = I + ITEMP
                  ICOL = I + MCOL
                  R(ILINE,ICOL) = SAVE
                  R(ICOL,ILINE) = SAVE
   40          CONTINUE
            LM1 = L
   45       CONTINUE
            IF (K .EQ. M) GO TO 60
            IF (IP(K) .EQ. 1) GO TO 60
            IUP = IP(K)
            LM1 = 1
            DO 55 L = 2,IUP
               MLAG = IKM + LM1
               CALL FTCRXY (X(1,K),X(1,M),N,Z0,Z0,MLAG,IRS,SAVE,IER)
               ITEMP = LM1 + MCOL
               J = MIN0(IP(M),IP(K)-LM1)
               DO 50 I = 1,J
                  ILINE = I + MLINE
                  ICOL = I + ITEMP
                  R(ILINE,ICOL) = SAVE
                  R(ICOL,ILINE) = SAVE
   50          CONTINUE
            LM1 = L
   55       CONTINUE
   60       MCOL = MCOL + IP(K)
   65    CONTINUE
         ITEMP = 1 - LAG(M)
         IUP = IP(M)
         DO 70 K = 1,IUP
            MLAG = ITEMP - K
            CALL FTCRXY (X(1,1),X(1,M),N,Z0,Z0,MLAG,IRS,SAVE,IER)
            T(MLINE+K) = SAVE
   70    CONTINUE
         MLINE = MLINE + IP(M)
   75 CONTINUE
C                            SOLVE LINEAR EQUATIONS
      CALL LEQT1F (R,I1,ISUM,IR,T,I0,WK(NS+1),IER)
      IF (IER .EQ. 129) GO TO 9000
C                                  CALCULATE PMAC
      ISUM = 0
      PMAC = 0.0
      IADD = MIN0(LAG(1)-1,ID(1))
      IF (IP(1).EQ.0)  IADD = ID(1)
      IADD = MAX0(IADD,0)
      IDIF = IADD
      DO 90 I = 1,NS
         IDIF = IDIF + MAX0(0,ID(I))
         IF (IP(I).EQ.0 .OR. ID(I).NE.0) GO TO 85
         SUM = 0.0
         K = IP(I)
         DO 80 J = 1,K
            SUM = SUM + T(J+ISUM)
   80    CONTINUE
         IF (I .EQ. 1) PMAC = WK(1)*(1.0 - SUM)
         IF (I .GT. 1) PMAC = PMAC - WK(I)*SUM
   85    ISUM = ISUM + IP(I)
   90 CONTINUE
C                                  UNDIFFERENCE PARAMETERS
      IF (IDIF .EQ. 0) GO TO 9005
      L = ISUM + 1
      K = L + IDIF
C                                  MAKE ROOM FOR UNDIFFERENCED PARMS.
      DO 95 I = 1,ISUM
         T(K-I) = T(L-I)
   95 CONTINUE
      ISUM = IADD
      DO 101 I = 1, ISUM
         T(I) = 0.0
  101 CONTINUE
C                                  COMPUTE UNDIFFERENCED PARAMETERS FOR
C                                  EACH TIME SERIES.
      DO 125 M = 1,NS
         IPP = IP(M)
         IF (IPP .EQ. 0) GO TO 125
         DO 100 I = 1,IPP
            T(I+ISUM) = T(I+IDIF)
  100    CONTINUE
         IDD = ID(M)
         IF (IDD .LE. 0) GO TO 120
         ITEMP = IPP + ISUM - 1
         DO 115 K = 1,IDD
            L = ITEMP + K
            T(L+1) = -T(L)
            LMI = L - ISUM
            IF (LMI .LE. 1) GO TO 115
            ISAVE = L + 2
            DO 105 J = 2,LMI
               I = ISAVE - J
               T(I) = T(I) - T(I-1)
  105       CONTINUE
  115    CONTINUE
  120    ISUM = ISUM + IPP + MAX0(0,IDD)
         IDIF = IDIF + IPP
  125 CONTINUE
C                                  ADJUST FOR THE UNDIFFERENCED
C                                  PREDICTED SERIES.
      IF (ID(1).EQ.0)  GO TO 9005
      WK(1) = 1.0
      ID1 = ID(1)
      DO 130 I = 1, ID1
        WK(I+1) = 0.0
  130 CONTINUE
      DO 140 I = 1, ID1
         DO 145 JM1 = 1, I
            J = I+2-JM1
            WK(J) = WK(J) - WK(J-1)
  145    CONTINUE
  140 CONTINUE
      DO 150 I = 1, ID1
         T(I) = T(I) + WK(I+1)
  150 CONTINUE
      GO TO 9005
 9000 CONTINUE
      CALL UERTST (IER,'FTWENX')
 9005 RETURN
      END
