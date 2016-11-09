C   IMSL ROUTINE NAME   - FTCRXY
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - CROSS-COVARIANCE OF TWO MUTUALLY
C                           STATIONARY TIME SERIES
C
C   USAGE               - CALL FTCRXY (X,Y,N,XBAR,YBAR,MLAG,IRS,C,IER)
C
C   ARGUMENTS    X      - INPUT VECTOR OF LENGTH N CONTAINING
C                           THE BASE TIME SERIES.
C                Y      - INPUT VECTOR OF LENGTH N CONTAINING
C                           THE CROSS TIME SERIES.
C                         NOTE THAT X AND Y MAY BE THE SAME TIME SERIES.
C                N      - INPUT LENGTH OF THE TIME SERIES.
C                XBAR   - INPUT MEAN OF TIME SERIES X.
C                YBAR   - INPUT MEAN OF TIME SERIES Y.
C                MLAG   - INPUT INDICATING THE CROSS-COVARIANCE IS CAL-
C                           CULATED WITH X LAGGING Y BY (-MLAG) UNITS.
C                           MLAG MAY BE POSITIVE, ZERO, OR NEGATIVE.
C                           N MINUS THE ABSOLUTE VALUE OF MLAG MUST BE
C                           POSITIVE.
C                IRS    - INPUT SCALING PARAMETER.
C                           IRS ALLOWS FLEXIBILITY AS A SCALING
C                           PARAMETER. IRS=N IS A COMMON CHOICE.
C                C      - OUTPUT CROSS-COVARIANCE AT LAG MLAG.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER=129 INDICATES THAT N MINUS THE ABSOLUTE
C                             VALUE OF MLAG IS NOT POSITIVE.
C
C   REQD. IMSL ROUTINES - UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE FTCRXY (X,Y,N,XBAR,YBAR,MLAG,IRS,C,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N,MLAG,IRS,IER
      REAL               X(N),Y(N),XBAR,YBAR,C
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,IL,IU
      DOUBLE PRECISION   TEMP,XN,XM,XB,YB,ZERO
      DATA               ZERO/0.0D0/
C                                  FIRST EXECUTABLE STATEMENT
      IER = 129
      IF (N-IABS(MLAG).LT.1) GO TO 9000
      IER = 0
C                                  CALCULATE LIMITS
      IL = 1
      IF (MLAG .GT. 0) IL = IL + MLAG
      IU = N
      IF (MLAG .LT. 0) IU = IU + MLAG
C                                  CALCULATE THE CROSS-COVARIANCE
      XB = XBAR
      YB = YBAR
      TEMP = ZERO
      DO 5 I = IL,IU
         TEMP = TEMP + (DBLE(X(I-MLAG))-XB)*(DBLE(Y(I))-YB)
    5 CONTINUE
      XN = N - IRS + 1
      XM = N - MLAG
      IF (MLAG .LT. 0) XM = N + MLAG
      C = XN * TEMP / XM
      GO TO 9005
 9000 CONTINUE
      CALL UERTST (IER,'FTCRXY')
 9005 RETURN
      END
