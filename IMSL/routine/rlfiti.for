C   IMSL ROUTINE NAME   - RLFITI
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - PURE REPLICATION ERROR DEGREES OF FREEDOM
C                           AND SUM OF SQUARES - IN CORE VERSION
C
C   USAGE               - CALL RLFITI (X,N,M,IX,Y,L,IY,SS,NDF,IER)
C
C   ARGUMENTS    X      - N BY M MATRIX OF INDEPENDENT VARIABLE
C                           SETTINGS. (INPUT)
C                N      - NUMBER OF ROWS IN X AND Y. (INPUT)
C                           N MUST BE LESS THAN OR EQUAL TO IX AND IY.
C                M      - NUMBER OF COLUMNS IN X. (INPUT)
C                IX     - ROW DIMENSION OF MATRIX X EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT
C                           IN THE CALLING PROGRAM. (INPUT)
C                Y      - N BY L MATRIX OF OBSERVATIONS ON L RESPONSE
C                           VARIABLES. (INPUT)
C                L      - NUMBER OF COLUMNS IN Y. (INPUT)
C                IY     - ROW DIMENSION OF MATRIX Y EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT
C                           IN THE CALLING PROGRAM. (INPUT)
C                SS     - VECTOR OF LENGTH L CONTAINING THE PURE
C                           ERROR SUM OF SQUARES FOR EACH RESPONSE.
C                           (OUTPUT)
C                NDF    - NUMBER OF DEGREES OF FREEDOM FOR THE PURE
C                           ERROR SUM OF SQUARES. (OUTPUT)
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         WARNING ERROR
C                           IER = 33 INDICATES THAT NO REPLICATES
C                             APPEARED IN THE X MATRIX
C                           IER = 34 INDICATES THAT THE PURE ERROR SUM
C                             OF SQUARES FOR ONE OR MORE OF THE
C                             RESPONSES IS EQUAL TO ZERO.
C
C   REQD. IMSL ROUTINES - SINGLE/UERTST,UGETIO
C                       - DOUBLE/UERTST,UGETIO,VXADD,VXMUL,VXSTO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE RLFITI (X,N,M,IX,Y,L,IY,SS,NDF,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N,M,IX,L,IY,NDF,IER
      REAL               X(IX,1),Y(IY,1),SS(L)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,ISTART,ICT,K,IEND,LL,J
      REAL               ZERO
      DOUBLE PRECISION   TEMP,D,DD,DZERO
      DATA               ZERO/0.0/,DZERO/0.0D0/
C                                  FIRST EXECUTABLE STATEMENT
      IER=0
      DO 5 I=1,L
         SS(I)=ZERO
    5 CONTINUE
      ISTART=1
      ICT=0
      NDF=0
C                                  SEARCH FOR GROUP OF REPLICATES
      DO 40 I=2,N
         DO 10 K=1,M
            IF (X(I,K) .NE. X(ISTART,K)) GO TO 15
   10    CONTINUE
C                                  A REPLICATE IS FOUND
         ICT=ICT+1
         IF (I .LT. N) GO TO 40
C                                  ROW IS NOT A REPLICATE
   15    IF (ICT .LE. 0) GO TO 35
         IEND=ISTART+ICT
         NDF=NDF+ICT
C                                  CALCULATE MEAN OF RESPONSES FOR
C                                  CURRENT GROUP OF REPLICATES
         DO 30 LL=1,L
            TEMP = DZERO
            DO 20 J=ISTART,IEND
               TEMP=TEMP+Y(J,LL)
   20       CONTINUE
            TEMP=TEMP/(ICT+1)
C                                  CALCULATE PURE ERROR SUM OF SQUARES
C                                  FOR EACH RESPONSE IN THE REPLICATE
C                                  GROUP
            DD = DZERO
            DO 25 J=ISTART,IEND
               D=Y(J,LL)-TEMP
               DD=DD+D*D
   25       CONTINUE
            SS(LL)=SS(LL)+DD
   30    CONTINUE
   35    ISTART=I
         ICT=0
   40 CONTINUE
C                                  NO REPLICATES APPEARED
      IF (NDF .NE. 0) GO TO 45
      IER=33
      GO TO 9000
C                                  CHECK FOR ZERO SUM OF SQUARES
   45 DO 50 I=1,L
         IF (SS(I) .NE. ZERO) GO TO 50
         IER = 34
         GO TO 9000
   50 CONTINUE
      GO TO 9005
 9000 CONTINUE
      CALL UERTST (IER,6HRLFITI)
 9005 RETURN
      END
