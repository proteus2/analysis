C   IMSL ROUTINE NAME   - RLFITO
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - PURE REPLICATION ERROR DEGREES OF FREEDOM
C                           AND SUM OF SQUARES - OUT OF CORE VERSION
C
C   USAGE               - CALL RLFITO (X,N,M,I,Y,L,T,W,S,ND,IER)
C
C   ARGUMENTS    X      - I-TH ROW OF THE N BY M DESIGN MATRIX. (INPUT)
C                           X IS A VECTOR OF LENGTH M.
C                N      - NUMBER OF ROWS IN THE DESIGN MATRIX AND THE
C                           RESPONSE MATRIX. (INPUT)
C                M      - NUMBER OF COLUMNS IN THE DESIGN MATRIX.
C                           (INPUT)
C                I      - NUMBER OF THE ROW OF THE DESIGN MATRIX AND
C                           OF THE RESPONSE MATRIX ENTERED TO RLFITO.
C                           (INPUT) I MUST BE GREATER THAN OR EQUAL TO
C                           ONE AND LESS THAN OR EQUAL TO N.
C                Y      - I-TH ROW OF THE N BY L RESPONSE MATRIX.
C                           (INPUT) Y IS A VECTOR OF LENGTH L.
C                L      - NUMBER OF COLUMNS IN THE RESPONSE MATRIX.
C                           (INPUT)
C                T      - WORK AREA OF LENGTH M. T MUST NOT BE ALTERED
C                           BY THE USER BETWEEN CALLS TO RLFITO.
C                W      - DOUBLE PRECISION WORK AREA OF LENGTH 3*L.
C                           W MUST NOT BE ALTERED BY THE USER BETWEEN
C                           CALLS TO RLFITO.
C                S      - VECTOR OF LENGTH L CONTAINING THE PURE ERROR
C                           SUM OF SQUARES FOR EACH RESPONSE, AFTER
C                           RLFITO IS CALLED FOR THE N-TH TIME. (OUTPUT)
C                ND     - NUMBER OF DEGREES OF FREEDOM FOR THE PURE
C                           ERROR SUM OF SQUARES, AFTER RLFITO IS
C                           CALLED FOR THE N-TH TIME. (OUTPUT)
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER=129 INDICATES THAT I IS GREATER THAN N
C                             OR I IS LESS THAN 1.
C                         WARNING ERROR
C                           IER=34 INDICATES THAT NO REPLICATES APPEARED
C
C   REQD. IMSL ROUTINES - UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE RLFITO (X,N,M,I,Y,L,T,W,S,ND,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N,M,I,L,ND,IER
      REAL               X(1),Y(1),T(1),S(1)
      DOUBLE PRECISION   W(3,1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            ICT,J
      DOUBLE PRECISION   Z
      DATA               ICT/0/
C                                  FIRST EXECUTABLE STATEMENT
      IER=0
      IF (I .LE. N .AND. I .GE. 1) GO TO 5
C                                  TERMINAL ERROR - I IS GREATER THAN
C                                  N OR LESS THAN 1
      IER=129
      GO TO 9000
    5 IF (I .NE. 1) GO TO 30
      ND=0
      DO 10 J=1,L
   10    S(J)=0.0
   15 DO 20 J=1,L
         W(1,J)=Y(J)
         W(2,J)=Y(J)
         W(3,J)=0.0D0
   20 CONTINUE
       DO 25 J=1,M
   25    T(J)=X(J)
      ICT=0
      IF (I .LT. N) GO TO 9005
      GO TO 60
C                                  SEARCH FOR REPLICATE
   30 DO 35 J=1,M
         IF (X(J) .NE. T(J)) GO TO 45
   35 CONTINUE
C                                  A REPLICATE IS FOUND
      DO 40 J=1,L
         W(2,J)=W(2,J)+Y(J)
         Z=Y(J)-W(1,J)
         W(3,J)=W(3,J)+Z*Z
   40 CONTINUE
      ICT=ICT+1
      IF (I .LT. N) GO TO 9005
C                                  ROW IS NOT A REPLICATE
   45 IF (ICT .LE. 0) GO TO 15
      ND=ND+ICT
      ICT=ICT+1
      DO 50 J=1,L
   50    W(2,J)=W(2,J)/ICT
      DO 55 J=1,L
         Z=W(2,J)-W(1,J)
         S(J)=S(J)+W(3,J)-ICT*Z*Z
   55 CONTINUE
      GO TO 15
   60 IF (ND .NE. 0) GO TO 9005
      IER=34
 9000 CONTINUE
      CALL UERTST (IER,6HRLFITO)
 9005 RETURN
      END
