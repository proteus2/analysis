C   IMSL ROUTINE NAME   - RLGQMO
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - CENTERING OF INDEPENDENT VARIABLE SETTINGS AND
C                           GENERATION OF UNCENTERED SQUARE AND CROSS
C                           PRODUCT TERMS - OUT OF CORE VERSION
C
C   USAGE               - CALL RLGQMO (X,N,M,I,XBAR,IER)
C
C   ARGUMENTS    X      - INPUT AND OUTPUT VECTOR OF LENGTH M*(M+3)/2.
C                         ON INPUT, X IS A VECTOR OF LENGTH M
C                           CONTAINING THE I-TH ROW OF THE N X M
C                           MATRIX OF INDEPENDENT VARIABLE SETTINGS.
C                         ON OUTPUT, X IS A VECTOR OF LENGTH M*(M+3)/2
C                           CONTAINING (AFTER THE I-TH CALL) THE
C                           CENTERED SETTING FOR EACH OF THE M INPUT
C                           VARIABLES AND THE SETTINGS FOR THE M*(M+1)/2
C                           GENERATED AND UNCENTERED SQUARE AND CROSS
C                           PRODUCT VARIABLES FOR THE I-TH ROW.
C                N      - INPUT NUMBER OF SETTINGS FOR EACH INDEPENDENT
C                           VARIABLE.
C                M      - INPUT LENGTH OF THE X VECTOR ON INPUT.
C                I      - INPUT NUMBER OF THE ROW OF THE MATRIX OF
C                           INDEPENDENT VARIABLE SETTINGS ENTERED TO
C                           RLGQMO. I MUST BE GREATER THAN OR EQUAL TO
C                           ONE AND LESS THAN OR EQUAL TO N.
C                XBAR   - INPUT AND OUTPUT VECTOR OF LENGTH M*(M+3)/2.
C                         ON INPUT, XBAR IS A VECTOR OF LENGTH M
C                           CONTAINING THE MEANS OF THE M INDEPENDENT
C                           VARIABLES.
C                         ON OUTPUT (AFTER THE N-TH CALL), XBAR IS
C                           VECTOR OF LENGTH M*(M+3)/2, WHOSE FIRST M
C                           COMPONENTS ARE UNCHANGED AND WHOSE
C                           REMAINING M*(M+1)/2 COMPONENTS CONTAIN THE
C                           MEANS OF THE GENERATED SQUARE AND CROSS
C                           PRODUCT VARIABLES.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER=129 INDICATES THAT I IS LESS THAN 1 OR
C                             GREATER THAN N.
C
C   REQD. IMSL ROUTINES - UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS      THE INPUT VARIABLES MAY BE CENTERED OR UNCENTERED.
C                THE INPUT MEANS MUST BE THE MEANS OF THESE VARIABLES.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE RLGQMO (X,N,M,I,XBAR,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N,M,I,IER
      REAL               X(1),XBAR(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            J,JJ,K,L,L1
C                                  FIRST EXECUTABLE STATEMENT
      IER=0
      IF (I .GE. 1 .AND. I .LE. N) GO TO 5
C                                  TERMINAL ERROR - I IS LESS THAN 1
C                                  OR GREATER THAN N
      IER=129
      GO TO 9000
    5 IF (I .NE. 1) GO TO 15
      J=M+1
      K=(M*(M+3))/2
         DO 10 JJ=J,K
   10    XBAR(JJ)=0.
   15    DO 20 J=1,M
            X(J)=X(J)-XBAR(J)
   20    CONTINUE
      L1=M+M
         DO 35 J=1,M
            DO 35 JJ=J,M
            IF (J .EQ. JJ) GO TO 25
            L1=L1+1
            L=L1
            GO TO 30
   25       L=M+J
   30       X(L)=X(J)*X(JJ)
   35       XBAR(L)=XBAR(L)+X(L)
      IF (I .NE. N) GO TO 9005
      L=M+1
      L1=(M*(M+3))/2
         DO 40 J=L,L1
   40    XBAR(J)=XBAR(J)/N
      GO TO 9005
 9000 CONTINUE
       CALL UERTST(IER,6HRLGQMO)
 9005 RETURN
      END
