C   IMSL ROUTINE NAME   - GTDDU
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JULY 1, 1983
C
C   PURPOSE             - D-SQUARE TALLY
C
C   USAGE               - CALL GTDDU (R,M,COUNT,K,IER)
C
C   ARGUMENTS    R      - INPUT VECTOR OF LENGTH M CONTAINING THE SE-
C                           QUENCE OF UNIFORM RANDOM NUMBERS.
C                M      - INPUT LENGTH OF SEQUENCE R IN CORE FOR THIS
C                           CALL.
C                COUNT  - INPUT/OUTPUT VECTOR OF LENGTH K WHICH CONTAINS
C                           TALLIES. PRIOR TO THE FIRST CALL TO GTDDU,
C                           COUNT SHOULD BE SET TO ZERO. IF MORE THAN
C                           ONE CALL IS NECESSARY TO EXHAUST THE SE-
C                           QUENCE OF RANDOM NUMBERS, THEN COUNT SHOULD
C                           NOT BE MODIFIED BETWEEN CALLS. COUNT CON-
C                           TAINS TALLIES OF DISTANCES (BETWEEN PAIRS OF
C                           POINTS) WHICH OCCURRED IN RESPECTIVE K CA-
C                           TEGORIES, USING R(I),R(I+1) AS THE FIRST
C                           POINT AND R(I+2),R(I+3) AS THE SECOND.
C                           I=1,5,9,...,J. J DOES NOT EXCEED M.
C                K      - INPUT NUMBER OF EQUIPROBABLE SUBDIVISIONS OF
C                           DISTANCE D (IN RANGE 0,2.0**0.5)
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER = 129 INDICATES M LESS THAN 4
C                             (INSUFFICIENT DATA)
C                           IER = 130 INDICATES K LESS THAN 2
C                             (INSUFFICIENT NUMBER OF SUBDIVISIONS OF
C                             DISTANCE)
C
C   REQD. IMSL ROUTINES - UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE GTDDU  (R,M,COUNT,K,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            M,K,IER
      REAL               R(M),COUNT(K)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,J
      REAL               XK,P,D,D1,PI,TWOPT6,PT3,PIM2
      DATA               PI/3.141593/
      DATA               TWOPT6/2.666667/
      DATA               PT3/.3333333/
      DATA               PIM2/1.141593/
C                                  FIRST EXECUTABLE STATEMENT
C                                  M LESS THAN 4  -TERMINAL ERROR
      IF (M .GE. 4) GO TO 10
      IER = 129
      GO TO 9000
C                                  K LESS THAN 2  -TERMINAL ERROR
   10 IF (K .GE. 2) GO TO 15
      IER = 130
      GO TO 9000
C                                  OPERATE ON R*S IN SETS OF 4
   15 IER = 0
      XK = K
      DO 30 I=1,M,4
         IF (I+3 .GT. M) GO TO 9005
C                                  COMPUTE SQUARE OF DISTANCE
         P = (R(I+2)-R(I))**2 + (R(I+3)-R(I+1))**2
         IF (P .GE. 1.0) GO TO 20
C                                  PROB FOR DIST LESS THAN 1.0
      P = P * PI - TWOPT6 * P * SQRT(P)+0.5*P*P
         GO TO 25
C                                  PROB FOR DIST NOT LESS THAN 1.0
   20    D = P - 1.0
         D1 = SQRT(D)
      P=PT3+PIM2*P+4.*D1+TWOPT6*D*D1-P*P*.5-4.*P*ARCOS(1./SQRT(P))
   25    J = K
C                                  CELL NO. FOR PROB LESS THAN 1.0
         IF (P .LT. 1.0) J=P*XK + 1.0
C                                  INCREMENT COUNT(J)
         COUNT(J) = COUNT(J) + 1.0
   30 CONTINUE
      GO TO 9005
 9000 CONTINUE
      CALL UERTST(IER,'GTDDU ')
 9005 RETURN
      END
