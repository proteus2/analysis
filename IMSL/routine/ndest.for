C   IMSL ROUTINE NAME   - NDEST
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1980
C
C   PURPOSE             - EVALUATE PROBABILITY DENSITY FUNCTION
C                           AT SPECIFIED POINTS
C
C   USAGE               - CALL NDEST (Y,N,F,M,IOPT,B,WK,EST,IER)
C
C   ARGUMENTS    Y      - INPUT VECTOR OF LENGTH N CONTAINING
C                           THE POINTS AT WHICH THE PROBABILITY
C                           DENSITY IS DESIRED.
C                N      - INPUT NUMBER OF POINTS AT WHICH DENSITY
C                           IS DESIRED.
C                F      - INPUT VECTOR OF LENGTH M CONTAINING THE
C                           GIVEN DENSITY FUNCTION VALUES.
C                M      - INPUT NUMBER OF ORDINATES SUPPLIED.  M
C                           MUST BE GREATER THAN 1 IF LINEAR
C                           INTERPOLATION IS DESIRED AND MUST BE
C                           GREATER THAN 3 IF QUASI-CUBIC
C                           INTERPOLATION IS DESIRED.
C                IOPT   - INPUT OPTION PARAMETER.  IF
C                           IOPT=1, LINEAR INTERPOLATION IS DESIRED
C                             AND THE DATA ARE EQUALLY SPACED.
C                           IOPT=2, LINEAR INTERPOLATION IS DESIRED
C                             AND THE DATA ARE NOT NECESSARILY
C                             EQUALLY SPACED.
C                           IOPT=3, QUASI-CUBIC HERMITE INTERPOLATION
C                             IS DESIRED AND THE DATA ARE EQUALLY
C                             SPACED.
C                           IOPT=4, QUASI-CUBIC HERMITE INTERPOLATION
C                             IS DESIRED AND THE DATA ARE NOT
C                             NECESSARILY EQUALLY SPACED.
C                B      - INPUT ARRAY OF LENGTH M IF IOPT=2,3 OR 4.
C                           OTHERWISE B IS OF LENGTH 2.  IF IOPT=2
C                           OR 4, B(I) IS THE ABSCISSA CORRESPONDING
C                           TO F(I).  THE ABSCISSAE MUST BE SPECIFIED
C                           IN INCREASING ORDER.  IF IOPT=1 OR 3,
C                           B(1) IS THE LOWER ENDPOINT OF THE SUPPORT
C                           OF THE DISTRIBUTION (CORRESPONDING TO
C                           F(1)) AND B(2) IS THE UPPER ENDPOINT
C                           (CORRESPONDING TO F(M)).  IF IOPT=3, THE
C                           EXTRA M-2 POSITIONS IN B ARE WORKSPACE.
C                WK     - WORK ARRAY.  IF IOPT=1 OR 2, ONLY WK(1)
C                           IS NEEDED.  OTHERWISE, WK MUST BE OF
C                           LENGTH 3*M.
C                EST    - OUTPUT ARRAY OF LENGTH N CONTAINING THE
C                           DENSITY FUNCTION ESTIMATES AT THE POINTS
C                           GIVEN IN Y.
C                IER    - ERROR PARAMETER.  (OUTPUT)
C                         TERMINAL ERROR
C                           IER=130 INDICATES M IS LESS THAN 4
C                             AND IOPT=3 OR 4, OR M IS LESS THAN
C                             2 AND IOPT=1 OR 2.
C                           IER=131 INDICATES INPUT ABSCISSAE ARE
C                             NOT IN ASCENDING ORDER.
C                           IER=132 INDICATES IOPT IS OUT OF RANGE.
C
C
C   REQD. IMSL ROUTINES - IQHSCU,UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C----------------------------------------------------------------------
C
      SUBROUTINE NDEST  (Y,N,F,M,IOPT,B,WK,EST,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N,M,IOPT,IC,IER
      REAL               Y(N),F(M),B(1),WK(1),EST(N)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,JPTR,J,M2
      REAL               BLEFT,BRIGHT,DEL,H,YI,ZERO
      DATA               ZERO/0.0/
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
      M2 = 2*M
      IF (M.GT.1) GO TO 5
      IER = 130
      GO TO 9000
    5 IF (IOPT.GE.1.AND.IOPT.LE.4) GO TO 10
      IER = 132
      GO TO 9000
   10 BLEFT = B(1)
      BRIGHT = B(2)
      IF (BRIGHT.GT.BLEFT) GO TO 15
      IER = 131
      GO TO 9000
   15 GO TO (20,35,55,55), IOPT
C                                  LINEAR CASE, EQUAL SPACING
   20 H = (BRIGHT-BLEFT)/(M-1)
      DO 30 I=1,N
         YI = Y(I)
         EST(I) = ZERO
         IF (YI.LT.BLEFT.OR.YI.GT.BRIGHT) GO TO 30
         IF (YI.LT.BRIGHT) GO TO 25
         EST(I) = F(M)
         GO TO 30
   25    DEL = YI-BLEFT
         J = DEL/H
         JPTR = J+1
         DEL = DEL-J*H
         EST(I) = F(JPTR)+(F(JPTR+1)-F(JPTR))*DEL/H
   30 CONTINUE
      GO TO 9005
   35 BRIGHT = B(M)
C                                  LINEAR CASE, UNEQUAL SPACING
      DO 40 J=2,M
         IF (B(J).GT.B(J-1)) GO TO 40
         IER = 131
         GO TO 9000
   40 CONTINUE
      DO 50 I=1,N
         YI = Y(I)
         EST(I) = ZERO
         IF (YI.LT.BLEFT.OR.YI.GT.BRIGHT) GO TO 50
         J = 1
   45    J = J+1
         IF (YI.GT.B(J)) GO TO 45
         EST(I) = F(J)-(F(J)-F(J-1))*(B(J)-YI)/(B(J)-B(J-1))
   50 CONTINUE
      GO TO 9005
C                                  QUASI-CUBIC CASE
   55 IF (M.GT.3) GO TO 60
      IER = 130
      GO TO 9000
   60 IF (IOPT.EQ.4) GO TO 70
C                                  EQUAL SPACING
      H = (BRIGHT-BLEFT)/(M-1)
      DO 65 J=2,M
   65 B(J) = B(J-1)+H
C                                  DO QUASI-CUBIC FIT
   70 CALL IQHSCU (B,F,M,WK,M,IER)
C                                  NOW DO INTERPOLATION
      IF (IOPT.EQ.4) GO TO 80
C                                  EQUAL SPACING
      DO 75 I=1,N
         YI = Y(I)
         EST(I) = ZERO
         IF (YI.LT.BLEFT.OR.YI.GT.BRIGHT) GO TO 75
         DEL = YI-BLEFT
         J = DEL/H
         JPTR = J+1
         DEL = DEL-J*H
         EST(I) = ((WK(M2+JPTR)*DEL+WK(M+JPTR))*DEL+WK(JPTR))*DEL
     1   +F(JPTR)
   75 CONTINUE
      B(2) = BRIGHT
      GO TO 9005
C                                  UNEQUAL SPACING
   80 BRIGHT = B(M)
      DO 90 I=1,N
         YI = Y(I)
         EST(I) = ZERO
         IF (YI.LT.BLEFT.OR.YI.GT.BRIGHT) GO TO 90
         J = 1
   85    J = J+1
         IF (YI.GT.B(J)) GO TO 85
         J = J-1
         DEL = YI-B(J)
         EST(I) = ((WK(M2+J)*DEL+WK(M+J))*DEL+WK(J))*DEL+F(J)
   90 CONTINUE
      GO TO 9005
 9000 CONTINUE
      CALL UERTST (IER,6HNDEST )
 9005 RETURN
      END

