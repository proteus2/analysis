C   IMSL ROUTINE NAME   - EBALAC
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1982
C
C   PURPOSE             - NUCLEUS CALLED ONLY BY IMSL ROUTINE EIGCC
C
C   REQD. IMSL ROUTINES - NONE REQUIRED
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE EBALAC (AR,AI,N,IA,K,L,D)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N,IA,K,L
      REAL               AR(IA,1),AI(IA,1),D(N)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            J,M,I,IEXC,L1,JJ
      REAL               RADIX,ZERO,ONE,PT95,B2,F,C,G,R,S,RRADIX,RB2
      LOGICAL            NOCONV
C                                  RADIX IS A MACHINE DEPENDENT
C                                    PARAMETER SPECIFYING THE BASE OF
C                                    THE MACHINE FLOATING POINT REPRE-
C                                    SENTATION
      DATA               RADIX/16.0/
      DATA               ZERO,ONE,PT95/0.0,1.0,0.95/
C                                  FIRST EXECUTABLE STATEMENT
      B2 = RADIX*RADIX
      RRADIX = ONE/RADIX
      RB2 = RRADIX*RRADIX
      K = 1
      L = N
      GO TO 30
C                                  IN-LINE PROCEDURE FOR ROW AND COLUMN
C                                    EXCHANGE
    5 D(M) = J
      IF (J .EQ. M) GO TO 20
      DO 10 I = 1,L
         F = AR(I,J)
         AR(I,J) = AR(I,M)
         AR(I,M) = F
         F = AI(I,J)
         AI(I,J) = AI(I,M)
         AI(I,M) = F
   10 CONTINUE
      DO 15 I=K,N
         F = AR(J,I)
         AR(J,I) = AR(M,I)
         AR(M,I) = F
         F = AI(J,I)
         AI(J,I) = AI(M,I)
         AI(M,I) = F
   15 CONTINUE
   20 GO TO (25,45), IEXC
C                                  SEARCH FOR ROWS ISOLATING AN
C                                    EIGENVALUE AND PUSH THEM DOWN
   25 IF (L .EQ. 1) GO TO 115
      L = L-1
C                                  DO J=L,1,-1
   30 L1 = L+1
      DO 40 JJ = 1,L
         J = L1-JJ
         DO 35 I = 1,L
            IF (I .EQ. J) GO TO 35
            IF (AR(J,I) .NE. ZERO .OR. AI(J,I) .NE. ZERO) GO TO 40
   35    CONTINUE
         M = L
         IEXC = 1
         GO TO 5
   40 CONTINUE
      GO TO 50
C                                  SEARCH FOR COLUMNS ISOLATING AN
C                                    EIGENVALUE AND PUSH THEM LEFT
   45 K = K+1
   50 DO 60 J = K,L
         DO 55 I = K,L
            IF (I .EQ. J) GO TO 55
            IF (AR(I,J) .NE. ZERO .OR. AI(I,J) .NE. ZERO) GO TO 60
   55    CONTINUE
         M = K
         IEXC = 2
         GO TO 5
   60 CONTINUE
C                                  BALANCE THE SUBMATRIX IN ROWS
C                                    K TO L
      DO 65 I = K,L
         D(I) = ONE
   65 CONTINUE
C                                  ITERATIVE LOOP FOR NORM REDUCTION
   70 NOCONV = .FALSE.
      DO 110 I = K,L
         C = ZERO
         R = ZERO
         DO 75 J = K,L
            IF (J .EQ. I) GO TO 75
            C = C+ABS(AR(J,I))+ABS(AI(J,I))
            R = R+ABS(AR(I,J))+ABS(AI(I,J))
   75    CONTINUE
         G = R*RRADIX
         F = ONE
         S = C+R
   80    IF (C .GE. G) GO TO 85
         F = F*RADIX
         C = C*B2
         GO TO 80
   85    G = R*RADIX
   90    IF (C .LT. G) GO TO 95
         F = F*RRADIX
         C = C*RB2
         GO TO 90
C                                  BALANCE
   95    IF ((C+R)/F .GE. PT95*S) GO TO 110
         G = ONE/F
         D(I) = D(I)*F
         NOCONV = .TRUE.
         DO 100 J = K,N
            AR(I,J) = AR(I,J)*G
            AI(I,J) = AI(I,J)*G
  100    CONTINUE
         DO 105 J = 1,L
            AR(J,I) = AR(J,I)*F
            AI(J,I) = AI(J,I)*F
  105    CONTINUE
  110 CONTINUE
      IF (NOCONV) GO TO 70
  115 RETURN
      END
