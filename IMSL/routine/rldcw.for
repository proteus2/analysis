C   IMSL ROUTINE NAME   - RLDCW
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - VARIANCES OF CODED ORTHOGONAL POLYNOMIAL
C                           REGRESSION COEFFICIENTS. FOR USAGE IN
C                           CONJUNCTION WITH IMSL ROUTINES RLFOTH AND
C                           RLFOTW, AND PROVIDED TO PREPARE INPUT FOR
C                           IMSL ROUTINE RLDCVA.
C
C   USAGE               - CALL RLDCW (SSE,X,W,N,ID,IOPT,A,B,V,P,IP,IER)
C
C   ARGUMENTS    SSE    - INPUT ERROR SUM OF SQUARES (OUTPUT ELEMENT
C                           S(MD+3) FROM RLFOTH OR RLFOTW).
C                X      - INPUT SCALED VECTOR OF LENGTH N CONTAINING
C                           INDEPENDENT VARIABLE SETTINGS (THE OUTPUT
C                           X FROM RLFOTH OR RLFOTW).
C                W      - INPUT VECTOR OF LENGTH N CONTAINING THE
C                           SQUARE ROOTS OF THE WEIGHTS AND REQUIRED
C                           ONLY IF IOPT IS NONZERO (THE OUTPUT W FROM
C                           RLFOTW). OTHERWISE, W IS USED AS A WORK
C                           AREA.
C                N      - INPUT LENGTH OF X. N MUST BE GREATER THAN 1.
C                ID     - INPUT DEGREE OF THE FITTED MODEL (THE OUTPUT
C                           ID FROM RLFOTH OR RLFOTW).
C                IOPT   - INPUT OPTION PARAMETER INDICATING THE SOURCE
C                           OF INPUT FOR THIS ROUTINE.
C                           IF IOPT IS ZERO, THE INPUT IS ACQUIRED
C                             FROM RLFOTH.
C                           IF IOPT IS NONZERO, THE INPUT SOURCE IS
C                             RLFOTW.
C                A      - INPUT VECTOR OF LENGTH ID CONTAINING
C                           ORTHOGONAL POLYNOMIAL GENERATION CONSTANTS
C                           (OUTPUT A FROM RLFOTW) AND USED ONLY IF
C                           IOPT IS NONZERO.
C                B      - INPUT VECTOR OF LENGTH ID CONTAINING
C                           ADDITIONAL CONSTANTS. SEE DESCRIPTION OF
C                           PARAMETER A ABOVE.
C                V      - OUTPUT VECTOR OF LENGTH ID+1 CONTAINING THE
C                           VARIANCES OF THE CODED INTERCEPT AND
C                           REGRESSION COEFFICIENTS REQUIRED AS INPUT
C                           TO RLDCVA.
C                P      - OUTPUT N BY ID MATRIX CONTAINING THE SETTINGS
C                           FOR THE ID ORTHOGONAL POLYNOMIALS.
C                IP     - INPUT ROW DIMENSION OF MATRIX P EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER = 129 INDICATES A TERMINAL ERROR
C                             OCCURRED IN IMSL ROUTINE RLPOL.
C                           IER = 130 INDICATES N WAS LESS THAN OR
C                             EQUAL TO 1.
C                         WARNING ERROR
C                           IER = 34 INDICATES A WARNING ERROR
C                             OCCURRED IN IMSL ROUTINE RLPOL.
C
C   REQD. IMSL ROUTINES - SINGLE/RLPOL,UERTST,UGETIO
C                       - DOUBLE/RLPOL,UERTST,UGETIO,VXADD,VXMUL,
C                           VXSTO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE RLDCW  (SSE,X,W,N,ID,IOPT,A,B,V,P,IP,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N,ID,IOPT,IP,IER
      REAL               SSE,X(N),W(N),A(ID),B(ID),V(1),P(IP,ID)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            NMID1,IDP1,I,JM1,JM2,J,JP1
      REAL               SM,SA,ZERO,ONE,SSQRD
      DOUBLE PRECISION   TEMP,DZERO,TOTW
      DATA               DZERO/0.0D0/
      DATA               ZERO/0.0/,ONE/1.0/
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
      IF (N.GT.1) GO TO 5
      IER = 130
      GO TO 9000
    5 NMID1 = N-ID-1
      IF (NMID1.NE.0) GO TO 15
      IDP1 = ID+1
      DO 10 I=1,IDP1
         V(I) = ZERO
   10 CONTINUE
      GO TO 9005
   15 SSQRD = SSE/NMID1
      IF (IOPT.EQ.0) GO TO 40
      TOTW = DZERO
      IF (ID.EQ.0) GO TO 30
      B(1) = ZERO
      DO 25 I=1,N
         TOTW = TOTW+DBLE(W(I))*DBLE(W(I))
         P(I,ID) = ONE
         JM1 = ID
         JM2 = ID
         DO 20 J=1,ID
            P(I,J) = (X(I)-A(J))*P(I,JM1)-B(J)*P(I,JM2)
            JM2 = JM1
            JM1 = J
   20    CONTINUE
   25 CONTINUE
      GO TO 50
   30 DO 35 I=1,N
         TOTW = TOTW+DBLE(W(I))*DBLE(W(I))
   35 CONTINUE
      GO TO 65
   40 TOTW = N
      DO 45 I=1,N
         W(I) = ONE
   45 CONTINUE
      IF (ID.EQ.0) GO TO 65
      CALL RLPOL (X,N,ID,SM,SA,V,V,P,IP,IER)
      IF (IER.LT.129) GO TO 50
      IER = 129
      GO TO 9000
   50 IF (ID.EQ.0) GO TO 65
      DO 60 J=1,ID
         TEMP = DZERO
         DO 55 I=1,N
            TEMP = TEMP+DBLE(P(I,J))*DBLE(P(I,J))*DBLE(W(I))*DBLE(W(I))
   55    CONTINUE
         JP1 = J+1
         V(JP1) = SSQRD/TEMP
   60 CONTINUE
   65 V(1) = SSQRD/TOTW
      IF (IER.EQ.0) GO TO 9005
      IER = 34
 9000 CONTINUE
      CALL UERTST (IER,6HRLDCW )
 9005 RETURN
      END
