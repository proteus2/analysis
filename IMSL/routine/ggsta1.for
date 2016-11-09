C   IMSL ROUTINE NAME   - GGSTA1
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1982
C
C   PURPOSE             - NUCLEUS CALLED ONLY BY IMSL ROUTINE GGSTA
C
C   USAGE               - FUNCTION GGSTA1(XARG)
C
C   ARGUMENTS    GGSTA1 - RESULTANT MODIFIED TANGENT FUNCTION
C                XARG   - DOUBLE PRECISION INPUT FOR FUNCTION
C
C   REQD. IMSL ROUTINES - NONE
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C
C-----------------------------------------------------------------------
C
      DOUBLE PRECISION   FUNCTION GGSTA1( XARG )
C                                  SPECIFICATIONS FOR ARGUMENTS
      DOUBLE PRECISION   XARG
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      DOUBLE PRECISION   PI,PIBY2,PIBY4,P0,P1,P2,Q0,Q1,Q2
      DOUBLE PRECISION   DZERO,DONE,TAN1,X,XX
      LOGICAL            NEG,INV
C                                  APPROX. 4283 FROM HART ET.AL.
C                                                   ( 1968, P.251 )
      DATA               P0,P1,P2,Q0,Q1,Q2
     1                   /.129221035D+3,-.887662377D+1,
     2                    .528644456D-1,.164529332D+3,
     3                   -.451320561D+2,1.0D0/
      DATA               PIBY4/.7853981633974483D0/
      DATA               PIBY2/1.570796326794897D0/
      DATA               PI/3.141592653589793D0/
      DATA               DZERO/0.0D0/,DONE/1.0D0/
C                                  FIRST EXECUTABLE STATEMENT
      X = DABS(XARG)
      IF(X .GT. PIBY4) GO TO 5
      X = X/PIBY4
C                                  CONVERT TO RANGE OF RATIONAL
      XX = X*X
      GGSTA1 = (P0+XX*(P1+XX*P2))/(PIBY4*(Q0+XX*(Q1+XX*Q2)))
      RETURN
C                                  PERFORM RANGE REDUCTION
C                                    ( IF NECESSARY )
    5 NEG = .FALSE.
      INV = .FALSE.
      NEG = XARG.LT.DZERO
      IF(X.LE.PIBY4) GO TO 15
      X = DMOD(X,PI)
      IF(X.LE.PIBY2) GO TO 10
      NEG = .NOT.NEG
      X = PI-X
   10 IF(X.LE.PIBY4) GO TO 15
      INV = .TRUE.
      X = PIBY2-X
   15 X = X/PIBY4
C                                  CONVERT TO RANGE OF RATIONAL
      XX = X*X
      TAN1 = X*(P0+XX*(P1+XX*P2))/(Q0+XX*(Q1+XX*Q2))
      IF(NEG) TAN1 = -TAN1
      IF(INV) TAN1 = DONE/TAN1
      GGSTA1 = TAN1/XARG
      RETURN
      END
