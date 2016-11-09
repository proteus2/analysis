C   IMSL ROUTINE NAME   - NDKER
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1981
C
C   PURPOSE             - NONPARAMETRIC PROBABILITY DENSITY FUNCTION
C                           (ONE DIMENSIONAL) ESTIMATION BY THE
C                           KERNEL METHOD.
C
C   USAGE               - CALL NDKER (X,N,H,DEL,XKER,B,NPT,F)
C
C   ARGUMENTS    X      - INPUT VECTOR OF LENGTH N CONTAINING THE
C                           RANDOM SAMPLE.
C                N      - INPUT SIZE OF RANDOM SAMPLE.
C                H      - INPUT WINDOW WIDTH.
C                DEL    - INPUT PARAMETER.  IF THE KERNEL HAS FINITE
C                           SUPPORT AND IF THE USER WISHES TO USE
C                           THIS PROPERTY FOR GREATER EFFICIENCY, DEL
C                           IS A POSITIVE NUMBER SUCH THAT XKER(Y)=0.0
C                           FOR ALL Y FOR WHICH ABS(Y).GT.DEL.
C                           IF NO SUCH DEL EXISTS OR IF THE USER DOES
C                           NOT WISH TO MAKE USE OF THIS PROPERTY,
C                           DEL SHOULD BE ASSIGNED ANY NONPOSITIVE
C                           VALUE.
C                XKER   - KERNEL FUNCTION.  XKER IS A USER-SUPPLIED
C                           EXTERNAL FUNCTION SUCH THAT XKER(Y) IS THE
C                           ORDINATE OF THE KERNEL AT THE POINT Y.
C                           XKER MUST APPEAR IN AN EXTERNAL STATEMENT
C                           IN THE CALLING PROGRAM.  XKER MUST BE TYPED
C                           APPROPRIATELY.  SEE PRECISION/HARDWARE.
C                B      - INPUT VECTOR OF LENGTH NPT CONTAINING THE
C                           ABSCISSAE AT WHICH DENSITY ESTIMATES ARE
C                           DESIRED.   B MUST BE IN ASCENDING ORDER
C                           IF DEL.GT.0.
C                NPT    - INPUT NUMBER OF POINTS AT WHICH IT IS
C                           DESIRED TO ESTIMATE THE DENSITY.
C                F      - OUTPUT VECTOR OF LENGTH NPT CONTAINING
C                           THE ESTIMATES OF THE DENSITY AT THE
C                           POINTS SPECIFIED IN B.
C
C
C   REQD. IMSL ROUTINES - SINGLE/NONE REQUIRED
C                       - DOUBLE/VXADD,VXSTO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE NDKER  (X,N,H,DEL,XKER,B,NPT,F)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N,NPT
      REAL               X(N),H,DEL,XKER,B(NPT),F(NPT)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            INDX,I,J,NPTP1
      REAL               BL,BU,DELM,FN,ONE,RN,X1,XX,Y
      DOUBLE PRECISION   FJ
      DATA               ONE /1.0/
C                                  FIRST EXECUTABLE STATEMENT
      RN = N
      FN = ONE/(RN*H)
      IF (DEL.GT.0.) GO TO 15
      DO 10 J=1,NPT
         XX = B(J)
         FJ = 0.D0
         DO 5 I=1,N
            Y = (XX-X(I))/H
            FJ = FJ+DBLE(XKER(Y))
    5    CONTINUE
         F(J) = FJ*FN
   10 CONTINUE
      GO TO 50
   15 CONTINUE
      DELM = -DEL
      DO 20 J=1,NPT
   20 F(J) = 0.0
      BL = B(1)-DEL*H
      BU = B(NPT)+DEL*H
      NPTP1 = NPT+1
      DO 40 I=1,N
         IF (X(I).LT.BL .OR. X(I).GT.BU) GO TO 40
         IF ((X(I)-BL).LT.(BU-X(I))) GO TO 30
         DO 25 J=1,NPT
            INDX = NPTP1-J
            XX = B(INDX)
            Y = (XX-X(I))/H
            IF (Y.GT.DEL) GO TO 25
            IF (Y.LT.DELM) GO TO 40
            F(INDX) = F(INDX)+XKER(Y)
   25    CONTINUE
         GO TO 40
   30    CONTINUE
         DO 35 J=1,NPT
            XX = B(J)
            Y = (XX-X(I))/H
            IF (Y.GT.DEL) GO TO 40
            IF (Y.LT.DELM) GO TO 35
            F(J) = F(J)+XKER(Y)
   35    CONTINUE
   40 CONTINUE
      DO 45 J=1,NPT
         F(J) = F(J)*FN
   45 CONTINUE
   50 RETURN
      END
