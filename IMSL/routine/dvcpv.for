C   IMSL ROUTINE NAME   - DVCPV
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1982
C
C   PURPOSE             - NUCLEUS CALLED ONLY BY IMSL SUBROUTINE DVCPR
C
C   REQD. IMSL ROUTINES - DVCPW,DVCPX,DVCPY,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE DVCPV (M,N,P,R,ALPHA,A1,B1,X,Y,LIN,A2,C2,DEL,JERROR,
     *                   IPRINT,EPS,IR,IC,UU,RES,MMAX,MTNMAX,NMAX,MMAX2,
     *                   F,HX,SK,GRADF,AUX,ICA,XAU,FCNI,FCNJ,FCNB)
C                                  SPECIFICATIONS FOR ARGUMENTS
C
      INTEGER            M,N,P,R,JERROR,IPRINT,MMAX,MTNMAX,NMAX,MMAX2,
     *                   IR(NMAX,MMAX),IC(NMAX,MMAX),ICA(MMAX),NIN,NOUT
      REAL               ALPHA(MMAX),A1(MMAX,MMAX),B1(MMAX,MMAX),X(NMAX)
     *                   ,Y(MTNMAX),A2(MTNMAX,MMAX),C2(MTNMAX,MMAX),
     *                   DEL(MMAX,MTNMAX),EPS,UU(MTNMAX),RES(MTNMAX),
     *                   F(MTNMAX),HX(NMAX),SK(MTNMAX),GRADF(MTNMAX),
     *                   AUX(MMAX2,MMAX),XAU(MTNMAX)
      LOGICAL            LIN
C                                  SPECIFICATIONS FOR COMMON /C1 /
      COMMON /C1/        EPSNU,CONT
      REAL               EPSNU
      LOGICAL            CONT
C                                  SPECIFICATIONS FOR COMMON /NEWT /
      COMMON /NEWT/      INWT,NU,CASI
      INTEGER            INWT,NU
      LOGICAL            CASI
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
C
      INTEGER            I1,I7,I8,I,J,K1JM,K1J,K1,MPNM,MPN,MP,P1
      REAL               DT,GNOR,PNOR,RABS1,RABS,REOLD,ROL1,ROL,SCPR,
     *                   SUM,TE,TMIN,TN,T,TPM8
      LOGICAL            SING,STEP
      EXTERNAL           FCNJ,FCNB
C                                  ERROR EXIT JERROR = 3 NON
C                                    CONVERGENCE. NEWTON ITERATION
C                                    STARTS
C                                  FIRST EXECUTABLE STATEMENT
      CALL UGETIO(1,NIN,NOUT)
      DT = 1.
      MP = M-P
      P1 = P+1
      T = 1.
      TE = 1.
      TPM8 = 2.**(-8)
      TMIN = TPM8
      STEP = .FALSE.
      SCPR = 0.
      MPN = M*N
      DO 5 I=1,MPN
    5 UU(I) = 0.
      MPNM = M*(N-1)
      INWT = 0
      ROL = 1.E20
      REOLD = 1.0E20
      IF (IPRINT.EQ.0) GO TO 20
      WRITE (NOUT,10) N
   10 FORMAT (7H NGRID=, I5)
      WRITE (NOUT,15) NU, EPS
   15 FORMAT (12H CORRECTION=, I4, 33H RESIDUAL FOR NEWTON SHOULD BE .L,
     *2HE., E14.7)
C                                  RESIDUAL COMPUTATION
   20 RABS = 0.
      DO 25 I7=1,N
         I8 = (I7-1)*M+1
         CALL FCNI(M,X(I7),Y(I8),F(I8))
   25 CONTINUE
      CALL FCNB(M,Y(1),Y(I8),ALPHA)
      IF (P.EQ.0) GO TO 35
      DO 30 I=1,P
         RES(I) = -ALPHA(I)
         RABS = RABS+RES(I)**2
   30 CONTINUE
C
   35 DO 45 I=2,N
         I1 = I-1
         K1 = I1*M
         DO 40 J=1,M
            K1J = K1+J
            K1JM = K1J-M
            RES(K1J-MP) = -Y(K1J)+Y(K1JM)+.5*HX(I1)*(F(K1J)+F(K1JM))
     *      +SK(K1JM)
            RABS = RABS+RES(K1J-MP)**2
   40    CONTINUE
   45 CONTINUE
C
      DO 50 J=P1,M
         RES(MPNM+J) = -ALPHA(J)
         RABS = RABS+ALPHA(J)**2
   50 CONTINUE
      RABS1 = SQRT(RABS)
      IF (IPRINT.EQ.0) GO TO 65
      WRITE (NOUT,55) INWT, RABS1
   55 FORMAT (18H NEWTON ITERATION=, I3, 20H MAX. ABS. RESIDUAL=, E14.7)
      WRITE (NOUT,60) (RES(J),J=1,MPN)
   60 FORMAT (11H RESIDUALS=/(1X, 6E12.3))
   65 IF (INWT.EQ.0 .AND. (EPSNU.EQ.1. .OR. EPSNU.EQ.0.)) GO TO 140
C
C                                  CHECK FOR CONVERGENCE
      IF (RABS1.GT.EPS) GO TO 70
      IF (EPSNU.LT.1.) RETURN
C                                  JACOBIAN IS KEPT CONSTANT UNTIL NEXT
C                                    MESH CHANGE. STEP AND ANGLE
C                                    CONTROL ARE SHORTCIRCUITED.
      CASI = .TRUE.
      RETURN
C                                  NEWTON EXIT.
C                                  CONVERGENCE OR TOO MANY ITERATIONS.
C                                  NEWTON TEST IN ORDER TO AVOID CYCLING
C
   70 IF ((REOLD-RABS.GE..5*T*SCPR .AND. INWT.LT.15) .OR. (INWT.EQ.1))
     *GO TO 140
      IF (INWT.EQ.15) GO TO 115
C                                  STEP CONTROL STARTS
      IF (.NOT.(CASI .OR. LIN)) GO TO 80
      IF (LIN) GO TO 75
      CASI = .FALSE.
      GO TO 145
   75 IF (RABS1.LE.100.*EPS) GO TO 120
      GO TO 135
   80 IF (STEP) GO TO 95
      STEP = .TRUE.
      IF (TE.LE.1.) GO TO 85
      TN = TE
      GO TO 90
   85 TN = 100.*TE
      IF (TN.GT.1.) GO TO 95
   90 TMIN = TN*TPM8
      GO TO 100
   95 TN = .5*T
      IF (TN.LT.TMIN) GO TO 115
      ROL = RABS
  100 DT = TN-T
      T = TN
      DO 105 I=1,MPN
  105 Y(I) = Y(I)+DT*UU(I)
      IF (IPRINT.EQ.0) GO TO 20
      WRITE (NOUT,110) T, TMIN
  110 FORMAT (20H STEP CONTROL=T - DT, 2E15.7)
      GO TO 20
C
  115 ROL1 = SQRT(ROL)
      IF (ROL1.GE.100.*EPS) GO TO 135
      RABS1 = ROL1
  120 DO 125 I=1,MPN
  125 Y(I) = Y(I)-DT*UU(I)
      IF (IPRINT.NE.0) WRITE (NOUT,130) NU, RABS1, EPS
  130 FORMAT (20H WARNING--CORRECTION, I4/26H RESIDUAL IN NEWTON COULD ,
     *18HONLY BE REDUCED TO, E12.2, 17H INSTEAD OF BELOW,
     *E12.2/56H IF UPON TERMINATION TOL HAS NOT BEEN REACHED THE PROBLE,
     *1HM, 32H REQUIRES A LONGER COMPUTER WORD)
C                                  NEWTON DID NOT QUITE REACH THE
C                                    TOLERANCE. FURTHER MESH
C                                    REFINEMENTS ARE NOT ALLOWED.
      JERROR = 4
      IF (EPSNU.LT.1.) RETURN
C                                  JACOBIAN IS KEPT CONSTANT UNTIL NEXT
C                                    MESH CHANGE. STEP AND ANGLE
C                                    CONTROL ARE SHORTCIRCUITED.
      CASI = .TRUE.
      RETURN
C                                  WE ASSUME DIVERGENCE AND RETURN
C                                    ERROR EXIT 3
  135 JERROR = 3
      RETURN
C
  140 IF (INWT.EQ.4) CASI = .FALSE.
  145 CALL DVCPW(M,N,P,R,X,Y,A1,B1,A2,C2,DEL,CASI,SING,IR,IC,UU,RES,LIN,
     *MMAX,MTNMAX,NMAX,MMAX2,HX,GRADF,AUX,ICA,XAU,FCNJ,FCNB)
      IF (SING .AND. IPRINT.NE.0) WRITE (NOUT,150)
  150 FORMAT (21H JACOBIAN IS SINGULAR)
      REOLD = RABS
      SCPR = 1.E-10
      IF (CASI .OR. LIN) GO TO 175
      TMIN = TPM8
      GNOR = 0.
      T = 1.
      STEP = .FALSE.
      SCPR = 0.
      PNOR = 0.
      DO 155 I=1,MPN
         PNOR = PNOR+UU(I)**2
         GNOR = GNOR+GRADF(I)**2
         SCPR = SCPR+GRADF(I)*UU(I)
  155 CONTINUE
      IF (PNOR.EQ.0.) GO TO 115
C                                  WE CHECK IF THE DIRECTION UU IS OF
C                                    DESCENT (AS IT SHOULD) AND ALSO IF
C                                    THE IDENTITY (GRADF,UU)=RABS IS
C                                    APPROXIMATELY VERIFIED. IF EITHER
C                                    ONE OF THESE CHECKS FAIL, WE TAKE
C                                    GRADF AS OUR NEXT SEARCH
C                                    DIRECTION, SINCE THE ABOVE
C                                    INDICATES THAT UU IS AN UNRELIABLE
C                                    DIRECTION DUE TO ILL-CONDITIONING
C                                    IN THE JACOBIAN.
      IF (SCPR) 165, 165, 160
  160 IF (ABS(SCPR-RABS).LE..1*RABS) GO TO 170
  165 SCPR = -1.
  170 TE = 1.
C                                  APPROXIMATE SOLUTION IS CORRECTED
  175 DO 180 I=1,MPN
         IF (SCPR.LE.0.) UU(I) = GRADF(I)
         Y(I) = Y(I)+UU(I)
  180 CONTINUE
      IF (SCPR.LE.0.) SCPR = GNOR
      T = 1.
      INWT = INWT+1
      GO TO 20
C
      END
