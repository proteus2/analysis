C   IMSL ROUTINE NAME   - DGRPS
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - NOVEMBER 1, 1984
C
C   PURPOSE             - NUCLEUS CALLED ONLY BY IMSL SUBROUTINE DGEAR
C
C   REQD. IMSL ROUTINES - LUDATF,LEQT1B,UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE DGRPS (FCN,FCNJ,Y,N0,CON,MITER,YMAX,SAVE1,SAVE2,PW,
     *                   EQUIL,IPIV,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N0,MITER,IPIV(1),IER
      REAL               Y(N0,1),CON,YMAX(1),SAVE1(1),SAVE2(1),PW(1),
     *                   EQUIL(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
C
      INTEGER            NC,MFC,KFLAG,JSTART,NQUSED,NSTEP,NFE,NJE,NPW,
     *                   NSQ,I,J1,J,NERROR,NSAVE1,NSAVE2,NEQUIL,NY,
     *                   IDUMMY(23),NLIM,II,IJ,LIM1,LIM2,NB,NLC,NUC,NWK
      REAL               SDUMMY(4)
      REAL               T,H,HMIN,HMAX,EPSC,UROUND,EPSJ,HUSED,D,R0,YJ,R,
     *                   D1,D2,WA,DUMMY(40)
      COMMON /DBAND/     NLC,NUC
      COMMON /GEAR/      T,H,HMIN,HMAX,EPSC,UROUND,EPSJ,HUSED,DUMMY,
     *                   SDUMMY,NC,MFC,KFLAG,JSTART,NSQ,NQUSED,NSTEP,
     *                   NFE,NJE,NPW,NERROR,NSAVE1,NSAVE2,NEQUIL,NY,
     *                   IDUMMY
C                                  THIS ROUTINE IS CALLED BY DGRST TO
C                                    COMPUTE AND PROCESS THE MATRIX P =
C                                    I - H*EL(1)*J , WHERE J IS AN
C                                    APPROXIMATION TO THE JACOBIAN. J
C                                    IS COMPUTED, EITHER BY THE USER-
C                                    SUPPLIED ROUTINE FCNJ IF MITER =
C                                    1, OR BY FINITE DIFFERENCING IF
C                                    MITER = 2. J IS STORED IN PW AND
C                                    REPLACED BY P, USING CON =
C                                    -H*EL(1). THEN P IS SUBJECTED TO
C                                    LU DECOMPOSITION IN PREPARATION
C                                    FOR LATER SOLUTION OF LINEAR
C                                    SYSTEMS WITH P AS COEFFICIENT
C                                    MATRIX. IN ADDITION TO VARIABLES
C                                    DESCRIBED PREVIOUSLY,
C                                    COMMUNICATION WITH DGRPS USES THE
C                                    FOLLOWING EPSJ = DSQRT(UROUND),
C                                    USED IN THE NUMERICAL JACOBIAN
C                                    INCREMENTS.
C
C                                  FIRST EXECUTABLE STATEMENT
      IF (NLC.EQ.-1) GO TO 45
C                                  BANDED JACOBIAN CASE
      NB = NLC+NUC+1
      NWK = NB*N0+1
      IF (MITER.EQ.2) GO TO 15
C                                  MITER = 1
      NLIM = NB*N0
      DO 5 I=1,NLIM
         PW(I) = 0.0
    5 CONTINUE
      CALL FCNJ(NC,T,Y,PW)
      DO 10 I=1,NLIM
         PW(I) = PW(I)*CON
   10 CONTINUE
      GO TO 35
C                                  MITER = 2
   15 D = 0.0
      DO 20 I=1,NC
   20 D = D+SAVE2(I)**2
      R0 = ABS(H)*SQRT(D)*1.0E+03*UROUND
      DO 30 J=1,NC
         YJ = Y(J,1)
         R = EPSJ*YMAX(J)
         R = AMAX1(R,R0)
         Y(J,1) = Y(J,1)+R
         D = CON/R
         CALL FCN(NC,T,Y,SAVE1)
         LIM1 = MAX0(1,J-NUC)
         LIM2 = MIN0(N0,J+NLC)
         DO 25 I=LIM1,LIM2
            IJ = (J-I+NLC)*N0+I
            PW(IJ) = (SAVE1(I)-SAVE2(I))*D
   25    CONTINUE
         Y(J,1) = YJ
   30 CONTINUE
C                                  ADD IDENTITY MATRIX.
   35 DO 40 I=1,NC
         II = NLC*N0+I
         PW(II) = PW(II)+1.0
   40 CONTINUE
C                                  DO LU DECOMPOSITION ON P
C
      CALL LEQT1B(PW,NC,NLC,NUC,N0,EQUIL,1,N0,1,PW(NWK),IER)
      RETURN
C                                  FULL JACOBIAN CASE
   45 IF (MITER.EQ.2) GO TO 55
C                                  MITER = 1
      CALL FCNJ(NC,T,Y,PW)
      DO 50 I=1,NSQ
   50 PW(I) = PW(I)*CON
      GO TO 75
C                                  MITER = 2
   55 D = 0.0
      DO 60 I=1,NC
   60 D = D+SAVE2(I)**2
      R0 = ABS(H)*SQRT(D)*1.0E+03*UROUND
      J1 = 0
      DO 70 J=1,NC
         YJ = Y(J,1)
         R = EPSJ*YMAX(J)
         R = AMAX1(R,R0)
         Y(J,1) = Y(J,1)+R
         D = CON/R
         CALL FCN(NC,T,Y,SAVE1)
         DO 65 I=1,NC
   65    PW(I+J1) = (SAVE1(I)-SAVE2(I))*D
         Y(J,1) = YJ
         J1 = J1+N0
   70 CONTINUE
C                                  ADD IDENTITY MATRIX.
   75 J = 1
      DO 80 I=1,NC
         PW(J) = PW(J)+1.0
         J = J+(N0+1)
   80 CONTINUE
C                                  DO LU DECOMPOSITION ON P.
C
      CALL LUDATF(PW,PW,NC,N0,0,D1,D2,IPIV,EQUIL,WA,IER)
      RETURN
      END
