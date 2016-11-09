C   IMSL ROUTINE NAME   - DPDEX
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1982
C
C   PURPOSE             - NUCLEUS CALLED ONLY BY IMSL SUBROUTINE DPDES
C
C   REQD. IMSL ROUTINES - DPDET,DPDEU,LEQT1B,UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE DPDEX (FCN,BNDRY,Y,N0,CON,MITER,YMAX,SAVE1,SAVE2,PW,
     *                   EQUIL,IPIV,XX,WORK,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N0,MITER,IPIV(1),IER
      REAL               Y(N0,1),CON,YMAX(1),SAVE1(1),SAVE2(1),PW(1),
     *                   EQUIL(1),XX(1),WORK(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
C
      INTEGER            NC,MFC,KFLAG,JSTART,NQUSED,NSTEP,NFE,NJE,NPW,
     *                   NSQ,I,J1,J,NERROR,NSAVE1,NSAVE2,NEQUIL,NY,IJ,
     *                   LIM1,LIM2,NB,NLIM,NWK,IDUMMY(23)
      REAL               SDUMMY(4)
      REAL               T,H,HMIN,HMAX,EPSC,UROUND,EPSJ,HUSED,D,R0,YJ,R,
     *                   D1,D2,WA,DUMMY(40)
      EXTERNAL           FCN,BNDRY
      COMMON /DBAND/     NLC,NUC
      INTEGER            NLC,NUC
      COMMON /GEAR/      T,H,HMIN,HMAX,EPSC,UROUND,EPSJ,HUSED,DUMMY,
     *                   SDUMMY,NC,MFC,KFLAG,JSTART,NSQ,NQUSED,NSTEP,
     *                   NFE,NJE,NPW,NERROR,NSAVE1,NSAVE2,NEQUIL,NY,
     *                   IDUMMY
C                                  THIS ROUTINE IS CALLED BY DGRST TO
C                                    COMPUTE AND PROCESS THE MATRIX P =
C                                    I - H*EL(1)*J , WHERE J IS AN
C                                    APPROXIMATION TO THE JACOBIAN. J
C                                    IS COMPUTED, EITHER BY THE USER-
C                                    SUPPLIED ROUTINE BNDRY IF MITER =
C                                    1, OR BY FINITE DIFFERENCING IF
C                                    MITER = 2. J IS STORED IN PW AND
C                                    REPLACED BY P, USING CON =
C                                    -H*EL(1). THEN P IS SUBJECTED TO
C                                    LU DECOMPOSITION IN PREPARATION
C                                    FOR LATER SOLUTION OF LINEAR
C                                    SYSTEMS WITH P AS COEFFICIENT
C                                    MATRIX. IN ADDITION TO VARIABLES
C                                    DESCRIBED PREVIOUSLY,
C                                    COMMUNICATION WITH DPDEX USES THE
C                                    FOLLOWING EPSJ = DSQRT(UROUND),
C                                    USED IN THE NUMERICAL JACOBIAN
C                                    INCREMENTS.
C
C                                  FIRST EXECUTABLE STATEMENT
C                                  BANDED JACOBIAN CASE
      NB = NLC+NUC+1
      NLIM = NB*N0
      NWK = NLIM+1
      DO 5 I=1,NLIM
    5 PW(I) = 0.0
C                                  MITER = 2
      DO 10 I=1,NC
   10 EQUIL(I) = SAVE2(I)
      CALL DPDET(NC,T,Y,SAVE1,XX,WORK,FCN,BNDRY,0)
      DO 15 I=1,NC
   15 SAVE2(I) = SAVE1(I)
      D = 0.0
      DO 20 I=1,NC
   20 D = D+SAVE2(I)**2
      R0 = ABS(H)*SQRT(D)*1.0E+03*UROUND
      DO 30 J=1,NC
         YJ = Y(J,1)
         R = EPSJ*YMAX(J)
         R = AMAX1(R,R0)
         Y(J,1) = Y(J,1)+R
         D = CON/R
         CALL DPDET(NC,T,Y,SAVE1,XX,WORK,FCN,BNDRY,J)
         LIM1 = MAX0(1,J-NUC)
         LIM2 = MIN0(N0,J+NLC)
         DO 25 I=LIM1,LIM2
            IJ = (J-I+NLC)*N0+I
            PW(IJ) = (SAVE1(I)-SAVE2(I))*D
   25    CONTINUE
         Y(J,1) = YJ
   30 CONTINUE
C                                  ADD A.
      DO 35 I=1,NLIM
         PW(I) = PW(I)+WORK(I)
   35 CONTINUE
      DO 40 I=1,NC
   40 SAVE2(I) = EQUIL(I)
C                                  DO LU DECOMPOSITION ON P
C
      CALL LEQT1B(PW,NC,NLC,NUC,N0,EQUIL,1,N0,1,PW(NWK),IER)
      RETURN
      END
