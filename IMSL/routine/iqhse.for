C   IMSL ROUTINE NAME   - IQHSE
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JULY 1, 1983
C
C   PURPOSE             - NUCLEUS CALLED ONLY BY IMSL SUBROUTINE IQHSCV
C
C   REQD. IMSL ROUTINES - SINGLE/LLSQF,UERTST,UGETIO,VBLA=SASUM,
C                           VBLA=SDOT,VBLA=SNRM2,VHS12
C                         DOUBLE/LLSQF,UERTST,UGETIO,VBLA=DASUM,
C                           VBLA=DDOT,VBLA=DNRM2,VHS12
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE IQHSE (NDP,XD,YD,ZD,PD)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            NDP
      REAL               XD(1),YD(1),ZD(1),PD(5,1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            IC,IM1,I,J,K,KBASIS,KER,NA,NC,NPTS,IWK(5)
      REAL               BB(5),B(12),R,D,C(12,5),RMD,TOL,REPS,EPS,
     *                   WKAREA(5),XDD,YDD,CLOSE(12),BIG,TEMP
      DATA               REPS/Z3C100000/
      DATA               NC,IC /5,12/
C                                  FIRST EXECUTABLE STATEMENT
      BIG = 0.0
      DO 10 I=2,NDP
         IM1 = I-1
         DO 5 J=1,IM1
            BIG = AMAX1(BIG,(XD(I)-XD(J))**2+(YD(I)-YD(J))**2)
    5    CONTINUE
   10 CONTINUE
      EPS = 100.*SQRT(BIG)*REPS
      DO 45 I=1,NDP
C                                  FIND RADIUS WHICH INCLUDES
C                                    10 NEIGHBORS
         DO 15 K=1,IC
            CLOSE(K) = BIG
   15    CONTINUE
         DO 25 J=1,NDP
            IF (I.EQ.J) GO TO 25
            XDD = XD(J)-XD(I)
            YDD = YD(J)-YD(I)
            D = XDD**2+YDD**2
            IF (D.GE.CLOSE(IC)) GO TO 25
            CLOSE(IC) = D
            DO 20 K=1,IC
               IF (CLOSE(K).LE.CLOSE(IC)) GO TO 20
               TEMP = CLOSE(K)
               CLOSE(K) = CLOSE(IC)
               CLOSE(IC) = TEMP
   20       CONTINUE
   25    CONTINUE
         R = SQRT(CLOSE(IC))
         NPTS = 0
         DO 30 J=1,NDP
            IF (I.EQ.J) GO TO 30
            XDD = XD(J)-XD(I)
            YDD = YD(J)-YD(I)
            D = XDD**2+YDD**2
            IF (D.GT.CLOSE(IC)) GO TO 30
            D = SQRT(D)
            RMD = R/D
            NPTS = NPTS+1
            C(NPTS,1) = RMD*XDD
            C(NPTS,2) = RMD*YDD
            C(NPTS,3) = C(NPTS,1)*XDD/2.
            C(NPTS,4) = C(NPTS,1)*YDD
            C(NPTS,5) = C(NPTS,2)*YDD/2.
            B(NPTS) = (ZD(J)-ZD(I))*RMD
            IF (NPTS.GE.IC) GO TO 35
   30    CONTINUE
   35    NA = NC
         IF (NPTS.LT.5) NA = 2
         KBASIS = 0
         TOL = EPS/R
C                                  LEAST-SQUARES FIT TO NEIGHBORS
C                                    USING QUADRATIC POLYNOMIAL
         CALL LLSQF(C,IC,NPTS,NA,B,TOL,KBASIS,BB,WKAREA,IWK,KER)
         IF (KER.GT.0) NA = 0
         DO 40 J=1,5
            PD(J,I) = 0.0
            IF (J.LE.NA) PD(J,I) = BB(J)
   40    CONTINUE
   45 CONTINUE
      RETURN
      END
