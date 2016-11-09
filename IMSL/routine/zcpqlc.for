C   IMSL ROUTINE NAME   - ZCPQLC
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - NUCLEUS CALLED ONLY BY IMSL SUBROUTINE
C                           ZCPOLY
C
C   REQD. IMSL ROUTINES - ZCPQLD,ZCPQLE,ZCPQLF,ZCPQLG,ZCPQLH,ZCPQLK,
C                           ZCPQLL,ZCPQLM
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE ZCPQLC (L2,ZR,ZI,CONV)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            L2
      REAL               ZR,ZI
      LOGICAL            CONV
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            NN,N,J,I
      REAL               PR(50),PI(50),HR(50),HI(50),QPR(50),QPI(50),
     1                   QHR(50),QHI(50),SHR(50),SHI(50),
     2                   SR,SI,TR,TI,PVR,PVI,ARE,RMRE,REPSR1,RINFP,
     3                   OTR,OTI,SVSR,SVSI,ZCPQLL,PT5
      LOGICAL            TEST,PASD,BOWL
      COMMON /ZCPQLN/    PR,PI,HR,HI,QPR,QPI,QHR,QHI,SHR,SHI,SR,SI,
     1                   TR,TI,PVR,PVI,ARE,RMRE,REPSR1,RINFP,NN
      DATA               PT5/0.5/
C                                  FIRST EXECUTABLE STATEMENT
      N = NN-1
C                                  COMPUTES L2 FIXED-SHIFT H
C                                    POLYNOMIALS AND TEST FOR
C                                    CONVERGENCE. INITIATES A
C                                    VARIABLE-SHIFT ITERATION AND
C                                    RETURN WITH THE APPROXIMATE ZERO
C                                    IF SUCCESSFUL.
C                                  L2 - LIMIT OF FIXED SHIFT STEPS
C                                  ZR,ZI - APPROXIMATE ZERO IF CONV IS
C                                    .TRUE. CONV - LOGICAL INDICATING
C                                    CONVERGENCE OF STAGE 3 ITERATION
C                                  EVALUATE P AT S
      CALL ZCPQLG (NN,SR,SI,PR,PI,QPR,QPI,PVR,PVI)
      TEST = .TRUE.
      PASD = .FALSE.
C                                  CALCULATE FIRST T = -P(S)/H(S)
      CALL ZCPQLE (BOWL)
C                                  MAIN LOOP FOR ONE SECOND STAGE STEP
      DO 25 J=1,L2
         OTR = TR
         OTI = TI
C                                  COMPUTE NEXT H POLYNOMIAL AND NEW T
         CALL ZCPQLF (BOWL)
         CALL ZCPQLE (BOWL)
         ZR = SR+TR
         ZI = SI+TI
C                                  TEST FOR CONVERGENCE UNLESS STAGE 3
C                                    HAS FAILED ONCE OR THIS IS THE
C                                    LAST H POLYNOMIAL
         IF (BOWL.OR..NOT.TEST.OR.J.EQ.L2) GO TO 25
         IF (ZCPQLL(TR-OTR,TI-OTI).GE.PT5*ZCPQLL(ZR,ZI)) GO TO 20
         IF (.NOT.PASD) GO TO 15
C                                  THE WEAK CONVERGENCE TEST HAS BEEN
C                                    PASSED TWICE, START THE THIRD
C                                    STAGE ITERATION, AFTER SAVING THE
C                                    CURRENT H POLYNOMIAL AND SHIFT.
         DO 5 I=1,N
            SHR(I) = HR(I)
            SHI(I) = HI(I)
    5    CONTINUE
         SVSR = SR
         SVSI = SI
         CALL ZCPQLD (10,ZR,ZI,CONV)
         IF (CONV) RETURN
C                                  THE ITERATION FAILED TO CONVERGE.
C                                    TURN OFF TESTING AND RESTORE
C                                    H,S,PV AND T.
         TEST = .FALSE.
         DO 10 I=1,N
            HR(I) = SHR(I)
            HI(I) = SHI(I)
   10    CONTINUE
         SR = SVSR
         SI = SVSI
         CALL ZCPQLG (NN,SR,SI,PR,PI,QPR,QPI,PVR,PVI)
         CALL ZCPQLE (BOWL)
         GO TO 25
   15    PASD = .TRUE.
         GO TO 25
   20    PASD = .FALSE.
   25 CONTINUE
C                                  ATTEMPT AN ITERATION WITH FINAL H
C                                    POLYNOMIAL FROM SECOND STAGE
      CALL ZCPQLD (10,ZR,ZI,CONV)
      RETURN
      END
