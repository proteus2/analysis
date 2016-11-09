C   IMSL ROUTINE NAME   - ZCPQLD
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - NUCLEUS CALLED ONLY BY IMSL SUBROUTINE
C                           ZCPOLY
C
C   REQD. IMSL ROUTINES - ZCPQLE,ZCPQLF,ZCPQLG,ZCPQLH,ZCPQLK,ZCPQLL,
C                           ZCPQLM
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE ZCPQLD (L3,ZR,ZI,CONV)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            L3
      REAL               ZR,ZI
      LOGICAL            CONV
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,NN,J
      REAL               PR(50),PI(50),HR(50),HI(50),QPR(50),QPI(50),
     1                   QHR(50),QHI(50),SHR(50),SHI(50),
     2                   SR,SI,TR,TI,PVR,PVI,ARE,RMRE,REPSR1,RINFP,
     3                   RMP,RMS,OMP,RELSTP,R1,R2,ZCPQLL,
     4                   ZCPQLH,TP,PT1,PTO5,ONE,TWENTY
      LOGICAL            B,BOWL
      COMMON /ZCPQLN/    PR,PI,HR,HI,QPR,QPI,QHR,QHI,SHR,SHI,SR,SI,
     1                   TR,TI,PVR,PVI,ARE,RMRE,REPSR1,RINFP,NN
      DATA               PT1,PTO5,ONE,TWENTY/0.1,0.05,1.0,20.0/
C                                  FIRST EXECUTABLE STATEMENT
      CONV = .FALSE.
      B = .FALSE.
      SR = ZR
      SI = ZI
C                                  CARRIES OUT THE THIRD STAGE
C                                    ITERATION.
C                                  L3 - LIMIT OF STEPS IN STAGE 3.
C                                  ZR,ZI - ON ENTRY CONTAINS THE
C                                    INITIAL ITERATE, IF THE ITERATION
C                                    CONVERGES IT CONTAINS THE FINAL
C                                    ITERATE ON EXIT
C                                  CONV - .TRUE. IF ITERATION CONVERGES
C                                  MAIN LOOP FOR STAGE THREE
      DO 30 I=1,L3
C                                  EVALUATE P AT S AND TEST FOR
C                                    CONVERGENCE
         CALL ZCPQLG (NN,SR,SI,PR,PI,QPR,QPI,PVR,PVI)
         RMP = ZCPQLL(PVR,PVI)
         RMS = ZCPQLL(SR,SI)
         IF (RMP.GT.TWENTY*ZCPQLH(NN,QPR,QPI,RMS,RMP,ARE,RMRE)) GO TO 5
C                                  POLYNOMIAL VALUE IS SMALLER IN VALUE
C                                    THAN A BOUND ON THE ERROR IN
C                                    EVALUATING P, TERMINATE THE
C                                    ITERATION
         CONV = .TRUE.
         ZR = SR
         ZI = SI
         RETURN
    5    IF (I.EQ.1) GO TO 20
         IF (B.OR.RMP.LT.OMP.OR.RELSTP.GE.PTO5) GO TO 15
C                                  ITERATION HAS STALLED. PROBABLY A
C                                    CLUSTER OF ZEROS. DO 5 FIXED SHIFT
C                                    STEPS INTO THE CLUSTER TO FORCE
C                                    ONE ZERO TO DOMINATE.
         TP = RELSTP
         B = .TRUE.
         IF (RELSTP.LT.REPSR1) TP = REPSR1
C1       R1 = DSQRT(TP)
         R1 = SQRT(TP)
         R2 = SR*(ONE+R1)-SI*R1
         SI = SR*R1+SI*(ONE+R1)
         SR = R2
         CALL ZCPQLG (NN,SR,SI,PR,PI,QPR,QPI,PVR,PVI)
         DO 10 J=1,5
            CALL ZCPQLE (BOWL)
            CALL ZCPQLF (BOWL)
   10    CONTINUE
         OMP = RINFP
         GO TO 25
C                                  EXIT IF POLYNOMIAL VALUE INCREASES
C                                    SIGNIFICANTLY
   15    IF (RMP*PT1.GT.OMP) RETURN
   20    OMP = RMP
C                                  CALCULATE NEXT ITERATE
   25    CALL ZCPQLE (BOWL)
         CALL ZCPQLF (BOWL)
         CALL ZCPQLE (BOWL)
         IF (BOWL) GO TO 30
         RELSTP = ZCPQLL(TR,TI)/ZCPQLL(SR,SI)
         SR = SR+TR
         SI = SI+TI
   30 CONTINUE
      RETURN
      END
