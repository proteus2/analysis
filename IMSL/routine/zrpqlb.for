C   IMSL ROUTINE NAME   - ZRPQLB
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - NUCLEUS CALLED ONLY BY IMSL SUBROUTINE
C                           ZRPOLY
C
C   REQD. IMSL ROUTINES - ZRPQLC,ZRPQLD,ZRPQLE,ZRPQLF,ZRPQLG,ZRPQLH,
C                           ZRPQLI
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE ZRPQLB (L2,NZ)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            L2,NZ
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            N,NN,J,ITYPE,I,IFLAG
      REAL               ARE,BETAS,BETAV,ETA,OSS,OTS,OTV,OVV,RMRE,SS,
     1                   TS,TSS,TV,TVV,VV
      REAL               P(101),QP(101),RK(101),QK(101),SVK(101)
      REAL               SR,SI,U,V,RA,RB,C,D,A1,A2,A3,
     1                   A6,A7,E,F,G,H,SZR,SZI,RLZR,RLZI,
     2                   SVU,SVV,UI,VI,S,ZERO
      LOGICAL            VPASS,SPASS,VTRY,STRY
      COMMON /ZRPQLJ/    P,QP,RK,QK,SVK,SR,SI,U,V,RA,RB,C,D,A1,A2,A3,A6,
     1                   A7,E,F,G,H,SZR,SZI,RLZR,RLZI,ETA,ARE,RMRE,N,NN
      DATA               ZERO/0.0/
C                                  FIRST EXECUTABLE STATEMENT
      NZ = 0
C                                  COMPUTES UP TO L2 FIXED SHIFT
C                                    K-POLYNOMIALS, TESTING FOR
C                                    CONVERGENCE IN THE LINEAR OR
C                                    QUADRATIC CASE. INITIATES ONE OF
C                                    THE VARIABLE SHIFT ITERATIONS AND
C                                    RETURNS WITH THE NUMBER OF ZEROS
C                                    FOUND.
C                                  L2 - LIMIT OF FIXED SHIFT STEPS
C                                  NZ -NUMBER OF ZEROS FOUND
      BETAV = .25
      BETAS = .25
      OSS = SR
      OVV = V
C                                  EVALUATE POLYNOMIAL BY SYNTHETIC
C                                    DIVISION
      CALL ZRPQLH (NN,U,V,P,QP,RA,RB)
      CALL ZRPQLE (ITYPE)
      DO 40 J=1,L2
C                                  CALCULATE NEXT K POLYNOMIAL AND
C                                    ESTIMATE V
         CALL ZRPQLF (ITYPE)
         CALL ZRPQLE (ITYPE)
         CALL ZRPQLG (ITYPE,UI,VI)
         VV = VI
C                                  ESTIMATE S
         SS = 0.
         IF (RK(N).NE.ZERO) SS = -P(NN)/RK(N)
         TV = 1.
         TS = 1.
         IF (J.EQ.1.OR.ITYPE.EQ.3) GO TO 35
C                                  COMPUTE RELATIVE MEASURES OF
C                                    CONVERGENCE OF S AND V SEQUENCES
         IF (VV.NE.0.) TV = ABS((VV-OVV)/VV)
         IF (SS.NE.0.) TS = ABS((SS-OSS)/SS)
C                                  IF DECREASING, MULTIPLY TWO MOST
C                                    RECENT CONVERGENCE MEASURES
         TVV = 1.
         IF (TV.LT.OTV) TVV = TV*OTV
         TSS = 1.
         IF (TS.LT.OTS) TSS = TS*OTS
C                                  COMPARE WITH CONVERGENCE CRITERIA
         VPASS = TVV.LT.BETAV
         SPASS = TSS.LT.BETAS
         IF (.NOT.(SPASS.OR.VPASS)) GO TO 35
C                                  AT LEAST ONE SEQUENCE HAS PASSED THE
C                                    CONVERGENCE TEST. STORE VARIABLES
C                                    BEFORE ITERATING
         SVU = U
         SVV = V
         DO 5 I=1,N
    5    SVK(I) = RK(I)
         S = SS
C                                  CHOOSE ITERATION ACCORDING TO THE
C                                    FASTEST CONVERGING SEQUENCE
         VTRY = .FALSE.
         STRY = .FALSE.
         IF (SPASS.AND.((.NOT.VPASS).OR.TSS.LT.TVV)) GO TO 20
   10    CALL ZRPQLC (UI,VI,NZ)
         IF (NZ.GT.0) RETURN
C                                  QUADRATIC ITERATION HAS FAILED. FLAG
C                                    THAT IT HAS BEEN TRIED AND
C                                    DECREASE THE CONVERGENCE
C                                    CRITERION.
         VTRY = .TRUE.
         BETAV = BETAV*.25
C                                  TRY LINEAR ITERATION IF IT HAS NOT
C                                    BEEN TRIED AND THE S SEQUENCE IS
C                                    CONVERGING
         IF (STRY.OR.(.NOT.SPASS)) GO TO 25
         DO 15 I=1,N
   15    RK(I) = SVK(I)
   20    CALL ZRPQLD (S,NZ,IFLAG)
         IF (NZ.GT.0) RETURN
C                                  LINEAR ITERATION HAS FAILED. FLAG
C                                    THAT IT HAS BEEN TRIED AND
C                                    DECREASE THE CONVERGENCE CRITERION
         STRY = .TRUE.
         BETAS = BETAS*.25
         IF (IFLAG.EQ.0) GO TO 25
C                                  IF LINEAR ITERATION SIGNALS AN
C                                    ALMOST DOUBLE REAL ZERO ATTEMPT
C                                    QUADRATIC INTERATION
         UI = -(S+S)
         VI = S*S
         GO TO 10
C                                  RESTORE VARIABLES
   25    U = SVU
         V = SVV
         DO 30 I=1,N
   30    RK(I) = SVK(I)
C                                  TRY QUADRATIC ITERATION IF IT HAS
C                                    NOT BEEN TRIED AND THE V SEQUENCE
C                                    IS CONVERGING
         IF (VPASS.AND.(.NOT.VTRY)) GO TO 10
C                                  RECOMPUTE QP AND SCALAR VALUES TO
C                                    CONTINUE THE SECOND STAGE
         CALL ZRPQLH (NN,U,V,P,QP,RA,RB)
         CALL ZRPQLE (ITYPE)
   35    OVV = VV
         OSS = SS
         OTV = TV
         OTS = TS
   40 CONTINUE
      RETURN
      END
