C   IMSL ROUTINE NAME   - DGRCS
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - NUCLEUS CALLED ONLY BY IMSL SUBROUTINE DGEAR
C
C   REQD. IMSL ROUTINES - NONE REQUIRED
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE DGRCS  (METH,NQ,EL,TQ,MAXDER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            METH,NQ,MAXDER
      REAL               TQ(1)
      REAL               EL(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            K
      REAL               PERTST(12,2,3)
      DATA               PERTST/1.,1.,2.,1.,.3158,.7407E-1,
     1                   .1391E-1,.2182E-2,.2945E-3,.3492E-4,
     2                   .3692E-5,.3524E-6,1.,1.,.5,.1667,
     3                   .4167E-1,7*1.,2.,12.,24.,37.89,
     4                   53.33,70.08,87.97,106.9,126.7,
     5                   147.4,168.8,191.0,2.0,4.5,7.333,
     6                   10.42,13.7,7*1.,12.0,24.0,37.89,
     7                   53.33,70.08,87.97,106.9,126.7,
     8                   147.4,168.8,191.0,1.,3.0,6.0,
     9                   9.167,12.5,8*1./
C                                  FIRST EXECUTABLE STATEMENT
      GO TO (5,10), METH
    5 MAXDER = 12
      GO TO (15,20,25,30,35,40,45,50,55,60,65,70), NQ
   10 MAXDER = 5
      GO TO (75,80,85,90,95), NQ
C                                  THE FOLLOWING COEFFICIENTS SHOULD BE
C                                    DEFINED TO MACHINE ACCURACY. FOR A
C                                    GIVEN ORDER NQ, THEY CAN BE
C                                    CALCULATED BY USE OF THE
C                                    GENERATING POLYNOMIAL L(T), WHOSE
C                                    COEFFICIENTS ARE EL(I).. L(T) =
C                                    EL(1) + EL(2)*T + ... +
C                                    EL(NQ+1)*T**NQ. FOR THE IMPLICIT
C                                    ADAMS METHODS, L(T) IS GIVEN BY
C                                    DL/DT = (T+1)*(T+2)* ...
C                                    *(T+NQ-1)/K, L(-1) = 0, WHERE K =
C                                    FACTORIAL(NQ-1). FOR THE GEAR
C                                    METHODS, L(T) = (T+1)*(T+2)* ...
C                                    *(T+NQ)/K, WHERE K =
C                                    FACTORIAL(NQ)*(1 + 1/2 + ... +
C                                    1/NQ). THE ORDER IN WHICH THE
C                                    GROUPS APPEAR BELOW IS.. IMPLICIT
C                                    ADAMS METHODS OF ORDERS 1 TO 12,
C                                    BACKWARD DIFFERENTIATION METHODS
C                                    OF ORDERS 1 TO 5.
   15 EL(1) = 1.0
      GO TO 100
   20 EL(1) = 0.5
      EL(3) = 0.5
      GO TO 100
   25 EL(1) = 4.166667E-01
      EL(3) = 0.75
      EL(4) = 1.666667E-01
      GO TO 100
   30 EL(1) = 0.375
      EL(3) = 9.166667E-01
      EL(4) = 3.333333E-01
      EL(5) = 4.166667E-02
      GO TO 100
   35 EL(1) = 3.486111E-01
      EL(3) = 1.041667E0
      EL(4) = 4.861111E-01
      EL(5) = 1.041667E-01
      EL(6) = 8.333333E-03
      GO TO 100
   40 EL(1) = 3.298611E-01
      EL(3) = 1.141667E+00
      EL(4) = 0.625E+00
      EL(5) = 1.770833E-01
      EL(6) = 0.025E+00
      EL(7) = 1.388889E-03
      GO TO 100
   45 EL(1) = 3.155919E-01
      EL(3) = 1.225E+00
      EL(4) = 7.518519E-01
      EL(5) = 2.552083E-01
      EL(6) = 4.861111E-02
      EL(7) = 4.861111E-03
      EL(8) = 1.984127E-04
      GO TO 100
   50 EL(1) = 3.042245E-01
      EL(3) = 1.296429E+00
      EL(4) = 8.685185E-01
      EL(5) = 3.357639E-01
      EL(6) = 7.777778E-02
      EL(7) = 1.064815E-02
      EL(8) = 7.936508E-04
      EL(9) = 2.480159E-05
      GO TO 100
   55 EL(1) = 2.948680E-01
      EL(3) = 1.358929E+00
      EL(4) = 9.765542E-01
      EL(5) = 4.171875E-01
      EL(6) = 1.113542E-01
      EL(7) = 0.01875E+00
      EL(8) = 1.934524E-03
      EL(9) = 1.116071E-04
      EL(10)= 2.755732E-06
      GO TO 100
   60 EL(1) = 2.869754E-01
      EL(3) = 1.414484E+00
      EL(4) = 1.077216E+00
      EL(5) = 4.985670E-01
      EL(6) = 1.484375E-01
      EL(7) = 2.906057E-02
      EL(8) = 3.720238E-03
      EL(9) = 2.996858E-04
      EL(10)= 1.377866E-05
      EL(11)= 2.755732E-07
      GO TO 100
   65 EL(1) = 2.801896E-01
      EL(3) = 1.464484E+00
      EL(4) = 1.171515E+00
      EL(5) = 5.793582E-01
      EL(6) = 1.883229E-01
      EL(7) = 4.143036E-02
      EL(8) = 6.211144E-03
      EL(9) = 6.252067E-04
      EL(10)= 4.041740E-05
      EL(11)= 1.515653E-06
      EL(12)= 2.505211E-08
      GO TO 100
   70 EL(1) = 2.742655E-01
      EL(3) = 1.509939E+00
      EL(4) = 1.260271E+00
      EL(5) = 6.592342E-01
      EL(6) = 2.304580E-01
      EL(7) = 5.569725E-02
      EL(8) = 9.439484E-03
      EL(9) = 1.119275E-03
      EL(10) = 9.093915E-05
      EL(11) = 4.822531E-06
      EL(12)= 1.503127E-07
      EL(13)= 2.087676E-09
      GO TO 100
C
   75 EL(1) = 1.0E+00
      GO TO 100
   80 EL(1) = 6.666667E-01
      EL(3) = 3.333333E-01
      GO TO 100
   85 EL(1) = 5.454545E-01
      EL(3) = EL(1)
      EL(4) = 9.090909E-02
      GO TO 100
   90 EL(1) = 0.48E+00
      EL(3) = 0.7E+00
      EL(4) = 0.2E+00
      EL(5) = 0.02E+00
      GO TO 100
   95 EL(1) = 4.379562E-01
      EL(3) = 8.211679E-01
      EL(4) = 3.102190E-01
      EL(5) = 5.474453E-02
      EL(6) = 3.649635E-03
C
  100 DO 105 K=1,3
         TQ(K) = PERTST(NQ,METH,K)
  105 CONTINUE
      TQ(4) = .5*TQ(2)/(NQ+2)
      RETURN
      END
