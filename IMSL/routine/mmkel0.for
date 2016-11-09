C   IMSL ROUTINE NAME   - MMKEL0
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - KELVIN FUNCTIONS OF THE FIRST KIND, (BER,BEI),
C                           AND OF THE SECOND KIND, (KER,KEI), OF ORDER
C                           ZERO
C
C   USAGE               - CALL MMKEL0 (X,BER,BEI,XKER,XKEI,IER)
C
C   ARGUMENTS    X      - INPUT ARGUMENT. IF X IS NEGATIVE, A WARNING
C                           ERROR IS PRODUCED AND VALUES OF POSITIVE
C                           MACHINE INFINITY WILL BE RETURNED FOR XKER
C                           AND XKEI.
C                BER    - OUTPUT VALUE OF BER X . (ORDER ZERO)
C                BEI    - OUTPUT VALUE OF BEI X . (ORDER ZERO)
C                XKER   - OUTPUT VALUE OF KER X . (ORDER ZERO)
C                XKEI   - OUTPUT VALUE OF KEI X . (ORDER ZERO)
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER = 129 INDICATES THAT THE ABSOLUTE VALUE
C                             OF X WAS GREATER THAN 119. BER AND BEI ARE
C                             SET TO ZERO. IF X IS NON-NEGATIVE, XKER
C                             AND XKEI ARE ALSO SET TO ZERO. OTHERWISE,
C                             XKER AND XKEI ARE SET TO POSITIVE MACHINE
C                             INFINITY.
C                         WARNING ERROR
C                           IER = 34 INDICATES THAT X IS NEGATIVE.
C                             XKER AND XKEI WILL BE RETURNED AS
C                             POSITIVE MACHINE INFINITY.
C
C   REQD. IMSL ROUTINES - UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE MMKEL0 (X,BER,BEI,XKER,XKEI,IER)
C
      DIMENSION          C1(9),C2(9),C3(9),C4(9),E1(9),E2(9)
      DOUBLE PRECISION   C1,C2,C3,C4,E1,E2,PIO8,RT2,XINF,
     *                   PI,EUL,TEN,ZERO,HALF,ONE,ARG,BER,BEI,B1,B2,B3,
     *                   B4,CON,DC,DCM,DE,DS,DSM,DSQ,PIO2,R1,R2,S,SM,T,
     *                   TM,TWOPI,X,XKER,XKEI,Z,ZI,ZIM,ZSQ,Z4,ZMAX
C
C                                  COEFFICIENTS FOR EVALUATION OF
C                                  BER-SUB-ZERO(X) FOR X GREATER THAN
C                                  0. AND LESS THAN OR EQUAL TO 10.
C
      DATA               C1(1)/5.16070465D-5/,C1(2)/-4.8987125727D-3/
      DATA               C1(3)/.25977730007D0/,C1(4)/-7.2422567278207D0/
      DATA               C1(5)/93.8596692971726D0/
      DATA               C1(6)/-470.9502795889968D0/
      DATA               C1(7)/678.1684027663091D0/
      DATA               C1(8)/-156.2499999995701D0/
      DATA               C1(9)/.9999999999974D0/
C
C                                  COEFFICIENTS FOR EVALUATION OF
C                                  BEI-SUB-ZERO(X) FOR X GREATER THAN 0.
C                                  AND X LESS THAN OR EQUAL TO 10.
C
      DATA               C2(1)/4.4913000D-6/,C2(2)/-5.444243175D-4/
      DATA               C2(3)/3.84288282734D-2/
      DATA               C2(4)/-1.4963342749742D0/
      DATA               C2(5)/28.9690338786499D0/
      DATA               C2(6)/-240.2807549442574D0/
      DATA               C2(7)/678.1684027769807D0/
      DATA               C2(8)/-434.0277777777479D0/
      DATA               C2(9)/24.9999999999998D0/
C
C                                  COEFFICIENTS FOR EVALUATION OF
C                                  KEI-SUB-ZERO(X) FOR X GREATER THAN
C                                  OR EQUAL TO ZERO AND X LESS THAN OR
C                                  EQUAL TO 10.
C
      DATA               C3(1)/1.54363047D-5/,C3(2)/-1.8064777860D-3/
      DATA               C3(3)/.1222087382192D0/
      DATA               C3(4)/-4.5187459132639D0/
      DATA               C3(5)/81.9524771606200D0/
      DATA               C3(6)/-623.0136717405201D0/
      DATA               C3(7)/1548.484519673099D0/
      DATA               C3(8)/-795.7175925924866D0/
      DATA               C3(9)/24.9999999999993D0/
C
C                                  COEFFICIENTS FOR EVALUATION OF
C                                  KER-SUB-ZERO(X) FOR X GREATER THAN OR
C                                  EQUAL TO ZERO AND X LESS THAN OR
C                                  EQUAL TO TEN
C
      DATA               C4(1)/1.2161109D-6/,C4(2)/-1.797627986D-4/
      DATA               C4(3)/1.59380149705D-2/
      DATA               C4(4)/-.8061529027876D0/
      DATA               C4(5)/21.2123451660231D0/
      DATA               C4(6)/-255.0971742710479D0/
      DATA               C4(7)/1153.828185281456D0/
      DATA               C4(8)/-1412.850839120364D0/
      DATA               C4(9)/234.375D0/
C
C                                  COEFFICIENTS FOR EVALUATION OF
C                                  AUXILIARY FUNCTIONS FOR X GREATER
C                                  THAN 10.
C
      DATA               E1(1)/4.92D-8/,E1(2)/1.452D-7/,E1(3)/1.35D-8/
      DATA               E1(4)/-1.6192D-6/,E1(5)/-1.12207D-5/
      DATA               E1(6)/-5.17869D-5/,E1(7)/7.0D-10/
      DATA               E1(8)/8.8388346D-3/,E1(9)/1.0D0/
      DATA               E2(1)/-2.43D-8/,E2(2)/7.5D-8/,E2(3)/5.929D-7/
      DATA               E2(4)/1.6431D-6/,E2(5)/-7.2D-9/
      DATA               E2(6)/-5.18006D-5/,E2(7)/-7.031241D-4/
      DATA               E2(8)/-8.8388340D-3/,E2(9)/0.0D0/
C
C                                  MISCELLANEOUS CONSTANTS
C
      DATA               PIO2/1.570796326794897D0/
      DATA               TWOPI/6.283185307179586D0/
      DATA               PIO8/.3926990816987242D0/
      DATA               RT2/.7071067811865475D0/
      DATA               XINF/0.7237005577332260D+76/
      DATA               PI/3.141592653589793D0/
      DATA               EUL/.5772156649015329D0/
      DATA               TEN/10.D0/,ZERO/0.D0/,HALF/.5D0/,ONE/1.D0/
      DATA               ZMAX/119.D0/
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
      Z = DABS(X)
      IF (Z .GT. TEN) GO TO 15
      IF (Z .EQ. ZERO) GO TO 10
C                                  CALCULATION OF FUNCTIONS FOR ABS(X)
C                                  LESS THAN 10.
      Z = Z/TEN
      ZSQ = Z*Z
      Z4 = ZSQ*ZSQ
      B1 = C1(1)
      B2 = C2(1)
      B3 = C3(1)
      B4 = C4(1)
      DO 5 I = 2,9
         B1 = B1*Z4+C1(I)
         B2 = B2*Z4+C2(I)
         B3 = B3*Z4+C3(I)
         B4 = B4*Z4+C4(I)
    5 CONTINUE
      BER = B1
      BEI = ZSQ*B2
      IF (X .LT. ZERO) GO TO 30
      R1 = ZSQ*B3
      R2 = Z4*B4
      CON = (DLOG(X*HALF)+EUL)
      XKEI = -PIO2*HALF*BER+(R1-BEI*CON)
      XKER = PIO2*HALF*BEI-(R2+BER*CON)
      GO TO 9005
C                                  X EQUAL 0. DEFAULT TO PROPER VALUES
   10 BER = ONE
      BEI = ZERO
      XKEI = -HALF*PIO2
      XKER = XINF
      GO TO 9005
C                                  X GREATER THAN 10. CALCULATE
C                                  AUXILIARY FUNCTIONS
   15 IF (Z .GT. ZMAX) GO TO 25
      ZI = TEN/Z
      ZIM = -ZI
      S = E1(1)
      SM = S
      T = E2(1)
      TM = T
      DO 20 I = 2,9
         S = S*ZI+E1(I)
         T = T*ZI+E2(I)
         SM = SM*ZIM+E1(I)
         TM = TM*ZIM+E2(I)
   20 CONTINUE
      ARG = Z*RT2
      DS = DSIN(ARG-PIO8)
      DC = DCOS(ARG-PIO8)
      DSM = DSIN(ARG+PIO8)
      DCM = DCOS(ARG+PIO8)
      DE = DEXP(ARG)
      DSQ = DSQRT(TWOPI*Z)
C                                  CALCULATE THE DESIRED FUNCTIONS
      BER = DE*(S*DC-T*DS)/DSQ
      BEI = DE*(T*DC+S*DS)/DSQ
      IF (X .LT. ZERO) GO TO 30
      XKEI = PI*(TM*DCM-SM*DSM)/(DE*DSQ)
      XKER = PI*(SM*DCM+TM*DSM)/(DE*DSQ)
      GO TO 9005
C                                  Z TOO LARGE.
   25 BER = ZERO
      BEI = ZERO
      IER = 129
      IF (X .LT. ZERO) GO TO 35
      XKEI = ZERO
      XKER = ZERO
      GO TO 9000
C                                  X LESS THAN 0. DEFAULT TO PROPER
C                                  VALUES
   30 IER = 34
   35 XKEI = XINF
      XKER = XINF
      IF (IER .EQ. 0) GO TO 9005
 9000 CONTINUE
      CALL UERTST (IER,6HMMKEL0)
 9005 RETURN
      END
