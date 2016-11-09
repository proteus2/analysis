C   IMSL ROUTINE NAME   - MMKELD
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1982
C
C   PURPOSE             - DERIVATIVES OF THE KELVIN FUNCTIONS (BER,BEI,
C                           KER, AND KEI) OF ORDER ZERO.
C
C   USAGE               - CALL MMKELD (X,BERP,BEIP,XKERP,XKEIP,IER)
C
C   ARGUMENTS    X      - INPUT ARGUMENT. IF X IS NEGATIVE, A WARNING
C                           ERROR IS PRODUCED AND VALUES OF POSITIVE
C                           MACHINE INFINITY WILL BE RETURNED FOR XKERP
C                           AND XKEIP.
C                BERP   - OUTPUT VALUE OF THE DERIVATIVE OF BER X .
C                           (ORDER ZERO)
C                BEIP   - OUTPUT VALUE OF THE DERIVATIVE OF BEI X .
C                           (ORDER ZERO)
C                XKERP  - OUTPUT VALUE OF THE DERIVATIVE OF KER X .
C                           (ORDER ZERO)
C                XKEIP  - OUTPUT VALUE OF THE DERIVATIVE OF KEI X .
C                           (ORDER ZERO)
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER = 129 INDICATES THAT THE ABSOLUTE VALUE
C                             OF X WAS GREATER THAN 119. BERP AND BEIP
C                             ARE SET TO ZERO. IF X IS NON-NEGATIVE,
C                             XKERP AND XKEIP ARE ALSO SET TO ZERO.
C                             OTHERWISE, XKERP AND XKEIP ARE SET TO
C                             POSITIVE MACHINE INFINITY.
C                         WARNING ERROR
C                           IER = 34 INDICATES THAT X IS NEGATIVE.
C                             XKERP AND XKEIP WILL BE RETURNED AS
C                             POSITIVE MACHINE INFINITY.
C
C   REQD. IMSL ROUTINES - MMKEL0,UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE MMKELD (X,BERP,BEIP,XKERP,XKEIP,IER)
C
      DIMENSION          D1(9),D2(9),D3(9),D4(9),E3(9),E4(9)
      DOUBLE PRECISION   ARG,BEI,BEIP,BER,BERP,B1,B2,B3,B4,CON,DC,DCM,
     *                   DE,DS,DSM,DSQ,D1,D2,D3,D4,EUL,E3,E4,PI,PIO2,
     *                   PIO8,RT2,R1P,R2P,TWOPI,U,UM,V,VM,X,XINF,XKEI,
     *                   XKEIP,XKER,XKERP,Z,ZI,ZIM,ZSQ,Z3,Z4,ZMAX
      DOUBLE PRECISION   TEN,ZERO,HALF
      DATA               TEN/10.D0/,ZERO/0.D0/,HALF/.5D0/
C
C                                  COEFFICIENTS FOR EVALUATION OF BERP-
C                                  SUB-ZERO(X) FOR X GREATER THAN 0. AND
C                                  LESS THAN OR EQUAL TO 10.
C
      DATA               D1(1)/-1.2506046D-6/,D1(2)/1.701453451D-4/
      DATA               D1(3)/-1.37246036190D-2/
      DATA               D1(4)/.6234726348243D0/
      DATA               D1(5)/-14.4845169498403D0/
      DATA               D1(6)/150.1754718432278D0/
      DATA               D1(7)/-565.1403356479486D0/
      DATA               D1(8)/542.5347222222147D0/
      DATA               D1(9)/-62.4999999999999D0/
C
C                                  COEFFICIENTS FOR EVALUATION OF BEIP-
C                                  SUB-ZERO(X) FOR X GREATER THAN 0.
C                                  AND LESS THAN OR EQUAL TO 10.
C
      DATA               D2(1)/1.52269884D-5/,D2(2)/-1.6331100837D-3/
      DATA               D2(3)/9.99147064932D-2/
      DATA               D2(4)/-3.2919352108579D0/
      DATA               D2(5)/52.1442608975905D0/
      DATA               D2(6)/-336.3930569023651D0/
      DATA               D2(7)/678.1684027747539D0/
      DATA               D2(8)/-260.4166666665533D0/
      DATA               D2(9)/4.9999999999993D0/
C
C                                  COEFFICIENTS FOR EVALUTION OF KEIP-
C                                  SUB-ZERO(X) FOR X GREATER THAN 0.
C                                  AND LESS THAN OR EQUAL TO 10.
C
      DATA               D3(1)/5.23294314D-5/
      DATA               D3(2)/-5.4188558408D-3/
      DATA               D3(3)/.3177418434686D0/
      DATA               D3(4)/-9.9412403209725D0/
      DATA               D3(5)/147.5144585913337D0/
      DATA               D3(6)/-872.2191403672455D0/
      DATA               D3(7)/1548.484519665204D0/
      DATA               D3(8)/-477.4305555551536D0/
      DATA               D3(9)/4.9999999999975D0/
C
C                                  COEFFICIENTS FOR EVALUATION OF KERP-
C                                  SUB-ZERO(X) FOR X GREATER THAN OR
C                                  EQUAL TO 0. AND LESS THAN OR EQUAL
C                                  TO 10.
C
      DATA               D4(1)/4.3682053D-6/,D4(2)/-5.752042283D-4/
      DATA               D4(3)/4.46263862145D-2/
      DATA               D4(4)/-1.9347669229237D0/
      DATA               D4(5)/42.4246903131088D0/
      DATA               D4(6)/-408.1554788292578D0/
      DATA               D4(7)/1384.593822337245D0/
      DATA               D4(8)/-1130.280671296269D0/
      DATA               D4(9)/93.7499999999998D0/
C
C
C                                  COEFFICIENTS FOR EVALUTION OF
C                                  AUXILIARY FUNCTIONS FOR X GREATER
C                                  THAN 10.
C
      DATA               E3(1)/-5.63D-8/,E3(2)/-1.671D-7/
      DATA               E3(3)/-1.47D-8/,E3(4)/1.9780D-6/
      DATA               E3(5)/1.44255D-5/,E3(6)/7.25024D-5/
      DATA               E3(7)/-8.0D-10/,E3(8)/-2.65165040D-2/
      DATA               E3(9)/1.0D0/
      DATA               E4(1)/-2.69D-8/,E4(2)/-8.83D-8/
      DATA               E4(3)/-6.992D-7/,E4(4)/-2.0042D-6/
      DATA               E4(5)/7.9D-9/,E4(6)/7.25179D-5/
      DATA               E4(7)/1.1718740D-3/,E4(8)/2.65165034D-2/
      DATA               E4(9)/0.0D0/
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
      DATA               ZMAX/119.D0/
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
      CALL MMKEL0(X,BER,BEI,XKER,XKEI,IER)
      Z = DABS(X)
      IF (Z .GT. TEN) GO TO 15
      IF (Z .EQ. ZERO) GO TO 10
C                                  CALCULATION OF FUNCTIONS FOR ABS(X)
C                                  LESS THAN 10.
      Z = Z/TEN
      ZSQ = Z*Z
      Z3 = ZSQ*Z
      Z4 = ZSQ*ZSQ
      B1 = D1(1)
      B2 = D2(1)
      B3 = D3(1)
      B4 = D4(1)
      DO 5 I = 2,9
         B1 = B1*Z4+D1(I)
         B2 = B2*Z4+D2(I)
         B3 = B3*Z4+D3(I)
         B4 = B4*Z4+D4(I)
    5 CONTINUE
      BERP = B1*Z3
      BEIP = Z*B2
      IF ( X .LT. ZERO) GO TO 30
      R1P = Z*B3
      R2P = Z3*B4
      CON = (DLOG(X*HALF) + EUL)
      V = DABS(X)
      XKEIP = -PIO2*HALF*BERP+(R1P-BEIP*CON-BEI/V)
      XKERP = PIO2*HALF*BEIP-(R2P+BERP*CON+BER/V)
      GO TO 9005
C                                  X EQUAL TO 0. DEFAULT TO PROPER
C                                  VALUES
   10 BERP = ZERO
      BEIP = ZERO
      XKEIP = ZERO
      XKERP = -XINF
      GO TO 9005
C                                  X GREATER THAN 10. CALCULATE
C                                  AUXILIARY FUNCTIONS
   15 IF (Z .GT. ZMAX) GO TO 25
      ZI = TEN/Z
      ZIM = -ZI
      U = E3(1)
      UM = U
      V = E4(1)
      VM = V
      DO 20 I = 2,9
         U = U*ZI+E3(I)
         V = V*ZI+E4(I)
         UM = UM*ZIM+E3(I)
         VM = VM*ZIM+E4(I)
   20 CONTINUE
      ARG = Z*RT2
      DS = DSIN(ARG-PIO8)
      DC = DCOS(ARG-PIO8)
      DSM = DSIN(ARG+PIO8)
      DCM = DCOS(ARG+PIO8)
      DE = DEXP(ARG)
      DSQ = DSQRT(TWOPI*Z)
C                                  CALCULATE THE DESIRED FUNCTIONS
      BERP = DE*(U*DCM-V*DSM)/DSQ
      BEIP = DE*(V*DCM+U*DSM)/DSQ
      IF (X .LT. ZERO) GO TO 30
      XKEIP = -PI*(VM*DC-UM*DS)/(DE*DSQ)
      XKERP = -PI*(UM*DC+VM*DS)/(DE*DSQ)
      GO TO 9005
C                                  Z TOO LARGE.
   25 BERP = ZERO
      BEIP = ZERO
      IER = 129
      IF (X .LT. ZERO) GO TO 35
      XKEIP = ZERO
      XKERP = ZERO
      GO TO 9000
C                                  X LESS THAN 0. DEFAULT TO PROPER
C                                  VALUES
   30 IER = 34
      BERP = -BERP
      BEIP = -BEIP
   35 XKERP = XINF
      XKEIP = XINF
      IF (IER .EQ. 0) GO TO 9005
 9000 CONTINUE
      CALL UERTST(IER,6HMMKELD)
 9005 RETURN
      END
