C   IMSL ROUTINE NAME   - MMDAS
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - DAWSON INTEGRAL
C
C   USAGE               - FUNCTION MMDAS (ARG,IER)
C
C   ARGUMENTS    MMDAS  - OUTPUT PARAMETER CONTAINING THE VALUE
C                           EXP(-ARG*ARG) * THE INTEGRAL (FROM 0 TO ARG)
C                           OF EXP(T*T) DT. MMDAS MUST BE TYPED APPRO-
C                           PRIATELY IN THE CALLING PROGRAM. (SEE THE
C                           PRECISION/HARDWARE SECTION.)
C                ARG    - INPUT ARGUMENT.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER = 129 INDICATES THAT UNDERFLOW WOULD
C                             HAVE OCCURRED. MMDAW IS SET TO ZERO.
C
C
C   REQD. IMSL ROUTINES - UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      DOUBLE PRECISION FUNCTION MMDAS (ARG,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            IER
      DOUBLE PRECISION   ARG
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            IEND,JEND,KEND,LEND,IENDP1,JENDP1,KENDP1,
     *                   LENDP1,I
      DOUBLE PRECISION   P1(10),Q1(10),P2(9),Q2(8),P3(10),Q3(9),P4(9),
     *                   Q4(8),FRAC,SUMP,SUMQ,W2,X,Y,ZCON,XLARGE,
     *                   XSMALL,ZERO,ONE,HALF,TWENT5,TWEL25,SIX25,XMAX
      DATA               ZERO/0.D0/,ONE/1.D0/,HALF/.5D0/,TWENT5/25.D0/,
     *                   TWEL25/12.25D0/,SIX25/6.25D0/
      DATA               XMAX/2.056880697D62/
      DATA               ZCON/Z4080000000000002/,
     *                   XLARGE/Z477FFBAC00000000/,
     *                   XSMALL/Z3A15798EE2308C3A/
C
C                                  COEFFICIENTS FOR R(9,9) APPROXIMATION
C                                  USED FOR ABS(ARG) LESS THAN 2.5
C
      DATA               P1(1)/Z40FFFFFFFFFFFFFC/,
     *                   P1(2)/ZC0238CDBE0EB32E0/,
     *                   P1(3)/Z3FC091A64C992A6E/,
     *                   P1(4)/ZBEBA6068CF41F366/,
     *                   P1(5)/Z3E1AAFC78FAC70DC/,
     *                   P1(6)/ZBCCFD51A51C04D47/,
     *                   P1(7)/Z3BF92DE05644969D/,
     *                   P1(8)/ZBA39EABA55EB3CF9/,
     *                   P1(9)/Z391CC398FF65D4B6/,
     *                   P1(10)/ZB72F539A0787AE4B/
      DATA               Q1(1)/Z4110000000000000/,
     *                   Q1(2)/Z40871DCEC9BF77CB/,
     *                   Q1(3)/Z4021D8B551AF9E36/,
     *                   Q1(4)/Z3F54F5DAB5A621CC/,
     *                   Q1(5)/Z3E9426C156134A4C/,
     *                   Q1(6)/Z3DBB9A0BEAAAA08D/,
     *                   Q1(7)/Z3CAE42CF1B91409C/,
     *                   Q1(8)/Z3B740A97E3A67939/,
     *                   Q1(9)/Z3A33398747CED1F8/,
     *                   Q1(10)/Z38BC4C9BB18B455D/
C
C                                  COEFFICIENTS FOR R(8,8) APPROXIMATION
C                                  IN J-FRACTION FORM, USED FOR 2.5
C                                  LESS THAN OR EQUAL TO ABS(ARG) LESS
C                                  THAN 3.5
C
      DATA               P2(1)/ZC11A9AD2BC45F509/,
     *                   P2(2)/ZC26BFF8D0B812B71/,
     *                   P2(3)/Z4260EC4F2724C319/,
     *                   P2(4)/Z414B41336DCA31B7/,
     *                   P2(5)/ZC1EA752D0A58DE37/,
     *                   P2(6)/Z415504B83C3A178A/,
     *                   P2(7)/Z416C2B41EF207D5B/,
     *                   P2(8)/Z4152628456C24702/,
     *                   P2(9)/Z407FEFDC294CFDAD/
      DATA               Q2(1)/Z40774476F5695B2F/,
     *                   Q2(2)/Z44293C87B6042CDA/,
     *                   Q2(3)/ZC1293A15BAF2EEB4/,
     *                   Q2(4)/Z42D178F9D0F274B8/,
     *                   Q2(5)/Z423927BA1B5C3478/,
     *                   Q2(6)/Z4311D1823BB82A4B/,
     *                   Q2(7)/ZC2ABD892165CB112/,
     *                   Q2(8)/Z404852419490BB14/
C
C                                  COEFFICIENTS FOR R(9,9) APPROXIMATION
C                                  IN J-FRACTION FORM, USED FOR ABS(ARG)
C                                  IN (3.5,5.0), INCLUSIVE.
C
      DATA               P3(1)/ZC148D3BE2BA2C06C/,
     *                   P3(2)/ZC212AA2A96672186/,
     *                   P3(3)/ZC175CF7D64052292/,
     *                   P3(4)/ZC242D739B0B60080/,
     *                   P3(5)/Z42307362CFFBF0DC/,
     *                   P3(6)/Z421AFAA396D7BF9D/,
     *                   P3(7)/ZC2218121571B8E90/,
     *                   P3(8)/Z41782781182AA26A/,
     *                   P3(9)/ZC117BFC9E99FFE2C/,
     *                   P3(10)/Z407FFFFCD3EDAF0B/
      DATA               Q3(1)/Z422CC8371A53A842/,
     *                   Q3(2)/Z4263DC582210F6D6/,
     *                   Q3(3)/Z41E061A33BE1B0E7/,
     *                   Q3(4)/Z43DA02D766BF95A9/,
     *                   Q3(5)/ZC19304F8D12A5519/,
     *                   Q3(6)/Z434D82F5C2A9FB2F/,
     *                   Q3(7)/ZC244CD705427E223/,
     *                   Q3(8)/ZC125001084E46DA6/,
     *                   Q3(9)/Z404002B8205F3193/
C
C                                  COEFFICIENTS FOR R(8,8) APPROXIMATION
C                                  IN J-FRACTION FORM, USED FOR ABS(ARG)
C                                  GREATER THAN 5.0
C
      DATA               P4(1)/ZC246C909F5E9317F/,
     *                   P4(2)/Z4217BD26686DF03A/,
     *                   P4(3)/ZC185DCC7C7CD5120/,
     *                   P4(4)/ZC21B5185BA4157FB/,
     *                   P4(5)/ZC15BE7F175370630/,
     *                   P4(6)/ZC166C905E0463FC5/,
     *                   P4(7)/ZC147FF94B5080098/,
     *                   P4(8)/ZC127FFFFFE5C6A8B/,
     *                   P4(9)/Z407FFFFFFFFFFE98/
      DATA               Q4(1)/Z4390FB8353EB91B5/,
     *                   Q4(2)/ZC17B94ABAB4ACC16/,
     *                   Q4(3)/Z42CADBB0B51C3CDC/,
     *                   Q4(4)/ZC21665D68C3B55FA/,
     *                   Q4(5)/ZC21000C5730A4065/,
     *                   Q4(6)/ZC1702668D9274223/,
     *                   Q4(7)/ZC1280001BEA98191/,
     *                   Q4(8)/Z40C000000008FAB1/
C                                  FIRST EXECUTABLE STATEMENT
      IEND = 8
      JEND = 9
      KEND = 8
      LEND = 9
      IENDP1 = IEND + 1
      JENDP1 = JEND + 1
      KENDP1 = KEND + 1
      LENDP1 = LEND + 1
      IER = 0
      X = ARG
      IF (DABS(X) .GT. XLARGE) GO TO 40
      IF (DABS(X) .LT. XSMALL) GO TO 45
      Y = X * X
      IF (Y .GE. SIX25) GO TO 10
C                                  ABS(X) IS LESS THAN 2.5
      SUMP = P1(LENDP1)*Y
      SUMQ = Q1(LENDP1)*Y
      DO 5 I=1,LEND
         II = LENDP1-I
         SUMP = (SUMP + P1(II))*Y
         SUMQ = (SUMQ + Q1(II))*Y
    5 CONTINUE
      MMDAS = X * SUMP / SUMQ
      GO TO 9005
   10 IF (Y .GE. TWEL25) GO TO 20
C                                  ABS(X) IS GREATER THAN OR EQUAL TO
C                                  2.5 AND LESS THAN 3.5
      FRAC = ZERO
      DO 15 I = 1,IEND
         FRAC = Q2(I) / (P2(I) + Y + FRAC)
   15 CONTINUE
      MMDAS = (P2(IENDP1) + FRAC) / X
      GO TO 9005
   20 IF (Y .GE. TWENT5) GO TO 30
C                                  ABS(X) IN GREATER THAN OR EQUAL TO
C                                  3.5 AND LESS THAN 5.
      FRAC = ZERO
      DO 25 I = 1,JEND
         FRAC = Q3(I) / (P3(I) + Y + FRAC)
   25 CONTINUE
      MMDAS = (P3(JENDP1) + FRAC) / X
      GO TO 9005
C                                  ABS(X) IN (5.0,XLARGE)
   30 W2 = ONE / X / X
      FRAC = ZERO
      DO 35 I = 1,KEND
         FRAC = Q4(I) / (P4(I) + Y + FRAC)
   35 CONTINUE
      FRAC = P4(KENDP1) + FRAC
      MMDAS = (ZCON + HALF * W2 * FRAC) / X
      GO TO 9005
C                                  ABS(X) IS GREATER THAN XLARGE
   40 IF (DABS(X) .GT. XMAX) GO TO 50
      MMDAS = HALF / X
      GO TO 9005
C                                  RETURN FOR SMALL X
   45 MMDAS = X
      GO TO 9005
   50 MMDAS = ZERO
      IER = 129
 9000 CONTINUE
      CALL UERTST(IER,6HMMDAS )
 9005 RETURN
      END
