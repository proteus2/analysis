C   IMSL ROUTINE NAME   - MMBSYN
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1981
C
C   PURPOSE             - BESSEL FUNCTION OF THE SECOND KIND OF NON-
C                           NEGATIVE REAL FRACTIONAL ORDER FOR REAL
C                           POSITIVE ARGUMENTS
C
C   USAGE               - CALL MMBSYN (ARG,ORDER,N,YN,IER)
C
C   ARGUMENTS    ARG    - INPUT ARGUMENT. ARG MUST BE TYPED APPRO-
C                           PRIATELY IN THE CALLING PROGRAM. (SEE THE
C                           PRECISION/HARDWARE SECTION.) ARG MUST BE
C                           GREATER THAN OR EQUAL TO XMIN AND LESS THAN
C                           OR EQUAL TO XMAX. XMIN IS OF THE ORDER OF
C                           10**(-38) AND XMAX IS OF THE ORDER OF 10**8.
C                           THE EXACT VALUES OF XMIN AND XMAX MAY ALLOW
C                           LARGER RANGES FOR ARG ON SOME COMPUTERS. SEE
C                           THE PROGRAMMING NOTES IN THE MANUAL FOR THE
C                           EXACT VALUES.
C                ORDER  - INPUT VALUE SPECIFYING THE DESIRED ORDER OF
C                           THE BESSEL FUNCTION. ORDER MUST BE TYPED
C                           APPROPRIATELY IN THE CALLING PROGRAM. (SEE
C                           THE PRECISION/HARDWARE SECTION.) ORDER MUST
C                           BE GREATER THAN OR EQUAL TO ZERO AND LESS
C                           THAN ONE.
C                N      - INPUT INTEGER VALUE SPECIFYING THE NUMBER
C                           OF FUNCTION VALUES TO BE COMPUTED.
C                YN     - OUTPUT VECTOR OF LENGTH N CONTAINING THE
C                           COMPUTED FUNCTION VALUES. YN MUST BE TYPED
C                           APPROPRIATELY IN THE CALLING PROGRAM. (SEE
C                           THE PRECISION/HARDWARE SECTION.) YN(1) WILL
C                           CONTAIN THE COMPUTED VALUE FOR THE INPUT
C                           ORDER, YN(2) WILL CONTAIN THE COMPUTED
C                           FUNCTION VALUE FOR ORDER + 1, YN(3) FOR
C                           ORDER + 2, ETC.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR.
C                           IER = 129 INDICATES THAT ORDER IS OUT OF
C                             RANGE OR ARG IS LESS THAN XMIN.
C                             YN(I), (I= 1,N) IS SET TO NEGATIVE
C                             MACHINE INFINITY.
C                           IER = 130 INDICATES THAT N IS LESS THAN OR
C                             EQUAL TO ZERO. YN(1) IS SET TO NEGATIVE
C                             MACHINE INFINITY.
C                           IER = 131 INDICATES THAT OVERFLOW WOULD
C                             HAVE OCCURED IN VECTOR YN, AT ELEMENT J.
C                             YN(I), (I = J,N) IS SET TO NEGATIVE
C                             MACHINE INFINITY.
C                           IER = 132 INDICATES THAT ARG IS GREATER
C                             THAN XMAX.  YN(I), (I=1,N) IS SET TO
C                             ZERO.
C
C   REQD. IMSL ROUTINES - UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE MMBSYN (ARG,ORDER,N,YN,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N,IER
      DOUBLE PRECISION   ARG,ORDER,YN(N)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,J,K,L,M,NMAX,NN
      INTEGER            IEND,MEND,M1,M2,M3,IENDP2,MENDM3
      DOUBLE PRECISION   CO,CUM,CUMP,C2,C2P,DELX,DEN,DENP,
     1                   DLNHAF,DLNX,DLNXNU,FACT(36),FOOLU,F1(17),
     2                   F2(16),F3(17),GAM,GNU,GNU1,GNU2,GP(7),GQ(7),
     3                   P(11,6),PI2,Q(10,6),
     4                   RTPI2,RTXNU,S,SCRAT(6),SUM,SUMP,T,TWOPI,
     5                   TWOPI1,TWOPI2,T1,T2,T3(6),V,VCOS,VPI,VSIN,
     6                   W(17),X,XC,XIN,XINF,XK,XM,XNU,XNUM,XNUMP,
     7                   Z,ZH,ZMAX,ZMAX1,ZMIN,ZSTOP,CON1,CON2
      DOUBLE PRECISION   P1(11),P4(11),P7(11),P10(11),P13(11),P16(11),
     *                   Q1(10),Q4(10),Q7(10),Q10(10),Q13(10),Q16(10)
      LOGICAL            NFLAG
      EQUIVALENCE        (P(1,1),P1(1)),(P(1,2),P4(1)),(P(1,3),P7(1)),
     *                   (P(1,4),P10(1)),(P(1,5),P13(1)),
     *                   (P(1,6),P16(1)),(Q(1,1),Q1(1)),
     *                   (Q(1,2),Q4(1)),(Q(1,3),Q7(1)),
     *                   (Q(1,4),Q10(1)),(Q(1,5),Q13(1)),
     *                   (Q(1,6),Q16(1))
C
C                                  MACHINE DEPENDENT CONSTANTS
C                                     PI2 = 2 / PI
C                                     RTPI2 = SQRT(2/PI)
C                                     DLNHAF = LOG(.5)
C                                     TWOPI = 2 * PI
C                                     TWOPI1 + TWOPI2 = 2*PI TO EXTRA
C                                     PRECISION
C                                     XINF = LARGEST POSITIVE MACHINE
C                                     NUMBER
C                                     ZMAX = 16**9, LARGEST ACCEPTABLE
C                                     ARGUMENT
C                                     ZMAX1 = SMALLEST FLOATING-POINT
C                                     CONSTANT WITH ENTIRELY INTEGER
C                                     REPRESENTATION
C                                     ZMIN = 6 * 16**-32, SMALLEST
C                                     ACCEPTABLE ARGUMENT
C
      DATA PI2/Z40A2F9836E4E4415/, RTPI2/Z40CC42299EA1B284/,
     1     DLNHAF/-.6931471805599453D0/, ZMAX1/Z4E10000000000000/,
     2     TWOPI1/Z416487ED00000000/, TWOPI2/Z3B5110B4611A6263/,
     3     ZMIN/Z2160000000000000/,XINF/Z7FFFFFFFFFFFFFFF/,
     4     GNU1/Z3FCD000000000000/,GNU2/Z40FF000000000000/,
     5     TWOPI/Z416487ED5110B461/,ZMAX/Z4A10000000000000/
C                                  COEFFICIENTS FOR APPROXIMATION OF
C                                    GAMMA(V) 1.0 .LE. V .LE. 2.0
      DATA CO/Z40E2DFC48DA77B56/,
     1     GP/ZC0C1B71B59A1A1F6, Z41B33F20CFA73CB3, Z4153CF867C239860,
     2        ZC23EBA40FFB0397B, Z43441182D7048BE6, Z43C3CDE7AC8F2232,
     3        Z42E8A532ACC72020/,
     4     GQ/Z441C1A16BED21CC5,ZC1A5004D879829C5,Z41E62A3573ECF95D,
     5        Z42C97F1D84DC37A0,ZC327558408F56C71,ZC358DA535E278586,
     6        Z4411F52476FDA8AB/
C                                  THE FOLLOWING DATA STATEMENTS
C                                    INITIALIZE COEFFICIENTS FOR
C                                    RATIONAL APPROXIMATIONS TO
C                                    Y-SUB-NU(X) FOR X AN INTEGER, 0.0
C                                    .LE. NU .LE. 2.0. THE ARRAYS PN
C                                    INITIALIZE NUMERATOR COEFFICIENTS
C                                    FOR X = N, WHILE THE ARRAYS QN
C                                    INITIALIZE THE CORRESPONDING
C                                    DENOMINATOR POLYNOMIALS.
C
      DATA P1/Z421DF771FDF58468,ZC350302BA4F574FE,Z442D2D9BBC3AA125,
     1      ZC4EB52A3F9CA322F,Z451CAE84505B3A92,Z45194A03EE053C23,
     2      ZC6167F120A225FDF,Z46230DEBE8D9359A,Z448419E8CE903348,
     3      ZC672907741B98B8A,Z458852079422A3BA/
      DATA P4/Z4252C89C7BCABBAA,ZC34F0B2980A4941A,Z4416AEE56B9BE524,
     1      Z4453105E54B444A5,ZC537FA0B12F3E253,Z452A3A14110D4AF0,
     2      Z46287D9C9971B330,ZC61F5C1D23471D76,ZC67F4A47B6B7AC6D,
     3      Z4537BF40333A8675,0.0D0/
      DATA P7/ZC2266B1D2BBF5926,Z431FD8FE1F5A8790,ZC285E3B3FE09999E,
     1      ZC45752916AF445A1,Z447227424EC3F6F5,Z454B8C435E8F4B90,
     2      ZC522F544296CF497,ZC5E866630484B65A,ZC4C9D228E5F12965,
     3      0.0D0,0.0D0/
      DATA P10/Z422E9F0632CCBC94,ZC326407BF11101E7,ZC336BFB5F7103498,
     1      Z446F3F7AE43BD48B,ZC34DE558BD996722,ZC55BD71CBFD45D29,
     2      ZC46B9D80F95EC535,Z46108DE43FF97828,Z45253F202E245922,
     3      0.0D0,0.0D0/
      DATA P13/ZC2401F652EFE4D2F,Z432A4DC22FA367D1,Z4392496D26E7287E,
     1      ZC47BB08EDC2A076C,ZC47A423959B3C4B4,Z4563CE7798F3092A,
     2      Z453360489B59D88F,ZC6118E06B3BAE715,ZC541BBC67992A414,
     3      0.0D0,0.0D0/
      DATA P16/Z4240DE542FF81F75,ZC32347A42380418B,ZC3C1A014FB333352,
     1      Z446A3D371A09DDBB,Z44CC212A0AD27BE2,ZC55608A2FEA54B03,
     2      ZC5505D457D109148,Z45F0B9873C87B9E5,Z45518D3814D1A39E,
     3      0.0D0,0.0D0/
      DATA Q1/ZC210B9ECF8D3CAF6,Z424ACB46A661F0AE,ZC2959ED723485A5A,
     1      Z43706736B79B4811,ZC4185928C3EA9825,ZC417C1DF8C50E3E8,
     2      ZC51D8AC294F51F19,Z4547BCD6F6F49E08,Z4610A72AFC14703D,
     3      Z4660895E406910E1/
      DATA Q4/Z41241A8926BA098C,Z424728CC949ABF92,ZC28DF5A3E74A38BF,
     1      ZC32613A09DFDD62D,ZC48168B7A5131D49,ZC52BB920FA3F1F69,
     2      ZC61484D8A3372D24,ZC63BE9860E900D94,ZC6CDAB5F516DC122,
     3      0.0D0/
      DATA Q7/Z414DE04EE2A528F3,Z42757DB2D1B07F46,Z4327589F668C8ECA,
     1      Z441AB1398E924673,Z447504BF77FB59E3,Z453011FBB00B1719,
     2      Z457E460AB23DD9FD,Z461E6160759CE749,0.0D0,0.0D0/
      DATA Q10/Z411CFAD0391802EA,Z426EF0263B6026E3,Z431717B9845B9AA3,
     1      Z441970D3A7E56845,Z44560A9E1F1C01D5,Z453381AF429B7EA9,
     2      Z45730FDB3FB6278D,Z4629D0BD5F42D177,0.0D0,0.0D0/
      DATA Q13/ZC024B83B10B74E56,Z426F56EF9A37EBE3,Z42BD6ECE60913004,
     1      Z4419C94F59DA7450,Z443CAA8E96F869E2,Z4537C79CD1A291A3,
     2      Z45623B273AD672A2,Z463487F4D98B31B8,0.0D0,0.0D0/
      DATA Q16/ZC115A6A304E7C3E0,Z426BD317F430F3FC,Z4240A3FD66EDDBAB,
     1      Z441846D5FE976786,Z4425B04E60D6A640,Z4534FD61C8C29C2A,
     2      Z4547D2F83F593541,Z463532C077DD35ED,0.0D0,0.0D0/
C
C                                  (Z .LT. 0.03125, NU .LT. .05 .OR. (1
C                                    - NU) .LT. 0.004) COEFFICIENTS FOR
C                                    POWER SERIES IN NU FOR T1, T2, T3,
C                                    RESPECTIVELY. T1, T2, T3 ARE USED
C                                    IN COEFFICIENTS FOR POWER SERIES
C                                    FOR W(Z) = (Z/2)**NU * Y-SUB-NU.
C
      DATA F1/Z3C584516DCC3F036,Z3C349A2BD3510A7F,Z3D312CBB2A2FB1EE,
     1        Z3D289F0E334F2351,ZBE1C6A9A69DF63E4,Z3E994D04CEB6B407,
     2        ZBE371ACB2FAAFAC1,ZBF564514CE83E1E6,Z401325920BBECA39,
     3        ZBF18BC2D963119F1,ZC06B00D27772FECF,Z40AB2C53740274B2,
     4        Z40A35D3B1C9C5EE4,ZC117090683A568FF,ZC0227399E50FDCE2,
     5        Z40517CC1B727220B,Z3C2B67790859DC55/
      DATA F2/Z3B99C72A9FA89929,Z3D180DD6520470F9,Z3D1601C0FD034655,
     1        Z3D96B2C61EC84773,ZBE217127E6F50C37,Z3E8C23F9D32DB02D,
     2        Z3E44C887195FF983,ZBF55C597F0E803B8,Z3FD59DADC096517C,
     3        Z4012BDFFE2483F74,ZC06D7AD5A4D964A8,Z40463DF6BEE5AE41,
     4        Z4114F226C3032B23,ZC0B249A1AE604540,ZC11921FB54442D18,
     5        Z405E124FA42E8A50/
      DATA F3/ZBC515C25745D59C5,Z3CA2A1FCFBD4BB1A,ZBD144A7D1AA4308F,
     1        Z3D289C7FA6A541A5,ZBD506B665E0466E4,Z3DA3AD90876199B2,
     2        ZBE132ED97D35A628,Z3E2BE3BF8C333779,ZBE3B822340F3CB76,
     3        Z3EE8AE8CAB930E7C,ZBD591C757C12BBE8,Z3F60CE3FAFE7732C,
     4        Z3F6A5C19A9255CC1,Z40218F4D09EC4196,Z40227399E50FDCE2,
     5        Z40517CC1B727220B,Z3C28B3F7C6C5E933/
C
C                                  1.0 / (FACTORIAL(N))
      DATA FACT/1.0D0,1.0D0,Z4080000000000000,Z402AAAAAAAAAAAAB,
     1      Z3FAAAAAAAAAAAAAB,Z3F22222222222222,Z3E5B05B05B05B05B,
     2      Z3DD00D00D00D00D0,Z3D1A01A01A01A01A,Z3C2E3BC74AAD8E67,
     3      Z3B49F93EDDE27D72,Z3A6B99159FD5138E,Z398F76C77FC6C4BE,
     4      Z38B092309D43684C,Z37C9CBA54603E4E9,Z36D73F9F399DC0F9,
     5      Z35D73F9F399DC0F9,Z34CA963B81856A53,Z33B413C31DCBECBC,
     6      Z3297A4DA340A0AB9,Z317950AE90080894,Z305C6E3BDB73D5C6,
     7      Z2F4338E5B6DFE14A,Z2E2EC368262C7034,Z2D1F2CF01972F578,
     8      Z2C13F3CCDD165FA9,Z2AC4742FE35272CD,Z29746AC70B733A8D,
     9      Z2842862898D42175,Z2724B3F31686B15B,Z2613932C5047D60E,
     *      Z24A1A6973C1FADE2,Z2350D34B9E0FD6F1,Z22273024A9BA1AA3,
     *      Z2112710231C0FD7A,Z1F86E2CE38B6C8F9/
C                                  FIRST EXECUTABLE STATEMENT
      ZSTOP = 17.5D0
      CON1 = 21.0D0
      CON2 = 25.0D0
      IEND = 16
      MEND = 14
      M1 = 13
      M2 = 9
      M3 = 4
      IENDP2 = IEND+2
      MENDM3 = MEND-3
      Z = ARG
      IER = 0
      V = ORDER
      NFLAG = (N.EQ.1)
      IOPT = 2
      IF (NFLAG) IOPT = 1
      IF (Z.GE.ZMIN) GO TO 5
C                                  ERROR RETURN FOR Z .LT. ZMIN
      IER = 129
    5 IF ((V.GE.0.0D0).AND.(V.LT.1.0D0)) GO TO 10
C
C                                  ERROR RETURN FOR V .LT. 0.0
      IER = 129
   10 IF (N.GT.0) GO TO 15
C                                  ERROR RETURN FOR N .LE. 0
      IER = 130
   15 IF (IER.EQ.0) GO TO 25
      J = N
      IF (IER.EQ.130) J = 1
      DO 20 I=1,J
         YN(I) = -XINF
   20 CONTINUE
      GO TO 9000
   25 IF (V.LT.ZMIN) V = 0.0D0
      IF (Z.LT..03125D0) GO TO 70
      IF (Z.GE.ZSTOP) GO TO 145
C                                  IF (1/32 .LE. Z .LT. 20.5), (X**NU)
C                                    * Y-SUB-NU(Z) IS APPROXIMATED BY A
C                                    SEQUENCE OF TAYLOR SERIES
C                                    EXPANSIONS, THE FIRST ABOUT A
C                                    GIVEN X, AND THE NEXT ABOUT X +
C                                    DELX. THE SEQUENCE IS CONTINUED
C                                    UNTIL Z IS WITHIN DELX OF THE
C                                    POINT OF EXPANSION.
      J = IDINT((Z+3.5D0)/3.0D0)
      X = DBLE(FLOAT(J))*3.0D0-2.0D0
      M = MAX0(9,12-J)
C                                  APPROXIMATION AT STARTING POINT X IS
C                                    OBTAINED AS A MINIMAX RATIONAL
      GNU = V+1.0D0
      XNUM = P(1,J)*V+P(2,J)
      XNUMP = XNUM+P(1,J)
      DEN = V+Q(1,J)
      DENP = DEN+1.0D0
      DO 30 I=3,M
         XNUM = XNUM*V+P(I,J)
         DEN = DEN*V+Q(I-1,J)
         XNUMP = XNUMP*GNU+P(I,J)
         DENP = DENP*GNU+Q(I-1,J)
   30 CONTINUE
   35 DELX = .0625D0
      IF (Z.GT.4.0D0) DELX = DELX+DELX
      T = Z-X
      IF (T.LT.0.0D0) DELX = -DELX
C                                  EXPAND IN TERMS OF THE AUXILIARY
C                                    FUNCTION W-SUB-NU(Z) = Z**V *
C                                    Y-SUB-NU(Z)
      M = MEND
      IF (Z.GT.4.0D0) M = MENDM3
      XNU = 1.0D0
      IF (X.NE.1.0D0) XNU = (X**V)*X
      W(1) = XNUMP/DENP*XNU
      W(2) = XNUM/DEN*XNU
      XNU = 2.0D0*GNU-1.0D0
      IF (Z.GE..5D0) GO TO 40
      T = Z+Z-1.5D0
      IF (Z.GE..25D0) GO TO 40
      T = 4.0D0*Z-2.0D0
      IF (Z.LT..125D0) T = 8.0D0*Z-2.5D0
   40 I = IDINT(T/DELX)+1
      C2 = 0.0D0
      C2P = 0.0D0
      CUM = W(1)
      CUMP = W(2)
      NMAX = M+3
      DO 65 NN=1,I
         IF (NN.LT.I) GO TO 45
         DELX = Z-X
         IF (DELX.EQ.0.0D0) GO TO 65
   45    IF (Z.GE.1.0D0) GO TO 50
         IF (NN.EQ.9) DELX = DMAX1(DELX,-.03125D0)
         IF (NN.EQ.17) DELX = DMAX1(DELX,-.015625D0)
         IF (NN.EQ.25) DELX = DMAX1(DELX,-.0078125D0)
C
C                                  CALCULATION OF DERIVATIVES OF W
   50    XIN = 1.0D0/X
         W(3) = XNU*XIN*W(2)-W(1)
         XK = 1.0D0
         S = XNU
         DO 55 K=4,NMAX
            S = S-1.0D0
            W(K) = (S*W(K-1)-XK*W(K-3))*XIN-W(K-2)
            XK = XK+1.0D0
   55    CONTINUE
         SUM = FACT(M+2)*W(M+2)*DELX
         SUMP = FACT(M+2)*W(M+3)*DELX
         J = M+1
         DO 60 K=1,M
            SUM = (SUM+FACT(J)*W(J))*DELX
            SUMP = (SUMP+FACT(J)*W(J+1))*DELX
            J = J-1
   60    CONTINUE
C                                  A SUMMING SCHEME SUGGESTED BY KAHAN
C                                    IS USED TO PRESERVE ACCURACY.
C                                    FOOLU IS AN EXTRA VARIABLE TO
C                                    NEGATE POSSIBLE OPTIMIZATION ERROR
         C2 = C2+SUM
         T = CUM+C2
         FOOLU = CUM-T
         C2 = C2+FOOLU
         CUM = T
         C2P = C2P+SUMP
         T = CUMP+C2P
         FOOLU = CUMP-T
         C2P = C2P+FOOLU
         CUMP = T
         X = X+DELX
         W(1) = CUM
         W(2) = CUMP
   65 CONTINUE
      XNU = (Z**V)*Z
      YN(1) = W(2)/XNU
      IF (NFLAG) GO TO 9005
      YN(2) = W(1)/XNU
      GO TO 170
C                                  FOR (Z .LT. .03125, .05 .LE. NU .LE.
C                                    .996), Y-SUB-NU = (J-SUB-NU *
C                                    COS(NU * PI) - J-SUB-(-NU) ) /
C                                    SIN(NU * PI)
   70 ZH = Z*0.5D0
      X = ZH*ZH
      IF ((V.LT.GNU1).OR.(V.GT.GNU2)) GO TO 95
      X = -X
      GNU = V
      VPI = V/PI2
      IF (V.GT.0.5D0) VPI = (1.0D0-V)/PI2
      VPI = VPI+VPI
      VSIN = DSIN(VPI)
      VCOS = DCOS(VPI)
      IF (V.GT.0.5D0) VCOS = -VCOS
      T1 = ZH**GNU
      T2 = 1.0D0/T1
      DO 90 I=1,2
         XC = T1
         DO 85 J=1,2
            T = GNU-0.5D0
            IF (J.EQ.2) T = T+2.0D0
            IF (T.GT.0.5D0) T = T-1.0D0
C                                  EVALUATION OF GAMMA(GNU + 1.0)
            XNUM = GP(1)*T
            DEN = T
            DO 75 L=2,7
               XNUM = (XNUM+GP(L))*T
               DEN = (DEN+GQ(L))*T
   75       CONTINUE
            DEN = DEN+GQ(1)
            GAM = CO+XNUM/DEN
            IF (J.LT.I) GAM = GAM*GNU
            IF (J.EQ.2) GAM = GAM/(GNU+1.0D0)
            IF (I*J.EQ.4) GAM = GAM/(GNU+2.0D0)
C
C                                  EVALUATION OF J-SUB-GNU
            SUM = FACT(5)
            XK = 4.0D0
            DO 80 K=1,4
               SUM = SUM*X/(GNU+XK)+FACT(5-K)
               XK = XK-1.0D0
   80       CONTINUE
            SCRAT(J) = XC*SUM/GAM
            IF(I*J.EQ.4) SCRAT(J) = SCRAT(J)/ZH
            GNU = -GNU
            XC = T2
   85    CONTINUE
         YN(I) = (VCOS*SCRAT(1)-SCRAT(2))/VSIN
         IF (NFLAG) GO TO 9005
         IF (I.EQ.2) GO TO 170
         VCOS = -VCOS
         VSIN = -VSIN
         GNU = V+1.0D0
         T1 = T1*ZH
   90 CONTINUE
      GO TO 170
C                                  FOR (Z .LT. .03125, NU .LT. .05 OR
C                                    NU .GT. .996) USE FULLERTON
C                                    EXPANSIONS OF W-SUB-NU(Z) = (Z /
C                                    2) ** V * Y-SUB-NU(Z)
   95 GNU = V
      IF (V.GT.GNU2) GNU = V - 1.0D0
      T1 = F1(17)
      T2 = 0.0D0
      T3(1) = F3(17)
      DO 100 I=1,16
         T1 = T1*GNU+F1(I)
         T2 = T2*GNU+F2(I)
         T3(1) = T3(1)*GNU+F3(I)
  100 CONTINUE
      T1 = T1/(1.0D0-GNU)
      T2 = T2/(1.0D0-GNU*GNU)
      T3(1) = T3(1)/(1.0D0+GNU)
      DLNX = 2.0D0*DLOG(ZH)
      DLNXNU = GNU*DLNX
      IF (DLNXNU.LT.DLNHAF) GO TO 110
C                                  IF X ** V .GT. 0.5 COMPUTE (X ** V
C                                    -1.0) / V USING POWER SERIES FOR E
C                                    ** (V * LOG(X))
      XNU = FACT(IENDP2)
      DO 105 I=1,IEND
         XNU = (XNU*DLNXNU)+FACT(IENDP2-I)
  105 CONTINUE
      XNU = DLNX*XNU
      SCRAT(1) = T2+XNU*T1
      XNU = GNU*XNU+1.0D0
      GO TO 115
  110 XNU = ZH**(GNU+GNU)
      SCRAT(1) = (XNU*T1-T3(1))/GNU
  115 XK = 1.0D0
C                                  SCRAT CONTAINS COEFFICIENTS FOR
C                                    EXPANSION W = SUM (SCRAT(K+1) *
C                                    (-X) ** K), K = 0 TO 4
      DO 120 I=2,6
         T = XK-GNU
         T3(I) = T3(I-1)/(T*XK)
         SCRAT(I) = (SCRAT(I-1)*(1.0D0-GNU/XK)/T-2.0D0*T3(I))/(XK+GNU)
         XK = XK+1.0D0
  120 CONTINUE
      XK = 5.0D0+GNU
      SUM = SCRAT(5)
      SUMP = XK*SCRAT(6)+T3(6)
      IF (X.LT.1.0D-30) X = 0.0D0
      DO 125 I=1,4
         K = 5-I
         XK = XK-1.0D0
         SUMP = (-SUMP)*X+XK*SCRAT(K+1)+T3(K+1)
         SUM = SUM*(-X)+SCRAT(K)
  125 CONTINUE
      SUMP = XNU*T1/ZH-SUMP*ZH
      RTXNU = DSQRT(XNU)
      XNU = GNU/ZH
      SUM = SUM/RTXNU
      IF (GNU.LT.0.0D0) GO TO 130
      YN(1) = SUM
      IF (NFLAG) GO TO 9005
      YN(2) = XNU*YN(1)-SUMP/RTXNU
      GO TO 170
C                                  IF GNU = 1.0 - V, CALCULATE
C                                    Y(Z,GNU+1)
  130 YN(1) = XNU*SUM-SUMP/RTXNU
      IF (NFLAG) GO TO 9005
      XNU = XNU+2.0D0/Z
C                                  CHECK FOR POTENTIAL OVERFLOW IN
C                                    COMPUTING SEQUENCE OF VALUES
      IF (XNU.GE.-XINF/YN(1)) GO TO 135
      YN(2) = XNU*YN(1)-SUM
      GO TO 170
  135 DO 140 I=2,N
         YN(I) = -XINF
  140 CONTINUE
      IER = 131
      GO TO 9000
C                                  USE AN ASYMPTOTIC SERIES FOR (Z .GE.
C                                    ZSTOP)
  145 IF (Z.LE.ZMAX) GO TO 155
      IER = 132
      DO 150 I=1,N
         YN(I) = 0.0D0
  150 CONTINUE
      GO TO 9000
  155 XC = RTPI2/DSQRT(Z)
      XIN = .015625D0/(Z*Z)
      M = 17
      IF (Z.GE.CON1) M = M1
      IF (Z.GE.CON2) M = M2
      IF (Z.GE.130.0D0) M = M3
      XM = 4.0D0*DBLE(FLOAT(M))
C                                  REDUCTION OF ARGUMENT FOR SIN AND
C                                    COS ROUTINES
      T = Z/TWOPI+ZMAX1
      T = T-ZMAX1
  158 X = ((Z-T*TWOPI1)-T*TWOPI2)-(V+0.5D0)/PI2
  159 VSIN = DSIN(X)
      VCOS = DCOS(X)
      GNU = V+V
      DO 165 I=1,2
         S = ((XM-1.0D0)-GNU)*((XM-1.0D0)+GNU)*XIN*0.5D0
         T = (GNU-(XM-3.0D0))*(GNU+(XM-3.0D0))
         XNUM = S*T*FACT(2*M+1)
         T1 = (GNU-(XM+1.0D0))*(GNU+(XM+1.0D0))
         DEN = S*T1*FACT(2*M+2)
         XK = XM
         K = M+M
         T1 = T
         DO 160 J=2,M
            XK = XK-4.0D0
            S = ((XK-1.0D0)-GNU)*((XK-1.0D0)+GNU)
            T = (GNU-(XK-3.0D0))*(GNU+(XK-3.0D0))
            XNUM = (XNUM+FACT(K-1))*S*T*XIN
            DEN = (DEN+FACT(K))*S*T1*XIN
            K = K-2
            T1 = T
  160    CONTINUE
         XNUM = XNUM+1.0D0
         DEN = (DEN+1.0D0)*(GNU*GNU-1.0D0)/(8.0D0*Z)
         YN(I) = XC*(XNUM*VSIN+DEN*VCOS)
         IF (NFLAG) GO TO 9005
         T = VSIN
         VSIN = -VCOS
         VCOS = T
         GNU = GNU+2.0D0
  165 CONTINUE
C                                  IF N .GT. 2, COMPUTE Y(Z,ORDER+I) I
C                                    = 2, N-1
  170 IF (N.LE.2) GO TO 9005
      GNU = V+V+2.0D0
      NN = N
      DO 185 J=3,NN
         IF (DABS(YN(J-1)).LT.1.0D0) GO TO 180
C                                  CHECK FOR POTENTIAL OVERFLOW IN
C                                    COMPUTING SEQUENCE OF VALUES
C
         IF (DABS(GNU/Z).LT.DABS(XINF/YN(J-1))) GO TO 180
         DO 175 L=J,N
            YN(L) = -XINF
  175    CONTINUE
         IER = 131
         GO TO 9000
  180    YN(J) = GNU*YN(J-1)/Z-YN(J-2)
         GNU = GNU+2.0D0
  185 CONTINUE
      GO TO 9005
 9000 CONTINUE
      CALL UERTST (IER,6HMMBSYN)
 9005 RETURN
      END
