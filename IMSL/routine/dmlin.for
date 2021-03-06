C   IMSL ROUTINE NAME   - DMLIN
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1982
C
C   PURPOSE             - NUMERICAL INTEGRATION OF A FUNCTION OF
C                           SEVERAL VARIABLES OVER A HYPER-RECTANGLE
C                           (GAUSSIAN METHOD)
C
C   USAGE               - FUNCTION DMLIN (F,A,B,N,MAXFCN,AERR,RERR,IER)
C
C   ARGUMENTS    DMLIN  - ESTIMATE OF THE INTEGRAL OF F
C                           OVER THE HYPER-RECTANGLE
C                               A(I).LE.X(I).LE.B(I), I=1...N
C                           (OUTPUT).
C                F      - A REAL FUNCTION SUBPROGRAM SUPPLIED BY THE
C                           USER. (INPUT) F DEFINES THE FUNCTION THAT
C                           IS TO BE INTEGRATED. F MUST BE DECLARED
C                           EXTERNAL IN THE CALLING PROGRAM AND MUST
C                           BE OF THE FORM
C                             REAL FUNCTION F(N,X)
C                             REAL X(N)
C                             F= ...
C                             RETURN
C                             END
C                           X MUST NOT BE ALTERED BY F.
C                A,B    - REAL VECTORS OF LENGTH N SUCH THAT THE
C                           HYPER-RECTANGLE IS DEFINED BY
C                                A(I).LE.X(I).LE.B(I), I=1...N
C                           (INPUT).
C                N      - NUMBER OF INDEPENDENT VARIABLES, OR
C                           DI             REAL X(N)
C                             F= ...
C                             RETURN
C                             END
C                           X MUST NOT BE ALTERED BY F.
C                A,B    - REAL VECTORS OF LENGTH N SUCH THAT THE
C                           HYPER-RECTANGLE IS DEFINED BY
C                                A(I).LE.X(I).LE.B(I), I=1...N
C                           (INPUT).
C                N      - NUMBER OF INDEPENDENT VARIABLES, OR
C                           DI         TERMINAL ERROR
C                           IER = 129 INDICATES THAT CONVERGENCE
C                             HAD NOT OCCURED AFTER MIN(MAXFCN,256**N)
C                             FUNCTION EVALUATIONS.
C                           IER = 130 INDICATES THAT N.GT.20
C
C   REQD. IMSL ROUTINES - UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS  1.  IF THE INTEGRAND IS RELATIVELY SMOOTH, AND/OR
C                FUNCTION EVALUATIONS ARE EXPENSIVE, DMLIN
C                MAY BE MORE EFFICIENT THAN OTHER METHODS,
C                EVEN FOR ONE AND TWO DIMENSIONAL PROBLEMS.
C                IF THE INTEGRAND IS BADLY BEHAVED, AN
C                ADAPTIVE SUBROUTINE MAY BE PREFERABLE FOR
C                ONE OR TWO DIMENSIONAL PROBLEMS.
C            2.  DMLIN TERMINATES WHEN EITHER THE ABSOLUTE
C                OR RELATIVE ERROR CRITERION APPEARS TO BE
C                SATISFIED.
C
C-----------------------------------------------------------------------
C
      REAL FUNCTION DMLIN (F,A,B,N,MAXFCN,AERR,RERR,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N,MAXFCN,IER
      REAL               A(N),B(N),AERR,F,RERR
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,IRES,J,M(20),NLIM,NN,NN2,NN2I,NNI
      REAL               BMA(20),BPA(20),BWT,RESULT(8),T(255),TINT,
     *                   TT(256),W(255),WT,WW(256),X(20)
      DATA T(1),T(2),T(3),T(4),T(5),T(6),T(7),T(8),T(9),T(10),T(11),
     *     T(12),T(13),T(14),T(15)/.5773503E0,
     *  .8611363E0,.3399810E0,.9602899E0,
     *  .7966665E0,.5255324E0,.1834346E0,
     *  .9894009E0,.9445750E0,.8656312E0,
     *  .7554044E0,.6178762E0,.4580168E0,
     *  .2816036E0,.9501251E-1/
      DATA T(16),T(17),T(18),T(19),T(20),T(21),T(22),T(23),T(24),T(25),
     *     T(26),T(27),T(28),T(29),T(30),T(31)/.9972639E0,
     *  .9856115E0,.9647623E0,.9349061E0,
     *  .8963212E0,.8493676E0,.7944838E0,
     *  .7321821E0,.6630443E0,.5877158E0,
     *  .5068999E0,.4213513E0,.3318686E0,
     *  .2392874E0,.1444720E0,.4830767E-1/
      DATA T(32),T(33),T(34),T(35),T(36),T(37),T(38),T(39),T(40),T(41),
     *     T(42),T(43),T(44),T(45),T(46),T(47)/.9993050E0,
     *  .9963401E0,.9910134E0,.9833363E0,
     *  .9733268E0,.9610088E0,.9464114E0,
     *  .9295692E0,.9105221E0,.8893154E0,
     *  .8659994E0,.8406293E0,.8132653E0,
     *  .7839724E0,.7528199E0,.7198819E0/
      DATA T(48),T(49),T(50),T(51),T(52),T(53),T(54),T(55),T(56),T(57),
     *     T(58),T(59),T(60),T(61),T(62),T(63)/.6852363E0,
     *  .6489655E0,.6111554E0,.5718956E0,
     *  .5312795E0,.4894031E0,.4463660E0,
     *  .4022702E0,.3572202E0,.3113229E0,
     *  .2646872E0,.2174236E0,.1696444E0,
     *  .1214628E0,.7299312E-1,.2435029E-1/
      DATA T(64),T(65),T(66),T(67),T(68),T(69),T(70),T(71),T(72),T(73),
     *     T(74),T(75),T(76),T(77),T(78),T(79)/.9998249E0,
     *  .9990775E0,.9977332E0,.9957928E0,
     *  .9932571E0,.9901278E0,.9864067E0,
     *  .9820961E0,.9771985E0,.9717168E0,
     *  .9656544E0,.9590148E0,.9518020E0,
     *  .9440203E0,.9356744E0,.9267693E0/
      DATA T(80),T(81),T(82),T(83),T(84),T(85),T(86),T(87),T(88),T(89),
     *     T(90),T(91),T(92),T(93),T(94),T(95),T(96),T(97)/
     *  .9173102E0,.9073029E0,.8967533E0,
     *  .8856677E0,.8740528E0,.8619155E0,
     *  .8492630E0,.8361029E0,.8224431E0,
     *  .8082918E0,.7936573E0,.7785485E0,
     *  .7629743E0,.7469442E0,.7304676E0,
     *  .7135544E0,.6962147E0,.6784589E0/
      DATA T(98),T(99),T(100),T(101),T(102),T(103),T(104),T(105),T(106),
     *     T(107),T(108),T(109),T(110),T(111),T(112),T(113),T(114)/
     *  .6602976E0,.6417417E0,.6228022E0,
     *  .6034905E0,.5838180E0,.5637966E0,
     *  .5434383E0,.5227552E0,.5017596E0,
     *  .4804641E0,.4588814E0,.4370245E0,
     *  .4149064E0,.3925403E0,.3699396E0,
     *  .3471177E0,.3240884E0/
      DATA T(115),T(116),T(117),T(118),T(119),T(120),T(121),T(122),
     *     T(123),T(124),T(125),T(126),T(127)/.3008654E0,
     *  .2774626E0,.2538940E0,.2301736E0,
     *  .2063156E0,.1823343E0,.1582440E0,
     *  .1340592E0,.1097942E0,
     *  .8546364E-1,.6108197E-1,.3666379E-1,
     *  .1222370E-1/
      DATA T(128),T(129),T(130),T(131),T(132),T(133),T(134),T(135),
     *     T(136),T(137),T(138),T(139),T(140),T(141),T(142)/
     *  .9999561E0,.9997684E0,.9994309E0,
     *  .9989435E0,.9983063E0,.9975193E0,
     *  .9965826E0,.9954965E0,.9942610E0,
     *  .9928763E0,.9913428E0,.9896605E0,
     *  .9878297E0,.9858508E0,.9837240E0/
      DATA T(143),T(144),T(145),T(146),T(147),T(148),T(149),T(150),
     *     T(151),T(152),T(153),T(154),T(155),T(156),T(157)/
     *  .9814496E0,.9790280E0,.9764595E0,
     *  .9737446E0,.9708836E0,.9678769E0,
     *  .9647251E0,.9614285E0,.9579877E0,
     *  .9544032E0,.9506755E0,.9468052E0,
     *  .9427929E0,.9386392E0, .9343446E0/
      DATA T(158),T(159),T(160),T(161),T(162),T(163),T(164),T(165),
     *     T(166),T(167),T(168),T(169),T(170),T(171),T(172)/
     *  .9299099E0,.9253357E0,.9206227E0,
     *  .9157716E0,.9107831E0,.9056580E0,
     *  .9003970E0,.8950010E0,.8894707E0,
     *  .8838069E0,.8780106E0,.8720826E0,
     *  .8660238E0,.8598350E0,.8535173E0/
      DATA T(173),T(174),T(175),T(176),T(177),T(178),T(179),T(180),
     *     T(181),T(182),T(183),T(184),T(185),T(186),T(187)/
     *  .8470715E0,.8404987E0,.8337997E0,
     *  .8269757E0,.8200277E0,.8129566E0,
     *  .8057636E0,.7984497E0,.7910160E0,
     *  .7834637E0,.7757938E0,.7680076E0,
     *  .7601062E0,.7520907E0,.7439624E0/
      DATA T(188),T(189),T(190),T(191),T(192),T(193),T(194),T(195),
     *     T(196),T(197),T(198),T(199),T(200),T(201),T(202)/
     *  .7357225E0,.7273723E0,.7189129E0,
     *  .7103457E0,.7016719E0,.6928929E0,
     *  .6840099E0,.6750243E0,.6659375E0,
     *  .6567508E0,.6474655E0,.6380831E0,
     *  .6286050E0,.6190327E0,.6093674E0/
      DATA T(203),T(204),T(205),T(206),T(207),T(208),T(209),T(210),
     *     T(211),T(212),T(213),T(214),T(215),T(216),T(217)/
     *  .5996107E0,.5897641E0,.5798290E0,
     *  .5698070E0,.5596994E0,.5495079E0,
     *  .5392340E0,.5288792E0,.5184450E0,
     *  .5079331E0,.4973450E0,.4866822E0,
     *  .4759465E0,.4651394E0,.4542624E0/
      DATA T(218),T(219),T(220),T(221),T(222),T(223),T(224),T(225),
     *     T(226),T(227),T(228),T(229),T(230),T(231),T(232)/
     *  .4433174E0,.4323058E0,.4212294E0,
     *  .4100898E0,.3988887E0,.3876278E0,
     *  .3763087E0,.3649331E0,.3535028E0,
     *  .3420195E0,.3304849E0,.3189007E0,
     *  .3072686E0,.2955905E0,.2838680E0/
      DATA T(233),T(234),T(235),T(236),T(237),T(238),T(239),T(240),
     *     T(241),T(242),T(243),T(244),T(245),T(246),T(247)/
     *  .2721029E0,.2602971E0,.2484521E0,
     *  .2365699E0,.2246523E0,.2127009E0,
     *  .2007176E0,.1887042E0,.1766625E0,
     *  .1645943E0,.1525014E0,.1403856E0,
     *  .1282488E0,.1160927E0,.1039192E0/
      DATA   T(248),T(249),T(250),T(251),T(252),T(253),T(254),T(255)/
     *  .9173013E-1,.7952729E-1,.6731252E-1,
     *  .5508766E-1,.4285453E-1,.3061497E-1,
     *  .1837082E-1,.6123912E-2/
      DATA W(1),W(2),W(3),W(4),W(5),W(6),W(7),W(8),W(9),W(10),W(11),
     *     W(12),W(13),W(14),W(15)/1.0E0,.3478548E0,
     *  .6521452E0,.1012285E0,.2223810E0,
     *  .3137066E0,.3626838E0,.2715246E-1,
     *  .6225352E-1,.9515851E-1,.1246290E0,
     *  .1495960E0,.1691565E0,.1826034E0,
     *  .1894506E0/
      DATA W(16),W(17),W(18),W(19),W(20),W(21),W(22),W(23),W(24),W(25),
     *     W(26),W(27),W(28),W(29),W(30),W(31)/.7018610E-2,
     *  .1627439E-1,.2539207E-1,.3427386E-1,
     *  .4283590E-1,.5099806E-1,.5868409E-1,
     *  .6582222E-1,.7234579E-1,.7819390E-1,
     *  .8331192E-1,.8765209E-1,.9117388E-1,
     *  .9384440E-1,.9563872E-1,.9654009E-1/
      DATA W(32),W(33),W(34),W(35),W(36),W(37),W(38),W(39),W(40),W(41),
     *     W(42),W(43),W(44),W(45),W(46),W(47)/.1783281E-2,
     *  .4147033E-2,.6504458E-2,.8846760E-2,
     *  .1116814E-1,.1346305E-1,.1572603E-1,
     *  .1795172E-1,.2013482E-1,.2227017E-1,
     *  .2435270E-1,.2637747E-1,.2833967E-1,
     *  .3023466E-1,.3205793E-1,.3380516E-1/
      DATA W(48),W(49),W(50),W(51),W(52),W(53),W(54),W(55),W(56),W(57),
     *     W(58),W(59),W(60),W(61),W(62),W(63)/ .3547221E-1,
     *  .3705513E-1,.3855015E-1,.3995374E-1,
     *  .4126256E-1,.4247352E-1,.4358372E-1,
     *  .4459056E-1,.4549163E-1, .4628480E-1,
     *  .4696818E-1,.4754017E-1,.4799939E-1,
     *  .4834476E-1,.4857547E-1,.4869096E-1/
      DATA W(64),W(65),W(66),W(67),W(68),W(69),W(70),W(71),W(72),W(73),
     *     W(74),W(75),W(76),W(77),W(78),W(79)/ .4493810E-3,
     *  .1045813E-2,.1642503E-2,.2238288E-2,
     *  .2832751E-2,.3425526E-2,.4016255E-2,
     *  .4604584E-2,.5190162E-2,.5772638E-2,
     *  .6351663E-2,.6926893E-2,.7497982E-2,
     *  .8064590E-2,.8626378E-2,.9183010E-2/
      DATA W(80),W(81),W(82),W(83),W(84),W(85),W(86),W(87),W(88),
     *   W(89),W(90),W(91),W(92),W(93),W(94),W(95)/.9734153E-2,
     *  .1027948E-1,.1081866E-1,.1135138E-1,
     *  .1187731E-1,.1239614E-1,.1290756E-1,
     *  .1341127E-1,.1390696E-1,.1439435E-1,
     *  .1487312E-1,.1534301E-1,.1580373E-1,
     *  .1625500E-1,.1669656E-1,.1712814E-1/
      DATA W(96),W(97),W(98),W(99),W(100),W(101),W(102),W(103),W(104),
     *     W(105),W(106),W(107),W(108),W(109),W(110),W(111),W(112)/
     *  .1754948E-1,.1796033E-1,.1836044E-1,
     *  .1874959E-1,.1912752E-1,.1949403E-1,
     *  .1984888E-1,.2019187E-1,.2052279E-1,
     *  .2084145E-1,.2114765E-1,.2144121E-1,
     *  .2172195E-1,.2198971E-1,.2224433E-1,
     *  .2248565E-1,.2271354E-1/
      DATA W(113),W(114),W(115),W(116),W(117),W(118),W(119),W(120),
     *     W(121),W(122),W(123),W(124),W(125),W(126),W(127)/
     *  .2292784E-1,.2312845E-1,.2331523E-1,
     *  .2348808E-1,.2364688E-1,.2379156E-1,
     *  .2392201E-1,.2403817E-1,.2413996E-1,
     *  .2422732E-1,.2430020E-1,.2435856E-1,
     *  .2440236E-1,.2443157E-1,.2444618E-1/
      DATA W(128),W(129),W(130),W(131),W(132),W(133),W(134),W(135),
     *     W(136),W(137),W(138),W(139),W(140),W(141),W(142)/
     *  .1127890E-3,.2625349E-3,.4124633E-3,
     *  .5623490E-3,.7121542E-3,.8618537E-3,
     *  .1011424E-2,.1160844E-2,.1310089E-2,
     *  .1459137E-2,.1607967E-2,.1756556E-2,
     *  .1904881E-2,.2052920E-2,.2200652E-2/
      DATA W(143),W(144),W(145),W(146),W(147),W(148),W(149),W(150),
     *     W(151),W(152),W(153),W(154),W(155),W(156),W(157)/
     *  .2348053E-2,.2495102E-2,.2641777E-2,
     *  .2788055E-2,.2933916E-2,.3079336E-2,
     *  .3224294E-2,.3368769E-2,.3512738E-2,
     *  .3656180E-2,.3799074E-2,.3941398E-2,
     *  .4083130E-2,.4224250E-2,.4364737E-2/
      DATA W(158),W(159),W(160),W(161),W(162),W(163),W(164),W(165),
     *     W(166),W(167),W(168),W(169),W(170),W(171),W(172)/
     *  .4504569E-2,.4643725E-2,.4782184E-2,
     *  .4919926E-2,.5056930E-2,.5193175E-2,
     *  .5328642E-2,.5463309E-2,.5597156E-2,
     *  .5730164E-2,.5862312E-2,.5993581E-2,
     *  .6123951E-2,.6253402E-2,.6381915E-2/
      DATA W(173),W(174),W(175),W(176),W(177),W(178),W(179),W(180),
     *     W(181),W(182),W(183),W(184),W(185),W(186),W(187)/
     *  .6509470E-2,.6636050E-2,.6761633E-2,
     *  .6886203E-2,.7009739E-2,.7132224E-2,
     *  .7253639E-2,.7373966E-2,.7493186E-2,
     *  .7611283E-2,.7728238E-2,.7844033E-2,
     *  .7958652E-2,.8072077E-2,.8184291E-2/
      DATA W(188),W(189),W(190),W(191),W(192),W(193),W(194),W(195),
     *     W(196),W(197),W(198),W(199),W(200),W(201),W(202)/
     *  .8295278E-2,.8405020E-2,.8513501E-2,
     *  .8620705E-2,.8726616E-2,.8831218E-2,
     *  .8934495E-2,.9036432E-2,.9137013E-2,
     *  .9236223E-2,.9334048E-2,.9430473E-2,
     *  .9525483E-2,.9619065E-2,.9711203E-2/
      DATA W(203),W(204),W(205),W(206),W(207),W(208),W(209),W(210),
     *     W(211),W(212),W(213),W(214),W(215),W(216),W(217)/
     *  .9801885E-2,.9891096E-2,.9978823E-2,
     *  .1006505E-1,.1014977E-1,.1023297E-1,
     *  .1031464E-1,.1039475E-1,.1047331E-1,
     *  .1055029E-1,.1062570E-1,.1069950E-1,
     *  .1077171E-1,.1084230E-1,.1091126E-1/
      DATA W(218),W(219),W(220),W(221),W(222),W(223),W(224),W(225),
     *     W(226),W(227),W(228),W(229),W(230),W(231),W(232)/
     *  .1097858E-1,.1104426E-1,.1110828E-1,
     *  .1117063E-1,.1123131E-1,.1129031E-1,
     *  .1134761E-1,.1140321E-1,.1145709E-1,
     *  .1150926E-1,.1155970E-1,.1160841E-1,
     *  .1165538E-1,.1170060E-1,.1174406E-1/
      DATA W(233),W(234),W(235),W(236),W(237),W(238),W(239),W(240),
     *     W(241),W(242),W(243),W(244),W(245),W(246),W(247)/
     *  .1178576E-1,.1182570E-1,.1186386E-1,
     *  .1190024E-1,.1193483E-1,.1196764E-1,
     *  .1199865E-1,.1202785E-1,.1205526E-1,
     *  .1208086E-1,.1210464E-1,.1212661E-1,
     *  .1214676E-1,.1216509E-1,.1218159E-1/
      DATA W(248),W(249),W(250),W(251),W(252),W(253),W(254),W(255)/
     *  .1219626E-1,.1220911E-1,.1222012E-1,
     *  .1222930E-1,.1223665E-1,.1224216E-1,
     *  .1224583E-1,.1224767E-1/
C                                  FIRST EXECUTABLE STATEMENT
      IER = 130
      IF (N.GT.20) GO TO 9005
      IER = 0
      DMLIN = 0.0
      BWT = 1.0
      DO 5 I=1,N
         BMA(I) = (B(I)-A(I))*0.5
         BPA(I) = (B(I)+A(I))*0.5
         BWT = BWT*BMA(I)
    5 CONTINUE
      NLIM = 8
C                                  USE 2**IRES POINT GAUSSIAN
C                                    FORMULA IN EACH DIMENSION
      DO 45 IRES=1,NLIM
         NN = 2**IRES
         IF ((2*NN)**N/(2**N-1).GT.MAXFCN) GO TO 50
         NN2 = NN/2
         DO 10 I=1,NN2
            NN2I = NN2-1+I
            NNI = NN+1-I
            TT(I) = -T(NN2I)
            TT(NNI) = -TT(I)
            WW(I) = W(NN2I)
            WW(NNI) = WW(I)
   10    CONTINUE
         DO 15 I=1,N
            M(I) = 1
   15    CONTINUE
         M(1) = 0
         TINT = 0.0
C                                  THE FOLLOWING IS EQUIVALENT
C                                    TO N NESTED DO-LOOPS
   20    M(1) = M(1)+1
         IF (N.EQ.1) GO TO 30
         DO 25 I=2,N
            IF (M(I-1).LE.NN) GO TO 25
            M(I-1) = 1
            M(I) = M(I)+1
   25    CONTINUE
   30    IF (M(N).GT.NN) GO TO 40
         WT = BWT
         DO 35 I=1,N
            J = M(I)
            WT = WT*WW(J)
            X(I) = BPA(I)+TT(J)*BMA(I)
   35    CONTINUE
         TINT = TINT+F(N,X)*WT
         GO TO 20
   40    RESULT(IRES) = TINT
         DMLIN = TINT
         IF (IRES.EQ.1) GO TO 45
C                                  COMPARE LAST TWO ESTIMATES
C                                    TO CHECK FOR CONVERGENCE
         IF (ABS(RESULT(IRES)-RESULT(IRES-1)).LE.
     *   AMAX1(AERR,RERR*ABS(RESULT(IRES)))) GO TO 9005
   45 CONTINUE
   50 IER = 129
 9005 IF (IER.GT.0) CALL UERTST(IER,6HDMLIN )
      RETURN
      END
