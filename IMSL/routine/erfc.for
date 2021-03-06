        FUNCTION ERFC(X)
        DIMENSION P1(4),Q1(4),P2(8),Q2(8),P3(5),Q3(5)
        LOGICAL NEGA
        DATA P1/242.667955203,21.979261628,6.996383489,-0.0356098437/
        DATA Q1/215.05887587,91.164905449,15.0827976304,1./
        DATA P2/300.45926102,451.918953712,339.208167343,
     &          152.989285047,43.162227222,7.21175825089,
     &          .5641955174789,-1.3686485738E-7/
        DATA Q2/3.004592609569E2,7.90950925328E2,9.3135409485E2,
     &          6.389802644656E2,2.77585444743E2,7.7000152935E1,
     &          12.782727319629,1./
        DATA P3/-2.996107077E-3,-4.947309106E-2,-2.2695659353E-1,
     &          -2.786613086E-1,-2.231924597E-2/
        DATA Q3/1.0620923053E-2,1.9130892611E-1,1.05167510707,
     &          1.98733201817,1./
        IF(X.LT.0) THEN
          NEGA=.TRUE.
          XTEMP=-X
        ELSE
          NEGA=.FALSE.
          XTEMP=X
        END IF
        PINV=0.56419
        IF (XTEMP.GE.0.46875) GO TO 100
        IF (XTEMP.EQ.0) THEN
            ERFC=1.
            GO TO 999
        END IF
        UP=0.
        DOWN=0.
        DO 10 I=1,4
        IM1=2*(I-1)
        UP=UP+P1(I)*(XTEMP**IM1)
        DOWN=DOWN+Q1(I)*(XTEMP**IM1)
10      CONTINUE
        ERFC=1.-XTEMP*UP/DOWN
        GO TO 999
100     IF(XTEMP.GE.4) GO TO 200
        UP=0.
        DOWN=0.
        DO 20 I=1,8
        IM1=I-1.
        UP=UP+P2(I)*(XTEMP**IM1)
        DOWN=DOWN+Q2(I)*(XTEMP**IM1)
20      CONTINUE
        ERFC=EXP(-XTEMP*XTEMP)*UP/DOWN
        GO TO 999
200     IF(XTEMP.GE.9) THEN
        ERFC=0.0
        GO TO 999
        END IF
        UP=0.
        DOWN=0.
        DO 30 I=1,5
        IM1=-2.*(I-1)
        UP=UP+P3(I)*(XTEMP**IM1)
        DOWN=DOWN+Q3(I)*(XTEMP**IM1)
30      CONTINUE
        ERFC=EXP(-XTEMP*XTEMP)/XTEMP*(PINV+UP/DOWN/XTEMP/XTEMP)
999     IF(NEGA) THEN
          ERFC=2.-ERFC
        END IF
        RETURN
        END
