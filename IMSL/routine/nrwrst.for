C   IMSL ROUTINE NAME   - NRWRST
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1981
C
C   PURPOSE             - WILCOXONS RANK-SUM TEST
C                           (MANN-WHITNEY TEST)
C
C   USAGE               - CALL NRWRST (X,M,N,EPS,IR,STAT,IER)
C
C   ARGUMENTS    X      - INPUT VECTOR OF LENGTH M+N
C                           CONTAINING THE SAMPLES FROM EACH OF TWO
C                           POPULATIONS. SAMPLE ONE IS CONTAINED IN THE
C                           ELEMENTS X(I), I=1,...,M AND SAMPLE TWO IS
C                           CONTAINED IN THE ELEMENTS X(J),
C                           J=M+1,...,M+N.
C                           ON OUTPUT, X CONTAINS THE SORTED COMBINED
C                           SAMPLE.
C                M      - INPUT SIZE OF SAMPLE ONE.
C                N      - INPUT SIZE OF SAMPLE TWO.
C                         NOTE THAT THE MINIMUM OF M AND N MUST BE
C                           GREATER THAN OR EQUAL TO ONE AND
C                           M+N MUST BE GREATER THAN OR EQUAL TO THREE.
C                EPS    - INPUT VALUE TO BE USED FOR DETECTING TIES
C                           IN THE COMBINED SAMPLES.
C                           IF EPS IS LESS THAN 10**(-5), THE EPSILON
C                           USED WILL BE 10**(-5). IF THE ELEMENTS OF
C                           THE SAMPLE ARE SUFFICIENTLY LARGE AS TO
C                           JUSTIFY THE NEED FOR A LARGER EPSILON, THEN
C                           EPS*X(I) IS THE EPSILON USED.
C                           SEE THE PROGRAMMING NOTES IN THE MANUAL
C                           DOCUMENT FOR FURTHER DETAILS.
C                IR     - WORK VECTOR OF LENGTH M+N.  ON OUTPUT, IR
C                           WILL CONTAIN THE ORIGINAL POSITIONS IN THE
C                           INPUT X VECTOR OF THE CORRESPONDING ELEMENTS
C                           IN THE SORTED X VECTOR, WITH THE SIGNS FOR
C                           THOSE DENOTING THE LARGER SAMPLE NEGATED.
C                STAT   - OUTPUT VECTOR OF LENGTH 6.
C                         STAT(1) CONTAINS W, THE RANK SUM STATISTIC,
C                           WITH TIES HANDLED SUCH THAT W IS AS SMALL
C                           AS IS POSSIBLE.
C                         STAT(2) CONTAINS 2*WBAR-W, WHERE WBAR IS THE
C                           EXPECTED VALUE OF W.
C                         STAT(3) CONTAINS THE PROBABILITY OF OBTAINING
C                           A STATISTIC LESS THAN OR EQUAL TO THE
C                           MINIMUM OF (W, 2*WBAR-W), IF THE POPULATIONS
C                           ARE IDENTICAL.
C                         STAT(4) CONTAINS THE SAME STATISTIC AS IN
C                           STAT(1) WITH THE EXCEPTION THAT TIES ARE
C                           HANDLED SUCH THAT W IS AS LARGE AS IS
C                           POSSIBLE.
C                         STAT(5) CONTAINS THE SAME STATISTIC AS IN
C                           STAT(2) WITH THE EXCEPTION THAT TIES ARE
C                           HANDLED SUCH THAT W IS AS LARGE AS IS
C                           POSSIBLE.
C                         STAT(6) CONTAINS THE SAME STATISTIC AS IN
C                           STAT(3) WITH THE EXCEPTION THAT TIES ARE
C                           HANDLED SUCH THAT W IS AS LARGE AS IS
C                           POSSIBLE.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER=129 INDICATES THAT THE MINIMUM OF M AND
C                             N (THE SAMPLE SIZES) IS LESS THAN 1
C                             OR M+N IS LESS THAN 3.
C                         WARNING ERROR
C                           IER=34 INDICATES THAT TIED OBSERVATIONS
C                             OCCURRED BETWEEN SAMPLES.  THIS ERROR
C                             OVERRIDES WARNING ERROR 35.
C                           IER=35 INDICATES THAT M AND N ARE LESS THAN
C                             25.
C
C   REQD. IMSL ROUTINES - MERRC=ERFC,UERTST,UGETIO,VSRTR
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS      THE MANN-WHITNEY U STATISTIC IS U=W-K*(K+1)/2 WHERE
C                K IS THE SMALLEST SAMPLE SIZE (MINIMUM OF M AND N)
C                AND W = STAT(1) OR STAT(4). TABLES OF CRITICAL VALUES
C                OF W OR 2*WBAR-W ARE AVAILABLE IN THE REFERENCES GIVEN
C                IN THE MANUAL DOCUMENT.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE NRWRST (X,M,N,EPS,IR,STAT,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            M,N,IR(1),IER
      REAL               X(1),EPS,STAT(6)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,IFLAG,IPLACE,ISTART,ISTOP,ISTRT,I32,JSTART,
     1                   JSTOP,KTB,LTB,MM,NSL,NSU,NTOD,NTOT
      REAL               EP1,EP11,EW,EWC,EWC1,EWC2,ONED12,RSQ1H,R1D2P,
     1                   SPI,VW,XMM,XNN,XTOT,ZI2,W(2),Z(2),EWC3(2),
     2                   EWC4(2)
      DATA               SPI/3.141593/
      DATA               ONED12/.8333333E-1/
      DATA               RSQ1H/.7071068/
C                                  FIRST EXECUTABLE STATEMENT
      R1D2P = 1.0/SQRT(2*SPI)
C                                  CHECK VALUES OF M AND N
      NTOT = M + N
      IF ((MIN0(M,N) .GE. 1) .AND. (NTOT .GE. 3)) GO TO 5
      IER = 129
      GO TO 9000
C                                  CHECK RANGE OF SEED
    5 IER = 0
      IF ((M .LT. 25) .AND. (N .LT. 25)) IER = 35
C                                  ASSIGN REAL POSITIONS FOR BOTH SAMP.
      DO 10 I=1,NTOT
         IR(I) = I
   10 CONTINUE
C                                  DETERMINE ENDPOINTS FOR NEGATIVE AND
C                                  POSITIVE INDICATIONS OF POSITION IN X
      IF (M .GE. N) GO TO 15
      ISTRT = M + 1
      ISTOP = M + N
      MM = M
      GO TO 20
   15 ISTRT = 1
      ISTOP = M
      MM = N
C                                  NEGATE IR FOR LARGER SAMPLE
   20 DO 25 I=ISTRT,ISTOP
         IR(I) = -I
   25 CONTINUE
      NTOD = NTOT - 1
C                                  SORT X VECTOR BY INVOKING VSRTR
      CALL VSRTR (X,NTOT,IR)
C                                  CHECK VALUE OF EPS
      EP1 = EPS
      IF (EPS .LT. 1.E-5) EP1 = 1.E-5
C                                  SEARCH FOR TIES
      IFLAG = 0
      IPLACE = 1
      NSL = 0
      NSU = 0
   30 ISTART = IPLACE
      JSTART = 0
      LTB = 0
      KTB = 0
      DO 55 I=ISTART,NTOD
         EP11 = EP1
         IF(X(I+1) .EQ. (X(I)+EP1)) EP11 = EP1*X(I)
         IF(X(I+1) .GT. X(I)+EP11) GO TO 35
         IER=34
         IF(JSTART .EQ. 0) JSTART = I
         IF(IR(I) .GT. 0) LTB = LTB+1
         KTB = KTB+1
         IF(I .NE. NTOD) GO TO 55
         IF(IR(I+1) .GT. 0) LTB = LTB+1
         IF(LTB .EQ. 0) GO TO 75
         KTB = KTB+1
         IFLAG = 1
         GO TO 60
   35    IF(I .EQ. 1) GO TO 45
         EP11 = EP1
         IF(X(I) .EQ. X(I-1)+EP1) EP11 = EP1*X(I-1)
         IF(X(I) .GT. X(I-1)+EP11) GO TO 45
         IF(IR(I) .GT. 0) LTB = LTB+1
         IF(LTB .NE. 0) GO TO 40
         JSTART = 0
         KTB = 0
         GO TO 50
   40    IPLACE = I+1
         KTB = KTB+1
         GO TO 60
   45    IF(IR(I) .LT. 0) GO TO 50
         NSL = NSL+I
         NSU = NSU+I
   50    IF(I .NE. NTOD) GO TO 55
         IF(IR(I+1) .LT. 0) GO TO 75
         NSL = NSL+I+1
         NSU = NSU+I+1
         GO TO 75
   55 CONTINUE
   60 JSTOP = JSTART+LTB-1
      DO 65 I = JSTART,JSTOP
         NSL = NSL+I
   65 CONTINUE
      ISTART = KTB-LTB+JSTART
      ISTOP = JSTART+KTB-1
      DO 70 I = ISTART,ISTOP
         NSU = NSU+I
   70 CONTINUE
      IF(IFLAG .EQ. 1) GO TO 75
      IF(IPLACE .LT. NTOT) GO TO 30
      IF(IR(IPLACE) .LT. 0) GO TO 75
      NSL = NSL+IPLACE
      NSU = NSU + IPLACE
C                                  CALCULATE W AND 2*WBAR-W
   75 XNN = NTOT - MM
      XMM = MM
      XTOT = NTOT
      EWC = XMM * (XTOT + 1.0)
      EWC1 = XNN * EWC
      STAT(1) = NSL
      STAT(4) = NSU
      STAT(2) = EWC - STAT(1)
      STAT(5) = EWC - STAT(4)
      EW = .5 * EWC
      VW = 1./SQRT(ONED12 * EWC1)
      EWC2 = (XMM **2 + XNN**2 + XMM * XNN + XTOT) * 0.05/EWC1
C                                  CALCULATE REJECTION PROBABILITIES
      I32 = 1
      DO 80 I = 1,2
         W(I) = AMIN1(STAT(I32), STAT(I32+1))
         Z(I) = W(I) - EW - 0.5
         IF ((W(I) - EW) .LT. 0.0) Z(I) = Z(I) + 1.0
         Z(I) = Z(I)*VW
         ZI2 = Z(I) * Z(I)
         EWC3(I) = .5*ERFC(-RSQ1H*Z(I))
         EWC4(I) = -R1D2P * (Z(I)*(ZI2 - 3.0))* EXP(-0.5*ZI2)
         STAT(I32+2) = EWC3(I) - EWC4(I) * EWC2
         I32 = I32+3
   80 CONTINUE
      IF (IER .EQ. 0) GO TO 9005
 9000 CONTINUE
      CALL UERTST(IER,6HNRWRST)
 9005 RETURN
      END
