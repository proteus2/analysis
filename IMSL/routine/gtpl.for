C   IMSL ROUTINE NAME   - GTPL
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - POKER TEST TALLY OF HAND TYPES AND STATISTICS
C
C   USAGE               - CALL GTPL (R,M,I1,I2,IQ,K2,K1,CSOBS,HSAVE,HE,
C                           HEH,H,IER)
C
C   ARGUMENTS    R      - INPUT RANDOM NUMBER VECTOR, M LONG, TO BE
C                           SUBJECTED TO THE POKER TEST
C                M      - INPUT LENGTH OF R
C                I1     - INPUT FIRST BIT IN R TO BE CONSIDERED. I1
C                           MUST BE GREATER THAN ZERO.
C                I2     - INPUT LAST BIT IN R TO BE CONSIDERED. I2
C                           MUST BE NO GREATER THAN THE NUMBER OF BITS
C                           IN A REAL WORD.
C                IQ     - INPUT NUMBER OF HANDS, OF THE M, WHICH ARE TO
C                           BE TALLIED BEFORE A CHI-SQUARE STATISTIC IS
C                           CALCULATED
C                K2     - INPUT NUMBER OF EQUIPROBABLE SUBDIVISIONS
C                           OF THE CHI-SQUARE DENSITY INTO WHICH THE
C                           STATISTICS CHI-SQUARE ARE TO BE TALLIED.
C                           THESE STATISTICS ARE CALCULATED FROM THE
C                           TALLIES OF IQ HANDS INTO K1 CATEGORIES.
C                K1     - OUTPUT, FOR USE IN GTPOK.
C                           K1= THE GREATEST INTEGER IN 1+(I2-I1+1)/2
C                CSOBS  - OUTPUT TALLYING VECTOR OF LENGTH K2.
C                           CSOBS MUST BE NULL ON THE FIRST ENTRY.
C                HSAVE  - OUTPUT VECTOR OF LENGTH K1 IN WHICH ALL
C                           TALLIES OF ALL HANDS ARE SAVED.  ON THE
C                           FIRST CALL HSAVE MUST BE A NULL VECTOR.
C                HE     - OUTPUT VECTOR OF LENGTH K1 CONTAINING THE
C                           EXPECTED VALUES OF PARTICULAR HAND TYPES.
C                HEH    - OUTPUT VECTOR SET TO HE*IQ OF LENGTH K1.
C                           HEH CONTAINS THE EXPECTED VALUES OF THE
C                           TALLIES OF THE IQ HANDS.
C                H      - WORK VECTOR OF LENGTH K1 USED TO HOLD EACH
C                           HAND TALLY FOR THE IQ HANDS BEFORE BEING
C                           MOVED TO HSAVE.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER = 129 INDICATES AN ERROR OCCURED IN
C                             GTPKP OR MDCH
C                           IER = 130 INDICATES IQ NOT A FACTOR OF M
C                           IER = 131 INDICATES I1 AND/OR I2 INCORRECTLY
C                             SPECIFIED.
C                         WARNING ERROR
C                           IER = 36 INDICATES AN EXPECTED VALUE LESS
C                             THAN FIVE WAS CALCULATED INTO HEH.
C                           IER = 37 INDICATES AN EXPECTED VALUE LESS
C                             THAN ONE WAS CALCULATED INTO HEH.  THE
C                             CHI-SQUARE TEST SHOULD BE SUSPECT.
C
C   REQD. IMSL ROUTINES - H32/GTPBC,GTPKP,MDCH,MDNOR,MERRC=ERFC,
C                           MGAMAD=DGAMMA,UERTST,UGETIO
C                       - H36,H48,H60/GTPBC,GTPKP,MDCH,MDNOR,
C                           MERRC=ERFC,MGAMA=GAMMA,UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE GTPL   (R,M,I1,I2,IQ,K2,K1,CSOBS,HSAVE,HE,HEH,H,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            M,I1,I2,IQ,K2,K1,IER
      REAL               R(1),CSOBS(1),HSAVE(1),HE(1),HEH(1),H(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            KT,I,IT,IEND,IST,J,JJ,KZ,II,KK,IET,JCT,J1,IBT
      REAL               TEMP1,HTEMP,XSQRD,H2M1,P1
      DATA               IBT/32/
C                                  FIRST EXECUTABLE STATEMENT
      IF(I1.LT.1 .OR. I1.GT.IBT .OR. I2.LE.I1 .OR. I2.GT.IBT) GO TO 50
      KT = I2-I1+1
      K1 = KT/2+1
C                                  OBTAIN PROBABILITIES OF HAND TYPES
      CALL GTPKP(KT,HE,IER)
      IF(IER.GT.127) GO TO 9000
C                                  OBTAIN EXPECTED NUMBER OF HAND TYPES
      DO 5 I=1,K1
         HEH(I) = HE(I)*IQ
         IF(HEH(I).LT.5.0) IER = 36
         IF(HEH(I).LT.1.0) IER = 37
    5 CONTINUE
      IT = M/IQ
      IF(IT*IQ .NE. M) GO TO 45
      IEND = IQ
      IST = 1
      DO 40 J=1,IT
C                                  ZERO TO HAND TYPE COUNTS
         DO 10 JJ=1,K1
   10    H(JJ) = 0.0
C                                  FIND ONE HAND TYPE
         DO 25 I=IST,IEND
            CALL GTPBC(I1,I2,R(I),KZ)
            DO 20 II=1,K1
               KK = II-1
               IF(KZ.EQ.KK) GO TO 15
               IF(KZ.NE.KT-KK) GO TO 20
   15          H(II) = H(II) + 1.
               GO TO 25
   20       CONTINUE
   25    CONTINUE
C                                  PREPARE FOR NEXT IQ HANDS
         IST = IST+IQ
         IEND = IEND+IQ
         TEMP1 = 0.0
C                                  CALCULATE ONE CHI-SQUARED STATISTIC
         DO 30 II=1,K1
            HTEMP = H(II)-HEH(II)
   30    TEMP1 = TEMP1+(HTEMP*HTEMP)/HEH(II)
         XSQRD = TEMP1
         H2M1 = K1-1
C                                  TALLY THE CHI-SQUARED STATISTIC
      CALL MDCH(XSQRD,H2M1,P1,IET)
      IF(IET.GT.127) GO TO 55
         JCT = 1+K2*P1
         IF(JCT.GT.K2) JCT = K2
         CSOBS(JCT) = CSOBS(JCT) + 1.0
C                                  RETAIN COUNTS OF HAND TYPES
         DO 35 J1=1,K1
   35    HSAVE(J1) = H(J1)+HSAVE(J1)
   40 CONTINUE
      IF(IER) 9000,9005,9000
   45 IER = 130
      GO TO 9000
   50 IER = 131
      GO TO 9000
   55 IER = IET
 9000 CONTINUE
      CALL UERTST(IER,'GTPL  ')
 9005 RETURN
      END
