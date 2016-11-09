C   IMSL ROUTINE NAME   - ANESTU
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1980
C
C   PURPOSE             - ANALYSIS OF COMPLETELY NESTED DESIGN DATA
C                           WITH UNEQUAL NUMBERS IN THE SUBCLASSES
C
C   USAGE               - CALL ANESTU (NF,NL,Y,GM,S,NDF,EMS,IWK,IER)
C
C   ARGUMENTS    NF     - INPUT NUMBER OF FACTORS (NUMBER OF DISTINCT
C                           SUBSCRIPTS IN THE MODEL). NF MUST BE GREATER
C                           THAN OR EQUAL TO TWO.
C                NL     - INPUT VECTOR CONTAINING THE NUMBER OF LEVELS
C                           OF EACH FACTOR AT EACH LEVEL OF THE FACTOR
C                           IN WHICH IT IS NESTED. LET NLT BE THE NUMBER
C                           OF ELEMENTS REQUIRED FOR THE ABOVE AND NLL
C                           BE THE NUMBER OF ELEMENTS REQUIRED FOR THE
C                           LEVEL NUMBERS OF THE LOWEST LEVEL FACTOR.
C                           NL MUST BE OF LENGTH 2*NLT-NLL WHERE
C                           NL(1),...,NL(NLT) ARE THE LEVEL NUMBERS
C                           AND NL(NLT+1),...,NL(2*NLT-NLL) IS WORK
C                           STORAGE. EACH OF THE FIRST NLT ELEMENTS OF
C                           NL MUST BE GREATER THAN OR EQUAL TO ONE.
C                           SEE REMARKS FOR FURTHER DETAILS.
C                Y      - INPUT VECTOR CONTAINING THE RESPONSES.
C                           Y IS OF LENGTH EQUAL TO THE TOTAL NUMBER
C                           OF RESPONSE VALUES. ELEMENTS OF Y ARE
C                           ORDERED IN THE SORT SEQUENCE OF THE MODEL
C                           SUBSCRIPTS. FOR EXAMPLE, FOR THE MODEL
C                           Y(IJK) = M+A(I)+B(IJ)+C(IJK) Y MUST BE IN
C                           THE SORT SEQUENCE OF THE MODEL SUBSCRIPT
C                           SET, IJK. ON OUTPUT THE CONTENTS OF Y HAVE
C                           BEEN DESTROYED.
C                GM     - OUTPUT GRAND MEAN.
C                S      - OUTPUT VECTOR OF LENGTH NF+1 CONTAINING THE
C                           SUM OF SQUARES FOR THE HIGHEST LEVEL FACTOR,
C                           THE NEXT HIGHEST LEVEL, AND SO FORTH.
C                           LOCATION NF+1 CONTAINS THE CORRECTED TOTAL
C                           SUM OF SQUARES.
C                NDF    - OUTPUT VECTOR OF LENGTH NF+1 CONTAINING
C                           DEGREES OF FREEDOM CORRESPONDING TO
C                           COMPONENTS OF THE S VECTOR.
C                EMS    - OUTPUT VECTOR OF LENGTH NF(NF+1)/2 CONTAINING
C                           EXPECTED MEAN SQUARE COEFFICIENTS.
C                           SEE PROGRAMMING NOTES IN THE MANUAL DOCUMENT
C                           FOR FURTHER DETAILS ON THE CONTENTS OF EMS.
C                IWK    - WORK VECTOR OF LENGTH 5*NF.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER=129 INDICATES THAT NF WAS SPECIFIED LESS
C                             THAN TWO.
C                           IER=130 INDICATES THAT AN ELEMENT OF NL WAS
C                             SPECIFIED LESS THAN ONE.
C                           IER=131 INDICATES THAT AN ELEMENT OF NDF WAS
C                             ZERO, ALSO IMPLYING THAT NL WAS
C                             INCORRECTLY SPECIFIED.
C
C   REQD. IMSL ROUTINES - SINGLE/UERTST,UGETIO
C                       - DOUBLE/UERTST,UGETIO,VXADD,VXMUL,VXSTO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS      CONSIDER THE EXAMPLE MODEL Y(IJK) = M+A(I)+B(IJ)+C(IJK)
C                WHERE I=1,2,...,IMAX, J=1,2,...,JMAX(I), AND
C                K=1,2,...,KMAX(IJ). THEN NL CONTAINS THE LEVEL NUMBERS
C                AS FOLLOWS
C                NL = (IMAX,JMAX(1),JMAX(2),...,JMAX(IMAX),KMAX(1,1),
C                      KMAX(1,2),...,KMAX(1,JMAX(1)),KMAX(2,1),
C                      KMAX(2,2),...,KMAX(2,JMAX(2)),...,KMAX(IMAX,1),
C                      KMAX(IMAX,2),...,KMAX(IMAX,JMAX(IMAX)),WORK
C                      STORAGE)
C
C-----------------------------------------------------------------------
C
      SUBROUTINE ANESTU (NF,NL,Y,GM,S,NDF,EMS,IWK,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            NF,NL(1),NDF(1),IWK(1),IER
      REAL               Y(1),GM,S(1),EMS(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,IC,ICT,IC1,IC2,ID1,IE,IEN,IEND,IE1,II,III,
     1                   IJ,IK,IL,INF0,INF1,INF2,IPT,IS,ISK,IST,IS1,
     2                   IY1,IY2,I1,I2,I3,J,JM1,J2,K,KNT,K1,K2,L,LL,
     3                   M,M1,N,NFMM,NFM1,NFM1M,NFM3,NFP1,NF2,NF2M,
     4                   NF2P,NF3,NI,N1
      REAL               TEMP,GNL
      DOUBLE PRECISION   GM1,GM2,GM3,G
      DOUBLE PRECISION   GM4
C                                  FIRST EXECUTABLE STATEMENT
      IF(NF .LT. 2) GO TO 185
      NFP1 = NF+1
      IER = 0
      NF3 = NF*3
      NF2 = NF+NF
      NF2P = NF2+1
      NF2M = NF2-1
      NFM1 = NF-1
      NFM3 = NF-3
C                                  COMPUTE POINTERS TO GROUPS.
      IPT = 1
      J = NFP1
      IWK(1) = 0
      IWK(2) = 1
      K = 2
      IST = 1
      IEND = 1
    5 K = K+1
      DO 10 I=IST,IEND
         IF(NL(I) .LT. 1) GO TO 190
   10 IPT = IPT+NL(I)
      IWK(K) = IPT
      IF(K .GE. NFP1) GO TO 15
      IST = IEND+1
      IEND = IWK(K)
      GO TO 5
   15 J = J-1
      K = K+1
      IF(IWK(J) .EQ. 0) GO TO 20
      IPT = IPT+IWK(J)-IWK(J-1)
      IWK(K) = IPT
      GO TO 15
C                                  COMPUTE LATTER PART OF NL VECTOR.
   20 M = NF+2
      ID1 = IWK(NFP1)
      J = NFM1
      IS1 = IWK(NF)+1
   25 IS = IWK(J)+1
      IEN = IWK(J+1)
      DO 35 I=IS,IEN
         IE1 = IS1+NL(I)-1
         ID1 = ID1+1
         NL(ID1) = 0
         DO 30 II=IS1,IE1
            NL(ID1) = NL(ID1)+NL(II)
   30    CONTINUE
         IS1 = IS1+NL(I)
   35 CONTINUE
      M = M+1
      J = J-1
      IF(J .NE. 0) GO TO 25
C                                  COMPUTE THE GRAND MEAN
      DO 40 I=1,NF
         IWK(M) = IWK(I+1)-IWK(I)
   40 M = M+1
      GM1 = 0.0D0
      J = IWK(NF2)
      J = NL(J)
      DO 45 I=1,J
   45 GM1 = GM1+DBLE(Y(I))
      GM = GM1/J
C                                  COMPUTE CORRECTED TOTAL SUM OF SQ.
      GM1 = 0.0D0
      DO 50 I=1,J
         GM2 = Y(I)-GM
   50 GM1 = GM1+(GM2*GM2)
      S(NFP1) = GM1
C                                  REPLACE ELEMENTS OF Y VECTOR WITH
C                                  APPROPRIATE TOTALS.
   55 IC = 1
      ICT = 0
      IC1 = 1
      DO 70 I=1,NF
         IJ = NF-ICT
         I1 = IWK(IJ)+1
         I2 = IWK(IJ+1)
         DO 65 II=I1,I2
            GM1 = 0.0D0
            I3 = NL(II)
            IC2 = IC1+I3-1
            DO 60 III=IC1,IC2
   60       GM1 = GM1+DBLE(Y(III))
            Y(IC) = GM1
            IC = IC+1
   65    IC1 = IC2+1
         ICT = ICT+1
         M = M-1
   70 IC1 = IC-IWK(M)
      IC = IC-1
C                                  COMPUTE FIRST NF ELEMENTS OF S VECTOR
      J2 = 1
      K2 = 1
      DO 85 I=1,NFM1
         INF1 = IWK(NF2M-I)
         INF2 = IWK(NF2-I)
         IY1 = IWK(NFP1)-IWK(I+2)
         IY2 = IWK(NFP1)-IWK(I+1)
         GM1 = 0.0D0
   75    IL = NL(K2)
         INF2 = INF2+1
         IY2 = IY2+1
         GM3 = 1.D0/NL(INF2)
         DO 80 J=1,IL
            INF1 = INF1+1
            IY1 = IY1+1
            GM2 = NL(INF1)
            GM4 = (DBLE(Y(IY1))/GM2) - (DBLE(Y(IY2))*GM3)
   80    GM1 = GM1 + (GM2*GM4*GM4)
         K2 = K2+1
         J2 = J2+IL
         IF(J2 - IWK(I+2)) 75,85,85
   85 S(I) = GM1
      S(NF) = S(NFP1)
      DO 90 I=1,NFM1
   90 S(NF) = S(NF)-S(I)
C                                  COMPUTE ELEMENTS OF NDF
      INF1 = NF+NFP1
      DO 95 I=1,NFM1
         NDF(I) = IWK(INF1+1)-IWK(INF1)
         IF(NDF(I) .EQ. 0) GO TO 195
   95 INF1 = INF1+1
      J = IWK(NF2)
      NDF(NFP1) = NL(J)-1
      NDF(NF) = NDF(NFP1)
      DO 100 I=1,NFM1
  100 NDF(NF) = NDF(NF)-NDF(I)
      IF(NDF(NF) .EQ. 0 .OR. NDF(NFP1) .EQ. 0) GO TO 195
C                                  COMPUTE EMS VECTOR.
      DO 105 IJ=1,NF2
  105 IWK(NF3+IJ) = IWK(IJ)
      KNT = 0
      IK = 0
      N = IWK(NF2)
      GNL = NL(N)
  110 KNT = KNT+1
      NFMM = NFM3-IK
      M1 = 0
      INF0 = IWK(NF3-IK)
      INF1 = IWK(NF+IK)
      DO 115 I=1,INF0
         INF1 = INF1+1
  115 M1 = M1+(NL(INF1)*NL(INF1))
      EMS(KNT) = M1/GNL
      IF(IK .GT. NFM3) GO TO 155
      K = 2
      IE = NF2+2
      LL = NF+IK
      NFM1M = NFM1-IK
      DO 150 I1 = 2,NFM1M
         KNT = KNT+1
         N = IWK(K)
         K = K+1
         IEND = IWK(IE)
         NI = IWK(NF2P-K) + 1
         GM1 = 0.0D0
         DO 140 I=1,IEND
            J = NL(N+I)
            ICT = 0
  120       M = 0
            IF(ICT .GE. NFMM) GO TO 130
            L = IWK(K+ICT)+1
            IWK(K+ICT) = IWK(K+ICT)+J
            DO 125 II=1,J
               M = M+NL(L)
  125       L = L+1
            ICT = ICT+1
            J = M
            GO TO 120
  130       L = IWK(LL)+1
            IWK(LL) = IWK(LL)+J
            DO 135 II=1,J
               M = M+(NL(L)*NL(L))
  135       L = L+1
            G = M
            GM1 = GM1+(G/NL(NI))
  140    NI = NI+1
         DO 145 IJ=1,NF2
  145    IWK(IJ) = IWK(NF3+IJ)
         IE = IE+1
         NFMM = NFMM-1
  150 EMS(KNT) = GM1
      IK = IK+1
      KNT = KNT+1
      EMS(KNT) = GNL
      GO TO 110
  155 J = (NF*(NF+1))/2
      EMS(KNT+1) = GNL
      JM1 = J-2
      K1 = 2
      K2 = 2
      DO 165 I=1,JM1
         J = J-1
         K1 = K1-1
         IF(K1 .GT. 0) GO TO 160
         K2 = K2+1
         K1 = K2
  160 EMS(J) = (EMS(J)-EMS(J-1))/NDF(K1)
  165 CONTINUE
      ISK = NF
      J = 1
      DO 170 I=1,NF
         EMS(J) = 1.0
         J = J+ISK
  170 ISK = ISK-1
      K = 3
      N1 = NF
      J = NF-2
      L = K+NFM1
  175 IF(J .LT. 1) GO TO 9005
      DO 180 I=1,J
         TEMP = EMS(K)
         EMS(K) = EMS(L)
         EMS(L) = TEMP
         K = K+1
  180 L = L+N1-I
      N1 = N1-1
      J = J-2
      L = K+NFP1
      K = K+3
      GO TO 175
  185 IER = 129
      GO TO 9000
  190 IER = 130
      GO TO 9000
  195 IER = 131
 9000 CONTINUE
      CALL UERTST(IER,'ANESTU')
 9005 RETURN
      END
