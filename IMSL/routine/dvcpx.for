C   IMSL ROUTINE NAME   - DVCPX
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1982
C
C   PURPOSE             - NUCLEUS CALLED ONLY BY IMSL SUBROUTINE DVCPR
C
C   REQD. IMSL ROUTINES - NONE
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE DVCPX (A,C,DEL,M,N,P,R,IR,IC,SING,ICA,AUX,MTNMAX,MMAX,
     *NMAX,MMAX2)
C                                  SPECIFICATIONS FOR ARGUMENTS
C
      INTEGER            M,N,P,R,MTNMAX,MMAX,NMAX,MMAX2,IR(NMAX,MMAX),
     *                   IC(NMAX,MMAX),ICA(MMAX)
      REAL               A(MTNMAX,MMAX),C(MTNMAX,MMAX),DEL(MMAX,MTNMAX),
     *                   AUX(MMAX2,MMAX)
      LOGICAL            SING
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
C
      INTEGER            I01,I0,I1S,I1,I2,I3,I4,I7,I8,II1,II2,II,INI,
     *                   IPIV,IPP,IP,IS1,ISM,ISP,IS,ITE,IXI,IXM,IXP,IX,
     *                   I,J1,JP,JSX,JS,JX,J,KS,M1,MP,N1,P1
      REAL               AK,TE,XMU,Y
C                                  FIRST EXECUTABLE STATEMENT
      SING = .FALSE.
      P1 = P+1
      MP = M+P
      N1 = N-1
C                                  MAIN LOOP
      DO 275 II=1,N1
         IX = (II-1)*M
         IXM = IX+M
         IXI = IX-M
         DO 10 I=1,M
            IC(II,I) = I
            IR(II,I) = I
            DO 5 J=1,M
    5       AUX(I,J) = A(IX+I,J)
   10    CONTINUE
         IF (P) 65, 65, 15
   15    DO 25 I=1,P
            DO 20 J=1,M
   20       AUX(I+M,J) = C(IX+I,J)
   25    CONTINUE
C                                  REDUCTION OF FIRST P ROWS OF A1(II)
         DO 60 I=1,P
            TE = 0.
            DO 30 J=I,M
               Y = ABS(AUX(I,J))
               IF (Y.LE.TE) GO TO 30
               TE = Y
               IPIV = J
   30       CONTINUE
            IF (TE.EQ.0) GO TO  335
C                                  COLUMN INTERCHANGES
            IF (IPIV.EQ.I) GO TO 45
            ITE = IC(II,I)
            IC(II,I) = IC(II,IPIV)
            IC(II,IPIV) = ITE
            DO 35 IS=1,MP
               TE = AUX(IS,I)
               AUX(IS,I) = AUX(IS,IPIV)
               AUX(IS,IPIV) = TE
   35       CONTINUE
            IF (II.EQ.1) GO TO 45
            I1 = IXI+P1
            I2 = IXI+M
            DO 40 IS=I1,I2
               TE = C(IS,I)
               C(IS,I) = C(IS,IPIV)
               C(IS,IPIV) = TE
   40       CONTINUE
C                                  FACTORIZATION
   45       AK = 1./AUX(I,I)
            I1 = I+1
            DO 55 IS=I1,MP
               XMU = AUX(IS,I)*AK
               AUX(IS,I) = XMU
               DO 50 JS=I1,M
   50          AUX(IS,JS) = AUX(IS,JS)-XMU*AUX(I,JS)
   55       CONTINUE
   60    CONTINUE
C                                  REDUCTION OF REMAINING PART OF
C                                    A1(II)
   65    DO 110 J=P1,M
            TE = 0.
            DO 70 I=J,MP
               Y = ABS(AUX(I,J))
               IF (Y.LE.TE) GO TO 70
               TE = Y
               IPP = I
   70       CONTINUE
            IF (TE.EQ.0.) GO TO  335
C                                  ROW INTERCHANGES
            IF (IPP.EQ.J) GO TO 95
            IPIV = IPP-P
            JP = J-P
            ITE = IR(II,JP)
            IR(II,JP) = IR(II,IPIV)
            IR(II,IPIV) = ITE
            DO 75 JS=1,M
               TE = AUX(J,JS)
               AUX(J,JS) = AUX(IPP,JS)
               AUX(IPP,JS) = TE
   75       CONTINUE
            II1 = IX+J
            II2 = IX+IPP
            IF (IPP.GT.M) GO TO 85
            DO 80 JS=1,M
               TE = C(II1,JS)
               C(II1,JS) = C(II2,JS)
               C(II2,JS) = TE
   80       CONTINUE
            GO TO 95
   85       DO 90 JS=1,M
               TE = C(II1,JS)
               C(II1,JS) = A(II2,JS)
               A(II2,JS) = TE
   90       CONTINUE
C                                  FACTORIZATION
   95       IF (J.EQ.M) GO TO 110
            AK = 1./AUX(J,J)
            J1 = J+1
            DO 105 IS=J1,MP
               XMU = AUX(IS,J)*AK
               AUX(IS,J) = XMU
               DO 100 JS=J1,M
  100          AUX(IS,JS) = AUX(IS,JS)-XMU*AUX(J,JS)
  105       CONTINUE
  110    CONTINUE
C                                  RESTORING B1(II+1) AND AL(II)
         IF (P.EQ.0) GO TO 125
         I1 = M-P+1
         DO 120 IS=I1,M
            IP = IR(II,IS)
            IF (IP.EQ.IS) GO TO 120
            ISP = IS+P+IXI
            IPP = IP+P+IX
            DO 115 JS=1,M
  115       C(ISP,JS) = A(IPP,JS)
  120    CONTINUE
  125    DO 135 I=1,M
            DO 130 J=1,M
  130       A(IX+I,J) = AUX(I,J)
  135    CONTINUE
C
         IS = IX+1
         ISP = IX+P
         DO 140 JS=1,M
  140    ICA(JS) = IC(II,JS)
         JS = 1
  145    IF (JS.GE.M) GO TO 180
         IP = ICA(JS)
         IF (IP.NE.JS) GO TO 150
         JS = JS+1
         GO TO 145
C
  150    INI = JS
         ICA(INI) = INI
C                                  COLUMN INTERCHANGES FOR B1(II+1)
         IF (P.EQ.0) GO TO 165
  155    DO 160 I1=IS,ISP
            TE = C(I1,INI)
            C(I1,INI) = C(I1,IP)
            C(I1,IP) = TE
  160    CONTINUE
  165    IF (R.EQ.0) GO TO 175
C                                  COLUMN INTERCHANGES OF D(II)
         IPP = IX+IP
         JX = IX+INI
         DO 170 IS1=1,R
            TE = DEL(IS1,IPP)
            DEL(IS1,IPP) = DEL(IS1,JX)
            DEL(IS1,JX) = TE
  170    CONTINUE
  175    INI = IP
         IP = ICA(IP)
         ICA(INI) = INI
         IF (IP.NE.JS) GO TO 155
         JS = JS+1
         GO TO 145
C
  180    IF (R.EQ.0) GO TO 230
C                                  SOLUTION OF DE(II)*AL(II)=-
C                                    DELTA(II)
         I1 = IX+1
         DO 210 IS=1,R
            DO 195 JS=I1,IXM
               TE = DEL(IS,JS)
               JSX = JS-IX
               IF (JS.EQ.I1) GO TO 190
               J1 = JS-1
               DO 185 KS=I1,J1
  185          TE = TE-DEL(IS,KS)*A(KS,JSX)
  190          DEL(IS,JS) = TE/A(JS,JSX)
  195       CONTINUE
            M1 = M-1
            DO 205 JS=1,M1
               J = IXM-JS
               TE = DEL(IS,J)
               J1 = J+1
               DO 200 KS=J1,IXM
  200          TE = TE-DEL(IS,KS)*A(KS,J-IX)
               DEL(IS,J) = TE
  205       CONTINUE
  210    CONTINUE
C                                  DE(II+1)=-DE(II)*GA(II)
         IXP = IX+P1
         DO 225 IS=1,R
            DO 220 JS=1,M
               TE = 0.
               DO 215 KS=IXP,IXM
  215          TE = TE-DEL(IS,KS)*C(KS,JS)
               DEL(IS,IXM+JS) = TE
  220       CONTINUE
  225    CONTINUE
C                                  SOLUTION OF BE(II+1)*AL(II)=B(II+1)
  230    I1 = IX+1
         IF (P.EQ.0) GO TO 275
         I2 = IX+P
         DO 270 IS=I1,I2
            DO 245 JS=1,M
               TE = C(IS,JS)
               IF (JS.EQ.1) GO TO 240
               J1 = JS-1
               DO 235 KS=1,J1
  235          TE = TE-C(IS,KS)*A(IX+KS,JS)
  240          C(IS,JS) = TE/A(IX+JS,JS)
  245       CONTINUE
            M1 = M-1
            DO 255 JS=1,M1
               J = M-JS
               TE = C(IS,J)
               J1 = J+1
               DO 250 KS=J1,M
  250          TE = TE-C(IS,KS)*A(IX+KS,J)
               C(IS,J) = TE
  255       CONTINUE
C                                  AL(II+1)=A(II+1)-BE(II+1)*GA(II)
            ISM = IS+M
            DO 265 JS=1,M
               TE = A(ISM,JS)
               DO 260 KS=P1,M
  260          TE = TE-C(IS,KS)*C(IX+KS,JS)
               A(ISM,JS) = TE
  265       CONTINUE
  270    CONTINUE
  275 CONTINUE
C                                  COMPLETE COMPUTATION OF ALFA(N)
      I2 = N1*M
      I1 = I2+P
      IF (R.EQ.0) GO TO 290
      DO 285 IS=1,R
         I1S = I1+IS
         DO 280 JS=1,M
  280    A(I1S,JS) = A(I1S,JS)+DEL(IS,I2+JS)
  285 CONTINUE
C                                  L U DECOMPOSITION OF ALFA(N) (COLUMN
C                                    PIVOTING)
  290 I4 = I1-M+1
      I1 = I2+1
      I3 = N*M-1
      I8 = I3+1
      DO 295 I=1,M
  295 IC(N,I) = I
      DO 330 I=I1,I3
         I0 = I-I2
         TE = 0.
         DO 300 J=I0,M
            Y = ABS(A(I,J))
            IF (Y.LE.TE) GO TO 300
            TE = Y
            IPIV = J
  300    CONTINUE
         IF (TE.EQ.0.) GO TO  335
C                                  COLUMN INTERCHANGES
         IF (IPIV.EQ.I0) GO TO 315
         ITE = IC(N,I0)
         IC(N,I0) = IC(N,IPIV)
         IC(N,IPIV) = ITE
         DO 305 IS=I1,I8
            TE = A(IS,I0)
            A(IS,I0) = A(IS,IPIV)
            A(IS,IPIV) = TE
  305    CONTINUE
         DO 310 IS=I4,I2
            TE = C(IS,I0)
            C(IS,I0) = C(IS,IPIV)
            C(IS,IPIV) = TE
  310    CONTINUE
C                                  FACTORIZATION
  315    AK = 1./A(I,I0)
         I01 = I0+1
         I7 = I+1
         DO 325 IS=I7,I8
            XMU = A(IS,I0)*AK
            A(IS,I0) = XMU
            DO 320 JS=I01,M
  320       A(IS,JS) = A(IS,JS)-XMU*A(I,JS)
  325    CONTINUE
  330 CONTINUE
C
      RETURN
  335 SING = .TRUE.
      RETURN
      END
