C   IMSL ROUTINE NAME   - ICSFKV
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1982
C
C   PURPOSE             - NUCLEUS CALLED ONLY BY IMSL SUBROUTINE
C                           ICSFKU
C
C   REQD. IMSL ROUTINES - NONE REQUIRED
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE ICSFKV (MODE,XK,NXK,X,NX,UERROR,TRPZWT,FCTL,FCT,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            MODE,NXK,NX,IER
      REAL               XK(NXK),X(NX),UERROR(NX),TRPZWT(NX),FCTL(NX),
     1                   FCT(NX,1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,IBOUND,ICLAST,ICUBE,ID,IFCT,IFCTL,ILAST,
     1                   ILM1,INSERT,INSIRT,INTERV,IORDER,ITRPZ,IUERR,
     2                   J,K,KLOC,KNOT,L,M,MKN,NJ
      REAL               BC,COEF,COEFL,ERBUT1,FXDKNT,VORD,VORDL,XI,XIL,
     1                   D(28),DIAG(28),DS,DUM1,DUM2,DX,DX2,ONE,Q,
     2                   THREE,TWO,XKNOT,ZERO,BEFORE,AFTER,HALF,HUNDRD
      LOGICAL            MODE3,REPEAT
      COMMON /ICSFK1/    COEFL(4,27),XIL(28),VORDL(28,2),BC(30),
     1                   COEF(4,381),XI(381),VORD(30,28,2),FXDKNT,
     2                   ERBUT1,KNOT,INTERV,ILAST,ICLAST,
     3                   IFCTL,IUERR,ITRPZ,ICUBE,IFCT,MKN,
     4                   INSERT,INSIRT(30),IORDER(28),MODE3
      DATA               ZERO/0.0/,ONE/1.0/,TWO/2.0/,THREE/3.0/,
     1                   HALF/.5/,HUNDRD/100./
C                                  FIRST EXECUTABLE STATEMENT
      XKNOT = XK(1)
      REPEAT = .FALSE.
      BEFORE = ZERO
      IF (MODE.EQ.3) GO TO 35
      IF (MODE.EQ.1) GO TO 15
      IF (MODE.EQ.2) GO TO 15
      ILAST = ILAST+1
      IF (ILAST.GT.4) GO TO 15
C                                  MODE = 0 ENTRY
C                                  CONSTRUCT BASIS FCT 1 TO 4
      ICLAST = ILAST
      IF (ILAST.GT.2) GO TO 10
      IF (ILAST.EQ.2) GO TO 5
C                                  FIRST BASIS FCT IS A CONSTANT
      XIL(1) = XK(1)
      XIL(2) = XK(NXK)
      INSERT = 2
      XI(1) = XIL(1)
      INSIRT(1) = 1
      KNOT = 2
      INTERV = 1
      IORDER(1) = 1
      IORDER(2) = 2
      VORDL(1,1) = ONE
      VORDL(2,1) = ONE
      VORDL(1,2) = ZERO
      VORDL(2,2) = ZERO
      GO TO 175
C                                  SECOND BASIS FCT IS A STRAIGHT LINE
    5 VORDL(2,2) = (VORDL(1,1)+VORDL(1,1))/(XIL(2)-XIL(1))
      VORDL(1,2) = -VORDL(2,2)
C                                  CONSTRUCT THIRD AND FOURTH BASIS FCT
C
   10 VORDL(2,1) = -VORDL(2,1)
      VORDL(2,2) = -VORDL(2,2)
      GO TO 110
   15 CONTINUE
C                                  MODE = 1 OR 2 ENTRY
C                                    INSERT NEW KNOT INTO XIL
      KNOT = KNOT+1
      ILAST = KNOT+2
      ICLAST = ICLAST+INTERV
      INTERV = KNOT-1
      DO 20 I=2,INTERV
         INSERT = I
         IF (XKNOT.LT.XIL(INSERT)) GO TO 25
   20 CONTINUE
   25 CONTINUE
      I = KNOT
      DO 30 L=INSERT,INTERV
         I = I-1
         XIL(I+1) = XIL(I)
         IORDER(I+1) = IORDER(I)
   30 CONTINUE
      IORDER(INSERT) = KNOT
C                                  MODE = 3 ENTRY
   35 CONTINUE
      IF (XKNOT.LE.XIL(INSERT-1)) GO TO 9000
      IF (XKNOT.GE.XIL(INSERT+1)) GO TO 9000
      XIL(INSERT) = XKNOT
C                                  UPDATE VORD FOR FCT M,
C                                    M=1,...,ILAST-1
      ILM1 = ILAST-1
      IBOUND = 0
      DO 50 I=1,ILM1
         ID = IBOUND+1
         IBOUND = ID+IDIM(I,4)
   40    IF (ID.EQ.IBOUND) GO TO 45
         IF (XKNOT.LT.XI(ID+1)) GO TO 45
         ID = ID+1
         GO TO 40
   45    DX = XKNOT-XI(ID)
         VORD(I,KNOT,1) = COEF(1,ID)+DX*(COEF(2,ID)+DX*(COEF(3,ID)+DX*
     1   COEF(4,ID)))
         VORD(I,KNOT,2) = COEF(2,ID)+DX*((COEF(3,ID)+COEF(3,ID))+DX*
     1   THREE*COEF(4,ID))
   50 CONTINUE
      IF (MODE.EQ.2) GO TO 70
   55 CONTINUE
C                                  MODE = 1
C                                  ADD BASIS FCT ILAST
C                                    = DIM(XIL(K),XKNOT)**3
      DO 65 K=1,KNOT
         DX = XIL(K)-XKNOT
         IF (DX.LE.ZERO) GO TO 60
         DX2 = DX*DX
         VORDL(K,1) = DX2*DX
         VORDL(K,2) = THREE*DX2
         GO TO 65
   60    VORDL(K,1) = ZERO
         VORDL(K,2) = ZERO
   65 CONTINUE
      GO TO 110
   70 CONTINUE
C                                  MODE = 2
C                                  REPLACE FCT ILAST BY INTERPOLATING
C                                    IT AT THE CURRENT SET OF KNOTS
C                                  IF FCT ILAST HAS NOT BEEN PREVIOUSLY
C                                    DEFINED PROCEED AS MODE=1
      IF (INSIRT(ILAST).EQ.0) GO TO 55
      VORDL(1,1) = VORD(ILAST,1,1)
      VORDL(1,2) = VORD(ILAST,1,2)
      VORDL(KNOT,1) = VORD(ILAST,2,1)
      VORDL(KNOT,2) = VORD(ILAST,2,2)
      IF (INTERV.LT.2) GO TO 110
      ID = ICLAST
      IBOUND = ICLAST+ILAST-4
      DO 85 K=2,INTERV
   75    IF (ID.EQ.IBOUND) GO TO 80
         IF (XIL(K).LT.XI(ID+1)) GO TO 80
         ID = ID+1
         GO TO 75
   80    DX = XIL(K)-XI(ID)
         VORDL(K,1) = COEF(1,ID)+DX*(COEF(2,ID)+DX*(COEF(3,ID)+DX*
     1   COEF(4,ID)))
   85 CONTINUE
C                                  INTERPOLATE
C                                    COMPUTE THE SLOPES VORDL(I,2),
C                                    I=2,...,KNOT-1 AT INTERIOR KNOTS OF
C                                    CUBIC SPLINE FOR GIVEN VALUES
C                                    VORDL(I,1),I=1,KNOT, AT ALL THE
C                                    KNOTS AND GIVEN BOUNDARY
C                                    DERIVATIVES
      D(1) = ZERO
      DIAG(1) = ONE
      DO 90 M=2,KNOT
         D(M) = XIL(M)-XIL(M-1)
         DIAG(M) = (VORDL(M,1)-VORDL(M-1,1))/D(M)
   90 CONTINUE
      DO 95 M=2,INTERV
         VORDL(M,2) = THREE*(D(M)*DIAG(M+1)+D(M+1)*DIAG(M))
         DIAG(M) = TWO*(D(M)+D(M+1))
   95 CONTINUE
      DO 100 M=2,INTERV
         Q = -D(M+1)/DIAG(M-1)
         DIAG(M) = DIAG(M)+Q*D(M-1)
         VORDL(M,2) = VORDL(M,2)+Q*VORDL(M-1,2)
  100 CONTINUE
      NJ = KNOT
      DO 105 M=2,INTERV
         NJ = NJ-1
         VORDL(NJ,2) = (VORDL(NJ,2)-D(NJ)*VORDL(NJ+1,2))/DIAG(NJ)
  105 CONTINUE
  110 CONTINUE
      IF (REPEAT) GO TO 120
C                                  UPDATE XI AND INSIRT
      ID = ICLAST-1
      DO 115 I=1,INTERV
         ID = ID+1
         XI(ID) = XIL(I)
  115 CONTINUE
      INSIRT(ILAST) = INSERT
C                                  EVALUATE COEFL AND FCTL
  120 CONTINUE
      DO 125 I=1,INTERV
         COEFL(1,I) = VORDL(I,1)
         COEFL(2,I) = VORDL(I,2)
         DX = XIL(I+1)-XIL(I)
         DUM1 = (VORDL(I+1,1)-VORDL(I,1))/DX
         DUM2 = VORDL(I,2)+VORDL(I+1,2)-DUM1-DUM1
         COEFL(3,I) = (DUM1-DUM2-VORDL(I,2))/DX
         COEFL(4,I) = DUM2/(DX*DX)
  125 CONTINUE
      J = 1
      DO 140 L=1,NX
  130    IF (J.EQ.INTERV) GO TO 135
         IF (X(L).LT.XIL(J+1)) GO TO 135
         J = J+1
         GO TO 130
  135    DX = X(L)-XIL(J)
         FCTL(L) = COEFL(1,J)+DX*(COEFL(2,J)+DX*(COEFL(3,J)+DX*
     1   COEFL(4,J)))
  140 CONTINUE
C                                  APPLY MODIFIED GRAM-SCHMIDT PROCESS
C                                    TO FCTL
      DS = ZERO
      DO 145 L=1,NX
  145 DS = DS+FCTL(L)*FCTL(L)*TRPZWT(L)
      IF (DS.LT.ZERO) DS=ZERO
      BEFORE = SQRT(DS)
      Q = ZERO
      IF (REPEAT.AND.ABS(BEFORE-ONE).GT.HALF) GO TO 205
      ILM1 = ILAST-1
      DO 170 J=1,ILM1
         DS = ZERO
         DO 150 L=1,NX
  150    DS = DS+FCT(L,J)*FCTL(L)*TRPZWT(L)
         DO 155 L=1,NX
  155    FCTL(L) = FCTL(L)-DS*FCT(L,J)
         DO 165 K=1,KNOT
            KLOC = IORDER(K)
            DO 160 I=1,2
  160       VORDL(K,I) = VORDL(K,I)-DS*VORD(J,KLOC,I)
  165    CONTINUE
  170 CONTINUE
  175 CONTINUE
C                                  EVALUATE NEW COEFL AND FCTL
      DO 180 I=1,INTERV
         COEFL(1,I) = VORDL(I,1)
         COEFL(2,I) = VORDL(I,2)
         DX = XIL(I+1)-XIL(I)
         DUM1 = (VORDL(I+1,1)-VORDL(I,1))/DX
         DUM2 = VORDL(I,2)+VORDL(I+1,2)-DUM1-DUM1
         COEFL(3,I) = (DUM1-DUM2-VORDL(I,2))/DX
         COEFL(4,I) = DUM2/(DX*DX)
  180 CONTINUE
      J = 1
      DO 195 L=1,NX
  185    IF (J.EQ.INTERV) GO TO 190
         IF (X(L).LT.XIL(J+1)) GO TO 190
         J = J+1
         GO TO 185
  190    DX = X(L)-XIL(J)
         FCTL(L) = COEFL(1,J)+DX*(COEFL(2,J)+DX*(COEFL(3,J)+DX*
     1   COEFL(4,J)))
  195 CONTINUE
C                                  COMPUTE THE COMPONENT BC(ILAST) OF
C                                    UERROR WRTO FCTL AND UPDATE THE
C                                    VARIOUS REPRESENTATIONS OF FCT
      DS = ZERO
      DO 200 L=1,NX
  200 DS = DS+FCTL(L)*FCTL(L)*TRPZWT(L)
      IF (DS.LT.ZERO) DS=ZERO
      AFTER = SQRT(DS)
      Q = ZERO
      IF (REPEAT.AND.ABS(AFTER-ONE).GT.HALF) GO TO 205
      DUM1 = AFTER + BEFORE
      IF (DUM1.NE.BEFORE) Q = ONE/AFTER
  205 CONTINUE
      DS = ZERO
      DO 210 L=1,NX
  210 DS = DS+FCTL(L)*UERROR(L)*TRPZWT(L)
      BC(ILAST) = DS*Q
      DO 220 K=1,KNOT
         KLOC = IORDER(K)
         DO 215 L=1,2
  215    VORD(ILAST,KLOC,L) = VORDL(K,L)*Q
  220 CONTINUE
      ID = ICLAST-1
      DO 230 I=1,INTERV
         ID = ID+1
         DO 225 L=1,4
  225    COEF(L,ID) = COEFL(L,I)*Q
  230 CONTINUE
      DO 235 L=1,NX
  235 FCT(L,ILAST) = FCTL(L)*Q
      IF (Q.EQ.ZERO) GO TO 9005
      IF (ILAST.EQ.1) GO TO 9005
      IF (HUNDRD*AFTER.GT.BEFORE) GO TO 9005
C                                  RE-ORTHOGONALIZE FCTL
      REPEAT = .TRUE.
      GO TO 70
 9000 CONTINUE
      IER = 129
 9005 RETURN
      END
