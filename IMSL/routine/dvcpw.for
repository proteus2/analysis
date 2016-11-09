C   IMSL ROUTINE NAME   - DVCPW
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1982
C
C   PURPOSE             - NUCLEUS CALLED ONLY BY IMSL SUBROUTINE DVCPR
C
C   REQD. IMSL ROUTINES - DVCPX,DVCPY
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE DVCPW (M,N,P,R,X,Y,A1,B1,A,C,DEL,CASI,SING,IR,IC,UU,
     *RES,LIN,MMAX,MTNMAX,NMAX,MMAX2,HX,GRADF,AUX,ICA,XAU,FCNJ,FCNB)
C
C                                  SPECIFICATIONS FOR ARGUMENTS
C
      INTEGER            M,N,P,R,MMAX,MTNMAX,NMAX,MMAX2,IR(NMAX,MMAX),
     *                   IC(NMAX,MMAX),ICA(1)
      REAL               X(NMAX),Y(MTNMAX),A1(MMAX,MMAX),B1(MMAX,MMAX),
     *                   A(MTNMAX,MMAX),C(MTNMAX,MMAX),DEL(MMAX,MTNMAX),
     *                   UU(MTNMAX),RES(MTNMAX),HX(NMAX),GRADF(MTNMAX),
     *                   AUX(MMAX2,MMAX),XAU(MTNMAX)
      LOGICAL            CASI,SING,LIN
C                                  SPECIFICATIONS FOR COMMON /NEWT /
      COMMON             /NEWT/ INWT,NU,CAS1
      INTEGER            INWT,NU
      LOGICAL            CAS1
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
C
      INTEGER            I1J,I1,I2J,I2,I3P,I3,I4,I7,I8,I9,II,IM,IP,IT,I,
     *                   J7,JJ,J,K,M1,MMP,N1,P1,PM,R1
      REAL               ECOND,H1,HAK,HBK,REPS,TE
      DATA               REPS/Z3C100000/
C                                  CASI = .TRUE. NO DECOMPOSITION
C                                  FIRST EXECUTABLE STATEMENT
      IF (CASI) GO TO 165
      P1 = P-1
      N1 = N-1
      M1 = M-P
      MMP = M1-1
      PM = P+1
      R1 = R-1
C                                  CONSTRUCTION OF A,C,DEL
      DO 15 I=1,N
         JJ = (I-1)*M
         I8 = (I-1)*M+1
         CALL FCNJ(M,X(I),Y(I8),DEL)
         DO 10 K=1,M
            DO 5 J=1,M
               C(JJ+J,K) = DEL(J,K)
    5       CONTINUE
   10    CONTINUE
   15 CONTINUE
      ECOND = SQRT(REPS)
      DO 30 I7=1,M
         HAK = AMAX1(Y(I7),0.1)*ECOND
         Y(I7) = Y(I7)+HAK
         CALL FCNB(M,Y(1),Y(I8),DEL(1,2))
         Y(I7) = Y(I7)-2.0*HAK
         CALL FCNB(M,Y(1),Y(I8),DEL(1,1))
         Y(I7) = Y(I7)+HAK
         DO 20 J7=1,M
            A1(J7,I7) = (DEL(J7,2)-DEL(J7,1))/(2.0*HAK)
   20    CONTINUE
         I9 = I8+I7-1
         HBK = AMAX1(Y(I9),0.1)*ECOND
         Y(I9) = Y(I9)+HBK
         CALL FCNB(M,Y(1),Y(I8),DEL(1,2))
         Y(I9) = Y(I9)-2.0*HBK
         CALL FCNB(M,Y(1),Y(I8),DEL(1,1))
         Y(I9) = Y(I9)+HBK
         DO 25 J7=1,M
            B1(J7,I7) = (DEL(J7,2)-DEL(J7,1))/(2.0*HBK)
   25    CONTINUE
   30 CONTINUE
      CALL FCNJ(M,X(N),Y(I8),DEL)
      IF (P.EQ.0) GO TO 45
      DO 40 I=1,P
         DO 35 J=1,M
   35    A(I,J) = A1(I,J)
   40 CONTINUE
   45 DO 85 II=1,N1
         H1 = .5*HX(II)
         I1 = (II-1)*M+1
         I2 = I1+MMP
         IT = 0
         DO 55 I=I1,I2
            IP = I+P
            IT = IT+1
            DO 50 J=1,M
               A(IP,J) = -H1*C(I,J)
               IF (IT.EQ.J) A(IP,J) = A(IP,J)-1.
   50       CONTINUE
   55    CONTINUE
         I3 = I2+M
         I4 = I1+P1
         IT = 0
         IF (P.EQ.0) GO TO 70
         DO 65 I=I1,I4
            IT = IT+1
            DO 60 J=1,M
               IM = I+M
               A(IM,J) = -H1*C(I3+IT,J)
               C(I,J) = -H1*C(I2+IT,J)
               IF (IT+M1.NE.J) GO TO 60
               A(IM,J) = A(IM,J)+1.
               C(I,J) = C(I,J)-1.
   60       CONTINUE
   65    CONTINUE
   70    IT = 0
         DO 80 I=I1,I2
            IP = I+P
            IT = IT+1
            IM = I+M
            DO 75 J=1,M
               C(IP,J) = -H1*C(IM,J)
               IF (IT.EQ.J) C(IP,J) = C(IP,J)+1.
   75       CONTINUE
   80    CONTINUE
   85 CONTINUE
C
      I1 = N1*M+PM
      I2 = N*M
      IT = P
      DO 95 I=I1,I2
         IT = IT+1
         DO 90 J=1,M
   90    A(I,J) = B1(IT,J)
   95 CONTINUE
      IF (R.EQ.0) GO TO 110
      DO 105 I=1,R
         DO 100 J=1,M
  100    DEL(I,J) = A1(I+P,J)
  105 CONTINUE
  110 IF (LIN) GO TO 160
C                                  COMPUTATION OF GRADIENT
      IT = 0
      DO 145 I=1,N
         I1 = (I-1)*M
         I2 = I1-M1
         I3 = I*M
         DO 140 J=1,M
            IT = IT+1
            TE = 0.
            DO 115 JJ=1,M
               I1J = I1+JJ
               TE = TE+A(I1J,J)*RES(I1J)
  115       CONTINUE
            IF (I.EQ.N) GO TO 125
            DO 120 JJ=1,P
  120       TE = TE+C(I1+JJ,J)*RES(I3+JJ)
  125       IF (I.EQ.1) GO TO 135
            DO 130 JJ=1,M1
               I2J = I2+JJ
               TE = TE+C(I2J,J)*RES(I2J)
  130       CONTINUE
  135       GRADF(I1+J) = TE
  140    CONTINUE
  145 CONTINUE
C
      IF (R.EQ.0) GO TO 160
      I3P = N1*M+P
      DO 155 J=1,M
         TE = 0.
         DO 150 I=1,R
  150    TE = TE+DEL(I,J)*RES(I3P+I)
         GRADF(J) = GRADF(J)+TE
  155 CONTINUE
  160 CALL DVCPX(A,C,DEL,M,N,P,R,IR,IC,SING,ICA,AUX,MTNMAX,MMAX,NMAX,
     *MMAX2)
      IF (SING) GO TO 170
  165 CALL DVCPY(A,C,DEL,RES,M,N,P,R,IR,IC,UU,MTNMAX,MMAX,NMAX,XAU)
  170 RETURN
      END
