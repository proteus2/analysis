C   IMSL ROUTINE NAME   - LLBQH
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1981
C
C   PURPOSE             - NUCLEUS CALLED ONLY BY IMSL ROUTINE LLBQF
C
C   PRECISION/HARDWARE  - SINGLE AND DOUBLE/H32
C                       - SINGLE/H36,H48,H60
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE LLBQH  (A,IA,M,N,M1,N1,B,QR,MPN,BASIC,X,
     1                   IPIVOT,RES,D,F,Y,IFAIL)
C
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            IA,M,N,M1,N1,MPN,IFAIL,IPIVOT(N)
      REAL               A(IA,N),B(M),QR(MPN,N),X(N),RES(M),D(N),F(MPN),
     1                   Y(N)
      LOGICAL            BASIC
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,J,K,MAXIT,N1P1,IP,IS
      REAL               C,RNDR1,RNDR2,RNDX1,RNDX2,RNR,RNX,RN
      DOUBLE PRECISION   DC
C                                  THIS SUBROUTINE USES THE
C                                    DECOMPOSITION STORED IN QR FOR THE
C                                    ITERATIVE REFINEMENT OF THE
C                                    SOLUTION CORRESPONDING TO THE
C                                    RIGHT HAND SIDE GIVEN IN B. IF THE
C                                    SOLUTION FAILS TO IMPROVE
C                                    SUFFICIENTLY, LLBQH RETURNS WITH
C                                    IFAIL SET TO IFAIL+1.
C                                    IF IFAIL .LE. 0,
C                                    THEN THE SUBROUTINE SOLVES THE
C                                    SYSTEM WITH THE (-IFAIL)-TH COLUMN
C                                    OF THE IDENTITY IN PLACE OF B.
C                                  FIRST EXECUTABLE STATEMENT
      N1P1 = N1+1
      MAXIT = 100
      DO 5 I=1,M
         F(I) = 0.0
         IF (IFAIL.GE.0) F(I) = B(I)
         IF (I.EQ.-IFAIL) F(I) = 1.0
         RES(I) = 0.0
    5 CONTINUE
      DO 10 IS=1,N
         X(IS) = 0.0
         J = M+IS
         F(J) = 0.0
         IF (M+IPIVOT(IS).EQ.-IFAIL) F(J) = 1.0
   10 CONTINUE
      DO 80 K=1,MAXIT
C                                  BEGIN K-TH ITERATION
         IF (K.EQ.1) GO TO 50
         IF (K.EQ.2) GO TO 15
         IF (((RNDX2.GE.RNDX1).OR.(RNX+RNDX2.EQ.RNX)).AND.
     1   ((RNDR2.GE.RNDR1).OR.(RNR+RNDR2.EQ.RNR))) GO TO 85
   15    RNDX1 = RNDX2
         RNDR1 = RNDR2
C                                  COMPUTE NEW RESIDUALS
         DO 20 I=1,M
   20    RES(I) = RES(I)+F(I)
         DO 25 IS=1,N
            J = M+IS
            IP = IPIVOT(IS)
            X(IP) = X(IP)+F(J)
            F(J) = 0.0
            C = 0.0
            IF (M+IP.EQ.-IFAIL) C = 1.0
            IF (IS.LE.N1) F(J) = -SDSDOT(M,-C,A(1,IP),1,RES(1),1)
   25    CONTINUE
         DO 45 I=1,M
            DC = 0.0D0
            IF (IFAIL.GE.0) DC = DBLE(B(I))
            IF (I.EQ.-IFAIL) DC = 1.0D0
            IF (I.GT.M1) DC = DC-DBLE(RES(I))
            IF ((N1.EQ.N).OR.BASIC) GO TO 35
            DO 30 IS=N1P1,N
               IP = IPIVOT(IS)
               DC = DC+DBLE(X(IP))*DBLE(QR(I,IS))
   30       CONTINUE
   35       DO 40 J=1,N
   40       DC = DC-DBLE(A(I,J))*DBLE(X(J))
            F(I) = DC
   45    CONTINUE
C                                  END NEW RESIDUALS
   50    CALL LLBQI (QR,MPN,M,N,M1,N1,D,F,Y)
         IF ((N1.EQ.N).OR.BASIC) GO TO 65
         DO 60 IS=N1P1,N
            C = -SDSDOT(IS,0.0,QR(M+1,IS),1,F(M+1),1)/D(IS)
            DO 55 J=1,IS
               I = M+J
   55       F(I) = F(I)+C*QR(I,IS)
   60    CONTINUE
   65    CONTINUE
         RNDX2 = 0.0
         DO 70 I=1,N
   70    RNDX2 = AMAX1(ABS(F(I+M)),RNDX2)
         RNDR2 = 0.0
         DO 75 I=1,M
   75    RNDR2 = AMAX1(ABS(F(I)),RNDR2)
         IF (K.NE.1) GO TO 80
         RNX = RNDX2
         RNR = RNDR2
   80 CONTINUE
C                                  END K-TH ITERATION
   85 CONTINUE
C                                  CONVERGENCE TEST
C
      RN=MAX0(10,N)
      RNDX1=RNX+RNDX2/RN
      RNDR1=RNR+RNDR2/RN
      IF ((RNDX1.NE.RNX).AND.(RNDR1.NE.RNR)) IFAIL=IFAIL+1
      RETURN
      END
