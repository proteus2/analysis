C   IMSL ROUTINE NAME   - LLBQG
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1982
C
C   PURPOSE             - NUCLEUS CALLED ONLY BY IMSL ROUTINE LLBQF
C
C   REQD. IMSL ROUTINES - SINGLE/VBLA=SDOT,VBLA=SDSDOT
C                       - DOUBLE/VBLA=DDOT
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE LLBQG  (QR,MPN,NP1,M,N,M1,N1,BASIC,TOL,IPIVOT,D)
C
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            MPN,NP1,M,N,M1,N1,IPIVOT(N)
      REAL               QR(MPN,NP1),TOL,D(N)
      LOGICAL            BASIC
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,J,KM1,K,MF,MH,MS,MV,IP,ISP1,IS
      REAL               C,DJ,DM,DSV,DS,PJ,PS,RSJ
      LOGICAL            FINIS,RANK
C                                  THIS SUBROUTINE USES THE MODIFIED
C                                    GRAM-SCHMIDT ALGORITHM WITH COLUMN
C                                    PIVOTING TO OBTAIN THE
C                                    DECOMPOSITION OF THE MATRIX STORED
C                                    IN QR NEEDED FOR ITERATIVE
C                                    REFINEMENT. N1 IS SET TO THE
C                                    NUMERICAL RANK OF QR. IF THE FIRST
C                                    M1 ROWS OF QR MODIFIED BY ROUNDING
C                                    ERRORS ARE LINEARLY DEPENDENT,
C                                    LLBQG RETURNS WITH N1 SET LESS
C                                    THAN M1. IF QR IS THE ZERO MATRIX,
C                                    LLBQG RETURNS WITH N1 SET TO 0.
C                                  FIRST EXECUTABLE STATEMENT
      MV = 1
      MH = M1
      N1 = N
      MS = M
      MF = 1
      FINIS=.FALSE.
      DO 5 J=1,N
         D(J) = SDOT(M,QR(1,J),1,QR(1,J),1)
    5 IPIVOT(J) = J
      DO 85 IS=1,N
C                                  BEGIN STEP NUMBER S OF THE
C                                    DECOMPOSITION
         K = M+IS
         IF (IS.NE.M1+1) GO TO 10
         MV = M1+1
         MH = M
   10    IF (FINIS) GO TO 25
         PS = 0.0
         IP = IS
C                                  BEGIN PIVOT SEARCH
         DO 15 J=IS,N
            DJ = SDOT(MH-MV+1,QR(MV,J),1,QR(MV,J),1)
            IF (DJ.EQ.0.0) GO TO 15
            PJ = -SDSDOT(MH-MV+1,0.0,QR(MV,J),1,QR(MV,NP1),1)**2/DJ
            IF (PS.LE.PJ) GO TO 15
            PS = PJ
            IP = J
   15    CONTINUE
         IF (IP.EQ.IS) GO TO 30
C                                  BEGIN COLUMN EXCHANGE
         DS = D(IP)
         D(IP) = D(IS)
         D(IS) = DS
         I = IPIVOT(IP)
         IPIVOT(IP) = IPIVOT(IS)
         IPIVOT(IS) = I
         KM1 = K-1
         DO 20 I=1,KM1
            C = QR(I,IP)
            QR(I,IP) = QR(I,IS)
            QR(I,IS) = C
   20    CONTINUE
C                                  END COLUMN EXCHANGE END PIVOT SEARCH
         GO TO 30
   25    MH = K-1
         MS = MH
   30    CONTINUE
         C = 0.0
         IF (FINIS) C = 1.0
         DS = SDSDOT(MH-MV+1,C,QR(MV,IS),1,QR(MV,IS),1)
         DSV = D(IS)
         D(IS) = DS
         IF (FINIS) GO TO 70
         IF (IS.GT.M) GO TO 35
         RANK = SQRT(DS).LE.TOL*SQRT(DSV)
         IF (.NOT.RANK) GO TO 70
   35    FINIS=.TRUE.
         N1 = IS-1
         IF ((N1.LT.M1).OR.(N1.EQ.0).OR.BASIC) GO TO 90
         MV = M+1
         MF = MV
         DO 65 IP=IS,N
            IF (M1.EQ.0) GO TO 55
            DO 40 I=1,M1
   40       QR(I,IP) = 0.0
            DO 50 J=1,M1
               C = SDSDOT(M,0.0,QR(1,J),1,QR(1,IP),1)/D(J)
               DO 45 I=1,M1
   45          QR(I,IP) = QR(I,IP)-C*QR(I,J)
   50       CONTINUE
   55       J = M+N1
   60       QR(J,IP) = -SDSDOT(M+N1-J+1,0.0,QR(J,J-M),MPN,QR(J,IP),1)
            J = J-1
            IF (J.GE.M+1) GO TO 60
   65    CONTINUE
         GO TO 25
   70    CONTINUE
         QR(K,IS) = -1.0
         ISP1 = IS+1
         IF (ISP1.GT.N) GO TO 85
C                                  BEGIN ORTHOGONALIZATION
         DO 80 J=ISP1,NP1
            RSJ = SDSDOT(MH-MV+1,0.0,QR(MV,J),1,QR(MV,IS),1)/DS
            QR(K,J) = RSJ
            DO 75 I=MF,MS
   75       QR(I,J) = QR(I,J)-RSJ*QR(I,IS)
   80    CONTINUE
C                                  END ORTHOGONALIZATION
   85 CONTINUE
C                                  END STEP S
   90 RETURN
      END
