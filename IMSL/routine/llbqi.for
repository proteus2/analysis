C   IMSL ROUTINE NAME   - LLBQI
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1980
C
C   PURPOSE             - NUCLEUS CALLED ONLY BY IMSL ROUTINE LLBQF
C
C   REQD. IMSL ROUTINES - SINGLE/VBLA=SDSDOT
C                       - DOUBLE/VBLA=DDOT
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE LLBQI (QR,MPN,M,N,M1,N1,D,F,Y)
C
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            MPN,M,N,M1,N1
      REAL               QR(MPN,N),D(N),F(MPN),Y(N)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,J,MH,MV,IS
      REAL               C
C                                  FIRST EXECUTABLE STATEMENT
      MV = 1
      MH = M1
      DO 15 IS=1,N1
         J = M+IS
         IF (IS.NE.M1+1) GO TO 5
         MV = M1+1
         MH = M
    5    Y(IS) = -SDSDOT(IS-1,-F(J),QR(M+1,IS),1,Y(1),1)
         C = 0.0
         IF (IS.GT.M1) C = Y(IS)
         C = SDSDOT(MH-MV+1,-C,QR(MV,IS),1,F(MV),1)/D(IS)
         F(J) = C
         DO 10 I=MV,M
   10    F(I) = F(I)-(C*QR(I,IS))
   15 CONTINUE
      IF (M1.EQ.0) GO TO 35
      DO 20 I=1,M1
   20 F(I) = 0.0
      DO 30 IS=1,M1
         C = SDSDOT(M,-Y(IS),QR(1,IS),1,F(1),1)/D(IS)
         DO 25 I=1,M1
   25    F(I) = F(I)-(C*QR(I,IS))
   30 CONTINUE
   35 CONTINUE
      J = M+N1
   40 F(J) = -SDSDOT(M+N1-J+1,0.0,QR(J,J-M),MPN,F(J),1)
      J = J-1
      IF (J.GE.M+1) GO TO 40
      RETURN
      END
