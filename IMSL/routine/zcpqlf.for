C   IMSL ROUTINE NAME   - ZCPQLF
C
C-----------------------------------------------------------------------
C
C   COMPUTER            - IBM/SINGLE
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - NUCLEUS CALLED ONLY BY IMSL SUBROUTINE
C                           ZCPOLY
C
C   REQD. IMSL ROUTINES - NONE REQUIRED
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE ZCPQLF (BOWL)
C                                  SPECIFICATIONS FOR ARGUMENTS
      LOGICAL            BOWL
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            NN,N,J
      REAL               PR(50),PI(50),HR(50),HI(50),QPR(50),QPI(50),
     1                   QHR(50),QHI(50),SHR(50),SHI(50),
     2                   SR,SI,TR,TI,PVR,PVI,ARE,RMRE,REPSR1,RINFP,
     3                   T1,T2,ZERO
      COMMON /ZCPQLN/    PR,PI,HR,HI,QPR,QPI,QHR,QHI,SHR,SHI,SR,SI,
     1                   TR,TI,PVR,PVI,ARE,RMRE,REPSR1,RINFP,NN
      DATA               ZERO/0.0/
C                                  FIRST EXECUTABLE STATEMENT
      N = NN-1
C                                  CALCULATES THE NEXT SHIFTED H
C                                    POLYNOMIAL
C                                  BOWL - LOGICAL, IF .TRUE. H(S) IS
C                                    ESSENTIALLY ZERO
      IF (BOWL) GO TO 10
      DO 5 J=2,N
         T1 = QHR(J-1)
         T2 = QHI(J-1)
         HR(J) = TR*T1-TI*T2+QPR(J)
         HI(J) = TR*T2+TI*T1+QPI(J)
    5 CONTINUE
      HR(1) = QPR(1)
      HI(1) = QPI(1)
      RETURN
C                                  IF H(S) IS ZERO REPLACE H WITH QH
   10 DO 15 J=2,N
         HR(J) = QHR(J-1)
         HI(J) = QHI(J-1)
   15 CONTINUE
      HR(1) = ZERO
      HI(1) = ZERO
      RETURN
      END
