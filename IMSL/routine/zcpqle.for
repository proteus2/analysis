C   IMSL ROUTINE NAME   - ZCPQLE
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - NUCLEUS CALLED ONLY BY IMSL SUBROUTINE
C                           ZCPOLY
C
C   REQD. IMSL ROUTINES - ZCPQLG,ZCPQLK,ZCPQLL,ZCPQLM
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE ZCPQLE (BOWL)
C                                  SPECIFICATIONS FOR ARGUMENTS
      LOGICAL            BOWL
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            NN,N
      REAL               PR(50),PI(50),HR(50),HI(50),QPR(50),QPI(50),
     1                   QHR(50),QHI(50),SHR(50),SHI(50),
     2                   SR,SI,TR,TI,PVR,PVI,ARE,RMRE,REPSR1,RINFP,
     3                   HVR,HVI,ZCPQLL,ZERO,TEN
      COMMON /ZCPQLN/    PR,PI,HR,HI,QPR,QPI,QHR,QHI,SHR,SHI,SR,SI,
     1                   TR,TI,PVR,PVI,ARE,RMRE,REPSR1,RINFP,NN
      DATA               ZERO,TEN/0.0,10.0/
C                                  FIRST EXECUTABLE STATEMENT
      N = NN-1
C                                  COMPUTES T = -P(S)/H(S).
C                                  BOWL - LOGICAL, SET TRUE IF H(S) IS
C                                    ESSENTIALLY ZERO.
C                                  EVALUATE H(S)
      CALL ZCPQLG (N,SR,SI,HR,HI,QHR,QHI,HVR,HVI)
      BOWL = ZCPQLL(HVR,HVI).LE.ARE*TEN*ZCPQLL(HR(N),HI(N))
      IF (BOWL) GO TO 5
      CALL ZCPQLK (-PVR,-PVI,HVR,HVI,TR,TI)
      RETURN
    5 TR = ZERO
      TI = ZERO
      RETURN
      END
