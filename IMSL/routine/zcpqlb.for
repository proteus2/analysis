C   IMSL ROUTINE NAME   - ZCPQLB
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - NUCLEUS CALLED ONLY BY IMSL SUBROUTINE
C                           ZCPOLY
C
C   REQD. IMSL ROUTINES - ZCPQLK,ZCPQLL,ZCPQLM
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE ZCPQLB (L1)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            L1
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            N,NN,NM1,I,JJ,J
      REAL               PR(50),PI(50),HR(50),HI(50),QPR(50),QPI(50),
     1                   QHR(50),QHI(50),SHR(50),SHI(50),
     2                   SR,SI,TR,TI,PVR,PVI,ARE,RMRE,REPSR1,RINFP,
     3                   XNI,T1,T2,ZCPQLL,ZERO,TEN,ONEDN,ONE
      COMMON /ZCPQLN/    PR,PI,HR,HI,QPR,QPI,QHR,QHI,SHR,SHI,SR,SI,
     1                   TR,TI,PVR,PVI,ARE,RMRE,REPSR1,RINFP,NN
      DATA               ZERO,TEN/0.0,10.0/,ONE/1.0/
C                                  FIRST EXECUTABLE STATEMENT
      N = NN-1
      NM1 = N-1
      ONEDN = ONE/N
C                                  COMPUTES THE DERIVATIVE POLYNOMIAL
C                                    AS THE INITIAL H POLYNOMIAL AND
C                                    COMPUTES L1 NO-SHIFT H
C                                    POLYNOMIALS.
      DO 5 I=1,N
         XNI = NN-I
         HR(I) = XNI*PR(I)*ONEDN
         HI(I) = XNI*PI(I)*ONEDN
    5 CONTINUE
      DO 25 JJ=1,L1
         IF (ZCPQLL(HR(N),HI(N)).LE.REPSR1*TEN*ZCPQLL(PR(N),PI(N)))
     1   GO TO 15
         CALL ZCPQLK (-PR(NN),-PI(NN),HR(N),HI(N),TR,TI)
         DO 10 I=1,NM1
            J = NN-I
            T1 = HR(J-1)
            T2 = HI(J-1)
            HR(J) = TR*T1-TI*T2+PR(J)
            HI(J) = TR*T2+TI*T1+PI(J)
   10    CONTINUE
         HR(1) = PR(1)
         HI(1) = PI(1)
         GO TO 25
C                                  IF THE CONSTANT TERM IS ESSENTIALLY
C                                    ZERO, SHIFT H COEFFICIENTS
   15    DO 20 I=1,NM1
            J = NN-I
            HR(J) = HR(J-1)
            HI(J) = HI(J-1)
   20    CONTINUE
         HR(1) = ZERO
         HI(1) = ZERO
   25 CONTINUE
      RETURN
      END
