C   IMSL ROUTINE NAME   - ICSSE
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1980
C
C   PURPOSE             - NUCLEUS CALLED ONLY BY IMSL SUBROUTINE ICSSCV
C
C   REQD. IMSL ROUTINES - NONE REQUIRED
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      REAL FUNCTION ICSSE (PLOG,A,X,IDUM,N,NMAX)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N,IDUM(1),NMAX
      REAL               PLOG,A(1),X(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,NL2
      REAL               AIP,F1,F2,P,T1,T2SQ,T2,T3,T4,VDERIV
      COMMON /VCOM/      VDERIV
C                                  CALCULATE THE CROSS VALIDATION
C                                    FUNCTION AND ITS DERIVATIVE
C                                  FIRST EXECUTABLE STATEMENT
      P = EXP(PLOG)
      T1 = 0.0
      T2 = 0.0
      T3 = 0.0
      T4 = 0.0
      NL2 = N-2
      DO 5 I=1,NL2
         AIP = 1.0/(A(I)+P)
         F1 = A(I)*AIP
         F2 = F1*X(I)
         T1 = T1+F2*F2
         T2 = T2+F1
         T3 = T3+F1*AIP
         T4 = T4+F2*F2*AIP
    5 CONTINUE
      T2SQ = T2*T2
      AIP = N
      ICSSE = AIP*T1/T2SQ
      VDERIV = 2.*AIP*P*(T1*T3-T2*T4)/T2SQ/T2
      RETURN
      END
