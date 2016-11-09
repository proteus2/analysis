C   IMSL ROUTINE NAME   - RLDCQM
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - DECODING OF A QUADRATIC REGRESSION MODEL
C
C   USAGE               - CALL RLDCQM (XBAR,M,BN)
C
C   ARGUMENTS    XBAR   - INPUT VECTOR OF LENGTH M*(M+3)/2 CONTAINING
C                           THE MEANS OF THE ORIGINAL INDEPENDENT
C                           VARIABLES AND THE MEANS OF THE GENERATED
C                           SQUARE AND CROSS PRODUCT VARIABLES.
C                           THE COMPONENTS OF XBAR ARE IN STANDARD
C                           ORDER. SEE PROGRAMMING NOTES IN THE
C                           MANUAL DOCUMENT FOR FURTHER DETAILS.
C                M      - INPUT NUMBER OF ORIGINAL INDEPENDENT
C                           VARIABLES.
C                BN     - INPUT AND OUTPUT VECTOR OF LENGTH
C                           M*(M+3)/2+1.
C                         ON INPUT, BN CONTAINS THE INTERCEPT AND
C                           REGRESSION COEFFICIENTS OF THE CENTERED
C                           MODEL IN STANDARD ORDER.
C                         ON OUTPUT, BN CONTAINS THE INTERCEPT AND
C                           REGRESSION COEFFICIENTS OF THE UNCENTERED
C                           MODEL IN STANDARD ORDER.
C                         SEE PROGRAMMING NOTES IN THE MANUAL
C                           DOCUMENT FOR FURTHER DETAILS.
C
C   REQD. IMSL ROUTINES - NONE REQUIRED
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE RLDCQM (XBAR,M,BN)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            M
      REAL               XBAR(1),BN(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            K,KK,J,I,L,M2,II,I1,IM1
      DOUBLE PRECISION   TEMP,T,TT
C                                  FIRST EXECUTABLE STATEMENT
      K=M+2
      KK=M*(M+3)/2+1
      TEMP=0.D0
      J=1
      DO 5 I=2,KK
         TEMP=TEMP+DBLE(BN(I))*DBLE(XBAR(J))
         J=I
    5 CONTINUE
      TEMP=BN(1)-TEMP
      T=0.D0
      DO 10 L=1,M
         TT=DBLE(XBAR(L))
         T=T+DBLE(BN(K))*TT*TT
         K=K+1
   10 CONTINUE
      TEMP=TEMP+T
      T=0.D0
      M2=M-1
      DO 15 L=1,M2
         II=L+1
      DO 15 J=II,M
         T=T+DBLE(BN(K))*DBLE(XBAR(L))*DBLE(XBAR(J))
   15    K=K+1
      BN(1)=TEMP+T
      K=M+M
      M2=K
      DO 35 I=1,M
         II=M2+I
         I1=I+1
         TEMP=XBAR(I)
         TT=BN(I1)-BN(M+I1)*(TEMP+TEMP)
         IF (I .EQ. 1) GO TO 25
         T=0.D0
         IM1=I-1
         DO 20 J=1,IM1
            T=T+DBLE(BN(II))*DBLE(XBAR(J))
            II=II+M-J-1
   20    CONTINUE
         IF (I .EQ. M) GO TO 35
         TT=TT-T
   25    T=0.D0
         DO 30 J=I1,M
   30       T=T+DBLE(XBAR(J))*DBLE(BN(K+J))
         K=K+M-I-1
   35    BN(I1)=TT-T
      RETURN
      END
