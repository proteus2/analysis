C   IMSL ROUTINE NAME   - FTTR
C
C-----------------------------------------------------------------------
C   LATEST REVISION     - JUNE 1, 1982
C
C   PURPOSE             - PARAMETER ESTIMATES FOR A UNIVARIATE
C                           TRANSFER FUNCTION MODEL
C
C   USAGE               - CALL FTTR (Y,X,N,T,IB,IND,A,IER)
C
C   ARGUMENTS    Y      - INPUT TIME SERIES OF LENGTH N REPRESENTING THE
C                           OUTPUT FROM A DYNAMIC SYSTEM. ON OUTPUT THE
C                           MEAN HAS BEEN REMOVED FROM Y AND IS STORED
C                           IN T(8). SEE DESCRIPTION OF T BELOW.
C                X      - INPUT TIME SERIES OF LENGTH N REPRESENTING THE
C                           INPUT TO A DYNAMIC SYSTEM. ON OUTPUT THE
C                           MEAN HAS BEEN REMOVED FROM X AND IS STORED
C                           IN T(9). SEE DESCRIPTION OF T BELOW.
C                N      - INPUT LENGTH OF Y AND X.
C                           N MUST BE GREATER THAN OR EQUAL TO 41.
C                T      - OUTPUT VECTOR OF LENGTH 10.
C                         T(1),T(2),...,T(7) CONTAIN THE ESTIMATES OF
C                           THE MODEL PARAMETERS. SEE THE PROGRAMMING
C                           NOTES IN THE MANUAL DOCUMENT FOR FURTHER
C                           DETAILS.
C                         T(8) CONTAINS THE MEAN OF SERIES Y.
C                         T(9) CONTAINS THE MEAN OF SERIES X.
C                         T(10) CONTAINS THE SUM OF SQUARES OF THE
C                           ESTIMATED WHITE NOISE SERIES. SEE
C                           SEE PROGRAMMING NOTES IN
C                           THE MANUAL DOCUMENT FOR FURTHER DETAILS.
C                IB     - OUTPUT LAG PARAMETER. SEE THE PROGRAMMING
C                           NOTES IN THE MANUAL DOCUMENT FOR FURTHER
C                           DETAILS.
C                IND    - INTEGER WORK VECTOR OF LENGTH 8.
C                A      - WORK VECTOR OF LENGTH 2*N+190.
C                         NOTE THAT A COMPLETE EXPLANATION OF THE MODEL
C                           AND THE RELATION OF THE OUTPUT PARAMETERS
C                           IN THAT MODEL IS GIVEN IN THE MANUAL
C                           DOCUMENT.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER=129 INDICATES THE NORMAL EQUATIONS IN
C                             FTWENX WERE SINGULAR.
C                           IER=130 INDICATES THE NORMAL EQUATIONS IN
C                             FTARPS WERE SINGULAR.
C                           IER=131 INDICATES INSUFFICIENT CORRELATION
C                             BETWEEN SERIES X AND Y, OR THAT
C                             IB IS GREATER THAN 20.
C                           IER=132 INDICATES N WAS LESS THAN 41.
C
C   REQD. IMSL ROUTINES - SINGLE/FTARPS,FTAUTO,FTCRXY,FTFREQ,FTMA,
C                           FTMA1,FTML,FTWENX,LEQT1F,LUDATN,LUELMN,
C                           UERSET,UERTST,UGETIO,VABMXF,VBLA=SNRM2,
C                           ZSPOW,ZSPWA,ZSPWB,ZSPWC,ZSPWD,ZSPWE,ZSPWF,
C                           ZSPWG
C                       - DOUBLE/FTARPS,FTAUTO,FTCRXY,FTFREQ,FTMA,
C                           FTMA1,FTML,FTWENX,LEQT1F,LUDATN,LUELMN,
C                           UERSET,UERTST,UGETIO,VABMXF,VBLA=DNRM2,
C                           VXADD,VXMUL,VXSTO,ZSPOW,ZSPWA,ZSPWB,ZSPWC,
C                           ZSPWD,ZSPWE,ZSPWF,ZSPWG
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS      IF IER=131 IS OBSERVED BECAUSE IB IS GREATER THAN
C                20, SHIFTING THE SERIES APPROPRIATELY AND PADDING
C                WITH THE MEAN VALUES IS ONE RECOURSE.
C
C----------------------------------------------------------------------
C
      SUBROUTINE FTTR (Y,X,N,T,IB,IND,A,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N,IB,IER,IND(1)
      REAL               Y(1),X(1),T(1),A(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,I1,I2,I5,IPN,ITA,ITABM,ITABP,ITAU,ITB,IUP,J,
     *                   JER,K,KN2129,KN2169,KP1,LEVEL,LEVOLD,N2,N2P128,
     *                   N2P129,N2P169,N2P48,N2P8,N2P88,N2P9,NM1
      REAL               PMAC,RATIO,WNV,Z(1)
      DOUBLE PRECISION   T1,T2,T3,T4,T5
      DATA               I1 /1/,I2 /2/,I5 /5/
C                                  FIRST EXECUTABLE STATEMENT
      LEVEL = 0
      CALL UERSET(LEVEL,LEVOLD)
      IF (N.GT.40) GO TO 5
      IER = 132
      GO TO 9000
C                                  USEFUL CONSTANTS
    5 IER = 0
      N2 = N+N
      N2P8 = N2+8
      N2P9 = N2+9
      N2P48 = N2+48
      N2P88 = N2+88
      N2P128 = N2+128
      N2P129 = N2+129
      N2P169 = N2+169
C                                  REMOVE MEANS FROM SERIES
      T1 = 0.D0
      T2 = 0.D0
      DO 10 I=1,N
         T1 = T1+X(I)
         T2 = T2+Y(I)
   10 CONTINUE
      T1 = T1/N
      T2 = T2/N
      T(8) = T2
      T(9) = T1
      DO 15 I=1,N
         X(I) = X(I)-T1
         Y(I) = Y(I)-T2
   15 CONTINUE
C                                  CALCULATE ARMA(3,2) MODEL FOR X
      IND(1) = N
      IND(2) = 3
      IND(3) = 2
      IND(4) = 0
      IND(5) = 30
      IND(6) = 6
      IND(7) = 0
      IND(8) = 5
      CALL FTML(X,IND,T,T(4),PMAC,WNV,A,A(11),IER)
      IF (IER.EQ.130) GO TO 9000
      IER = 0
C                                  PREWHITEN X SERIES
C                                  TRANSFORM Y SERIES
      DO 20 I=1,3
         A(I) = 0.0
         A(I+N) = 0.0
   20 CONTINUE
      DO 35 I=4,N
         T1 = X(I)
         T2 = Y(I)
         DO 25 J=1,3
            T3 = T(J)
            T1 = T1-T3*X(I-J)
            T2 = T2-T3*Y(I-J)
   25    CONTINUE
         IPN = I+N
         DO 30 J=1,2
            T3 = T(J+3)
            T1 = T1+T3*A(I-J)
            T2 = T2+T3*A(IPN-J)
   30    CONTINUE
         A(I) = T1
         A(IPN) = T2
   35 CONTINUE
C                                  CALCULATE CROSS COVARINACES
      IND(1) = 0
      IND(2) = N
      IND(3) = -1
      IND(4) = 20
      IND(5) = 0
      IND(6) = 0
      A(N2+1) = 1.0
      A(N2+2) = 0.0
      CALL FTFREQ(A,IND,A(N2+1),A(N2+3),A(N2+9),Z,Z,A(N2+89),Z,Z,Z,Z,
     *JER)
      IND(1) = 1
      IND(4) = 40
      CALL FTFREQ(A,IND,A(N2+1),A(N2+3),A(N2+9),Z,Z,A(N2+89),Z,Z,Z,Z,
     *JER)
C                                  NORMALIZE COVARIANCES
      T1 = 1.0D0/A(N2+4)
      T2 = 1.0D0/A(N2+6)
      T3 = DSQRT(T1*T2)
      A(N2P129) = A(N2P129)*T3
      DO 40 I=1,40
         IF (I.LE.20) A(N2P8+I) = A(N2P8+I)*T1
         A(N2P48+I) = A(N2P48+I)*T2
         A(N2P88+I) = A(N2P88+I)*T3
         A(N2P129+I) = A(N2P129+I)*T3
   40 CONTINUE
C                                  VARIANCES OF IMPULSE PARAMETERS
      DO 45 I=1,20
         J = 21-I
         A(N2P9+J) = A(N2P8+J)
   45 CONTINUE
      A(N2P9) = 1.0
      A(N2P48) = 1.0
      K = 0
      DO 55 KP1=1,21
         KN2129 = K+N2P129
         KN2169 = KP1+N2P169
         A(KN2169) = 0.0
         DO 50 I=1,41
            ITAU = I-21
            ITA = IABS(ITAU)+N2P9
            ITB = ITA+39
            ITABP = KN2129+ITAU
            ITABM = KN2129-ITAU
            A(KN2169) = A(KN2169)+(A(ITA)*A(ITB)+A(ITABP)*A(ITABM)
     *      +A(KN2129)**2*(A(ITAU+N2P129)**2+0.5*(A(ITA)**2+A(ITB)**2))
     *      -(A(KN2129)+A(KN2129))*(A(ITA)*A(ITABP)+A(-ITAU+N2P129)*
     *      A(K+ITB)))
   50    CONTINUE
         A(KN2169) = A(KN2169)/(N-K)
         K = KP1
   55 CONTINUE
C                                  FIND FIRST SIGNIFICANT IMPULSE PARA-
C                                  METER AT 2*STANDARD DEVIATION LEVEL
      RATIO = A(N2+6)/A(N2+4)
      IB = 21
      DO 60 I=1,21
         IF (A(I+N2P128)**2*RATIO/A(I+N2P169).GE.4.0) IB = MIN0(I-1,IB)
   60 CONTINUE
      IF (IB.EQ.21) IER = 131
      IF (IER.EQ.131) GO TO 9000
C                                  CALCULATE MODEL PARAMETER ESTIMATES
C                                  OF Y LAGGING X BY IB TIME UNITS
      NM1 = N-1
      DO 65 I=2,N
         A(I) = Y(I)
         A(I+NM1) = X(I)
   65 CONTINUE
      A(1) = Y(1)
      A(N2) = 0.0
      IND(1) = 2
      IND(2) = 3
      IND(3) = 1
      IND(4) = IB+1
      IND(5) = -1
      IND(6) = -1
      CALL FTWENX(A,N,I2,N,IND,IND(3),IND(5),A(N2+1),I5,T,PMAC,A(N2+26),
     *IER)
      IF (IER.NE.0) GO TO 9000
C                                  CALCULATE RESIDUAL NOISE SERIES
      IUP = IB+3
      A(IUP-1) = 0.0
      A(IUP-2) = 0.0
      T1 = T(1)
      T2 = T(2)
      T3 = T(3)
      T4 = T(4)
      T5 = T(5)
      DO 70 I=IUP,N
         K = I-IB
         A(I) = Y(I)-T1*(DBLE(Y(I-1))-A(I-1))-T2*(DBLE(Y(I-2))-A(I-2))
     *   -T3*X(K)-T4*X(K-1)-T5*X(K-2)
   70 CONTINUE
      CALL FTWENX(A(IUP),N,I1,N-IUP+1,IND,IND(3),IND(5),A(N+1),I2,T(6),
     *PMAC,A(N+5),IER)
      IF (IER.NE.0) GO TO 9000
C                                  SUM OF SQUARES OF WHITE NOISE
      T1 = T(6)
      T2 = T(7)
      T3 = 0.0D0
      DO 75 I=IUP,N
         T3 = T3+(A(I)-T1*A(I-1)-T2*A(I-2))**2
   75 CONTINUE
      T(10) = T3
      CALL UERSET(LEVOLD,LEVEL)
      IF (IER.EQ.0) GO TO 9005
 9000 CONTINUE
      CALL UERSET(LEVOLD,LEVEL)
      CALL UERTST(IER,'FTTR  ')
 9005 CONTINUE
      RETURN
      END
