C   IMSL ROUTINE NAME   - MSENO
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1982
C
C   PURPOSE             - EXPECTED VALUES OF NORMAL ORDER STATISTICS
C
C   USAGE               - CALL MSENO (IFIRST,ILAST,N,EX,IER)
C
C   ARGUMENTS    IFIRST - INPUT.  FIRST (SMALLEST) ORDER STATISTIC
C                           WHOSE EXPECTED VALUE IS TO BE COMPUTED.
C                ILAST  - INPUT.  LAST ORDER STATISTIC WHOSE EXPECTED
C                           VALUE IS TO BE COMPUTED.
C                N      - INPUT.  SAMPLE SIZE.
C                EX     - OUTPUT.  VECTOR OF LENGTH (ILAST + 1 - IFIRST)
C                           CONTAINING THE EXPECTED VALUES OF THE IFIRST
C                           THROUGH THE ILAST ORDER STATISTICS FROM A
C                           STANDARD NORMAL DISTRIBUTION.
C                IER    - ERROR PARAMETER.  (OUTPUT)
C                         WARNING WITH FIX ERROR
C                           IER = 65  INDICATES IFIRST IS LESS THAN 1
C                             OR ILAST IS GREATER THAN N.  THE COM-
C                             PUTATIONS PROCEED AS IF THE PARAMETER(S)
C                             OUT OF RANGE HAD BEEN AT THE APPROPRIATE
C                             LIMIT(S) OF THE RANGE.
C                         TERMINAL ERROR
C                           IER = 129  INDICATES IFIRST IS GREATER
C                             THAN ILAST.
C
C   REQD. IMSL ROUTINES - H32,H36/MLGAMD=DLGAMA,UERTST,UGETIO
C                         H48,H60/MLGAMA=ALGAMA,UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE MSENO (IFIRST,ILAST,N,EX,IER)
C
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            IFIRST,ILAST,N,IER
      REAL               EX(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,IST,IOE,J,JM,JP,J2
      INTEGER            K,KK,L,LST,M,M1,M2
      REAL               ZERO,MSENP
      DATA               ZERO /0.0/
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
      IST = IFIRST
      LST = ILAST
      IF(IFIRST.GT.ILAST) GO TO 95
      IF(ILAST.LE.N) GO TO 5
      IF(IFIRST.GT.N) IST = N
      IER = 65
      LST = N
    5 IF(IFIRST.GE.1) GO TO 10
      IF(ILAST.LT.1) LST = 1
      IER = 65
      IST = 1
   10 IF(IER.EQ.0) GO TO 15
      CALL UERTST(IER,6HMSENO )
   15 CONTINUE
      IOE = MOD(N,2)
      J2 = LST-IST+1
      M = (N+1)/2
      IF(IST-M) 20, 80, 85
   20 IF(LST-M) 70, 65, 25
C                                  COMPUTE EXPECTATIONS FOR CASES WHEN
C                                    ORDER STATISTICS ARE FROM BOTH
C                                    NEGATIVE AND POSITIVE SIDES OF
C                                    THE NORMAL (0,1) DISTRIBUTION.
   25 M1 = M
      M2 = M+1
      JM = 1
      IF(IOE.EQ.0) GO TO 30
      EX(M-IST+1) = ZERO
      JM = 2
      M1 = M-1
   30 IF(N-IST-LST) 50, 35, 35
C                                  COMPUTE EXPECTATIONS FOR NEGATIVE
C                                    SIDE OF DISTRIBUTION AND FILL
C                                    POSITIVE SIDE FROM COMPUTED EX.
   35 J = 0
      DO 40  I=IST,M1
         J = J+1
         EX(J) = -MSENP(I,N)
   40 CONTINUE
      J = 0
      K = M2-IST+1
      KK = K-JM
      DO 45  I=M2,LST
         EX(K+J) = -EX(KK-J)
         J = J+1
   45 CONTINUE
      RETURN
C                                  COMPUTE EXPECTATIONS FOR POSITIVE
C                                    SIDE OF DISTRIBUTION AND FILL
C                                    NEGATIVE SIDE FROM COMPUTED EX.
   50 CONTINUE
      J = J2
      K = N+1-LST
      DO 55  I=K,M1
         EX(J) = MSENP(I,N)
         J = J-1
   55 CONTINUE
      J = 1
      JP = (M1-IST)*2+3+IOE
      DO 60  I=IST,M1
         EX(J) = -EX(JP-J)
         J = J+1
   60 CONTINUE
      RETURN
C                                  COMPUTE EXPECTATIONS FOR CASES WHEN
C                                    ORDER STATISTICS ARE FROM ONLY
C                                    THE NEGATIVE SIDE OF THE NORMAL
C                                    (0,1) DISTRIBUTION.
   65 IF(IOE.EQ.0) GO TO 70
      EX(J2) = ZERO
      LST = LST-1
   70 J = 0
      DO 75  I=IST,LST
         J = J+1
         EX(J) = -MSENP(I,N)
   75 CONTINUE
      RETURN
C                                  COMPUTE EXPECTATIONS FOR CASES WHEN
C                                    ORDER STATISTICS ARE FROM ONLY
C                                    THE POSITIVE SIDE OF THE NORMAL
C                                    (0,1) DISTRIBUTION.
   80 IF(IOE.EQ.0) GO TO 85
      EX(1) = ZERO
      IST = IST+1
   85 J = J2
      K = N+1-LST
      L = N+1-IST
      DO 90  I=K,L
         EX(J) = MSENP(I,N)
         J = J-1
   90 CONTINUE
      RETURN
   95 IER = 129
      CALL UERTST(IER,6HMSENO )
      RETURN
      END
