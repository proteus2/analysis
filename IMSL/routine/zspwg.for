C   IMSL ROUTINE NAME   - ZSPWG
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1982
C
C   PURPOSE             - NUCLEUS CALLED ONLY BY IMSL SUBROUTINE ZSPOW
C
C   REQD. IMSL ROUTINES - SINGLE/VBLA=SNRM2
C                       - DOUBLE/VBLA=DNRM2
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE ZSPWG (M,N,A,LDA,PIVOT,IPVT,LIPVT,RDIAG,ACNORM,WA)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            M,N,LDA,LIPVT,IPVT(LIPVT)
      REAL               A(LDA,N),RDIAG(N),ACNORM(N),WA(N)
      LOGICAL            PIVOT
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,JP1,J,KMAX,K,MINMN
      REAL               AJNORM,EPSMCH,ONE,P05,SPMPAR,SUM,TEMP,ZERO
      REAL               SNRM2
      DATA               SPMPAR /Z3C100000/
      DATA               ONE,P05,ZERO /1.0E0,5.0E-2,0.0E0/
C                                  EPSMCH IS THE MACHINE PRECISION.
C                                  FIRST EXECUTABLE STATEMENT
      EPSMCH = SPMPAR
C                                  COMPUTE THE INITIAL COLUMN NORMS AND
C                                  INITIALIZE SEVERAL ARRAYS.
      DO 5 J=1,N
         ACNORM(J) = SNRM2(M,A(1,J),1)
         RDIAG(J) = ACNORM(J)
         WA(J) = RDIAG(J)
         IF (PIVOT) IPVT(J) = J
    5 CONTINUE
C                                  REDUCE A TO R WITH HOUSEHOLDER
C                                  TRANSFORMATIONS.
      MINMN = MIN0(M,N)
      DO 55 J=1,MINMN
         IF (.NOT.PIVOT) GO TO 20
C                                  BRING THE COLUMN OF LARGEST NORM INTO
C                                  THE PIVOT POSITION.
         KMAX = J
         DO 10 K=J,N
            IF (RDIAG(K).GT.RDIAG(KMAX)) KMAX = K
   10    CONTINUE
         IF (KMAX.EQ.J) GO TO 20
         DO 15 I=1,M
            TEMP = A(I,J)
            A(I,J) = A(I,KMAX)
            A(I,KMAX) = TEMP
   15    CONTINUE
         RDIAG(KMAX) = RDIAG(J)
         WA(KMAX) = WA(J)
         K = IPVT(J)
         IPVT(J) = IPVT(KMAX)
         IPVT(KMAX) = K
   20    CONTINUE
C                                  COMPUTE THE HOUSEHOLDER
C                                  TRANSFORMATION TO REDUCE THE J-TH
C                                  COLUMN OF A TO A MULTIPLE OF THE J-TH
C                                  UNIT VECTOR.
         AJNORM = SNRM2(M-J+1,A(J,J),1)
         IF (AJNORM.EQ.ZERO) GO TO 50
         IF (A(J,J).LT.ZERO) AJNORM = -AJNORM
         DO 25 I=J,M
            A(I,J) = A(I,J)/AJNORM
   25    CONTINUE
         A(J,J) = A(J,J)+ONE
C                                  APPLY THE TRANSFORMATION TO THE
C                                  REMAINING COLUMNS AND UPDATE THE
C                                  NORMS.
         JP1 = J+1
         IF (N.LT.JP1) GO TO 50
         DO 45 K=JP1,N
            SUM = ZERO
            DO 30 I=J,M
               SUM = SUM+A(I,J)*A(I,K)
   30       CONTINUE
            TEMP = SUM/A(J,J)
            DO 35 I=J,M
               A(I,K) = A(I,K)-TEMP*A(I,J)
   35       CONTINUE
            IF (.NOT.PIVOT .OR. RDIAG(K).EQ.ZERO) GO TO 40
            TEMP = A(J,K)/RDIAG(K)
            RDIAG(K) = RDIAG(K)*SQRT(AMAX1(ZERO,ONE-TEMP**2))
            IF (P05*(RDIAG(K)/WA(K))**2.GT.EPSMCH) GO TO 40
            RDIAG(K) = SNRM2(M-J,A(JP1,K),1)
            WA(K) = RDIAG(K)
   40       CONTINUE
   45    CONTINUE
   50    CONTINUE
         RDIAG(J) = -AJNORM
   55 CONTINUE
      RETURN
      END
