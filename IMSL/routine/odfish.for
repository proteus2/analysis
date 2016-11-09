C   IMSL ROUTINE NAME   - ODFISH
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - NOVEMBER 1, 1979
C
C   PURPOSE             - LINEAR DISCRIMINANT ANALYSIS METHOD OF FISHER
C                           FOR REDUCING THE NUMBER OF VARIABLES
C
C   USAGE               - CALL ODFISH (X,IX,NG,NV,ND,XM,IXM,NNV,E,C,IC,
C                                      SW,SB,ISB,EX,CX,ICX,IS,IER)
C
C   ARGUMENTS    X      - INPUT/OUTPUT MATRIX OF DIMENSION NV BY THE
C                           SUM OF THE ND(I). ON INPUT, X CONTAINS THE
C                           GROUP DATA.
C                         ON OUTPUT, EACH COLUMN DATA VECTOR IS REPLACED
C                           BY THE TRANSFORMED NNV BY 1 DATA VECTOR.
C                           (SEE REMARK 1).
C                IX     - INPUT ROW DIMENSION OF MATRIX X EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM.
C                NG     - INPUT NUMBER OF DISTINCT GROUPS.
C                NV     - INPUT NUMBER OF DATA VARIABLES.
C                ND     - INPUT VECTOR OF LENGTH NG.  ND(I) CONTAINS
C                           THE NUMBER OF DATA VECTORS FROM THE I-TH
C                           GROUP IN THE DATA MATRIX X.
C                XM     - OUTPUT MATRIX OF DIMENSION NV BY NG+1.
C                           THE I-TH COLUMN OF XM CONTAINS THE
C                           MEAN VECTOR FOR THE I-TH GROUP.  COLUMN
C                           NG+1 CONTAINS THE OVERALL MEAN VECTOR.
C                IXM    - INPUT ROW DIMENSION OF MATRIX XM EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM.
C                NNV    - OUTPUT NUMBER OF NEW VARIABLES COMPUTED.
C                           NNV IS GIVEN BY THE MINIMUM OF NV AND NG-1.
C                E      - OUTPUT VECTOR OF LENGTH NV CONTAINING THE
C                           EIGENVALUES IN THE FIRST NNV LOCATIONS IN
C                           DESCENDING ORDER. THE REMAINDER OF E IS
C                           USED AS WORK AREA. (SEE REMARK 2).
C                C      - OUTPUT MATRIX OF DIMENSION NV BY NNV.  THE
C                           I-TH COLUMN OF C CONTAINS THE VECTOR FOR
C                           COMPUTING THE I-TH NEW LINEAR VARIABLE.
C                           (SEE REMARK 1).
C                IC     - INPUT ROW DIMENSION OF MATRIX C EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM.
C                SW     - WORK VECTOR OF LENGTH ((NV+1)*NV)/2 OR 4,
C                           WHICHEVER IS GREATER. (SEE REMARK 3).
C                SB     - WORK MATRIX OF DIMENSION NV BY NV.
C                           (SEE REMARK 3).
C                ISB    - INPUT ROW DIMENSION OF MATRIX SB EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM.
C                EX     - COMPLEX WORK VECTOR OF LENGTH NV CONTAINING
C                           THE UNSORTED EIGENVALUES.  (SEE REMARK 4).
C                         NOTE - THE ROUTINE TREATS EX AS A REAL VECTOR
C                           OF LENGTH 2*NV. AN APPROPRIATE EQUIVALENCE
C                           STATEMENT MAY BE REQUIRED. SEE DOCUMENT
C                           EXAMPLE.
C                CX     - COMPLEX WORK MATRIX OF DIMENSION NV BY NV
C                           CONTAINING THE EIGENVECTORS ASSOCIATED WITH
C                           EX.  (SEE REMARK 4).
C                         NOTE - THE ROUTINE TREATS CX AS A REAL VECTOR
C                           OF LENGTH 2*NV*NV. AN APPROPRIATE EQUI-
C                           VALENCE STATEMENT MAY BE REQUIRED. SEE
C                           DOCUMENT EXAMPLE.
C                ICX    - INPUT ROW DIMENSION OF MATRIX CX EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM.
C                IS     - INTEGER WORK VECTOR OF LENGTH NV.
C                IER    - OUTPUT ERROR PARAMETER
C                         TERMINAL ERROR
C                           IER = 129 INDICATES THAT THE COVARIANCE
C                             MATRIX IS SINGULAR.  CHECK FOR
C                             REDUNDANCIES IN THE ORIGINAL VARIABLES.
C                           IER = 130 INDICATES THAT THE EIGENVALUES
C                             COULD NOT BE COMPUTED.  REDUCE THE
C                             NUMBER OF VARIABLES AND TRY AGAIN.
C                           IER = 131 INDICATES ONE OF IX,IXM,IC,ISB, OR
C                             ICX IS TOO SMALL.  REDIMENSION MATRICES.
C
C
C   REQD. IMSL ROUTINES - SINGLE/EBALAF,EBBCKF,EHBCKF,EHESSF,EIGRF,
C                           EQRH3F,LEQT1P,LUDECP,LUELMP,UERTST,UGETIO,
C                           VSRTR
C                       - DOUBLE/EBALAF,EBBCKF,EHBCKF,EHESSF,EIGRF,
C                           EQRH3F,LEQT1P,LUDECP,LUELMP,UERTST,UGETIO,
C                           VSRTRD
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS  1.  ON INPUT, THE MATRIX X CONTAINS THE ND(1)+ND(2)+...
C                ...+ND(NG) SAMPLE (COLUMN) VECTORS. THE FIRST ND(1)
C                COLUMNS OF MATRIX X SHOULD CONTAIN THE DATA FOR GROUP
C                ONE. THE NEXT ND(2) COLUMNS OF MATRIX X SHOULD CONTAIN
C                THE DATA FOR GROUP TWO, AND SO ON.
C                ON OUTPUT, EACH ORIGINAL OBSERVATION X(I) (WHERE X(I)
C                REPRESENTS A COLUMN VECTOR OF INPUT MATRIX X) IS
C                REPLACED BY ITS TRANSFORMED OBSERVATION ( A NNV BY 1
C                DATA VECTOR). THE TRANSFORMED OBSERVATION IS GIVEN BY
C                THE FORMULA (C TRANSPOSE) * X(I) WHERE C IS THE OUTPUT
C                MATRIX C.
C            2.  THE COLUMN VECTORS IN THE TRANSFORMATION MATRIX C ARE
C                IN ORDER OF DECREASING IMPORTANCE, THAT IS, THE EIGEN-
C                VECTOR ASSOCIATED WITH THE LARGEST EIGENVALUE IS IN
C                THE FIRST COLUMN OF C, AND SO ON.
C            3.  IF THE WITHIN AND BETWEEN SUMS OF SQUARES MATRICES, SW
C                AND SB, ARE DESIRED, THEY MAY BE PRINTED BY MODIFYING
C                THE SOURCE CODE AT THE DESIGNATED COMMENT CARD.
C            4.  ONE SHOULD PRUDENTLY INSPECT THE COMPLEX EIGENVALUE
C                VECTOR AND EIGENVECTOR MATRIX AND CHECK THAT THERE ARE
C                ONLY NG-1 OR NV NONZERO EIGENVALUES AND THAT THE IMAG-
C                INARY PARTS OF THE EIGENVECTORS ARE TRULY ZERO.
C
C-----------------------------------------------------------------------
C
        SUBROUTINE ODFISH (X,IX,NG,NV,ND,XM,IXM,NNV,E,C,IC,
     *                    SW,SB,ISB,EX,CX,ICX,IS,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            IX,NG,NV,ND(1),IXM,NNV,IC,ISB,ICX,IS(1),IER
      REAL               X(IX,1),XM(IXM,1),E(1),C(IC,1),SW(1),SB(ISB,1)
      REAL               EX(2,1),CX(2,ICX,1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,IDGT,IJOB,J,K,KHI,KLO,L,N,NGP1
      REAL               D1,D2,XMI,XMJ
      DOUBLE PRECISION   DSUM,DTOT,DZERO
      DATA               DZERO/0.0D0/
C                       ERROR CHECKS
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
      IF (IX.GE.NV.AND.IXM.GE.NV.AND.IC.GE.NV.AND.ISB.GE.NV.AND.ICX.GE.N
     1V) GO TO 5
      IER = 131
      GO TO 9000
C                       COMPUTE MEAN VECTOR MATRIX
    5 NGP1 = NG+1
      DO 20 I=1,NV
         DTOT = DZERO
         KHI = 0
         DO 15 J=1,NG
            KLO = KHI+1
            KHI = KLO+ND(J)-1
            DSUM = DZERO
            DO 10 K=KLO,KHI
               DSUM = DSUM+X(I,K)
   10       CONTINUE
            XM(I,J) = DSUM/ND(J)
            DTOT = DTOT+DSUM
   15    CONTINUE
         XM(I,NGP1) = DTOT/KHI
   20 CONTINUE
      N = KHI
C                       COMPUTE POOLED SUMS OF SQUARES MATRIX
      IR = 0
      DO 40 I=1,NV
         DO 35 J=1,I
            DSUM = DZERO
            KHI = 0
            DO 30 L=1,NG
               KLO = KHI+1
               KHI = KLO+ND(L)-1
               XMI = XM(I,L)
               XMJ = XM(J,L)
               DO 25 K=KLO,KHI
                  DSUM = DSUM+(X(I,K)-XMI)*(X(J,K)-XMJ)
   25          CONTINUE
   30       CONTINUE
            IR = IR+1
            SW(IR) = DSUM
   35    CONTINUE
   40 CONTINUE
C                       PRINT WITHIN SUMS OF SQUARES MATRIX HERE
C                       COMPUTE BETWEEN SUMS OF SQUARES MATRIX
      DO 55 I=1,NV
         DO 50 J=1,I
            XMI = XM(I,NGP1)
            XMJ = XM(J,NGP1)
            DSUM = DZERO
            DO 45 K=1,NG
               DSUM = DSUM+ND(K)*(XM(I,K)-XMI)*(XM(J,K)-XMJ)
   45       CONTINUE
            SB(I,J) = DSUM
            SB(J,I) = DSUM
   50    CONTINUE
   55 CONTINUE
C                       PRINT BETWEEN SUMS OF SQUARES MATRIX HERE
C                       FORM MATRIX TO GET EIGENVALUES AND VECTORS
      IDGT = 0
      CALL LEQT1P (SW,NV,NV,SB,ISB,IDGT,D1,D2,IER)
      IF (IER.NE.0) GO TO 9000
C                       COMPUTE EIGEN INFORMATION OF SB
      IJOB = 1
      CALL EIGRF (SB,NV,ISB,IJOB,EX,CX,ICX,SW,IER)
      IF (IER.EQ.0) GO TO 60
      IER = 130
      GO TO 9000
C                       EXTRACT REAL EIGENSTRUCTURE
   60 DO 65 I=1,NV
         E(I) = -EX(1,I)
         IS(I) = I
   65 CONTINUE
C                       SORTED INTO DESCENDING ORDER (BY NEGATIVE)
      CALL VSRTR (E,NV,IS)
C                       THE NUMBER OF RETAINED EIGENVALUES = NNV
      NNV = MIN0(NV,NG-1)
      DO 75 J=1,NNV
         E(J) = -E(J)
         K = IS(J)
         DO 70 I=1,NV
            C(I,J) = CX(1,I,K)
   70    CONTINUE
   75 CONTINUE
C                       COMPUTE FISHER PROJECTION VARIABLES
      DO 95 J=1,N
         DO 85 K=1,NNV
            DSUM = DZERO
            DO 80 I=1,NV
               DSUM = DSUM+X(I,J)*C(I,K)
   80       CONTINUE
            SW(K) = DSUM
   85    CONTINUE
         DO 90 K=1,NNV
            X(K,J) = SW(K)
   90    CONTINUE
   95 CONTINUE
      GO TO 9005
 9000 CONTINUE
      CALL UERTST (IER,6HODFISH)
 9005 RETURN
      END
