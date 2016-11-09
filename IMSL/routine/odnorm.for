
 
C   IMSL ROUTINE NAME   - ODNORM
C
C-----------------------------------------------------------------------
C
C   COMPUTER            - IBM/SINGLE
C
C   LATEST REVISION     - JUNE 1, 1982
C
C   PURPOSE             - MULTIVARIATE NORMAL LINEAR DISCRIMINANT
C                           ANALYSIS AMONG SEVERAL KNOWN GROUPS
C
C   USAGE               - CALL ODNORM (X,IX,NG,NV,ND,P,
C                                      XM,IXM,S,W,IW,D,ID,IER)
C
C   ARGUMENTS    X      - INPUT MATRIX OF DIMENSION NV BY THE SUM OF
C                           ND(I) CONTAINING THE GROUP DATA.
C                           (SEE REMARK 1).
C                IX     - INPUT ROW DIMENSION OF MATRIX X EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM.
C                NG     - INPUT NUMBER OF DISTINCT GROUPS.
C                NV     - INPUT NUMBER OF DATA VARIABLES.
C                ND     - INPUT VECTOR OF LENGTH NG.  ND(I) CONTAINS
C                           THE NUMBER OF DATA VECTORS FROM THE I-TH
C                           GROUP IN THE DATA MATRIX X.
C                P      - INPUT/OUTPUT VECTOR OF LENGTH NG. ON INPUT,
C                           P(I) CONTAINS THE BAYESIAN PRIOR PROBABILITY
C                           THAT A RANDOM OBSERVATION BELONGS TO THE
C                           I-TH GROUP. (SEE REMARK 2).
C                         ON OUTPUT, P CONTAINS THE CONSTANTS OF THE
C                           NG LINEAR DISCRIMINANT FUNCTIONS.
C                XM     - OUTPUT MATRIX OF DIMENSION NV BY NG.
C                           THE I-TH COLUMN OF XM CONTAINS THE
C                           MEAN VECTOR FOR THE I-TH GROUP.
C                IXM    - INPUT ROW DIMENSION OF MATRIX XM EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM.
C                S      - WORK VECTOR OF LENGTH ((NV+1)*NV)/2.
C                W      - OUTPUT MATRIX OF DIMENSION NV BY NG.  THE
C                           I-TH COLUMN OF W CONTAINS THE VECTOR FOR
C                           THE I-TH LINEAR DISCRIMINANT FUNCTION.
C                IW     - INPUT ROW DIMENSION OF MATRIX W EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM.
C                D      - OUTPUT PERFORMANCE MATRIX OF DIMENSION NG BY
C                           NG. (SEE REMARK 3).
C                ID     - INPUT ROW DIMENSION OF MATRIX D EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM.
C                IER    - OUTPUT ERROR PARAMETER
C                         TERMINAL ERROR
C                           IER = 129 INDICATES THAT THE COVARIANCE
C                             MATRIX IS SINGULAR. CHECK THE VARIABLES
C                             FOR REDUNDANCIES.
C                           IER = 130 INDICATES ONE OF IX,IXM,IW, OR ID
C                             IS TOO SMALL.  REDIMENSION MATRICES.
C
C   PRECISION/HARDWARE  - SINGLE AND DOUBLE/H32
C                       - SINGLE/H36,H48,H60
C
C   REQD. IMSL ROUTINES - LEQT1P,LUDECP,LUELMP,UERTST,UGETIO
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
C            2.  COMMON CHOICES FOR THE BAYESIAN PRIOR PROBABILITIES
C                P(I) ARE 1/NG, ND(I)/(ND(1)+ND(2)+...+ND(NG)), OR
C                SUBJECTIVE ESTIMATES. THESE SHOULD SUM TO ONE.
C            3.  THE RESULTING BAYESIAN CLASSIFIER IS TESTED ON THE
C                ORIGINAL DATA X. THE OUTPUT MATRIX D SUMMARIZES THE
C                PERFORMANCE. D(I,J) IS THE FRACTION OF SAMPLES IN THE
C                I-TH GROUP THAT WERE CLASSIFIED IN THE J-TH GROUP.
C                D(I,I) IS FRACTION CLASSIFIED FOR GROUP I. THE ACTUAL
C                PERFORMANCE OF THE BAYESIAN CLASSIFIER ON UNLABELED
C                OBSERVATIONS WILL BE SLIGHTLY WORSE.
C
C   COPYRIGHT           - 1978 BY IMSL, INC. ALL RIGHTS RESERVED.
C
C   WARRANTY            - IMSL WARRANTS ONLY THAT IMSL TESTING HAS BEEN
C                           APPLIED TO THIS CODE. NO OTHER WARRANTY,
C                           EXPRESSED OR IMPLIED, IS APPLICABLE.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE ODNORM (X,IX,NG,NV,ND,P,XM,IXM,S,W,IW,D,ID,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            IX,NG,NV,ND(1),IXM,IW,ID,IER
      REAL               X(IX,1),P(1),XM(IXM,1),S(1),W(IW,1),D(ID,1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,IDGT,IR,J,K,KHI,KLO,L,LMAX,N,NDI,NDJ,NMNG
      REAL               D1,D2,XMAX,XMI,XMJ
      DOUBLE PRECISION   DSUM,DZERO
      DATA               DZERO/0.0D0/
C
C                       CHECK SUFFICIENT MATRIX DIMENSIONS
C
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
      IF (IX.GE.NV.AND.IXM.GE.NV.AND.IW.GE.NV.AND.ID.GE.NG) GO TO 5
      IER = 130
      GO TO 9000
C                       COMPUTE MATRIX OF MEAN VECTORS
    5 KHI = 0
      DO 20 J=1,NG
         NDJ = ND(J)
         KLO = KHI+1
         KHI = KLO+NDJ-1
         DO 15 I=1,NV
            DSUM = DZERO
            DO 10 K=KLO,KHI
               DSUM = DSUM+X(I,K)
   10       CONTINUE
            XM(I,J) = DSUM/NDJ
   15    CONTINUE
   20 CONTINUE
      N = KHI
C                       COMPUTE POOLED COVARIANCE MATRIX
      NMNG = N-NG
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
            S(IR) = DSUM/NMNG
   35    CONTINUE
   40 CONTINUE
      DO 50 J=1,NG
         DO 45 I=1,NV
            W(I,J) = XM(I,J)
   45    CONTINUE
   50 CONTINUE
C                       SOLVE FOR DISCRIMINANT WEIGHT VECTORS
      IDGT = 0
      CALL LEQT1P (S,NG,NV,W,IW,IDGT,D1,D2,IER)
      IF (IER.NE.0) GO TO 9000
C                       COMPUTE CONSTANT WEIGHTS FOR EACH GROUP
      DO 60 I=1,NG
         DSUM = DZERO
         DO 55 J=1,NV
            DSUM = DSUM+XM(J,I)*W(J,I)
   55    CONTINUE
         P(I) = -DSUM/2+ALOG(P(I))
   60 CONTINUE
C                       COMPUTE PERFORMANCE MATRIX
      KHI = 0
      DO 90 I=1,NG
         NDI = ND(I)
         DO 65 J=1,NG
            D(I,J) = 0
   65    CONTINUE
         KLO = KHI+1
         KHI = KLO+NDI-1
         DO 80 K=KLO,KHI
            XMAX = -1.E20
            DO 75 L=1,NG
               DSUM = DZERO
               DO 70 J=1,NV
                  DSUM = DSUM+W(J,L)*X(J,K)
   70          CONTINUE
               DSUM = DSUM+P(L)
               IF (DSUM.LT.XMAX) GO TO 75
               XMAX = DSUM
               LMAX = L
   75       CONTINUE
            D(I,LMAX) = D(I,LMAX)+1
   80    CONTINUE
         DO 85 J=1,NG
            D(I,J) = D(I,J)/NDI
   85    CONTINUE
   90 CONTINUE
      GO TO 9005
 9000 CONTINUE
      CALL UERTST (IER,6HODNORM)
 9005 RETURN
      END
 
R; T=0.05/0.53 22:54:33
INUE
 