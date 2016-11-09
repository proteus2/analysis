C   IMSL ROUTINE NAME   - ACTRST
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1980
C
C   PURPOSE             - CONTRAST ESTIMATES AND SUMS OF SQUARES
C
C   USAGE               - CALL ACTRST (T,NR,N,ID,P,IP,Q,SQ)
C
C   ARGUMENTS    T      - INPUT VECTOR OF LENGTH N CONTAINING THE SAMPLE
C                           MEANS FOR THE LEVELS OF A CLASSIFICATION
C                           VARIABLE.
C                NR     - INPUT NUMBER OF RESPONSES IN EACH SAMPLE MEAN.
C                N      - INPUT NUMBER OF SAMPLE MEANS INVOLVED IN THE
C                           CONTRASTS.
C                ID     - INPUT NUMBER OF CONTRASTS TO BE ESTIMATED.
C                P      - INPUT N BY ID MATRIX CONTAINING THE CONTRAST
C                           COEFFICIENTS FOR EACH OF THE ID CONTRASTS.
C                IP     - INPUT ROW DIMENSION OF MATRIX P EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM.
C                Q      - OUTPUT VECTOR OF LENGTH ID CONTAINING CONTRAST
C                           ESTIMATES.
C                SQ     - OUTPUT VECTOR OF LENGTH ID CONTAINING THE SUMS
C                           OF SQUARES ASSOCIATED WITH EACH CONTRAST.
C
C   REQD. IMSL ROUTINES - SINGLE/NONE REQUIRED
C                       - DOUBLE/VXADD,VXMUL,VXSTO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS  1.  THE ORDERING OF COMPONENTS IN VECTORS Q AND SQ
C                CORRESPONDS TO THE ORDERING OF CONTRAST COEFFICIENTS
C                IN THE MATRIX P.
C            2.  WHEN CONTRAST COEFFICIENTS (LINEAR,QUADRATIC,CUBIC,
C                AND SO FORTH) ARE DESIRED THROUGH DEGREE ID, IMSL
C                ROUTINE RLPOL MAY BE USED TO GENERATE THE MATRIX P.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE ACTRST (T,NR,N,ID,P,IP,Q,SQ)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            NR,N,ID,IP
      REAL               T(N),P(IP,ID),Q(ID),SQ(ID)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,J
      REAL               XNR
      DOUBLE PRECISION   S,V
C                                  FIRST EXECUTABLE STATEMENT
      XNR = NR
      DO 10 I = 1,ID
C                                  COMPUTE ESTIMATE OF THE I-TH CONTRAST
         S = 0.0D0
         DO 5 J = 1,N
            S = S+DBLE(P(J,I))*DBLE(T(J))
    5    CONTINUE
         Q(I) = S
   10 CONTINUE
      DO 20 I = 1,ID
C                                  COMPUTE SUM OF SQUARES ATTRIBUTABLE
C                                  TO THE I-TH CONTRAST
         S = 0.0D0
         DO 15 J = 1,N
            V = DBLE(P(J,I))
            S = S+V*V
   15    CONTINUE
         SQ(I) = Q(I)*Q(I)*XNR/S
   20 CONTINUE
      RETURN
      END
