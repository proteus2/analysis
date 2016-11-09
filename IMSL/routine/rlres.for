C   IMSL ROUTINE NAME   - RLRES
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - PERFORM A RESIDUAL ANALYSIS FOR A FITTED
C                           REGRESSION MODEL
C
C   USAGE               - CALL RLRES (XY,IX,MM,N,IH,M,BETA,SDR,RES,
C                           IR,IER)
C
C   ARGUMENTS    XY     - INPUT N BY MM MATRIX OF N DATA POINTS USED
C                           TO FIT A MULTIPLE REGRESSION MODEL.
C                           TYPICAL USAGE WOULD IMPLY THAT ALL MM
C                           VARIABLES ENTERED THE MODEL WITH THE
C                           DEPENDENT (RESPONSE) VARIABLE SETTINGS
C                           APPEARING IN COLUMN MM.  IN THIS CASE
C                           M = MM, AND THE FIRST M ELEMENTS OF IH ARE
C                           SET TO 1,2,...,M.  MORE GENERALLY, THE
C                           FITTED MODEL MAY BE BASED ON A REORDERED
C                           SUBSET OF THE VARIABLES STORED IN XY.
C                           SEE PARAMETER IH.
C                IX     - INPUT ROW DIMENSION OF THE MATRIX XY EXACTLY
C                           AS SPECIFIED IN THE DIMENSION STATEMENT IN
C                           THE CALLING PROGRAM.
C                MM     - INPUT NUMBER OF COLUMNS IN XY (NUMBER OF
C                           VARIABLES IN THE ORIGINAL DATA MATRIX).
C                N      - INPUT NUMBER OF ROWS IN XY.
C                           N MUST BE GREATER THAN OR EQUAL TO 1.
C                IH     - INPUT NON-NEGATIVE VECTOR OF LENGTH M + MM
C                           INDICATING WHICH VARIABLES (AND THEIR ORDER)
C                           WERE USED IN DETERMINING THE PARAMETER
C                           ESTIMATES IN BETA.  FOR TYPICAL USAGE
C                           IH IS SPECIFIED AS STATED IN THE DESCRIPTION
C                           OF PARAMETER XY.  HOWEVER, IF IMSL ROUTINE
C                           RLSUM IS USED TO SELECT A SUBSET PROBLEM
C                           FOR MODEL FITTING, IH SHOULD BE SPECIFIED
C                           AS IS INPUT IH IN THE RLSUM DOCUMENT.
C                           (IN THIS CASE IH(I) = J IMPLIES VARIABLE J
C                           OF XY APPEARED AS VARIABLE I IN THE FITTED
C                           MODEL, AND IH(M) ALWAYS SPECIFIES THE
C                           RESPONSE VARIABLE.)  THE REMAINING
C                           LOCATIONS ARE WORK STORAGE.
C                M      - INPUT NUMBER, NOT EXCEEDING MM, OF
C                           INDEPENDENT AND DEPENDENT VARIABLES ACTUALLY
C                           USED IN FITTING THE MODEL.
C                BETA   - INPUT VECTOR OF LENGTH M CONTAINING REGRESSION
C                           COEFFICIENT ESTIMATES, IN THE ORDER
C                           SPECIFIED IN THE FITTED MODEL, IN THE FIRST
C                           M - 1 LOCATIONS, AND THE INTERCEPT ESTIMATE
C                           IN BETA(M).
C                SDR    - INPUT STANDARD DEVIATION OF RESIDUALS (SQUARE
C                           ROOT OF THE ERROR MEAN SQUARE--SEE OUTPUT
C                           PARAMETER ANOVA(12) IN RLMUL).
C                RES    - OUTPUT N BY 4 MATRIX CONTAINING THE RESIDUAL
C                           ANALYSIS CORRESPONDING TO XY.  COLUMNS 1,
C                           2, 3, AND 4 CONTAIN THE OBSERVED RESPONSES,
C                           PREDICTED RESPONSES, RESIDUALS, AND
C                           STANDARDIZED RESIDUALS, RESPECTIVELY.
C                IR     - INPUT ROW DIMENSION OF THE MATRIX RES EXACTLY
C                           AS SPECIFIED IN THE DIMENSION STATEMENT IN
C                           THE CALLING PROGRAM.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER=129 INDICATES THAT THE FIRST M ELEMENTS
C                             OF VECTOR IH DID NOT CONSTITUTE A SUBSET
C                             OF THE INTEGERS (1,2,..,MM).
C                           IER=130 INDICATES THAT M EXCEEDED MM OR THAT
C                             N WAS LESS THAN 1.
C
C   REQD. IMSL ROUTINES - SINGLE/UERTST,UGETIO
C                       - DOUBLE/UERTST,UGETIO,VXADD,VXMUL,VXSTO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS  1.  TYPICAL USAGE OF RLRES WOULD OCCUR SUBSEQUENT TO A
C                CALL TO IMSL ROUTINE RLMUL.
C            2.  FOLLOWING A CALL (OR CALLS) TO RLRES, THE USER MAY
C                WISH TO USE ROUTINES IN IMSL CHAPTER U TO PLOT THE
C                STANDARDIZED RESIDUALS OR TO PRINT THE MATRIX RES.
C            3.  TO CONSERVE CORE, THE ENTIRE DATA MATRIX NEED NOT BE
C                IN CORE AT ONCE. MULTIPLE CALLS MAY BE MADE TO RLRES
C                SUCH THAT XY IS A SUBSET OF THE FULL DATA MATRIX
C                USED IN FITTING THE MODEL.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE RLRES  (XY,IX,MM,N,IH,M,BETA,SDR,RES,IR,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            IX,MM,N,IH(1),M,IR,IER
      REAL               XY(IX,MM),BETA(M),SDR,RES(IR,4)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,ITEMP,J,JJ,K,L,M1
      DOUBLE PRECISION   STAT
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
      IF(M.LE.MM.AND.N.GE.1) GO TO 5
      IER = 130
      GO TO 9000
C                                  TERMINAL ERROR 2
    5 DO 10 I=1,MM
   10 IH(M+I) = I
      DO 25 I = 1,M
         DO 15 J=I,MM
            JJ = J
            IF(IH(M+J).EQ.IH(I)) GO TO 20
   15    CONTINUE
         IER = 129
         GO TO 9000
C                                  TERMINAL ERROR 1
   20    ITEMP = IH(M+I)
         IH(M+I) = IH(M+JJ)
         IH(M+JJ) = ITEMP
   25 CONTINUE
      L = IH(M)
      DO 35 I=1,N
         STAT = BETA(M)
         RES(I,1) = XY(I,L)
         M1 = M-1
         DO 30 J=1,M1
            K = IH(J)
            STAT = STAT + DBLE(BETA(J))*DBLE(XY(I,K))
   30    CONTINUE
         RES(I,2) = STAT
         RES(I,3) = RES(I,1) - STAT
         RES(I,4) = RES(I,3)/SDR
   35 CONTINUE
      GO TO 9005
 9000 CONTINUE
      CALL UERTST(IER,6HRLRES )
 9005 RETURN
      END
