C   IMSL ROUTINE NAME   - USPDF
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1982
C
C   PURPOSE             - PLOT OF TWO SAMPLE CUMULATIVE PROBABILITY
C                           DISTRIBUTION FUNCTIONS AGAINST THEIR SPECTRA
C
C   USAGE               - CALL USPDF (X,N,M,W,IW,IR)
C
C   ARGUMENTS    X      - VECTOR CONTAINING SAMPLE ONE FOLLOWED
C                           BY SAMPLE TWO. (INPUT/OUTPUT)
C                           X HAS LENGTH N+M.
C                         ON OUTPUT, X IS DESTROYED.
C                N      - SIZE OF SAMPLE ONE. (INPUT)
C                M      - SIZE OF SAMPLE TWO. (INPUT)
C                W      - (N+M) BY 2 MATRIX, USED AS WORK STORAGE.
C                IW     - ROW DIMENSION OF MATRIX W EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM. (INPUT)
C                IR     - ARRAY OF SIZE N+M USED AS WORK STORAGE.
C
C   REQD. IMSL ROUTINES - UERTST,UGETIO,USPLO,VSRTR
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS  1.  OUTPUT IS WRITTEN TO THE UNIT SPECIFIED BY IMSL
C                ROUTINE UGETIO. SEE THE UGETIO DOCUMENT FOR DETAILS.
C            2.  THE PLOTTED FUNCTIONS ARE STEP FUNCTIONS. THE POINTS
C                PLOTTED SHOULD BE CONNECTED BY THE USER STARTING
C                AT THE LOWER LEFT CORNER OF THE GRAPH.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE USPDF  (X,N,M,W,IW,IR)
C
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N,M,IW,IR(1)
      REAL               X(1),W(IW,1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,IER,J,K,K1,NM,N1
      REAL               XO(2),XI(2),RANGE(4)
      DATA               RANGE/3*0.0,1.0/
C                                  FIRST EXECUTABLE STATEMENT
      CALL UGETIO(1,NIN,NOUT)
C                                  ASSIGN SAMPLE NUMBER TO ELEMENTS IN X
      DO 5 I = 1,N
    5 IR(I) = 1
      NM = N+M
      N1 = N+1
      DO 10 I = N1,NM
   10 IR(I) = 2
C                                  SORT THE SAMPLES
      CALL VSRTR (X,NM,IR)
C                                  CALCULATE THE CUMULATIVE PDF
      XI(1) = 1.0/N
      XI(2) = 1.0/M
      XO(1) = 0.
      XO(2) = 0.
      DO 15 I = 1,NM
         J = IR(I)
         XO(J) = XO(J)+XI(J)
         W(I,1) = XO(1)
   15 W(I,2) = XO(2)
C                                  ELIMINATE MULTIPLE POINTS
      K = 1
      J = 1
      DO 20 I = 2,NM
         IF (X(J) .EQ. X(I)) GO TO 20
         K1 = I-1
         X(K) = X(K1)
         W(K,1) = W(K1,1)
         W(K,2) = W(K1,2)
         K = K+1
         J = I
   20 CONTINUE
      X(K) = X(NM)
      W(K,1) = 1.
      W(K,2) = 1.
      NM = NM+1
      CALL USPLO (X,W,IW,K,2,1,
     * 52HCUMULATIVE SAMPLE PROBABILITY DISTRIBUTION FUNCTIONS,52,
     * 13HSAMPLE VALUES,13,11HPROBABILITY,11,RANGE,2H12,1,IER)
      WRITE (NOUT,25)
   25 FORMAT(/25X,30H SAMPLE 1 = 1     SAMPLE 2 = 2)
      RETURN
      END
