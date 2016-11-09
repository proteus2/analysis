C   IMSL ROUTINE NAME   - USPC
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - NOVEMBER 1, 1984
C
C   PURPOSE             - PRINT A SAMPLE CDF, A THEORETICAL CDF AND
C                           CONFIDENCE BAND INFORMATION. PLOT THESE
C                           ON OPTION.
C
C   USAGE               - CALL USPC (CDF,X,N,N12,N95,IP,IC,W)
C
C   ARGUMENTS    CDF    - USER SUPPLIED PROBABILITY DISTRIBUTION
C                           FUNCTION CONSISTING OF 2 ARGUMENTS (X,P).
C                           X IS THE SAMPLE POINT AND P IS THE
C                           RESULTING THEORETICAL PROBABILITY AT THE
C                           POINT X (INTEGERAL OF THE DENSITY TO X).
C                           CDF MUST APPEAR IN AN EXTERNAL STATEMENT
C                           IN THE CALLING PROGRAM. (INPUT)
C                X      - INPUT VECTOR OF LENGTH N CONTAINING SAMPLE OF
C                           SIZE N. X MUST BE SORTED INTO ASCENDING
C                           ORDER PRIOR TO ENTERING USPC.
C                         ON OUTPUT, X IS DESTROYED. (INPUT/OUTPUT)
C                N      - SIZE OF SAMPLE. (INPUT)
C                N12    - CONFIDENCE BAND OPTION. (INPUT)
C                           IF N12 = 2,TWO-SIDED CONFIDENCE BAND
C                             INFORMATION IS DESIRED.
C                           IF N12 = 1, POSITIVE ONE-SIDED CONFIDENCE
C                             BAND INFORMATION IS DESIRED.
C                           OTHERWISE, NEGATIVE ONE-SIDED CONFIDENCE
C                             BAND INFORMATION IS DESIRED.
C                N95    - CONFIDENCE BAND OPTION. (INPUT)
C                           IF N95 = 95,THE 95 PER CENT BAND IS DESIRED.
C                           OTHERWISE, THE 99 PER CENT BAND IS DESIRED.
C                IP     - PLOTTING OPTION. (INPUT)
C                           IF IP = 1, A PLOT OF THE SAMPLE VS THE
C                             THEORETICAL CDF IS DESIRED.
C                           OTHERWISE, ONLY PRINTING OCCURS.
C                IC     - CONFIDENCE BAND OPTION. (INPUT)
C                           IF IC = 1, CONFIDENCE BAND INFORMATION
C                             WILL ALSO BE PLOTTED.
C                W      - WORK AREA OF SIZE 4*N.  SEE REMARKS.
C
C   REQD. IMSL ROUTINES - UERTST,UGETIO,USPLO,USPKD
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS  1.  IF A GRAPH IS REQUESTED, THE USER MUST COMPLETE THE
C                GRAPH BY CONNECTING THE POINTS AND REMEMBERING THAT
C                SAMPLE CDFS ARE STEP FUNCTIONS.
C            2.  W MAY BE DIMENSIONED 3N INSTEAD OF 4N FOR A ONE-SIDED
C                CONFIDENCE BAND.
C            3.  CONFIDENCE BANDS ARE PLOTTED AROUND THE SAMPLE CDF.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE USPC   (CDF,X,N,N12,N95,IP,IC,W)
C
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N,N12,N95,IP,IC
      REAL               X(1),W(N,1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,IER,IJ,IX1,IX2,J,K,KEND
      REAL               RANGE(4),D,SN,T,XINC,XO
      DATA               RANGE/4*0.0/
C                                  FIRST EXECUTABLE STATEMENT
      CALL UGETIO(1,NIN,NOUT)
      IJ = N95
      IF (IJ .NE. 95) IJ = 99
      SN=N
      XO=0.
      XINC = 1./SN
      SN = SQRT(XINC)
      IF (N12 .EQ. 2) GO TO 5
C                                  ONE-SIDED BAND
C                                  99 PER CENT BAND
C
      D = 1.517428*SN
C                                  95 PER CENT BAND
C
      IF (N95 .EQ. 95) D = 1.223849*SN
C                                  NEGATIVE BAND
C
      IF (N12 .NE. 1) D = -D
      KEND = 3
      GO TO 10
C                                  TWO-SIDED BAND - 99 PER CENT
    5 D = -1.627628*SN
C                                  95 PERCENT
      IF (N95 .EQ. 95) D = -1.358100*SN
      KEND = 4
   10 J = 1
      K = 1
      DO 20 I = 2,N
         XO = XO+XINC
         IF (X(J) .EQ. X(I)) GO TO 20
C                                  ELIMINATE MULTIPLE POINTS
         X(K) = X(J)
C                                  CALCULATE SAMPLE CDF
         W(K,1) = XO
C                                  CALCULATE THEORETICAL CDF
         CALL CDF (X(K),W(K,2))
         T = W(K,1)+D
         IF (T .GT. 1.) T = 1.
         IF (T .LT. 0.) T = 0.
         W(K,3) = T
         IF (KEND .EQ. 3) GO TO 15
         T = W(K,1)-D
         IF (T .GT. 1.) T = 1.
         IF (T .LT. 0.) T = 0.
         W(K,4) = T
   15    K = K+1
   20 J = I
      X(K) = X(N)
      W(K,1) = 1.0
      CALL CDF(X(K),W(K,2))
      T = W(K,1)+D
      IF (T .GT. 1.) T = 1.
      IF (T .LT. 0.) T = 0.
      W(K,3) = T
      IF (KEND .EQ. 3) GO TO 25
      T = W(K,1)-D
      IF (T .GT. 1.) T = 1.
      IF (T .LT. 0.) T = 0.
      W(K,4) = T
C                                  PRINT TABLE
   25 KEND = 1
      J = 50
   30 J = AMIN0 (J,K)
      WRITE (NOUT,50) IJ
C                                  ONE-SIDED BAND
C                                  POSITIVE BAND
      IX1 = 1
      IX2 = 3
      IF (N12 .EQ. 1) GO TO 35
C                                  NEGATIVE BAND
      IX1 = 3
      IX2 = 1
C                                  TWO-SIDED BAND
      IF (N12 .EQ. 2) IX2 = 4
   35 DO 40 I = KEND,J
         WRITE (NOUT,55) X(I),W(I,1),W(I,2),W(I,IX1),W(I,IX2)
   40 CONTINUE
      IF (J .GE. K) GO TO 45
      KEND = KEND+50
      J = J+50
      GO TO 30
   45 IF (IP .NE. 1) GO TO 80
      KEND = 2
      IF (IC .EQ. 1) KEND = 3
      IF (N12 .EQ. 2 .AND. IC .EQ. 1) KEND = 4
      J = KEND+1
   50 FORMAT(/1H1,13X,1HX,18X,5HFN(X),15X,4HF(X),17X,I2,14H PER CENT BAN
     *D)
   55 FORMAT(3E20.6,7X,1H(,E13.6,1H,,E13.6,1H))
      CALL USPLO (X,W,N,K,KEND,1,
     * 38HCUMULATIVE SAMPLE AND THEORETICAL CDFS,38,13HSAMPLE VALUES,
     * 13,11HPROBABILITY,11,RANGE,10H1234567890,1,IER)
      WRITE (NOUT,60)
   60 FORMAT(/25X,15H SAMPLE CDF = 1,5X,20H THEORETICAL CDF = 2)
      IF (IC .NE. 1) GO TO 80
      IF (N12 .NE. 2) GO TO 70
      WRITE (NOUT,65)
   65 FORMAT(/25X,27H CONFIDENCE BANDS = 3 AND 4)
      GO TO 80
   70 WRITE (NOUT,75)
   75 FORMAT(/25X,20H CONFIDENCE BAND = 3)
   80 RETURN
      END
