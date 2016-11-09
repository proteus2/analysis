C   IMSL ROUTINE NAME   - RLSTP
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - REGRESSION MODEL SELECTION USING A FORWARD
C                           STEPWISE ALGORITHM WITH RESULTS AVAILABLE
C                           AFTER EACH STEP
C
C   USAGE               - CALL RLSTP (A,M,N,ALFAI,ALFAO,JX,IH,B,IOPT,
C                           IER)
C
C   ARGUMENTS    A      - ON INPUT, A IS AN (M+1) BY (M+1) SYMMETRIC
C                           MATRIX STORED IN SYMMETRIC STORAGE MODE.
C                           IT CONTAINS THE CORRECTED SUMS OF SQUARES
C                           AND CROSS PRODUCTS OF THE INDEPENDENT
C                           AND DEPENDENT VARIABLES. A VECTOR OF LENGTH
C                           (M+1)*(M+2)/2 IS REQUIRED.
C                         ON OUTPUT, A CONTAINS THE INVERSE OF THE
C                           SUBMATRIX OF THE CORRECTED SUMS OF SQUARES
C                           AND CROSS PRODUCTS MATRIX, THE
C                           REGRESSION COEFFICIENTS AND THE ERROR SUM
C                           OF SQUARES, FOR THE VARIABLES IN REGRESSION.
C                           THE REMAINING COMPONENTS OF A ARE GENERALLY
C                           NOT OF INTEREST. SEE REMARKS.
C                M      - INPUT NUMBER OF INDEPENDENT VARIABLES UNDER
C                           CONSIDERATION.
C                N      - INPUT NUMBER OF DATA POINTS.
C                ALFAI  - INPUT SIGNIFICANCE LEVEL FOR ENTERING
C                           VARIABLES. ALFAI MUST BE IN THE INCLUSIVE
C                           RANGE (0,.95). THE CHOICE 0.05 IS A COMMON
C                           ONE.
C                ALFAO  - INPUT SIGNIFICANCE LEVEL FOR DELETING
C                           VARIABLES. ALFAO MUST BE IN THE INCLUSIVE
C                           RANGE (0,.95), AND ALFAO MUST NOT BE LESS
C                           THAN ALFAI.
C                JX     - INPUT VECTOR OF LENGTH M CONTAINING ZEROS AND
C                           ONES USED TO CONTROL THE FORCING OF
C                           VARIABLES INTO THE MODEL.
C                         IF FORCING IS DESIRED, JX(I) SHOULD BE SET TO
C                           ONE IF THE I-TH INDEPENDENT VARIABLE WILL
C                           BE FORCED INTO THE MODEL,
C                           WHERE I = 1,2,...,M.
C                         IF FORCING IS NOT DESIRED, ALL ELEMENTS
C                           OF JX SHOULD BE SET TO ZERO.
C                IH     - INPUT/OUTPUT VECTOR OF LENGTH M+1.
C                         ON THE FIRST ENTRY TO RLSTP, IH(1)
C                           MUST BE EQUAL TO ZERO.
C                         ON OUTPUT, IH IS A VECTOR OF PLUS AND MINUS
C                           ONES. IH(I)=-1, IMPLIES THAT THE I-TH
C                           INDEPENDENT VARIABLE IS NOT IN THE MODEL.
C                           IH(I)=1 IMPLIES THAT THE I-TH INDEPENDENT
C                           VARIABLE IS IN THE MODEL. I=1,2,...,M.
C                B      - OUTPUT VECTOR OF LENGTH M+1 CONTAINING THE
C                           COEFFICIENTS OF THE FITTED MODEL. B(I) IS
C                           THE COEFFICIENT FOR THE I-TH INDEPENDENT
C                           VARIABLE WHEN IH(I) = 1. I = 1,2,...,M.
C                           OTHERWISE B(I) = 0.0.
C                IOPT   - ON INPUT, IOPT IS AN OPTION PARAMETER THAT
C                           DETERMINES WHETHER OR NOT RESULTS FROM
C                           RLSTP ARE RETURNED TO THE CALLING
C                           PROGRAM AFTER EACH STEP.
C                         IF IOPT IS EQUAL TO ONE, RESULTS ARE GIVEN
C                           AFTER EACH STEP.
C                         IF IOPT IS EQUAL TO ZERO, RESULTS ARE
C                           RETURNED ONLY AFTER THE FINAL MODEL HAS
C                           BEEN SELECTED.
C                         ON OUTPUT, IOPT IS UNCHANGED UNTIL THE FINAL
C                           MODEL HAS BEEN SELECTED, AT WHICH TIME
C                           IOPT IS EQUAL TO -1.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER=129 INDICATES THAT AN ERROR OCCURRED IN
C                             IMSL ROUTINE MDFD.
C                           IER=130 INDICATES THAT CYCLING (VARIABLES
C                             WERE CHECKED FOR ADDITION OR DELETION 2*M
C                             TIMES) OCCURRED.
C                           IER=131 INDICATES THAT THE SIGNIFICANCE
C                             LEVELS ALFAI AND ALFAO WERE SPECIFIED
C                             INCORRECTLY.
C                           IER=132 INDICATES THAT AN ATTEMPT TO ENTER
C                             A VARIABLE INTO THE MODEL CAUSED THE
C                             INVOLVED SUBMATRIX OF THE INPUT MATRIX A
C                             TO BE ALGORITHMICALLY SINGULAR. THE
C                             CORRESPONDING ELEMENT IN IH IS SET TO 2.
C                         WARNING ERROR
C                           IER=37 INDICATES THAT THE LAST ELEMENT OF
C                             INPUT A WAS ZERO (I.E. THE RESPONSE
C                             VARIABLE WAS CONSTANT) OR A PERFECT FIT
C                             TO THE DATA WAS OBTAINED, POSSIBLY
C                             OCCURRING PRIOR TO NORMAL (AS EXPECTED BY
C                             THE USER) TERMINATION OF RLSTP. OUTPUT
C                             INFORMATION FOR THE SELECTED MODEL IS
C                             CORRECT.
C
C   REQD. IMSL ROUTINES - MDFD,MERRC=ERFC,UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS  1.  AT EACH STEP, AFTER AN INDEPENDENT VARIABLE HAS
C                BEEN ADDED OR DELETED, A((M+1)*(M+2)/2) CONTAINS
C                THE ERROR SUM OF SQUARES CORRESPONDING TO THE SET
C                OF INDEPENDENT VARIABLES IN THE MODEL AT THAT TIME.
C            2.  AT EACH STEP, AFTER AN INDEPENDENT VARIABLE HAS
C                BEEN ADDED OR DELETED, THE INVERSE OF THE SUBMATRIX
C                OF THE CORRECTED SUMS OF SQUARES AND CROSS PRODUCTS
C                MATRIX, CORRESPONDING TO THE INDEPENDENT VARIABLES
C                PRESENTLY IN THE MODEL, IS AVAILABLE. THE ELEMENTS
C                OF THE INVERSE ARE SCATTERED THROUGHOUT A. THE
C                IMSL SUBROUTINE RLSUBM IS PROVIDED TO RETRIEVE AND
C                ORGANIZE THE INVERSE ELEMENTS IN A SUBMATRIX STORED
C                IN SYMMETRIC STORAGE MODE. THE ORDER OF THE SUB-
C                MATRIX WILL NEVER EXCEED THE NUMBER OF INDEPENDENT
C                VARIABLES IN THE MODEL. THE OUTPUT MATRIX FROM
C                RLSUBM, (CALLED S IN RLSUBM) MAY BE DIMENSIONED
C                ACCORDINGLY. WHEN ONLY FINAL RESULTS FROM RLSTP
C                ARE DESIRED, STORAGE REQUIREMENTS ARE MINIMIZED
C                IF THE MATRIX, A, IS REPLACED BY THE SUBMATRIX.
C            3.  WHEN RESULTS ARE TO BE RETURNED TO THE CALLING
C                PROGRAM AFTER EACH ADDITION OR DELETION OF AN
C                INDEPENDENT VARIABLE, IH(1) MUST BE SET EQUAL TO
C                ZERO FOR THE INITIAL CALL AND THE INITIAL CALL
C                ONLY. PARAMETERS IN THE CALLING SEQUENCE MUST NOT
C                BE MODIFIED BETWEEN THE INTERMEDIATE CALLS.
C            4.  WHEN RLSTP ATTEMPTS TO ENTER HIGHLY CORRELATED
C                VARIABLES INTO THE MODEL, AN OVERFLOW CAN OCCUR.
C                UNDER USUAL USAGE CONDITIONS THE ALGORITHM WILL
C                NOT ATTEMPT TO ENTER A VARIABLE THAT IS HIGHLY
C                CORRELATED WITH OTHER VARIABLES ALREADY IN THE
C                MODEL. HOWEVER, THIS SITUATION CAN OCCUR WHEN
C                HIGHLY CORRELATED VARIABLES ARE FORCED INTO THE
C                MODEL VIA THE USAGE OF THE JX VECTOR OR WHEN THE
C                SIGNIFICANCE LEVEL FOR ENTERING VARIABLES, ALFAI,
C                IS SET AT A HIGH VALUE. WHEN THE OVERFLOW OCCURS,
C                A RUN RETURNING RESULTS AFTER EVERY STEP IS
C                HELPFUL IN TAKING CORRECTIVE ACTION.
C            5.  TERMINAL ERROR IER=132 CAN OCCUR, FOR EXAMPLE WHEN
C                AN ATTEMPT IS MADE TO ENTER THE SECOND OF TWO
C                HIGHLY CORRELATED INDEPENDENT VARIABLES. TWO
C                SITUATIONS CAN MAKE THIS ERROR POSSIBLE.
C                A. ALFAI SET TO 1.0, WHICH EFFECTIVELY FORCES ALL
C                   VARIABLES INTO THE MODEL
C                B. FORCING TWO VERY HIGHLY CORRELATED VARIABLES
C                   INTO THE MODEL VIA THE USAGE OF THE JX VECTOR
C            6.  WHEN RLSTP IS USED WITH IOPT EQUAL TO ONE, IT MAY
C                BE NECESSARY FOR THE USER TO INCLUDE A LABELLED
C                COMMON BLOCK IN HIS/HER MAIN PROGRAM TO PROTECT
C                THE VALUES OF VARIABLES INTERNAL TO SUBROUTINE
C                RLSTP. THE COMMON BLOCK HAS THE FORM
C                   COMMON  /RLST1/DDD,III(12)
C                WHERE DDD IS TYPED DOUBLE PRECISION AND III IS
C                TYPED INTEGER IN THE CALLING PROGRAM.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE RLSTP  (A,M,N,ALFAI,ALFAO,JX,IH,B,IOPT,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            M,N,IH(1),JX(1),IOPT,IER
      REAL               A(1),ALFAI,ALFAO,B(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            JXSUM,I,IB,ICOUNT,IJ,ISW,ISW2,I1,J,K,KMAX,KX,
     1                   L,LD,LG,LL,LM1,LP1,LR,L1,MI,M1,M2,NDF,NERDF,
     2                   NV,NVV,N1,N2
      REAL               F,P
      REAL               ZERO,D,T,HUND,EPS,ZMAX,ZMAX1,ONE,PNF
      DOUBLE PRECISION   X,Z,R,Y,RR,DZERO
      COMMON             /RLST1/Z,K,LR,M1,ICOUNT,N2,M2,N1,MI,NV,LL,NVV,
     1                   NDF
      DATA               ZERO/0./,HUND/100./,ONE/1./,PNF/.95/
      DATA               EPS/Z3C100000/
      DATA               DZERO/0.0D0/
C                                  FIRST EXECUTABLE STATEMENT
      JXSUM = 0
      DO 1 I=1,M
         JXSUM = JXSUM+JX(I)
    1 CONTINUE
      IF (JXSUM .NE. 0) GO TO 2
C                                  NO FORCING
      ISW = 0
      NV = 1
      NVV = M+1
      GO TO 10
C                                  FORCING
    2 ISW = 1
      IF (IH(1).NE.0) GO TO 95
      LL = 1
      NV = 0
      DO 5 I=1,M
         IF (JX(I).EQ.1) NV = NV+1
    5 CONTINUE
      NVV = NV
   10 IF (IH(1).NE.0) GO TO 95
C                                  FIRST ENTRY
      IER = 0
      M1 = M+1
      K = ((M+2)*M1)/2
      LR = (M1*M)/2
      ICOUNT = 0
      N2 = N-2
      M2 = M+M
      N1 = N-1
      MI = MIN0(2,M)
      NDF = 0
      ZMAX = ZERO
      KMAX = 0
      DO 15 I=1,M1
         B(I) = ZERO
         KMAX = KMAX+I
         IH(I) = -1
         ZMAX1 = A(KMAX)
         IF (ZMAX1 .GT. ZMAX .AND. I .NE. M1) ZMAX = ZMAX1
   15 CONTINUE
      IF (ALFAI.GT.ALFAO) GO TO 16
      IF (ALFAI.GE.ZERO .AND. ALFAI.LE.PNF .AND. ALFAO.GE.ZERO .AND.
     2 ALFAO.LE.PNF) GO TO 20
C                                  TERMINAL ERROR - SIGNIFICANCE
C                                  LEVELS SPECIFIED INCORRECTLY
   16 CONTINUE
      IOPT = -1
      IER = 131
      GO TO 9000
C                                  FIND BEST VARIABLE TO ADD TO THE
C                                  MODEL
   20 IF (NDF.GE.M.OR.NDF.GE.N2) GO TO 115
      X = DZERO
      Z = DBLE(A(K))
      IF (Z .GT. DZERO) GO TO 23
C                                  PERFECT FIT
      IER = 37
      GO TO 115
   23 CONTINUE
      ICOUNT = ICOUNT+1
      L1 = 0
      IF (ISW.EQ.0.OR.NV.EQ.0) GO TO 30
C                                  THESE VARIABLES MUST BE FORCED INTO
C                                  THE MODEL
      LG = LL
      DO 25 I=LG,M
         L1 = L1+I
         IF (JX(I).NE.1) GO TO 25
         IF (IH(I).EQ.1) GO TO 25
         RR = DBLE(A(LR+I))
         R = DZERO
         IF (A(L1).GT.ZERO.AND.Z.GT.DZERO) R=RR*RR/(Z*A(L1))
         L = I
         LL = I+1
         X = R
         NV = NV-1
         KX = 1
         GO TO 40
   25 CONTINUE
   30 L1 = 0
      DO 35 I=1,M
         L1 = L1+I
         IF (IH(I).EQ.1) GO TO 35
         RR = DBLE(A(LR+I))
         R = DZERO
         IF (A(L1).GT.ZERO.AND.Z.GT.DZERO) R=RR*RR/(Z*A(L1))
         IF (R.LT.X) GO TO 35
         L = I
         X = R
   35 CONTINUE
C                                  FIND ERROR DEGREES OF FREEDOM AND
C                                  CALCULATE F PROBABILITY
      KX = 1
      NERDF = N2-NDF
      IF (X .GT. ONE-EPS) GO TO 40
      F = NERDF*X/(1.D0-X)
      CALL MDFD (F,1,NERDF,P,IER)
C                                  CONSOLIDATE ONE ERROR FROM MDFD
      IF (IER.GT.128) IER=129
      IF (IER.GE.128) GO TO 115
      P = 1.-P
      IF (P.GT.ALFAI) GO TO 115
C                                  ADD THE VARIABLE TO THE MODEL
C                                  PERFORM THE JORDAN REDUCTION
   40 NDF = NDF+KX
      LM1 = L-1
      LP1 = L+1
      L1 = (L*LM1)/2
      LD = L1+L
C                                  COEFFICIENT MATRIX IS ALGORITHMICALLY
C                                  SINGULAR WHEN VARIABLE L IS ADDED TO
C                                  THE MODEL
      IF (A(LD).EQ.ZERO) GO TO 110
      X = 1.D0/A(LD)
      B(L) = X
      Y = IH(L)*X
      IF (L.EQ.1) GO TO 50
      DO 45 I=1,LM1
         B(I) = A(L1+I)*X
   45 CONTINUE
   50 J = L+LD
      DO 55 I=LP1,M1
         B(I) = A(J)*Y*IH(I)
         J = J+I
   55 CONTINUE
      I1 = 0
      DO 70 I=1,M1
         IF (I.EQ.L) GO TO 65
         IB = (I*I1)/2
         R = A(IB+L)
         IF (I.LT.L) R = A(L1+I)*IH(L)*IH(I)
         DO 60 J=1,I
            IF (J.EQ.L) GO TO 60
            IJ = IB+J
            D = A(IJ)-B(J)*R
            T = D/HUND+A(IJ)
            IF (T.EQ.A(IJ)) D = ZERO
            A(IJ) = D
   60    CONTINUE
         IF (A(IJ).LT.ZERO) A(IJ) = ZERO
   65    I1 = I
   70 CONTINUE
      IF (L.EQ.1) GO TO 80
      DO 75 I=1,LM1
         A(L1+I) = B(I)
         B(I) = ZERO
   75 CONTINUE
   80 J = L+L1
      A(J) = X
      B(L) = ZERO
      J = J+L
      DO 85 I=LP1,M1
         A(J) = -A(J)*X
         B(I) = ZERO
         J = I+J
   85 CONTINUE
      IH(L) = -IH(L)
      IF (ICOUNT.LE.M2) GO TO 90
      IER = 130
      GO TO 120
   90 IF (IOPT.NE.0) GO TO 120
   95 IF (ISW.EQ.1.AND.NV.NE.0) GO TO 20
      IF (ISW.EQ.1.AND.NVV.EQ.M) GO TO 115
      IF (ICOUNT.LT.MI.OR.MOD(ICOUNT,2).EQ.0) GO TO 20
      IF (NDF.LE.M.AND.M.LE.2) GO TO 20
C                                  DETERMINE WHETHER A VARIABLE SHOULD
C                                  BE DELETED
      Z = DBLE(A(K))
      IF (Z .GT. DZERO) GO TO 96
C                                  PERFECT FIT. DO NOT ATTEMPT TO
C                                    DELETE ANY VARIABLES FROM THE
C                                    MODEL.
      IER = 37
      GO TO 115
   96 ISW2 = 0
      ICOUNT = ICOUNT+1
      L1 = 0
      KX = -1
      J = 0
      L = 0
      DO 105 I=1,M
         L1 = L1+I
         IF (ISW.EQ.1) J = JX(I)
         IF (IH(I).EQ.-1.OR.J.EQ.1) GO TO 105
         R = DBLE(A(LR+I))
         R = R*R/A(L1)
         IF (ISW2.EQ.0) GO TO 100
         IF (R.GT.X) GO TO 105
  100    ISW2 = 1
         X = R
         L = I
  105 CONTINUE
      IF (L .EQ. 0) GO TO 20
      NERDF = N1-NDF
      F = X*NERDF/Z
      CALL MDFD (F,1,NERDF,P,IER)
C                                  CONSOLIDATE ERROR FROM MDFD
      IF (IER.GT.128) IER=129
      IF (IER.GE.128) GO TO 115
      P = 1.-P
      IF (P.GT.ALFAO) GO TO 40
      GO TO 20
  110 IER = 132
      IH(L) = 2
  115 IOPT = -1
      IF (A(K) .LE. ZERO .AND. IER .EQ. 0) IER = 37
  120 DO 125 I=1,M1
         IF (IH(I).EQ.1) B(I) = -A(LR+I)
  125 CONTINUE
      IF (IER.EQ.0) GO TO 9005
 9000 CONTINUE
      CALL UERTST (IER,6HRLSTP )
 9005 RETURN
      END
