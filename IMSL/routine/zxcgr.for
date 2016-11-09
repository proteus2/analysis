C   IMSL ROUTINE NAME   - ZXCGR
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - NOVEMBER 1, 1979
C
C   PURPOSE             - A CONJUGATE GRADIENT ALGORITHM FOR FINDING
C                           THE MINIMUM OF A FUNCTION OF N VARIABLES
C
C   USAGE               - CALL ZXCGR (FUNCT,N,ACC,MAXFN,DFPRED,X,G,F,W,
C                           IER)
C
C   ARGUMENTS    FUNCT  - A USER SUPPLIED SUBROUTINE WHICH CALCULATES
C                           THE OBJECTIVE FUNCTION AND ITS GRADIENT
C                           FOR GIVEN PARAMETER VALUES
C                           X(1),X(2),...,X(N).
C                           THE CALLING SEQUENCE HAS THE FOLLOWING FORM
C                           CALL FUNCT (N,X,F,G)
C                           WHERE X AND G ARE VECTORS OF LENGTH N.
C                           THE SCALAR F IS FOR THE OBJECTIVE FUNCTION.
C                           G(1), G(2), ..., G(N) ARE FOR THE COMPONENTS
C                           OF THE GRADIENT OF F.
C                           FUNCT MUST APPEAR IN AN EXTERNAL STATEMENT
C                           IN THE CALLING PROGRAM. FUNCT MUST NOT
C                           ALTER THE VALUES OF X(I),I=1,...,N OR N.
C                N      - THE NUMBER OF PARAMETERS OF THE OBJECTIVE
C                           FUNCTION. (INPUT) (I.E.,THE LENGTH OF X)
C                ACC    - CONVERGENCE CRITERION. (INPUT)
C                           THE CALCULATION ENDS WHEN THE SUM OF SQUARES
C                           OF THE COMPONENTS OF G IS LESS THAN ACC.
C                MAXFN  - MAXIMUM NUMBER OF FUNCTION EVALUATIONS (I.E.,
C                           CALLS TO SUBROUTINE FUNCT) ALLOWED. (INPUT)
C                           IF MAXFN IS SET TO ZERO, THEN THERE IS
C                           NO RESTRICTION ON THE NUMBER OF FUNCTION
C                           EVALUATIONS.
C                DFPRED - A ROUGH ESTIMATE OF THE EXPECTED REDUCTION
C                           IN F, WHICH IS USED TO DETERMINE THE SIZE
C                           OF THE INITIAL CHANGE TO X. (INPUT)
C                           NOTE THAT DFPRED IS THE EXPECTED REDUCTION
C                           ITSELF, AND DOES NOT DEPEND ON ANY RATIOS.
C                           A BAD VALUE OF DFPRED CAUSES AN ERROR
C                           MESSAGE, WITH IER=129, AND A RETURN ON THE
C                           FIRST ITERATION. (SEE THE DESCRIPTION OF
C                           IER BELOW)
C                X      - VECTOR OF LENGTH N CONTAINING PARAMETER
C                           VALUES.
C                         ON INPUT, X MUST CONTAIN THE INITIAL
C                           PARAMETER ESTIMATES.
C                         ON OUTPUT, X CONTAINS THE FINAL PARAMETER
C                           ESTIMATES AS DETERMINED BY ZXCGR.
C                G      - A VECTOR OF LENGTH N CONTAINING THE
C                           COMPONENTS OF THE GRADIENT OF F AT THE
C                           FINAL PARAMETER ESTIMATES. (OUTPUT)
C                F      - A SCALAR CONTAINING THE VALUE OF THE FUNCTION
C                           AT THE FINAL PARAMETER ESTIMATES. (OUTPUT)
C                W      - WORK VECTOR OF LENGTH 6*N.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         IER = 0 IMPLIES THAT CONVERGENCE WAS
C                           ACHIEVED AND NO ERRORS OCCURRED.
C                         TERMINAL ERROR
C                           IER = 129 IMPLIES THAT THE LINE SEARCH OF
C                             AN INTEGRATION WAS ABANDONED. THIS
C                             ERROR MAY BE CAUSED BY AN ERROR IN THE
C                             GRADIENT.
C                           IER = 130 IMPLIES THAT THE CALCULATION
C                             CANNOT CONTINUE BECAUSE THE SEARCH
C                             DIRECTION IS UPHILL.
C                           IER = 131 IMPLIES THAT THE ITERATION WAS
C                             TERMINATED BECAUSE MAXFN WAS EXCEEDED.
C                           IER = 132 IMPLIES THAT THE CALCULATION
C                             WAS TERMINATED BECAUSE TWO CONSECUTIVE
C                             ITERATIONS FAILED TO REDUCE F.
C
C   REQD. IMSL ROUTINES - UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS  1.  THE ROUTINE INCLUDES NO THOROUGH CHECKS ON THE PART
C                OF THE USER PROGRAM THAT CALCULATES THE DERIVATIVES
C                OF THE OBJECTIVE FUNCTION. THEREFORE, BECAUSE
C                DERIVATIVE CALCULATION IS A FREQUENT SOURCE OF
C                ERROR, THE USER SHOULD VERIFY INDEPENDENTLY THE
C                CORRECTNESS OF THE DERIVATIVES THAT ARE GIVEN TO
C                THE ROUTINE.
C            2.  BECAUSE OF THE CLOSE RELATION BETWEEN THE CONJUGATE
C                GRADIENT METHOD AND THE METHOD OF STEEPEST DESCENTS,
C                IT IS VERY HELPFUL TO CHOOSE THE SCALE OF THE
C                VARIABLES IN A WAY THAT BALANCES THE MAGNITUDES OF
C                THE COMPONENTS OF A TYPICAL DERIVATE VECTOR. IT
C                CAN BE PARTICULARLY INEFFICIENT IF A FEW COMPONENTS
C                OF THE GRADIENT ARE MUCH LARGER THAN THE REST.
C            3.  IF THE VALUE OF THE PARAMETER ACC IN THE ARGUMENT
C                LIST OF THE ROUTINE IS SET TO ZERO, THEN THE
C                SUBROUTINE WILL CONTINUE ITS CALCULATION UNTIL IT
C                STOPS REDUCING THE OBJECTIVE FUNCTION. IN THIS CASE
C                THE USUAL BEHAVIOUR IS THAT CHANGES IN THE
C                OBJECTIVE FUNCTION BECOME DOMINATED BY COMPUTER
C                ROUNDING ERRORS BEFORE PRECISION IS LOST IN THE
C                GRADIENT VECTOR. THEREFORE, BECAUSE THE POINT OF
C                VIEW HAS BEEN TAKEN THAT THE USER REQUIRES THE
C                LEAST POSSIBLE VALUE OF THE FUNCTION, A VALUE OF
C                THE OBJECTIVE FUNCTION THAT IS SMALL DUE TO
C                COMPUTER ROUNDING ERRORS CAN PREVENT FURTHER
C                PROGRESS. HENCE THE PRECISION IN THE FINAL VALUES
C                OF THE VARIABLES MAY BE ONLY ABOUT HALF THE
C                NUMBER OF SIGNIFICANT DIGITS IN THE COMPUTER
C                ARITHMETIC, BUT THE LEAST VALUE OF F IS USUALLY
C                FOUND TO QUITE HIGH ACCURACY.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE ZXCGR  (FUNCT,N,ACC,MAXFN,DFPRED,X,G,F,W,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N,MAXFN,IER
      REAL               ACC,DFPRED,X(N),G(N),F,W(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            MAXLIN,MXFCON,I,IGINIT,IGOPT,IRETRY,IRSDG,
     1                   IRSDX,ITERC,ITERFM,ITERRS,IXOPT,NCALLS,NFBEG,
     2                   NFOPT
      REAL               BETA,DDSPLN,DFPR,FCH,FINIT,FMIN,GAMDEN,GAMA,
     1                   GINIT,GMIN,GNEW,GSPLN,GSQRD,SBOUND,STEP,STEPCH,
     2                   STMIN,SUM,WORK
      DATA               MAXLIN/5/,MXFCON/2/
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
C                                  THE WORKING SPACE ARRAY IS SPLIT
C                                    INTO SIX VECTORS OF LENGTH N. THE
C                                    FIRST PART IS USED FOR THE SEARCH
C                                    DIRECTION OF AN ITERATION. THE
C                                    SECOND AND THIRD PARTS CONTAIN THE
C                                    INFORMATION THAT IS REQUIRED BY
C                                    THE CONJUGACY CONDITIONS OF THE
C                                    RESTART PROCEDURE. THE FOURTH PART
C                                    CONTAINS THE GRADIENT AT THE START
C                                    OF AN ITERATION. THE FIFTH PART
C                                    CONTAINS THE PARAMETERS THAT GIVE
C                                    THE LEAST CALCULATED VALUE OF F.
C                                    THE SIXTH PART CONTAINS THE
C                                    GRADIENT VECTOR WHERE F IS LEAST.
      IRSDX = N
      IRSDG = IRSDX+N
      IGINIT = IRSDG+N
      IXOPT = IGINIT+N
      IGOPT = IXOPT+N
C                                  SET SOME PARAMETERS TO BEGIN THE
C                                    CALCULATION. ITERC AND
C                                    NCALLS COUNT THE NUMBER OF
C                                    ITERATIONS AND CALLS OF FUNCT.
C                                    ITERFM IS THE NUMBER OF THE MOST
C                                    RECENT ITERATION THAT DECREASES F.
      ITERC = 0
      NCALLS = 0
      ITERFM = ITERC
C                                  CALL SUBROUTINE FUNCT. LET THE
C                                    INITIAL SEARCH DIRECTION BE MINUS
C                                    THE GRADIENT VECTOR. USUALLY THE
C                                    PARAMETER ITERRS GIVES THE
C                                    ITERATION NUMBER OF THE MOST
C                                    RECENT RESTART, BUT IT IS SET TO
C                                    ZERO WHEN THE STEEPEST DESCENT
C                                    DIRECTION IS USED.
    5 NCALLS = NCALLS+1
      CALL FUNCT (N,X,F,G)
      IF (NCALLS.GE.2) GO TO 20
   10 DO 15 I=1,N
   15 W(I) = -G(I)
      ITERRS = 0
      IF (ITERC.GT.0) GO TO 80
C                                  SET SUM TO G SQUARED. GMIN AND GNEW
C                                    ARE THE OLD AND THE NEW
C                                    DIRECTIONAL DERIVATIVES ALONG THE
C                                    CURRENT SEARCH DIRECTION. LET FCH
C                                    BE THE DIFFERENCE BETWEEN F AND
C                                    THE PREVIOUS BEST VALUE OF THE
C                                    OBJECTIVE FUNCTION.
   20 GNEW = 0.0
      SUM = 0.0
      DO 25 I=1,N
         GNEW = GNEW+W(I)*G(I)
   25 SUM = SUM+G(I)**2
      IF (NCALLS.EQ.1) GO TO 35
      FCH = F-FMIN
C                                  STORE THE VALUES OF X, F AND G, IF
C                                    THEY ARE THE BEST THAT HAVE BEEN
C                                    CALCULATED SO FAR, AND NOTE G
C                                    SQUARED AND THE VALUE OF NCALLS.
C                                    TEST FOR CONVERGENCE.
      IF (FCH) 35,30,50
   30 IF (GNEW/GMIN.LT.-1.0) GO TO 45
   35 FMIN = F
      GSQRD = SUM
      NFOPT = NCALLS
      DO 40 I=1,N
         W(IXOPT+I) = X(I)
   40 W(IGOPT+I) = G(I)
   45 IF (SUM.LE.ACC) GO TO 9005
C                                  TEST IF THE VALUE OF MAXFN ALLOWS
C                                    ANOTHER CALL OF FUNCT.
   50 IF (NCALLS.NE.MAXFN) GO TO 55
      IER = 131
      GO TO 9000
   55 IF (NCALLS.GT.1) GO TO 100
C                                  SET DFPR TO THE ESTIMATE OF THE
C                                    REDUCTION IN F GIVEN IN THE
C                                    ARGUMENT LIST, IN ORDER THAT THE
C                                    INITIAL CHANGE TO THE PARAMETERS
C                                    IS OF A SUITABLE SIZE. THE VALUE
C                                    OF STMIN IS USUALLY THE
C                                    STEP-LENGTH OF THE MOST RECENT
C                                    LINE SEARCH THAT GIVES THE LEAST
C                                    CALCULATED VALUE OF F.
      DFPR = DFPRED
      STMIN = DFPRED/GSQRD
C                                  BEGIN THE ITERATION
   80 ITERC = ITERC+1
C                                  STORE THE INITIAL FUNCTION VALUE AND
C                                    GRADIENT, CALCULATE THE INITIAL
C                                    DIRECTIONAL DERIVATIVE, AND BRANCH
C                                    IF ITS VALUE IS NOT NEGATIVE. SET
C                                    SBOUND TO MINUS ONE TO INDICATE
C                                    THAT A BOUND ON THE STEP IS NOT
C                                    KNOWN YET, AND SET NFBEG TO THE
C                                    CURRENT VALUE OF NCALLS. THE
C                                    PARAMETER IRETRY SHOWS THE NUMBER
C                                    OF ATTEMPTS AT SATISFYING THE BETA
C                                    CONDITION.
      FINIT = F
      GINIT = 0.0
      DO 85 I=1,N
         W(IGINIT+I) = G(I)
   85 GINIT = GINIT+W(I)*G(I)
      IF (GINIT.GE.0.0) GO TO 165
      GMIN = GINIT
      SBOUND = -1.0
      NFBEG = NCALLS
      IRETRY = -1
C                                  SET STEPCH SO THAT THE INITIAL
C                                    STEP-LENGTH IS CONSISTENT WITH THE
C                                    PREDICTED REDUCTION IN F, SUBJECT
C                                    TO THE CONDITION THAT IT DOES NOT
C                                    EXCEED THE STEP-LENGTH OF THE
C                                    PREVIOUS ITERATION. LET STMIN BE
C                                    THE STEP TO THE LEAST CALCULATED
C                                    VALUE OF F.
      STEPCH = AMIN1(STMIN,ABS(DFPR/GINIT))
      STMIN = 0.0
C                                  CALL SUBROUTINE FUNCT AT THE VALUE
C                                    OF X THAT IS DEFINED BY THE NEW
C                                    CHANGE TO THE STEP-LENGTH, AND LET
C                                    THE NEW STEP-LENGTH BE STEP. THE
C                                    VARIABLE WORK IS USED AS WORK
C                                    SPACE.
   90 STEP = STMIN+STEPCH
      WORK = 0.0
      DO 95 I=1,N
         X(I) = W(IXOPT+I)+STEPCH*W(I)
   95 WORK = AMAX1(WORK,ABS(X(I)-W(IXOPT+I)))
      IF (WORK.GT.0.0) GO TO 5
C                                  TERMINATE THE LINE SEARCH IF STEPCH
C                                    IS EFFECTIVELY ZERO.
      IF (NCALLS.GT.NFBEG+1) GO TO 115
      IF (ABS(GMIN/GINIT)-0.2) 170,170,115
C                                  LET SPLN BE THE QUADRATIC SPLINE
C                                    THAT INTERPOLATES THE CALCULATED
C                                    FUNCTION VALUES AND DIRECTIONAL
C                                    DERIVATIVES AT THE POINTS STMIN
C                                    AND STEP OF THE LINE SEARCH, WHERE
C                                    THE KNOT OF THE SPLINE IS AT
C                                    0.5*(STMIN+STEP). REVISE STMIN,
C                                    GMIN AND SBOUND, AND SET DDSPLN TO
C                                    THE SECOND DERIVATIVE OF SPLN AT
C                                    THE NEW STMIN. HOWEVER, IF FCH IS
C                                    ZERO, IT IS ASSUMED THAT THE
C                                    MAXIMUM ACCURACY IS ALMOST
C                                    ACHIEVED, SO DDSPLN IS CALCULATED
C                                    USING ONLY THE CHANGE IN THE
C                                    GRADIENT.
  100 WORK = (FCH+FCH)/STEPCH-GNEW-GMIN
      DDSPLN = (GNEW-GMIN)/STEPCH
      IF (NCALLS.GT.NFOPT) SBOUND = STEP
      IF (NCALLS.GT.NFOPT) GO TO 105
      IF (GMIN*GNEW.LE.0.0) SBOUND = STMIN
      STMIN = STEP
      GMIN = GNEW
      STEPCH = -STEPCH
  105 IF (FCH.NE.0.0) DDSPLN = DDSPLN+(WORK+WORK)/STEPCH
C
C                                  TEST FOR CONVERGENCE OF THE LINE
C                                    SEARCH, BUT FORCE AT LEAST TWO
C                                    STEPS TO BE TAKEN IN ORDER NOT TO
C                                    LOSE QUADRATIC TERMINATION.
      IF (GMIN.EQ.0.0) GO TO 170
      IF (NCALLS.LE.NFBEG+1) GO TO 120
      IF (ABS(GMIN/GINIT).LE.0.2) GO TO 170
C                                  APPLY THE TEST THAT DEPENDS ON THE
C                                    PARAMETER MAXLIN.
  110 IF (NCALLS.LT.NFOPT+MAXLIN) GO TO 120
  115 IER = 129
      GO TO 170
C                                  SET STEPCH TO THE GREATEST CHANGE TO
C                                    THE CURRENT VALUE OF STMIN THAT IS
C                                    ALLOWED BY THE BOUND ON THE LINE
C                                    SEARCH. SET GSPLN TO THE GRADIENT
C                                    OF THE QUADRATIC SPLINE AT
C                                    (STMIN+STEPCH). HENCE CALCULATE
C                                    THE VALUE OF STEPCH THAT MINIMIZES
C                                    THE SPLINE FUNCTION, AND THEN
C                                    OBTAIN THE NEW FUNCTION AND
C                                    GRADIENT VECTOR, FOR THIS VALUE OF
C                                    THE CHANGE TO THE STEP-LENGTH.
  120 STEPCH = 0.5*(SBOUND-STMIN)
      IF (SBOUND.LT.-0.5) STEPCH = 9.0*STMIN
      GSPLN = GMIN+STEPCH*DDSPLN
      IF (GMIN*GSPLN.LT.0.0) STEPCH = STEPCH*GMIN/(GMIN-GSPLN)
      GO TO 90
C                                  CALCULATE THE VALUE OF BETA THAT
C                                    OCCURS IN THE NEW SEARCH
C                                    DIRECTION.
  125 SUM = 0.0
      DO 130 I=1,N
  130 SUM = SUM+G(I)*W(IGINIT+I)
      BETA = (GSQRD-SUM)/(GMIN-GINIT)
C                                  TEST THAT THE NEW SEARCH DIRECTION
C                                    CAN BE MADE DOWNHILL. IF IT
C                                    CANNOT, THEN MAKE ONE ATTEMPT TO
C                                    IMPROVE THE ACCURACY OF THE LINE
C                                    SEARCH.
      IF (ABS(BETA*GMIN).LE.0.2*GSQRD) GO TO 135
      IRETRY = IRETRY+1
      IF (IRETRY.LE.0) GO TO 110
C                                  APPLY THE TEST THAT DEPENDS ON THE
C                                    PARAMETER MXFCON.
C                                    SET DFPR TO THE PREDICTED
C                                    REDUCTION IN F ON THE NEXT
C                                    ITERATION.
  135 IF (F.LT.FINIT) ITERFM = ITERC
      IF (ITERC.LT.ITERFM+MXFCON) GO TO 140
      IER = 132
      GO TO 9000
  140 DFPR = STMIN*GINIT
C                                  BRANCH IF A RESTART PROCEDURE IS
C                                    REQUIRED DUE TO THE ITERATION
C                                    NUMBER OR DUE TO THE SCALAR
C                                    PRODUCT OF CONSECUTIVE GRADIENTS.
      IF (IRETRY.GT.0) GO TO 10
      IF (ITERRS.EQ.0) GO TO 155
      IF (ITERC-ITERRS.GE.N) GO TO 155
      IF (ABS(SUM).GE.0.2*GSQRD) GO TO 155
C                                  CALCULATE THE VALUE OF GAMA THAT
C                                    OCCURS IN THE NEW SEARCH
C                                    DIRECTION, AND SET SUM TO A SCALAR
C                                    PRODUCT FOR THE TEST BELOW. THE
C                                    VALUE OF GAMDEN IS SET BY THE
C                                    RESTART PROCEDURE.
      GAMA = 0.0
      SUM = 0.0
      DO 145 I=1,N
         GAMA = GAMA+G(I)*W(IRSDG+I)
  145 SUM = SUM+G(I)*W(IRSDX+I)
      GAMA = GAMA/GAMDEN
C                                  RESTART IF THE NEW SEARCH DIRECTION
C                                    IS NOT SUFFICIENTLY DOWNHILL.
C
      IF (ABS(BETA*GMIN+GAMA*SUM).GE.0.2*GSQRD) GO TO 155
C
C                                  CALCULATE THE NEW SEARCH DIRECTION.
      DO 150 I=1,N
  150 W(I) = -G(I)+BETA*W(I)+GAMA*W(IRSDX+I)
      GO TO 80
C                                  APPLY THE RESTART PROCEDURE.
  155 GAMDEN = GMIN-GINIT
      DO 160 I=1,N
         W(IRSDX+I) = W(I)
         W(IRSDG+I) = G(I)-W(IGINIT+I)
  160 W(I) = -G(I)+BETA*W(I)
      ITERRS = ITERC
      GO TO 80
C                                  SET IER TO INDICATE THAT THE SEARCH
C                                    DIRECTION IS UPHILL.
  165 IER = 130
C                                  ENSURE THAT F, X AND G ARE OPTIMAL.
  170 IF (NCALLS.EQ.NFOPT) GO TO 180
      F = FMIN
      DO 175 I=1,N
         X(I) = W(IXOPT+I)
  175 G(I) = W(IGOPT+I)
  180 IF (IER.EQ.0) GO TO 125
 9000 CONTINUE
      CALL UERTST (IER,6HZXCGR )
 9005 RETURN
      END
