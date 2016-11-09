C   IMSL ROUTINE NAME   - OCLINK
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1982
C
C   PURPOSE             - PERFORM A SINGLE-LINKAGE OR COMPLETE-LINKAGE
C                           HIERARCHICAL CLUSTER ANALYSIS GIVEN A
C                           SIMILARITY MATRIX
C
C   USAGE               - CALL OCLINK (ND,IOPT,XSIM,CLEVEL,ICLSON,
C                                     ICRSON,IPTR,IER)
C
C   ARGUMENTS    ND     - INPUT NUMBER OF DATA POINTS TO BE CLUSTERED.
C                           ND-1 CLUSTERS ARE FORMED NUMBERED CONSEC-
C                           UTIVELY FROM ND+1 TO ND+(ND-1).  ND MUST
C                           BE GREATER THAN 2. (SEE REMARKS)
C                IOPT   - INPUT OPTIONS VECTOR OF LENGTH 2.
C                           IOPT(1)=0 IMPLIES SINGLE-LINKAGE DESIRED.
C                             OTHERWISE, COMPLETE-LINKAGE IS PERFORMED.
C                           IOPT(2)=0 IMPLIES THE SIMILARITIES ARE
C                             DISTANCE-LIKE (I.E.  SMALLER IMPLIES
C                             CLOSER).  OTHERWISE, THE SIMILARITIES ARE
C                             ASSUMED TO BE POSITIVE CORRELATION-LIKE
C                             SIMILARITIES. (SEE REMARKS).
C                XSIM   - INPUT/OUTPUT VECTOR OF LENGTH ((ND+1)*ND)/2
C                           CONTAINING THE SIMILARITY MATRIX IN SYM-
C                           METRIC STORAGE MODE.  FOR I GREATER THAN J,
C                           XSIM(((I-1)*I)/2 + J) CONTAINS THE SIMILAR-
C                           ITY OF THE I-TH AND J-TH DATA POINTS.
C                           ON INPUT, THE DIAGONAL ELEMENTS (SIMILARITY
C                           TO ITSELF) ARE ARBITRARY AND DO NOT NEED TO
C                           BE DEFINED.  XSIM IS DESTROYED ON OUTPUT.
C                CLEVEL - OUTPUT VECTOR OF LENGTH ND-1.  CLEVEL(K)
C                           CONTAINS THE SIMILARITY LEVEL AT WHICH
C                           CLUSTER ND+K WAS FORMED.
C                ICLSON - OUTPUT VECTOR OF LEFTSONS OF LENGTH ND-1.
C                           CLUSTER NUMBER ND+K WAS FORMED BY MERGING
C                           CLUSTER ICLSON(K) WITH CLUSTER ICRSON(K).
C                ICRSON - OUTPUT VECTOR OF RIGHTSONS OF LENGTH ND-1.
C                           THE RIGHTSON OF CLUSTER ND+K IS CONTAINED
C                           IN ICRSON(K).
C                IPTR   - WORK VECTOR OF LENGTH ND.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER = 129 INDICATES ND WAS LESS THAN 3.
C
C   REQD. IMSL ROUTINES - UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS  1.  THE DATA CLUSTERS ARE NUMBERED 1 TO ND. THE ND-1
C                CLUSTERS FORMED BY MERGING ARE NUMBERED ND+1 TO
C                ND+(ND-1) AND DECREASE IN SIMILARITY, MAKING IT EASY
C                TO IDENTIFY THE MOST SIMILAR CLUSTERS.
C            2.  SIMILARITIES GENERALLY SHOULD BE NONNEGATIVE. RAW
C                CORRELATIONS TAKE ON VALUES R WHERE R LIES IN THE
C                CLOSED INTERVAL (-1, 1) AND SHOULD BE MADE POSITIVE.
C                IF R = 1 AND R = -1 BOTH MEAN HIGH SIMILARITY, THEN THE
C                TRANSFORMATIONS, RR EQUAL TO THE ABSOLUTE VALUE OF R,
C                OR, RR EQUAL TO R SQUARED, ARE APPROPRIATE. IF R = -1
C                REPRESENTS VERY LOW SIMILARITY, THEN THE TRANSFORMA-
C                TION, RR = 1-R, BECOMES A DISTANCE-LIKE SIMILARITY.
C            3.  NOTE THAT THE ORIGINAL DATA MATRIX OF THE USER DOES
C                NOT ENTER OCLINK, ALLOWING THE USER TO DEFINE, VIA
C                XSIM, WHATEVER MEASURE OF SIMILARITY SEEMS MOST
C                APPROPRIATE.
C            4.  A COMPUTER PRINTING OF A BINARY TREE IS AVAILABLE
C                THROUGH IMSL ROUTINE USTREE.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE OCLINK (ND,IOPT,XSIM,CLEVEL,ICLSON,ICRSON,IPTR,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            ND,IOPT(2),ICLSON(1),ICRSON(1),IPTR(ND),IER
      REAL               XSIM(1),CLEVEL(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            IOPT1,ISGN,IR,IM1,I,J,INUM,JSAVE,ICL,ISTAR,
     *                   JSTAR,ISM1,JSM1,ISP1,JSP1,ITYPE,JTYPE,IRR,IP,
     *                   JP,IUP,ILO,IS,KSTAR
      REAL               ZERO,TENTH,ONEP1,XBIG,XMIN,TEMP
      DATA               ZERO,TENTH,ONEP1,XBIG/0.0,0.1,1.1,1.E9/
C                                  ERROR CHECK
C                                  FIRST EXECUTABLE STATEMENT
      IF (ND.LT.3) GO TO 135
      IER = 0
      IOPT1 = IOPT(1)
C                                  PUT CLUSTER DESIGNATION ON DIAGONAL
C                                  CHANGE SIGNS OF SIMILARITIES
      ISGN = 1
      IF (IOPT(2).NE.0) ISGN = -ISGN
      XSIM(1) = ONEP1
      IR = 1
      IM1 = 1
      DO 15 I=2,ND
         IF (ISGN.EQ.1) GO TO 10
         DO 5 J=1,IM1
    5    XSIM(IR+J) = -XSIM(IR+J)
   10    IR = IR+I
         IM1 = I
   15 XSIM(IR) = FLOAT(I)+TENTH
C                                  IPTR(I) CONTAINS THE COLUMN WITH
C                                    THE SMALLEST SIMILARITY IN ROW I
C                                    (IPTR(1) UNDEFINED-ROW 1 IS NULL)
      IPTR(2) = 1
C                                  INUM POINTS TO ROW WITH SMALLEST
C                                    SIMILARITY.
      INUM = 2
      XMIN = XSIM(2)
      IR = 3
      IM1 = 2
      DO 25 I=3,ND
         TEMP = XSIM(IR+1)
         JSAVE = 1
         DO 20 J=2,IM1
            IF (TEMP.LT.XSIM(IR+J)) GO TO 20
            TEMP = XSIM(IR+J)
            JSAVE = J
   20    CONTINUE
         IM1 = I
         IR = IR+I
         IPTR(I) = JSAVE
         IF (XMIN.LT.TEMP) GO TO 25
         INUM = I
         XMIN = TEMP
   25 CONTINUE
      ICL = 0
C                                  MERGE THE TWO CLUSTERS CLOSEST
C                                    TOGETHER.  THEY ARE IN ROWS
C                                    ISTAR AND JSTAR.
   30 ICL = ICL+1
      ISTAR = INUM
      JSTAR = IPTR(INUM)
      ISM1 = ((ISTAR-1)*ISTAR)/2
      JSM1 = ((JSTAR-1)*JSTAR)/2
      ISP1 = ((ISTAR+1)*ISTAR)/2
      JSP1 = ((JSTAR+1)*JSTAR)/2
C                                  UPDATE LEFT AND RIGHT SON ARRAYS
C                                  RECOVER CURRENT CLUSTER NUMBER
C                                    ON DIAGONAL.
      ITYPE = XSIM(ISP1)
      JTYPE = XSIM(JSP1)
      IF (ITYPE.GT.ND.AND.JTYPE.GT.ND) GO TO 35
C                                  AT LEAST ONE SON IS A PURE DATA
C                                    CLUSTER.  MAKE RIGHTSON DATA ALWAYS
      ICLSON(ICL) = MAX0(ITYPE,JTYPE)
      ICRSON(ICL) = MIN0(ITYPE,JTYPE)
      GO TO 45
C                                  BOTH SONS ARE HIGH LEVEL CLUSTERS.
C                                  MAKE LEFTSON MORE SIMILAR OF TWO
   35 IF (CLEVEL(JTYPE-ND).LT.CLEVEL(ITYPE-ND)) GO TO 40
      ICLSON(ICL) = ITYPE
      ICRSON(ICL) = JTYPE
      GO TO 45
   40 ICLSON(ICL) = JTYPE
      ICRSON(ICL) = ITYPE
C                                  RECORD LEVEL OF MERGING
   45 CLEVEL(ICL) = XSIM(ISM1+JSTAR)
      IF (ICL.GE.ND-1) GO TO 125
C                                  MARK MERGED CLUSTERS-DELETE ROW ISTAR
      XSIM(JSP1) = FLOAT(ND+ICL)+TENTH
      XSIM(ISP1) = -TENTH
C                                  UPDATE SIMILARITY MATRIX
      IR = 0
      DO 75 I=1,ND
         IRR = IR
         IR = IR+I
C                                  SKIP CURRENT AND DELETED CLUSTERS
         IF (XSIM(IR).LT.ZERO) GO TO 75
         IF (I.EQ.JSTAR) GO TO 75
C                                  FIND XSIM(ISTAR,I) AND XSIM(JSTAR,I)
         IF (I.GT.ISTAR) GO TO 50
         IP = ISM1+I
         GO TO 55
   50    IP = IRR+ISTAR
         GO TO 60
   55    IF (I.GT.JSTAR) GO TO 60
         JP = JSM1+I
         GO TO 65
   60    JP = IRR+JSTAR
C                                  SINGLE/COMPLETE LINKAGE OPTION
   65    IF (IOPT1.EQ.0) GO TO 70
         XSIM(JP) = AMAX1(XSIM(IP),XSIM(JP))
         GO TO 75
   70    XSIM(JP) = AMIN1(XSIM(IP),XSIM(JP))
   75 CONTINUE
C                                  UPDATE IPTR ARRAY, FINDING NEW ROW
C                                    MINIMUMS AND OVERALL MINIMUM ALSO
      XMIN = XBIG
C                                  FIRST DO THE NEW CLUSTER ROW, JSTAR.
C                                  IF JSTAR IS LESS THAN 3, THEN
C                                    IPTR ARRAY IS OK.
      IF (JSTAR.LT.3) GO TO 85
      TEMP = XBIG
      IUP = JSTAR-1
      IR = 0
      DO 80 J=1,IUP
         IR = IR+J
         IF (XSIM(IR).LT.ZERO) GO TO 80
         JP = JSM1+J
         IF (TEMP.LT.XSIM(JP)) GO TO 80
         TEMP = XSIM(JP)
         IPTR(JSTAR) = J
   80 CONTINUE
C                                  FIND OVERALL MIN OF 1ST JSTAR ROWS
   85 IF (JSTAR.EQ.1) GO TO 95
      IR = 1
      DO 90 I=2,JSTAR
         IRR = IR
         IR = IR+I
         IF (XSIM(IR).LT.ZERO) GO TO 90
         IP = IRR+IPTR(I)
         IF (XMIN.LT.XSIM(IP)) GO TO 90
         XMIN = XSIM(IP)
         INUM = I
   90 CONTINUE
C                                  UPDATE SIMILARITIES IN ROWS AFTER
C                                    JSTAR.
   95 ILO = JSTAR+1
      IR = ((ILO-1)*ILO)/2
C                                  SINCE JSTAR .LT. ISTAR, ILO .LE. ND
C                                    IS ASSURED
      IM1 = ILO-1
      DO 120 I=ILO,ND
         IS = IR
         IR = IR+I
         IF (XSIM(IR).LT.ZERO) GO TO 120
         KSTAR = IPTR(I)
C                                  IF KSTAR=ISTAR, THE LAST MINIMUM
C                                    WAS DELETED-SEARCH WHOLE ROW
         IF (KSTAR.EQ.ISTAR) GO TO 105
         IF (IOPT1.EQ.0) GO TO 100
C                                  COMPLETE LINKAGE OPTION.  VALUE IN
C                                    COLUMN JSTAR WAS INCREASED.  IF
C                                    KSTAR=JSTAR, SEARCH WHOLE ROW.
C                                    OTHERWISE IPTR AND MIN ARE OK.
         IF (KSTAR.EQ.JSTAR) GO TO 105
         GO TO 115
C                                  SINGLE LINKAGE OPTION.  VALUE IN
C                                    COLUMN JSTAR WAS DECREASED.  IF
C                                    KDTAR=JSTAR, IPTR IS OK. OTHERWISE,
C                                    COMPARE OLD MIN KSTAR WITH JSTAR.
  100    IF (KSTAR.EQ.JSTAR) GO TO 115
         IP = IS+KSTAR
         JP = IS+JSTAR
         IF (XSIM(IP).LT.XSIM(JP)) GO TO 115
         IPTR(I) = JSTAR
         GO TO 115
C                                  SEARCH ENTIRE ROW FOR MINIMUM.
  105    TEMP = XBIG
         IRR = 0
         DO 110 J=1,IM1
            IRR = IRR+J
            IF (XSIM(IRR).LT.ZERO) GO TO 110
            JP = IS+J
            IF (TEMP.LT.XSIM(JP)) GO TO 110
            TEMP = XSIM(JP)
            IPTR(I) = J
  110    CONTINUE
C                                  LOOK FOR OVERALL MINIMUM.
  115    IP = IS+IPTR(I)
         IF (XMIN.LT.XSIM(IP)) GO TO 120
         XMIN = XSIM(IP)
         INUM = I
  120 IM1 = I
      GO TO 30
C                                  CHECK CLEVEL USING IOPT AND ISGN.
  125 IF (ISGN.EQ.1) GO TO 9005
      DO 130 I=1,ICL
  130 CLEVEL(I) = -CLEVEL(I)
      GO TO 9005
  135 IER = 129
 9000 CONTINUE
      CALL UERTST(IER,6HOCLINK)
 9005 RETURN
      END
