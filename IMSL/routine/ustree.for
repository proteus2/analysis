C   IMSL ROUTINE NAME   - USTREE
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1982
C
C   PURPOSE             - PRINT A BINARY TREE (WHICH MAY REPRESENT THE
C                           OUTPUT OF A CLUSTERING ALGORITHM IN
C                           CHAPTER O)
C
C   USAGE               - CALL USTREE (ND,ICLSON,ICRSON,CLEVEL,IND,
C                           XSIM,IOUT,CLVLSK,NCLRST,LEFTRT,STARST,IER)
C
C   ARGUMENTS    ND     - INPUT NUMBER OF DATA POINTS OR TERMINAL NODES
C                           (WHICH MUST BE NUMBERED FROM 1 TO ND).
C                           THE ND-1 NON-TERMINAL NODES MUST BE NUMBERED
C                           FROM ND+1 TO ND+(ND-1).
C                           ND IN THE RANGE 2,3,...,500 IS ALLOWED.
C                ICLSON - INPUT VECTOR OF LENGTH ND.  THE FIRST ND-1
C                           LOCATIONS CONTAIN THE LEFTSON NODES.
C                           LOCATION ND IS USED AS WORK STORAGE.
C                           NODE NUMBER ND+K HAS LEFTSON NODE GIVEN BY
C                           ICLSON(K) AND RIGHTSON NODE GIVEN BY
C                           ICRSON(K), FOR K = 1 TO ND-1.
C                ICRSON - INPUT VECTOR OF RIGHTSON NODES OF LENGTH
C                           ND-1.  SEE ICLSON DESCRIPTION ABOVE.
C                CLEVEL - INPUT VECTOR OF LENGTH ND.  THE FIRST ND-1
C                           LOCATIONS CONTAIN THE SIMILARITY LEVELS.
C                           LOCATION ND IS USED AS WORK STORAGE.
C                           NODE ND+K IS PLOTTED ON A VERTICAL SCALE AT
C                           THE VALUE CLEVEL(K), FOR K=1 TO ND-1.
C                IND    - INPUT VECTOR OF LENGTH 4. IND(I) CONTAINS WHEN
C                           I=1, THE HEAD NODE OF THE SUBTREE TO BE
C                             PRINTED.  MUST BE BETWEEN ND AND 2*ND,
C                             EXCLUSIVELY. SEE REMARKS.
C                           I=2, NUMBER OF PRINTABLE SPACES PER LINE ON
C                             THE OUTPUT (PRINTER) DEVICE. IND(2) MUST
C                             EXCEED 4.
C                           I=3, NUMBER OF HORIZONTAL SLICES OF TREE
C                             DESIRED TO PROVIDE THE NECESSARY DETAIL.
C                             IND(3) MUST BE POSITIVE.
C                           I=4, NUMBER OF FILLER LINES PRINTED BETWEEN
C                             NODE LINES (MUST BE NONNEGATIVE, 1 IS
C                             USUALLY SUFFICIENT).
C                XSIM   - INPUT VECTOR OF LENGTH 2 CONTAINING THE
C                           INTERVAL ON THE VERTICAL SCALE USED TO PLOT
C                           THE TREE (SEE VECTOR CLEVEL). LEVEL XSIM(1)
C                           IS WHERE THE TERMINAL NODES (FOR IMSL
C                           CLUSTERING ROUTINES, THE DATA POINTS) ARE
C                           PRINTED. THE OTHER INTERVAL ENDPOINT XSIM(2)
C                           SHOULD INCLUDE THE LEVEL FOR THE HEAD
C                           NODE (IND(1)). SEE REMARKS.
C                IOUT   - WORK VECTOR OF LENGTH IND(2)-4.
C                CLVLSK - WORK VECTOR OF LENGTH ND.
C                NCLRST - WORK VECTOR OF LENGTH ND.
C                LEFTRT - WORK VECTOR OF LENGTH ND.
C                STARST - WORK VECTOR OF LENGTH ND.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER=129 INDICATES ONE OF IND(1), IND(2),
C                             IND(3), OR ND WAS INCORRECTLY SPECIFIED.
C                           IER=130 INDICATES A REVERSAL WAS FOUND.
C                             CHECK INPUT ARRAYS. SEE REMARKS.
C                           IER=131 INDICATES ALGORITHM DETECTED
C                             ABNORMALITY PROBABLY DUE TO INCORRECT NODE
C                             NUMBERING.
C
C   REQD. IMSL ROUTINES - UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS  1.  THE TREE MAY BE TOO LARGE TO FIT IN IND(2) SPACES
C                REPRESENTING THE INTERVAL (XSIM(1),XSIM(2)). IF SO,
C                USTREE CAN PRINT THE TREE IN IND(3) SLICES OF THIS
C                INTERVAL WHICH MAY BE CUT AND TAPED TOGETHER.
C            2.  TO PRINT THE ENTIRE TREE FROM IMSL SUBROUTINE OCLINK,
C                THE HEAD NODE IND(1) = ND+(ND-1).
C            3.  OUTPUT IS WRITTEN TO THE UNIT SPECIFIED BY IMSL
C                ROUTINE UGETIO. SEE THE UGETIO DOCUMENT FOR DETAILS.
C            4.  REVERSALS (IER=130) MAY OCCUR IN TWO WAYS. FIRST, TWO
C                NODES MAY BE JOINED AT A LOWER LEVEL (CLOSER TO
C                XSIM(1)). SECOND, THE LEVEL OF THE HEAD NODE MAY LIE
C                ABOVE THE INTERVAL (XSIM(1),XSIM(2)).
C            5.  FOR PROPER DISPLAY, THE TREE CREATED BY USTREE SHOULD
C                BE TURNED TO AN UPRIGHT POSITION.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE USTREE (ND,ICLSON,ICRSON,CLEVEL,IND,XSIM,IOUT,CLVLSK,
     *                   NCLRST,LEFTRT,STARST,IER)
C
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            ND,IER,ICLSON(1),ICRSON(1),IND(4),IOUT(1),
     *                   NCLRST(1),LEFTRT(1)
      REAL               CLEVEL(1),XSIM(2),CLVLSK(1),STARST(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,IBACK,IBLANK,IFMT(9),IHEAD,II4A,II4B,IOPT,
     *                   IPLUS,IRET,ISERCH,ISTAR,ITOP,ITOPST,ITOPSV,J,K,
     *                   K1,K1P,K2,KTHRED,L,N,NC,NDATA,NIN,NODE,NOUT,NP,
     *                   NPAGES,NREPS,NSP4,NSPACE,NUM(10)
      REAL               CLVL,CLVLM1,CONS,EPS,RANGE,RHISIM,RLOSIM,RSIGN,
     *                   RSIGNE,XHISIM,XLOSIM
      DATA               NUM /1H0,1H1,1H2,1H3,1H4,1H5,1H6,1H7,1H8,1H9/
      DATA               ISTAR /1H*/,IBLANK /1H /,IPLUS /1H+/
      DATA               IFMT /1H(,1H ,1H ,2HX,,1H ,1H ,1H ,2HA1,1H)/
      DATA               II4A /2HI4/,II4B /2H,1/
C                                  FIRST EXECUTABLE STATEMENT
      CALL UGETIO(1,NIN,NOUT)
      NC = ND-1
      IHEAD = IND(1)
      NSPACE = IND(2)
      NPAGES = IND(3)
      NREPS = IND(4)
      XHISIM = XSIM(1)
      XLOSIM = XSIM(2)
      IER = 0
C                                  ERROR CHECKS
      IF (NSPACE.LE.4 .OR. ND.GT.500 .OR. IHEAD.LE.ND .OR.
     *IHEAD.GE.ND+ND .OR. NPAGES.LE.0 .OR. ND.LT.2) GO TO 200
C                                  SOME SIGN INFORMATION FOR SCALES
      RSIGN = 1.0
      RSIGNE = 1.0001
      IF (XHISIM.LT.XLOSIM) GO TO 5
      RSIGN = -1.0
      RSIGNE = -0.9999
    5 NSP4 = NSPACE-4
C                                  INITIALIZE PAGE SPLITTING OPTIONS
      K = NSP4
      K1 = 1
      IRET = 2
      GO TO 115
   10 DO 15 I=1,J
   15 IFMT(4+I) = IOUT(I)
      RANGE = (XLOSIM-XHISIM)/NPAGES
      NP = 0
   20 NP = NP+1
      IF (NP.GT.NPAGES) GO TO 9005
      RHISIM = XLOSIM-RANGE*NP
      RLOSIM = RHISIM+RANGE
      WRITE (NOUT,25) RHISIM, RLOSIM
   25 FORMAT (22H1SIMILARITY RANGE FROM, F10.4, 3H TO, F10.4/)
      IFMT(2) = NUM(6)
      WRITE (NOUT,IFMT) (IPLUS,I=1,NSP4)
C                                  INITIALIZE STACKS AND OTHER VARIABLES
      NODE = IHEAD-ND
      ITOP = 1
      CLVLSK(ITOP) = XLOSIM
      NCLRST(ITOP) = ND
      ITOPST = 0
C                                  SET DUMMY CLUSTER LEFTSON TO IHEAD
      ICLSON(ND) = IHEAD
      CLEVEL(ND) = XLOSIM
C                                  INITIALLY BLANK OUT IOUT
      CONS = NSP4/(RLOSIM-RHISIM)
      DO 30 I=1,NSP4
   30 IOUT(I) = IBLANK
C                                  BEGIN TREE WALK AND OUTPUT
C                                  STACK CURRENT NODE FOR NEXT WALK
   35 ITOP = ITOP+1
      CLVLSK(ITOP) = CLEVEL(NODE)
C                                  CHECK FOR REVERSALS
C                                  RESIGN GIVES VERTICAL SCALE DIRECTION
      IF (RSIGN*CLVLSK(ITOP).GT.RSIGNE*CLVLSK(ITOP-1)) GO TO 195
      NCLRST(ITOP) = NODE
C                                  SEARCH ND ICLSONS AND ND-1 ICRSONS
C                                  FOR VALUE = NODE+ND.  ISERCH = 0 FOR
C                                  LEFTSON, ISERCH = 1 FOR RIGHTSON
      K = NODE+ND
      ISERCH = 0
      DO 40 I=1,ND
         IF (IABS(ICLSON(I)).EQ.K) GO TO 50
   40 CONTINUE
      ISERCH = 1
      DO 45 I=1,NC
         IF (ICRSON(I).EQ.K) GO TO 50
   45 CONTINUE
      GO TO 190
   50 LEFTRT(ITOP) = ISERCH
C                                  NATURE OF LEFTSON OF CURRENT NODE
      IF (ICLSON(NODE).LE.ND) GO TO 55
C                                  LEFTSON IS CLUSTER - MOVE AND STACK
      NODE = ICLSON(NODE)-ND
      GO TO 35
C                                  LEFTSON IS DATA. PRINT DATA, CLUSTER
   55 NDATA = ICLSON(NODE)
      GO TO 65
C                                  RIGHTSON DATA. PRINT DATA, CLUSTER
   60 NDATA = ICRSON(NODE)
C                                  PRINT DATA LINE UP TO LEVEL-(ITOP)
   65 CLVL = CLVLSK(ITOP)
      K = 1.0001+(CLVL-RHISIM)*CONS
      IF (K.LT.1) GO TO 75
      K = MIN0(K,NSP4)
      DO 70 I=1,K
   70 IOUT(I) = ISTAR
   75 IF (NP.NE.NPAGES) GO TO 80
      IFMT(2) = II4A
      IFMT(3) = II4B
      WRITE (NOUT,IFMT) NDATA, (IOUT(I),I=1,NSP4)
      IFMT(3) = IBLANK
      GO TO 85
   80 IFMT(2) = NUM(6)
      WRITE (NOUT,IFMT) (IOUT(I),I=1,NSP4)
C                                  ADJUST STAR-STACK
C                                  FIND THREAD THROUGH MARK ON LEFTSONS
C                                  CHECK IF DATA IS LEFTSON - (RETURN)
   85 KTHRED = 0
      IF (NDATA.EQ.ICLSON(NODE)) GO TO 95
C                                  DATA IS RIGHTSON
C                                  FIND UNPRINTED CLUSTER IN STACK
      K = ITOP+1
   90 K = K-1
      IF (K.LT.1) GO TO 190
      L = NCLRST(K)
      IF (ICLSON(L).LT.0) GO TO 90
      KTHRED = L
      ITOPSV = K
C                                  PRINT FILLER LINE
   95 IBACK = 1
      CLVL = CLVLSK(ITOP)
      IOPT = 0
      IF (KTHRED.EQ.0) GO TO 145
      IF (RSIGN*CLEVEL(KTHRED).GT.RSIGN*CLEVEL(IHEAD-ND)) GO TO 205
      IOPT = 1
      GO TO 145
C                                  FOLLOW THREAD TO NODE - POP STACKS
C                                  ONLY IF THE THREAD IS NOT ZERO
  100 IF (KTHRED.EQ.0) GO TO 105
      NODE = KTHRED
C                                  POP STACK DOWN TO CLUSTER LEVEL NODE
      ITOP = ITOPSV
C                                  PRINT CLUSTER NODE BETWEEN LEVELS
  105 CLVL = CLVLSK(ITOP)
      CLVLM1 = CLVLSK(ITOP-1)
      K1 = 1.0001+(CLVL-RHISIM)*CONS
      K2 = 1.0001+(CLVLM1-RHISIM)*CONS
      IF (K1.GT.NSP4 .OR. K2.LT.1) GO TO 135
      K1P = MAX0(K1,1)
      K2 = MIN0(K2,NSP4)
      DO 110 I=K1P,K2
  110 IOUT(I) = ISTAR
      IRET = 1
      K = NODE+ND
  115 J = 0
      L = 1000
      DO 130 I=1,3
         L = L/10
         N = K/L
         IF (N.NE.0) GO TO 120
         IF (J.EQ.0) GO TO 130
  120    IF (K1+J.LT.1) GO TO 125
         IOUT(K1+J) = NUM(N+1)
  125    J = J+1
         K = K-N*L
  130 CONTINUE
      IF (IRET.EQ.2) GO TO 10
  135 IFMT(2) = NUM(6)
      WRITE (NOUT,IFMT) (IOUT(I),I=1,NSP4)
C                                  AFTER CLUSTER PRINTED, MARK LEFTSON
      ICLSON(NODE) = -ICLSON(NODE)
C                                  ADJUST STAR-STACK
      IBACK = 2
      CLVL = CLVLSK(ITOP-1)
      IOPT = LEFTRT(ITOP)
      GO TO 145
C                                  EXAMINE RIGHTSON OF CLUSTER
  140 IF (ICRSON(NODE).LE.ND) GO TO 60
      NODE = ICRSON(NODE)-ND
      GO TO 35
C                                  ADD OR DELETE LEVEL FROM *-STACK
C                                  A BLANK SPOT IS -1.E10
  145 IF (IOPT.EQ.1) GO TO 160
C                                  ADD CLVL TO *-STACK
C                                  DO NOT STACK XLOSIM
      IF (RSIGN*CLVL.GE.XLOSIM/RSIGNE) GO TO 170
      IF (ITOPST.EQ.0) GO TO 155
      DO 150 I=1,ITOPST
         IF (STARST(I).GT.-1.E9) GO TO 150
         STARST(I) = CLVL
         GO TO 170
  150 CONTINUE
  155 ITOPST = ITOPST+1
      STARST(ITOPST) = CLVL
      GO TO 170
C                                  DELETE CLVL FROM *-STACK
  160 EPS = 1.E-5*AMAX1(ABS(CLVL),1.E-3)
      DO 165 I=1,ITOPST
C                                  USE RELATIVE COMPARISON TO DELETE
         IF (ABS(CLVL-STARST(I)).GT.EPS) GO TO 165
         STARST(I) = -1.E10
         GO TO 170
  165 CONTINUE
      GO TO 190
C                                  PRINT FILLER LINE
C                                  BLANKS OUT AND PUTS IN STAR STACK
  170 DO 175 I=1,NSP4
  175 IOUT(I) = IBLANK
      DO 180 I=1,ITOPST
         IF (STARST(I).LT.-1.E9) GO TO 180
         K = 1.0001+(STARST(I)-RHISIM)*CONS
         IF (K.LT.1 .OR. K.GT.NSP4) GO TO 180
         IOUT(K) = ISTAR
  180 CONTINUE
      IFMT(2) = NUM(6)
      IF (NREPS.EQ.0) GO TO (100, 140), IBACK
      DO 185 I=1,NREPS
  185 WRITE (NOUT,IFMT) (IOUT(J),J=1,NSP4)
      GO TO (100, 140), IBACK
C                                  TERMINAL ERROR ENCOUNTERED
  190 IER = IER+1
  195 IER = IER+1
  200 IER = IER+129
  205 DO 210 I=1,ND
  210 ICLSON(I) = IABS(ICLSON(I))
      IF (IER.NE.0) GO TO 9000
      IFMT(2) = NUM(6)
      WRITE (NOUT,IFMT) (IPLUS,I=1,NSP4)
      GO TO 20
 9000 CONTINUE
      CALL UERTST(IER,6HUSTREE)
 9005 RETURN
      END
