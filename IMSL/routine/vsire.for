C   IMSL ROUTINE NAME   - VSIRE
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1982
C
C   PURPOSE             - NUCLEUS CALLED BY IMSL ROUTINES WHICH SORT
C                           A MATRIX WITH RESPECT TO VECTORS.
C
C   REQD. IMSL ROUTINES - UERTST,UGETIO,VSCMPE,VBLA=SCOPY
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE VSIRE (A,IA,NR,NC,IOPT,NK,WK,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            IA,NR,NC,IOPT,NK,IER
      REAL               A(IA,1),WK(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,IJ,IL(21),INCA,INCJC,INCW,INKA,INKW,IROW1,
     *                   IROW2,IROWA,IROWB,IU(21),J,JCOL1,JCOL2,JCOLA,
     *                   JCOLB,K,KC,KOMPAR,KR,L,M,NC1,NK1,NR1
      REAL               DR,R,R1,UR
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
      NC1 = NC
      NR1 = NR
      NK1 = NK
C                                  CHECK INPUT. IF CONDITIONS ARE NOT
C                                    MET, SET IER AND RETURN WITHOUT
C                                    SORTING.
      IF (.NOT.(NR.LE.0 .OR. NC.LE.0 .OR. NK.EQ.0)) GO TO 5
      IER = 129
      GO TO 9000
C                                  CHECK INPUT. IF CONDITIONS ARE NOT
C                                    MET, WARN, FIX, AND CONTINUE.
    5 IF (NR1.LE.IA) GO TO 10
      NR1 = IA
      IER = 65
   10 IF (IABS(NK1).LE.NC1) GO TO 15
      IF (NK1.LT.0) NK1 = -NC1
      IF (NK1.GT.0) NK1 = NC1
      IER = 65
C
   15 IROW1 = 1
      IROW2 = NR1
      JCOL1 = 1
      JCOL2 = NK1
      INCJC = 1
      IF (JCOL2.LT.0) JCOL1 = -NK1
      IF (JCOL2.LT.0) INCJC = -1
      IF (JCOL2.LT.0) JCOL2 = 1
      IROWA = MIN0(IROW1,IROW2)
      IROWB = MAX0(IROW1,IROW2)
      JCOLA = MIN0(JCOL1,JCOL2)
      JCOLB = MAX0(JCOL1,JCOL2)
      KC = MAX0((JCOL2-JCOL1+INCJC)/INCJC,0)
C
C                                  INITIALIZATION OF POINTERS FOR VBLA.
C                                    FOR THE PARTICULAR CASE WHERE KC=1
C                                    AND INCJC .LT. 0., OFFSET ROW
C                                    MUST BE SET DIFFERENTLY FOR THE
C                                    FIRST COMPARISON ELEMENT DUE TO THE
C                                    SPECIAL HANDLING OF NEGATIVE
C                                    INCREMENTS BY VBLAX=VSCMPE.
C
      IF (KC.EQ.1 .AND. INCJC.LT.0) JCOLA = JCOLB
      M = 1
      I = IROWA
      J = IROWB
      INCW = 1
      INKW = INCJC
      INCA = IA
      INKA = INCJC*IA
      R = 0.375
      R1 = 0.5898437
      UR = 3.90625E-02
      DR = 0.21875
   20 CONTINUE
      IF (I.EQ.J) GO TO 60
      IF (R.GT.R1) GO TO 25
      R = R+UR
      GO TO 30
   25 CONTINUE
      R = R-DR
   30 CONTINUE
      K = I
C                                  SELECT A CENTRAL ROW OF THE
C                                    MATRIX A AND SAVE IT IN VECTOR WK.
      IJ = I+(J-I)*R
      CALL SCOPY(NC1,A(IJ,1),INCA,WK(1),INCW)
C                                  IF 1ST ROW OF MATRIX IS GREATER
C                                    THAN VECTOR WK INTERCHANGE WITH WK.
      CALL VSCMPE(IOPT,KC,A(I,JCOLA),INKA,WK(JCOLA),INKW,KOMPAR)
      IF (KOMPAR.EQ.-1 .OR. KOMPAR.EQ.0) GO TO 35
      CALL SCOPY(NC1,A(I,1),INCA,A(IJ,1),INCA)
      CALL SCOPY(NC1,WK(1),INCW,A(I,1),INCA)
      CALL SCOPY(NC1,A(IJ,1),INCA,WK(1),INCW)
   35 CONTINUE
      L = J
C                                  IF LAST ROW OF MATRIX IS LESS
C                                    THAN WK INTERCHANGE WITH WK.
      CALL VSCMPE(IOPT,KC,A(J,JCOLA),INKA,WK(JCOLA),INKW,KOMPAR)
      IF (KOMPAR.EQ.1 .OR. KOMPAR.EQ.0) GO TO 45
      CALL SCOPY(NC1,A(J,1),INCA,A(IJ,1),INCA)
      CALL SCOPY(NC1,WK(1),INCW,A(J,1),INCA)
      CALL SCOPY(NC1,A(IJ,1),INCA,WK(1),INCW)
C                                  IF 1ST ROW OF MATRIX IS GREATER
C                                    THAN VECTOR WK INTERCHANGE WITH WK.
      CALL VSCMPE(IOPT,KC,A(I,JCOLA),INKA,WK(JCOLA),INKW,KOMPAR)
      IF (KOMPAR.EQ.-1 .OR. KOMPAR.EQ.0) GO TO 45
      CALL SCOPY(NC1,A(I,1),INCA,A(IJ,1),INCA)
      CALL SCOPY(NC1,WK(1),INCW,A(I,1),INCA)
      CALL SCOPY(NC1,A(IJ,1),INCA,WK(1),INCW)
      GO TO 45
   40 CONTINUE
      CALL VSCMPE(IOPT,KC,A(L,JCOLA),INKA,A(K,JCOLA),INKA,KOMPAR)
      IF (KOMPAR.EQ.0) GO TO 45
      CALL SCOPY(NC1,A(L,1),INCA,WK(NC1+1),INCW)
      CALL SCOPY(NC1,A(K,1),INCA,A(L,1),INCA)
      CALL SCOPY(NC1,WK(NC1+1),INCW,A(K,1),INCA)
C                                  FIND A ROW IN THE SECOND HALF OF
C                                    THE MATRIX WHICH IS SMALLER THAN
C                                    WK.
   45 CONTINUE
      L = L-1
      CALL VSCMPE(IOPT,KC,A(L,JCOLA),INKA,WK(JCOLA),INKW,KOMPAR)
      IF (KOMPAR.EQ.1) GO TO 45
C                                  FIND A ROW IN THE SECOND HALF OF
C                                    THE MATRIX WHICH IS GREATER THAN
C                                    WK.
   50 CONTINUE
      K = K+1
      CALL VSCMPE(IOPT,KC,A(K,JCOLA),INKA,WK(JCOLA),INKW,KOMPAR)
      IF (KOMPAR.EQ.-1) GO TO 50
C                                  INTERCHANGE THESE ROWS
      IF (K.LE.L) GO TO 40
C                                  SAVE UPPER AND LOWER SUBSCRIPTS OF
C                                    THE MATRIX A SECTION REMAINING TO
C                                    BE SORTED.
      IF ((L-I).LE.(J-K)) GO TO 55
      IL(M) = I
      IU(M) = L
      I = K
      M = M+1
      GO TO 65
   55 CONTINUE
      IL(M) = K
      IU(M) = J
      J = L
      M = M+1
      GO TO 65
C
C                                  BEGIN AGAIN ON ANOTHER PORTION OF
C                                    THE UNSORTED MATRIX.
   60 CONTINUE
      M = M-1
      IF (M.EQ.0) GO TO 9000
      I = IL(M)
      J = IU(M)
   65 CONTINUE
      IF ((J-I).GE.11) GO TO 30
      IF (I.EQ.IROWA) GO TO 20
      I = I-1
   70 CONTINUE
      I = I+1
      IF (I.EQ.J) GO TO 60
      CALL SCOPY(NC1,A(I+1,1),INCA,WK(1),INCW)
      CALL VSCMPE(IOPT,KC,A(I,JCOLA),INKA,WK(JCOLA),INKW,KOMPAR)
      IF (KOMPAR.EQ.-1 .OR. KOMPAR.EQ.0) GO TO 70
      K = I
   75 CONTINUE
      CALL SCOPY(NC1,A(K,1),INCA,A(K+1,1),INCA)
      K = K-1
      CALL VSCMPE(IOPT,KC,WK(JCOLA),INKW,A(K,JCOLA),INKA,KOMPAR)
      IF (KOMPAR.EQ.-1) GO TO 75
      CALL SCOPY(NC1,WK(1),INCW,A(K+1,1),INCA)
      GO TO 70
 9000 CONTINUE
      IF (IER.NE.0) CALL UERTST(IER,6HVSIRE )
 9005 RETURN
      END
