C   IMSL ROUTINE NAME   - VSORA
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1982
C
C   PURPOSE             - SORTING OF COLUMNS OF A REAL MATRIX INTO
C                           ASCENDING ORDER OF KEYS IN ROWS.
C
C   USAGE               - CALL VSORA (A,IA,NR,NC,NK,WK,IER)
C
C   ARGUMENTS    A      - INPUT/OUTPUT NR BY NC MATRIX
C                           ON INPUT, A IS THE MATRIX TO BE SORTED.
C                           ON OUTPUT, THE FIRST NC COLUMNS OF A ARE
C                             SORTED IN ASCENDING ORDER OF KEYS IN
C                             THE FIRST NK ROWS.
C                IA     - ROW DIMENSION OF MATRIX A EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM. (INPUT)
C                NR     - INPUT NUMBER OF ROWS OF A.
C                NC     - INPUT NUMBER OF COLUMNS OF A.
C                NK     - INPUT NUMBER OF ROWS TO BE USED AS KEYS.
C                           IF NK IS POSITIVE, THE FIRST ROW IS THE
C                             MOST SIGNIFICANT KEY AND ROW NK IS THE
C                             LEAST SIGNIFICANT.
C                           IF NK IS NEGATIVE, ROW ABS(NK) IS THE MOST
C                             SIGNIFICANT KEY, AND THE FIRST ROW IS THE
C                             LEAST SIGNIFICANT.
C                WK     - WORK VECTOR OF LENGTH (2*NR).
C                IER    - OUTPUT ERROR PARAMETER.
C                           TERMINAL ERROR
C                             IER = 129 INDICATES THAT NR OR NC IS
C                             LESS THAN 0 OR THAT NK = 0. NO SORTING
C                             IS PERFORMED.
C                           WARNING WITH FIX ERROR
C                             IER = 65 INDICATES THAT ABS(NK) IS LARGER
C                             THAN NR, OR THAT NR AND/OR ABS(NK) IS
C                             LARGER THAN IA. IF NR IS LARGER THAN IA,
C                             NR IS LIMITED TO IA ROWS. IF ABS(NK) IS
C                             LARGER THAN NR, NK IS LIMITED TO NR KEYS.
C
C   REQD. IMSL ROUTINES - UERTST,UGETIO,VSCMP,VBLA=SCOPY
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS      IF NK IS POSITIVE, THE COLUMNS OF A ARE SORTED SO THAT
C                FOR K=1,2,...,NC-1,
C                  1.  A(1,K).LE.A(1,K+1), AND
C                  2.  FOR J=2,3,...,NK, IF A(I,K).EQ.A(I,K+1)
C                      FOR I=1,2,...,J-1, THEN A(J,K).LE.A(J,K+1).
C                IF NK IS NEGATIVE, THE COLUMNS OF A ARE SORTED SO THAT
C                FOR K=1,2,...,NC-1,
C                  1.  A(-NK,K).LE.A(-NK,K+1), AND
C                  2.  FOR J=1,2,...,-NK-1, IF A(I,K).EQ.A(I,K+1)
C                      FOR I=J+1,J+2,...,-NK, THEN A(J,K).LE.A(J,K+1).
C
C-----------------------------------------------------------------------
C
      SUBROUTINE VSORA (A,IA,NR,NC,NK,WK,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            IA,NR,NC,NK,IER
      REAL               A(IA,1),WK(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,IJ,IL(21),INCA,INCIR,INCW,INKA,INKW,IROW1,
     *                   IROW2,IROWA,IROWB,IU(21),J,JCOL1,JCOL2,JCOLA,
     *                   JCOLB,K,KOMPAR,KR,L,M,NK1,NR1
      REAL               DR,R,R1,UR
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
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
   10 IF (IABS(NK1).LE.NR1) GO TO 15
      IF (NK1.LT.0) NK1 = -NR1
      IF (NK1.GT.0) NK1 = NR1
      IER = 65
C
   15 JCOL1 = 1
      JCOL2 = NC
      IROW1 = 1
      IROW2 = NK1
      INCIR = 1
      IF (IROW2.LT.0) IROW1 = -NK1
      IF (IROW2.LT.0) INCIR = -1
      IF (IROW2.LT.0) IROW2 = 1
      JCOLA = MIN0(JCOL1,JCOL2)
      JCOLB = MAX0(JCOL1,JCOL2)
      IROWA = MIN0(IROW1,IROW2)
      IROWB = MAX0(IROW1,IROW2)
      KR = MAX0((IROW2-IROW1+INCIR)/INCIR,0)
C                                  INITIALIZATION OF POINTERS FOR VBLA.
C                                    FOR THE PARTICULAR CASE WHERE KR=1
C                                    AND INCIR .LT. 0, OFFSET ROW MUST
C                                    BE SET FOR THE FIRST COMPARISON
C                                    ELEMENT DUE TO THE SPECIAL
C                                    HANDLING OF NEGATIVE INCREMENTS BY
C                                    VBLAX=VSCMP.
      IF (KR.EQ.1 .AND. INCIR.LT.0) IROWA = IROWB
      M = 1
      I = JCOLA
      J = JCOLB
      INCW = 1
      INKW = INCIR
      INCA = 1
      INKA = INCIR
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
C                                  SELECT A CENTRAL COLUMN OF A SAVE IT
C                                    IN VECTOR WK.
      IJ = I+(J-I)*R
      CALL SCOPY(NR1,A(1,IJ),INCA,WK(1),INCW)
C                                  IF 1ST COLUMN OF MATRIX IS GREATER
C                                    THAN VECTOR WK, INTERCHANGE WITH
C                                    WK
      CALL VSCMP(KR,A(IROWA,I),INKA,WK(IROWA),INKW,KOMPAR)
      IF (KOMPAR.EQ.-1 .OR. KOMPAR.EQ.0) GO TO 35
      CALL SCOPY(NR1,A(1,I),INCA,A(1,IJ),INCA)
      CALL SCOPY(NR1,WK(1),INCW,A(1,I),INCA)
      CALL SCOPY(NR1,A(1,IJ),INCA,WK(1),INCW)
   35 CONTINUE
      L = J
C                                  IF LAST COLUMN OF MATRIX IS LESS
C                                    THAN VECTOR WK, INTERCHANGE WITH
C                                    WK
      CALL VSCMP(KR,A(IROWA,J),INKA,WK(IROWA),INKW,KOMPAR)
      IF (KOMPAR.EQ.1 .OR. KOMPAR.EQ.0) GO TO 45
      CALL SCOPY(NR1,A(1,J),INCA,A(1,IJ),INCA)
      CALL SCOPY(NR1,WK(1),INCW,A(1,J),INCA)
      CALL SCOPY(NR1,A(1,IJ),INCA,WK(1),INCW)
C                                  IF 1ST COLUMN OF MATRIX IS GREATER
C                                    THAN VECTOR WK, INTERCHANGE WITH
C                                    WK
      CALL VSCMP(KR,A(IROWA,I),INKA,WK(IROWA),INKW,KOMPAR)
      IF (KOMPAR.EQ.-1 .OR. KOMPAR.EQ.0) GO TO 45
      CALL SCOPY(NR1,A(1,I),INCA,A(1,IJ),INCA)
      CALL SCOPY(NR1,WK(1),INCW,A(1,I),INCA)
      CALL SCOPY(NR1,A(1,IJ),INCA,WK(1),INCW)
      GO TO 45
   40 CONTINUE
      CALL VSCMP(KR,A(IROWA,L),INKA,A(IROWA,K),INKA,KOMPAR)
      IF (KOMPAR.EQ.0) GO TO 45
      CALL SCOPY(NR1,A(1,L),INCA,WK(NR1+1),INCW)
      CALL SCOPY(NR1,A(1,K),INCA,A(1,L),INCA)
      CALL SCOPY(NR1,WK(NR1+1),INCW,A(1,K),INCA)
C                                  FIND A COLUMN IN THE SECOND HALF OF
C                                    THE MATRIX WHICH IS LESS THAN WK.
   45 CONTINUE
      L = L-1
      CALL VSCMP(KR,A(IROWA,L),INKA,WK(IROWA),INKW,KOMPAR)
      IF (KOMPAR.EQ.1) GO TO 45
C                                  FIND A COLUMN IN THE SECOND HALF OF
C                                    THE MATRIX WHICH IS GREATER THAN
C                                    WK.
   50 CONTINUE
      K = K+1
      CALL VSCMP(KR,A(IROWA,K),INKA,WK(IROWA),INKW,KOMPAR)
      IF (KOMPAR.EQ.-1) GO TO 50
C                                  INTERCHANGE THESE COLUMNS
      IF (K.LE.L) GO TO 40
C                                  SAVE UPPER AND LOWER SUBSCRIPTS OF
C                                    SECTION OF A REMAINING TO BE
C                                    SORTED.
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
C                                  BEGIN AGAIN ON ANOTHER PORTION OF
C                                    THE UNSORTED MATRIX.
   60 CONTINUE
      M = M-1
      IF (M.EQ.0) GO TO 9000
      I = IL(M)
      J = IU(M)
   65 CONTINUE
      IF ((J-I).GE.11) GO TO 30
      IF (I.EQ.JCOLA) GO TO 20
      I = I-1
   70 CONTINUE
      I = I+1
      IF (I.EQ.J) GO TO 60
      CALL SCOPY(NR1,A(1,I+1),INCA,WK(1),INCW)
      CALL VSCMP(KR,A(IROWA,I),INKA,WK(IROWA),INKW,KOMPAR)
      IF (KOMPAR.EQ.-1 .OR. KOMPAR.EQ.0) GO TO 70
      K = I
   75 CONTINUE
      CALL SCOPY(NR1,A(1,K),INCA,A(1,K+1),INCA)
      K = K-1
      CALL VSCMP(KR,WK(IROWA),INKW,A(IROWA,K),INKA,KOMPAR)
      IF (KOMPAR.EQ.-1) GO TO 75
      CALL SCOPY(NR1,WK(1),INCW,A(1,K+1),INCA)
      GO TO 70
 9000 CONTINUE
      IF (IER.NE.0) CALL UERTST(IER,6HVSORA )
 9005 RETURN
      END
