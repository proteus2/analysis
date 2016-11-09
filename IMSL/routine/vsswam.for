C   IMSL ROUTINE NAME   - VSSWAM
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - NOVEMBER 1, 1984
C
C   PURPOSE             - NUCLEUS CALLED BY IMSL ROUTINES WHICH SORT
C                           A MATRIX WITH RESPECT TO VECTORS.
C
C   REQD. IMSL ROUTINES - VBLA=SSWAP
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE VSSWAM(A,IA,NR,NC,IOP,INCR,INCC)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            IA,NR,NC,IOP,INCR,INCC
      REAL               A(IA,1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,INCA,IX,IY,M,MP1,NH,NS
C                                  FIRST EXECUTABLE STATEMENT
      IF (NR.LE.0 .OR. NC.LE.0) RETURN
      IF (IOP.EQ.1) GO TO 55
C                                  SWAP COLUMNS
      IF (INCC.EQ.INCR) IF (INCC-1) 5, 20, 45
C                                  CODE FOR UNEQUAL OR NONPOSITIVE
C                                    INCREMENTS.
    5 IX = 1
      NH = NC/2
      IF (INCC.LT.0) IX = (-NC+1)*INCC + 1
      IF (NH.EQ.0) GO TO 15
      DO 10 I = 1, NH
        CALL SSWAP(NR,A(1,IX),INCR,A(1,NC-IX+1),INCR)
        IX = IX + INCC
   10 CONTINUE
   15 RETURN
C                                  CODE FOR BOTH INCREMENTS EQUAL TO 1
C                                    CLEAN-UP LOOP SO REMAINING NUMBER
C                                    OF COLUMNS A MULTIPLE OF 3.
   20 NH = NC/2
      M = NH - (NH/3)*3
      IF (M.EQ.0) GO TO 30
      DO 25 I = 1, M
        CALL SSWAP(NR,A(1,I),INCR,A(1,NC-I+1),INCR)
   25 CONTINUE
      IF (NH.LT.3) RETURN
   30 MP1 = M + 1
      IF (NH.EQ.0) GO TO 40
      DO 35 I = MP1, NH, 3
        CALL SSWAP(NR,A(1,I),INCR,A(1,NC-I+1),INCR)
        CALL SSWAP(NR,A(1,I+1),INCR,A(1,NC-I),INCR)
        CALL SSWAP(NR,A(1,I+2),INCR,A(1,NC-I-1),INCR)
   35 CONTINUE
   40 RETURN
   45 CONTINUE
C                                  CODE FOR EQUAL, POSITIVE, NONUNIT
C                                    INCREMENTS.
      NH = NC/2
      NS = NH*INCC
      IF (NH.EQ.0) GO TO 55
      DO 50 I = 1, NS, INCC
        CALL SSWAP(NR,A(1,I),INCR,A(1,NC-I+1),INCR)
   50 CONTINUE
C                                  SWAP ROWS
   55 INCA = INCC*IA
      IF (INCR.EQ.INCC) IF (INCR-1) 60, 75, 100
C                                  CODE FOR UNEQUAL OR NONPOSITIVE
C                                    INCREMENTS.
   60 IX = 1
      NH = NR/2
      IF (INCR.LT.0) IX = (-NR+1)*INCR + 1
      IF (NH.EQ.0) GO TO 70
      DO 65 I = 1, NH
        CALL SSWAP(NC,A(IX,1),INCA,A(NR-IX+1,1),INCA)
        IX = IX + INCR
   65 CONTINUE
   70 RETURN
C                                  CODE FOR BOTH INCREMENTS EQUAL TO 1
C                                    CLEAN-UP LOOP SO REMAINING NUMBER
C                                    OF COLUMNS A MULTIPLE OF 3.
   75 NH = NR/2
      M = NH - (NH/3)*3
      IF (M.EQ.0) GO TO 85
      DO 80 I = 1, M
        CALL SSWAP(NC,A(I,1),INCA,A(NR-I+1,1),INCA)
   80 CONTINUE
      IF (NH.LT.3) RETURN
   85 MP1 = M + 1
      IF (NH.EQ.0) GO TO 95
      DO 90 I = MP1, NH, 3
        CALL SSWAP(NC,A(I,1),INCA,A(NR-I+1,1),INCA)
        CALL SSWAP(NC,A(I+1,1),INCA,A(NR-I,1),INCA)
        CALL SSWAP(NC,A(I+2,1),INCA,A(NR-I-1,1),INCA)
   90 CONTINUE
   95 RETURN
  100 CONTINUE
C                                  CODE FOR EQUAL, POSITIVE, NONUNIT
C                                    INCREMENTS.
      NH = NR/2
      NS = NH*INCR
      IF (NH.EQ.0) GO TO 110
      DO 105 I = 1, NS, INCR
        CALL SSWAP(NC,A(I,1),INCA,A(NR-I+1,1),INCA)
  105 CONTINUE
  110 RETURN
      END
