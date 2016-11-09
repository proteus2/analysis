C   IMSL ROUTINE NAME   - VABMXS
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - MAXIMUM ABSOLUTE VALUE OF THE ELEMENTS OF A
C                           ROW OR COLUMN OF A MATRIX STORED IN
C                           SYMMETRIC STORAGE MODE
C
C   USAGE               - CALL VABMXS (V,L,IRO,J,VMAX)
C
C   ARGUMENTS    V      - NAME OF THE MATRIX (ADDRESS OF THE FIRST
C                           ELEMENT). (INPUT)
C                L      - ORDER OF THE MATRIX. (INPUT)
C                IRO    - ROW OR COLUMN NUMBER (IN RELATION TO THE
C                           MATRIX AS IF IN FULL STORAGE MODE) OF THE
C                           ELEMENTS TO BE TESTED. (INPUT)
C                J      - THE ELEMENT FOUND TO HAVE THE MAXIMUM
C                           ABSOLUTE VALUE WOULD BE V(IRO,J) = V(J,IRO)
C                           IF THE MATRIX WERE STORED IN FULL STORAGE
C                           MODE. (OUTPUT)
C                VMAX   - MAXIMUM ABSOLUTE VALUE OF THE ELEMENTS OF THE
C                           ROW OR COLUMN. (OUTPUT)
C
C   REQD. IMSL ROUTINES - NONE REQUIRED
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE VABMXS (V,L,IRO,J,VMAX)
C
      REAL               V(1),T,VMAX,ZERO
      DATA               ZERO/0.0/
C                                  FIRST EXECUTABLE STATEMENT
      VMAX = ZERO
      IRO1=IRO-1
      IBEG=1+(IRO*IRO1)/2
      J=IBEG
      IEND=IBEG+IRO1
         DO 5 I=IBEG,IEND
         T=ABS(V(I))
         IF (T .LE. VMAX)  GO TO 5
         VMAX=T
         J=I
    5    CONTINUE
      J=J-IBEG+1
      K=IEND+IRO
      IEND=L-IRO
      IF (IEND .EQ. 0) GO TO 15
         DO 10 I=1,IEND
         T=ABS(V(K))
         IF (T .LE. VMAX) GO TO 10
         VMAX=T
         J=IRO+I
         K=K+IRO+I
   10    CONTINUE
   15 RETURN
      END
