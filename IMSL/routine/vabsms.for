C   IMSL ROUTINE NAME   - VABSMS
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - SUM OF THE ABSOLUTE VALUES OF THE ELEMENTS OF
C                           A ROW (OR COLUMN) OF A MATRIX STORED IN
C                           SYMMETRIC STORAGE MODE
C
C   USAGE               - CALL VABSMS (V,L,IRO,VSUM)
C
C   ARGUMENTS    V      - NAME OF THE MATRIX (ADDRESS OF THE FIRST
C                           ELEMENT). (INPUT)
C                L      - ORDER OF THE MATRIX. (INPUT)
C                IRO    - ROW OR COLUMN NUMBER (IN RELATION TO THE
C                           MATRIX) OF THE ELEMENTS TO BE SUMMED.
C                           (INPUT)
C                VSUM   - SUM OF THE ABSOLUTE VALUES. (OUTPUT)
C
C   REQD. IMSL ROUTINES - SINGLE/NONE REQUIRED
C                       - DOUBLE/VXADD,VXSTO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE VABSMS (V,L,IRO,VSUM)
C
      DIMENSION          V(1)
      DOUBLE PRECISION   TEMP
C                                  FIRST EXECUTABLE STATEMENT
      IRO1=IRO-1
      IBEG=1+(IRO*IRO1)/2
      IEND=IBEG+IRO1
      TEMP = 0.D0
         DO 5 I=IBEG,IEND
         TEMP=TEMP+ABS(V(I))
    5    CONTINUE
      K=IEND+IRO
      IEND=L-IRO
      IF (IEND .EQ. 0) GO TO 15
         DO 10 I=1,IEND
         TEMP=TEMP+ABS(V(K))
         K=K+IRO+I
   10    CONTINUE
   15 VSUM=TEMP
      RETURN
      END
