C   IMSL ROUTINE NAME   - VSARX
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1982
C
C   PURPOSE             - NUCLEUS CALLED BY IMSL ROUTINES WHICH SORT
C                           A MATRIX WITH RESPECT TO VECTORS.
C
C   REQD. IMSL ROUTINES - VBLA=SCOPY,VNABSX
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE VSARX (A,IA,NR,NC,IOP,IR,WK)
C
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            IA,NR,NC,IOP,IR(1)
      REAL               A(IA,1),WK(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            IPTR,K,L
C                                  FIRST EXECUTABLE STATEMENT
      IPTR = 1
      IF (IOP.EQ.1) GO TO 30
C                                  SORT A BY COLUMNS (IOP = 0)
C                                  CHECK IF ALL COLUMNS ARE SORTED
    5 IF (IPTR.GE.NC) GO TO 55
C                                  CHECK IF COLUMN IPTR HAS BEEN SORTED
      IF (IR(IPTR).GT.0) GO TO 15
   10 IPTR = IPTR+1
      GO TO 5
C                                  CHECK IF COLUMN IPTR NEED BE MOVED
   15 IF (IR(IPTR).EQ.IPTR) GO TO 10
      K = IPTR
C                                  STORE COLUMN IPTR IN TEMPORARY VECTOR
      CALL SCOPY(NR,A(1,K),1,WK,1)
   20 L = IR(K)
C                                  CHECK IF TEMPORARY VECTOR NEEDED HERE
      IF (L.EQ.IPTR) GO TO 25
C                                  INSERT COLUMN L INTO COLUMN K
      CALL SCOPY(NR,A(1,L),1,A(1,K),1)
C                                  MARK COLUMN K AS ALREADY SORTED
      IR(K) = -IR(K)
      K = L
      GO TO 20
C                                  INSERT TEMPORARY VECTOR IN COLUMN K
   25 CALL SCOPY(NR,WK,1,A(1,K),1)
      IR(K) = -IR(K)
      GO TO 10
C                                  SORT A BY ROWS (IOP = 1)
   30 IF (IPTR.GE.NR) GO TO 55
      IF (IR(IPTR).GT.0) GO TO 40
   35 IPTR = IPTR+1
      GO TO 30
   40 IF (IR(IPTR).EQ.IPTR) GO TO 35
      K = IPTR
      CALL SCOPY(NC,A(K,1),IA,WK,1)
   45 L = IR(K)
      IF (L.EQ.IPTR) GO TO 50
      CALL SCOPY(NC,A(L,1),IA,A(K,1),IA)
      IR(K) = -IR(K)
      K = L
      GO TO 45
   50 CALL SCOPY(NC,WK,1,A(K,1),IA)
      IR(K) = -IR(K)
      GO TO 35
   55 IF (IOP.EQ.0) CALL VNABSX(NC,IR,1)
      IF (IOP.EQ.1) CALL VNABSX(NR,IR,1)
      RETURN
      END
