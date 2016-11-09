C   IMSL ROUTINE NAME   - VNPINV
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1982
C
C   PURPOSE             - NUCLEUS CALLED BY IMSL ROUTINES WHICH SORT
C                           A MATRIX WITH RESPECT TO VECTORS.
C
C   REQD. IMSL ROUTINES - VNINI,VNSCPY,VSRTP
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE VNPINV (IR,LA,WK)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            LA,IR(1)
      REAL               WK(1)
C                                  FIRST EXECUTABLE STATEMENT
      IF (LA.LE.0) RETURN
      CALL VNSCPY(LA,IR,1,WK,1)
      CALL VNINI(LA,IR,1,1,1)
      CALL VSRTP(WK,LA,IR)
      RETURN
      END
