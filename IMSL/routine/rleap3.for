C   IMSL ROUTINE NAME   - RLEAP3
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - NUCLEUS CALLED ONLY BY IMSL SUBROUTINE
C                           RLEAP
C
C   REQD. IMSL ROUTINES - NONE REQUIRED
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE RLEAP3 (RSS,CAB,KO,CL,RM,N,NS)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            KO,N,NS
      REAL               RSS,CAB,CL(NS,KO),RM(NS,1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            L
C                                  FIRST EXECUTABLE STATEMENT
      DO 5 L = 1,KO
         IF (CAB .EQ. CL(N,L)) RETURN
    5 CONTINUE
      L = 0
   10 L = L+1
      IF (RSS .GT. RM(N,L+1)) GO TO 15
      RM(N,L) = RM(N,L+1)
      CL(N,L) = CL(N,L+1)
      GO TO 10
   15 RM(N,L) = RSS
      CL(N,L) = CAB
      RETURN
      END
