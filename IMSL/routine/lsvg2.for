C   IMSL ROUTINE NAME   - LSVG2
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JULY 1, 1983
C
C   PURPOSE             - NUCLEUS CALLED ONLY BY IMSL ROUTINE LSVDB
C
C   REQD. IMSL ROUTINES - NONE REQUIRED
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE LSVG2  (CS,SN,X,Y)
C
C                                  SPECIFICATIONS FOR ARGUMENTS
      REAL               CS,SN,X,Y
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      REAL               XR
C                                  FIRST EXECUTABLE STATEMENT
      XR=CS*X+SN*Y
      Y=-SN*X+CS*Y
      X=XR
      RETURN
      END
