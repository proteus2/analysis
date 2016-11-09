C   IMSL ROUTINE NAME   - DBLND
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1982
C
C   PURPOSE             - NUCLEUS CALLED ONLY BY IMSL SUBROUTINE DBLIN
C
C   REQD. IMSL ROUTINES - NONE REQUIRED
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      REAL FUNCTION DBLND (F,Y)
C                                  SPECIFICATIONS FOR ARGUMENTS
      REAL               F,Y
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      REAL               XTEMP,YTEMP
      REAL               X1,AY1,BY1,AERR1
      COMMON /DBLIC/     X1,AY1,BY1,AERR1
C                                  FIRST EXECUTABLE STATEMENT
      XTEMP = X1
      YTEMP = Y
      DBLND = F(XTEMP,YTEMP)
      RETURN
      END
