C   IMSL ROUTINE NAME   - DBLNB
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1982
C
C   PURPOSE             - NUCLEUS CALLED ONLY BY IMSL SUBROUTINE DBLIN
C
C   REQD. IMSL ROUTINES - DBLNC,DBLND,UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      REAL FUNCTION DBLNB (F,X,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            IER
      REAL               F,X
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            JER
      REAL               X1,AY1,BY1,AERR1,DBLNC,ERR,ZERO
      EXTERNAL           F
      COMMON /DBLIC/     X1,AY1,BY1,AERR1
      DATA               ZERO/0.0/
C                                  FIRST EXECUTABLE STATEMENT
      X1 = X
      DBLNB = DBLNC(F,AY1,BY1,AERR1,ZERO,ERR,JER)
      IER = MAX0(IER,JER)
      RETURN
      END
