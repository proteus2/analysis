C   IMSL ROUTINE NAME   - VSCMPE
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1982
C
C   PURPOSE             - NUCLEUS CALLED BY IMSL ROUTINES WHICH SORT
C                           A MATRIX WITH RESPECT TO VECTORS.
C
C   REQD. IMSL ROUTINES - VBLA=ISAMAX,VBLA=SASUM,VBLA=SNRM2
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE VSCMPE (IOPT,N,SX,INCX,SY,INCY,KOMPAR)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            IOPT,N,INCX,INCY,KOMPAR
      REAL               SX(1),SY(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,I1,I2,IX,IY,M,MP1
C                                  SPECIFICATIONS FOR FUNCTIONS
      INTEGER            ISAMAX
      REAL               SASUM
      REAL               SNRM2
C                                  FIRST EXECUTABLE STATEMENT
      KOMPAR = -999
      IF (N.LE.0) GO TO 50
      KOMPAR = 0
      IF (IOPT.EQ.2) GO TO 20
      IF (IOPT.EQ.3) GO TO 35
C                                  CODE FOR L-2 NORM COMPARISON
    5 IF (SNRM2(N,SX,INCX)-SNRM2(N,SY,INCY)) 10, 50, 15
   10 KOMPAR = -1
      GO TO 50
   15 KOMPAR = 1
      GO TO 50
C                                  CODE FOR L-1 NORM COMPARISON
   20 IF (SASUM(N,SX,INCX)-SASUM(N,SY,INCY)) 25, 50, 30
   25 KOMPAR = -1
      GO TO 50
   30 KOMPAR = 1
      GO TO 50
C                                  CODE FOR L-INFINITY NORM COMPARISON
   35 I1 = 1+(ISAMAX(N,SX,INCX)-1)*INCX
      I2 = 1+(ISAMAX(N,SY,INCY)-1)*INCY
      IF (ABS(SX(I1))-ABS(SY(I2))) 40, 50, 45
   40 KOMPAR = -1
      GO TO 50
   45 KOMPAR = 1
   50 RETURN
      END
