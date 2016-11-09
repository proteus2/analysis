C   IMSL ROUTINE NAME   - VDCMPE
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1982
C
C   PURPOSE             - NUCLEUS CALLED BY IMSL ROUTINES WHICH SORT
C                           A MATRIX WITH RESPECT TO VECTORS.
C
C   REQD. IMSL ROUTINES - VBLA=IDAMAX,VBLA=DASUM,VBLA=DNRM2
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE VDCMPE (IOPT,N,SX,INCX,SY,INCY,KOMPAR)
C                                  SPECIFICATIONS FOR FUNCTIONS
C                                  SPECIFICATIONS FOR FUNCTIONS
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            IOPT,N,INCX,INCY,KOMPAR
      DOUBLE PRECISION   SX(1),SY(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,I1,I2,IX,IY,M,MP1
C                                  SPECIFICATIONS FOR FUNCTIONS
      DOUBLE PRECISION   DASUM
      DOUBLE PRECISION   DNRM2
      INTEGER            IDAMAX
C                                  FIRST EXECUTABLE STATEMENT
      KOMPAR = -999
      IF (N.LE.0) GO TO 45
      KOMPAR = 0
      IF (IOPT.EQ.2) GO TO 15
      IF (IOPT.EQ.3) GO TO 30
C                                  CODE FOR L-2 NORM COMPARISON
      IF (DNRM2(N,SX,INCX)-DNRM2(N,SY,INCY)) 5, 45, 10
    5 KOMPAR = -1
      GO TO 45
   10 KOMPAR = 1
      GO TO 45
C                                  CODE FOR L-1 NORM COMPARISON
   15 IF (DASUM(N,SX,INCX)-DASUM(N,SY,INCY)) 20, 45, 25
   20 KOMPAR = -1
      GO TO 45
   25 KOMPAR = 1
      GO TO 45
C                                  CODE FOR L-INFINITY NORM COMPARISON
   30 I1 = 1+(IDAMAX(N,SX,INCX)-1)*INCX
      I2 = 1+(IDAMAX(N,SY,INCY)-1)*INCY
      IF (DABS(SX(I1))-DABS(SY(I2))) 35, 45, 40
   35 KOMPAR = -1
      GO TO 45
   40 KOMPAR = 1
   45 RETURN
      END
