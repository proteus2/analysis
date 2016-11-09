C   IMSL ROUTINE NAME   - VSCMPA
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1982
C
C   PURPOSE             - NUCLEUS CALLED BY IMSL ROUTINES WHICH SORT
C                           A MATRIX WITH RESPECT TO VECTORS.
C
C   REQD. IMSL ROUTINES - NONE REQUIRED
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE VSCMPA (N,SX,INCX,SY,INCY,KOMPAR)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N,INCX,INCY,KOMPAR
      REAL               SX(1),SY(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,IX,IY,M,MP1
C                                  FIRST EXECUTABLE STATEMENT
      KOMPAR = -999
      IF (N.LE.0) GO TO 165
      KOMPAR = 0
      IF (INCX.EQ.1 .AND. INCY.EQ.1) GO TO 25
C                                  CODE FOR UNEQUAL INCREMENTS OR
C                                    INCREMENTS NOT EQUAL TO 1.
      IX = 1
      IY = 1
      IF (INCX.LT.0) IX = (-N+1)*INCX+1
      IF (INCY.LT.0) IY = (-N+1)*INCY+1
      DO 20 I=1,N
         IF (ABS(SX(IX))-ABS(SY(IY))) 5, 15, 10
    5    KOMPAR = -1
         GO TO 165
   10    KOMPAR = 1
         GO TO 165
   15    IX = IX+INCX
         IY = IY+INCY
   20 CONTINUE
      GO TO 165
C                                  CODE FOR BOTH INCREMENTS EQUAL TO 1.
C                                    CLEAN UP LOOP.
   25 M = MOD(N,7)
      IF (M.EQ.0) GO TO 50
      DO 40 I=1,M
         IF (ABS(SX(I))-ABS(SY(I))) 30, 40, 35
   30    KOMPAR = -1
         GO TO 45
   35    KOMPAR = 1
         GO TO 45
   40 CONTINUE
   45 IF (N.LT.7) GO TO 165
   50 MP1 = M+1
      DO 160 I=MP1,N,7
         IF (ABS(SX(I))-ABS(SY(I))) 55, 65, 60
   55    KOMPAR = -1
         GO TO 165
   60    KOMPAR = 1
         GO TO 165
   65    IF (ABS(SX(I+1))-ABS(SY(I+1))) 70, 80, 75
   70    KOMPAR = -1
         GO TO 165
   75    KOMPAR = 1
         GO TO 165
   80    IF (ABS(SX(I+2))-ABS(SY(I+2))) 85, 95, 90
   85    KOMPAR = -1
         GO TO 165
   90    KOMPAR = 1
         GO TO 165
   95    IF (ABS(SX(I+3))-ABS(SY(I+3))) 100, 110, 105
  100    KOMPAR = -1
         GO TO 165
  105    KOMPAR = 1
         GO TO 165
  110    IF (ABS(SX(I+4))-ABS(SY(I+4))) 115, 125, 120
  115    KOMPAR = -1
         GO TO 165
  120    KOMPAR = 1
         GO TO 165
  125    IF (ABS(SX(I+5))-ABS(SY(I+5))) 130, 140, 135
  130    KOMPAR = -1
         GO TO 165
  135    KOMPAR = 1
         GO TO 165
  140    IF (ABS(SX(I+6))-ABS(SY(I+6))) 145, 155, 150
  145    KOMPAR = -1
         GO TO 165
  150    KOMPAR = 1
         GO TO 165
  155    CONTINUE
  160 CONTINUE
  165 RETURN
      END
