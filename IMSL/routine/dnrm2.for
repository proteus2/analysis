C   IMSL ROUTINE NAME   - VBLA=DNRM2
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - COMPUTE THE EUCLIDEAN LENGTH OR L2 NORM
C                           OF A DOUBLE PRECISION VECTOR
C
C   USAGE               - FUNCTION DNRM2 (N,DX,INCX)
C
C   ARGUMENTS    DNRM2  - DOUBLE PRECISION SQUARE ROOT OF THE SUM FROM
C                           I=1 TO N OF X(I)**2. (OUTPUT)
C                           X(I) REFERS TO A SPECIFIC ELEMENT OF DX.
C                           SEE INCX ARGUMENT DESCRIPTION.
C                N      - LENGTH OF VECTOR X. (INPUT)
C                DX     - DOUBLE PRECISION VECTOR OF LENGTH N*INCX.
C                           (INPUT)
C                INCX   - DISPLACEMENT BETWEEN ELEMENTS OF DX. (INPUT)
C                           X(I) IS DEFINED TO BE DX(1+(I-1)*INCX).
C                           INCX MUST BE GREATER THAN ZERO.
C
C   REQD. IMSL ROUTINES - NONE REQUIRED
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      DOUBLE PRECISION FUNCTION DNRM2 (N,DX,INCX)
C
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N,INCX
      DOUBLE PRECISION   DX(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,J,NEXT,NN
      DOUBLE PRECISION   CUTLO,CUTHI,SUM,XMAX,ZERO,ONE,HITEST
      DATA               ZERO, ONE /0.0D0, 1.0D0/
      DATA               CUTLO, CUTHI / 8.232D-11,  1.304D19 /
C                                  FIRST EXECUTABLE STATEMENT
      IF (N.GT.0) GO TO 5
      DNRM2 = ZERO
      GO TO 70
C
    5 ASSIGN 15 TO NEXT
      SUM = ZERO
      NN = N*INCX
C                                  BEGIN MAIN LOOP
      I = 1
   10 GO TO NEXT, (15,20,35,40)
   15 IF (DABS(DX(I)).GT.CUTLO) GO TO 55
      ASSIGN 20 TO NEXT
      XMAX = ZERO
C                                  PHASE 1. SUM IS ZERO
   20 IF (DX(I).EQ.ZERO) GO TO 65
      IF (DABS(DX(I)).GT.CUTLO) GO TO 55
C                                  PREPARE FOR PHASE 2.
      ASSIGN 35 TO NEXT
      GO TO 30
C                                  PREPARE FOR PHASE 4.
   25 I = J
      ASSIGN 40 TO NEXT
      SUM = (SUM/DX(I))/DX(I)
   30 XMAX = DABS(DX(I))
      GO TO 45
C                                  PHASE 2. SUM IS SMALL. SCALE TO
C                                    AVOID DESTRUCTIVE UNDERFLOW.
   35 IF (DABS(DX(I)).GT.CUTLO) GO TO 50
C                                  COMMON CODE FOR PHASES 2 AND 4. IN
C                                    PHASE 4 SUM IS LARGE. SCALE TO
C                                    AVOID OVERFLOW.
   40 IF (DABS(DX(I)).LE.XMAX) GO TO 45
      SUM = ONE+SUM*(XMAX/DX(I))**2
      XMAX = DABS(DX(I))
      GO TO 65
C
   45 SUM = SUM+(DX(I)/XMAX)**2
      GO TO 65
C                                  PREPARE FOR PHASE 3.
   50 SUM = (SUM*XMAX)*XMAX
C                                  FOR REAL OR D.P. SET HITEST =
C                                    CUTHI/N FOR COMPLEX SET HITEST =
C                                    CUTHI/(2*N)
   55 HITEST = CUTHI/FLOAT(N)
C                                  PHASE 3. SUM IS MID-RANGE. NO
C                                    SCALING.
      DO 60 J=I,NN,INCX
         IF (DABS(DX(J)).GE.HITEST) GO TO 25
   60 SUM = SUM+DX(J)**2
      DNRM2 = DSQRT(SUM)
      GO TO 70
C
   65 CONTINUE
      I = I+INCX
      IF (I.LE.NN) GO TO 10
C                                  END OF MAIN LOOP. COMPUTE SQUARE
C                                    ROOT AND ADJUST FOR SCALING.
      DNRM2 = XMAX*DSQRT(SUM)
   70 CONTINUE
      RETURN
      END
