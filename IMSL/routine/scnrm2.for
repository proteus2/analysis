C   IMSL ROUTINE NAME   - VBLA=SCNRM2
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - COMPUTE THE EUCLIDEAN LENGTH OR L2 NORM
C                           OF A COMPLEX VECTOR
C
C   USAGE               - FUNCTION SCNRM2 (N,CX,INCX)
C
C   ARGUMENTS    SCNRM2 - SQUARE ROOT OF THE SUM FROM I=1 TO N OF
C                           X(I)**2. (OUTPUT)
C                           X(I) REFERS TO A SPECIFIC ELEMENT OF CX.
C                           SEE INCX ARGUMENT DESCRIPTION.
C                N      - LENGTH OF VECTOR X. (INPUT)
C                CX     - COMPLEX VECTOR OF LENGTH N*INCX. (INPUT)
C                INCX   - DISPLACEMENT BETWEEN ELEMENTS OF CX. (INPUT)
C                           X(I) IS DEFINED TO BE CX(1+(I-1)*INCX).
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
      REAL FUNCTION SCNRM2 (N,CX,INCX)
C
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N,INCX
      COMPLEX            CX(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      LOGICAL            IMAG,SCALE
      INTEGER            I,NEXT,NN
      REAL               CUTLO,CUTHI,HITEST,SUM,XMAX,ABSX,ZERO,ONE
      DATA               ZERO, ONE /0.0E0, 1.0E0/
      DATA               CUTLO, CUTHI / 4.441E-16,  1.304E19 /
C                                  FIRST EXECUTABLE STATEMENT
      IF (N.GT.0) GO TO 5
      SCNRM2 = ZERO
      GO TO 70
C
    5 ASSIGN 10 TO NEXT
      SUM = ZERO
      NN = N*INCX
C                                  BEGIN MAIN LOOP
      DO 65 I=1,NN,INCX
         ABSX = ABS(REAL(CX(I)))
         IMAG = .FALSE.
         GO TO NEXT, (10,15,30,55,35)
   10    IF (ABSX.GT.CUTLO) GO TO 50
         ASSIGN 15 TO NEXT
         SCALE = .FALSE.
C                                  PHASE 1. SUM IS ZERO
   15    IF (ABSX.EQ.ZERO) GO TO 60
         IF (ABSX.GT.CUTLO) GO TO 50
C                                  PREPARE FOR PHASE 2.
         ASSIGN 30 TO NEXT
         GO TO 25
C                                  PREPARE FOR PHASE 4.
   20    ASSIGN 35 TO NEXT
         SUM = (SUM/ABSX)/ABSX
   25    SCALE = .TRUE.
         XMAX = ABSX
         GO TO 40
C                                  PHASE 2. SUM IS SMALL. SCALE TO
C                                    AVOID DESTRUCTIVE UNDERFLOW.
   30    IF (ABSX.GT.CUTLO) GO TO 45
C                                  COMMON CODE FOR PHASES 2 AND 4. IN
C                                    PHASE 4 SUM IS LARGE. SCALE TO
C                                    AVOID OVERFLOW.
   35    IF (ABSX.LE.XMAX) GO TO 40
         SUM = ONE+SUM*(XMAX/ABSX)**2
         XMAX = ABSX
         GO TO 60
C
   40    SUM = SUM+(ABSX/XMAX)**2
         GO TO 60
C                                  PREPARE FOR PHASE 3.
   45    SUM = (SUM*XMAX)*XMAX
C
   50    ASSIGN 55 TO NEXT
         SCALE = .FALSE.
C                                  FOR REAL OR D.P. SET HITEST =
C                                    CUTHI/N FOR COMPLEX SET HITEST =
C                                    CUTHI/(2*N)
         HITEST = CUTHI/FLOAT(N)
C                                  PHASE 3. SUM IS MID-RANGE. NO
C                                    SCALING.
   55    IF (ABSX.GE.HITEST) GO TO 20
         SUM = SUM+ABSX**2
   60    CONTINUE
C                                  CONTROL SELECTION OF REAL AND
C                                    IMAGINARY PARTS.
         IF (IMAG) GO TO 65
         ABSX = ABS(AIMAG(CX(I)))
         IMAG = .TRUE.
         GO TO NEXT, (15,30,55,35)
C
   65 CONTINUE
C                                  END OF MAIN LOOP. COMPUTE SQUARE
C                                    ROOT AND ADJUST FOR SCALING.
      SCNRM2 = SQRT(SUM)
      IF (SCALE) SCNRM2 = SCNRM2*XMAX
   70 CONTINUE
      RETURN
      END
