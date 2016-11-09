C   IMSL ROUTINE NAME   - ICSSG
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1982
C
C   PURPOSE             - NUCLEUS CALLED ONLY BY IMSL SUBROUTINE ICSSCV
C
C   REQD. IMSL ROUTINES - ICSSF,ICSSH
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE ICSSG  (Z,Y,L,V,N,H,RO,AWK)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N
      REAL               Z(1),Y(1),L(1),V(1),H,RO,AWK(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      REAL               DELTA,ERR,GF1,GF2,GF3,GF4,R1,R2,R3,R4,TAU,
     *                   TOL,RATIO
C                                  TOL = MACHINE PRECISION PARAMETER
      DATA               TOL/Z3C100000/
      DATA               RATIO/2.0/
      DATA               TAU/1.618034/
C                                  CALCULATE THE MINIMUM OF THE
C                                    CROSS VALIDATION FUNCTION.
C                                  FIRST EXECUTABLE STATEMENT
      R1 = 1./(H*H*H)
      R2 = RATIO*R1
      CALL ICSSH (Z,Y,L,V,N,H,R2,GF2,AWK)
    5 CALL ICSSH (Z,Y,L,V,N,H,R1,GF1,AWK)
      IF (GF1.GT.GF2) GO TO 10
      R2 = R1
      GF2 = GF1
      R1 = R1/RATIO
      GO TO 5
   10 R3 = RATIO*R2
   15 CALL ICSSH (Z,Y,L,V,N,H,R3,GF3,AWK)
      IF (GF3.GT.GF2) GO TO 20
      R2 = R3
      GF2 = GF3
      R3 = RATIO*R3
      GO TO 15
   20 R2 = R3
      GF2 = GF3
      DELTA = (R2-R1)/TAU
      R4 = R1+DELTA
      R3 = R2-DELTA
      CALL ICSSH (Z,Y,L,V,N,H,R3,GF3,AWK)
      CALL ICSSH (Z,Y,L,V,N,H,R4,GF4,AWK)
   25 IF (GF3-GF4) 30,30,35
   30 R2 = R4
      GF2 = GF4
      R4 = R3
      GF4 = GF3
      DELTA = DELTA/TAU
      R3 = R2-DELTA
      CALL ICSSH (Z,Y,L,V,N,H,R3,GF3,AWK)
      GO TO 40
   35 R1 = R3
      GF1 = GF3
      R3 = R4
      GF3 = GF4
      DELTA = DELTA/TAU
      R4 = R1+DELTA
      CALL ICSSH (Z,Y,L,V,N,H,R4,GF4,AWK)
   40 ERR = (R2-R1)/(R1+R2)
      IF (ERR.GT.AMAX1(100.0*TOL,1.0E-6)) GO TO 25
      RO = (R1+R2)*.5
      RETURN
      END
