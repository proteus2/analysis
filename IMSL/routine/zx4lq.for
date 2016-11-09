C   IMSL ROUTINE NAME   - ZX4LQ
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1980
C
C   PURPOSE             - NUCLEUS CALLED ONLY BY IMSL SUBROUTINE ZX4LP
C
C   REQD. IMSL ROUTINES - SINGLE/VBLA=SROTM,VBLA=SROTMG
C                       - DOUBLE/VBLA=DROTM,VBLA=DROTMG
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE ZX4LQ (B,M,J,MDB,D,IFLAG)
C
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            M,J,MDB,IFLAG
      REAL               B(MDB,MDB),D(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            MP2
      REAL               SPARAM(5),ZERO
C                                  FIRST EXECUTABLE STATEMENT
C                                  ELIMINATE ELEMENTS 1,...,J OF COL.
C                                    M+2 TO FORM A LOWER TRAPAZOID OF
C                                    DIMENSION M+2 BY J. A POSITIVE
C                                    VALUE OF IFLAG INDICATES THAT A
C                                    ZERO SCALE FACTOR, D(IFLAG)=0.,
C                                    WAS SEEN.
C
      IFLAG = 0
      MP2 = M+2
      ZERO = 0.
      IF (J.LE.0.OR.J.GE.MP2) RETURN
      DO 5 I=1,J
         CALL SROTMG (D(I),D(MP2),B(I,I),B(I,MP2),SPARAM)
         IF (D(I).EQ.ZERO) IFLAG = I
         CALL SROTM (MP2-I,B(I+1,I),1,B(I+1,MP2),1,SPARAM)
    5 CONTINUE
      RETURN
      END
