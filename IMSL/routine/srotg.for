C   IMSL ROUTINE NAME   - VBLA=SROTG
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - NOVEMBER 1, 1984
C
C   PURPOSE             - CONSTRUCT GIVENS PLANE ROTATION
C                           (SINGLE PRECISION)
C
C   USAGE               - CALL SROTG (SA,SB,SC,SS)
C
C   ARGUMENTS    SA     - FIRST ELEMENT OF VECTOR. (INPUT/OUTPUT)
C                         ON OUTPUT, R=(+/-)SQRT(SA**2 + SB**2)
C                           OVERWRITES SA.
C                SB     - SECOND ELEMENT OF VECTOR. (INPUT/OUTPUT)
C                         ON OUTPUT, Z OVERWRITES SB.
C                           Z IS DEFINED TO BE..
C                           SS,     IF ABS(SA).GT.ABS(SB)
C                           1.0/SC, IF ABS(SB).GE.ABS(SA) AND SC.NE.0.0
C                           1.,     IF SC.EQ.0.0.
C                SC     - ELEMENT OF OUTPUT TRANSFORMATION MATRIX.
C                           SEE REMARKS.
C                SS     - ELEMENT OF OUTPUT TRANSFORMATION MATRIX.
C                           SEE REMARKS.
C
C   REQD. IMSL ROUTINES - NONE REQUIRED
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS      SROTG CONSTRUCTS THE GIVENS TRANSFORMATION
C                    ( SC  SS )
C                G = (        ) ,  SC**2 + SS**2 = 1 ,
C                    (-SS  SC )
C                WHICH ZEROS THE SECOND ELEMENT OF (SA,SB)**T.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE SROTG  (SA,SB,SC,SS)
C
C                                  SPECIFICATIONS FOR ARGUMENTS
      REAL               SA,SB,SC,SS
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      REAL               R,U,V
C                                  FIRST EXECUTABLE STATEMENT
      IF (ABS(SA).LE.ABS(SB)) GO TO 5
C                                  HERE ABS(SA) .GT. ABS(SB)
      U = SA+SA
      V = SB/U
C                                  NOTE THAT U AND R HAVE THE SIGN OF
C                                    SA
      R = SQRT(.25+V**2)*U
C                                  NOTE THAT SC IS POSITIVE
      SC = SA/R
      SS = V*(SC+SC)
      SB = SS
      SA = R
      RETURN
C                                  HERE ABS(SA) .LE. ABS(SB)
    5 IF (SB.EQ.0.) GO TO 15
      U = SB+SB
      V = SA/U
C                                  NOTE THAT U AND R HAVE THE SIGN OF
C                                    SB (R IS IMMEDIATELY STORED IN SA)
      SA = SQRT(.25+V**2)*U
C                                  NOTE THAT SS IS POSITIVE
      SS = SB/SA
      SC = V*(SS+SS)
      IF (SC.EQ.0.) GO TO 10
      SB = 1./SC
      RETURN
   10 SB = 1.
      RETURN
C                                  HERE SA = SB = 0.
   15 SC = 1.
      SS = 0.
      SA = 0.
      SB = 0.
      RETURN
C
      END
