C   IMSL ROUTINE NAME   - VBLA=DROTG
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - NOVEMBER 1, 1984
C
C   PURPOSE             - CONSTRUCT GIVENS PLANE ROTATION
C                           (DOUBLE PRECISION)
C
C   USAGE               - CALL DROTG (DA,DB,DC,DS)
C
C   ARGUMENTS    DA     - FIRST ELEMENT OF DOUBLE PRECISION VECTOR.
C                           (INPUT/OUTPUT)
C                         ON OUTPUT, R=(+/-)DSQRT(DA**2 + DB**2)
C                           OVERWRITES DA.
C                DB     - SECOND ELEMENT OF DOUBLE PRECISION VECTOR.
C                           (INPUT/OUTPUT)
C                         ON OUTPUT, Z OVERWRITES DB.
C                           Z IS DEFINED TO BE..
C                           DS,     IF DABS(DA).GT.DABS(DB)
C                           1.0D0/DC, IF DABS(DB).GE.DABS(DA) AND
C                                     DC.NE.0.0D0
C                           1.0D0,    IF DC.EQ.0.0D0.
C                DC     - DOUBLE PRECISION ELEMENT OF OUTPUT
C                           TRANSFORMATION MATRIX. SEE REMARKS.
C                DS     - DOUBLE PRECISION ELEMENT OF OUTPUT
C                           TRANSFORMATION MATRIX. SEE REMARKS.
C
C   REQD. IMSL ROUTINES - NONE REQUIRED
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS      DROTG CONSTRUCTS THE GIVENS TRANSFORMATION
C                    ( DC  DS )
C                G = (        ) ,  DC**2 + DS**2 = 1 ,
C                    (-DS  DC )
C                WHICH ZEROS THE SECOND ELEMENT OF (DA,DB)**T.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE DROTG  (DA,DB,DC,DS)
C
C                                  SPECIFICATIONS FOR ARGUMENTS
      DOUBLE PRECISION   DA,DB,DC,DS
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      DOUBLE PRECISION   U,V,R
C                                  FIRST EXECUTABLE STATEMENT
      IF (DABS(DA).LE.DABS(DB)) GO TO 5
C                                  HERE DABS(DA) .GT. DABS(DB)
      U = DA+DA
      V = DB/U
C                                  NOTE THAT U AND R HAVE THE SIGN OF
C                                    DA
      R = DSQRT(.25D0+V**2)*U
C                                  NOTE THAT DC IS POSITIVE
      DC = DA/R
      DS = V*(DC+DC)
      DB = DS
      DA = R
      RETURN
C                                  HERE DABS(DA) .LE. DABS(DB)
    5 IF (DB.EQ.0.D0) GO TO 15
      U = DB+DB
      V = DA/U
C                                  NOTE THAT U AND R HAVE THE SIGN OF
C                                    DB (R IS IMMEDIATELY STORED IN DA)
      DA = DSQRT(.25D0+V**2)*U
C                                  NOTE THAT DS IS POSITIVE
      DS = DB/DA
      DC = V*(DS+DS)
      IF (DC.EQ.0.D0) GO TO 10
      DB = 1.D0/DC
      RETURN
   10 DB = 1.D0
      RETURN
C                                  HERE DA = DB = 0.D0
   15 DC = 1.D0
      DS = 0.D0
      DA = 0.D0
      DB = 0.D0
      RETURN
C
      END
