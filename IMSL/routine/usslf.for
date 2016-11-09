C   IMSL ROUTINE NAME   - USSLF
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1981
C
C   PURPOSE             - PRINT A STEM AND LEAF DISPLAY
C
C   USAGE               - CALL USSLF (X,N,IUNIT,MAXL)
C
C   ARGUMENTS    X      - INPUT VECTOR OF LENGTH N CONTAINING THE
C                           DATA TO BE DISPLAYED.
C                         ON OUTPUT, X WILL BE SORTED.
C                N      - THE NUMBER OF OBSERVATIONS IN X. (INPUT)
C                IUNIT  - SIZE OF THE INCREMENTS ON THE STEM. (INPUT)
C                           IF IUNIT IS SET SO SMALL THAT THE LENGTH OF
C                           THE STEM IS MORE THAN 60 LINES, USSLF WILL
C                           USE AN IUNIT SUCH THAT THE STEM WILL BE NO
C                           LONGER THAN 60 LINES. HOWEVER, IF THE USER
C                           SPECIFIES IUNIT AS A NEGATIVE INTEGER,
C                           USSLF WILL USE THE ABSOLUTE VALUE OF IUNIT,
C                           EVEN IF THE STEM WOULD BECOME VERY LONG.
C                           A COMMON VALUE FOR IUNIT IS 10.
C                MAXL   - MAXIMUM DISPLAY WIDTH. (INPUT)
C                           MAXL MUST BE 80 OR 129.
C
C   REQD. IMSL ROUTINES - UGETIO,VSRTA
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE USSLF  (X,N,IUNIT,MAXL)
C
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N,IUNIT,MAXL
      REAL               X(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            ICOUNT,IDASH,II,ITEMP(120),IX,I,J,LABEL,LINES,
     1                   MAXI,MII,MPI,NEG,NIN,NOUT,NUNIT
      DATA               IDASH/1H-/
C                                  FIRST EXECUTABLE STATEMENT
      CALL VSRTA (X,N)
      NEG = 0
      IF (X(1)+.5.LT.0.) NEG = 1
      MII = 0
      MAXI = 120
      IF (MAXL.NE.129) MAXI = 70
      NUNIT=IUNIT
      IF (NUNIT.LT.0) GO TO 10
    5 LINES = (X(N)-X(1))/NUNIT+1
      IF (LINES.LT.60) GO TO 15
      NUNIT = NUNIT*2
      GO TO 5
   10 NUNIT = IABS(NUNIT)
   15 CONTINUE
      ICOUNT = INT(X(1)/NUNIT)
      LABEL = ICOUNT*NUNIT/10
      II = 0
      CALL UGETIO (1,NIN,NOUT)
      DO 50 I=1,N
C                                  GET NEXT DATA POINT
         IX = X(I)
         IF (X(I).LT.0..AND.X(I).LE.IX-.5) IX = IX-1
         IF (X(I).GE.0..AND.X(I).GE.IX+.5) IX = IX+1
         IF (IX.GE.0.AND.NEG.EQ.0) GO TO 25
   20    IF (IX.EQ.0) GO TO 30
         IF (IX.LE.NUNIT*ICOUNT) GO TO 45
         GO TO 30
   25    IF (IX.LT.NUNIT*(ICOUNT+1)) GO TO 45
   30    IF (LABEL.EQ.0.AND.NEG.EQ.1) GO TO 40
         IF (II.EQ.0) WRITE (NOUT,55) LABEL
         IF (II.GT.0.AND.II.LT.MAXI) WRITE (NOUT,60) LABEL,(ITEMP(J),J
     1   =1,II)
         IF (II.EQ.70) WRITE (NOUT,65) LABEL,(ITEMP(J),J=1,II)
         IF (II.EQ.120) WRITE (NOUT,70) LABEL,(ITEMP(J),J=1,II)
         II = 0
   35    ICOUNT = ICOUNT+1
         LABEL = ICOUNT*NUNIT/10
         IF (IX.LT.0) GO TO 20
         GO TO 25
   40    IF (II.EQ.0) WRITE (NOUT,80)
         IF (II.GT.0.AND.II.LT.MAXI) WRITE (NOUT,85) (ITEMP(J),J=1,II)
         IF (II.EQ.70) WRITE (NOUT,90) (ITEMP(J),J=1,II)
         IF (II.EQ.120) WRITE (NOUT,95) (ITEMP(J),J=1,II)
         II = 0
         IF (LABEL.NE.0.OR.ICOUNT.NE.0) GO TO 35
         NEG = 0
         GO TO 25
   45    IF (II.EQ.MAXI) GO TO 50
         II = II+1
         IF (II.GT.MII) MII = II
         IF (IX.LT.0) IX = IABS(IX)
         ITEMP(II) = MOD(IX,10)
   50 CONTINUE
      IF (II.LT.MAXI) WRITE (NOUT,60) LABEL,(ITEMP(J),J=1,II)
      IF (II.EQ.70) WRITE (NOUT,65) LABEL,(ITEMP(J),J=1,II)
      IF (II.EQ.120) WRITE (NOUT,70) LABEL,(ITEMP(J),J=1,II)
C                                  PRINT DASHED COUNTING LINE
      MPI = MII/10+1
      IF (MII.EQ.70.OR.MII.EQ.120) MPI = MII/10
      WRITE (NOUT,75) (IDASH,I=1,MPI)
   55 FORMAT (1X,I5)
   60 FORMAT (1X,I5,1X,120I1)
   65 FORMAT (1X,I5,1X,70I1,1X,1H+)
   70 FORMAT (1X,I5,1X,120I1,1X,1H+)
   75 FORMAT (6X,1H+,12(A1,8(1H-),1H+))
   80 FORMAT (4X,1H-,1H0)
   85 FORMAT (4X,1H-,1H0,1X,120I1)
   90 FORMAT (4X,1H-,1H0,1X,70I1,1X,1H+)
   95 FORMAT (4X,1H-,1H0,1X,120I1,1X,1H+)
      RETURN
      END
