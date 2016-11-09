C   IMSL ROUTINE NAME   - USWFM
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1982
C
C   PURPOSE             - PRINT A MATRIX STORED IN FULL STORAGE MODE.
C
C   USAGE               - CALL USWFM (ITITLE,NC,A,IA,N,M,IOPT)
C
C   ARGUMENTS    ITITLE - A CHARACTER STRING TO PROVIDE A TITLE. (INPUT)
C                           THE LENGTH OF ITITLE MUST NOT EXCEED 20.
C                NC     - LENGTH OF ITITLE.  (INPUT)
C                           IF NC.EQ.0 THE TITLE  IS NOT PRINTED.
C                A      - NAME OF THE MATRIX TO BE PRINTED. (INPUT)
C                IA     - ROW DIMENSION OF MATRIX A EXACTLY AS SPECIFIED
C                           IN THE DIMENSION STATEMENT IN THE CALLING
C                           PROGRAM. (INPUT)
C                N      - NUMBER OF ROWS OF MATRIX A TO BE PRINTED.
C                           (INPUT)
C                M      - NUMBER OF COLUMNS OF MATRIX A TO BE PRINTED.
C                           (INPUT)
C                IOPT   - OPTION INDICATING THE FORMAT STATEMENT TO BE
C                           USED. (INPUT)
C
C                           OPTIONS FOR    OPTIONS FOR
C                           129 COLUMNS    80 COLUMNS     FORMAT
C                           -----------    -----------    ------
C
C                                1              2         F18.5
C                                3              4         E15.6
C                                5              6         E25.ISIG
C                         NOTE - ISIG IS INTENDED TO GIVE NEARLY FULL
C                         PRECISION REPRESENTATION OF THE MATRIX
C                         ELEMENTS. SEE REMARKS SECTION.
C                         IF IOPT IS NOT IN THE RANGE 1 TO 6, 6 IS USED.
C
C   REQD. IMSL ROUTINES - UGETIO,USPKD
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS  1.  IF AN ASTERISK APPEARS IN THE HIGH ORDER POSITION OF
C                OF THE FIELD, THEN THE NUMBER TO BE PRINTED EXCEEDS
C                THE MAGNITUDE ALLOWED BY THE SELECTED FORMAT.
C            2.  ISIG IS DEFINED TO BE THE NUMBER OF DIGITS IN DECIMAL
C                CONSTANTS. MORE SPECIFICALLY,
C
C                  ISIG     PRECISION     HARDWARE
C                  ----     ---------     --------
C                    7        SINGLE         H32
C                   16        DOUBLE         H32
C                    9        SINGLE         H36
C                   11        SINGLE         H48
C                   14        SINGLE         H60
C
C            3.  OUTPUT IS WRITTEN TO THE UNIT SPECIFIED BY IMSL
C                ROUTINE UGETIO. SEE THE UGETIO DOCUMENT FOR DETAILS.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE USWFM (ITITLE,NC,A,IA,N,M,IOPT)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            NC,IA,N,M,IOPT
      REAL               A(IA,M)
      INTEGER            ITITLE(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,IOPT0,IOPT1,IPICK(6),J,JTITLE(20),K,L,LBEG,
     *                   LEND,MM,MTIME,NCA,NCAMTB,NCMAX,NIN,NOUT
      DATA               IPICK /6,3,8,5,5,3/
      DATA               NCMAX /20/
C                                  WRITE MATRIX TITLE
C                                  SET IOPT0. IF IOPT IS OUT OF RANGE
C                                  USE 6
C                                  FIRST EXECUTABLE STATEMENT
      IOPT0 = 6
      IF (IOPT.GE.1 .AND. IOPT.LE.6) IOPT0 = IOPT
C                                  WRITE TITLE
      NCA = MIN0(NC,NCMAX)
      CALL USPKD(ITITLE,NCA,JTITLE,NCAMTB)
      CALL UGETIO(1,NIN,NOUT)
      IF (NCA.GT.0) WRITE (NOUT,145) (JTITLE(I),I=1,NCA)
C                                  WRITE COLUMN NUMBERS
      GO TO (5, 10, 15, 20, 25, 30), IOPT0
    5 WRITE (NOUT,85) (I,I=1,M)
      GO TO 35
   10 WRITE (NOUT,90) (I,I=1,M)
      GO TO 35
   15 WRITE (NOUT,95) (I,I=1,M)
      GO TO 35
   20 WRITE (NOUT,100) (I,I=1,M)
      GO TO 35
   25 WRITE (NOUT,135) (I,I=1,M)
      GO TO 35
   30 WRITE (NOUT,140) (I,I=1,M)
   35 IOPT1 = IPICK(IOPT0)
      MTIME = M/IOPT1
      IF (MTIME*IOPT1.EQ.M) MTIME = MTIME-1
C                                  WRITE THE MATRIX
      DO 80 I=1,N
         LEND = IOPT1
         L = MIN0(IOPT1,M)
C                                  WRITE THE FIRST LINE OF THE ROW
         GO TO (40, 40, 45, 45, 50, 50), IOPT0
   40    WRITE (NOUT,105) I, (A(I,J),J=1,L)
         GO TO 55
   45    WRITE (NOUT,110) I, (A(I,J),J=1,L)
         GO TO 55
   50    WRITE (NOUT,120) I, (A(I,J),J=1,L)
   55    MM = M-L
C                                  WRITE THE REMAINING LINES
         IF (MTIME.EQ.0) GO TO 80
         DO 75 K=1,MTIME
            LBEG = LEND+1
            L = MIN0(IOPT1,MM)
            MM = MM-L
            LEND = LEND+L
            GO TO (60, 60, 65, 65, 70, 70), IOPT0
   60       WRITE (NOUT,115) (A(I,J),J=LBEG,LEND)
            GO TO 75
   65       WRITE (NOUT,130) (A(I,J),J=LBEG,LEND)
            GO TO 75
   70       WRITE (NOUT,125) (A(I,J),J=LBEG,LEND)
   75    CONTINUE
   80 CONTINUE
   85 FORMAT (14X, I4, 15X, I4, 15X, I4, 15X, I4, 15X, I4, 15X, I4)
   90 FORMAT (14X, I4, 15X, I4, 15X, I4)
   95 FORMAT (10X, I4, 11X, I4, 11X, I4, 11X, I4, 11X, I4, 11X, I4,
     *11X, I4, 11X, I4)
  100 FORMAT (10X, I4, 11X, I4, 11X, I4, 11X, I4, 11X, I4)
  105 FORMAT (/I4, 6(1X, F18.5))
  110 FORMAT (/I4, 8E15.6)
  115 FORMAT (5X, F18.5, 1X, F18.5, 1X, F18.5, 1X, F18.5, 1X, F18.5,
     *1X, F18.5)
  120 FORMAT (/I4, 5E25.7)
  125 FORMAT (4X, 5E25.7)
  130 FORMAT (4X, 8E15.6)
  135 FORMAT (14X, I4, 21X, I4, 21X, I4, 21X, I4, 21X, I4)
  140 FORMAT (14X, I4, 21X, I4, 21X, I4)
  145 FORMAT (1X, 20A1)
      RETURN
      END
