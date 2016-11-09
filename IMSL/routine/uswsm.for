C   IMSL ROUTINE NAME   - USWSM
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1982
C
C   PURPOSE             - PRINT A MATRIX STORED IN SYMMETRIC STORAGE
C                           MODE.
C
C   USAGE               - CALL USWSM (ITITLE,NC,A,M,IOPT)
C
C   ARGUMENTS    ITITLE - A CHARACTER STRING TO PROVIDE A TITLE. (INPUT)
C                           THE LENGTH OF ITITLE MUST NOT EXCEED 20.
C                NC     - LENGTH OF ITITLE.  (INPUT)
C                           IF NC.EQ.0 THE TITLE  IS NOT PRINTED.
C                A      - NAME OF THE MATRIX TO BE PRINTED. (INPUT)
C                M      - NUMBER OF ROWS (OR COLUMNS) OF MATRIX A TO
C                           BE PRINTED. (INPUT)
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
      SUBROUTINE USWSM (ITITLE,NC,A,M,IOPT)
C
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            NC,M,IOPT
      REAL               A(1)
      INTEGER            ITITLE(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,IOPT0,IOPT1,IPICK(6),J,JTITLE(20),LADD,LEND,
     *                   LI,LMAX,LSTART,LSTOP,NCA,NCAMTB,NCMAX,NIN,
     *                   NOUT
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
      IF (NCA.GT.0) WRITE (NOUT,155) (JTITLE(I),I=1,NCA)
      LI = 1
C                                  WRITE COLUMN NUMBERS
      GO TO (5, 10, 15, 20, 25, 30), IOPT0
    5 WRITE (NOUT,95) (I,I=1,M)
      GO TO 35
   10 WRITE (NOUT,100) (I,I=1,M)
      GO TO 35
   15 WRITE (NOUT,105) (I,I=1,M)
      GO TO 35
   20 WRITE (NOUT,110) (I,I=1,M)
      GO TO 35
   25 WRITE (NOUT,145) (I,I=1,M)
      GO TO 35
   30 WRITE (NOUT,150) (I,I=1,M)
   35 IOPT1 = IPICK(IOPT0)
C                                  WRITE THE MATRIX
      DO 90 I=1,M
         LMAX = LI+I-1
         LEND = LMAX
         IF (I.GT.IOPT1) LEND = LI+IOPT1-1
C                                  WRITE THE FIRST LINE OF THE ROW
         GO TO (40, 40, 45, 45, 50, 50), IOPT0
   40    WRITE (NOUT,115) I, (A(J),J=LI,LEND)
         GO TO 55
   45    WRITE (NOUT,120) I, (A(J),J=LI,LEND)
         GO TO 55
   50    WRITE (NOUT,130) I, (A(J),J=LI,LEND)
   55    CONTINUE
C                                  WRITE THE REMAINING LINES
         IF (I.LE.IOPT1) GO TO 85
         LSTOP = LEND
   60    LSTART = LSTOP+1
         LADD = MIN0(IOPT1-1,I-IOPT1-1)
         LSTOP = MIN0(LSTART+LADD,LMAX)
         GO TO (65, 65, 70, 70, 75, 75), IOPT0
   65    WRITE (NOUT,125) (A(J),J=LSTART,LSTOP)
         GO TO 80
   70    WRITE (NOUT,140) (A(J),J=LSTART,LSTOP)
         GO TO 80
   75    WRITE (NOUT,135) (A(J),J=LSTART,LSTOP)
   80    CONTINUE
         IF (LSTOP.NE.LMAX) GO TO 60
   85    LI = LI+I
   90 CONTINUE
   95 FORMAT (14X, I4, 15X, I4, 15X, I4, 15X, I4, 15X, I4, 15X, I4)
  100 FORMAT (14X, I4, 15X, I4, 15X, I4)
  105 FORMAT (10X, I4, 11X, I4, 11X, I4, 11X, I4, 11X, I4, 11X, I4,
     *11X, I4, 11X, I4)
  110 FORMAT (10X, I4, 11X, I4, 11X, I4, 11X, I4, 11X, I4)
  115 FORMAT (/I4, 6(1X, F18.5))
  120 FORMAT (/I4, 8E15.6)
  125 FORMAT (5X, F18.5, 1X, F18.5, 1X, F18.5, 1X, F18.5, 1X, F18.5,
     *1X, F18.5)
  130 FORMAT (/I4, 5E25.7)
  135 FORMAT (4X, 5E25.7)
  140 FORMAT (4X, 8E15.6)
  145 FORMAT (14X, I4, 21X, I4, 21X, I4, 21X, I4, 21X, I4)
  150 FORMAT (14X, I4, 21X, I4, 21X, I4)
  155 FORMAT (1X, 20A1)
      RETURN
      END
