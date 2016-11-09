C   IMSL ROUTINE NAME   - USWBS
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - NOVEMBER 1, 1984
C
C   PURPOSE             - PRINT A MATRIX STORED IN BAND SYMMETRIC
C                           STORAGE MODE.
C
C   USAGE               - CALL USWBS (ITITLE,NC,A,IA,M,NLC,WK,IOPT)
C
C   ARGUMENTS    ITITLE - A CHARACTER STRING TO PROVIDE A TITLE. (INPUT)
C                           THE LENGTH OF ITITLE MUST NOT EXCEED 20.
C                NC     - LENGTH OF ITITLE.  (INPUT)
C                           IF NC.EQ.0 THE TITLE  IS NOT PRINTED.
C                A      - THE MATRIX TO BE PRINTED. (INPUT)
C                           A IS STORED IN BAND SYMMETRIC STORAGE MODE
C                           AND HAS DIMENSION M BY (2*NLC)+1.
C                IA     - ROW DIMENSION OF MATRIX A EXACTLY AS SPECIFIED
C                           IN THE DIMENSION STATEMENT IN THE CALLING
C                           PROGRAM. (INPUT)
C                M      - NUMBER OF ROWS OF MATRIX A TO BE PRINTED.
C                           (INPUT)
C                NLC    - NUMBER OF UPPER OR LOWER CODIAGONALS IN
C                           MATRIX A. (INPUT)
C                WK     - WORK AREA VECTOR OF LENGTH M.
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
      SUBROUTINE USWBS (ITITLE,NC,A,IA,M,NLC,WK,IOPT)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            NC,IA,M,NLC,IOPT
      REAL               A(IA,1),WK(M)
      INTEGER            ITITLE(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,ICNT,II,IK,IOPT0,IOPT1,IPICK(6),J,JBEG,
     *                   JTITLE(20),K,KK,L,LBEG,LEND,MM,MTIME,N2,NCA,
     *                   NCAMTB,NCC,NCMAX,NIN,NLCP1,NLCP2,NOUT
      REAL               ZERO
      DATA               IPICK /6,3,8,5,5,3/
      DATA               NCMAX /20/
      DATA               ZERO /0.0/
C                                  WRITE MATRIX TITLE
C                                  FIRST EXECUTABLE STATEMENT
      NLCP1 = NLC+1
      NLCP2 = NLC+2
      ICNT = 1
C                                  SET IOPT0. IF IOPT IS OUT OF RANGE
C                                  USE 6
      IOPT0 = 6
      IF (IOPT.GE.1 .AND. IOPT.LE.6) IOPT0 = IOPT
C                                  WRITE TITLE
      NCA = MIN0(NC,NCMAX)
      CALL USPKD(ITITLE,NCA,JTITLE,NCAMTB)
      CALL UGETIO(1,NIN,NOUT)
      IF (NCA.GT.0) WRITE (NOUT,160) (JTITLE(I),I=1,NCA)
C                                  WRITE COLUMN NUMBERS
      GO TO (5, 10, 15, 20, 25, 30), IOPT0
    5 WRITE (NOUT,100) (I,I=1,M)
      GO TO 35
   10 WRITE (NOUT,105) (I,I=1,M)
      GO TO 35
   15 WRITE (NOUT,110) (I,I=1,M)
      GO TO 35
   20 WRITE (NOUT,115) (I,I=1,M)
      GO TO 35
   25 WRITE (NOUT,150) (I,I=1,M)
      GO TO 35
   30 WRITE (NOUT,155) (I,I=1,M)
   35 IOPT1 = IPICK(IOPT0)
C                                  WRITE THE MATRIX
      DO 95 I=1,M
         JBEG = MAX0(NLCP2-I,1)
         KK = 1
         IF (I.LE.NLCP1) GO TO 45
         DO 40 IK=1,ICNT
            WK(KK) = ZERO
            KK = KK+1
   40    CONTINUE
         ICNT = ICNT+1
   45    DO 50 IK=JBEG,NLCP1
            WK(KK) = A(I,IK)
            KK = KK+1
   50    CONTINUE
         LEND = IOPT1
         L = MIN0(IOPT1,M)
         II = MIN0(L,I)
C                                  WRITE THE FIRST LINE OF THE ROW
         GO TO (55, 55, 60, 60, 65, 65), IOPT0
   55    WRITE (NOUT,120) I, (WK(J),J=1,II)
         GO TO 70
   60    WRITE (NOUT,125) I, (WK(J),J=1,II)
         GO TO 70
   65    WRITE (NOUT,135) I, (WK(J),J=1,II)
   70    MM = I-L
C                                  WRITE THE REMAINING LINES
         IF (I.LE.IOPT1) GO TO 95
         MTIME = I/IOPT1
         IF (MTIME*IOPT1.EQ.I) MTIME = MTIME-1
         DO 90 K=1,MTIME
            LBEG = LEND+1
            L = MIN0(IOPT1,MM)
            MM = MM-L
            LEND = LEND+L
            GO TO (75, 75, 80, 80, 85, 85), IOPT0
   75       WRITE (NOUT,130) (WK(J),J=LBEG,LEND)
            GO TO 90
   80       WRITE (NOUT,145) (WK(J),J=LBEG,LEND)
            GO TO 90
   85       WRITE (NOUT,140) (WK(J),J=LBEG,LEND)
   90    CONTINUE
   95 CONTINUE
  100 FORMAT (14X, I4, 15X, I4, 15X, I4, 15X, I4, 15X, I4, 15X, I4)
  105 FORMAT (14X, I4, 15X, I4, 15X, I4)
  110 FORMAT (10X, I4, 11X, I4, 11X, I4, 11X, I4, 11X, I4, 11X, I4,
     *11X, I4, 11X, I4)
  115 FORMAT (10X, I4, 11X, I4, 11X, I4, 11X, I4, 11X, I4)
  120 FORMAT (/I4, 6(1X, F18.5))
  125 FORMAT (/I4, 8E15.6)
  130 FORMAT (5X, F18.5, 1X, F18.5, 1X, F18.5, 1X, F18.5, 1X, F18.5,
     *1X, F18.5)
  135 FORMAT (/I4, 5E25.7)
  140 FORMAT (4X, 5E25.7)
  145 FORMAT (4X, 8E15.6)
  150 FORMAT (14X, I4, 21X, I4, 21X, I4, 21X, I4, 21X, I4)
  155 FORMAT (14X, I4, 21X, I4, 21X, I4)
  160 FORMAT (1X, 20A1)
      RETURN
      END
