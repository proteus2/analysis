C   IMSL ROUTINE NAME   - USWFV
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JULY 1, 1983
C
C   PURPOSE             - PRINT A VECTOR.
C
C   USAGE               - CALL USWFV (ITITLE,NC,A,M,INC,IOPT)
C
C   ARGUMENTS    ITITLE - A CHARACTER STRING TO PROVIDE A TITLE. (INPUT)
C                           THE LENGTH OF ITITLE MUST NOT EXCEED 20.
C                NC     - LENGTH OF ITITLE.  (INPUT)
C                           IF NC.EQ.0 THE TITLE  IS NOT PRINTED.
C                A      - THE VECTOR TO BE PRINTED. (INPUT)
C                M      - NUMBER OF ELEMENTS OF VECTOR A TO BE PRINTED.
C                           (INPUT)
C                INC    - DISPLACEMENT BETWEEN ELEMENTS OF THE VECTOR
C                           TO BE PRINTED. (INPUT) USWFV PRINTS
C                           A(1+(I-1)*INC) FOR I=1,...,M. SEE REMARKS.
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
C            4.  IF A IS A 1-DIMENSIONAL ARRAY, THEN INC=1. FORTRAN
C                STORES MATRICES COLUMNWISE. HENCE, IF A IS A 2-
C                DIMENSIONAL ARRAY AND A ROW OF A IS TO BE PRINTED,
C                INC= THE ROW DIMENSION OF A EXACTLY AS SPECIFIED IN
C                THE CALLING PROGRAM. IF A IS A 2-DIMENSIONAL ARRAY AND
C                A COLUMN OF A IS TO BE PRINTED, INC=1.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE USWFV (ITITLE,NC,A,M,INC,IOPT)
C
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            NC,M,INC,IOPT
      REAL               A(1)
      INTEGER            ITITLE(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,IBEG2,IEND,IEND2,IOPT0,IOPT1,IPICK(6),J,
     *                   JTITLE(20),K,L,LBEG,LEND,MM,MTIME,NCA,NCAMTB,
     *                   NCMAX,NIN,NOUT
      DATA               IPICK /6,3,8,5,5,3/
      DATA               NCMAX /20/
C                                  WRITE VECTOR TITLE
C                                  SET IOPT0. IF IOPT IS OUT OF RANGE
C                                  USE 6
C                                  FIRST EXECUTABLE STATEMENT
      IOPT0 = 6
      IF (IOPT.GE.1 .AND. IOPT.LE.6) IOPT0 = IOPT
C                                  WRITE TITLE
      NCA = MIN0(NC,NCMAX)
      CALL USPKD(ITITLE,NCA,JTITLE,NCAMTB)
      CALL UGETIO(1,NIN,NOUT)
      IF (NCA.GT.0) WRITE (NOUT,80) (JTITLE(I),I=1,NCA)
      IOPT1 = IPICK(IOPT0)
      MTIME = M/IOPT1
      IF (MTIME*IOPT1.EQ.M) MTIME = MTIME-1
C                                  WRITE THE VECTOR
      LEND = IOPT1
      L = MIN0(IOPT1,M)
C                                  INITIALIZE INCTMP AND CHECK INC,
C                                    USING 1 IF INC IS LESS THAN 1
      INCTMP = INC
      IF (INC.LT.1) INCTMP = 1
C                                  WRITE THE FIRST LINE OF THE VECTOR
      IEND = 1+(L-1)*INCTMP
      GO TO (5, 5, 10, 10, 15, 15), IOPT0
    5 WRITE (NOUT,50) (A(J),J=1,IEND,INCTMP)
      GO TO 20
   10 WRITE (NOUT,55) (A(J),J=1,IEND,INCTMP)
      GO TO 20
   15 WRITE (NOUT,65) (A(J),J=1,IEND,INCTMP)
   20 MM = M-L
C                                  WRITE THE REMAINING LINES
      IF (MTIME.EQ.0) GO TO 45
      DO 40 K=1,MTIME
         LBEG = LEND+1
         L = MIN0(IOPT1,MM)
         MM = MM-L
         LEND = LEND+L
         IBEG2 = 1+(LBEG-1)*INCTMP
         IEND2 = 1+(LEND-1)*INCTMP
         GO TO (25, 25, 30, 30, 35, 35), IOPT0
   25    WRITE (NOUT,60) (A(J),J=IBEG2,IEND2,INCTMP)
         GO TO 40
   30    WRITE (NOUT,75) (A(J),J=IBEG2,IEND2,INCTMP)
         GO TO 40
   35    WRITE (NOUT,70) (A(J),J=IBEG2,IEND2,INCTMP)
   40 CONTINUE
   45 CONTINUE
   50 FORMAT (/, 4X, 6(1X, F18.5))
   55 FORMAT (/, 4X, 8E15.6)
   60 FORMAT (5X, F18.5, 1X, F18.5, 1X, F18.5, 1X, F18.5, 1X, F18.5,
     *1X, F18.5)
   65 FORMAT (/, 4X, 5E25.7)
   70 FORMAT (4X, 5E25.7)
   75 FORMAT (4X, 8E15.6)
   80 FORMAT (1X, 20A1)
      RETURN
      END
