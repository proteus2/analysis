C   IMSL ROUTINE NAME   - USWCV
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JULY 1, 1983
C
C   PURPOSE             - PRINT A COMPLEX VECTOR
C
C   USAGE               - CALL USWCV (ITITLE,NC,A,M,INC,IOPT)
C
C   ARGUMENTS    ITITLE - A CHARACTER STRING TO PROVIDE A TITLE. (INPUT)
C                           THE LENGTH OF ITITLE MUST NOT EXCEED 20.
C                NC     - LENGTH OF ITITLE.  (INPUT)
C                           IF NC.EQ.0 THE TITLE  IS NOT PRINTED.
C                A      - COMPLEX VECTOR TO BE PRINTED. (INPUT)
C                M      - NUMBER OF ELEMENTS OF VECTOR A TO BE PRINTED.
C                           (INPUT)
C                INC    - DISPLACEMENT BETWEEN ELEMENTS OF THE VECTOR
C                           TO BE PRINTED. (INPUT) USWCV PRINTS
C                           A(1+(I-1)*INC) FOR I=1,...,M. SEE REMARKS.
C                IOPT   - OPTION INDICATING THE FORMAT STATEMENT TO BE
C                           USED. (INPUT)
C
C                           OPTIONS FOR    OPTIONS FOR
C                           129 COLUMNS    80 COLUMNS     FORMAT
C                           -----------    -----------    ------
C
C                                1              2         F12.4
C                                3              4         E12.4
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
C                STORES MATRICES COLUMN WISE. HENCE, IF A IS A 2-
C                DIMENSIONAL ARRAY AND A ROW OF A IS TO BE PRINTED,
C                INC= THE ROW DIMENSION OF A EXACTLY AS SPECIFIED IN
C                THE CALLING PROGRAM. IF A IS A 2-DIMENSIONAL ARRAY AND
C                A COLUMN OF A IS TO BE PRINTED, INC=1.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE USWCV (ITITLE,NC,A,M,INC,IOPT)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            NC,M,INC,IOPT
      COMPLEX            A(M)
      INTEGER            ITITLE(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,IBEG2,IEND,IEND2,IOPT0,IOPT1,J,K,L,LBEG,
     *                   LEND,MM,MTIME,NCA,NCAMTB,NIN,NOUT,JTITLE(20)
      INTEGER            ILP,IRP,NCMAX,IPICK(6)
      DATA               IPICK /4,2,4,2,2,1/
      DATA               NCMAX /20/
      DATA               ILP /1H(/,IRP /1H)/
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
      IF (NCA.GT.0) WRITE (NOUT,140) (JTITLE(I),I=1,NCA)
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
C                                  WRITE THE FIRST LINE OF THE ROW
      IEND = 1+(L-1)*INCTMP
      GO TO (5, 10, 15, 20, 25, 30), IOPT0
    5 WRITE (NOUT,80) (ILP,A(J),IRP,J=1,IEND,INCTMP)
      GO TO 35
   10 WRITE (NOUT,85) (ILP,A(J),IRP,J=1,IEND,INCTMP)
      GO TO 35
   15 WRITE (NOUT,130) (ILP,A(J),IRP,J=1,IEND,INCTMP)
      GO TO 35
   20 WRITE (NOUT,135) (ILP,A(J),IRP,J=1,IEND,INCTMP)
      GO TO 35
   25 WRITE (NOUT,100) (ILP,A(J),IRP,J=1,IEND,INCTMP)
      GO TO 35
   30 WRITE (NOUT,105) (ILP,A(J),IRP,J=1,IEND,INCTMP)
   35 MM = M-L
C                                  WRITE THE REMAINING LINES
      IF (MTIME.EQ.0) GO TO 75
      DO 70 K=1,MTIME
         LBEG = LEND+1
         L = MIN0(IOPT1,MM)
         MM = MM-L
         LEND = LEND+L
         IBEG2 = 1+(LBEG-1)*INCTMP
         IEND2 = 1+(LEND-1)*INCTMP
         GO TO (40, 45, 50, 55, 60, 65), IOPT0
   40    WRITE (NOUT,90) (ILP,A(J),IRP,J=IBEG2,IEND2,INCTMP)
         GO TO 70
   45    WRITE (NOUT,95) (ILP,A(J),IRP,J=IBEG2,IEND2,INCTMP)
         GO TO 70
   50    WRITE (NOUT,120) (ILP,A(J),IRP,J=IBEG2,IEND2,INCTMP)
         GO TO 70
   55    WRITE (NOUT,125) (ILP,A(J),IRP,J=IBEG2,IEND2,INCTMP)
         GO TO 70
   60    WRITE (NOUT,110) (ILP,A(J),IRP,J=IBEG2,IEND2,INCTMP)
         GO TO 70
   65    WRITE (NOUT,115) (ILP,A(J),IRP,J=IBEG2,IEND2,INCTMP)
   70 CONTINUE
   75 CONTINUE
   80 FORMAT (/5X, 4(A1, F12.4, 1H,, F12.4, A1, 2X))
   85 FORMAT (/5X, 2(A1, F12.4, 1H,, F12.4, A1, 2X))
   90 FORMAT (5X, 4(A1, F12.4, 1H,, F12.4, A1, 2X))
   95 FORMAT (5X, 2(A1, F12.4, 1H,, F12.4, A1, 2X))
  100 FORMAT (/5X, 2(A1, E25.7, 1H,, E25.7, A1, 2X))
  105 FORMAT (/5X, A1, E25.7, 1H,, E25.7, A1)
  110 FORMAT (5X, 2(A1, E25.7, 1H,, E25.7, A1, 2X))
  115 FORMAT (5X, A1, E25.7, 1H,, E25.7, A1)
  120 FORMAT (5X, 4(A1, E12.4, 1H,, E12.4, A1, 2X))
  125 FORMAT (5X, 2(A1, E12.4, 1H,, E12.4, A1, 2X))
  130 FORMAT (/5X, 4(A1, E12.4, 1H,, E12.4, A1, 2X))
  135 FORMAT (/5X, 2(A1, E12.4, 1H,, E12.4, A1, 2X))
  140 FORMAT (1X, 20A1)
      RETURN
      END
