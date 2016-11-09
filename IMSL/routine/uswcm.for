C   IMSL ROUTINE NAME   - USWCM
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1982
C
C   PURPOSE             - PRINT A COMPLEX MATRIX STORED IN FULL
C                           STORAGE MODE
C
C   USAGE               - CALL USWCM (ITITLE,NC,A,IA,N,M,IOPT)
C
C   ARGUMENTS    ITITLE - A CHARACTER STRING TO PROVIDE A TITLE. (INPUT)
C                           THE LENGTH OF ITITLE MUST NOT EXCEED 20.
C                NC     - LENGTH OF ITITLE.  (INPUT)
C                           IF NC.EQ.0 THE TITLE  IS NOT PRINTED.
C                A      - COMPLEX MATRIX TO BE PRINTED. (INPUT)
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
C   REMARKS  1.  IF AN ASTERISK APPEARS IN THE HIGH ORDER POSITION
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
      SUBROUTINE USWCM (ITITLE,NC,A,IA,N,M,IOPT)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            NC,IA,N,M,IOPT
      COMPLEX            A(IA,M)
      INTEGER            ITITLE(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,ILP,IOPT0,IOPT1,IPICK(6),IRP,J,JTITLE(20),K,
     *                   L,LBEG,LEND,MM,MTIME,NCA,NCAMTB,NCMAX,NIN,
     *                   NOUT
      DATA               IPICK /4,2,4,2,2,1/
      DATA               NCMAX /20/
      DATA               ILP /1H(/,IRP /1H)/
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
      IF (NCA.GT.0) WRITE (NOUT,205) (JTITLE(I),I=1,NCA)
C                                  WRITE COLUMN NUMBERS
      GO TO (5, 10, 15, 20, 25, 30), IOPT0
    5 WRITE (NOUT,115) (I,I=1,M)
      GO TO 35
   10 WRITE (NOUT,120) (I,I=1,M)
      GO TO 35
   15 WRITE (NOUT,125) (I,I=1,M)
      GO TO 35
   20 WRITE (NOUT,130) (I,I=1,M)
      GO TO 35
   25 WRITE (NOUT,195) (I,I=1,M)
      GO TO 35
   30 WRITE (NOUT,200) (I,I=1,M)
   35 IOPT1 = IPICK(IOPT0)
      MTIME = M/IOPT1
      IF (MTIME*IOPT1.EQ.M) MTIME = MTIME-1
C                                  WRITE THE MATRIX
      DO 110 I=1,N
         LEND = IOPT1
         L = MIN0(IOPT1,M)
C                                  WRITE THE FIRST LINE OF THE ROW
         GO TO (40, 45, 50, 55, 60, 65), IOPT0
   40    WRITE (NOUT,135) I, (ILP,A(I,J),IRP,J=1,L)
         GO TO 70
   45    WRITE (NOUT,140) I, (ILP,A(I,J),IRP,J=1,L)
         GO TO 70
   50    WRITE (NOUT,185) I, (ILP,A(I,J),IRP,J=1,L)
         GO TO 70
   55    WRITE (NOUT,190) I, (ILP,A(I,J),IRP,J=1,L)
         GO TO 70
   60    WRITE (NOUT,155) I, (ILP,A(I,J),IRP,J=1,L)
         GO TO 70
   65    WRITE (NOUT,160) I, (ILP,A(I,J),IRP,J=1,L)
   70    MM = M-L
C                                  WRITE THE REMAINING LINES
         IF (MTIME.EQ.0) GO TO 110
         DO 105 K=1,MTIME
            LBEG = LEND+1
            L = MIN0(IOPT1,MM)
            MM = MM-L
            LEND = LEND+L
            GO TO (75, 80, 85, 90, 95, 100), IOPT0
   75       WRITE (NOUT,145) (ILP,A(I,J),IRP,J=LBEG,LEND)
            GO TO 105
   80       WRITE (NOUT,150) (ILP,A(I,J),IRP,J=LBEG,LEND)
            GO TO 105
   85       WRITE (NOUT,175) (ILP,A(I,J),IRP,J=LBEG,LEND)
            GO TO 105
   90       WRITE (NOUT,180) (ILP,A(I,J),IRP,J=LBEG,LEND)
            GO TO 105
   95       WRITE (NOUT,165) (ILP,A(I,J),IRP,J=LBEG,LEND)
            GO TO 105
  100       WRITE (NOUT,170) (ILP,A(I,J),IRP,J=LBEG,LEND)
  105    CONTINUE
  110 CONTINUE
  115 FORMAT (15X, I4, 25X, I4, 25X, I4, 25X, I4)
  120 FORMAT (15X, I4, 25X, I4)
  125 FORMAT (15X, I4, 25X, I4, 25X, I4, 25X, I4)
  130 FORMAT (15X, I4, 25X, I4)
  135 FORMAT (/I4, 1X, 4(A1, F12.4, 1H,, F12.4, A1, 2X))
  140 FORMAT (/I4, 1X, 2(A1, F12.4, 1H,, F12.4, A1, 2X))
  145 FORMAT (5X, 4(A1, F12.4, 1H,, F12.4, A1, 2X))
  150 FORMAT (5X, 2(A1, F12.4, 1H,, F12.4, A1, 2X))
  155 FORMAT (/I4, 1X, 2(A1, E25.7, 1H,, E25.7, A1, 2X))
  160 FORMAT (/I4, 1X, A1, E25.7, 1H,, E25.7, A1)
  165 FORMAT (5X, 2(A1, E25.7, 1H,, E25.7, A1, 2X))
  170 FORMAT (5X, A1, E25.7, 1H,, E25.7, A1)
  175 FORMAT (5X, 4(A1, E12.4, 1H,, E12.4, A1, 2X))
  180 FORMAT (5X, 2(A1, E12.4, 1H,, E12.4, A1, 2X))
  185 FORMAT (/I4, 1X, 4(A1, E12.4, 1H,, E12.4, A1, 2X))
  190 FORMAT (/I4, 1X, 2(A1, E12.4, 1H,, E12.4, A1, 2X))
  195 FORMAT (28X, I4, 51X, I4)
  200 FORMAT (28X, I4)
  205 FORMAT (1X, 20A1)
      RETURN
      END
