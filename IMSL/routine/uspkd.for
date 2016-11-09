C   IMSL ROUTINE NAME   - USPKD
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - NOVEMBER 1, 1984
C
C   PURPOSE             - NUCLEUS CALLED BY IMSL ROUTINES THAT HAVE
C                           CHARACTER STRING ARGUMENTS
C
C   USAGE               - CALL USPKD  (PACKED,NCHARS,UNPAKD,NCHMTB)
C
C   ARGUMENTS    PACKED - CHARACTER STRING TO BE UNPACKED.(INPUT)
C                NCHARS - LENGTH OF PACKED. (INPUT)  SEE REMARKS.
C                UNPAKD - INTEGER ARRAY TO RECEIVE THE UNPACKED
C                         REPRESENTATION OF THE STRING. (OUTPUT)
C                NCHMTB - NCHARS MINUS TRAILING BLANKS. (OUTPUT)
C
C   REQD. IMSL ROUTINES - NONE
C
C   REMARKS  1.  USPKD UNPACKS A CHARACTER STRING INTO AN INTEGER ARRAY
C                IN (A1) FORMAT.
C            2.  UP TO 129 CHARACTERS MAY BE USED.  ANY IN EXCESS OF
C                THAT ARE IGNORED.
C
C-----------------------------------------------------------------------
      SUBROUTINE USPKD  (PACKED,NCHARS,UNPAKD,NCHMTB)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            NC,NCHARS,NCHMTB
C
      LOGICAL*1          UNPAKD(1),PACKED(1),LBYTE,LBLANK
      INTEGER*2          IBYTE,IBLANK
      EQUIVALENCE (LBYTE,IBYTE)
      DATA               LBLANK /1H /
      DATA               IBYTE /1H /
      DATA               IBLANK /1H /
C                                  INITIALIZE NCHMTB
      NCHMTB = 0
C                                  RETURN IF NCHARS IS LE ZERO
      IF(NCHARS.LE.0) RETURN
C                                  SET NC=NUMBER OF CHARS TO BE DECODED
      NC = MIN0 (129,NCHARS)
      NWORDS = NC*4
      J = 1
      DO 110 I = 1,NWORDS,4
      UNPAKD(I) = PACKED(J)
      UNPAKD(I+1) = LBLANK
      UNPAKD(I+2) = LBLANK
      UNPAKD(I+3) = LBLANK
  110 J = J+1
C                                  CHECK UNPAKD ARRAY AND SET NCHMTB
C                                  BASED ON TRAILING BLANKS FOUND
      DO 200 N = 1,NWORDS,4
         NN = NWORDS - N - 2
         LBYTE = UNPAKD(NN)
         IF(IBYTE .NE. IBLANK) GO TO 210
  200 CONTINUE
      NN = 0
  210 NCHMTB = (NN + 3) / 4
      RETURN
      END

