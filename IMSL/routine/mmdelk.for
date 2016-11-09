C   IMSL ROUTINE NAME   - MMDELK
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - COMPLETE ELLIPTIC INTEGRAL OF THE FIRST KIND
C
C   USAGE               - FUNCTION MMDELK (IOPT,ARG,IER)
C
C   ARGUMENTS    MMDELK - OUTPUT VALUE OF THE INTEGRAL. MMDELK MUST BE
C                           TYPED APPROPRIATELY IN THE CALLING PROGRAM.
C                           (SEE PRECISION/HARDWARE SECTION.)
C                IOPT   - INPUT OPTION.
C                         IF IOPT = 1, THE INTEGRAL (FROM 0 TO PI/2) OF
C                           1./SQRT(1.-ARG*(SIN(PHI))**2) D(PHI) WILL
C                           BE EVALUATED. ARG MUST BE GREATER THAN OR
C                           EQUAL TO 0.0 AND LESS THAN 1.
C                         IF IOPT = 2, THE INTEGRAL (FROM 0 TO PI/2) OF
C                           1./SQRT(1.-ARG**2*(SIN(PHI))**2) D(PHI) WILL
C                           BE EVALUATED. THE ABSOLUTE VALUE OF ARG MUST
C                           BE LESS THAN 1.
C                         IF IOPT = 3, THE INTEGRAL (FROM 0 TO PI/2) OF
C                           1./SQRT(1-(1-ARG)(SIN(PHI)**2)) D(PHI) WILL
C                           BE EVALUATED. ARG MUST BE GREATER THAN 0.0
C                           AND LESS THAN OR EQUAL TO 1.
C                ARG    - INPUT PARAMETER. SEE IOPT DESCRIPTION.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL
C                           IER = 129 INDICATES THAT IOPT IS LESS THAN
C                             1 OR GREATER THAN 3. MMDELK IS SET TO
C                             MACHINE INFINITY.
C                           IER = 130 INDICATES THAT ARG IS OUT OF
C                             RANGE. MMDELK IS SET TO MACHINE INFINITY.
C
C   REQD. IMSL ROUTINES - UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      DOUBLE PRECISION FUNCTION MMDELK (IOPT,ARG,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            IOPT,IER
      DOUBLE PRECISION   ARG
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I
      DOUBLE PRECISION   CAY,CAYSQ,ZERO,ONE,EPS,SUMA,SUMB,ETA,XINF
      DOUBLE PRECISION   A(11),B(11)
      DATA               ZERO/0.0D0/,ONE/1.0D0/
      DATA               XINF/Z7FFFFFFFFFFFFFFF/
      DATA               EPS/Z3110000000000000/
      DATA               A(1)/Z3D92136ADBAAF320/,A(2)/Z3E968323C78B29EC/
      DATA               A(3)/Z3F20C7C205DDCBDF/,A(4)/Z3F285759D56E0A36/
      DATA               A(5)/Z3F1C0C8DA6A336AF/,A(6)/Z3F194FCF1C842636/
      DATA               A(7)/Z3F2400C6796EBAAD/,A(8)/Z3F3D2FA479710958/
      DATA               A(9)/Z3F7E816C52ADCD51/,
     *                   A(10)/Z4018B90BFBE9E073/,
     *                   A(11)/Z41162E42FEFA39F0/
      DATA               B(1)/Z3D1F249BC455C765/,B(2)/Z3E3C651F04C5E0BA/
      DATA               B(3)/Z3F187817B9824765/,B(4)/Z3F3F9D5F662AB575/
      DATA               B(5)/Z3F62066E10DB00C3/,B(6)/Z3F7B643194EADC14/
      DATA               B(7)/Z3F9919669DC978A3/,B(8)/Z3FC7FFE9B1E21075/
      DATA               B(9)/Z4011FFFFFEE111D4/
      DATA               B(10)/Z401FFFFFFFFFE621/
      DATA               B(11)/Z4080000000000001/
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
      CAYSQ = ARG
      IF (IOPT .LT. 1 .OR. IOPT .GT. 3) GO TO 40
      GO TO (5,25,30),IOPT
C                                  IOPT = 1
    5 IF (CAYSQ .LT. ZERO) GO TO 35
      ETA = ONE - CAYSQ
   10 IF (ETA .LE. ZERO) GO TO 35
      IF (ETA .LT. EPS) GO TO 20
      SUMA = A(1)
      SUMB = B(1)
      DO 15 I = 2,11
         SUMA = SUMA * ETA + A(I)
         SUMB = SUMB * ETA + B(I)
   15 CONTINUE
      MMDELK = SUMA - DLOG(ETA) * SUMB
      GO TO 9005
C                                  RETURN FOR SMALL ARGUMENT
   20 MMDELK = A(11) - DLOG(ETA) * B(11)
      GO TO 9005
C                                  IOPT = 2, CALCULATION OF 1-CAY**2 IS
C                                  REFORMULATED TO PRESERVE ACCURACY
   25 CAY = DABS(CAYSQ)
      ETA = ONE - CAY
      ETA = ETA + CAY * ETA
      GO TO 10
C                                  IOPT = 3
   30 IF (CAYSQ .GT. ONE) GO TO 35
      ETA = CAYSQ
      GO TO 10
C                                  TERMINAL ERROR - ARG IS OUT OF
C                                  RANGE
   35 IER = 130
      MMDELK = XINF
      GO TO 9000
C                                  TERMINAL ERROR - IOPT IS OUT OF RANGE
   40 MMDELK = XINF
      IER = 129
 9000 CONTINUE
      CALL UERTST(IER,6HMMDELK)
 9005 RETURN
      END
