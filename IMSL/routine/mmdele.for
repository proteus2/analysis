C   IMSL ROUTINE NAME   - MMDELE
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - NOVEMBER 1, 1979
C
C   PURPOSE             - COMPLETE ELLIPTIC INTEGRAL OF THE SECOND KIND
C
C   USAGE               - FUNCTION MMDELE (IOPT,ARG,IER)
C
C   ARGUMENTS    MMDELE - OUTPUT VALUE OF THE INTEGRAL. MMDELE MUST BE
C                           TYPED APPROPRIATELY IN THE CALLING PROGRAM.
C                           (SEE PRECISION/HARDWARE SECTION.)
C                IOPT   - INPUT OPTION.
C                         IF IOPT = 1, THE INTEGRAL (FROM 0 TO PI/2) OF
C                           SQRT(1-ARG*(SIN(PHI))**2) D(PHI) WILL BE
C                           EVALUATED. ARG MUST BE GREATER THAN OR
C                           EQUAL TO 0.0 AND LESS THAN OR EQUAL TO 1.
C                         IF IOPT = 2, THE INTEGRAL (FROM 0 TO PI/2) OF
C                           SQRT(1-ARG**2*(SIN(PHI))**2) D(PHI) WILL BE
C                           EVALUATED. THE ABSOLUTE VALUE OF ARG MUST
C                           BE LESS THAN OR EQUAL TO 1.
C                         IF IOPT = 3, THE INTEGRAL (FROM 0 TO PI/2) OF
C                           SQRT(1-(1-ARG)(SIN(PHI))**2) D(PHI) WILL BE
C                           EVALUATED. THE ABSOLUTE VALUE OF ARG MUST
C                           BE LESS THAN OR EQUAL TO 1.
C                ARG    - INPUT PARAMETER. SEE DESCRIPTION OF IOPT.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER = 129 INDICATES THAT IOPT IS LESS THAN
C                             1 OR GREATER THAN 3. MMDELE IS SET TO
C                             MACHINE INFINITY.
C                           IER = 130 INDICATES THAT ARG IS OUT OF
C                             RANGE. MMDELE IS SET TO MACHINE INFINITY.
C
C   REQD. IMSL ROUTINES - UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      DOUBLE PRECISION FUNCTION MMDELE (IOPT,ARG,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            IOPT,IER
      DOUBLE PRECISION   ARG
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I
      DOUBLE PRECISION   CAY,CAYSQ,EPS,ETA,ONE,ONEP,SUMA,SUMB,XINF,ZERO
      DOUBLE PRECISION   A(10),B(10)
      DATA               ZERO/0.D0/,ONEP/1.D0/
      DATA               XINF/Z7FFFFFFFFFFFFFFF/,
     *                   EPS/Z3210000000000000/
      DATA               A(1)/Z3D9CBA0842E0E54E/,A(2)/Z3EA1C69D130B9BD4/
      DATA               A(3)/Z3F23621029FBF8E3/,A(4)/Z3F2C1DD4540DD62D/
      DATA               A(5)/Z3F200849BAEEE829/,A(6)/Z3F1F1C085DC262C4/
      DATA               A(7)/Z3F2F6399B5A9EBA6/,A(8)/Z3F596C5052952415/
      DATA               A(9)/Z3FE8AC904E91BFF1/,
     *                   A(10)/Z40717217F7D2D940/
      DATA               B(1)/Z3D216823413270CB/,B(2)/Z3E40DEA596008FA8/
      DATA               B(3)/Z3F1A589833B593F2/,B(4)/Z3F44D44BAB6BEF41/
      DATA               B(5)/Z3F6B1708F175E9BA/,B(6)/Z3F892137826299E3/
      DATA               B(7)/Z3FAEF8F46DAE9025/,B(8)/Z3FEFFFE82D86CCB2/
      DATA               B(9)/Z4017FFFFFECD75A9/,
     *                   B(10)/Z403FFFFFFFFFE45A/
      DATA               ONE/Z4110000000000001/
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
      CAYSQ = ARG
      IF (IOPT .LT. 1 .OR. IOPT .GT. 3) GO TO 40
      GO TO (5,25,30),IOPT
C                                  IOPT = 1
    5 IF (CAYSQ .LT. ZERO) GO TO 35
      ETA = ONEP - CAYSQ
   10 IF (ETA .LT. ZERO) GO TO 35
      IF (ETA .LT. EPS) GO TO 20
      SUMA = ZERO
      SUMB = ZERO
      DO 15 I = 1,10
         SUMA = (SUMA + A(I)) * ETA
         SUMB = (SUMB + B(I)) * ETA
   15 CONTINUE
      MMDELE = SUMA - DLOG(ETA) * SUMB
      MMDELE = MMDELE + ONE
      GO TO 9005
   20 MMDELE = ONEP
      GO TO 9005
C                                  IOPT = 2
C                                  CALCULATION OF 1-ARG**2 IS REFORMU-
C                                  LATED TO PRESERVE ACCURACY
   25 CAY = DABS(CAYSQ)
      ETA = ONEP - CAY
      ETA = ETA + CAY * ETA
      GO TO 10
C                                  IOPT = 3
   30 IF (CAYSQ .GT. ONEP) GO TO 35
      ETA = CAYSQ
      GO TO 10
C                                  TERMINAL ERROR - ARG IS OUT OF RANGE
   35 IER = 130
      MMDELE = XINF
      GO TO 9000
C                                  TERMINAL ERROR - IOPT IS OUT OF
C                                  RANGE
   40 MMDELE = XINF
      IER = 129
 9000 CONTINUE
      CALL UERTST (IER,6HMMDELE)
 9005 RETURN
      END
