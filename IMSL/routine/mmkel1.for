C   IMSL ROUTINE NAME   - MMKEL1
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - KELVIN FUNCTIONS OF THE FIRST KIND, (BER,BEI),
C                           AND OF THE SECOND KIND, (KER,KEI), OF ORDER
C                           ONE
C
C   USAGE               - CALL MMKEL1 (X,BER1,BEI1,XKER1,XKEI1,IER)
C
C   ARGUMENTS    X      - INPUT ARGUMENT. IF X IS NEGATIVE, A WARNING
C                           ERROR IS PRODUCED AND VALUES OF POSITIVE
C                           MACHINE INFINITY WILL BE RETURNED FOR XKER1
C                           AND XKEI1.
C                BER1   - OUTPUT VALUE OF BER X . (ORDER ONE)
C                BEI1   - OUTPUT VALUE OF BEI X . (ORDER ONE)
C                XKER1  - OUTPUT VALUE OF KER X . (ORDER ONE)
C                XKEI1  - OUTPUT VALUE OF KEI X . (ORDER ONE)
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER = 129 INDICATES THAT THE ABSOLUTE VALUE
C                             OF X WAS GREATER THAN 119. BER1 AND BEI1
C                             ARE SET TO ZERO. IF X IS NON-NEGATIVE,
C                             XKER1 AND XKEI1 ARE ALSO SET TO ZERO.
C                             OTHERWISE, XKER1 AND XKEI1 ARE SET TO
C                             POSITIVE MACHINE INFINITY.
C                         WARNING ERROR
C                           IER = 34 INDICATES THAT X IS NEGATIVE.
C                             XKER1 AND XKEI1 WILL BE RETURNED AS
C                             POSITIVE MACHINE INFINITY.
C
C   REQD. IMSL ROUTINES - MMKELD,MMKEL0,UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE MMKEL1 (X,BER1,BEI1,XKER1,XKEI1,IER)
C
      DOUBLE PRECISION   BER1,BEI1,BERP,BEIP,RT2,X,XINF,XKEIP,XKEI1,
     *                   XKERP,XKER1,ZERO,ZMAX
      DATA               XINF/0.7237005577332260D+76/,ZERO/0.0D0/
      DATA               RT2/.7071067811865475D0/
      DATA               ZMAX/119.D0/
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
      IF (X .EQ. ZERO) GO TO 15
      IF (DABS(X) .GT. ZMAX) GO TO 10
      CALL MMKELD(X,BERP,BEIP,XKERP,XKEIP,IER)
      BEI1 = (BERP+BEIP) *RT2
      BER1 = (BERP - BEIP) * RT2
      IF (X .LT. ZERO) GO TO 5
      XKEI1 = (XKERP + XKEIP) * RT2
      XKER1 = (XKERP - XKEIP) * RT2
      GO TO 9005
C                                  ARGUMENT IS NEGATIVE
    5 XKER1 = XINF
      XKEI1 = XINF
      IER = 34
      GO TO 9000
   10 BEI1 = ZERO
      BER1 = ZERO
      XKER1 = ZERO
      XKEI1 = ZERO
      IER = 129
      IF (X .GT. ZERO) GO TO 9000
      XKER1 = XINF
      XKEI1 = XINF
      GO TO 9000
C                                  ARGUMENT IS 0.0
   15 BEI1 = ZERO
      BER1 = ZERO
      XKER1 = -XINF
      XKEI1 = -XINF
      GO TO 9005
 9000 CONTINUE
      CALL UERTST(IER,6HMMKEL1)
 9005 RETURN
      END
