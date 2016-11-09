C   IMSL ROUTINE NAME   - GTMNT
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - MOMENTS AND STANDARDIZED MOMENTS OF UNIFORM
C                           RANDOM NUMBERS
C
C   USAGE               - CALL GTMNT (RBAR,R,N,NR,STAT,IER)
C
C   ARGUMENTS    RBAR   - INPUT MEAN OF THE SEQUENCE. RBAR IS USED ONLY
C                           ON INITIAL ENTRY TO GTMNT AND MUST HAVE BEEN
C                           OBTAINED PRIOR TO ENTRY
C                           (SEE BASIC STATISTICS CHAPTER)
C                R      - INPUT VECTOR CONTAINING SEQUENCE OF UNIFORM
C                           (0,1) RANDOM NUMBERS.
C                N      - INPUT TOTAL LENGTH OF SEQUENCE R.
C                NR     - INPUT LENGTH OF SEQUENCE R IN CORE ON THIS
C                           CALL.
C                STAT   - OUTPUT VECTOR OF LENGTH 4. THE I-TH ELEMENT
C                           OF STAT CONTAINS, WHEN
C                           I=1, THE STANDARDIZED MEAN
C                         NOTE - ON ALL ENTRIES EXCEPT THE FINAL ENTRY,
C                             STAT(1) MUST BE ZERO. ON THE FINAL ENTRY,
C                             STAT(1) MUST BE NON-ZERO. (SEE REMARKS)
C                           I=2, THE ESTIMATED VARIANCE. ON THE FIRST
C                             ENTRY TO GTMNT, STAT(2) MUST BE ZERO.
C                           I=3, THE STANDARDIZED ESTIMATE OF THE
C                             VARIANCE
C                           I=4, THE ESTIMATED THIRD MOMENT. ON THE
C                             FIRST ENTRY TO GTMNT, STAT(4) MUST BE
C                             ZERO.
C                             (EXPECTED VALUE = 0.25)
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER = 129 INDICATES THAT THE INPUT NR IS
C                             LESS THAN 1 OR THE TOTAL, N, IS LESS
C                             THAN 2.
C
C   REQD. IMSL ROUTINES - UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS      STAT(1) IS USED AS A SWITCH TO CAUSE FINAL CALCULATIONS
C                TO BE PERFORMED ON THE FINAL CALL TO GTMNT WHEN MUL-
C                TIPLE CALLS ARE REQUIRED TO EXHAUST THE SEQUENCE OF
C                RANDOM NUMBERS. THE ELEMENTS OF STAT ARE NOT CALCULATED
C                UNTIL ALL THE ELEMENTS OF R HAVE BEEN ENTERED. THE
C                CALLING PROGRAM MUST ASSIGN THE VALUE ZERO TO STAT(1)
C                PRIOR TO ALL INVOCATIONS EXCEPT THE FINAL CALL. AT THAT
C                TIME STAT(1) MUST HAVE BEEN ASSIGNED A NON-ZERO VALUE
C                SO THAT THE ROUTINE MAY BE TRIGGERED TO COMPUTE THE
C                FINAL RESULTS. IF ONLY ONE CALL TO THE ROUTINE IS
C                NECESSARY (THE TOTAL SEQUENCE R IS IN MEMORY), STAT(1)
C                MUST BE NON-ZERO ON THAT CALL TO GTMNT.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE GTMNT  (RBAR,R,N,NR,STAT,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N,NR,IER
      REAL               RBAR,R(1),STAT(4)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            NR1,I
      REAL               XNR,CON1,PT083,XMU,XSIG,FIVFIV,TWOPT7,ATEPT3
      DOUBLE PRECISION   SUM2,SUM4
      DATA               CON1/3.464102/
      DATA               PT083/.8333333E-1/
      DATA               FIVFIV/5.555556/
      DATA               TWOPT7/2.777778/
      DATA               ATEPT3/8.333333/
      DATA               NR1/0/
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
    5 SUM2 = 0.D0
      SUM4 = 0.D0
      IF (NR .LT. 1) GO TO 15
      DO 10  I=1,NR
         SUM2 = SUM2 + (DBLE(R(I))-DBLE(RBAR))**2
         SUM4 = SUM4 + DBLE(R(I))**3
   10 CONTINUE
      STAT(2) = STAT(2) + SUM2
      STAT(4) = STAT(4) + SUM4
C                                  CHECK FOR FURTHER ENTRIES
      IF (STAT(1) .EQ. 0.0) GO TO 9005
      NR1 = N
      IF (NR1 .LE. 1) GO TO 15
      XNR = NR1
C                                  CALCULATE STD. 1ST MOMENT
      STAT(1) = CON1 * SQRT(XNR) * (RBAR-0.5)
C                                  CALCULATE 2ND MOMENT
      STAT(2) = STAT(2)/(XNR-1.0)
C                                  CALCULATE STD. 2ND MOMENT
      XMU = (PT083 * (XNR-1.0))/XNR
      XSIG = (FIVFIV + TWOPT7/XNR - ATEPT3 /
     1   (XNR*XNR)) * (1.E-3/XNR)
      XSIG = SQRT(XSIG)
      STAT(3) = (STAT(2)-XMU)/XSIG
C                                  CALCULATE 3RD MOMENT
      STAT(4) = STAT(4)/XNR
      GO TO 9005
   15 IER = 129
 9000 CONTINUE
      CALL UERTST(IER,'HGTMNT ')
 9005 RETURN
      END
