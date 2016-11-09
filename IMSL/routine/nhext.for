C   IMSL ROUTINE NAME   - NHEXT
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1981
C
C   PURPOSE             - FISHERS EXACT METHOD FOR 2 BY 2 TABLES
C
C   USAGE               - CALL NHEXT (ITAB,IT,P,IER)
C
C   ARGUMENTS    ITAB   - INPUT/OUTPUT 3 BY 3 MATRIX.
C                         ON INPUT, ITAB(I,J) CONTAINS THE INPUT TABLE
C                           OF CATEGORIZATION J OF THE OBSERVATIONS
C                           INTO SAMPLE I, FOR I=1,2 AND J=1,2.
C                         ON OUTPUT, THE THIRD ROW OF ITAB CONTAINS THE
C                           COLUMN TOTALS, THE THIRD COLUMN OF ITAB
C                           CONTAINS THE ROW TOTALS, AND ITAB(3,3)
C                           CONTAINS THE TABLE TOTAL.
C                           THE REST OF ITAB HAS BEEN
C                           MODIFIED SO AS TO CONTAIN THE OPPOSITE TABLE
C                           (EXCEPT IN CERTAIN CASES INVOLVING NEGATIVE
C                           INTEGERS). SEE THE ALGORITHM SECTION IN THE
C                           MANUAL DOCUMENT FOR FURTHER DETAILS.
C                IT     - INPUT ROW DIMENSION OF THE MATRIX ITAB
C                           EXACTLY AS SPECIFIED IN THE DIMENSION
C                           STATEMENT IN THE CALLING PROGRAM.
C                P      - OUTPUT VECTOR OF LENGTH 2.
C                         P(1) CONTAINS THE PROBABILITY FOR A ONE
C                           TAILED TEST, OF ACHIEVING A TABLE AS
C                           EXTREME OR MORE EXTREME THAN ITAB.
C                         P(2) CONTAINS THE PROBABILITY OF AN OUTCOME IN
C                           THE OPPOSITE DIRECTION AS EXTREME OR MORE
C                           EXTREME THAN THE ADJUSTED ITAB.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER=129 MEANS THAT SOME ROW OR COLUMN TOTAL
C                             IN ITAB IS ZERO.
C                           IER=130 MEANS THAT AN ERROR OCCURRED IN
C                             SUBROUTINE MDHYP.
C
C   REQD. IMSL ROUTINES - MDHYP,UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS  1.  OBSERVATIONS CAN BE CATEGORIZED USING IMSL ROUTINES
C                BDCOU1 OR BDCOU2.
C            2.  THE ASSUMPTIONS ON WHICH THIS TEST IS BASED ARE
C                (A)  EACH SAMPLED POPULATION IS INFINITE
C                (B)  EACH UNIT IN EACH POPULATION IS A SUCCESS OR
C                     FAILURE, BUT NOT BOTH (THE MEMBERS OF EACH
C                     DICHOTOMY ARE MUTUALLY EXCLUSIVE AND EXHAUSTIVE)
C                (C)  OBSERVATIONS ARE INDEPENDENT
C
C-----------------------------------------------------------------------
C
      SUBROUTINE NHEXT  (ITAB,IT,P,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            IT,IER,ITAB(IT,3)
      REAL               P(2)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            JER,JX,K,L,N
      REAL               XK,RTAB1,RTAB2
      DOUBLE PRECISION   PEQK,PLEK
C                                  FIRST EXECUTABLE STATEMENT
      IER=0
      JER=0
C                                  COMPUTE ROW, COLUMN, AND TABLE TOTALS
      ITAB(1,3)=ITAB(1,1)+ITAB(1,2)
      ITAB(2,3)=ITAB(2,1)+ITAB(2,2)
      ITAB(3,1)=ITAB(1,1)+ITAB(2,1)
      ITAB(3,2)=ITAB(1,2)+ITAB(2,2)
      ITAB(3,3)=ITAB(3,1)+ITAB(3,2)
      IF (ITAB(1,3) .EQ. 0 .OR. ITAB(2,3) .EQ. 0) GO TO 25
      IF (ITAB(3,1) .EQ. 0 .OR. ITAB(3,2) .EQ. 0) GO TO 25
      L=1
C                                  FIND THE HYPERGEOMETRIC CALL
C                                  PARAMETERS
    5 JX = 1
      RTAB1 = FLOAT(ITAB(1,1))*FLOAT(ITAB(2,2))
      RTAB2 = FLOAT(ITAB(2,1))*FLOAT(ITAB(1,2))
      IF(RTAB1 - RTAB2) 20,10,15
   10 IF(MIN0(ITAB(1,1),ITAB(2,2)).LT. MIN0(ITAB(2,1),ITAB(1,2)))GO TO
     *   20
   15 JX = 2
C                                  COMPUTE PROBABILITY OF THIS OR MORE
C                                  EXTREME TABLE
   20 CALL MDHYP (ITAB(1,JX),ITAB(1,3),ITAB(3,3),ITAB(3,JX),PEQK,PLEK,
     *JER)
      IF (JER .NE. 0) GO TO 30
      P(L)=PLEK
C                                  L=2 FOR THE OPPOSITE EXTREME TABLE
      IF (L.EQ.2) GO TO 9005
C                                  SHIFT TO THE OTHER DIAGONAL
      N = JX
      JX = 3-JX
      P(2) = 0.0
C                                  SET UP THE OPPOSITE EXTREME TABLE
      XK = 2.*FLOAT(ITAB(1,3))*FLOAT(ITAB(3,JX))-
     1     FLOAT(ITAB(3,3))*FLOAT(ITAB(1,JX))
      IF(XK.LT.0.0) GO TO 9005
      K = XK/FLOAT(ITAB(3,3))
      ITAB(1,JX) = K
      ITAB(2,JX) = ITAB(3,JX)-K
      ITAB(1,N) = ITAB(1,3) - K
      ITAB(2,N) = ITAB(3,N)-ITAB(1,N)
      L = 2
      IF(ITAB(2,N))9005,5,5
   25 IER = 129
      GO TO 9000
   30 IER = 130
 9000 CONTINUE
      CALL UERTST(IER,6HNHEXT )
 9005 RETURN
      END
