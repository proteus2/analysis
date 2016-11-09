C   IMSL ROUTINE NAME   - ANESTE
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - ANALYSIS OF COMPLETELY NESTED DESIGN DATA
C                           WITH EQUAL NUMBERS IN THE SUBCLASSES
C
C   USAGE               - CALL ANESTE (NF,NL,Y,S,NDF,IER)
C
C   ARGUMENTS    NF     - INPUT NUMBER OF FACTORS (NUMBER OF DISTINCT
C                           SUBSCRIPTS IN THE MODEL). NF MUST GREATER
C                           THAN OR EQUAL TO TWO.
C                NL     - INPUT VECTOR OF LENGTH NF CONTAINING THE
C                           NUMBER OF LEVELS FOR EACH OF THE NF FACTORS.
C                           NL(I) MUST BE GREATER THAN OR EQUAL TO
C                           TWO FOR I=1,2,...,NF.
C                Y      - INPUT VECTOR OF LENGTH 2**NF+(NL(1)+1)*
C                           (NL(2)+1)*...*(NL(NF)+1).
C                           THE FIRST 2**NF LOCATIONS ARE WORK STORAGE.
C                           Y CONTAINS THE RESPONSES IN LOCATIONS
C                           2**NF+1, 2*NF+2, ... ,
C                           2**NF+NL(1)*NL(2)*...*NL(NF).
C                           THE REMAINING COMPONENTS ARE WORK STORAGE.
C                S      - OUTPUT VECTOR OF LENGTH 2**NF. THE FIRST NF
C                           LOCATIONS CONTAIN THE SUM OF SQUARES
C                           FOR EACH FACTOR NESTED IN THE NEXT
C                           HIGHER LEVEL FACTOR. THE NF+1 LOCATION
C                           CONTAINS THE CORRECTED TOTAL SUM OF SQUARES.
C                           THE REMAINING COMPONENTS OF S ARE
C                           WORK STORAGE.
C                NDF    - OUTPUT VECTOR OF LENGTH 2**NF.
C                           THE FIRST NF LOCATIONS CONTAIN THE
C                           DEGREES OF FREEDOM CORRESPONDING
C                           TO COMPONENTS OF THE S VECTOR. THE NF+1
C                           LOCATION CONTAINS THE CORRECTED TOTAL
C                           DEGREES OF FREEDOM. THE REMAINING COMPONENTS
C                           OF NDF ARE WORK STORAGE.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER=129 INDICATES THAT NF WAS SPECIFIED
C                             LESS THAN 2.
C                           IER=130 INDICATES THAT SOME NL(I) WAS
C                             SPECIFIED LESS THAN TWO FOR I=1,2,...,NF.
C
C   REQD. IMSL ROUTINES - SINGLE/AFACN,UERTST,UGETIO
C                       - DOUBLE/AFACN,UERTST,UGETIO,VXADD,VXMUL,
C                           VXSTO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE ANESTE (NF,NL,Y,S,NDF,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            NF,NL(NF),NDF(1),IER
      REAL               Y(1),S(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,INC,INX,IOPT,J,L,NFM1,NSUM
      REAL               SUM
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
C                                  CALL AFACN TO COMPUTE SUM
C                                    OF SQUARES AND DEGREES OF
C                                    FREEDOM FOR EACH FACTOR
      IOPT = 1
      CALL AFACN (IOPT,NF,NL,Y,S,NDF,IER)
C                                  CHECK ERROR INDICATOR
      IF (IER .EQ. 0) GO TO 5
      GO TO 9000
C                                  SET STARTING LOCATION AND DISTANCE
C                                  BETWEEN ITEMS
    5 INX = 1
      INC = 2
      L = 2**NF
C                                  SUM THE SUMS OF SQUARES AND
C                                  THE DEGREES OF FREEDOM
      DO  15  I=1,NF
         SUM = 0.0
         NSUM = 0
         DO 10  J=INX,L,INC
            SUM = SUM + S(J)
            NSUM = NSUM + NDF(J)
   10    CONTINUE
         S(INX) = SUM
         NDF(INX) = NSUM
         INX = INC
         INC = INC * 2
   15 CONTINUE
C                                  MOVE THE OUTPUT INTO SEQUENTIAL
C                                  LOCATIONS
      NFM1 = NF-1
      DO  20  I=2,NFM1
         J = 2**I
         S(I+1) = S(J)
         NDF(I+1) = NDF(J)
   20 CONTINUE
      S(NF+1) = S(L)
      NDF(NF+1) = NDF(L)
      GO TO 9005
 9000 CONTINUE
      CALL UERTST (IER,'ANESTE')
 9005 RETURN
      END
