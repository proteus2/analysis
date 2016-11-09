C   IMSL ROUTINE NAME   - VDCPS
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - DECOMPOSE AN INTEGER INTO ITS PRIME FACTORS
C
C   USAGE               - CALL VDCPS (N,NPF,IPF,IEXP,IPWR)
C
C   ARGUMENTS    N      - INTEGER NUMBER TO BE FACTORED. (INPUT)
C                NPF    - INTEGER NUMBER OF DIFFERENT PRIME FACTORS OF
C                           THE ABSOLUTE VALUE OF N. (OUTPUT)
C                           IF N IS -1, 0, OR 1, NPF IS SET TO 0.
C                IPF    - INTEGER VECTOR OF LENGTH 13. (OUTPUT) IPF(I)
C                           CONTAINS THE PRIME FACTORS OF THE ABSOLUTE
C                           VALUE OF N, FOR I=1,2,...,NPF. THE
C                           REMAINING 13-NPF LOCATIONS OF IPF ARE NOT
C                           USED.
C                IEXP   - INTEGER VECTOR OF LENGTH 13. (OUTPUT) IEXP(I)
C                           CONTAINS THE EXPONENT OF IPF(I), FOR
C                           I=1,2,...,NPF. THE REMAINING 13-NPF
C                           LOCATIONS OF IEXP ARE NOT USED.
C                IPWR   - INTEGER VECTOR OF LENGTH 13. (OUTPUT) IPWR(I)
C                           CONTAINS THE QUANTITY IPF(I)**IEXP(I), FOR
C                           I=1,2,...,NPF. THE REMAINING 13-NPF
C                           LOCATIONS OF IPWR ARE NOT USED.
C
C   REQD. IMSL ROUTINES - NONE REQUIRED
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS      THE OUTPUT FROM VDCPS SHOULD BE INTERPRETED IN THE
C                FOLLOWING MANNER:
C                ABS(N) = IPF(1)**IEXP(1) * IPF(2)**IEXP(2) * ...
C                * IPF(NPF)**IEXP(NPF)
C
C-----------------------------------------------------------------------
C
      SUBROUTINE VDCPS  (N,NPF,IPF,IEXP,IPWR)
C
      DIMENSION          IPF(13),IEXP(13),IPWR(13)
C                                  FIRST EXECUTABLE STATEMENT
      I = 0
      NF = IABS(N)
C                                  IF N IS -1, 0, OR 1, SET NPF TO 0
      IF (NF .LE. 1) GO TO 30
      IFC = 0
      ID = 2
    5 IQ = NF/ID
      IF (NF .NE. ID*IQ) GO TO 15
C                                  ID IS A FACTOR OF NF
      NF = IQ
      IF (ID .LE. IFC) GO TO 10
C                                  ID IS A NEW DIVISOR
      I = I+1
C                                  UPDATE THE PRIME FACTOR VECTOR,
C                                    THE POWER VECTOR, AND SET THE
C                                    EXPONENT TO 1
      IPF(I) = ID
      IPWR(I) = ID
      IEXP(I) = 1
      IFC = ID
      GO TO 5
C                                  ID FACTORS N MORE THAN ONCE. UPDATE
C                                    THE POWER VECTOR AND THE EXPONENT
C                                    VECTOR.
   10 IPWR(I) = ID*IPWR(I)
      IEXP(I) = IEXP(I)+1
      GO TO 5
C                                  ID IS NOT A FACTOR OF NF. IF NF IS
C                                    NOT PRIME, CONTINUE TO FACTOR.
   15 IF (IQ .LE. ID) GO TO 25
C                                  UPDATE THE DIVISOR
      IF (ID .GT. 2) GO TO 20
      ID = 3
      GO TO 5
   20 ID = ID+2
      GO TO 5
C                                  NF IS PRIME OR 1
   25 IF (NF .LE. 1) GO TO 30
C                                  NF IS PRIME. UPDATE THE PRIME
C                                    FACTOR VECTOR, THE POWER VECTOR,
C                                    AND THE EXPONENT VECTOR.
      I = I+1
      IPF(I) = NF
      IPWR(I) = NF
      IEXP(I) = 1
   30 NPF = I
      RETURN
      END
