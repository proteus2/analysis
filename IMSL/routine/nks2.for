C   IMSL ROUTINE NAME   - NKS2
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - KOLMOGOROV-SMIRNOV TWO-SAMPLE TEST
C
C   USAGE               - CALL NKS2 (F,N,G,M,PDIF,IER)
C
C   ARGUMENTS    F      - INPUT VECTOR OF LENGTH N CONTAINING SAMPLE
C                           ONE SORTED IN ASCENDING ORDER. F SHOULD BE
C                           THE LARGER (IN TERMS OF THE NUMBER OF
C                           ELEMENTS) OF THE TWO SAMPLES.
C                N      - INPUT NUMBER OF ELEMENTS IN SAMPLE ONE.
C                           N SHOULD BE GREATER THAN OR EQUAL TO M
C                           (SEE DESCRIPTION OF M BELOW).
C                G      - INPUT VECTOR OF LENGTH M CONTAINING SAMPLE
C                           TWO SORTED IN ASCENDING ORDER. G SHOULD BE
C                           THE SMALLER (IN TERMS OF THE NUMBER OF
C                           ELEMENTS) OF THE TWO SAMPLES.
C                M      - INPUT NUMBER OF ELEMENTS IN SAMPLE TWO.
C                           M SHOULD BE LESS THAN OR EQUAL TO N.
C                PDIF   - OUTPUT VECTOR OF LENGTH 6.
C                         PDIF(1) CONTAINS DMN, THE MAXIMUM ABSOLUTE
C                           DIFFERENCE, OVER THE SAMPLES, OF F-G.
C                         PDIF(2) CONTAINS DMN+, THE MAXIMUM DIFFERENCE
C                           OF F-G.
C                         PDIF(3) CONTAINS DMN-, THE MINIMUM DIFFERENCE
C                           OF F-G.
C                         PDIF(4) CONTAINS Z, THE SAMPLE VALUE OF THE
C                           KOLMOGOROV-SMIRNOV STATISTIC, USED TO
C                           EVALUATE PROBABILITIES.
C                         PDIF(5) CONTAINS THE PROBABILITY OF THE
C                           STATISTIC EXCEEDING Z IF THE HYPOTHESIS OF
C                           EQUALITY IS TRUE AND THE ALTERNATIVE IS ONE
C                           SIDED.
C                         PDIF(6) CONTAINS THE PROBABILITY OF THE
C                           STATISTIC EXCEEDING Z IF THE HYPOTHESIS OF
C                           EQUALITY IS TRUE AND THE ALTERNATIVE
C                           IS TWO SIDED.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         WARNING ERROR
C                           IER=33 IMPLIES M/N SMALL, N LESS THAN 50
C                           IER=34 IMPLIES M LESS THAN 80, M/N SMALL
C                           IER=35 IMPLIES N LESS THAN 100 AND A
C                             MULTIPLE OF M
C                           IER=36 IMPLIES N LESS THAN 50, AND NOT A
C                             MULTIPLE OF M
C                           IER=37 IMPLIES N IN (50,100) AND NOT A
C                             MULTIPLE OF M
C
C   REQD. IMSL ROUTINES - MDSMR,UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS  1.  THE USER MAY WISH TO USE ROUTINES IN CHAPTER V OF THE
C                IMSL LIBRARY TO SORT THE OBSERVATIONS IN ASCENDING
C                ORDER.
C            2.  IF ANY OF THE WARNING ERRORS OCCUR IN NKS2, THE USER
C                SHOULD CONSULT THE PROGRAMMING NOTES IN THE MANUAL
C                DOCUMENT FOR FURTHER INFORMATION.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE NKS2   (F,N,G,M,PDIF,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N,M,IER
      REAL               F(N),G(M),PDIF(6)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,J,K,L,NL100,NL50
      REAL               CC,E,FM,FMD,FN,FND,TOTRDS
      DATA               TOTRDS/.6666667/
C                                  FIRST EXECUTABLE STATEMENT
      IER=0
      FN=N
      FM=M
      FND=1./FN
      FMD=1./FM
      PDIF(2)=0.0
      PDIF(3)=0.0
      I=0
      J=0
      K=0
C                                  SMIRNOV CONTINUITY CORRECTION
      CC=SQRT(FND)*TOTRDS
      NL50=1
      NL100=1
      L=0
C                                  CHECK FOR TIES
    5 IF(F(I+1)-G(J+1))10,15,80
   10 K=1
      GO TO 20
   15 K=0
   20 I=I+1
      IF(I.GE.N)GO TO 90
C                                  CHECK FOR TIES
   25 IF(F(I+1).LE.F(I))GO TO 20
   30 IF(K.EQ.0)GO TO 80
   35 E=I*FND-J*FMD
C                                  DMN+
      PDIF(2)=AMAX1(PDIF(2),E)
C                                  DMN-
      PDIF(3)=AMIN1(PDIF(3),E)
      IF(L.EQ.0)GO TO 5
C                                  DMN
   40 PDIF(1)=AMAX1(PDIF(2),ABS(PDIF(3)))
C                                  SET WARNINGS, CHOOSE APPROXIMATION TO
C                                  PROBABILITIES
      IF(N.LT.50)NL50=0
      IF(N.LT.100)NL100=0
      IF(FM*FND.LT..1)GO TO 65
      IF(N.EQ.M*(N/M))GO TO 60
C                                  ADJUST SMIRNOV CONTINUITY CORRECTION
      CC=.6*CC
      IF(NL100.EQ.0)GO TO 55
C                                  CALCULATE THE STATISTIC Z
   45 PDIF(4)=PDIF(1)*SQRT((FN*FM)/(FN+FM))+CC
C                                  OBTAIN PROBABILITIES
   50 CALL MDSMR(PDIF(4),PDIF(5),PDIF(6))
      PDIF(5)=1.-PDIF(5)
      PDIF(6)=1.-PDIF(6)
C                                  CHECK ERROR INDICATOR
      IF(IER)9000,9005,9000
C                                  SET ERROR INDICATORS
   55 IER=37
      IF(NL50.EQ.0)IER=36
      GO TO 45
   60 IF(NL100.EQ.0)IER=35
      GO TO 45
   65 IF(NL50.NE.0)GO TO 75
      IER=33
C                                  DM
   70 PDIF(4)=PDIF(1)*SQRT(FM)
      GO TO 50
   75 IF(M.GE.80)GO TO 95
      IER=34
      GO TO 70
   80 J=J+1
      IF(J.GE.M)GO TO 85
C                                  CHECK FOR TIES
      IF(G(J+1).LE.G(J))GO TO 80
      GO TO 35
   85 L=1
      GO TO 35
   90 L=1
      GO TO 30
   95 IF(NL100.NE.0)GO TO 45
      GO TO 70
 9000 CONTINUE
      CALL UERTST(IER,6HNKS2  )
 9005 RETURN
      END
