C   IMSL ROUTINE NAME   - AMEANS
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - PREPARATION OF A SET OF UNBALANCED DATA FOR
C                           ANALYSIS BY THE METHOD OF UNWEIGHTED MEANS
C
C   USAGE               - CALL AMEANS (Y,N,K,YM,HN,SS,IER)
C
C   ARGUMENTS    Y      - INPUT VECTOR OF LENGTH N(1)+N(2)+...+N(K)
C                           CONTAINING THE OBSERVATIONS.
C                           ALL OBSERVATIONS IN A PARTICULAR CELL
C                           MUST APPEAR CONSECUTIVELY IN VECTOR Y.
C                N      - INPUT VECTOR OF LENGTH K CONTAINING THE
C                           NUMBER OF OBSERVATIONS PER CELL.
C                           N(I) MUST BE GREATER THAN ZERO FOR
C                           I=1,2,...,K.
C                K      - INPUT NUMBER OF CELLS.
C                YM     - OUTPUT VECTOR OF LENGTH K CONTAINING THE CELL
C                           MEANS. THE ORDER OF THE MEANS IN YM
C                           CORRESPONDS TO THE ORDER OF THE OBSERVATIONS
C                           IN Y.
C                HN     - OUTPUT HARMONIC MEAN OF CELL SAMPLE SIZES.
C                SS     - OUTPUT WITHIN CELL ERROR SUM OF SQUARES.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER=129 INDICATES THAT SOME N(I) FOR I=1,2,
C                             ...,K WAS SPECIFIED LESS THAN 1.
C
C   REQD. IMSL ROUTINES - SINGLE/UERTST,UGETIO
C                       - DOUBLE/UERTST,UGETIO,VXADD,VXSTO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE AMEANS (Y,N,K,YM,HN,SS,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N(1),K,IER
      REAL               Y(1),YM(1),HN,SS
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,J,NBEG,NEND
      DOUBLE PRECISION   TEMP,S,FN
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
C                                  CALCULATE MEANS
      DO 5 I = 1,K
         IF (N(I) .LE. 0) GO TO 25
    5 CONTINUE
      NBEG = 1
      HN = 0.0
      SS = 0.0
      NEND = 0
      DO 20 I = 1,K
         NEND = NEND+N(I)
         FN = 1.0D0/N(I)
         TEMP = 0.0D0
         DO 10 J = NBEG,NEND
            TEMP = TEMP+Y(J)
   10    CONTINUE
         TEMP = TEMP*FN
         YM(I) = TEMP
C                                  CALCULATE WITHIN CELL ERROR SUM OF
C                                  SQUARES
         S = 0.0D0
         DO 15 J = NBEG,NEND
            S = Y(J)-TEMP
            SS = SS+S*S
   15    CONTINUE
         NBEG = NEND+1
C                                  CALCULATE HARMONIC MEAN OF CELL
C                                  SAMPLE SIZES
         HN = HN+FN
   20 CONTINUE
      HN = K/HN
      GO TO 9005
C                                  SOME N(I) IS LESS THAN 1
   25 IER = 129
 9000 CONTINUE
      CALL UERTST (IER,'AMEANS'
 9005 RETURN
      END
