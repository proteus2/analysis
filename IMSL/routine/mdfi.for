C   IMSL ROUTINE NAME   - MDFI
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1982
C
C   PURPOSE             - INVERSE F PROBABILITY DISTRIBUTION FUNCTION
C
C   USAGE               - CALL MDFI (P,D1,D2,X,IER)
C
C   ARGUMENTS    P      - INPUT PROBABILITY IN THE EXCLUSIVE RANGE (0,1)
C                D1     - INPUT NUMERATOR DEGREES OF FREEDOM OF THE F
C                           DISTRIBUTION
C                D2     - INPUT DENOMINATOR DEGREES OF FREEDOM OF THE F
C                           DISTRIBUTION
C                X      - OUTPUT VALUE SUCH THAT THE PROBABILITY THAT
C                           A RANDOM VARIABLE DISTRIBUTED F(D1,D2) IS
C                           LESS THAN OR EQUAL TO X IS P.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER = 129 MEANS AN ERROR OCCURRED IN MDBETI
C                           IER = 130 MEANS P WAS NOT IN THE EXCLUSIVE
C                             RANGE (0,1)
C                         WARNING ERROR
C                           IER = 35 MEANS THE COMPUTED VALUE OF X
C                             WOULD HAVE PRODUCED AN OVERFLOW. X IS SET
C                             TO MACHINE INFINITY.
C
C   REQD. IMSL ROUTINES - H32/MDBETA,MDBETI,MLGAMD=DLGAMA,UERTST,UGETIO
C                       - H36,H48,H60/MDBETA,MDBETI,MLGAMA=ALGAMA,
C                           UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE MDFI(P,D1,D2,X,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            IER
      REAL               P,D1,D2,X
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      REAL               A,B,RINFP,RDELP
      DATA               RINFP/Z7FFFFFFF/
      DATA               RDELP/Z40FFFFFF/
C                                  FIRST EXECUTABLE STATEMENT
      IF (P.LE.0..OR.P.GE.1.) GO TO 10
      IER = 0
C                                  REPARAMETERIZE FOR BETA USAGE
      A = .5*D1
      B = .5*D2
      IF (P.GT.0.5) GO TO 5
      CALL MDBETI (P,A,B,X,IER)
      IF (IER.NE.0) GO TO 15
      IF (X.GE.RDELP) GO TO 20
C                                  MODIFY BETA OUTPUT FOR X FROM
C                                  F DISTRIBUTION
      X = D2*X/(D1*(1.D0-DBLE(X)))
      GO TO 9005
    5 PP = 1.0-P
      CALL MDBETI (PP,B,A,X,IER)
      IF (IER.NE.0) GO TO 15
      IF (X.EQ.0.0) GO TO 20
      X = (1.0/X-1.0)*D2/D1
      GO TO 9005
   10 IER = 130
      GO TO 9000
   15 IER = 129
      GO TO 9000
   20 IER = 35
      X = RINFP
 9000 CONTINUE
      CALL UERTST (IER,6HMDFI  )
 9005 RETURN
      END
