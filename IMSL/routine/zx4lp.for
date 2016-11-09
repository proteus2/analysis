C   IMSL ROUTINE NAME   - ZX4LP
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1980
C
C   PURPOSE             - SOLVE THE LINEAR PROGRAMMING PROBLEM VIA THE
C                           REVISED SIMPLEX ALGORITHM (ALTERNATE EASY TO
C                           USE VERSION)
C
C   USAGE               - CALL ZX4LP (A,IA,B,C,N,M1,M2,S,PSOL,DSOL,RW,
C                           IW,IER)
C
C   ARGUMENTS    A      - MATRIX OF DIMENSION M1+M2+2 BY N+M1+2
C                           CONTAINING THE COEFFICIENTS OF THE M1
C                           INEQUALITY CONSTRAINTS IN THE FIRST M1
C                           ROWS FOLLOWED BY THE COEFFICIENTS OF
C                           THE M2 EQUALITY CONSTRAINTS. (INPUT)
C                           THE LAST 2 ROWS AND THE LAST M1+2
C                           COLUMNS OF A ARE USED ONLY AS WORKING
C                           STORAGE.
C                IA     - ROW DIMENSION OF MATRIX A EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT
C                           IN THE CALLING PROGRAM. (INPUT) IA
C                           MUST BE AT LEAST M1+M2+2.
C                B      - VECTOR OF LENGTH M1+M2 CONTAINING THE
C                           RIGHT HAND SIDES OF THE INEQUALITY
C                           CONSTRAINTS IN ITS FIRST M1 LOCATIONS
C                           FOLLOWED BY THE M2 RIGHT HAND SIDES
C                           OF THE EQUALITY CONSTRAINTS. (INPUT)
C                C      - VECTOR OF LENGTH N CONTAINING THE
C                           COEFFICIENTS OF THE OBJECTIVE FUNCTION.
C                           (INPUT)
C                N      - NUMBER OF UNKNOWNS IN THE MODEL. (INPUT)
C                M1     - NUMBER OF INEQUALITY CONSTRAINTS. (INPUT)
C                M2     - NUMBER OF EQUALITY CONSTRAINTS. (INPUT)
C                S      - VALUE OF THE OBJECTIVE FUNCTION. (OUTPUT)
C                PSOL   - VECTOR OF LENGTH N CONTAINING THE PRIMAL
C                           SOLUTION. (OUTPUT)
C                DSOL   - VECTOR OF LENGTH M1+M2 CONTAINING THE
C                           DUAL SOLUTION. (OUTPUT) SEE REMARKS.
C                RW     - REAL WORK VECTOR OF LENGTH
C                           (M1+M2+2)*(M1+M2+4)+2*(N+M1).
C                IW     - INTEGER WORK VECTOR OF LENGTH
C                           2*(N+M1)+3.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER = 129 INDICATES THAT IA IS LESS
C                             THAN M1+M2+2.
C                           IER = 130 INDICATES THAT M2 IS GREATER
C                             THAN N.
C                           IER = 131 INDICATES THAT AN EXCESSIVE
C                             NUMBER OF ITERATIONS WERE DONE.
C                           IER = 132 INDICATES THAT REDUNDANCIES
C                             ARE PRESENT IN THE CONSTRAINTS.
C                           IER = 133 INDICATES THAT B=0.
C                           IER = 134 INDICATES THAT THE OBJECTIVE
C                             FUNCTION IS UNBOUNDED.
C                           IER = 135 INDICATES THAT THE CONSTRAINTS
C                             ARE INFEASIBLE.
C                           IER = 136 INDICATES THAT C(TRANSPOSE)*PSOL
C                             IS NOT EQUAL TO B(TRANSPOSE)*DSOL, OR
C                             PSOL OR DSOL DO NOT SATISFY THE
C                             CONSTRAINTS.
C
C   REQD. IMSL ROUTINES - SINGLE/UERTST,UGETIO,VBLA=ISAMAX,VBLA=SASUM,
C                           VBLA=SAXPY,VBLA=SCOPY,VBLA=SDOT,VBLA=SROTM,
C                           VBLA=SROTMG,VBLA=SSCAL,VBLA=SSWAP,ZX4LQ,
C                           ZX4LR
C                       - DOUBLE/UERTST,UGETIO,VBLA=IDAMAX,VBLA=DASUM,
C                           VBLA=DAXPY,VBLA=DCOPY,VBLA=DDOT,VBLA=DROTM,
C                           VBLA=DROTMG,VBLA=DSCAL,VBLA=DSWAP,ZX4LQ,
C                           ZX4LR
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE ZX4LP (A,IA,B,C,N,M1,M2,S,PSOL,DSOL,RW,IW,IER)
C
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            IA,N,M1,M2,IW(1),IER
      REAL               A(IA,1),B(1),C(1),S,PSOL(1),DSOL(1),RW(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            L1,L2,L3,MLP,MODE,NLP
      REAL               ASUM,EPS,PRGOPT(1),REPS,SD,SUMP
C                                  SEPS=MACHINE PRECISION
      DATA REPS/Z3C100000/
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
      EPS = SQRT(REPS)
      NLP = N+M1
      MLP = M1+M2
      IF (M1.EQ.0) GO TO 15
      DO 10 J=1,M1
         DO 5 I=1,MLP
            A(I,J+N) = 0.0
    5    CONTINUE
         A(J,J+N) = 1.0
         A(MLP+1,J+N) = 0.0
   10 CONTINUE
   15 DO 20 J=1,N
         A(MLP+1,J) = -C(J)
   20 CONTINUE
      DO 25 I=1,MLP
         A(I,NLP+1) = B(I)
   25 CONTINUE
      A(MLP+1,NLP+1) = 0.0
      DO 30 I=1,NLP
         IW(I) = 0
   30 CONTINUE
      PRGOPT(1) = 1
      L1 = MLP+3
      L2 = L1+NLP
      L3 = L2+MLP
C                                  CALL ZX4LR TO SOLVE LP PROBLEM
C
      CALL ZX4LR (A,IA,MLP,NLP,PRGOPT,RW,DSOL,MODE,RW(NLP+1),RW(NLP
     1+L1),RW(NLP+L2),RW(NLP+L3),MLP+2,IW,IW(NLP+3),IER)
      DO 35 I=1,N
         PSOL(I) = RW(I)
   35 CONTINUE
      DO 40 J=1,MLP
         DSOL(J) = -DSOL(J)
   40 CONTINUE
C                                  CHECK TO SEE IF PSOL AND DSOL
C                                  SATISFY THE CONSTRAINTS AND CHECK IF
C                                  C(TRANSPOSE)*PSOL=B(TRANSPOSE)*DSOL
      S = SDOT(N,PSOL,1,C,1)
      IF (IER.GE.128) GO TO 9005
      SD = SDOT(MLP,DSOL,1,B,1)
      ASUM = ABS(S)
      DO 45 I=1,MLP
         ASUM = ASUM+ABS(B(I))*ABS(DSOL(I))
   45 CONTINUE
      SUMP = ABS(S-SD)
      IF (ASUM.NE.0.0) SUMP = SUMP/ASUM
      IF (SUMP.LE.EPS) GO TO 50
      IER = 136
      GO TO 9000
   50 DO 60 I=1,MLP
         ASUM = ABS(B(I))
         DO 55 J=1,N
            ASUM = ASUM+ABS(A(I,J))*ABS(PSOL(J))
   55    CONTINUE
         SUMP = SDOT(N,A(I,1),IA,PSOL,1)-B(I)
         IF (I.GT.M1) SUMP = ABS(SUMP)
         IF (ASUM.NE.0.0) SUMP = SUMP/ASUM
         IF (SUMP.LE.EPS) GO TO 60
         IER = 136
         GO TO 9000
   60 CONTINUE
      DO 70 J=1,N
         ASUM = ABS(C(J))
         DO 65 I=1,MLP
            ASUM = ASUM+ABS(A(I,J))*ABS(DSOL(I))
   65    CONTINUE
         SUMP = C(J)-SDOT(MLP,A(1,J),1,DSOL,1)
         IF (ASUM.NE.0.0) SUMP = SUMP/ASUM
         IF (SUMP.LE.EPS) GO TO 70
         IER = 136
         GO TO 9000
   70 CONTINUE
      IER = 0
 9000 IF (IER.EQ.0) GO TO 9005
      CALL UERTST (IER,6HZX4LP )
 9005 RETURN
      END
