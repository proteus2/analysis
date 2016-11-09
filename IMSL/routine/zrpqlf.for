C   IMSL ROUTINE NAME   - ZRPQLF
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - NUCLEUS CALLED ONLY BY IMSL SUBROUTINE
C                           ZRPOLY
C
C   REQD. IMSL ROUTINES - NONE REQUIRED
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE ZRPQLF (ITYPE)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            ITYPE
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            N,NN,I
      REAL               ARE,ETA,RMRE
      REAL               P(101),QP(101),RK(101),QK(101),SVK(101)
      REAL               SR,SI,U,V,RA,RB,C,D,A1,A2,A3,
     1                   A6,A7,E,F,G,H,SZR,SZI,RLZR,RLZI,TEMP,ZERO
      COMMON /ZRPQLJ/    P,QP,RK,QK,SVK,SR,SI,U,V,RA,RB,C,D,A1,A2,A3,A6,
     1                   A7,E,F,G,H,SZR,SZI,RLZR,RLZI,ETA,ARE,RMRE,N,NN
      DATA               ZERO/0.0/
C                                  COMPUTES THE NEXT K POLYNOMIALS
C                                    USING SCALARS COMPUTED IN ZRPQLE
C                                  FIRST EXECUTABLE STATEMENT
      IF (ITYPE.EQ.3) GO TO 20
      TEMP = RA
      IF (ITYPE.EQ.1) TEMP = RB
      IF (ABS(A1).GT.ABS(TEMP)*ETA*10.) GO TO 10
C                                  IF A1 IS NEARLY ZERO THEN USE A
C                                    SPECIAL FORM OF THE RECURRENCE
      RK(1) = ZERO
      RK(2) = -A7*QP(1)
      DO 5 I=3,N
    5 RK(I) = A3*QK(I-2)-A7*QP(I-1)
      RETURN
C                                  USE SCALED FORM OF THE RECURRENCE
   10 A7 = A7/A1
      A3 = A3/A1
      RK(1) = QP(1)
      RK(2) = QP(2)-A7*QP(1)
      DO 15 I=3,N
   15 RK(I) = A3*QK(I-2)-A7*QP(I-1)+QP(I)
      RETURN
C                                  USE UNSCALED FORM OF THE RECURRENCE
C                                    IF TYPE IS 3
   20 RK(1) = ZERO
      RK(2) = ZERO
      DO 25 I=3,N
   25 RK(I) = QK(I-2)
      RETURN
      END
