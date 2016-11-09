C   IMSL ROUTINE NAME   - RLDOPM
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - COEFFICIENT DECODER FOR AN ORTHOGONAL
C                           POLYNOMIAL REGRESSION MODEL
C
C   USAGE               - CALL RLDOPM (C,ID,A,B,T)
C
C   ARGUMENTS    C      - VECTOR OF LENGTH ID+3.
C                         ON INPUT, THE FIRST ID+1 COMPONENTS OF C
C                           CONTAIN THE REGRESSION COEFFICIENTS OF THE
C                           FITTED POLYNOMIAL IN ASCENDING DEGREE ORDER.
C                           C(ID+2) AND C(ID+3) CONTAIN THE SCALING
C                           CONSTANTS A AND B, RESPECTIVELY, USED IN
C                           THE TRANSFORMATION SX = A*X+B.
C                         ON OUTPUT, THE FIRST ID+1 COMPONENTS OF C
C                           CONTAIN THE CONSTANT TERM AND COEFFICIENTS
C                           OF THE DECODED MODEL.
C                ID     - INPUT DEGREE OF THE ORTHOGONAL POLYNOMIAL
C                           MODEL TO BE DECODED.
C                A      - INPUT VECTOR OF LENGTH ID CONTAINING THE
C                           CONSTANTS USED IN GENERATING THE ORTHOGONAL
C                           POLYNOMIALS.
C                B      - INPUT VECTOR OF LENGTH ID CONTAINING
C                           ADDITIONAL CONSTANTS USED IN GENERATING
C                           THE ORTHOGONAL POLYNOMIALS.
C                T      - DOUBLE PRECISION WORK VECTOR OF LENGTH
C                           4*(ID+1).
C
C   REQD. IMSL ROUTINES - NONE REQUIRED
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS      RLDOPM IS DESIGNED TO BE USED IN CONJUNCTION WITH
C                IMSL SUBROUTINES RLFOTH OR RLFOTW. INPUT ARGUMENTS
C                A, B, AND C ARE COMPUTED BY RLFOTH OR RLFOTW.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE RLDOPM (C,ID,A,B,T)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            ID
      REAL               C(1),A(ID),B(ID)
      DOUBLE PRECISION   T(4,1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            ID1,I,K,J,II,II1
      DOUBLE PRECISION   ALPHA,BETA,S,SSS,SS,SKL,SL,TT,ZERO,ONE
      DATA               ZERO/0.0D0/,ONE/1.0D0/
C                                  FIRST EXECUTABLE STATEMENT
      ID1=ID+1
C                                  PERFORM PARTIAL DECODING - COMPUTE
C                                  COEFFICIENTS FOR THE MODEL EXPRESSED
C                                  AS A POLYNOMIAL FUNCTION OF THE
C                                  SCALED INDEPENDENT VARIABLE
      DO 10 I=1,ID1
         T(4,I)=C(I)
         DO 5 K=1,3
            T(K,I)=ZERO
    5    CONTINUE
   10 CONTINUE
      DO 25 I=2,ID1
         T(2,I)=ONE
         K=I-1
         ALPHA=A(K)
         BETA=B(K)
         DO 15 II=2,I
            II1=II-1
            T(3,II)=T(2,II1)-T(2,II)*ALPHA-BETA*T(1,II)
            T(4,II1)=T(4,II1)+C(I)*T(3,II)
   15    CONTINUE
         IF (I .EQ. ID1) GO TO 30
         DO 20 II=1,I
            T(1,II)=T(2,II)
            T(2,II)=T(3,II)
   20    CONTINUE
   25 CONTINUE
C                                  FINISH DECODING - COMPUTE
C                                  COEFFICIENTS FOR THE MODEL EXPRESSED
C                                  AS A POLYNOMIAL FUNCTION OF THE
C                                  ORIGINAL INDEPENDENT VARIABLE
   30 BETA=C(ID+2)
      SSS=ONE
      SS=C(ID+3)
      S=ONE
      SL=ONE
      SKL=ONE
      DO 50 J=1,ID1
         TT=SSS*T(4,J)
         K=J+1
         IF (K .GT. ID1) GO TO 45
         IF (J .EQ. 1) GO TO 35
         SKL=J
         SL=ONE
   35    DO 40 I=K,ID1
            S=SS*S*SKL/SL
            TT=TT+T(4,I)*S
            IF (J .EQ. 1) GO TO 40
            SKL=SKL+ONE
            SL=SL+ONE
   40    CONTINUE
   45    C(J)=TT
         SSS=SSS*BETA
         S=SSS
   50 CONTINUE
      RETURN
      END
