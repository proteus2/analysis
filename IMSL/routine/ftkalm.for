C   IMSL ROUTINE NAME   - FTKALM
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1982
C
C   PURPOSE             - KALMAN FILTERING
C
C   USAGE               - CALL FTKALM(K,X,H,G,Y,S,Q,R,P,IN,IS,IL,N,M1,L,
C                           T1,T2,IT,T3,IER)
C
C   ARGUMENTS    K      - INPUT STEP COUNTER. K=0,1,2,...
C                           WHEN K IS EQUAL TO ZERO, VECTOR X SHOULD
C                           CONTAIN THE PRIOR ESTIMATE OF THE MEAN OF X,
C                           AND THE PROGRAM CALCULATES THE ESTIMATED
C                           VARIANCE OF X AS P=GQG-TRANSPOSE AT STEP 0.
C                X      - INPUT/OUTPUT VECTOR OF LENGTH N. ON INPUT,
C                           X IS THE STATE VECTOR AT STEP K, AND ON
C                           OUTPUT, X CONTAINS THE ESTIMATED STATE
C                           VECTOR AT STEP K+1.
C                H      - INPUT MATRIX OF DIMENSION N BY N. H IS THE
C                           TRANSITION MATRIX AT STEP K.
C                G      - INPUT MATRIX OF DIMENSION N BY L AT STEP K.
C                Y      - INPUT OBSERVATION VECTOR OF LENGTH M1 AT
C                           STEP K+1.
C                S      - INPUT MATRIX OF DIMENSION M1 BY N AT STEP K+1.
C                Q      - INPUT COVARIANCE MATRIX OF DIMENSION L BY L
C                           AT STEP K.
C                R      - INPUT COVARIANCE MATRIX OF DIMENSION
C                           M1 BY M1 AT STEP K+1.
C                P      - INPUT/OUTPUT MATRIX OF DIMENSION N BY N.
C                           ON INPUT, P IS THE VARIANCE MATRIX OF X
C                           AT STEP K.  ON OUTPUT, P IS THE ESTIMATED
C                           VARIANCE MATRIX OF X AT STEP K+1.
C                IN     - INPUT ROW DIMENSION OF THE MATRICES H,G, AND P
C                           EXACTLY AS SPECIFIED IN THE DIMENSION STATE-
C                           MENT IN THE CALLING PROGRAM.
C                IS     - INPUT ROW DIMENSION OF THE MATRICES S AND R
C                           EXACTLY AS SPECIFIED IN THE DIMENSION STATE-
C                           MENT IN THE CALLING PROGRAM.
C                IL     - INPUT ROW DIMENSION OF THE MATRIX Q EXACTLY
C                           AS SPECIFIED IN THE DIMENSION STATEMENT IN
C                           THE CALLING PROGRAM.
C                N      - INPUT SCALAR. SEE DESCRIPTIONS OF X,H,G,M,P.
C                           N MUST BE GREATER THAN 0.
C                M1     - INPUT SCALAR. SEE DESCRIPTIONS OF Y,S,R,T3.
C                           M1 MUST BE GREATER THAN 0.
C                L      - INPUT SCALAR. SEE DESCRIPTIONS OF G,Q.
C                           L MUST BE GREATER THAN 0.
C                T1     - WORK MATRIX OF DIMENSION NM BY NM, WHERE
C                           NM IS THE MAXIMUM OF N AND M1.
C                T2     - WORK MATRIX OF DIMENSION NM BY NML, WHERE
C                           NML IS THE MAXIMUM OF NM AND L.
C                IT     - ROW DIMENSION OF THE MATRICES T1 AND T2
C                           EXACTLY AS SPECIFIED IN THE DIMENSION
C                           STATEMENT IN THE CALLING PROGRAM.
C                T3     - WORK VECTOR OF LENGTH M1.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER=129 INDICATES ONE OF IN, IS, IL, OR IT
C                             IS TOO SMALL, OR THAT ONE OF N, M1,
C                             OR L IS NOT A POSITIVE INTEGER.
C                           IER=130 INDICATES AN ERROR OCCURRED IN
C                             IMSL ROUTINE LEQT1F.
C
C   REQD. IMSL ROUTINES - SINGLE/LEQT1F,LUDATN,LUELMN,UERTST,UGETIO,
C                           VMULFF,VMULFP
C                       - DOUBLE/LEQT1F,LUDATN,LUELMN,UERTST,UGETIO,
C                           VMULFF,VMULFP,VXADD,VXMUL,VXSTO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE FTKALM (K,X,H,G,Y,S,Q,R,P,IN,IS,IL,N,M1,L,T1,T2,IT,
     1                   T3,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            K,IN,IS,IL,N,M1,L,IT,IER
      REAL               X(1),H(IN,1),G(IN,1),Y(1),S(IS,1),Q(IL,1),
     1                   R(IS,1),P(IN,1),T1(IT,1),T2(IT,1),T3(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,I0,I1,J
      DATA               I0/0/,I1/1/
C                                  FIRST EXECUTABLE STATEMENT
      IF (IN.GE.N .AND. IS.GE.M1 .AND. IL.GE.L
     *   .AND. (IT.GE.N .OR. IT.GE.M1) .AND. N.GT.0
     *   .AND. M1.GT.0 .AND. L.GT.0) GO TO 5
      IER = 129
      GO TO 9000
    5 IER = 0
C                                  CALCULATE P IF K = ZERO
      IF (K .NE. 0) GO TO 10
      CALL VMULFF ( G, Q, N, L, L,IN,IL,T2,IT,IER)
      CALL VMULFP (T2, G, N, L, N,IT,IN, P,IN,IER)
C                                  CALCULATE X-PRIME AT STEP K+1
   10 CALL VMULFF ( H, X, N, N,I1,IN,IN,T1,IT,IER)
      DO 15 I = 1,N
         X(I) = T1(I,1)
   15 CONTINUE
C                                  CALCULATE P-PRIME AT STEP K+1
      CALL VMULFF ( H, P, N, N, N,IN,IN,T1,IT,IER)
      CALL VMULFP (T1, H, N, N, N,IT,IN, P,IN,IER)
      CALL VMULFF ( G, Q, N, L, L,IN,IL,T2,IT,IER)
      CALL VMULFP (T2, G, N, L, N,IT,IN,T1,IT,IER)
      DO 25 I = 1,N
         DO 20 J = 1,N
            P(I,J) = P(I,J) + T1(I,J)
   20    CONTINUE
   25 CONTINUE
C                                  CALCULATE MATRIX K AT STEP K+1
      CALL VMULFP ( S, P,M1, N, N,IS,IN,T2,IT,IER)
      CALL VMULFP (T2, S,M1, N,M1,IT,IS,T1,IT,IER)
      DO 35 I = 1,M1
         DO 30 J = 1,M1
            T1(I,J) = T1(I,J) + R(J,I)
   30    CONTINUE
   35 CONTINUE
      CALL LEQT1F(T1, N,M1,IT,T2,I0,T3,IER)
      IF (IER .EQ. 0) GO TO 40
      IER = 130
      GO TO 9000
   40 DO 50 I = 1,M1
         DO 45 J = 1,N
            T1(J,I) = T2(I,J)
   45    CONTINUE
   50 CONTINUE
C                                  CALCULATE X-HAT AT STEP K+1
   55 CALL VMULFF ( S, X,M1, N,I1,IS,IN,T3,IS,IER)
      DO 60 I = 1,M1
         T3(I) = T3(I) - Y(I)
   60 CONTINUE
      CALL VMULFF (T1,T3, N,M1,I1,IT,IS,T2,IT,IER)
      DO 65 I = 1,N
         X(I) = X(I) - T2(I,1)
   65 CONTINUE
C                                  CALCULATE P AT STEP K+1
      CALL VMULFF (T1, S, N,M1, N,IT,IS,T2,IT,IER)
      CALL VMULFF (T2, P, N, N, N,IT,IN,T1,IT,IER)
      DO 75 I = 1,N
         DO 70 J = 1,N
            P(I,J) = P(I,J) - T1(I,J)
   70    CONTINUE
   75 CONTINUE
      GO TO 9005
 9000 CONTINUE
      CALL UERTST (IER,'FTKALM')
 9005 RETURN
      END
