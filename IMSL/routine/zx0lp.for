C   IMSL ROUTINE NAME   - ZX0LP
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - SOLVE THE LINEAR PROGRAMMING PROBLEM
C                           (PHASE ONE OR PHASE TWO) VIA THE REVISED
C                           SIMPLEX ALGORITHM
C
C   USAGE               - CALL ZX0LP (IPHASE,C,D,ICOLMS,ROW,K,M,N,ITMAX,
C                           LIC,IR,COPI,IDES,X,WA,IER)
C
C   ARGUMENTS    IPHASE - INPUT OPTION PARAMETER.
C                         IPHASE = 1 IMPLIES SOLVE THE PHASE ONE LINEAR
C                           PROGRAMMING PROBLEM
C                         IPHASE = 2 IMPLIES SOLVE THE PHASE TWO LINEAR
C                           PROGRAMMING PROBLEM
C                C      - THE (M+1) BY N MATRIX C CONTAINS THE FIRST N
C                           COLUMNS OF THE SIMPLEX TABLEAU. (INPUT)
C                           FOR IPHASE=1,
C                           C(1,I)=-C(2,I)-C(3,I)-...-C(M+1,I),I=1,...,N
C                           C(2,I)= THE NEGATIVE OF THE I-TH COST COEFF-
C                           ICIENT.
C                           C(J,I)= THE I-TH COEFFICIENT OF CONSTRAINT
C                                   J-2 TIMES THE SIGN OF THE RIGHT
C                                   HAND SIDE OF CONSTRAINT J-2.
C                                   J=3,4,...,M+1
C                           FOR IPHASE=2,
C                           C(1,I)= THE NEGATIVE OF THE I-TH COST COEFF-
C                           ICIENT
C                           C(J,I)= THE I-TH COEFFICIENT OF CONSTRAINT
C                                   J-1. J=2,3,...,M+1.
C                           MATRIX C IS INPUT. NOTE THAT, FOR IPHASE=1,
C                           THE NUMBER OF CONSTRAINTS IS M-1 (THE COST
C                           EQUATION BEING CONSIDERED THE OTHER
C                           CONSTRAINT)
C                D      - THE VECTOR  D  OF LENGTH (M+1) CONTAINS THE
C                           LAST COLUMN OF THE SIMPLEX TABLEAU. (INPUT)
C                           FOR IPHASE=1,
C                           D(1)=D(2)=0.0
C                           D(J)= THE VALUE OF THE RIGHT HAND
C                           SIDE OF CONSTRAINT J-2.  J=3,4,...,M+1
C                           FOR IPHASE=2,
C                           D(1)=0.0
C                           D(J)= THE VALUE OF THE RIGHT HAND SIDE
C                           OF CONSTRAINT J-1.  J=2,3,...,M+1
C                ICOLMS - INPUT PARAMETERS. ICOLMS, ROW, AND K TO-
C                ROW        GETHER DESCRIBE THE REST OF THE SIMPLEX TAB-
C                K          LEAU (I.E. COLUMNS N+1,N+2,...,N+LIC).
C                           (INPUT) LET THESE COLUMNS BE DENOTED BY IC.
C                           THUS IC IS AN (M+1) BY LIC MATRIX. ALL
C                           ELEMENTS OF IC ARE ZERO EXCEPT THE
C                           FOLLOWING,
C                           IC(K,I)=ROW(I) I=1,...,LIC
C                           IC(ABS(ICOLMS(I)),I)=1.0*SIGN(ICOLMS(I))
C                           I=1,2,...,LIC FOR ABS(ICOLMS(I)).NE.K.
C                           USUALLY, THE USER WILL SET THE INPUT
C                           PARAMETER K EQUAL TO 1 SO THAT THE INPUT
C                           VECTOR ROW SPECIFIES THE FIRST ROW OF IC.
C                           HOWEVER, IN SOME SPECIAL CASES THE
C                           USER MAY WISH TO SPECIFY ANOTHER
C                           ROW.  THE INPUT VECTOR ICOLMS OF
C                           LENGTH LIC SPECIFIES WHICH ELEMENT IN EACH
C                           COLUMN OF IC (EXCLUSIVE OF ROW K ) IS PLUS
C                           OR MINUS ONE. ICOLMS SPECIFIES THE PLACEMENT
C                           OF COST, SLACK, AND ARTIFICIAL VARIABLES.
C                           THE COMPLETE SIMPLEX TABLEAU IS THUS GIVEN
C                           BY THE CONCATENATION OF C, IC, AND D.
C                M      - FOR IPHASE=1, THE INPUT PARAMETER M EQUALS THE
C                           NUMBER OF CONSTRAINTS PLUS ONE (THE COST
C                           EQUATION). FOR IPHASE=2, M EQUALS THE NUMBER
C                           OF CONSTRAINTS.
C                N      - THE INPUT PARAMETER N SPECIFIES THE NUMBER OF
C                           VARIABLES (NOT INCLUDING COST, SLACK OR
C                           ARTIFICIAL VARIABLES)
C                ITMAX  - MAXIMUM ALLOWABLE NUMBER OF ITERATIONS OF THE
C                           MODIFIED SIMPLEX ALGORITHM. CHARACTERISTIC-
C                           ALLY THIS ALGORITHM REQUIRES ABOUT 5*N
C                           ITERATIONS TO CONVERGE TO THE OPTIMAL SOLU-
C                           TION. IT IS NOT DIFFICULT, HOWEVER, TO
C                           CONSTRUCT EXAMPLES WHICH REQUIRE MANY MORE
C                           ITERATIONS. (INPUT)
C                LIC    - THE LENGTH OF INPUT VECTORS ROW AND ICOLMS.
C                           (INPUT)
C                           FOR IPHASE=1, LIC IS THE NUMBER OF SLACK
C                           VARIABLES PLUS THE NUMBER OF ARTIFICIAL
C                           VARIABLES PLUS TWO (PHASE ONE COST AND PHASE
C                           TWO COST).
C                           FOR IPHASE=2, LIC IS THE NUMBER OF SLACK
C                           VARIABLES PLUS ONE (PHASE TWO COST).
C                IR     - ROW DIMENSION OF THE MATRICES C AND COPI
C                           EXACTLY AS SPECIFIED IN THE DIMENSION
C                           STATEMENT IN THE CALLING PROGRAM. (INPUT)
C                COPI   - AN (M+1) BY (M+1) INPUT MATRIX WHICH IS THE
C                           INVERSE OF THE STARTING BASIS. THE STARTING
C                           BASIS CONSISTS OF (M+1) COLUMNS FROM (C,IC)
C                           AND IS FURTHER DESCRIBED BY PARAMETER IDES.
C                           USUALLY ALL COLUMNS OF THE STARTING BASIS
C                           ARE CONTAINED IN IC. ALSO, USUALLY ALL COL-
C                           UMNS OF THE STARTING BASIS ARE THEMSELVES
C                           COLUMNS OF THE (M+1) BY (M+1) IDENTITY MAT-
C                           RIX. THUS, USUALLY, COPI WILL BE THE IDEN-
C                           TITY MATRIX.  IN SOME APPLICATIONS HOWEVER,
C                           IT IS CONVENIENT TO USE OTHER STARTING
C                           BASES.  THE SUBROUTINE UPDATES
C                           COPI DURING EACH ITERATION AS THE
C                           TRIAL BASIS CHANGES.
C                IDES   - AN INPUT VECTOR OF LENGTH (M+1) WHICH CONTAINS
C                           THE STARTING BASIS DESCRIPTION. IF COLUMN
C                           IDES(I) OF (C,IC) IS PREMULTIPLIED BY COPI
C                           THE RESULT SHOULD BE THE I-TH COLUMN OF THE
C                           (M+1) BY (M+1) IDENTITY).  THE SUBROUTINE
C                           UPDATES IDES DURING EACH ITERATION AS THE
C                           TRIAL BASIS CHANGES.
C                           FOR IPHASE=1, IDES(1) SHOULD POINT TO THE
C                           PHASE ONE COST VARIABLE AND IDES(2) SHOULD
C                           POINT TO THE PHASE TWO COST VARIABLE.
C                           FOR IPHASE=2, IDES(1) SHOULD POINT TO THE
C                           (PHASE TWO) COST VARIABLE.
C                           FOR IPHASE=1, THE ROUTINE WILL LEAVE THE
C                           STARTING IDES(1) AND IDES(2) IN THE BASIS AT
C                           AT ALL TIMES.
C                           FOR IPHASE=2, THE ROUTINE WILL LEAVE THE
C                           STARTING IDES(1) IN THE BASIS AT ALL TIMES.
C                X      - AN OUTPUT VECTOR OF LENGTH (M+1) WHICH
C                           CONTAINS THE SOLUTION. THERE ARE N+LIC
C                           VARIABLES IN THE SIMPLEX TABLEAU. THE FIRST
C                           N VARIABLES ARE THE USERS ORIGINAL VARIABLES
C                           THE (N+1) VARIABLE IS THE (PHASE TWO) COST.
C                           THE NEXT FEW VARIABLES ARE THE SLACK
C                           VARIABLES (IF ANY). FOR IPHASE=2, THERE ARE
C                           NO MORE VARIABLES.  FOR IPHASE=1, THE SLACK
C                           VARIABLES (IF ANY) ARE FOLLOWED BY THE PHASE
C                           ONE (ARTIFICIAL) COST. THIS IN TURN IS
C                           FOLLOWED BY ONE OR MORE ARTIFICIAL VARIABLES
C                           IN THE FINAL SOLUTION. ALL THE VARIABLES ARE
C                           ZERO EXCEPT VARIABLE(IDES(I))=X(I), I=1,
C                           ...,M+1. THE OUTPUT VECTOR X ALSO HAS THE
C                           PROPERTY THAT IT IS EQUAL TO THE INPUT
C                           VECTOR D PRE-MULTIPLIED BY THE OUTPUT
C                           MATRIX COPI.
C                WA     - WORK VECTOR OF LENGTH (M+1)
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER = 129 INDICATES ORIGINAL RIGHT HAND SIDE
C                             CONTAINS NEGATIVE ELEMENTS
C                           IER = 130 INDICATES STARTING SOLUTION NOT
C                             FEASIBLE
C                           IER = 131 INDICATES COST CRITERION HAS
C                             UNBOUNDED VALUES
C                           IER = 132 INDICATES THE MAXIMUM NUMBER OF
C                             ITERATIONS WAS REACHED
C
C   REQD. IMSL ROUTINES - UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE ZX0LP  (IPHASE,C,D,ICOLMS,ROW,K,M,N,ITMAX,LIC,IR,COPI,
     1                   IDES,X,WA,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            IPHASE,ICOLMS(1),K,M,N,ITMAX,LIC,IR,IDES(1),IER
      REAL               C(IR,1),D(1),ROW(1),COPI(IR,1),X(1),WA(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,ISW,MP1,J,L,IP,ITER,IQ,II,IN
      REAL               R1,R2,D1,ZERO,ONE,EPS,D2,D3,T,RMT
      DATA               EPS/1.0E-5/
C                                  EPS IS USED IN TESTS FOR ZERO
C                                    IF ABS(T) .LE. EPS, THEN T IS
C                                    CONSIDERED TO BE ZERO
      DATA               ZERO/0.0/,ONE/1.0/
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
      IF (IPHASE .EQ. 1) ISW = 3
      IF (IPHASE .NE. 1) ISW = 2
C                                  MAKE SURE RHS IS NON-NEGATIVE
      ITER = 0
      MP1 = M+1
      RMT = MP1*50
      DO 10 I=ISW,MP1
         IF (D(I).GE.ZERO) GO TO 10
         IER = 129
         GO TO 9000
   10 CONTINUE
C                                  MAKE SURE STARTING SOLN. IS POSITIVE
      DO 20 I=1,MP1
         D1 = ZERO
         D2 = ZERO
         DO 15 J=1,MP1
            T = COPI(I,J)*D(J)
            D1 = D1+T
            D2 = D2+ABS(T)
   15    CONTINUE
         X(I) = D1
         D2 = D2*RMT
         IF (D2+ABS(D1) .EQ. D2) D1 = ZERO
         IF (I.LT.ISW) GO TO 20
         IF (D1.LT.-EPS) IER = 130
   20 CONTINUE
      IF (IER.NE.0) GO TO 9000
C                                  FIND NEXT PIVOT COLUMN (IQ)
   25 IQ = 0
      D3 = ZERO
      DO 26 I=1,MP1
   26 D3 = D3+ABS(COPI(1,I))
      R1 = -EPS
      DO 35 J=1,N
         D1 = ZERO
         D2 = ZERO
         DO 30 I=1,MP1
            T = COPI(1,I)*C(I,J)
            D1 = D1+T
            D2 = D2+ABS(T)
   30    CONTINUE
         D2 = D2*RMT
         IF (D2+ABS(D1) .EQ. D2) GO TO 35
         IF (D1.GE.R1) GO TO 35
         R1 = D1
         IQ = J
   35 CONTINUE
      DO 40 L=1,LIC
         II = ICOLMS(L)
         J = IABS(II)
         D1 = COPI(1,J)
         IF (II.LT.0) D1 = -D1
         IF (J.EQ.K) D1 = ZERO
         D2 = ONE
         IF (J.EQ.K) D2 = ZERO
         D1 = D1+COPI(1,K)*ROW(L)
         D2 = D2+ABS(ROW(L))
         D2 = D2*D3*RMT
         IF (D2+ABS(D1) .EQ. D2) GO TO 40
         IF (D1.GE.R1) GO TO 40
         R1 = D1
         IQ = L+N
   40 CONTINUE
C                                  FIND NEXT PIVOT ROW (IP)
      R1 = -ONE
      IP = 0
      DO 65 I=1,MP1
         D1 = ZERO
         D3 = ZERO
         DO 45 J=1,MP1
            D1 = D1+COPI(I,J)*D(J)
            D3 = D3+ABS(COPI(I,J))
   45    CONTINUE
         X(I) = D1
         IF (IQ.GT.N) GO TO 55
         IF (IQ.EQ.0) GO TO 65
         D1 = ZERO
         D2 = ZERO
         DO 50 J=1,MP1
            T = COPI(I,J)*C(J,IQ)
            D1 = D1+T
            D2 = D2+ABS(T)
   50    CONTINUE
         GO TO 60
   55    IN = IQ-N
         II = ICOLMS(IN)
         J = IABS(II)
         D1 = COPI(I,J)
         IF (II.LT.0) D1 = -D1
         IF (J.EQ.K) D1 = ZERO
         D2 = ONE
         IF (J.EQ.K) D2 = ZERO
         D1 = D1+COPI(I,K)*ROW(IN)
         D2 = D2+ABS(ROW(IN))
         D2 = D2*D3
   60    WA(I) = D1
         D2 = D2*RMT
         IF (D2+ABS(D1) .EQ. D2) GO TO 65
         IF (D1.LE.EPS) GO TO 65
         IF (I.LT.ISW) GO TO 65
         R2 = X(I)/D1
         IF (R2.GE.R1.AND.R1.GE.ZERO) GO TO 65
         R1 = R2
         IF (R1.LT.ZERO) R1 = ZERO
         IP = I
   65 CONTINUE
C                                  DETERMINE IF PRESENT SOLN. IS OPTIMAL
      IF (IQ.EQ.0) GO TO 9005
      DO 70 J=1,MP1
         IF (IDES(J).EQ.IQ) GO TO 9005
   70 CONTINUE
      IF (ITER.GE.ITMAX) GO TO 95
      IF (IP.NE.0) GO TO 75
C                                  UNBOUNDED FEASIBLE SOLUTION
C                                    PSOL(IDES(I)) = X(I)-TH*WA(I)
C                                    PSOL(IQ) = TH FOR TH .GT. 0.0
C                                    THE OBJECTIVE FUNCTION TENDS TO
C                                    INFINITY AS TH TENDS TO INFINITY
      IER = 131
      GO TO 9000
C                                  FORM E(IP,IQ)*COPI
   75 D1 = WA(IP)
      DO 80 J=1,MP1
         COPI(IP,J) = COPI(IP,J)/D1
   80 CONTINUE
      DO 90 I=1,MP1
         IF (I.EQ.IP) GO TO 90
         D1 = WA(I)
         DO 85 J=1,MP1
            COPI(I,J) = COPI(I,J)-COPI(IP,J)*D1
   85    CONTINUE
   90 CONTINUE
      IDES(IP) = IQ
      ITER = ITER+1
      GO TO 25
   95 IER = 132
 9000 CONTINUE
      CALL UERTST (IER,6HZX0LP )
 9005 RETURN
      END
