C   IMSL ROUTINE NAME   - ZCPOLY
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - ZEROS OF A POLYNOMIAL WITH COMPLEX
C                           COEFFICIENTS (JENKINS-TRAUB)
C
C   USAGE               - CALL ZCPOLY (A,NDEG,Z,IER)
C
C   ARGUMENTS    A      - INPUT COMPLEX VECTOR OF LENGTH NDEG+1
C                           CONTAINING THE COEFFICIENTS IN ORDER OF
C                           DECREASING POWERS OF THE VARIABLE.
C                         NOTE - THE ROUTINE TREATS A AS A REAL VECTOR
C                           OF LENGTH 2*(NDEG+1). AN APPROPRIATE
C                           EQUIVALENCE STATEMENT MAY BE REQUIRED.
C                           SEE DOCUMENT EXAMPLE.
C                NDEG   - INPUT INTEGER DEGREE OF THE POLYNOMIAL.
C                           NDEG MUST BE GREATER THAN 0 AND LESS
C                           THAN 50.
C                Z      - OUTPUT COMPLEX VECTOR OF LENGTH NDEG
C                           CONTAINING THE COMPUTED ROOTS OF THE
C                           POLYNOMIAL.
C                         NOTE - THE ROUTINE TREATS Z AS A REAL VECTOR
C                           OF LENGTH 2*NDEG. AN APPROPRIATE
C                           EQUIVALENCE STATEMENT MAY BE REQUIRED.
C                           SEE DOCUMENT EXAMPLE.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER=129, INDICATES THAT THE DEGREE OF THE
C                             POLYNOMIAL IS GREATER THAN 49 OR LESS
C                             THAN 1.
C                           IER=130, INDICATES THAT THE LEADING
C                             COEFFICIENT IS ZERO.
C                           IER=131, INDICATES THAT ZCPOLY FOUND FEWER
C                             THAN NDEG ZEROS. IF ONLY M ZEROS ARE
C                             FOUND, Z(J),J=M+1,...,NDEG ARE SET TO
C                             POSITIVE MACHINE INFINITY.
C
C   REQD. IMSL ROUTINES - UERTST,UGETIO,ZCPQLB,ZCPQLC,ZCPQLD,ZCPQLE,
C                           ZCPQLF,ZCPQLG,ZCPQLH,ZCPQLI,ZCPQLJ,ZCPQLK,
C                           ZCPQLL,ZCPQLM
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE ZCPOLY (A,NDEG,Z,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            NDEG,IER
      REAL               A(1),Z(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,ICNT1,ICNT2,II,INX,INXI,J,NN,NN2,NPI,N1,N2
      REAL               PR(50),PI(50),HR(50),HI(50),QPR(50),QPI(50),
     1                   QHR(50),QHI(50),SHR(50),SHI(50)
      REAL               SR,SI,TR,TI,PVR,PVI,ARE,RMRE,REPSR1,RINFP,
     1                   XX,YY,COSR,SINR,REPSP,RADIX,XXX,ZR,ZI,BND,
     2                   ZCPQLL,ZCPQLJ,ZCPQLI,ZERO,ONE,TWO,RSQ2
      LOGICAL            CONV
      COMMON /ZCPQLN/    PR,PI,HR,HI,QPR,QPI,QHR,QHI,SHR,SHI,SR,SI,
     1                   TR,TI,PVR,PVI,ARE,RMRE,REPSR1,RINFP,NN
      DATA               ZERO,ONE,TWO/0.0,1.0,2.0/
      DATA               RSQ2/1.414214/
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
C                                  INITIALIZATION OF CONSTANTS
      IF (NDEG .GT. 49 .OR. NDEG .LT. 1) GO TO 80
      CALL ZCPQLM (REPSR1,RINFP,REPSP,RADIX)
      ARE = REPSR1
      RMRE = TWO*RSQ2*REPSR1
      XX = .7071068
      YY = -XX
      COSR = -.06975647
      SINR = .9975641
      NN = NDEG+1
C                                  ALGORITHM FAILS IF THE LEADING
C                                    COEFFICIENT IS ZERO.
      IF (A(1).NE.ZERO.OR.A(2).NE.ZERO) GO TO 5
      IER = 130
      GO TO 9000
C                                  REMOVE THE ZEROS AT THE ORIGIN IF
C                                    ANY
    5 NN2 = NN+NN
      IF (A(NN2-1).NE.ZERO.OR.A(NN2).NE.ZERO) GO TO 10
      INX = NDEG-NN+2
      INXI = INX+NDEG
      Z(INXI) = ZERO
      Z(INX) = ZERO
      NN = NN-1
      IF (NN .EQ. 1) GO TO 9005
      GO TO 5
C                                  MAKE A COPY OF THE COEFFICIENTS
   10 DO 15 I=1,NN
         II = I+I
         PR(I) = A(II-1)
         PI(I) = A(II)
         SHR(I) = ZCPQLL(PR(I),PI(I))
   15 CONTINUE
C                                  SCALE THE POLYNOMIAL
      BND = ZCPQLJ(NN,SHR,REPSR1,RINFP,REPSP,RADIX)
      IF (BND.EQ.ONE) GO TO 25
      DO 20 I=1,NN
         PR(I) = BND*PR(I)
         PI(I) = BND*PI(I)
   20 CONTINUE
C                                  START THE ALGORITHM FOR ONE ZERO
   25 IF (NN.GT.2) GO TO 30
C                                  CALCULATE THE FINAL ZERO AND RETURN
      CALL ZCPQLK (-PR(2),-PI(2),PR(1),PI(1),Z(NDEG),Z(NDEG+NDEG))
      GO TO 60
C                                  CALCULATE BND, A LOWER BOUND ON THE
C                                    MODULUS OF THE ZEROS
   30 DO 35 I=1,NN
         SHR(I) = ZCPQLL(PR(I),PI(I))
   35 CONTINUE
      BND = ZCPQLI(NN,SHR,SHI)
C                                  OUTER LOOP TO CONTROL 2 MAJOR PASSES
C                                    WITH DIFFERENT SEQUENCES OF
C                                    SHIFTS.
      DO 55 ICNT1=1,2
C                                  FIRST STAGE CALCULATION, NO SHIFT
         CALL ZCPQLB (5)
C                                  INNER LOOP TO SELECT A SHIFT
         DO 50 ICNT2=1,9
C                                  SHIFT IS CHOSEN WITH MODULUS BND AND
C                                    AMPLITUDE ROTATED BY 94 DEGREES
C                                    FROM THE PREVIOUS SHIFT
            XXX = COSR*XX-SINR*YY
            YY = SINR*XX+COSR*YY
            XX = XXX
            SR = BND*XX
            SI = BND*YY
C                                  SECOND STAGE CALCULATION, FIXED
C                                    SHIFT.
            CALL ZCPQLC (10*ICNT2,ZR,ZI,CONV)
            IF (.NOT.CONV) GO TO 45
C                                  THE SECOND STAGE JUMPS DIRECTLY TO
C                                    THE THIRD STAGE ITERATION. IF
C                                    SUCCESSFUL THE ZERO IS STORED AND
C                                    THE POLYNOMIAL DEFLATED.
            INX = NDEG+2-NN
            INXI = INX+NDEG
            Z(INX) = ZR
            Z(INXI) = ZI
            NN = NN-1
            DO 40 I=1,NN
               PR(I) = QPR(I)
               PI(I) = QPI(I)
   40       CONTINUE
            GO TO 25
   45       CONTINUE
C                                  IF THE ITERATION IS UNSUCCESSFUL
C                                    ANOTHER SHIFT IS CHOSEN.
   50    CONTINUE
C                                  IF 9 SHIFTS FAIL, THE OUTER LOOP IS
C                                    REPEATED WITH ANOTHER SEQUENCE OF
C                                    SHIFTS.
   55 CONTINUE
C                                  THE ZEROFINDER HAS FAILED ON TWO
C                                    MAJOR PASSES. RETURN EMPTY HANDED.
C
      IER = 131
C                                  CONVERT ZEROS(Z) IN COMPLEX FORM
   60 DO 65 I=1,NDEG
         NPI=NDEG+I
         PI(I) = Z(NPI)
   65 CONTINUE
      N2 = NDEG+NDEG
      J = NDEG
      DO 70 I=1,NDEG
         Z(N2-1) = Z(J)
         Z(N2) = PI(J)
         N2 = N2-2
         J = J-1
   70 CONTINUE
      IF (IER .EQ. 0) GO TO 9005
C                                  SET UNFOUND ZEROS TO MACHINE INFINITY
      N2 = 2*(NDEG-NN)+3
      N1 = NN-1
      DO 75 I=1,N1
         Z(N2) = RINFP
         Z(N2+1) = RINFP
         N2 = N2+2
   75 CONTINUE
      GO TO 9000
   80 IER = 129
 9000 CONTINUE
      CALL UERTST (IER,6HZCPOLY)
 9005 RETURN
      END
