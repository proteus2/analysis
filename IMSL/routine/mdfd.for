C   IMSL ROUTINE NAME   - MDFD
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1981
C
C   PURPOSE             - F PROBABILITY DISTRIBUTION FUNCTION
C
C   USAGE               - CALL MDFD (F,N1,N2,P,IER)
C
C   ARGUMENTS    F      - INPUT CONSTANT TO WHICH INTEGRATION IS
C                           PERFORMED. F MUST BE GREATER THAN OR EQUAL
C                           TO ZERO.
C                N1     - INPUT FIRST DEGREE OF FREEDOM. A POSITIVE
C                           INTEGER.
C                N2     - INPUT SECOND DEGREE OF FREEDOM, A POSITIVE
C                           INTEGER.
C                P      - OUTPUT PROBABILITY THAT A RANDOM VARIABLE
C                           FOLLOWING THE F DISTRIBUTION WITH DEGREES
C                           OF FREEDOM N1 AND N2 WILL BE LESS THAN OR
C                           EQUAL TO INPUT F.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                           TERMINAL ERROR
C                             IER = 129 INDICATES EITHER N1 OR N2 IS
C                               LESS THAN ONE OR N1+N2 IS GREATER THAN
C                               20,000.  P IS SET TO POSITIVE MACHINE
C                               INFINITY.
C                             IER = 130 INDICATES F IS LESS THAN ZERO.
C                               P IS SET TO POSITIVE MACHINE INFINITY.
C
C   REQD. IMSL ROUTINES - MERRC=ERFC,UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE MDFD   (F,N1,N2,P,IER)
C
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N1,N2,IER
      REAL               F,P
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I1,I2P,I2,I,L2,MNM,MXM
      REAL               ACONS,A,BIGE,B,CBR1,CBR2,C,DPL,DP,F1F,F1,F2P,
     1                   F2,R1D3,R2D9,R2DPI,RINFP,STS,S,TEMP1,TEMP,
     2                   THETA,VP,X1,X2,XI
      DATA               RINFP/Z7FFFFFFF/
      DATA               R2D9/.2222222/,R1D3/.3333333/
C                                  R2DPI = 2/PI
      DATA               R2DPI/.6366198/
      DATA               BIGE/174.6/,ACONS/1.E50/
C                                  FIRST EXECUTABLE STATEMENT
C                                  TEST FOR INVALID INPUT
      MXM = MAX0(N1,N2)
      MNM = MIN0(N1,N2)
      IF (MNM.LT.1.OR.MXM.GT.(20000-MNM)) GO TO 100
      IF (F.LT.0.0) GO TO 105
      IER = 0
      IF (F.EQ.0.0) GO TO 115
      F1 = N1
      F2 = N2
      DP = 0.0
      VP = F1+F2-2.0
      F1F = F1*F
      F2P = F2+F1F
      X1 = F2/F2P
      X2 = 1.0-X1
      IF (X2.EQ.0.0) GO TO 115
      IF ((N1/2)*2-N1.EQ.0.AND.N1.LE.500) GO TO 5
      IF ((N2/2)*2-N2.EQ.0.AND.N2.LE.500) GO TO 30
      IF (N1+N2.LE.500) GO TO 55
      F1 = R2D9/F1
      F2 = R2D9/F2
      CBR1 = R1D3*ALOG(F)
      IF (ABS(CBR1).GT.BIGE) GO TO 120
      CBR1 = EXP(CBR1)
      CBR2 = CBR1*CBR1
      S = (CBR1*(1.0-F2)-1.0+F1)/SQRT(F1+CBR2*F2)
      P=.7071068
      P = .5*ERFC(-P*S)
      GO TO 95
C                                  N1 IS EVEN AND LESS THAN 500
    5 TEMP1 = 0.
      TEMP = .5*F2*ALOG(X1)
      IF (N1.EQ.2) GO TO 25
      I1 = N1-2
      XI = F1
      DO 10 I2=2,I1,2
         L2 = I2
         XI = XI-2.
         VP = VP-2.
         DP = X2*VP/XI*(1.+DP)
         IF (DP.GT.ACONS) GO TO 15
   10 CONTINUE
      GO TO 25
   15 IF (L2.GE.I1) GO TO 25
      DPL = ALOG(DP)
      I2P = L2+2
      XI = F1-I2P
      DO 20 I2=I2P,I1,2
         VP = VP-2.
         DPL = DPL+ALOG(X2*VP/XI)
         XI = XI-2.
   20 CONTINUE
      TEMP = TEMP+DPL
      IF (ABS(TEMP).LE.BIGE) TEMP1 = EXP(TEMP)
      P = 1.-TEMP1
      GO TO 95
   25 IF (ABS(TEMP).LE.BIGE) TEMP1 = EXP(TEMP)
      P = 1.0-TEMP1*(1.0+DP)
      GO TO 95
C                                  N2 IS EVEN AND LESS THAN 500
   30 TEMP1 = 0.
      TEMP = .5*F1*ALOG(X2)
      IF (N2.EQ.2) GO TO 50
      I1 = N2-2
      XI = F2
      DO 35 I2=2,I1,2
         L2 = I2
         XI = XI-2.
         VP = VP-2.
         DP = X1*VP/XI*(1.+DP)
         IF (DP.GT.ACONS) GO TO 40
   35 CONTINUE
      GO TO 50
   40 IF (L2.GE.I1) GO TO 50
      DPL = ALOG(DP)
      I2P = L2+2
      XI = F2-I2P
      DO 45 I2=I2P,I1,2
         VP = VP-2.
         DPL = DPL+ALOG(X1*VP/XI)
         XI = XI-2.
   45 CONTINUE
      TEMP = TEMP+DPL
      IF (ABS(TEMP).LE.BIGE) TEMP1 = EXP(TEMP)
      P = TEMP1
      GO TO 95
   50 IF (ABS(TEMP).LE.BIGE) TEMP1 = EXP(TEMP)
      P = TEMP1*(1.+DP)
      GO TO 95
C                                  SUM OF DFS ARE LE 500 AND ODD
   55 DP = SQRT(F1F/F2)
      THETA = ATAN(DP)
      STS = F1F/F2P
      A = 0.0
      B = 0.0
      IF (N2.EQ.1) GO TO 70
      IF (N2.EQ.3) GO TO 65
      I1 = N2-3
      XI = F2
      DO 60 I2=2,I1,2
         XI = XI-2.
         A = X1*(XI-1.0)/XI*(1.0+A)
   60 CONTINUE
   65 A = X1*DP*(1.0+A)
   70 A = A+THETA
      IF (N1.EQ.1) GO TO 90
      IF (N1.EQ.3) GO TO 80
      I1 = N1-3
      XI = F1
      DO 75 I2=2,I1,2
         XI = XI-2.
         VP = VP-2.
         B = STS*VP/XI*(1.0+B)
   75 CONTINUE
   80 B = DP*X1*(1.0+B)
      IF (N2.EQ.1) GO TO 90
      I2 = N2/2
      C = 1.0
      DO 85 I=1,I2
         B = B*X1*C/(C-0.5)
         C = C+1.0
   85 CONTINUE
   90 P = R2DPI*(A-B)
   95 IF (P.LT.0.0) P = 0.0
      IF (P.GT.1.0) P = 1.0
      GO TO 9005
  100 IER = 129
      GO TO 110
  105 IER = 130
  110 P = RINFP
      GO TO 9000
  115 P = 0.0
      GO TO 9005
  120 P = .5
      GO TO 9005
 9000 CONTINUE
      CALL UERTST (IER,6HMDFD  )
 9005 RETURN
      END
