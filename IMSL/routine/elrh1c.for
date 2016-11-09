C   IMSL ROUTINE NAME   - ELRH1C
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1982
C
C   PURPOSE             - NUCLEUS CALLED ONLY BY IMSL ROUTINE EIGCC
C
C   REQD. IMSL ROUTINES - UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE ELRH1C (HR,HI,K,L,N,IH,WR,WI,INFER,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            K,L,N,IH,INFER,IER
      REAL               HR(IH,1),HI(IH,1),WR(1),WI(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,NN,ITS,NM1,NPL,LL,M,MM1,NMJ,J,M1,MM,MP1,
     *                   IM1,JM1
      REAL               T1(2),T2(2),T3(2),ZERO,ONE,TWO,RDELP,TR,TI,
     *                   SR,SI,XR,XI,YR,YI,ZR,ZI
      COMPLEX            X,Y,Z
      EQUIVALENCE        (X,T1(1)),(T1(1),XR),(T1(2),XI),
     1                   (Y,T2(1)),(T2(1),YR),(T2(2),YI),
     2                   (Z,T3(1)),(T3(1),ZR),(T3(2),ZI)
      DATA               ZERO,ONE,TWO/0.0,1.0,2.0/
      DATA               RDELP/Z3C100000/
C                                  FIRST EXECUTABLE STATEMENT
      INFER=0
      IER=0
      DO 5 I=1,N
         IF (I .GE. K .AND. I .LE. L) GO TO 5
         WR(I)=HR(I,I)
         WI(I)=HI(I,I)
    5 CONTINUE
      NN=L
      TR=ZERO
      TI=ZERO
C                                  SEARCH FOR NEXT EIGENVALUE
   10 IF (NN .LT. K) GO TO 9005
      ITS=0
      NM1=NN-1
      IF (NN .EQ. K) GO TO 25
C                                  LOOK FOR SINGLE SMALL SUB-DIAGONAL
C                                    ELEMENTS
   15 NPL=NN+K
      DO 20 LL=K,NM1
         M=NPL-LL
         MM1=M-1
         IF ( ABS(HR(M,MM1)) +  ABS(HI(M,MM1)) .LE.
     1    RDELP*( ABS(HR(MM1,MM1)) +  ABS(HI(MM1,MM1)) +
     2       ABS(HR(M,M)) +  ABS(HI(M,M)))) GO TO 30
   20 CONTINUE
   25 M=K
   30 IF (M .EQ. NN) GO TO 110
      IF (ITS .EQ. 30) GO TO 115
C                                  FORM SHIFT
      IF (ITS .EQ. 10 .OR. ITS .EQ. 20) GO TO 35
      SR=HR(NN,NN)
      SI=HI(NN,NN)
      XR=HR(NM1,NN)*HR(NN,NM1)-HI(NM1,NN)*HI(NN,NM1)
      XI=HR(NM1,NN)*HI(NN,NM1)+HI(NM1,NN)*HR(NN,NM1)
      IF (XR .EQ. ZERO .AND. XI .EQ. ZERO) GO TO 40
      YR=(HR(NM1,NM1)-SR)/TWO
      YI=(HI(NM1,NM1)-SI)/TWO
      Z= CSQRT( CMPLX(YR**2-YI**2+XR,TWO*YR*YI+XI))
      IF (YR*ZR+YI*ZI .LT. ZERO) Z=-Z
      X=X/(Y+Z)
      SR=SR-XR
      SI=SI-XI
      GO TO 40
   35 SR= ABS(HR(NN,NM1))+ ABS(HR(NM1,NN-2))
      SI= ABS(HI(NN,NM1))+ ABS(HI(NM1,NN-2))
   40 DO 45 I=K,NN
         HR(I,I)=HR(I,I)-SR
         HI(I,I)=HI(I,I)-SI
   45 CONTINUE
      TR=TR+SR
      TI=TI+SI
      ITS=ITS+1
C                                  LOOK FOR TWO CONSECUTIVE SMALL
C                                    SUB-DIAGONAL ELEMENTS
      XR= ABS(HR(NM1,NM1))+ ABS(HI(NM1,NM1))
      YR= ABS(HR(NN,NM1))+ ABS(HI(NN,NM1))
      ZR= ABS(HR(NN,NN))+ ABS(HI(NN,NN))
      NMJ=NM1-M
      IF (NMJ .EQ. 0) GO TO 55
C                                  FOR MM=NN-1 STEP -1 UNTIL M+1 DO
      DO 50 J=1,NMJ
         MM=NN-J
         M1=MM-1
         YI=YR
         YR= ABS(HR(MM,M1))+ ABS(HI(MM,M1))
         XI=ZR
         ZR=XR
         XR= ABS(HR(M1,M1))+ ABS(HI(M1,M1))
         IF (YR.LE.RDELP*ZR/YI*(ZR+XR+XI)) GO TO 60
   50 CONTINUE
   55 MM=M
C                                  TRIANGULAR DECOMPOSITION
   60 MP1=MM+1
      DO 85 I=MP1,NN
         IM1=I-1
         XR=HR(IM1,IM1)
         XI=HI(IM1,IM1)
         YR=HR(I,IM1)
         YI=HI(I,IM1)
         IF ( ABS(XR)+ ABS(XI) .GE.  ABS(YR)+ ABS(YI)) GO TO 70
C                                  INTERCHANGE ROWS OF HR AND HI
         DO 65 J=IM1,NN
            ZR=HR(IM1,J)
            HR(IM1,J)=HR(I,J)
            HR(I,J)=ZR
            ZI=HI(IM1,J)
            HI(IM1,J)=HI(I,J)
            HI(I,J)=ZI
   65    CONTINUE
         Z=X/Y
         WR(I)=ONE
         GO TO 75
   70    Z=Y/X
         WR(I)=-ONE
   75    HR(I,IM1)=ZR
         HI(I,IM1)=ZI
         DO 80 J=I,NN
            HR(I,J)=HR(I,J)-ZR*HR(IM1,J)+ZI*HI(IM1,J)
            HI(I,J)=HI(I,J)-ZR*HI(IM1,J)-ZI*HR(IM1,J)
   80    CONTINUE
   85 CONTINUE
C                                  COMPOSITION
      DO 105 J=MP1,NN
         JM1=J-1
         XR=HR(J,JM1)
         XI=HI(J,JM1)
         HR(J,JM1)=ZERO
         HI(J,JM1)=ZERO
C                                  INTERCHANGE COLUMNS OF HR AND HI IF
C                                    NECESSARY
         IF (WR(J) .LE. ZERO) GO TO 95
         DO 90 I=M,J
            ZR=HR(I,JM1)
            HR(I,JM1)=HR(I,J)
            HR(I,J)=ZR
            ZI=HI(I,JM1)
            HI(I,JM1)=HI(I,J)
            HI(I,J)=ZI
   90    CONTINUE
   95    DO 100 I=M,J
            HR(I,JM1)=HR(I,JM1)+XR*HR(I,J)-XI*HI(I,J)
            HI(I,JM1)=HI(I,JM1)+XR*HI(I,J)+XI*HR(I,J)
  100    CONTINUE
  105 CONTINUE
      GO TO 15
C                                  A ROOT FOUND
  110 WR(NN)=HR(NN,NN)+TR
      WI(NN)=HI(NN,NN)+TI
      NN=NM1
      GO TO 10
C                                  SET ERROR-NO CONVERGENCE TO AN
C                                    EIGENVALUE AFTER 30 ITERATIONS
  115 INFER=NN
      IER=129
 9000 CONTINUE
      CALL UERTST (IER,6HELRH1C)
 9005 RETURN
      END
