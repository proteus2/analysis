C   IMSL ROUTINE NAME   - BEGRPS
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1982
C
C   PURPOSE             - MOMENTS ESTIMATION FOR GROUPED DATA WITH
C                           AND WITHOUT SHEPPARDS CORRECTIONS
C
C   USAGE               - CALL BEGRPS (N,C,CI,U,UC,IER)
C
C   ARGUMENTS    N      - NUMBER OF GROUPS OF DATA. (INPUT)
C                C      - VECTOR OF LENGTH N CONTAINING FREQUENCIES
C                           WITHIN GROUPS. (INPUT)
C                CI     - VECTOR OF LENGTH TWO CONTAINING THE CENTER OF
C                           THE LOWEST CLASS INTERVAL, CI(1), AND THE
C                           CLASS WIDTH, CI(2). (INPUT)
C                U      - VECTOR OF LENGTH FOUR CONTAINING UNCORRECTED
C                           CENTRAL (EXCEPTING THE MEAN) MOMENTS.
C                           (OUTPUT)
C                UC     - VECTOR OF LENGTH FOUR CONTAINING CORRECTED
C                           CENTRAL MOMENTS. (OUTPUT)
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER = 129 INDICATES THAT N IS LESS THAN 1.
C                           IER = 130 INDICATES THAT THE SUM OF THE
C                             FREQUENCIES IS WITHIN 10**(-6) OF ZERO.
C
C   REQD. IMSL ROUTINES - UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE BEGRPS (N,C,CI,U,UC,IER)
C
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N,IER
      REAL               C(N),CI(2),U(4),UC(4)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,J
      REAL               ASIX,ZERO
      DOUBLE PRECISION   DBH,DZERO,E,H,ONE,P5,SXTY,T,T1,Z
      DATA               ZERO /0.0/,P5 /0.5D0/,ASIX /.1666667/,
     *                   SXTY /.1166667D0/,ONE /1.0D0/,DZERO /0.D0/
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
      IF (N.LE.0) GO TO 25
      T = DZERO
      U(1) = ZERO
      H = CI(2)
      DBH = DBLE(CI(1))-H
      DO 5 I=1,N
         T = T+DBLE(C(I))
    5 CONTINUE
      IF ((DABS(T)).GT.1.D-6) GO TO 10
      IER = 130
      GO TO 9000
   10 T1 = ONE/T
      DO 20 J=1,4
         E = DBH-DBLE(U(1))
         Z = DZERO
         DO 15 I=1,N
            E = E+H
            Z = Z+DBLE(C(I))*E**J
   15    CONTINUE
         U(J) = Z*T1
   20 CONTINUE
      E = P5*H*H
      UC(1) = U(1)
      UC(2) = U(2)-SNGL(E)*ASIX
      UC(3) = U(3)-U(1)*SNGL(P5*E)
      UC(4) = U(4)-U(2)*SNGL(E)+SNGL(E*E*SXTY)
      GO TO 9005
   25 IER = 129
 9000 CONTINUE
      CALL UERTST(IER,'BEGRPS')
 9005 RETURN
      END
