C   IMSL ROUTINE NAME   - RLEAP1
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - NUCLEUS CALLED ONLY BY IMSL SUBROUTINE
C                           RLEAP
C
C   REQD. IMSL ROUTINES - NONE REQUIRED
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE RLEAP1 (NC,LB,L,IPI,MV,RS,BND,ILI,JC,ID,XI,II,NI,ND,KZ)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            LB,L,MV,JC,II,ND,KZ,ID(ND),NC(ND,ND),
     1                   ILI(ND,ND),IPI(ND),NI(ND)
      REAL               RS,BND,XI(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            J,I1,IS,KZZ,LBB,IN,MM,IP,IJ,I2,KA,KN,I3
      REAL               B
C                                  FIRST EXECUTABLE STATEMENT
    5 I1 = IPI(JC)
C                                  FIND SOURCE MATRIX
      IF (LB .LE. NI(I1)) GO TO 10
      JC = JC-1
      GO TO 5
C                                  ADJUST FOR PREVIOUS PIVOTS
   10 KZZ = (KZ*(KZ-1))/2
      LBB = (LB*(LB-1))/2
      DO 25 J = JC,MV
         IN = IPI(J)
         L = ILI(IN,LB)
         MM = ID(IN)
         IF (J .EQ. MV) GO TO 30
         IS = IPI(J+1)
         IP = ILI(IN,IS-1)
         I1 = MAX0(IP,L)
         IJ = MIN0(IP,L)
         I1 = MM+(I1*(I1-1))/2+IJ
         I2 = MM+(IP*(IP+1))/2
         B = XI(I1)/XI(I2)
         KA = IS
   15    IF (KA .GT. LB) GO TO 20
         KN = ILI(IN,KA)
         I1 = ID(IS)+LBB+KA
         I2 = MAX0(KN,L)
         IJ = MIN0(KN,L)
         I2 = MM+(I2*(I2-1))/2+IJ
         I3 = MAX0(KN,IP)
         IJ = MIN0(KN,IP)
         I3 = MM+(I3*(I3-1))/2+IJ
         XI(I1) = XI(I2)-B*XI(I3)
         KA = KA+1
         GO TO 15
   20    I1 = ID(IS)+KZZ+LB
         I2 = MM+KZZ+L
         I3 = MM+KZZ+IP
         XI(I1) = XI(I2)-B*XI(I3)
         NI(IS) = LB
         ILI(IS,LB) = LB
         IF (II .EQ. 0) NC(IS,LB) = NC(IN,L)
   25 CONTINUE
C                                  CURRENT PIVOT
   30 I1 = MM+KZZ+L
      I2 = MM+(L*(L+1))/2
      RS = BND-XI(I1)*XI(I1)/XI(I2)
      RETURN
      END
