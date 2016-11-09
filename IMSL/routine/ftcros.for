C   IMSL ROUTINE NAME   - FTCROS
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - MEANS, VARIANCES, CROSS-COVARIANCES, AND
C                           CROSS-CORRELATIONS FOR TWO MUTUALLY
C                           STATIONARY N CHANNEL TIME SERIES
C
C   USAGE               - CALL FTCROS (XY,NC,ISW,EMUSIG,ACV,IA,IB,AC,
C                           IC,ID,IER)
C
C   ARGUMENTS    XY     - INPUT VECTOR OF LENGTH NC(1)*NC(3)+NC(2)*NC(3)
C                           CONTAINING THE TWO TIME SERIES, BOTH IN
C                           TRACE MODE. TRACE MODE IMPLIES THAT THE
C                           ORDERING OF INPUT IS CHANNEL ONE OF THE
C                           FIRST SERIES FOLLOWED BY CHANNEL TWO OF
C                           THE FIRST SERIES, AND SO FORTH FOR THE
C                           FIRST SERIES. THE SECOND SERIES FOLLOWS
C                           THE FIRST SERIES IN XY IN THE SAME MANNER.
C                NC     - INPUT VECTOR OF LENGTH 4.
C                         NC(1) CONTAINS THE LENGTH OF THE FIRST
C                           TIME SERIES. NC(1) MUST BE GREATER THAN OR
C                           EQUAL TO TWO.
C                         NC(2) CONTAINS THE LENGTH OF THE SECOND
C                           TIME SERIES. NC(2) MUST BE GREATER THAN OR
C                           EQUAL TO TWO.
C                         NC(3) CONTAINS THE NUMBER OF CHANNELS IN
C                           EACH TIME SERIES. NC(3) MUST BE GREATER
C                           THAN OR EQUAL TO ONE.
C                         NC(4) CONTAINS THE MAXIMUM NUMBER OF TIME
C                           LAG UNITS FOR THE CROSS-COVARIANCES AND
C                           CROSS-CORRELATIONS. NOT REQUIRED
C                           IF ISW IS EQUAL TO ONE. NC(4) MUST BE
C                           GREATER THAN OR EQUAL TO ZERO AND
C                           LESS THAN THE MINIMUM OF NC(1) AND NC(2).
C                ISW    - INPUT OPTION PARAMETER.
C                         IF ISW = 1, THE ROUTINE COMPUTES THE MEAN AND
C                           VARIANCE FOR EACH CHANNEL OF EACH TIME
C                           SERIES.
C                         IF ISW = 2, THE ROUTINE COMPUTES
C                           CROSS-COVARIANCES.
C                         IF ISW = 3, THE ROUTINE COMPUTES THE MEAN AND
C                           VARIANCE FOR EACH CHANNEL OF EACH TIME
C                           SERIES AND CROSS-COVARIANCES.
C                         IF ISW = 4, THE ROUTINE COMPUTES
C                           CROSS-COVARIANCES AND CROSS-CORRELATIONS.
C                         IF ISW = 5, THE ROUTINE COMPUTES THE MEAN AND
C                           VARIANCE FOR EACH CHANNEL OF EACH TIME
C                           SERIES, CROSS-COVARIANCES, AND CROSS-
C                           CORRELATIONS.
C                EMUSIG - INPUT/OUTPUT VECTOR OF LENGTH 4*NC(3)
C                           CONTAINING THE MEANS FOR THE NC(3) CHANNELS
C                           OF BOTH TIME SERIES IN THE FIRST 2*NC(3)
C                           LOCATIONS.  THE REMAINING LOCATIONS CONTAIN
C                           THE CORRESPONDING VARIANCES.  EMUSIG IS
C                           INPUT IF ISW IS EVEN.  OTHERWISE IT IS
C                           OUTPUT.
C                ACV    - OUTPUT ARRAY OF DIMENSION NC(3) BY NC(3) BY
C                           (NC(4)+1) CONTAINING THE CROSS-COVARIANCES
C                           IN MULTIPLEXED MODE. SEE REMARKS.
C                           ACV IS DEFINED ONLY WHEN ISW IS 2,3,4,OR 5.
C                IA     - INPUT FIRST DIMENSION OF THE ARRAY ACV
C                           EXACTLY AS SPECIFIED IN THE DIMENSION
C                           STATEMENT IN THE CALLING PROGRAM.
C                IB     - INPUT SECOND DIMENSION OF THE ARRAY ACV
C                           EXACTLY AS SPECIFIED IN THE DIMENSION
C                           STATEMENT IN THE CALLING PROGRAM.
C                AC     - OUTPUT ARRAY OF DIMENSION NC(3) BY NC(3) BY
C                           (NC(4)+1) CONTAINING THE CROSS-CORRELATIONS
C                           IN MULTIPLEXED MODE. SEE REMARKS.
C                           AC IS DEFINED ONLY WHEN ISW EQUALS 4 OR 5.
C                IC     - INPUT FIRST DIMENSION OF THE ARRAY AC
C                           EXACTLY AS SPECIFIED IN THE DIMENSION
C                           STATEMENT IN THE CALLING PROGRAM.
C                ID     - INPUT SECOND DIMENSION OF THE ARRAY AC
C                           EXACTLY AS SPECIFIED IN THE DIMENSION
C                           STATEMENT IN THE CALLING PROGRAM.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER=129 INDICATES ISW IS NOT 1,2,3,4, OR 5.
C                           IER=130 INDICATES NC(1) OR NC(2) ARE
C                             LESS THAN 2.
C                           IER=131 INDICATES NC(3) IS LESS THAN 1 OR
C                             THAT NC(4) IS LESS THAN ZERO WHEN ISW
C                             IS 2,3,4, OR 5.
C                           IER=132 INDICATES THE SMALLER OF NC(1) AND
C                             NC(2) IS LESS THAN OR EQUAL TO NC(4) WHEN
C                             ISW IS 2,3,4, OR 5.
C
C   REQD. IMSL ROUTINES - SINGLE/UERTST,UGETIO
C                       - DOUBLE/UERTST,UGETIO,VXADD,VXMUL,VXSTO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS      THE ORDERING OF ELEMENTS IN OUTPUT ARRAYS ACV AND
C                AC IS IMPLIED BY THE FORMULA GIVEN FOR ACV IN THE
C                ALGORITHM SECTION OF THE MANUAL DOCUMENT. FOR EXAMPLE,
C                THE FIRST PLANE OF THE RESPECTIVE THREE DIMENSIONAL
C                ARRAYS CONTAINS THE ZERO TIME LAG RESULTS.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE FTCROS (XY,NC,ISW,EMUSIG,ACV,IA,IB,AC,IC,ID,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            NC(4),ISW,IA,IB,IC,ID,IER
      REAL               XY(1),EMUSIG(1),ACV(IA,IB,1),AC(IC,ID,1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,IE,IEO,IM,IP,IS,ISE,IS1,IT,J,K,KK,KK2,KK3,
     1                   KX,KX1,KY,KY1,NC1,NC2,NC3,NC4,NC4P
      DOUBLE PRECISION   TEMP,GP
C                                  FIRST EXECUTABLE STATEMENT
      IF(ISW .LT. 1 .OR. ISW .GT. 5) GO TO 75
      IEO = ((ISW/2)*2)-ISW
      IF(NC(1) .GT. 1 .AND. NC(2) .GT. 1) GO TO 5
      IER = 130
      GO TO 9000
    5 IER = 0
      NC1 = NC(1)
      NC2 = NC(2)
      NC3 = NC(3)
      NC4 = NC(4)
      IF(ISW .EQ. 1) GO TO 15
      IF(NC3 .GT. 0 .AND. NC4 .GE. 0) GO TO 10
      IER = 131
      GO TO 9000
   10 IF(IEO .EQ. 0) GO TO 50
C                                  COMPUTE THE MEAN VALUES
   15 IE = 0
      ISE = 0
      DO 30 K=1,2
         IS1 = ISE+1
         ISE = ISE+NC3
         DO 25 I=IS1,ISE
            TEMP = 0.0D0
            IS = IE+1
            IE = IE+NC(K)
            DO 20 J=IS,IE
   20       TEMP = TEMP+DBLE(XY(J))
            EMUSIG(I) = TEMP/NC(K)
   25    CONTINUE
   30 CONTINUE
C                                  COMPUTE THE VARIANCES
      IE = 0
      IM = ISE
      ISE = 0
      DO 45 K=1,2
         IS1 = ISE+1
         ISE = ISE+NC3
         DO 40 I=IS1,ISE
            TEMP = 0.0D0
            IS = IE+1
            IE = IE+NC(K)
            DO 35 J=IS,IE
   35       TEMP = TEMP+(DBLE(XY(J))-DBLE(EMUSIG(I)))**2
            IM = IM+1
   40    EMUSIG(IM) = TEMP/NC(K)
   45 CONTINUE
C                                  COMPUTE THE CROSS-COVARIANCES AND
C                                  THE CROSS-CORRELATIONS.
      IF(ISW .EQ. 1) GO TO 9005
   50 KK = NC1*NC3
      IP = MIN0(NC1,NC2)
      IF(IP .LE. NC4) GO TO 80
      GP = IP
      NC4P = NC4+1
      KK2 = 2*NC3
      KK3 = 3*NC3
      IE = NC1*NC3-NC2
      DO 70 I=1,NC4P
         KY = IE
         IS = IP
         IF(NC2 .LT. NC1) IS = MIN0(IP+I-1,NC1)
         DO 65 J=1,NC3
            KY = KY+NC2
            KX = -NC1
            DO 60 K=1,NC3
               KX = KX+NC1
               TEMP = 0.0D0
               KY1 = KY
               DO 55 IT=I,IS
                  KX1 = KX+IT
                  KY1 = KY1+1
                  TEMP = TEMP+((DBLE(XY(KX1))-DBLE(EMUSIG(K)))*
     *              (DBLE(XY(KY1))-DBLE(EMUSIG(NC3+J))))
   55          CONTINUE
               ACV(K,J,I) = TEMP/GP
               IF(ISW .LT. 4) GO TO 60
               AC(K,J,I) = ACV(K,J,I)/SQRT(EMUSIG(KK2+K)*EMUSIG(KK3+J))
   60       CONTINUE
   65    CONTINUE
   70 CONTINUE
      GO TO 9005
   75 IER = 129
      GO TO 9000
   80 IER = 132
 9000 CONTINUE
      CALL UERTST(IER,'FTCROS')
 9005 CONTINUE
      RETURN
      END
