C   IMSL ROUTINE NAME   - RLPRDI
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - CONFIDENCE INTERVALS FOR THE TRUE RESPONSE
C                           AND FOR THE AVERAGE OF A SET OF FUTURE
C                           OBSERVATIONS ON THE RESPONSE - IN CORE
C                           VERSION
C
C   USAGE               - CALL RLPRDI (Y,V,N,TS,NR,C,IC)
C
C   ARGUMENTS    Y      - VECTOR OF LENGTH N CONTAINING
C                           THE  PREDICTED RESPONSES. (INPUT)
C                V      - VECTOR OF LENGTH N CONTAINING THE ESTIMATED
C                           VARIANCES OF THE PREDICTED RESPONSES,
C                           DIVIDED BY THE ERROR MEAN SQUARE. (INPUT)
C                           ALL ELEMENTS IN V MUST BE NON-NEGATIVE.
C                N      - NUMBER OF POINTS IN THE DESIGN SPACE AT WHICH
C                           CONFIDENCE INTERVALS ARE DESIRED. (INPUT)
C                TS     - THE PRODUCT OF THE APPROPRIATE STUDENTS T
C                           VALUE AND THE SQUARE ROOT OF THE ERROR MEAN
C                           SQUARE. (INPUT)
C                NR     - NUMBER OF FUTURE OBSERVATIONS FOR WHICH A
C                           CONFIDENCE INTERVAL IS DESIRED ON THE
C                           AVERAGE OF THESE OBSERVATIONS. (INPUT)
C                C      - N BY 4 MATRIX OF CONFIDENCE INTERVALS AT EACH
C                           OF THE N DESIGN POINTS. (OUTPUT)
C                           COLUMNS 1 AND 2 CONTAIN THE LOWER AND UPPER
C                           LIMITS, RESPECTIVELY, FOR THE TRUE MEAN.
C                           COLUMNS 3 AND 4 CONTAIN THE LOWER AND UPPER
C                           LIMITS, RESPECTIVELY, FOR THE MEAN OF NR
C                           FUTURE OBSERVATIONS.
C                IC     - ROW DIMENSION OF MATRIX C EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN
C                           THE CALLING PROGRAM. (INPUT)
C
C   REQD. IMSL ROUTINES - NONE REQUIRED
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE RLPRDI (Y,V,N,TS,NR,C,IC)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N,NR,IC
      REAL               Y(1),V(1),TS,C(IC,4)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I
      REAL               D1,D2,Z,XNR,ONE
      DATA               ONE/1.0/
C                                  FIRST EXECUTABLE STATEMENT
      XNR=ONE/NR
      DO 5 I=1,N
         D1=TS*SQRT(V(I))
         D2=TS*SQRT(XNR+V(I))
         Z=Y(I)
         C(I,1)=Z-D1
         C(I,2)=Z+D1
         C(I,3)=Z-D2
         C(I,4)=Z+D2
    5 CONTINUE
      RETURN
      END
