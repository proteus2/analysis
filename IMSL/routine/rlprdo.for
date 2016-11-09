C   IMSL ROUTINE NAME   - RLPRDO
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - CONFIDENCE INTERVALS FOR THE TRUE RESPONSE
C                           AND FOR THE AVERAGE OF A SET OF FUTURE
C                           OBSERVATIONS ON THE RESPONSE - OUT OF CORE
C                           VERSION
C
C   USAGE               - CALL RLPRDO (Y,V,TS,NR,C)
C
C   ARGUMENTS    Y      - PREDICTED RESPONSE VALUE. (INPUT)
C                V      - ESTIMATED VARIANCE OF PREDICTED RESPONSE
C                           DIVIDED BY THE ERROR MEAN SQUARE. (INPUT)
C                           V MUST BE NON-NEGATIVE.
C                TS     - THE PRODUCT OF THE APPROPRIATE STUDENTS T
C                           VALUE AND THE SQUARE ROOT OF THE ERROR
C                           MEAN SQUARE. (INPUT)
C                NR     - NUMBER OF FUTURE OBSERVATIONS FOR WHICH A
C                           CONFIDENCE INTERVAL IS DESIRED ON THE
C                           AVERAGE OF THE OBSERVATIONS. (INPUT)
C                C      - VECTOR OF LENGTH 4 CONTAINING CONFIDENCE
C                           INTERVALS AT THE DESIGN POINT. (OUTPUT)
C                           ELEMENTS 1 AND 2 CONTAIN THE LOWER AND UPPER
C                           LIMITS, RESPECTIVELY, FOR THE TRUE MEAN.
C                           ELEMENTS 3 AND 4 CONTAIN THE LOWER AND UPPER
C                           LIMITS, RESPECTIVELY, FOR THE MEAN OF NR
C                           FUTURE OBSERVATIONS.
C
C   REQD. IMSL ROUTINES - NONE REQUIRED
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE RLPRDO (Y,V,TS,NR,C)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            NR
      REAL               Y,V,TS,C(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      REAL               D1,D2
C                                  FIRST EXECUTABLE STATEMENT
      D1=TS*SQRT(V)
      D2=TS*SQRT(1./NR+V)
      C(1)=Y-D1
      C(2)=Y+D1
      C(3)=Y-D2
      C(4)=Y+D2
      RETURN
      END
