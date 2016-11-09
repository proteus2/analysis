C   IMSL ROUTINE NAME   - VABSMF
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - SUM OF THE ABSOLUTE VALUES OF THE ELEMENTS OF
C                           A VECTOR OR A SUBSET OF A VECTOR
C
C   USAGE               - CALL VABSMF (V,L,INC,VSUM)
C
C   ARGUMENTS    V      - NAME OF THE VECTOR (ADDRESS OF THE FIRST
C                           ELEMENT). (INPUT)
C                L      - NUMBER OF ELEMENTS TO BE SUMMED. (INPUT)
C                INC    - INCREMENT BETWEEN SUCCESSIVE ELEMENTS OF THE
C                           SUBSET VECTOR. (INPUT)
C                VSUM   - SUM OF THE ABSOLUTE VALUES. (OUTPUT)
C
C   REQD. IMSL ROUTINES - NONE REQUIRED
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE VABSMF (V,L,INC,VSUM)
C
      REAL               V(1),VSUM
C                                  FIRST EXECUTABLE STATEMENT
      VSUM = 0.0
      J=1+(L-1)*INC
      DO 5 I=1,J,INC
    5 VSUM=VSUM+ABS(V(I))
      RETURN
      END
