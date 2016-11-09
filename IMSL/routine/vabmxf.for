C   IMSL ROUTINE NAME   - VABMXF
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - MAXIMUM ABSOLUTE VALUE OF THE ELEMENTS OF A
C                           VECTOR OR A SUBSET OF THE ELEMENTS OF A
C                           VECTOR
C
C   USAGE               - CALL VABMXF (V,L,INC,J,VMAX)
C
C   ARGUMENTS    V      - THE VECTOR (ADDRESS OF THE FIRST ELEMENT).
C                           (INPUT)
C                L      - LENGTH OF THE SUBSET VECTOR. (INPUT)
C                INC    - INCREMENT BETWEEN SUCCESSIVE ELEMENTS OF THE
C                           SUBSET VECTOR. (INPUT)
C                J      - INDEX IN THE SUBSET VECTOR OF THE ELEMENT
C                           HAVING THE MAXIMUM ABSOLUTE VALUE. (OUTPUT)
C                VMAX   - MAXIMUM ABSOLUTE VALUE. (OUTPUT)
C
C   REQD. IMSL ROUTINES - NONE REQUIRED
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE VABMXF (V,L,INC,J,VMAX)
C
      DIMENSION          V(1)
      DATA               ZERO/0.0/
C                                  FIRST EXECUTABLE STATEMENT
      VMAX = ZERO
      J=1
      IEND=1+(L-1)*INC
         DO 5 I=1,IEND,INC
         T=ABS(V(I))
         IF (T .LE. VMAX) GO TO 5
         VMAX=T
         J=I
    5    CONTINUE
      J=(J-1)/INC+1
      RETURN
      END
