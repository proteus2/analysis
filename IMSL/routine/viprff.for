C   IMSL ROUTINE NAME   - VIPRFF
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - NOVEMBER 1, 1984
C
C   PURPOSE             - VECTOR INNER PRODUCT OF TWO VECTORS OR
C                           SUBSETS OF TWO VECTORS
C
C   USAGE               - CALL VIPRFF (X,Y,L,IX,IY,XYIP)
C
C   ARGUMENTS    X      - FIRST VECTOR (ADDRESS OF THE FIRST ELEMENT).
C                           (INPUT)
C                Y      - SECOND VECTOR (ADDRESS OF THE FIRST ELEMENT).
C                           (INPUT)
C                L      - NUMBER OF ELEMENTS OF VECTOR X OR Y INVOLVED
C                           IN THE INNER PRODUCT. (INPUT)
C                IX     - INCREMENT BETWEEN SUCCESSIVE ELEMENTS OF X.
C                           (INPUT)
C                IY     - INCREMENT BETWEEN SUCCESSIVE ELEMENTS OF Y.
C                           (INPUT)
C                XYIP   - INNER PRODUCT. (OUTPUT)  XYIP = X(1)*Y(1)+
C                           X(1+IX)*Y(1+IY)+...+X(1+L*IX)*Y(1+L*IY).
C
C   REQD. IMSL ROUTINES - SINGLE/VBLA=SDOT
C                       - DOUBLE/VBLA=DDOT
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE VIPRFF (X,Y,L,IX,IY,XYIP)
C                                 SPECIFICATIONS FOR ARGUMENTS
      REAL               X(1),XYIP,Y(1)
      INTEGER            IX,IY,L
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      REAL               SDOT
C                                  FIRST EXECUTABLE STATEMENT
      XYIP = SDOT(L,X,IX,Y,IY)
      RETURN
      END
