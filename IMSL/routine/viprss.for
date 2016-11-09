C   IMSL ROUTINE NAME   - VIPRSS
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - VECTOR INNER PRODUCT OF TWO VECTORS EACH OF
C                           WHICH IS PART OF SOME MATRIX STORED IN
C                           SYMMETRIC MODE
C
C   USAGE               - CALL VIPRSS (X,Y,L,IX,IY,XYIP)
C
C   ARGUMENTS    X      - NAME OF THE FIRST MATRIX (ADDRESS OF THE
C                           FIRST ELEMENT). (INPUT)
C                Y      - NAME OF THE SECOND MATRIX (ADDRESS OF THE
C                           FIRST ELEMENT). (INPUT)
C                L      - NUMBER OF ELEMENTS OF VECTOR X OR Y INVOLVED
C                           IN THE INNER PRODUCT. (INPUT)
C                IX     - X IS ROW (OR COLUMN) NUMBER IX OF A SYMMETRIC
C                           MATRIX. (INPUT)
C                IY     - Y IS ROW (OR COLUMN) NUMBER IY OF A SYMMETRIC
C                           MATRIX. (INPUT)
C                XYIP   - RESULTANT INNER PRODUCT. (OUTPUT)
C
C   REQD. IMSL ROUTINES - SINGLE/NONE REQUIRED
C                       - DOUBLE/VXADD,VXMUL,VXSTO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE VIPRSS (X,Y,L,IX,IY,XYIP)
C
      DIMENSION          X(1),Y(1)
      DOUBLE PRECISION   TEMP
C                                  FIRST EXECUTABLE STATEMENT
      I=(IX*(IX-1))/2+1
      J=(IY*(IY-1))/2+1
      TEMP = 0.D0
         DO 20 K=1,L
         TEMP=DBLE(X(I))*DBLE(Y(J))+TEMP
         IF(K .LT. IX) GO TO 5
         I=I+K
         GO TO 10
    5    I=I+1
   10    IF (K .LT. IY) GO TO 15
         J=J+K
         GO TO 20
   15    J=J+1
   20    CONTINUE
      XYIP=TEMP
      RETURN
      END
