C   IMSL ROUTINE NAME   - VNRMS1
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - 1-NORM OF MATRICES (SYMMETRIC STORAGE MODE)
C
C   USAGE               - CALL VNRMS1 (A,N,XNRMA)
C
C   ARGUMENTS    A      - N BY N MATRIX (STORED IN SYMMETRIC STORAGE
C                           MODE AS A VECTOR) FOR WHICH THE 1-NORM WILL
C                           BE COMPUTED. (INPUT)
C                N      - ORDER OF THE MATRIX STORED IN VECTOR A.
C                           (INPUT)
C                XNRMA  - VALUE OF THE 1-NORM OF A. (OUTPUT)
C
C
C   REQD. IMSL ROUTINES - SINGLE/VABSMS
C                       - DOUBLE/VABSMS,VXADD,VXSTO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE VNRMS1 (A,N,XNRMA)
C
      DIMENSION          A(1)
C                                  FIRST EXECUTABLE STATEMENT
      XNRMA=0.0
         DO 5 J=1,N
         CALL VABSMS(A(1),N,J,VSUM)
         IF(VSUM.GE.XNRMA) XNRMA=VSUM
    5    CONTINUE
      RETURN
      END
