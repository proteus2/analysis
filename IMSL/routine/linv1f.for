C   IMSL ROUTINE NAME   - LINV1F
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1982
C
C   PURPOSE             - INVERSION OF A MATRIX - FULL STORAGE MODE -
C                           SPACE ECONOMIZER SOLUTION
C
C   USAGE               - CALL LINV1F (A,N,IA,AINV,IDGT,WKAREA,IER)
C
C   ARGUMENTS    A      - INPUT MATRIX OF DIMENSION N BY N CONTAINING
C                           THE MATRIX TO BE INVERTED.
C                         ON OUTPUT, A IS REPLACED BY THE LU
C                           DECOMPOSITION OF A ROWWISE PERMUTATION OF A.
C                N      - ORDER OF A. (INPUT)
C                IA     - ROW DIMENSION OF MATRICES A AND AINV EXACTLY
C                           AS SPECIFIED IN THE DIMENSION STATEMENT IN
C                           THE CALLING PROGRAM. (INPUT)
C                AINV   - OUTPUT MATRIX OF DIMENSION N BY N CONTAINING
C                           THE INVERSE OF A. A AND AINV MUST OCCUPY
C                           SEPARATE CORE LOCATIONS.
C                IDGT   - INPUT OPTION.
C                         IF IDGT IS GREATER THAN 0, THE ELEMENTS OF A
C                           ARE ASSUMED TO BE CORRECT TO IDGT DECIMAL
C                           DIGITS AND THE ROUTINE PERFORMS AN ACCURACY
C                           TEST.
C                         IF IDGT EQUALS ZERO, THE ACCURACY TEST IS
C                           BYPASSED.
C                WKAREA - WORK AREA OF DIMENSION GREATER THAN OR EQUAL
C                           TO N.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER=129 INDICATES THAT A IS
C                             ALGORITHMICALLY SINGULAR. (SEE THE
C                             CHAPTER L PRELUDE.)
C                         WARNING ERROR
C                           IER=34 INDICATES THAT THE ACCURACY TEST
C                             FAILED. THE COMPUTED SOLUTION MAY BE IN
C                             ERROR BY MORE THAN CAN BE ACCOUNTED FOR
C                             BY THE UNCERTAINTY OF THE DATA. THIS
C                             WARNING CAN BE PRODUCED ONLY IF IDGT IS
C                             GREATER THAN 0 ON INPUT. SEE CHAPTER L
C                             PRELUDE FOR FURTHER DISCUSSION.
C
C   REQD. IMSL ROUTINES - LEQT1F,LUDATN,LUELMN,UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE LINV1F (A,N,IA,AINV,IDGT,WKAREA,IER)
C
      REAL               A(IA,N),AINV(IA,N),WKAREA(1),ZERO,ONE
      DATA               ZERO/0.0/,ONE/1.0/
C                                  FIRST EXECUTABLE STATEMENT
      IER=0
      DO 10 I=1,N
         DO 5 J=1,N
            AINV(I,J) = ZERO
    5    CONTINUE
         AINV(I,I) = ONE
   10 CONTINUE
      CALL LEQT1F (A,N,N,IA,AINV,IDGT,WKAREA,IER)
      IF (IER .EQ. 0) GO TO 9005
 9000 CONTINUE
      CALL UERTST (IER,6HLINV1F)
 9005 RETURN
      END
