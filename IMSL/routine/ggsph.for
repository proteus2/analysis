C   IMSL ROUTINE NAME   - GGSPH
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JUNE 1, 1980
C
C   PURPOSE             - GENERATION OF UNIFORM RANDOM DEVIATES FROM
C                           THE SURFACE OF THE UNIT SPHERE IN 3 OR 4
C                           SPACE
C
C   USAGE               - CALL GGSPH (DSEED,NR,IOPT,IZ,Z,IER)
C
C   ARGUMENTS    DSEED  - INPUT/OUTPUT DOUBLE PRECISION VARIABLE
C                           ASSIGNED AN INTEGER VALUE IN THE
C                           EXCLUSIVE RANGE (1.D0, 2147483647.D0).
C                           DSEED IS REPLACED BY A NEW VALUE TO BE
C                           USED IN A SUBSEQUENT CALL.
C                NR     - INPUT.  NUMBER OF POINTS ON THE (3 OR 4)
C                           SPHERE TO BE GENERATED.
C                IOPT   - INPUT.  OPTION PARAMETER.
C                           IOPT = 3 CAUSES GENERATION OF POINTS ON THE
C                             3-SPHERE.
C                           IOPT = 4 CAUSES 4 SPHERE POINT GENERATION.
C                           IF IOPT IS NOT 3 OR 4, A WARNING IS GIVEN
C                             (SEE DESCRIPTION OF PARAMETER IER) AND
C                             3-SPHERE POINTS ARE GENERATED.
C                IZ     - INPUT ROW DIMENSION OF MATRIX Z EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
C                           CALLING PROGRAM.
C                Z      - OUTPUT MATRIX OF DIMENSION (NR,3) OR (NR,4).
C                           Z(1,1), Z(1,2), Z(1,3) AND SUBSEQUENT GROUPS
C                           OF 3 ELEMENTS OF Z, CONTAIN THE COORDINATES
C                           OF THE VARIOUS POINTS IN 3 SPACE, IF IOPT IS
C                           NOT EQUAL TO 4.
C                           SIMILARLY, Z(1,1), Z(1,2), Z(1,3), Z(1,4)
C                           AND SUBSEQUENT GROUPS OF 4 ELEMENTS OF Z,
C                           HOLD THE COORDINATES OF THE POINTS IN 4
C                           SPACE IF IOPT EQUALS 4.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                           WARNING ERROR
C                             IER = 33 INDICATES IOPT WAS NOT EQUAL TO
C                               3 OR 4.  3-SPHERE POINTS WERE GENERATED.
C
C   REQD. IMSL ROUTINES - GGUBS,UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE GGSPH   (DSEED,NR,IOPT,IZ,Z,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            NR,IOPT,IZ,IER
      REAL               Z(IZ,1)
      DOUBLE PRECISION   DSEED
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            J,K,K1
      REAL               R(4),U(4),S(4),S2
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
      DO 20 J=1,NR
         K = 1
         GO TO 10
    5    K = 3
   10    CALL GGUBS(DSEED,2,R(K))
         K1 = K+1
         U(K) = R(K)+R(K)-1.0
         U(K1) = R(K1)+R(K1)-1.0
         S(K) = U(K)*U(K)+U(K1)*U(K1)
         IF(S(K) .GE. 1.) GO TO 10
         IF(K .EQ. 3) GO TO 15
         IF(IOPT .EQ. 4) GO TO 5
         S2 = SQRT(1.0-S(1))
         Z(J,1) = (U(1)+U(1))*S2
         Z(J,2) = (U(2)+U(2))*S2
         Z(J,3) = 1.0-S(1)-S(1)
         GO TO 20
   15    Z(J,1) = U(1)
         Z(J,2) = U(2)
         S2 = SQRT((1.0-S(1))/S(3))
         Z(J,3) = U(3)*S2
         Z(J,4) = U(4)*S2
   20 CONTINUE
      IF(IOPT .EQ. 3 .OR. IOPT .EQ. 4) GO TO 9005
 9000 CONTINUE
      IER = 33
      CALL UERTST(IER,'GGSPH ')
 9005 RETURN
      END
