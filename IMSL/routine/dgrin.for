C   IMSL ROUTINE NAME   - DGRIN
C
C-----------------------------------------------------------------------
C
C   LATEST REVISION     - JANUARY 1, 1978
C
C   PURPOSE             - NUCLEUS CALLED ONLY BY IMSL SUBROUTINE DGEAR
C
C   REQD. IMSL ROUTINES - NONE REQUIRED
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C-----------------------------------------------------------------------
C
      SUBROUTINE DGRIN  (TOUT,Y,N0,Y0)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N0
      REAL               TOUT,Y0(N0),Y(N0,1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            NC,MFC,KFLAG,I,L,J,JSTART,NSQ,NQUSED,NSTEP,
     1                   NFE,NJE,NPW,NERROR,NSAVE1,NSAVE2,NEQUIL,NY,
     2                   IDUMMY(23)
      REAL               SDUMMY(4)
      REAL               T,H,HMIN,HMAX,EPSC,UROUND,EPSJ,HUSED,S,S1,
     1                   DUMMY(40)
      COMMON /GEAR/      T,H,HMIN,HMAX,EPSC,UROUND,EPSJ,HUSED,DUMMY,
     1                   SDUMMY,NC,MFC,KFLAG,JSTART,NSQ,NQUSED,NSTEP,
     2                   NFE,NJE,NPW,NERROR,NSAVE1,NSAVE2,NEQUIL,NY,
     3                   IDUMMY
C                                  FIRST EXECUTABLE STATEMENT
      DO 5 I = 1,NC
         Y0(I) = Y(I,1)
    5 CONTINUE
C                                  THIS SUBROUTINE COMPUTES INTERPOLATED
C                                    VALUES OF THE DEPENDENT VARIABLE
C                                    Y AND STORES THEM IN Y0. THE
C                                    INTERPOLATION IS TO THE
C                                    POINT T = TOUT, AND USES THE
C                                    NORDSIECK HISTORY ARRAY Y, AS
C                                    FOLLOWS..
C                                               NQ
C                                    Y0(I)  =  SUM  Y(I,J+1)*S**J ,
C                                              J=0
C                                    WHERE S = -(T-TOUT)/H.
      L = JSTART + 1
      S = (TOUT - T)/H
      S1 = 1.0
      DO 15 J = 2,L
         S1 = S1*S
         DO 10 I = 1,NC
            Y0(I) = Y0(I) + S1*Y(I,J)
   10    CONTINUE
   15 CONTINUE
      RETURN
      END
