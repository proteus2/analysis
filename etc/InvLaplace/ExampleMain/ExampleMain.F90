PROGRAM ExampleMain
!
!  This program should be compiled using gfortran 4.8 or higher
!  (ifort 17 or higher) because ATAN for complex argument and 
!  Bessel function BESJ0 are not implemented in lower-version compilers.
!  
   USE Base,       ONLY: i4, r4, r8, r16
   USE InvLaplace, ONLY: PI, LINVSIDI
!
   IMPLICIT NONE
!
   REAL(r8), DIMENSION(:), ALLOCATABLE :: t, y
   REAL(r8) :: a0, b0, h, c, ftrue, feval, fdiff, maxdiff, tx
   REAL(r8) :: time1, time2
   INTEGER(i4) :: iexample, i, m, p, ierr
!
   WRITE(6,'(A)') 'Inversion of the Laplace Transform for given values of'
   WRITE(6,'(A)') 'the transformed function for complex arguments'
   WRITE(6,'(A)') '(Bromwich Integral, Algorithm of Sidi)'
!
   DO
!
      CALL CHOICE
      IF (iexample == 0) EXIT
!
      WRITE(6,'(A)') 'Enter p, c (> conv. absc.), a0, b0, m: '
      READ(*,*) p, c, a0, b0, m
!
      ALLOCATE(t(1:m))
      ALLOCATE(y(1:m))
!
      h = (b0 - a0)/(m - 1)
!
      DO i = 1, m
         tx = a0 + (i-1)*h
         t(i) = tx
      END DO 
!
      CALL CPU_TIME(time1)
      CALL LINVSIDI(FS, p, c, m, t, y, ierr)
      CALL CPU_TIME(time2)
!
      maxdiff = 0._r8
!
      WRITE(6,'(/A/)') '     t         Result                 True Solution           Error  '
      DO i = 1, m
         ftrue = FT(t(i))
         feval = y(i)
         fdiff = ABS(feval - ftrue)
         IF (maxdiff < fdiff) maxdiff = fdiff
         WRITE(6,'(E11.4,1X,E23.15,1X,E23.15,1X,E11.4)') t(i), feval, ftrue, fdiff
      END DO
      WRITE(6,'(A,E11.4)') 'Max. abs. Error: ', maxdiff
      WRITE(6,'(A,F15.4)') 'Computing time: ', time2-time1
! 
      IF (ALLOCATED(t)) DEALLOCATE(t)
      IF (ALLOCATED(y)) DEALLOCATE(y)
!
   END DO
!
CONTAINS
!
!-------------------------------------------------------------------------------
!
!  SUBROUTINE CHOICE
!
!-------------------------------------------------------------------------------
!
   SUBROUTINE CHOICE
!
      IMPLICIT NONE
!
      WRITE(6,'(A)') 'Choice of Example'
      WRITE(6,'(A)') '( 0) End'
      WRITE(6,'(A)') '( 1) F(s) = 1/(s + 0.5)                f(t) = exp(-t/2)'
      WRITE(6,'(A)') '( 2) F(s) = (s+0.2)/((s+0.2)^2 + 1)    f(t) = exp(-0.2*t)*cos(t)'
      WRITE(6,'(A)') '( 3) F(s) = 1/s                        f(t) = 1'
      WRITE(6,'(A)') '( 4) F(s) = 1/(s*s)                    f(t) = t'
      WRITE(6,'(A)') '( 5) F(s) = 1/(s*s*s)                  f(t) = t*t/2'
      WRITE(6,'(A)') '( 6) F(s) = 1/sqrt(s)                  f(t) = 1.0/sqrt(pi*t)'
      WRITE(6,'(A)') '( 7) F(s) = 1/(s*s + 1)                f(t) = sin(t)'
      WRITE(6,'(A)') '( 8) F(s) = 1/sqrt(s*s + 1)            f(t) = J0(t)'
      WRITE(6,'(A)') '( 9) F(s) = exp(-1/s)/s                f(t) = J0(2*srqt(t))'
      WRITE(6,'(A)') '(10) F(s) = 1/(s*sqrt(s))              f(t) = 2*sqrt(t/pi)'
      WRITE(6,'(A)') '(11) F(s) = log((s*s + 1)/(s*s + 4))   f(t) = 2*(cos(2*t) - cos(t))/t'
      WRITE(6,'(A)') '(12) F(s) = atan(1/s)                  f(t) = sin(t)/t'
      WRITE(6,'(A)') '(13) F(s) = exp(-s)/s                  f(t) = u(t-1)'
      WRITE(6,'(A)') '(14) F(s) = exp(-s)/(s*(s+1))          f(t) = 0 (t<=1), 1-exp(-(t-1)), (1<=t)'
      WRITE(6,'(A)') '(15) F(s) = 1/(s^2 + s + 1)            f(t) = 2/sqrt(3)*exp(-t/2)*sin(t*sqrt(3)/2)'
      WRITE(6,'(A)') '(16) F(s) = log(s)/s                   f(t) = -gamma - log(t)'
      WRITE(6,'(A)') '(17) F(s) = 1/(s - 1)                  f(t) = exp(t)'
      WRITE(6,'(A)') '(18) F(s) = exp(-4*sqrt(s))            f(t) = 2*exp(-4/t)/sqrt(pi*t^3)'
      WRITE(6,'(A)') '(19) F(s) = 1/(s^3 - 8)                f(t) = 1/12*exp(-t)(exp(3*t) - cos(sqrt(3)*t)'
      WRITE(6,'(A)') '                                                           - sqrt(3)*sin(sqrt(3)*t))'
      WRITE(6,'(A)') '(20) F(s) = 1/(s*(1+exp(s)))           f(t) = 0 (2k < t < 2k+1), 1 (2k+1 < t < 2k+2)'
      WRITE(6,'(A)') '(21) F(s) = 1/(s*(s+1))*(1+exp(-pi*s)) f(t) = |sin(t)|'
      WRITE(6,'(A)') '            /(1-exp(-pi*s))'
      WRITE(6,'(A)') '(22) F(s) = 1/((s-1)*sqrt(s))          f(t) = exp(t)*erf(sqrt(t))'
      WRITE(6,'(A)') '(23) F(s) = exp(-5*sqrt(s))/s          f(t) = 1-erf(5/(2*sqrt(t)))'
!
      DO
         WRITE(6,'(A)') 'Choose : '
         READ(*,*) iexample
         IF (iexample >= 0 .AND. iexample <= 23) EXIT
      END DO
!
      RETURN
   END SUBROUTINE CHOICE  
! 
!-------------------------------------------------------------------------------
!
!  FUNCTION FS
!
!-------------------------------------------------------------------------------
!
   FUNCTION FS(s) RESULT(ftval)
!
      IMPLICIT NONE
!
      COMPLEX(r8), INTENT(IN) :: s
!
      COMPLEX(r8) :: ftval, ai
!!!
!
      ai = (0._r8, 1._r8)
!
      SELECT CASE (iexample)
         CASE ( 1)
            ftval = 1.0_r8 / (s + 0.5_r8)
         CASE ( 2) 
            ftval = (s + 0.2_r8)/((s + 0.2_r8)*(s + 0.2_r8) + 1.0_r8)
         CASE ( 3)
            ftval = 1._r8 / s
         CASE ( 4)
            ftval = 1._r8 / (s*s)
         CASE ( 5)
            ftval = 1._r8 / (s*s*s)
         CASE ( 6)
            ftval = 1._r8 / SQRT(s)
         CASE ( 7)
            ftval = 1._r8 / (s*s + 1._r8)
         CASE ( 8)
            ftval = 1._r8 / SQRT(s*s + 1._r8)
         CASE ( 9)
            ftval = EXP(-1._r8 / s) / s
         CASE (10)
            ftval = 1._r8 / (s * SQRT(s))
         CASE (11)
            ftval = LOG((s*s + 1._r8) / (s*s + 4._r8))
         CASE (12)
            ftval = ATAN(1._r8 / s)
         CASE (13)
            ftval = EXP(-s) / s
         CASE (14)
            ftval = EXP(-s) / (s*(s+1._r8))
         CASE (15)
            ftval = 1._r8 / (s*s + s + 1._r8)
         CASE (16)
            ftval = LOG(s) / s
         CASE (17)
            ftval = 1._r8 / (s - 1._r8)
         CASE (18)
            ftval = EXP(-4._r8 * SQRT(s))
         CASE (19)
            ftval = 1._r8 / (s*s*s - 8._r8)
         CASE (20)
            ftval = 1._r8 / (s*(1._r8 + EXP(s)))
         CASE (21)
            ftval = 1._r8/(s*s + 1._r8)*(1._r8 + EXP(-PI*s))/(1.0 - EXP(-PI*s))
         CASE (22)
            ftval = 1._r8/((s-1._r8)*SQRT(s))
         CASE (23)
            ftval = EXP(-5._r8*SQRT(s))/s
         CASE DEFAULT
            ftval = 0._r8 + ai*0._r8
      END SELECT
!
      RETURN
   END FUNCTION FS
!
!-------------------------------------------------------------------------------
!
!  FUNCTION FT
!
!-------------------------------------------------------------------------------
!
   FUNCTION FT(t) RESULT(ftval)
!
      IMPLICIT NONE
!
      REAL(r8), INTENT(IN) :: t
!
      COMPLEX(r8) :: ftval, gam
!!!
!
      gam = 0.577215664901532_r8
!
      SELECT CASE (iexample)
         CASE ( 1)
            ftval = EXP(-t*0.5_r8)
         CASE ( 2)
            ftval = EXP(-0.2_r8*t) * COS(t)
         CASE ( 3)
            ftval = 1._r8
         CASE ( 4)
            ftval = t
         CASE ( 5)
            ftval = t * t / 2._r8
         CASE ( 6)
            ftval = 1._r8 / SQRT(PI*t)
         CASE ( 7)
            ftval = SIN(t)
         CASE ( 8)
            ftval = BESJ0(t)
         CASE ( 9)
            ftval = BESJ0(2._r8*SQRT(t))
         CASE (10)
            ftval = 2._r8*SQRT(t/PI);
         CASE (11) 
            ftval = FTN(t, 11)
         CASE (12)
            ftval = FTN(t, 12)
         CASE (13)
            ftval = U(t-1)
         CASE (14)
            ftval = FTN(t, 14)
         CASE (15)
            ftval = 2._r8/SQRT(3._r8)*EXP(-t/2._r8)*SIN(t*SQRT(3._r8)/2._r8)
         CASE (16)
            ftval = -gam-LOG(t)
         CASE (17)
            ftval = EXP(t)
         CASE (18)
            ftval = FTN(t, 18)
         CASE (19)
            ftval = 1._r8/12._r8*EXP(-t)*  &
                 (EXP(3._r8*t)-COS(SQRT(3._r8)*t)-SQRT(3._r8)*SIN(SQRT(3._r8)*t))
         CASE (20)
            ftval = SQWAVE(t)
         CASE (21)
            ftval = ABS(SIN(t))
         CASE (22)
            ftval = EXP(t)*ERF(SQRT(t))
         CASE (23)
            ftval = 1._r8 - ERF(5._r8/(2._r8*SQRT(t)))
         CASE DEFAULT
            ftval = 0._r8
      END SELECT
!
      RETURN
   END FUNCTION FT

!-------------------------------------------------------------------------------
!
!  FUNCTION FTN
!
!-------------------------------------------------------------------------------
!
   FUNCTION FTN(t, icase) RESULT(fval)
!
      IMPLICIT NONE
!
      REAL(r8),    INTENT(IN) :: t
      INTEGER(i4), INTENT(IN) :: icase
!
      REAL(r8) :: fval
!!!
!
      IF (icase == 11) THEN
         IF (t == 0._r8) THEN
            fval = 0._r8
         ELSE
            fval = 2._r8*(COS(2._r8*t) - COS(t))/t
         END IF
      ELSE IF (icase == 12) THEN
         IF (t == 0._r8) THEN
            fval = 1._r8
         ELSE
            fval = SIN(t) / t
         END IF
      ELSE IF (icase == 14) THEN
         IF (t <= 1._r8) THEN
            fval = 0._r8
         ELSE
            fval = 1._r8 - EXP(-(t-1._r8))
         END IF
      ELSE IF (icase == 18) THEN
         IF (ABS(t) <= 0.001_r8) THEN
            fval = 0._r8
         ELSE
            fval = 2._r8*exp(-4._r8/t) / SQRT(PI*t*t*t)
         END IF
      END IF
!
      RETURN
   END FUNCTION FTN
!
!-------------------------------------------------------------------------------
!
!  FUNCTION U
!
!-------------------------------------------------------------------------------
!
   FUNCTION U(t) RESULT(fval)
!
      IMPLICIT NONE
!
      REAL(r8), INTENT(IN) :: t
!
      REAL(r8) :: weps, fval
!!!
!
      weps = 1.0e-10_r8
!
      IF (t < -weps) THEN
         fval = 0._r8
      ELSE IF (t > weps) THEN
         fval = 1._r8
      ELSE
         fval = .5_r8
      END IF
!
      RETURN
   END FUNCTION U
!
!-------------------------------------------------------------------------------
!
!  FUNCTION SQWAVE
!
!-------------------------------------------------------------------------------
!
   FUNCTION SQWAVE(t) RESULT(fval)
!
      IMPLICIT NONE
!
      REAL(r8), INTENT(IN) :: t
!
      REAL(r8) :: fval, tfl, eps
      INTEGER(i4) :: r, tf
!!!
!
      tfl = FLOOR(t)
      eps = 1.0e-7_r8
!
      tf = INT(tfl)
      r = MOD(tf, 2)
!
      IF (ABS(t - tfl) <= eps) THEN
         fval = 0.5_r8
      ELSE IF (r == 0) THEN
         fval = 0.0_r8
      ELSE
         fval = 1.0_r8
      END IF
!
      RETURN
   END FUNCTION SQWAVE
! 
END
