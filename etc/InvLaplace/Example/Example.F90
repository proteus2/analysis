PROGRAM Example
!
!  This program should be compiled using gfortran 4.8 or higher
!  (ifort 17 or higher) because ATAN for complex argument and 
!  Bessel function BESJ0 are not implemented in lower-version compilers.
!  
   USE Base,       ONLY: i4, r4, r8, r16
   USE InvLaplace, ONLY: PI, LINVSIDI
   USE InvLaplaceExtern, ONLY: iexample, CHOICE, FS, FT
!
   IMPLICIT NONE
!
   REAL(r8), DIMENSION(:), ALLOCATABLE :: t, y
   REAL(r8) :: a0, b0, h, c, ftrue, feval, fdiff, maxdiff, tx
   REAL(r8) :: time1, time2
   INTEGER(i4) :: i, m, p, ierr
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
      CALL LINVSIDI(p, c, m, t, y, ierr)
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
END
