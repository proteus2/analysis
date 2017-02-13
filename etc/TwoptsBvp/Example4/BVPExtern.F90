MODULE BVPExtern
!
   USE Base,      ONLY: i4, r8
   USE BVPShared, ONLY: pdebug, use_c, comp_c,  &
                        uval0, nminit, iprint, idum,  &
                        MTLOAD
!
   IMPLICIT NONE
!
   PRIVATE
   PUBLIC  :: INITU, FSUB, DFSUB, GSUB, DGSUB
!
CONTAINS
!
   SUBROUTINE INITU(ncomp, nmsh, xx, nudim, u,rpar,ipar)
!
      USE Base,      ONLY: i4, r8
      USE BVPShared, ONLY: pdebug, use_c, comp_c,  &
                           uval0, nminit, iprint, idum,  &
                           MTLOAD
!
      IMPLICIT NONE
!
      INTEGER(i4), INTENT(IN) :: ncomp, nudim
      INTEGER(i4), INTENT(INOUT) :: nmsh
      REAL(r8),    INTENT(INOUT) :: xx(*), u(ncomp,*), rpar(*)
      INTEGER(i4), INTENT(INOUT) :: ipar(*)
!
      INTEGER(i4) :: iprob
!!!
!
!     This routine must be provided to reset u after re-meshing 
!     for linear problems or for nonlinear problems
!     when interpolation of the old solution is not used.
!
      IF (iprint /= -1) WRITE(6,99) uval0
99    FORMAT(1H ,'initu, uval0',1pd15.5)
!
      IF (iprob == 33) THEN
         CALL MTLOAD33(ncomp, nmsh, xx, nudim, u)
      ELSE
!        This version sets all elements of u to the constant uval0.
         CALL MTLOAD(ncomp, nmsh, uval0, nudim, u)
      END IF
!
      RETURN
   END SUBROUTINE INITU
!!!
!!!
   SUBROUTINE MTLOAD33(nrow, ncol, xx, nrowx, xmat)
!
      IMPLICIT NONE
!
      INTEGER(i4), INTENT(IN) :: nrow, nrowx
      INTEGER(i4), INTENT(INOUT) :: ncol
      REAL(r8),    INTENT(INOUT) :: xx(*), xmat(nrowx,ncol) 
!
      INTEGER(i4) :: i, j
!
!     for problem 33 the initial condition is set to 
!     xmat(1,:) = 2 x  - 1
!     xmat(2,:) = 2
!     xmat(3:end,:) 0
!!!
!
      IF (nrow <= 0 .OR. ncol <= 0) RETURN
!
      DO j = 1, ncol
         DO i = 1, nrow
!     SCMODIFIED: from
!     xmat(i,j) = 2.D0 * xx(i-1)-1.D0
            IF (i == 1) THEN
               xmat(i,j) = 2.D0 * xx(j)-1.D0
            ELSE IF (i == 2) THEN
               xmat(i,j) = 2.D0
            ELSE
               xmat(i,j) = 0.D0
            END IF
         END DO
      END DO
!
      RETURN
   END SUBROUTINE MTLOAD33
!!!
!!!
   SUBROUTINE FSUB(ncomp, x, z, f, rpar, ipar)
!
      USE Base, ONLY: i4, r8
!
      IMPLICIT NONE
!
      INTEGER(i4), INTENT(IN)    :: ncomp
      REAL(r8),    INTENT(IN)    :: x, z(*)
      REAL(r8),    INTENT(INOUT) :: f(*), rpar(*)
      INTEGER(i4), INTENT(INOUT) :: ipar(*)
!
      REAL(r8) :: eps, pi, pix, ax, apx, ga
      INTEGER(i4) :: iprob
!!!
!
      iprob = ipar(1)
      eps   = rpar(1)
      pi    = rpar(2)
!
      IF (iprob == 1) THEN
         f(1) = z(2)
         f(2) = z(1) / eps
      ELSE IF (iprob == 2) THEN
         f(1) = z(2)
         f(2) = z(2) / eps
      ELSE IF (iprob == 3) THEN
         f(1) = z(2)
         f(2) = (-(2.0D0+DCOS(pi*x))*z(2)+z(1)-(1.0d0+eps*pi*pi)*dcos(pi*x)  &
                                    -(2.0d0+DCOS(pi*x))*pi*sin(pi*x))/eps
      ELSE IF (iprob == 4) THEN
         f(1) = z(2)
         f(2) = ((1.0D0+eps)*z(1)-z(2)) / eps
      ELSE IF (iprob == 5) THEN
         pix  = pi * x
         f(1) = z(2)
         f(2) = (z(1)+x*z(2)-(1.0D0+eps*pi*pi)*dcos(pix)+x*pi*dsin(pix))/eps
      ELSE IF (iprob == 6) THEN
         pix = pi * x
         f(1) = z(2)
         f(2) = (-x*z(2)-eps*pi*pi*DCOS(pix)-pi*X*DSIN(pix))/eps
      ELSE IF (iprob == 7) THEN
         pix  = pi * x
         f(1) = z(2)
         f(2) = (-x*z(2)+z(1)-(1.0D0+eps*pi*pi)*DCOS(pix)-pix*DSIN(pix))/eps
      ELSE IF (iprob == 8) THEN
         f(1) = z(2)
         f(2) = -z(2)/eps
      ELSE IF (iprob == 9) THEN
         f(1) = z(2)
         f(2) = -(4.0D0*x*z(2)+2.0D0*z(1))/(eps+x*x)
      ELSE IF (iprob == 10) THEN
         f(1) = z(2)
         f(2) = -x*z(2)/eps
      ELSE IF (iprob == 11) THEN
         pix = pi*x
         f(1) = z(2)
         f(2) = (z(1)-eps*pi*pi*DCOS(pix)-DCOS(pix))/(eps)
      ELSE IF (iprob == 12) THEN
         pix = pi*x
         f(1) = z(2)
         f(2) = (z(1)-eps*pi*pi*DCOS(pix)-DCOS(pix))/(eps)
      ELSE IF (iprob == 13) THEN
         pix = pi*x
         f(1) =  z(2)
         f(2) =  (z(1)-eps*pi*pi*DCOS(pix)-DCOS(pix))/(eps)
      ELSE IF (iprob == 14) THEN
         pix = pi * x
         f(1) = z(2)
         f(2) = (z(1)-eps*pi*pi*DCOS(pix)-DCOS(pix))/(eps)
      ELSE IF (iprob == 15) THEN
         f(1) = z(2)
         f(2) = x * z(1)/eps
      ELSE IF (iprob == 16) THEN
         f(1) = z(2)
         f(2) = -z(1)*pi*pi/(4.0D0*eps*eps)
      ELSE IF (iprob == 17) THEN
         f(1) = z(2)
         f(2) = -3.0D0*eps*z(1)/(eps+X*X)**2
      ELSE IF (iprob == 18) THEN
         f(1) = z(2)
         f(2) = -z(2)/eps
      ELSE IF (iprob == 19) THEN
         pix = pi*x
         f(1) = z(2)
         f(2) = ((pi/2.D0)*DSIN(pix/2.D0)*DEXP(2.D0*z(1))-DEXP(z(1))*z(2))/eps
      ELSE IF (iprob == 20) THEN
         f(1) = z(2)
         f(2) = (1.D0-z(2)*z(2))/eps
      ELSE IF (iprob == 21) THEN
         f(1) = z(2)
         f(2) = (z(1)*(1.D0+z(1))-DEXP((-2.D0*X)/eps))/(eps*eps)
      ELSE IF (iprob == 22) THEN
         f(1) = z(2)
         f(2) = -(z(2)+z(1)*z(1))/eps
      ELSE IF (iprob == 23) THEN
         f(1) = z(2)
         f(2) = DSINH(z(1)/eps)/eps
      ELSE IF (iprob == 24) THEN
         ax = 1.D0+x**2
         apx = 2.D0*x
         ga = 1.4D0
         f(1) = z(2)
         f(2) = (((1.D0+ga)/2.D0-eps*apx)*z(1)*z(2)-z(2)/z(1)-  &
                   (apx/ax)*(1.D0-(ga-1.D0)*z(1)**2/2.D0))/(eps*ax*z(1))
      ELSE IF (iprob.GE.25.AND.iprob.LE.30) THEN
         f(1) = z(2)
         f(2) = (z(1)*(1.D0-z(2)))/eps
      ELSE IF (iprob == 31) THEN
         f(1) = DSIN(z(2))
         f(2) = z(3)
         f(3) = -z(4)/eps
         f(4) = ((z(1)-1.D0)*DCOS(z(2))-z(3)*(1.D0/DCOS(z(2))+  &
                                                eps*z(4)*DTAN(z(2))))/eps
      ELSE IF (iprob == 32) THEN
         f(1) = z(2)
         f(2) = z(3)
         f(3) = z(4)
         f(4) = (z(2)*z(3)-z(1)*z(4))/eps
      ELSE IF (iprob == 33) THEN
         f(1) = z(2)
         f(2) = (z(1)*z(4) - z(3)*z(2))/eps
         f(3) = z(4)
         f(4) = z(5)
         f(5) = z(6)
         f(6) = (-z(3)*z(6) - z(1)*z(2))/eps
      ELSE IF (iprob == 34) THEN
         f(1) = z(2)
         f(2) = -eps*DEXP(z(1))
      ELSE IF (iprob == 35) THEN
         f(1) = z(2)
         f(2) = X*z(2)/eps-z(1)/eps
      ENDIF
!
      RETURN
   END SUBROUTINE FSUB
!!!
!!!
   SUBROUTINE DFSUB(ncomp, x, z, df, rpar, ipar)
!
      USE Base, ONLY: i4, r8
!
      IMPLICIT NONE
!
      INTEGER(i4), INTENT(IN)    :: ncomp
      REAL(r8),    INTENT(IN)    :: x, z(*)
      REAL(r8),    INTENT(INOUT) :: df(ncomp,*), rpar(*)
      INTEGER(i4), INTENT(INOUT) :: ipar(*)
!
      REAL(r8) :: eps, pi, pix, ax, apx, ga, z2cs, z2sc
      INTEGER(i4) :: iprob
!!!
!
      iprob = ipar(1)
      eps   = rpar(1)
      pi    = rpar(2)
!
      IF (iprob == 1) THEN
         df(1,1) = 0.0D0
         df(1,2) = 1.0D0
         df(2,1) = 1.D0/eps
         df(2,2) = 0.D0
      ELSE IF (iprob == 2) THEN
         df(1,1) = 0.D0
         df(1,2) = 1.0D0
         df(2,1) = 0.0D0
         df(2,2) = 1.0D0/eps
      ELSE IF (iprob == 3) THEN
         df(1,1) = 0.0D0
         df(1,2) = 1.0D0
         df(2,1) = 1.0D0/eps
         df(2,2) = -(2.0D0+DCOS(pi*x))/eps
      ELSE IF (iprob == 4) THEN
         df(1,1) = 0.0D0
         df(1,2) = 1.0D0
         df(2,1) = (1.0D0+eps)/eps
         df(2,2) = -1.0D0/eps
      ELSE IF (iprob == 5) THEN
         df(1,1) = 0.0D0
         df(1,2) = 1.0D0
         df(2,1) = 1.0D0/eps
         df(2,2) = x/eps
      ELSE IF (iprob == 6) THEN
         df(1,1) = 0.0D0
         df(1,2) = 1.0D0
         df(2,1) = 0.0D0
         df(2,2) = -x/eps
      ELSE IF (iprob == 7) THEN
         df(1,1) = 0.0D0
         df(1,2) = 1.0D0
         df(2,1) = 1.0D0/eps
         df(2,2) = -x/eps
      ELSE IF (iprob == 8) THEN
         df(1,1) = 0.0D0
         df(1,2) = 1.0D0
         df(2,1) = 0.0D0
         df(2,2) = -1.0D0/eps
      ELSE IF (iprob == 9) THEN
         df(1,1) = 0.0D0
         df(1,2) = 1.0D0
         df(2,1) = -2.0D0/(eps+x*x)
         df(2,2) = -4.0D0*x/(eps+x*x)
      ELSE IF (iprob == 10) THEN
         df(1,1) = 0.0D0
         df(1,2) = 1.0D0
         df(2,1) = 0.0D0
         df(2,2) = -x/eps
      ELSE IF (iprob == 11) THEN
         df(1,1) = 0.0D0
         df(1,2) = 1.0D0
         df(2,1) = 1.0D0/(eps)
         df(2,2) = 0.0D0
      ELSE IF (iprob == 12) THEN
         df(1,1) = 0.0D0
         df(1,2) = 1.0D0
         df(2,1) = 1.0D0/(eps)
         df(2,2) = 0.0D0
      ELSE IF (iprob == 13) THEN
         df(1,1) = 0.0D0
         df(1,2) = 1.0D0
         df(2,1) = 1.0D0/(eps)
         df(2,2) = 0.0D0
      ELSE IF (iprob == 14) THEN
         df(1,1) = 0.0D0
         df(1,2) = 1.0D0
         df(2,1) = 1.0D0/(eps)
         df(2,2) = 0.0D0
      ELSE IF (iprob == 15) THEN
         df(1,1) = 0.0D0
         df(1,2) = 1.0D0
         df(2,1) = X/eps
         df(2,2) = 0.0D0
      ELSE IF (iprob == 16) THEN
         df(1,1) = 0.0D0
         df(1,2) = 1.0D0
         df(2,1) = -pi*pi/(4.0D0*eps*eps)
         df(2,2) = 0.0D0
      ELSE IF (iprob == 17) THEN
         df(1,1) = 0.0D0
         df(1,2) = 1.0D0
         df(2,1) = -3.0D0*eps/(eps+x*x)**2
         df(2,2) = 0.0D0
      ELSE IF (iprob == 18) THEN
         df(1,1) = 0.0D0
         df(1,2) = 1.0D0
         df(2,1) = 0.0D0
         df(2,2) = -1.D0/eps
      ELSE IF (iprob == 19) THEN
         df(1,1) = 0.0D0
         df(1,2) = 1.0D0
         df(2,1) = (pi*DSIN(pi*x/2.D0)*DEXP(2.D0*z(1))-DEXP(z(1))*z(2))/eps
         df(2,2) = -DEXP(z(1))/eps
      ELSE IF (iprob == 20) THEN
         df(1,1) = 0.0D0
         df(1,2) = 1.0D0
         df(2,1) = 0.0D0
         df(2,2) = -2.D0*z(2)/eps
      ELSE IF (iprob == 21) THEN
         df(1,1) = 0.0D0
         df(1,2) = 1.0D0
         df(2,1) = (1.D0+2.D0*z(1))/(eps*eps)
         df(2,2) = 0.0D0
      ELSE IF (iprob == 22) THEN
         df(1,1) = 0.0D0
         df(1,2) = 1.0D0
         df(2,1) = -2.D0*z(1)/eps
         df(2,2) = -1.D0/eps
      ELSE IF (iprob == 23) THEN
         df(1,1) = 0.0D0
         df(1,2) = 1.0D0
         df(2,1) = DCOSH(z(1)/eps)/(eps*eps)
         df(2,2) = 0.D0
      ELSE IF (iprob == 24) THEN
         ax = 1.D0+x**2
         apx = 2.D0*x
         ga = 1.4D0
         df(1,1) = 0.0D0
         df(1,2) = 1.0D0
         df(2,1) = ((2.D0*z(2))/z(1)**3+apx/(ax*z(1)**2)+   &
                    (apx*(ga-1.D0))/(2.D0*ax))/(eps*ax)
         df(2,2) = (((1.D0+ga)/2.D0-eps*apx)*z(1)-1.D0/z(1))/(eps*ax*z(1))
      ELSE IF (iprob >= 25 .AND. iprob <= 30) THEN
         df(1,1) = 0.0D0
         df(1,2) = 1.0D0
         df(2,1) = (1.D0-z(2))/eps
         df(2,2) = -z(1)/eps
      ELSE IF (iprob >= 25 .AND. iprob <= 30) THEN
         df(1,1) = 0.0D0
         df(1,2) = 1.0D0
         df(2,1) = (1.D0-z(2))/eps
         df(2,2) = -z(1)/eps
      ELSE IF (iprob == 31) THEN
         z2cs = DCOS(z(2))
         df(1,1) = 0.0D0
         df(1,2) = z2cs
         df(1,3) = 0.0D0
         df(1,4) = 0.0D0
         df(2,1) = 0.0D0
         df(2,2) = 0.0D0
         df(2,3) = 1.0D0
         df(2,4) = 0.0D0
         df(3,1) = 0.0D0
         df(3,2) = 0.0D0
         df(3,3) = 0.0D0
         df(3,4) = -1.0D0/eps
         df(4,1) = z2cs/eps
         z2sc = 1.D0/z2cs
         df(4,2) = (-(z(1)-1.D0)*DSIN(z(2))-z(3)*z2sc*(DTAN(z(2))+  &
                                               eps*z(4)*z2sc))/eps
         df(4,3) = -(z2sc+eps*z(4)*DTAN(z(2)))/eps
         df(4,4) = (-z(3)*eps*DTAN(z(2)))/eps
      ELSE IF (iprob == 32) THEN
         df(1,1) = 0.0D0
         df(1,2) = 1.D0
         df(1,3) = 0.0D0
         df(1,4) = 0.0D0
         df(2,1) = 0.0D0
         df(2,2) = 0.0D0
         df(2,3) = 1.0D0
         df(2,4) = 0.0D0
         df(3,1) = 0.0D0
         df(3,2) = 0.0D0
         df(3,3) = 0.0D0
         df(3,4) = 1.0D0
         df(4,1) = -z(4)/eps
         df(4,2) = z(3)/eps
         df(4,3) = z(2)/eps
         df(4,4) = -z(1)/eps
      ELSE IF (iprob == 33) THEN
         df(1,2) = 1.0D0
         df(2,1) = z(4)/eps
         df(2,2) = -z(3)/eps
         df(2,3) = -z(2)/eps
         df(2,4) = z(1)/eps
         df(3,4) = 1.0D0
         df(4,5) = 1.0D0
         df(5,6) = 1.0D0
         df(6,1) = -z(2)/eps
         df(6,2) = -z(1)/eps
         df(6,3) = -z(6)/eps
         df(6,6) = -z(3)/eps
      ELSE IF (iprob == 34) THEN
         df(1,1) = 0
         df(1,2) = 1
         df(2,1) = -eps*DEXP(z(1))
         df(2,2) = 0
      ELSE IF (iprob == 35) THEN
         df(1,1) = 0.0D0
         df(1,2) = 1.0D0
         df(2,1) = -1.0D0/eps
         df(2,2) = x/eps
      END IF
!
      RETURN
   END SUBROUTINE DFSUB
!!!
!!!
   SUBROUTINE GSUB(i, ncomp, z, g, rpar, ipar)
!
      USE Base, ONLY: i4, r8
!
      IMPLICIT NONE
!
      INTEGER(i4), INTENT(IN)    :: i, ncomp
      REAL(r8),    INTENT(IN)    :: z(*)
      REAL(r8),    INTENT(INOUT) :: g, rpar(*)
      INTEGER(i4), INTENT(INOUT) :: ipar(*)
!
      REAL(r8) :: eps, pi, pix, ax, apx, ga, z2cs, z2sc, x
      INTEGER(i4) :: iprob
!!!
!
      iprob = ipar(1)
      eps   = rpar(1)
      pi    = rpar(2)
!
      IF (iprob == 1) THEN
         IF (i == 1) g = z(1)-1.D0
         IF (i == 2) g = z(1)
      ELSE IF (iprob == 2) THEN
         IF (i == 1) g = z(1)-1.D0
         IF (i == 2) g = z(1)
      ELSE IF (iprob == 3) THEN
         g = z(1)+1.0D0
      ELSE IF (iprob == 4) THEN
         IF (i == 1) g = z(1)-1.0D0-DEXP(-2.0D0)
         IF (i == 2) g = z(1)-1.0D0-DEXP(-2.0D0*(1.0D0+eps)/eps)
      ELSE IF (iprob == 5) THEN
         g = z(1)+1.0D0
      ELSE IF (iprob == 6) THEN
         IF (i == 1) g = z(1)+2.0D0
         IF (i == 2) g = z(1)
      ELSE IF (iprob == 7) THEN
         IF (i == 1) g = z(1)+1.0D0
         IF (i == 2) g = z(1)-1.0D0
      ELSE IF (iprob == 8) THEN
         IF (i == 1) g = z(1)-1.0D0
         IF (i == 2) g = z(1)-2.0D0
      ELSE IF (iprob == 9) THEN
         g = z(1)-1.0D0/(1.0D0+EPS)
      ELSE IF (iprob == 10) THEN
         IF (i == 1) g = z(1)-0.0D0
         IF (i == 2) g = z(1)-2.0D0
      ELSE IF (iprob == 11) THEN
         g = z(1)+1.0D0
      ELSE IF (iprob == 12) THEN
         IF (i == 1) g = z(1)+1.0D0
         IF (i == 2) g = z(1)
      ELSE IF (iprob == 13) THEN
         IF (i == 1) g = z(1)
         IF (i == 2) g = z(1)+1.0D0
      ELSE IF (iprob == 14) THEN
         g = z(1)
      ELSE IF (iprob == 15) THEN
         g = z(1)-1.0D0
      ELSE IF (iprob == 16) THEN
         IF (i == 1) g = z(1)
         IF (i == 2) g = z(1)-DSIN(PI/(2.0D0*eps))
      ELSE IF (iprob == 17) THEN
         IF (i == 1) g = z(1)+0.1D0/DSQRT(eps+0.01D0)
         IF (i == 2) g = z(1)-0.1D0/DSQRT(eps+0.01D0)
      ELSE IF (iprob == 18) THEN
         IF (i == 1) g = z(1)-1.D0
         IF (i == 2) g = z(1)-DEXP(-1.D0/(4.D0*EPS))
      ELSE IF (iprob == 19) THEN
         IF (i == 1) g = z(1)
         IF (i == 2) g = z(1)
      ELSE IF (iprob == 20) THEN
         IF (i == 1) THEN
            x = -0.745D0/eps
            g = z(1)-1.D0-eps*(-X+DLOG((DEXP(2.D0*X)+1.D0)/2.D0))
         ELSE
            x = 0.255D0/eps
            g = z(1)-1.D0-eps*(X+DLOG((DEXP(-2.D0*X)+1.D0)/2.D0))
         ENDIF
      ELSE IF (iprob == 21) THEN
         IF (i == 1) g = z(1)-1.D0
         IF (i == 2) g = z(1)-DEXP(-1.D0/eps)
      ELSE IF (iprob == 22) THEN
         IF (i == 1) g = z(1)
         IF (i == 2) g = z(1)-0.5D0
      ELSE IF (iprob == 23) THEN
         IF (i == 1) g = z(1)
         IF (i == 2) g = z(1)-1.D0
      ELSE IF (iprob == 24) THEN
         IF (i == 1) g = z(1)-0.9129D0
         IF (i == 2) g = z(1)-0.375D0
      ELSE IF (iprob == 25) THEN
         IF (i == 1) g = z(1)+1.D0/3.D0
         IF (i == 2) g = z(1)-1.D0/3.D0
      ELSE IF (iprob == 26) THEN
         IF (i == 1) g = z(1)-1.D0
         IF (i == 2) g = z(1)+1.D0/3.D0
      ELSE IF (iprob == 27) THEN
         IF (i == 1) g = z(1)-1.D0
         IF (i == 2) g = z(1)-1.D0/3.D0
      ELSE IF (iprob == 28) THEN
         IF (i == 1) g = z(1)-1.D0
         IF (i == 2) g = z(1)-3.D0/2.D0
      ELSE IF (iprob == 29) THEN
         IF (i == 1) g = z(1)
         IF (i == 2) g = z(1)-3.D0/2.D0
      ELSE IF (iprob == 30) THEN
         IF (i == 1) g = z(1)+7.D0/6.D0
         IF (i == 2) g = z(1)-3.D0/2.D0
      ELSE IF (iprob == 31) THEN
         IF (i == 1) g = z(1)
         IF (i == 2) g = z(3)
         IF (i == 3) g = z(1)
         IF (i == 4) g = z(3)
      ELSE IF (iprob == 32) THEN
         IF (i == 1) g = z(1)
         IF (i == 2) g = z(2)
         IF (i == 3) g = z(1)-1.D0
         IF (i == 4) g = z(2)
      ELSE IF (iprob == 33) THEN
         IF (i == 1) g = z(1)+1.D0
         IF (i == 2) g = z(3)
         IF (i == 3) g = z(4)
         IF (i == 4) g = z(1)-1.D0
         IF (i == 5) g = z(3)
         IF (i == 6) g = z(4)
      ELSE IF (iprob == 34) THEN
         g = z(1)
      ELSE IF (iprob == 35) THEN
         IF (i == 1) g = z(1)-1.0d0
         IF (i == 2) g = z(1)-2.0d0
      ENDIF
!
      RETURN
   END SUBROUTINE GSUB
!!!
!!!
   SUBROUTINE DGSUB(i, ncomp, z, dg, rpar, ipar)
!
      USE Base, ONLY: i4, r8
!
      IMPLICIT NONE
!
      INTEGER(i4), INTENT(IN)    :: i, ncomp
      REAL(r8),    INTENT(IN)    :: z(*)
      REAL(r8),    INTENT(INOUT) :: dg(*), rpar(*)
      INTEGER(i4), INTENT(INOUT) :: ipar(*)
!
      REAL(r8) :: eps, pi, pix, ax, apx, ga, z2cs, z2sc
      INTEGER(i4) :: iprob
!!!
!
      iprob = ipar(1)
      eps   = rpar(1)
      pi    = rpar(2)
!
      dg(1) = 1.D0
      dg(2) = 0.D0
      IF (iprob == 31) THEN
         dg(1) = 0.D0
         dg(2) = 0.D0
         dg(3) = 0.D0
         dg(4) = 0.D0
         IF (i == 1 .OR. i == 3) dg(1) = 1.D0
         IF (i == 2 .OR. i == 4) dg(3) = 1.D0
      ELSE IF (iprob == 32) THEN
         dg(1) = 0.D0
         dg(2) = 0.D0
         dg(3) = 0.D0
         dg(4) = 0.D0
         IF (i == 1 .OR. i == 3) dg(1) = 1.D0
         IF (i == 2 .OR. i == 4) dg(2) = 1.D0
      ELSE IF (iprob == 33) THEN
         dg(1) = 0.D0
         dg(2) = 0.D0
         dg(3) = 0.D0
         dg(4) = 0.D0
         dg(5) = 0.D0
         dg(6) = 0.D0
         IF (i == 1 .OR. i == 4) dg(1) = 1.D0
         IF (i == 2 .OR. i == 5) dg(3) = 1.D0
         IF (i == 3 .OR. i == 6) dg(4) = 1.D0
      ENDIF
!
      RETURN
   END SUBROUTINE DGSUB
!
END
