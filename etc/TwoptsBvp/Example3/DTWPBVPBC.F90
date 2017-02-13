! Driver for TWPBVPC with all the test problems
!
!  version of July 2006
! 
!  subroutine  used by TWPBVPL
!     fsub, dfsub, gsub, dgsub, initu
!
! 
PROGRAM  DTWPBVPBC
!
   USE Base, ONLY: i4, r8
   USE BVPC, ONLY: TWPBVPC_INIT
!
   IMPLICIT NONE
!
   REAL(r8) :: es
   INTEGER(i4) :: iprob, nmsh
   LOGICAL :: bsucc
!
   WRITE(*,'(/,A,/)') 'PROGRAM DTWPBVPBC '
!      
   WRITE(*,'(A)') 'Input the problem number (0 to 35): '
   WRITE(*,'(A)') '(0 = batch process - save to ''batchc.res'')' 
   READ(*,*) iprob 
! 
   CALL TWPBVPC_INIT
!
   IF (iprob == 0) THEN
      CALL RUNBATCH()
      STOP
   ENDIF
!      
   IF ((iprob > 35) .OR. (iprob < 0)) THEN
      WRITE(*,*) ' The problem number is incorrect'
      STOP
   ENDIF   
!    
   CALL RUNPROB(iprob, 1.0D-8, 1.0D-8, bsucc, nmsh, .FALSE., es)
!
CONTAINS
!
   SUBROUTINE RUNBATCH()
!
      IMPLICIT NONE
!
      REAL(r8) :: tols(3), eps, es
      INTEGER(i4) :: iprob, j, nmsh
      LOGICAL :: bsucc
!
      OPEN(UNIT=26, FILE='batchc.res', STATUS='UNKNOWN')
!
      tols(1) = 1.0D-4
      tols(2) = 1.0D-6
      tols(3) = 1.0D-8
!
      DO iprob = 1, 35
         WRITE(26,*) 'NP'
         DO j = 1, 3
            eps = 1.0D0
            DO WHILE (.TRUE.)
               CALL RUNPROB(iprob, eps, tols(j), bsucc, nmsh, .TRUE., es)
               IF (bsucc) THEN
2609              FORMAT(1x, I2, A, I10, A, D10.3, A, D10.3, A, D10.3)
                  WRITE(26,2609) iprob, ',', nmsh, ',', eps, ',', tols(j), ',', es
               ELSE
                  EXIT
               END IF
               eps = eps/1.0D1
               IF (eps < 1.0D-6) EXIT
            END DO
         END DO
      END DO
!
      RETURN
   END SUBROUTINE RUNBATCH
!!!
!!!
   SUBROUTINE RUNPROB(iprob, epsin, etol, bsucc, nmsh, batch, es)
!
      USE BVPShared, ONLY: pdebug, use_c, comp_c, uval0,  &
                           nminit, iprint, idum
      USE BVPC,      ONLY: TWPBVPC
!
      IMPLICIT NONE
!
      INTEGER(i4), INTENT(IN)    :: iprob
      INTEGER(i4), INTENT(INOUT) :: nmsh
      REAL(r8),    INTENT(IN)    :: epsin, etol
      REAL(r8),    INTENT(INOUT) :: es
      LOGICAL,     INTENT(IN)    :: batch
      LOGICAL,     INTENT(INOUT) :: bsucc
!
      INTEGER(i4), PARAMETER :: LWRKFL = 2500000
      INTEGER(i4), PARAMETER :: LWRKIN = 500000
      INTEGER(i4), PARAMETER :: NUDIM = 6
      INTEGER(i4), PARAMETER :: NXXDIM = 100000
      INTEGER(i4), PARAMETER :: LISERIES = 100
!
      INTEGER(i4) :: ntol, ncomp, nlbc, nfxpnt, nmax, ind
      INTEGER(i4) :: iflbvp, iseries, indnms, i, nmshacc, k, j
! 
      REAL(r8) :: pi, eps, aleft, aright, ckappa1, gamma1
      REAL(r8) :: ckappa, xpt, errmax, xptold, exx
!
      REAL(r8) :: wrk(LWRKFL)
      REAL(r8) :: u(NUDIM,NXXDIM), xx(NXXDIM)
      REAL(r8) :: uacc(NUDIM,NXXDIM), xxacc(NXXDIM)      
      REAL(r8) :: tol(NUDIM), fixpnt(1), fixpntacc(NXXDIM)
      REAL(r8) :: rpar(2)
      INTEGER(i4) :: iwrk(LWRKIN)
      INTEGER(i4) :: ltol(NUDIM)
      INTEGER(i4) :: ipar(1)
      LOGICAL :: linear, giveu, givmsh
!!!
!
      bsucc = .FALSE.
      pi = 4.0D0*DATAN(1.0D0)
! 
!     If you do not like to use the conditioning in the mesh selection
!     use_c  = .false.
!     If you do not like to compute the conditioning parameters
!     comp_c = .false.     
!
      use_c  =  .true.
      comp_c =  .true.
!
      eps = epsin
!     take the reciprocal if we are dealing with problem 34
      IF (iprob == 34) eps = 1.D0 / epsin
      tol(1) = etol
!
      iprint = -1
!
      IF (.NOT. batch) THEN      
         WRITE(*,'(A)') 'Input epsilon '
         READ(*,*) eps
         WRITE(*,'(A)') 'Input tolerance '
         READ(*,*) TOL(1)
      END IF
!
      DO ind = 1, LWRKFL
         wrk(ind) = 0d0
      ENDDO
      DO ind = 1, LWRKIN
         iwrk(ind) = 0
      ENDDO
!
      WRITE(6,*) 'Problem =', iprob
      WRITE(6,*) 'Tol = ',tol(1),' eps = ',eps
!
      IF (iprob == 33) THEN
         ntol = 6
      ELSE
         ntol = 2
      END IF 
!
      ipar(1) = iprob
      rpar(1) = eps
      rpar(2) = pi
!
      CALL INITIAL(ncomp, aleft, aright, ntol, tol, ltol, rpar, ipar)
!
      IF (iprob <= 30 .OR. iprob >= 34) THEN
         nlbc = 1
      ELSE IF (iprob <= 32) THEN
         nlbc = 2
      ELSE IF (iprob == 33) THEN
         nlbc = 3
      END IF
!
      nmsh   = 0
      nfxpnt = 0
! The initial approximation is equal to uval0      
      uval0  = 0.D0     
      IF (iprob == 24) uval0 = 0.5D0
      IF (iprob == 25) uval0 = 0.5D0
      IF (iprob == 27) uval0 = 0.5D0
      IF (iprob == 28) uval0 = 0.5D0
      IF (iprob == 29) uval0 = 0.5D0
! 
      IF (iprob < 19 .OR. iprob > 34) THEN
         linear = .TRUE.
      ELSE
         linear = .FALSE.
      ENDIF
!
      giveu  = .FALSE.
      givmsh = .FALSE.
      pdebug = .FALSE.
!
!     main subroutine
!
      CALL TWPBVPC(ncomp, nlbc, aleft, aright, nfxpnt, fixpnt, ntol,  &
                   ltol, tol, linear, givmsh, giveu, nmsh, NXXDIM,  &
                   xx, NUDIM, u, nmax, LWRKFL, wrk, LWRKIN, iwrk,  &
                   ckappa1, gamma1, ckappa, rpar, ipar, iflbvp)
!
!     SCMODIFIED     
!     CHECK THE NUMBER OF MESHPOINTS
!
      IF (nmsh > NXXDIM) THEN
         WRITE(*,*) 'Assertion Failure, max meshpoints bypassed!'
         WRITE(*,*) 'Returning to avoid corruption...'
         RETURN
      END IF
!
!     NOW WE WANT TO MEASURE THE DIFFERENCE BETWEEN SOLVED MESH
!     AND A "MORE ACCURATE" MESH      
!
      DO i = 2, nmsh-1
         fixpntacc(2*(i-1)-1) = (xx(i-1)+xx(i))/2.0D0
         fixpntacc(2*(i-1)) = xx(i)
      END DO
      fixpntacc(2*(nmsh-2)+1) = (xx(nmsh-1)+xx(nmsh))/2.0
      nfxpnt = 2*(nmsh-2)+1
      nmshacc = 0
! 
      tol(1) = eps/1.0D2
! 
      DO ind = 1, LWRKFL
         wrk(ind) = 0d0
      ENDDO
      DO ind = 1, LWRKIN
         iwrk(ind) = 0
      ENDDO
!
      CALL TWPBVPC(ncomp, nlbc, aleft, aright, nfxpnt, fixpntacc, ntol,  &
                   ltol, tol, linear, givmsh, giveu, nmshacc, NXXDIM,  &
                   xxacc, NUDIM, uacc, nmax, LWRKFL, wrk, LWRKIN, iwrk,  &
                   ckappa1, gamma1, ckappa, rpar, ipar, iflbvp)
!
      IF (nmshacc > NXXDIM) THEN
         WRITE(*,*) 'Assertion Failure, max meshpoints bypassed!'
         WRITE(*,*) 'Returning to avoid corruption...'
         RETURN
      END IF     
      IF (nmshacc == nmsh) THEN
         WRITE(*,*) 'Assertion Failure, ACC mesh not good enough!'
         RETURN
      END IF
      IF (nmshacc == 0) THEN
         WRITE(*,*) 'Computation of ''accurate'' mesh failed.'
         RETURN
      END IF
!
!.... CHECK ERRORS 
!
      xpt = aleft
      es = TNORM(u, uacc, 1, 1, nlbc*2, NUDIM)
      es = es + TNORM(u, uacc, nmsh, nmshacc, nlbc*2, NUDIM)
!
!     LOOP OVER THE OTHER MESH POINTS      
!
      k = 2
      DO i = 2, nmsh-1
         DO WHILE(DABS(xx(i)-xxacc(k)) > 1.0D-10)
            k = k + 1
            IF (k >= nmshacc) THEN
               WRITE(*, *) 'Assertion failure, could not check all points!'
               RETURN
            END IF            
        END DO
        es= es + TNORM(u, uacc, i, k, nlbc*2, NUDIM)
      END DO
!      
      es = es/nmsh
      errmax = 1.0D-16
!      
      IF (.NOT. batch) THEN
         IF (iprob <= 21 .OR. iprob == 35) THEN
            WRITE(*,*) '    X             U(X)       EXACT      ERROR  ' 
         ELSE
            WRITE(*,*) '    X             U(X)         '
         END IF 
      END IF
!
      DO j = 1, nmsh
         xptold = xpt
         xpt = xx(j)
         IF (iprob <= 21 .OR. iprob == 35) THEN
            IF (j == 1 .OR. j == nmsh) THEN 
               exx = u(1,j)
            ELSE  
               exx = EXACT(xpt, rpar, ipar)
            END IF   
            errmax = DMAX1(errmax,DABS((u(1,j)-exx)/(1d0+DABS(exx))))
            IF (.NOT. batch) THEN
               WRITE(6,109) xpt,u(1,j),exx,DABS((u(1,j)-exx))/(1d0+DABS(exx))
            END IF
109         FORMAT(1x,  4(D12.3))
         ELSE 
108         FORMAT(1x,  2(D12.3))   
            IF (.NOT. batch) THEN
               WRITE(6,108) xpt,u(1,j)
            END IF 
            errmax = 0.D0
         END IF
      END DO
!
      IF (.NOT. batch) THEN  
         WRITE(*,'(/,/)')
         IF (iprob <= 21 .OR. iprob == 35) THEN
            IF (comp_c) THEN
               WRITE(6,*) '   eps      nmsh    error      kappa1      gamma1      kappa   '
               WRITE(6,103) eps, nmsh, errmax, ckappa1, gamma1, ckappa
            ELSE
               WRITE(6,*) '   eps     nmsh      error   '
               WRITE(6,103) eps, nmsh, errmax
            END IF
         ELSE
            IF (comp_c) THEN
               WRITE(6,*) '   eps     nmsh     kappa1      gamma1       kappa  '
               WRITE(6,103) eps, nmsh,  ckappa1, gamma1, ckappa
            ELSE
               WRITE(6,*) '   eps     nmsh '
               WRITE(6,103) eps, nmsh
            END IF
         END IF
         WRITE(6,*)'------------------------------------------------'
      END IF
!
100   FORMAT(/,/)
103   FORMAT(1x, D10.3, I8, 5(D12.3))
! 
      IF (iflbvp == 0) bsucc = .TRUE.
      IF (.not. batch) WRITE(*,*) 'Actual accuracy = ', es
!
      RETURN
   END SUBROUTINE RUNPROB
!!!
!!!
   FUNCTION TNORM(a, b, i1, i2, ncmp, nudim) RESULT(retval)
!
      IMPLICIT NONE
!
      REAL(r8),    INTENT(IN) :: a(nudim,*), b(nudim,*)
      INTEGER(i4), INTENT(IN) :: i1, i2, ncmp, nudim
!
      REAL(r8) :: den, retval
      INTEGER(i4) :: i
!!!
!
      retval = 0.D0
      DO i = 1, ncmp
         den = DABS(b(i,i2))
         den = MAX(den, 1.0D0)
         retval = retval + (a(i,i1)-b(i,i2))*(a(i,i1)-b(i,i2))/(den*den)
      END DO
! 
      retval = DSQRT(retval)
!
      RETURN      
   END FUNCTION TNORM
!!!
!!!
   SUBROUTINE INITIAL(ncomp,aleft,aright,ntol,tol,ltol,rpar,ipar)
!
      IMPLICIT NONE
!
      INTEGER(i4), INTENT(INOUT) :: ncomp, ntol
      REAL(r8),    INTENT(INOUT) :: aleft, aright
      REAL(r8),    INTENT(INOUT) :: tol(*), rpar(*)
      INTEGER(i4), INTENT(INOUT) :: ltol(*), ipar(*)
!
      REAL(r8) :: eps, pi, xx
      INTEGER(i4) :: iprob
      INTEGER(i4) :: j
!!!
!
      iprob = ipar(1)
      eps   = rpar(1)
      pi    = rpar(2)
!
      DO j = 1, ntol
         ltol(j) = j
         tol(j) = tol(1)
      END DO
!
      IF (iprob <= 30 .OR. iprob >= 34) THEN
         ncomp = 2
      ELSE IF (iprob <= 32) THEN
         ncomp = 4
      ELSE 
         ncomp = 6
      ENDIF
!
      IF (iprob ==  1 .OR. iprob ==  2 .OR. iprob == 8 .OR.  &
          iprob == 16 .OR. (iprob > 18 .AND. iprob <= 34)) THEN
         aleft  = 0.D0
         aright = 1.0D0
      ELSE IF (iprob == 17) THEN
         aleft = -0.1D0
         aright = 0.1D0
      ELSE IF (iprob == 18) THEN
         aleft = 0.D0
         aright = 0.25D0
      ELSE 
         aleft = -1.0D0
         aright = 1.0D0
      END IF
!
      RETURN
   END SUBROUTINE INITIAL
!!!
!!!
   FUNCTION EXACT(x, rpar, ipar) RESULT(retval)
!
      IMPLICIT NONE
!
      REAL(r8) :: x, rpar(*)
      INTEGER(i4) :: ipar(*) 
!
      REAL(r8) :: retval
      REAL(r8) :: eps, pi, pix, sqep, aa, bb, cc, dd, za, z1, xx
      INTEGER(i4) :: iprob
!!!
!
      iprob = ipar(1)
      eps   = rpar(1)
      pi    = rpar(2)
! 
      IF (iprob == 1) THEN
         retval = (DEXP(-x/DSQRT(eps))-DEXP(-(2.D0-x)/DSQRT(eps)))  &
                                    /(1.D0-DEXP(-2.D0/DSQRT(eps)))
      ELSE IF (iprob == 2) THEN
         retval = (1.D0-DEXP((x-1.D0)/eps))/(1.D0-DEXP(-1.D0/eps))
      ELSE IF (iprob == 3) THEN
         retval = DCOS(pi*x)
      ELSE IF (iprob == 4) THEN
         aa = 0.D0
         bb = (1.0D0+eps)*(1.0D0+x)/eps
         IF (bb < 100.0D0) aa=DEXP(-bb)
         retval = DEXP(x-1.0D0)+aa
      ELSE IF (iprob == 5) THEN
         retval = DCOS(pi*x)
      ELSE IF (iprob == 6) THEN
         pix = pi*x
         sqep = DSQRT(2.D0*eps)
         retval = DCOS(pix)+ERF(x/sqep)/ERF(1.d0/sqep)
      ELSE IF (iprob == 7) THEN
         cc = 0.0D0
         dd = 0.0D0
         aa = x*x/(2.0D0*eps)
         bb = 1.0D0/(2.0D0*eps)
         IF (aa < 100.0D0) cc = DEXP(-aa)
         IF (bb < 100.0D0) dd = DEXP(-bb)
         retval = DCOS(pi*x)+x+(x*ERF(x/SQRT(2.0D0*eps))  &
                                       +SQRT(2.0D0*eps/pi)*cc)/  &
                  (ERF(1.0D0/SQRT(2.0D0*eps))+SQRT(2.0D0*eps/pi)*dd)
      ELSE IF (iprob == 8) THEN
         aa = 0.0D0
         cc = 0.0D0
         bb = x/eps
         dd = 1.0D0/eps
         IF (bb < 100.0D0) aa=DEXP(-bb)
         IF (dd < 100.0D0) cc=DEXP(-dd)
         retval = (2.0D0-cc-aa)/(1.0D0-cc)
      ELSE IF (iprob == 9) THEN
         retval = 1.0D0/(x*x+eps)
      ELSE IF (iprob == 10) THEN
         retval = 1.0D0+ERF(x/SQRT(2.0D0*eps))/ERF(1.0D0/SQRT(2.0D0*eps))
      ELSE IF (iprob == 11) THEN
         retval = DCOS(pi*x)
      ELSE IF (iprob == 12) THEN
         retval = DCOS(pi*x)+DEXP(-(1.0D0-x)/DSQRT(eps))
      ELSE IF (iprob == 13) THEN
         retval = DCOS(PI*X)+DEXP(-(1.0D0+X)/DSQRT(eps))
      ELSE IF (iprob == 14) THEN
         retval = DCOS(PI*X)+DEXP(-(1.0D0+X)/DSQRT(eps))+  &
                                        DEXP(-(1.0D0-X)/DSQRT(eps))
      ELSE IF (iprob == 15) THEN
         aa = 0.0D0
         bb = 0.0D0
         IF ((1.0D0-x)/eps < 100.0D0) bb = DEXP((x-1.0D0)/eps)
         IF ((1.0D0+x)/eps < 100.0D0) aa = DEXP(-(x+1.0D0)/eps)
         za = 1.0D0/(eps**(1.0D0/3.0D0))*x
         z1 = 1.0D0/(eps**(1.0D0/3.0D0))
      ELSE IF (iprob == 16) THEN
         retval = DSIN(pi*x/(2.0D0*eps))
      ELSE IF (iprob == 17) THEN
         retval = X/SQRT(eps+x*x)
      ELSE IF (iprob == 18) THEN
         retval = DEXP(-X/eps)
      ELSE IF (iprob == 19) THEN
         retval = -DLOG((1.D0+DCOS(pi*x/2.D0))*  &
                                  (1.D0-0.5D0*DEXP(-x/(2.D0*eps))))
      ELSE IF (iprob == 20) THEN
         xx = (x-0.745D0)/eps
         IF (xx > 0.D0) THEN
            retval = 1.D0+eps*(xx+DLOG((1.D0+DEXP(-2.D0*xx))/2.D0))
         ELSE
            retval = 1.D0+eps*(-xx+DLOG((1.D0+DEXP(2.D0*xx))/2.D0))
         END IF
      ELSE IF (iprob == 21) THEN
         retval = DEXP(-x/eps)
      ELSE IF (iprob == 35) THEN
         IF (x == -1.0D+0) THEN
            retval = 1.0D+0
         ELSE IF (x == 1.D+0) THEN
            retval = 2.0D+0
         ELSE IF (x <= -0.9999D+0) THEN
            retval = -0.5D+0 + 3.0D+0/2.0D+0*EXP(-(x+1)/eps)
         ELSE IF (x >= 0.9999D+0) THEN
            retval = 0.5D+0 + 3.0D+0/2.0D+0*EXP(-(1-x)/eps)
         ELSE
            retval = x/2.0D+0
         END IF 
      END IF
!
      RETURN
   END FUNCTION EXACT
!
END
