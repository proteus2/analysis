MODULE BVPC
!
   USE Base,      ONLY: i4, r4, r8
   USE BVPShared, ONLY: flmin, flmax, epsmch,  &
                        pdebug, use_c, comp_c,  &
                        uval0, nminit, iprint, idum,  &
                        alp1, alp2, alp3,  &
                        bet0, bet2, bet3, bet4,  &
                        a1, b1, c1, d1, e2, f2, c2, d2, b2, a2,  &
                        p3, q3, e3, f3, c3, d3, a3, b3,  &
                        a4, p4, x4, e4, c4, &
                        a5, b5, c5, d5, e5, f5, a6, b6, c6,  &
                        ABDNRM, DASUM, DONEST,  &
                        SPRT, MPRT, COLROW, INVERSE, CRDCMP, CRSLVE,  &
                        LUFAC, LUSOL, DCOPY, DAXPY, DDOT,  &
                        DSCAL, DSWAP, IDAMAX, DLOAD,  &
                        MAXPY, MATCOP, MTLOAD, MSSQ, DSSQ
   USE BVPExtern, ONLY: INITU, FSUB, DFSUB, GSUB, DGSUB
!
   IMPLICIT NONE
!
   PRIVATE
   PUBLIC :: TWPBVPC_INIT, TWPBVPC
!
CONTAINS
!
   SUBROUTINE TWPBVPC_INIT
!
      USE BVPShared, ONLY: pdebug, use_c, comp_c, uval0,  &
                           nminit, iprint, idum
!
      IMPLICIT NONE
!
!     BLOCK DATA
!
      nminit = 16
      pdebug = .false.
      iprint = -1
      uval0 = 0.0_r8
      use_c = .true.
      comp_c = .true.
!
      RETURN
   END SUBROUTINE TWPBVPC_INIT
!
!-------------------------------------------------------------------------------
!
!  SUBROUTINE TWPBVPC
!
!  Code converted using TO_F90 by Alan Miller
!  Date: 2017-02-06  Time: 18:03:10
!
!  The subroutine TWPBVPC is intended to solve two-point boundary value problems
!
!  References:
!     Cash, J. R., and Mazzia, F., 2005:
!        A new mesh selection algorithm, based on conditioning, for two-point
!        boundary value codes. J. Comput. Appl. Math. 184, 2, 362--381.
!
!  Revision History:
!
!  03JUL2006: Added rpar and ipar in the functions
!             DFSUB, DGSUB, FSUB, GSUB
!             changed the name of the variable double in ddouble
!
!  20SEP2004: This is a modified version of twpbvp that uses the conditioning
!             in the mesh selection.
!
!  New subroutines not included in the old version:
!     CONDESTIM, ESTIMKAPPA, MONCOND, SELCONDMSH, SELCONDERRMSH, SMPSELCONDMSH
!
!  Updates subroutines:
!     BVPSOL, CONV4, FAIL4, FAIL6, CONV8, FAIL8, NEWTEQ, MSHREF, WTCHDG
!
!  Auxiliary function not used by the old version:
!     DONEST, ABDNRM
!
!  The common block algprs contains two more variable
!  logical pdebug, use_c, comp_c
!  common/algprs/ nminit, pdebug, iprint, idum, uval0, use_c, comp_c
!
!  The starting subroutine is twpbvp by J.R. Cash and M.H. Wright
!
!-------------------------------------------------------------------------------
!
   SUBROUTINE TWPBVPC(ncomp, nlbc, aleft, aright,  &
                      nfxpnt, fixpnt, ntol, ltol, tol,  &
                      linear, givmsh, giveu, nmsh,  &
                      nxxdim, xx, nudim, u, nmax,  &
                      lwrkfl, wrk, lwrkin, iwrk,  &
                      ckappa1, gamma1, ckappa, rpar, ipar, iflbvp)
!
      IMPLICIT NONE
!
      INTEGER(i4), INTENT(IN)    :: ncomp
      INTEGER(i4), INTENT(INOUT) :: nlbc
      REAL(r8),    INTENT(IN)    :: aleft
      REAL(r8),    INTENT(IN)    :: aright
      INTEGER(i4), INTENT(IN)    :: nfxpnt
      REAL(r8),    INTENT(IN)    :: fixpnt(*)
      INTEGER(i4), INTENT(IN)    :: ntol
      INTEGER(i4), INTENT(IN)    :: ltol(*)
      REAL(r8),    INTENT(IN)    :: tol(*)
      LOGICAL,     INTENT(IN)    :: linear
      LOGICAL,     INTENT(IN)    :: givmsh
      LOGICAL,     INTENT(IN)    :: giveu
      INTEGER(i4), INTENT(INOUT) :: nmsh
      INTEGER(i4), INTENT(IN)    :: nxxdim
      REAL(r8),    INTENT(INOUT) :: xx(*)
      INTEGER(i4), INTENT(IN)    :: nudim
      REAL(r8),    INTENT(INOUT) :: u(nudim,*)
      INTEGER(i4), INTENT(INOUT) :: nmax
      INTEGER(i4), INTENT(IN)    :: lwrkfl
      REAL(r8),    INTENT(INOUT) :: wrk(lwrkfl)
      INTEGER(i4), INTENT(IN)    :: lwrkin
      INTEGER(i4), INTENT(INOUT) :: iwrk(lwrkin)
      REAL(r8),    INTENT(INOUT) :: ckappa1
      REAL(r8),    INTENT(INOUT) :: gamma1
      REAL(r8),    INTENT(INOUT) :: ckappa
      REAL(r8),    INTENT(INOUT) :: rpar(*)
      INTEGER(i4), INTENT(INOUT) :: ipar(*)
      INTEGER(i4), INTENT(INOUT) :: iflbvp
!
      REAL(r8), PARAMETER :: zero = 0.0D+0
      INTEGER(i4) :: isp, iden, nmax1, nmax2
      INTEGER(i4) :: irhs, lrhs, itpblk, ltpblk, ibtblk, lbtblk, iajac, lajac
      INTEGER(i4) :: ibhold, lbhold, ichold, lchold, ifval, lfval, idef, ldef
      INTEGER(i4) :: idefex, ldefex, idef6, ldef6, idefim, ldefim, idef8, ldef8
      INTEGER(i4) :: iusve, lusve, iuold, luold, itmrhs, ltmrhs, irhtri, lrhtri
      INTEGER(i4) :: idelu, ldelu, ixmer, lxmer, irerr, lrerr, iutri, lutri
      INTEGER(i4) :: iermx, lermx, irtdc, lrtdc, ixxold, lxxold, iuint, luint
      INTEGER(i4) :: iftmp, lftmp, idgtm, ldgtm, idftm1, ldftm1, idftm2, ldftm2
      INTEGER(i4) :: itmp, ltmp, idsq, ldsq, idexr, ldexr, ietst6, letst6 
      INTEGER(i4) :: ietst8, letst8, iamg, lamg, ic1, lc1, iwrkrhs, lwrkrhs
      INTEGER(i4) :: ilast, iiref, liref, iihcom, lihcom, iipvbk, lipvbk
      INTEGER(i4) :: iipvlu, lipvlu, iisign, lisign
      INTEGER(i4) :: i

!     Check for invalid input parameters.  If any parameters are
!     invalid, exit with the flag iflbvp set to -1.
!
      iflbvp = -1
      IF (ncomp <= 0)  RETURN
      IF (nlbc < 0 .OR. nlbc > ncomp) RETURN
      IF (aleft >= aright) RETURN
!
      IF (nfxpnt < 0)  RETURN
      IF (givmsh .AND. nmsh < nfxpnt+2) RETURN
      IF (givmsh .AND. xx(1) /= aleft) RETURN
!     SCMODIFIED add an extra condition to avoid accessing xx(0)
      IF (nmsh > 0) THEN
         IF (givmsh .AND. xx(nmsh) /= aright) RETURN
      END IF
      IF (nfxpnt > 0) THEN
         IF (fixpnt(1) <= aleft) RETURN
         IF (fixpnt(nfxpnt) >= aright) RETURN
         DO i = 1, nfxpnt-1
            IF (fixpnt(i+1) <= fixpnt(i)) RETURN
         END DO
      END IF
!
      IF (ntol < 1) RETURN
      DO i = 1, ntol
         IF (ltol(i) < 0 .OR. ltol(i) > ncomp) RETURN
         IF (tol(i) <= zero) RETURN
      END DO
!
      IF (giveu .AND. .NOT. givmsh) RETURN
      IF (use_c .AND. .NOT. comp_c) RETURN
      IF (nudim <= 0) RETURN
      IF (lwrkfl <= 0 .OR. lwrkin <= 0) RETURN
!
!     Calculate maximum number of mesh points possible with the
!     given floating-point and integer workspace.
!
      isp = lwrkfl - 2*ntol - 22*ncomp - 6*ncomp*ncomp
      iden = 6*ncomp*ncomp + 22*ncomp + 3
      nmax1 = isp/iden
!
      isp = lwrkin - 2*ncomp
      nmax2 = isp/(2*ncomp+3)
!
      nmax = MIN(nmax1,nmax2)
!     nmax from workspace
      nmax = MIN(nmax, nxxdim)
!     nmax from size of u and xx
!
      IF (iprint >= 0) WRITE(6,901) nmax
901   FORMAT(1H ,'nmax from workspace =',i8)
!
      IF (nmax <= 1) RETURN
!
!     Partition floating point workspace.
!
      irhs = 1
      lrhs = ncomp*nmax
!
      itpblk = irhs + lrhs
      ltpblk = ncomp*nlbc
!
      ibtblk = itpblk + ltpblk
      lbtblk = ncomp*(ncomp - nlbc)
!
      iajac = ibtblk + lbtblk
      lajac = 2*ncomp*ncomp*nmax
!
      ibhold = iajac + lajac
      lbhold = ncomp*ncomp*nmax
!
      ichold = ibhold + lbhold
      lchold = ncomp*ncomp*nmax
!
      ifval = ichold + lchold
      lfval = ncomp*nmax
!
      idef = ifval + lfval
      ldef = ncomp*(nmax-1)
!
      idefex = idef + ldef
      ldefex = ncomp*(nmax-1)
!
!     def6 uses the same space as defexp
!
      idef6 = idefex
      ldef6 = ncomp*(nmax-1)
!
      idefim = idef6 + ldef6
      ldefim = ncomp*(nmax-1)
!
!     def8 uses the same space as defimp
!
      idef8 = idefim
      ldef8 = ncomp*(nmax-1)
!
      iusve = idef8 + ldef8
      lusve = ncomp*nmax
!
      iuold = iusve + lusve
      luold = ncomp*nmax
!
      itmrhs = iuold + luold
      ltmrhs = ncomp*nmax
!
      irhtri = itmrhs + ltmrhs
      lrhtri = ncomp*nmax
!
      idelu = irhtri + lrhtri
      ldelu = ncomp*nmax
!
      ixmer = idelu + ldelu
      lxmer = ncomp*nmax
!
!     rerr occupies the same space as xmerit
!
      irerr = ixmer
      lrerr = ncomp*nmax
!
      iutri = irerr + lrerr
      lutri = ncomp*nmax
!
      iermx = iutri + lutri
      lermx = nmax
!
      irtdc = iermx + lermx
      lrtdc = nmax
!
      ixxold = irtdc + lrtdc
      lxxold = nmax
!
      iuint = ixxold + lxxold
      luint = ncomp
!
      iftmp = iuint + luint
      lftmp = ncomp
!
      idgtm = iftmp + lftmp
      ldgtm = ncomp
!
      idftm1 = idgtm + ldgtm
      ldftm1 = ncomp*ncomp
!
      idftm2 = idftm1 + ldftm1
      ldftm2 = ncomp*ncomp
!
      itmp = idftm2 + ldftm2
      ltmp = ncomp*8
!
      idsq = itmp + ltmp
      ldsq = ncomp*ncomp
!
      idexr = idsq + ldsq
      ldexr = ncomp
!
      ietst6 = idexr + ldexr
      letst6 = ntol
!
      ietst8 = ietst6 + letst6
      letst8 = ntol
!
      iamg = ietst8 + letst8
      lamg = ncomp*nmax
!
      ic1 = iamg + lamg
      lc1 = ncomp*ncomp*nmax
!
      iwrkrhs = ic1+lc1
      lwrkrhs = ncomp*nmax
!
      ilast = iwrkrhs +  lwrkrhs
!
      IF (iprint == 1) WRITE(6,903) ilast
903   FORMAT(1H ,'ilast',i10)
!
!     Partition integer workspace.
!
      iiref = 1
      liref = nmax
!
      iihcom = iiref + liref
      lihcom = nmax
!
      iipvbk = iihcom + lihcom
      lipvbk = ncomp*nmax
!
      iipvlu = iipvbk + lipvbk
      lipvlu = ncomp
!
      iisign = iipvlu + lipvlu
      lisign = ncomp*nmax
!
      CALL BVPSOL(ncomp, nmsh, nlbc, aleft, aright,  &
                  nfxpnt, fixpnt, ntol, ltol, tol, nmax,   &
                  linear, giveu, givmsh, xx, nudim, u,  &
                  wrk(idefex), wrk(idefim), wrk(idef),  &
                  wrk(idelu), wrk(irhs), wrk(ifval),  &
                  wrk(itpblk), wrk(ibtblk), wrk(iajac), wrk(ibhold),  &
                  wrk(ichold), iwrk(iipvbk), iwrk(iipvlu), iwrk(iisign),  &
                  wrk(iuint), wrk(iftmp), wrk(itmrhs),  &
                  wrk(idftm1), wrk(idftm2), wrk(idgtm),  &
                  wrk(iutri), wrk(irhtri), wrk(ixmer),  &
                  wrk(ixxold), wrk(iuold), wrk(iusve),  &
                  wrk(itmp), wrk(idsq), wrk(idexr), wrk(irtdc),  &
                  wrk(irerr), wrk(ietst6), wrk(ietst8), wrk(iermx),  &
                  iwrk(iihcom), iwrk(iiref), wrk(idef6), wrk(idef8),  &
                  iflbvp,  &
                  wrk(iamg), wrk(ic1), wrk(iwrkrhs),  &
                  ckappa1, gamma1, ckappa, rpar, ipar)
!
      RETURN
   END SUBROUTINE TWPBVPC
!!!
!!!
   SUBROUTINE BVPSOL(ncomp, nmsh, nlbc, aleft, aright,  &
                     nfxpnt, fixpnt, ntol, ltol, tol, nmax,  &
                     linear, giveu, givmsh, xx, nudim, u,  &
                     defexp, defimp, def,  &
                     delu, rhs, fval,  &
                     topblk, botblk, ajac, bhold,  &
                     chold, ipvblk, ipivlu, isign,  &
                     uint, ftmp, tmprhs,  &
                     dftmp1, dftmp2, dgtm,  &
                     utrial, rhstri, xmerit,  &
                     xxold, uold, usave,  &
                     tmp, dsq, dexr, ratdc,  &
                     rerr, etest6, etest8, ermx,  &
                     ihcomp, irefin, def6, def8,  &
                     iflbvp,  &
                     amg, c1, wrkrhs,  &
                     ckappa1, gamma1, ckappa, rpar, ipar)
!
      USE BVPShared, ONLY: flmin, flmax, epsmch,  &
                           pdebug, use_c, comp_c, uval0, nminit, iprint, idum
!
      IMPLICIT NONE
!
      INTEGER(i4), INTENT(IN)    :: ncomp
      INTEGER(i4), INTENT(INOUT) :: nmsh
      INTEGER(i4), INTENT(IN)    :: nlbc
      REAL(r8),    INTENT(IN)    :: aleft
      REAL(r8),    INTENT(IN)    :: aright
      INTEGER(i4), INTENT(IN)    :: nfxpnt
      REAL(r8),    INTENT(IN)    :: fixpnt(*)
      INTEGER(i4), INTENT(IN)    :: ntol
      INTEGER(i4), INTENT(IN)    :: ltol(ntol)
      REAL(r8),    INTENT(IN)    :: tol(ntol)
      INTEGER(i4), INTENT(INOUT) :: nmax
      LOGICAL,     INTENT(IN)    :: linear
      LOGICAL,     INTENT(IN)    :: giveu
      LOGICAL,     INTENT(IN)    :: givmsh
      REAL(r8),    INTENT(INOUT) :: xx(*)
      INTEGER(i4), INTENT(IN)    :: nudim
      REAL(r8),    INTENT(INOUT) :: u(nudim,*)
      REAL(r8),    INTENT(INOUT) :: defexp(ncomp,*)
      REAL(r8),    INTENT(INOUT) :: defimp(ncomp,*)
      REAL(r8),    INTENT(INOUT) :: def(ncomp,*)
      REAL(r8),    INTENT(INOUT) :: delu(ncomp,*)
      REAL(r8),    INTENT(INOUT) :: rhs(*)
      REAL(r8),    INTENT(INOUT) :: fval(ncomp,*)
      REAL(r8),    INTENT(INOUT) :: topblk(nlbc,*)
      REAL(r8),    INTENT(INOUT) :: botblk(ncomp-nlbc,*)
      REAL(r8),    INTENT(INOUT) :: ajac(ncomp,2*ncomp,*)
      REAL(r8),    INTENT(INOUT) :: bhold(ncomp,ncomp,*)
      REAL(r8),    INTENT(INOUT) :: chold(ncomp,ncomp,*)
      INTEGER(i4), INTENT(INOUT) :: ipvblk(*)
      INTEGER(i4), INTENT(INOUT) :: ipivlu(*)
      INTEGER(i4), INTENT(INOUT) :: isign(*)
      REAL(r8),    INTENT(INOUT) :: uint(ncomp)
      REAL(r8),    INTENT(INOUT) :: ftmp(ncomp)
      REAL(r8),    INTENT(INOUT) :: tmprhs(*)
      REAL(r8),    INTENT(INOUT) :: dftmp1(ncomp,ncomp)
      REAL(r8),    INTENT(INOUT) :: dftmp2(ncomp,ncomp)
      REAL(r8),    INTENT(INOUT) :: dgtm(ncomp)
      REAL(r8),    INTENT(INOUT) :: utrial(ncomp,*)
      REAL(r8),    INTENT(INOUT) :: rhstri(*)
      REAL(r8),    INTENT(INOUT) :: xmerit(ncomp,*)
      REAL(r8),    INTENT(INOUT) :: xxold(*)
      REAL(r8),    INTENT(INOUT) :: uold(ncomp,*)
      REAL(r8),    INTENT(INOUT) :: usave(ncomp,*)
      REAL(r8),    INTENT(INOUT) :: tmp(ncomp,8)
      REAL(r8),    INTENT(INOUT) :: dsq(ncomp,ncomp)
      REAL(r8),    INTENT(INOUT) :: dexr(ncomp)
      REAL(r8),    INTENT(INOUT) :: ratdc(*)
      REAL(r8),    INTENT(INOUT) :: rerr(ncomp,*)
      REAL(r8),    INTENT(INOUT) :: etest6(*)
      REAL(r8),    INTENT(INOUT) :: etest8(*)
      REAL(r8),    INTENT(INOUT) :: ermx(*)
      INTEGER(i4), INTENT(INOUT) :: ihcomp(*)
      INTEGER(i4), INTENT(INOUT) :: irefin(*)
      REAL(r8),    INTENT(INOUT) :: def6(ncomp,*)
      REAL(r8),    INTENT(INOUT) :: def8(ncomp,*)
      INTEGER(i4), INTENT(INOUT) :: iflbvp
      REAL(r8),    INTENT(INOUT) :: amg(*)
      REAL(r8),    INTENT(INOUT) :: c1(ncomp,*)
      REAL(r8),    INTENT(INOUT) :: wrkrhs(*)
      REAL(r8),    INTENT(INOUT) :: ckappa1
      REAL(r8),    INTENT(INOUT) :: gamma1
      REAL(r8),    INTENT(INOUT) :: ckappa
      REAL(r8),    INTENT(INOUT) :: rpar(*)
      INTEGER(i4), INTENT(INOUT) :: ipar(*)
!
      REAL(r8),    PARAMETER :: zero = 0.0D+0
      REAL(r8),    PARAMETER :: one = 1.0D+0
      REAL(r8),    PARAMETER :: third = 0.33D+0
      REAL(r8),    PARAMETER :: fourth = 0.25D+0
      REAL(r8),    PARAMETER :: quan6 = 0.1D+0
      INTEGER(i4), PARAMETER :: itcondmax = 10
!
      REAL(r8)    :: gamma1old, ckappa1old, tolmin, rnsq, err6
      REAL(r8)    :: fxfct = 10.0D+0
      LOGICAL     :: ddouble
      LOGICAL     :: smooth, succes, strctr, trst6, reaft6
      LOGICAL     :: onto6, onto8, ludone, rhsgiv
      LOGICAL     :: first4, first8
      LOGICAL     :: frscal
      LOGICAL     :: stab_kappa, stab_gamma, stab_cond, stiff_cond
      LOGICAL, SAVE :: mchset = .TRUE.
      LOGICAL     :: maxmsh = .FALSE.
      INTEGER(i4) :: nmold, numbig, nummed, itcond, iorder, iflnwt, itnwt
      INTEGER(i4) :: ninter, n, iter
      INTEGER(i4) :: i
!!!
!
      frscal = .true.
      IF (mchset) THEN
         flmin = TINY(1._r8)  ! D1MACH(1)
         flmax = HUGE(1._r8)  ! D1MACH(2)
         epsmch = EPSILON(1._r8) / RADIX(1._r8)  !  D1MACH(3)
         IF (pdebug) WRITE(6,901) epsmch
         mchset = .false.
      END IF
!
!     The routine stcons calculates integration constants stored in
!     labeled common consts.
!
      CALL STCONS
!
!     Set up arrays for the error tests.
!
      IF (.NOT. linear) THEN
         CALL DLOAD(ntol, one, etest6, 1)
      ELSE
         DO i = 1, ntol
            etest6(i) = one/MAX(quan6, tol(i)**third)
         END DO
      END IF
!
      nmold  = 1
      smooth = .false.
      strctr = .false.
      trst6  = .true.
      reaft6 = .false.
      numbig = 0
      nummed = 0
      first4 = .true.
      first8 = .true.
      onto6  = .false.
!
!     Initialize parameter for the conditioning estimation
!
      itcond  = 0
      gamma1old  = flmax
      gamma1     = flmax
      ckappa1old  = flmax
      ckappa1     = flmax
      ckappa      = flmax
      stiff_cond = .false.
      stab_cond  = .false.
      tolmin = flmax
      DO i = 1, ntol
         tolmin = MIN(tol(i),tolmin)
      END DO
!
!     If givmsh is .true., the initial number of mesh points must be
!     provided by the user in nmsh, and the mesh points must be
!     contained in the array xx (of dimension nmsh).
!     Otherwise, nmsh is set to its default value, and a
!     uniform initial mesh is created.
!
      IF (.NOT. giveu .AND. .NOT. givmsh) THEN
         nmsh = nminit
         IF (nmsh < nfxpnt+2) nmsh = nfxpnt + 2
         CALL UNIMSH(nmsh, aleft, aright, nfxpnt, fixpnt, xx)
      END IF
      IF (pdebug) THEN
         WRITE(6,902)
         CALL SPRT(nmsh, xx)
      END IF
!
      IF (.NOT. giveu) CALL INITU(ncomp, nmsh, xx, nudim, u, rpar,ipar)
!
!     top of logic for 4th order solution ****
!
400   CONTINUE
      IF (iprint == 1) WRITE(6,903) nmsh
!
!     Set the def (deferred correction) array to zero.
! 
      CALL MTLOAD(ncomp, nmsh-1, zero, ncomp, def)
      iorder = 4
!
!     The routine fneval calls fsub at the mesh xx and the
!     solution u, and saves the values in the array fval.
!
      CALL FNEVAL(ncomp, nmsh, xx, nudim, u, fval, rpar, ipar)
!
!     Try to compute a 4th order solution by solving a system of nonlinear
!     equations.
!
      IF (linear) THEN
!
         ludone = .false.
         CALL LINEQ(ncomp, nmsh, nlbc, ludone, xx, nudim, u, def,  &
                    delu, rhs, fval, uint, ftmp, dftmp1, dftmp2, dgtm,  &
                    tmprhs, ajac, topblk, botblk, bhold, chold, ipvblk,  &
                    iflnwt, rpar, ipar)
!
!        Call fneval to evaluate the fval array at the new solution u.
!        (Such a call is not necessary for the nonlinear case because
!        fval is called within newteq for the new u.)
!
         CALL FNEVAL(ncomp, nmsh, xx, nudim, u, fval, rpar, ipar)
!
      ELSE
! 
         rhsgiv = .false.
         CALL NEWTEQ(ncomp, nmsh, nlbc, rhsgiv, ntol, ltol, tol,  &
                     xx, nudim, u, def, delu, rhs, fval,  &
                     utrial, rhstri, uint,  &
                     ftmp, dftmp1, dftmp2, dgtm, tmprhs, xmerit,  &
                     ajac, topblk, botblk, bhold, chold, ipvblk,  &
                     itnwt, iflnwt, rpar, ipar, frscal)
!
      END IF   ! if (linear) then
!
!     these flags are used in the mesh selection strategy
!
      IF (iflnwt == 0) THEN
!
!        COMPUTE ESTIMATIONS OF CONDITION NUMBERS KAPPA1 and GAMMA1
!        related to perturbation to the boundary conditions
!
         n = nmsh * ncomp
         ninter = nmsh - 1
!
         IF (comp_c) THEN
!
            gamma1old = gamma1
            ckappa1old = ckappa1
!
            CALL CONDESTIM(aleft, aright, nmsh, ncomp, n, xx, topblk,  &
                           nlbc, ncomp, ajac, ncomp, 2*ncomp, ninter, &
                           botblk, ncomp-nlbc, ipvblk, isign,  &
                           amg, c1 , wrkrhs, ckappa1, gamma1)
!
            IF (iprint >= 0) THEN
!     
!              COMPUTE ESTIMATION OF THE CONDITION NUMBER KAPPA
! 
               CALL ESTIMKAPPA(nmsh, ncomp, n, xx, topblk, nlbc, ncomp, ajac,  &
                               ncomp, 2*ncomp, ninter, botblk, ncomp-nlbc, ipvblk,  &
                               isign, c1, wrkrhs, ckappa)
!
            END IF
!    
            IF (iprint >= 0) THEN
               WRITE(6,1001) ckappa1/gamma1
               WRITE(6,1002) gamma1
               WRITE(6,1003) ckappa1
               WRITE(6,1004) ckappa
            END IF
!
            stab_kappa = ABS(ckappa1old-ckappa1)/(1D0+ckappa1) < 5D-2  &
                         .AND. ckappa1 < flmax .AND. gamma1 < flmax
! 
            stab_gamma = ABS(gamma1old-gamma1)/(1D0+gamma1) < 5D-2  &
                          .AND. gamma1 < flmax .AND. ckappa1 < flmax
! 
            stab_cond = stab_kappa .AND. stab_gamma
!
            stiff_cond = (( (ckappa1/gamma1 >= 10D0  )))
! 
            IF (iprint == 1) THEN
               WRITE(6,*) 'stab_kappa = ', stab_kappa
               WRITE(6,*) 'stab_gamma = ', stab_gamma
               WRITE(6,*) 'stiff_cond = ', stiff_cond
            END IF
! 
         END IF
!        endif if (comp_c)
!
         CALL CONV4(ncomp, nmsh, ntol, ltol, tol, linear, nmax,  &
                    xx, nudim, u, defexp, defimp, def, fval,  &
                    tmp, bhold, chold, dsq, dexr, usave, ratdc,  &
                    rerr, ipivlu, nmold, xxold,  &
                    smooth, reaft6, onto6, strctr, trst6, ddouble,  &
                    maxmsh, succes, first4,  &
                    amg, stab_cond, ckappa, stiff_cond, wrkrhs,  &
                    nfxpnt, fixpnt, irefin, rpar, ipar)
!
         IF (pdebug .AND. .NOT. onto6) WRITE (6,904)
! 
      ELSE
!
         IF (comp_c) THEN
!
            IF (iflnwt /= -1) THEN
!
               gamma1old = gamma1
               ckappa1old = ckappa1
               n =nmsh*ncomp
               ninter=nmsh-1
!
               CALL CONDESTIM(aleft, aright, nmsh, ncomp, n, xx, topblk,  &
                              nlbc, ncomp, ajac, ncomp, 2*ncomp, ninter,  &
                              botblk, ncomp-nlbc, ipvblk, isign,  &
                              amg, c1, wrkrhs, ckappa1, gamma1)
!
               IF (iprint >= 0) THEN
                  CALL ESTIMKAPPA(nmsh, ncomp, n, xx, topblk, nlbc, ncomp, ajac,  &
                                  ncomp, 2*ncomp, ninter, botblk, ncomp-nlbc, ipvblk,  &
                                  isign, c1, wrkrhs, ckappa)
               END IF
!
            END IF
! 
            IF (iprint >= 0) THEN
               WRITE(6,1001) ckappa1/gamma1
               WRITE(6,1002) gamma1
               WRITE(6,1003) ckappa1
               WRITE(6,1004) ckappa
            END IF
! 
            stab_kappa = ABS(ckappa1old-ckappa1)/(1D0+ckappa1) < 5D-2  &
                         .AND. ckappa1 < flmax .AND. gamma1 < flmax
! 
            stab_gamma = ABS(gamma1old-gamma1)/(1D0+gamma1) < 5D-2  &
                         .AND. gamma1 < flmax .AND. ckappa1 < flmax
! 
            stab_cond = stab_kappa .AND. stab_gamma
!
            stiff_cond = (( (ckappa1/gamma1 >= 5D0  )))
!
            IF (iprint == 1) THEN
               WRITE(6,*) 'stab_kappa = ',stab_kappa
               WRITE(6,*) 'stab_gamma = ', stab_gamma
               WRITE(6,*) 'stiff_cond = ', stiff_cond
            END IF
!
         END IF
!        end if if(comp_c)
!
         succes = .false.
         onto6 = .false.
         reaft6 = .false.
!
         CALL FAIL4(ncomp, nmsh, nlbc, ntol, ltol,  &
                    xx, nudim, u, rhs, linear, nmax, nmold, xxold, uold, ratdc,  &
                    iorder, iflnwt, itnwt, ddouble , maxmsh,  &
                    numbig, nummed, wrkrhs, amg, stab_cond, stiff_cond,  &
                    nfxpnt, fixpnt, irefin, itcond, itcondmax, rpar, ipar)
!
      END IF  !   if (iflnwt == 0) then
!
      IF (succes) THEN
         IF (comp_c) THEN
            CALL ESTIMKAPPA(nmsh, ncomp, n, xx, topblk, nlbc, ncomp, ajac,  &
                            ncomp, 2*ncomp, ninter, botblk, ncomp-nlbc, ipvblk,  &
                            isign, c1, wrkrhs, ckappa)
         END IF
         IF (iprint /= -1 .AND. comp_c) THEN
            IF (ckappa >= tolmin/epsmch) WRITE(6,1005)
         END IF
         iflbvp = 0
         RETURN
      ELSE IF (maxmsh) THEN
         GO TO 900
      ELSE IF (.NOT. onto6)  THEN
         GO TO 400
      END IF
!
!     To reach here, onto6 must be .true.
!     logic for 6th order ****
!
      IF (iprint == 1) WRITE(6,905)
!
!     Save the 4th order solution on this mesh in uold.
!
      CALL MATCOP(nudim, ncomp, ncomp, nmsh, u, uold)
!
!     Copy the current mesh into the xxold array.
!
      nmold = nmsh
      CALL DCOPY(nmold, xx, 1, xxold, 1)
!
      iorder = 6
!
      IF (linear) THEN
!
         CALL LINEQ(ncomp, nmsh, nlbc, ludone, xx, nudim, u, def,  &
                    delu, rhs, fval, uint, ftmp, dftmp1, dftmp2, dgtm,  &
                    tmprhs, ajac, topblk, botblk, bhold, chold, ipvblk,  &
                    iflnwt,rpar,ipar)
!
      ELSE
!
         CALL FIXJAC(ncomp, nmsh, nlbc, iorder, ntol, ltol, tol, xx,  &
                     nudim, u, def, def, delu, rhs, fval, utrial, rhstri,  &
                     rnsq, uint, ftmp, tmprhs, ajac, topblk, botblk, ipvblk,  &
                     iflnwt,rpar,ipar)
!
!        If the fixed Jacobian iterations fail but rnsq is small,
!        try a Newton procedure.  Set rhsgiv to indicate that
!        the right-hand side and fval have already been evaluated
!        at the given u.
! 
         IF (iflnwt == -3 .AND. rnsq < fxfct*epsmch) THEN
!
            rhsgiv = .true.
!    
            CALL NEWTEQ(ncomp, nmsh, nlbc, rhsgiv, ntol, ltol, tol,  &
                        xx, nudim, u, def, delu, rhs, fval,  &
                        utrial, rhstri, uint,  &
                        ftmp, dftmp1, dftmp2, dgtm, tmprhs, xmerit,  &
                        ajac, topblk, botblk, bhold, chold, ipvblk,  &
                        iter, iflnwt, rpar, ipar, frscal)
! 
         END IF
!
      END IF
!
      IF (iflnwt == 0) THEN
!
         CALL CONV6(ncomp, nmsh, ntol, ltol, tol, nudim, u, uold,  &
                    etest6, err6, trst6, onto8, reaft6, succes)
!
      ELSE
! 
         onto8 = .false.
!  
         CALL FAIL6(ncomp, nmsh, nlbc, ntol, ltol, tol, nfxpnt, fixpnt,  &
                    iorder, nmax, xx, nudim, u, rhs, usave,  &
                    xxold, uold, nmold, ihcomp, irefin,  &
                    rerr, ermx, ratdc, reaft6, ddouble, succes, maxmsh,  &
                    numbig, nummed, wrkrhs, amg, stab_cond,  &
                    ckappa1, gamma1, ckappa, stiff_cond, itcond, itcondmax)
!
      END IF
!
      IF (succes) THEN
         IF (comp_c) THEN
            CALL ESTIMKAPPA(nmsh, ncomp, n, xx, topblk, nlbc, ncomp, ajac,  &
                            ncomp, 2*ncomp, ninter, botblk, ncomp-nlbc, ipvblk,  &
                            isign, c1, wrkrhs, ckappa)
         END IF
         IF (iprint /= -1 .AND. comp_c) THEN
            IF (ckappa >= tolmin/epsmch) WRITE(6,1005)
         END IF
         iflbvp = 0
         RETURN
      ELSE IF (maxmsh) THEN
         GO TO 900
      ELSE IF (.NOT. onto8) THEN
         GO TO 400
      END IF
!
!     Logic for trying to calculate 8th order solution
!
      IF (iprint == 1) WRITE(6,906)
!
      CALL MATCOP(nudim, ncomp, ncomp, nmsh, u, uold)
!     Copy the current mesh into the xxold array.
      nmold = nmsh
      CALL DCOPY(nmold, xx, 1, xxold, 1)
!
!     Save the old deferred correction vector def in def6.
!
      CALL MATCOP(ncomp, ncomp, ncomp, nmsh-1, def, def6)
!
!     For linear problems, calculate the fval array for the
!     new solution u.
!
      IF (linear) CALL FNEVAL(ncomp, nmsh, xx, nudim, u, fval, rpar,ipar)
!
!     Calculate 8th order deferred corrections (the array def8).
!
      CALL DF8CAL(ncomp, nmsh, xx, nudim, u, fval, def8, tmp, rpar,ipar)
!
!     For linear problems, the def array is the def8 array.
!     For nonlinear problems, add the def8 array to the
!     already-calculated def array.
!
      IF (linear) THEN
         CALL MATCOP(ncomp, ncomp, ncomp, nmsh-1, def8, def)
      ELSE
         CALL MAXPY (ncomp, nmsh-1, one, def8, ncomp, def)
      END IF
!
      iorder = 8
!
      IF (linear) THEN
!
         CALL LINEQ(ncomp, nmsh, nlbc, ludone, xx, nudim, u, def, delu,  &
                    rhs, fval, uint, ftmp, dftmp1, dftmp2, dgtm, tmprhs,  &
                    ajac, topblk, botblk, bhold, chold, ipvblk,  &
                    iflnwt,rpar,ipar)
!
      ELSE
! 
         CALL FIXJAC(ncomp, nmsh, nlbc, iorder, ntol, ltol, tol, xx,  &
                     nudim, u, def, def8, delu, rhs, fval, utrial, rhstri,  &
                     rnsq, uint, ftmp, tmprhs, ajac, topblk, botblk, ipvblk,  &
                     iflnwt, rpar, ipar)
! 
!        If the fixed Jacobian iterations fail but rnsq is small,
!        try a Newton procedure.  Set rhsgiv to indicate that
!        the right-hand side and fval have already been evaluated
!        at the given u.
!
         IF (iflnwt == -3 .AND. rnsq < fxfct*epsmch) THEN
!
            rhsgiv = .true.
! 
            CALL NEWTEQ(ncomp, nmsh, nlbc, rhsgiv, ntol, ltol, tol,  &
                        xx, nudim, u, def, delu, rhs, fval,  &
                        utrial, rhstri, uint,  &
                        ftmp, dftmp1, dftmp2, dgtm, tmprhs, xmerit,  &
                        ajac, topblk, botblk, bhold, chold, ipvblk,  &
                        iter, iflnwt, rpar, ipar, frscal)
!
         END IF
!
      END IF
!
      IF (iflnwt == 0) THEN
!
         CALL CONV8(ncomp, nmsh, ntol, ltol, tol, nfxpnt, fixpnt, linear,  &
                    nmax, xx, nudim, u, def, def6, def8, uold,  &
                    ihcomp, irefin, ermx, err6, etest8, strctr, ddouble,  &
                    nmold, xxold, maxmsh, succes, first8, wrkrhs,  &
                    amg, stab_cond, ckappa1, gamma1, ckappa, stiff_cond,  &
                    rpar, ipar)
! 
      ELSE
! 
         succes = .false.
         CALL FAIL8(ncomp, nmsh, nfxpnt, fixpnt, nmax, ntol, ltol, tol, nmold,  &
                    xx, nudim, u, def6, xxold, uold, ihcomp, irefin, ermx, ddouble , maxmsh,  &
                    wrkrhs,amg, stiff_cond, stab_cond)
! 
      END IF
!
      IF (maxmsh) THEN
         GO TO 900
      ELSE IF (.NOT. succes) THEN
         GO TO 400
      END IF
!
!     Successful termination.
!
      IF (comp_c) THEN
         CALL ESTIMKAPPA(nmsh, ncomp, n, xx, topblk, nlbc, ncomp, ajac,  &
                         ncomp, 2*ncomp, ninter, botblk, ncomp-nlbc, ipvblk,  &
                         isign,c1,wrkrhs,ckappa)
      END IF
      IF (iprint /= -1 .AND. comp_c ) THEN
         IF ( ckappa >= tolmin/epsmch) WRITE(6,1005)
      END IF
      iflbvp = 0
!
      RETURN
900   CONTINUE
!
!     Error exit---too many mesh points.
!
      iflbvp = 1
      WRITE(6,*) 'Terminated, too many mesh points'
      IF (iprint /= -1 .AND. comp_c ) THEN
         IF (linear .AND. ckappa >= tolmin/epsmch) WRITE(6,1006)
         IF (.NOT.linear .AND. ckappa >= tolmin/epsmch) WRITE(6,1007)
      END IF
!
901   FORMAT(1H ,'epsmch',1PE10.3)
902   FORMAT(1H ,'initial mesh')
903   FORMAT(1H ,'start 4th order, nmsh',i5)
904   FORMAT(1H ,'do not go on to 6th')
905   FORMAT(1H ,'start 6th order')
906   FORMAT(1H ,'start 8th order')
1001  FORMAT(1H ,'stiffness = ',1PE11.3)
1002  FORMAT(1H ,'gamma1    = ',1PE11.3)
1003  FORMAT(1H ,'kappa1    = ',1PE11.3)
1004  FORMAT(1H ,'kappa     = ',1PE11.3)
1005  FORMAT(1H ,'The problem is ill-conditioned, ',  &
      ' the solution could be inaccurate')
1006  FORMAT(1H ,'The problem is ill-conditioned,',  &
      ' try with a less stringent tolerance')
1007  FORMAT(1H ,'The problem is ill-conditioned,',  &
      ' try with a less stringent tolerance', ' or with a different initial guess' )
!
      RETURN
   END SUBROUTINE BVPSOL
!!!
!!!
   SUBROUTINE CONDESTIM(aleft, aright, nmsh, ncomp, n, xx, topblk,  &
                        nrwtop, novrlp, array, nrwblk, nclblk, nbloks, &
                        botblk, nrwbot, ipvcd, isign,  &
                        omg, c1, work, kppa, gamma)
!
!     COMPUTES THE FIRST AND LAST BLOCK COLUMN OF INVERSE MATRIX AND
!     ESTIMATE THE CONDITION NUMBERS KPPA, GAMMA
!
      IMPLICIT NONE
!
      REAL(r8),    INTENT(IN)    :: aleft
      REAL(r8),    INTENT(IN)    :: aright
      INTEGER(i4), INTENT(INOUT) :: nmsh
      INTEGER(i4), INTENT(IN)    :: ncomp
      INTEGER(i4), INTENT(IN)    :: n
      REAL(r8),    INTENT(INOUT) :: xx(*)
      REAL(r8),    INTENT(INOUT) :: topblk(nrwtop,*)
      INTEGER(i4), INTENT(IN)    :: nrwtop
      INTEGER(i4), INTENT(IN)    :: novrlp
      REAL(r8),    INTENT(INOUT) :: array(nrwblk,nclblk,*)
      INTEGER(i4), INTENT(IN)    :: nrwblk
      INTEGER(i4), INTENT(IN)    :: nclblk
      INTEGER(i4), INTENT(IN)    :: nbloks
      REAL(r8),    INTENT(INOUT) :: botblk(nrwbot,*)
      INTEGER(i4), INTENT(IN)    :: nrwbot
      INTEGER(i4), INTENT(INOUT) :: ipvcd(*)
      INTEGER(i4), INTENT(INOUT) :: isign(*)
      REAL(r8),    INTENT(INOUT) :: omg(*)
      REAL(r8),    INTENT(INOUT) :: c1(ncomp,*)
      REAL(r8),    INTENT(INOUT) :: work(*)
      REAL(r8),    INTENT(INOUT) :: kppa
      REAL(r8),    INTENT(INOUT) :: gamma
!
      REAL(r8) :: gamma1, minmg, csum
      INTEGER(i4) :: idmx, idmn, idomg, job
      INTEGER(i4) :: k, i, j, l
!!!
!
      minmg = 0.0D+0
      kppa = 0.0D+0
!
      DO k = 1, ncomp
!
         DO l = 1, n
            work(l) = 0.0D0
            IF (k <= nrwtop) THEN
               work(k) = 1.0D0
            ELSE
               work(n-ncomp+k) = 1.0D0
            END IF
         END DO
!  
         CALL CRSLVE(topblk, nrwtop, novrlp, array, nrwblk, nclblk, nbloks,  &
                     botblk, nrwbot, ipvcd, work, 0)
!
         DO l = 1, n
            c1(k,l)=work(l)
         END DO
!
      END DO
!
!     ESTIMATION OF CONDITION NUMBER BY BRUGNANO,TRIGIANTE & MAZZIA
!     infinity-norm
!
      l = 1
      DO j = 1, n-ncomp+1, ncomp
         omg(l) = 0D0
         DO i = j,j+ncomp-1
            csum = 0D0
            DO k = 1, ncomp
               csum = csum + DABS(c1(k,i))
            END DO
            omg(l) = DMAX1(omg(l),csum)
         END DO
         l = l + 1
      END DO
!
!     1-norm
!
      gamma = 0.0D0
!
      DO i = 2, nbloks+1
         IF (omg(i) > omg(i-1)) THEN
            gamma1=omg(i)* (xx(i)-xx(i-1))
         ELSE
            gamma1=omg(i-1)* (xx(i)-xx(i-1))
         END IF
         gamma = gamma + gamma1
      END DO
      idomg = IDAMAX(nbloks+1,omg,1)
      gamma = gamma/(aright-aleft)
!
      idmx = 1
      idmn = 1
      minmg = omg(1)
      kppa = omg(1)
!
      DO j = 2, nbloks+1
         IF (omg(j) > kppa) THEN
            kppa = omg(j)
            idmx = j
         END IF
         IF (omg(j) < minmg) THEN
            minmg = omg(j)
            idmn = j
         END IF
      END DO
!
      RETURN
   END SUBROUTINE CONDESTIM
!!!
!!!
   SUBROUTINE ESTIMKAPPA(nmsh, ncomp, n, xx, topblk, nrwtop, novrlp, array,  &
                         nrwblk, nclblk, nbloks, botblk, nrwbot, ipvcd,  &
                         isign, c1, work, ckappa)
!
!     ESTIMATE THE CONDITION NUMBER  kappa
!
      IMPLICIT NONE
!
      INTEGER(i4), INTENT(INOUT) :: nmsh
      INTEGER(i4), INTENT(IN)    :: ncomp
      INTEGER(i4), INTENT(IN)    :: n
      REAL(r8),    INTENT(INOUT) :: xx(*)
      REAL(r8),    INTENT(INOUT) :: topblk(nrwtop,*)
      INTEGER(i4), INTENT(IN)    :: nrwtop
      INTEGER(i4), INTENT(IN)    :: novrlp
      REAL(r8),    INTENT(INOUT) :: array(nrwblk,nclblk,*)
      INTEGER(i4), INTENT(IN)    :: nrwblk
      INTEGER(i4), INTENT(IN)    :: nclblk
      INTEGER(i4), INTENT(IN)    :: nbloks
      REAL(r8),    INTENT(INOUT) :: botblk(nrwbot,*)
      INTEGER(i4), INTENT(IN)    :: nrwbot
      INTEGER(i4), INTENT(INOUT) :: ipvcd(*)
      INTEGER(i4), INTENT(INOUT) :: isign(*)
      REAL(r8),    INTENT(INOUT) :: c1(ncomp,*)
      REAL(r8),    INTENT(INOUT) :: work(*)
      REAL(r8),    INTENT(INOUT) :: ckappa
!
      INTEGER(i4) :: job
      INTEGER(i4) :: kase, isolve
      INTEGER(i4) :: k, i, j, l, jj
!
!     DO THE CONDITION NUMBER ESTIMATION BY HIGHAM:
!     (DONE IN COLROW) infinity-norm
!
      isolve = 0
      kase = 0
!
55    CALL DONEST(n, c1, work, isign, ckappa, kase)
!
      IF (kase /= 0) THEN
!
         isolve = isolve+1
         IF (kase == 1) THEN
            job = 1
         ELSE
            job = 0
         END IF
!
         IF (job == 0) THEN
            DO i = 1, nbloks
               DO j = 1, nrwblk
                  work((i-1)*nrwblk+nrwtop+j) =  &
                  (xx(i+1)-xx(i))*work((i-1)*nrwblk+nrwtop+j)
               END DO
            END DO
         END IF
!
         CALL CRSLVE(topblk, nrwtop, novrlp, array, nrwblk, nclblk, nbloks,  &
                     botblk, nrwbot, ipvcd, work, job)
!
         IF (job == 1) THEN
            DO i = 1, nbloks
               DO j = 1, nrwblk
                  work((i-1)*nrwblk+nrwtop+j) =  &
                  (xx(i+1)-xx(i))*work((i-1)*nrwblk+nrwtop+j)
               END DO
            END DO
         END IF
!
         GO TO 55
      END IF
!
      RETURN
   END SUBROUTINE ESTIMKAPPA
!!!
!!!
   SUBROUTINE CONV4(ncomp, nmsh, ntol, ltol, tol, linear, nmax,  &
                    xx, nudim, u, defexp, defimp, def, fval,  &
                    tmp, bhold, chold, dsq, dexr, usave, ratdc,  &
                    rerr, ipivot, nmold, xxold,  &
                    smooth, reaft6, onto6, strctr, trst6, ddouble ,  &
                    maxmsh, succes, first4,  &
                    amg, stab_cond, ckappa, stiff_cond, r4,  &
                    nfxpnt, fixpnt, irefin, rpar, ipar)
!
      USE BVPShared, ONLY: pdebug, use_c, comp_c, uval0, nminit,  &
                           iprint, idum,  &
                           flmin, flmax, epsmch
      USE BVPExtern, ONLY: INITU, FSUB
!
      IMPLICIT NONE
!
      INTEGER(i4), INTENT(IN)    :: ncomp
      INTEGER(i4), INTENT(INOUT) :: nmsh
      INTEGER(i4), INTENT(IN)    :: ntol
      INTEGER(i4), INTENT(IN)    :: ltol(ntol)
      REAL(r8),    INTENT(IN)    :: tol(ntol)
      LOGICAL,     INTENT(IN)    :: linear
      INTEGER(i4), INTENT(INOUT) :: nmax
      REAL(r8),    INTENT(INOUT) :: xx(*)
      INTEGER(i4), INTENT(IN)    :: nudim
      REAL(r8),    INTENT(INOUT) :: u(nudim,*)
      REAL(r8),    INTENT(INOUT) :: defexp(ncomp,*)
      REAL(r8),    INTENT(INOUT) :: defimp(ncomp,*)
      REAL(r8),    INTENT(INOUT) :: def(ncomp,*)
      REAL(r8),    INTENT(INOUT) :: fval(ncomp,*)
      REAL(r8),    INTENT(INOUT) :: tmp(ncomp,4)
      REAL(r8),    INTENT(INOUT) :: bhold(ncomp,ncomp,*)
      REAL(r8),    INTENT(INOUT) :: chold(ncomp, ncomp,*)
      REAL(r8),    INTENT(INOUT) :: dsq(ncomp,ncomp)
      REAL(r8),    INTENT(INOUT) :: dexr(ncomp)
      REAL(r8),    INTENT(INOUT) :: usave(ncomp,*)
      REAL(r8),    INTENT(INOUT) :: ratdc(*)
      REAL(r8),    INTENT(INOUT) :: rerr(ncomp,*)
      INTEGER(i4), INTENT(INOUT) :: ipivot(*)
      INTEGER(i4), INTENT(INOUT) :: nmold
      REAL(r8),    INTENT(INOUT) :: xxold(*)
      LOGICAL,     INTENT(INOUT) :: smooth
      LOGICAL,     INTENT(INOUT) :: reaft6
      LOGICAL,     INTENT(INOUT) :: onto6
      LOGICAL,     INTENT(INOUT) :: strctr
      LOGICAL,     INTENT(INOUT) :: trst6
      LOGICAL,     INTENT(INOUT) :: ddouble
      LOGICAL,     INTENT(INOUT) :: maxmsh
      LOGICAL,     INTENT(INOUT) :: succes
      LOGICAL,     INTENT(INOUT) :: first4
      REAL(r8),    INTENT(INOUT) :: amg(*)
      LOGICAL,     INTENT(INOUT) :: stab_cond
      REAL(r8),    INTENT(INOUT) :: ckappa
      LOGICAL,     INTENT(INOUT) :: stiff_cond
      REAL(r8),    INTENT(INOUT) :: r4(*)
      INTEGER(i4), INTENT(IN)    :: nfxpnt
      REAL(r8),    INTENT(IN)    :: fixpnt(*)
      INTEGER(i4), INTENT(INOUT) :: irefin(*)
      REAL(r8),    INTENT(INOUT) :: rpar(*)
      INTEGER(i4), INTENT(INOUT) :: ipar(*)
!
      REAL(r8), PARAMETER :: hundth = 1.0D+5
      REAL(r8), PARAMETER :: rerfct = 1.5D+0
      REAL(r8), PARAMETER :: power = 1.0D+0/6.0D+0
      REAL(r8), PARAMETER :: one = 1.0D+0
      REAL(r8), PARAMETER :: zero = 0.0D+0
      REAL(r8), PARAMETER :: huge = 1.0D+30
!
      REAL(r8), SAVE :: dfold, oldrt1
      REAL(r8) :: dfctol, dfexmx, derivm, dfimmx, rat1, rat2, remax, drat
      INTEGER(i4) :: incmp, inmsh, intol, itlmx, numadd
      LOGICAL, SAVE :: savedu, reposs
      LOGICAL :: callrt, oscchk, adjrer
!
!     The Newton iteration converged for a 4th order solution.
!
      IF (iprint == 1) WRITE(6,901)
!
      IF (first4) THEN
         dfold = zero
         oldrt1 = huge
         savedu = .false.
         reposs = .false.
         first4 = .false.
      END IF
!
      succes = .false.
      maxmsh = .false.
!
!     Compute the explicit deferred correction array defexp.
!     The parameter fval is an ncomp by nmsh array, calculated
!     by a previous call to fneval for this mesh and solution.
!
      CALL DFEXCL(ncomp, nmsh, xx, nudim, u, defexp, fval, tmp, rpar, ipar)
!
      CALL MATCOP(ncomp, ncomp, ncomp, nmsh-1, defexp, def)
!
!     If the problem has been declared as smooth, or we might wish
!     to perform Richardson extrapolation after trying to calculate a
!     6th order solution, just go on to the 6th order logic (don't
!     bother to calculate implicit deferred corrections).
!
      IF ( smooth .OR. reaft6) THEN
         IF (smooth .AND. pdebug) WRITE(6,902)
         IF (reaft6 .AND. pdebug) WRITE(6,903)
         onto6 = .true.
         RETURN
      END IF
!
!     Compute the cheap implicit deferred correction array defimp.
!     The array chold must be unchanged from a call to jaccal
!     for this mesh and solution.
!     The temporary arrays dsq and dexr are calculated inside dfimcl.
!
      CALL DFIMCL(ncomp, nmsh, defexp, chold, dsq, dexr, ipivot, defimp)
!
!     Call dccal to calculate: dfexmx, the maximum-magnitude element
!     of defexp in components for which tolerances are specified;
!     incmp, inmsh, and intol, the indices of the component,
!     mesh interval, and tolerance of dfexmx; derivm, the
!     maximum-magnitude element of fval(incmp,*) for all mesh
!     points; dfimmx, the maximum-magnitude element of defimp in
!     component incmp; the ratios rat1 and rat2; and the array
!     ratdc (used in osc).
!
      dfctol = hundth*epsmch
      CALL DCCAL(ncomp, nmsh, ntol, ltol, defexp, defimp, dfctol, fval,  &
                 ratdc, dfexmx, incmp, inmsh, intol, derivm, dfimmx,  &
                 rat1, rat2)
!
!     decid4 sets logical flags that determine the next step of the
!     algorithm.
!
      CALL DECID4(linear, rat1, rat2, dfexmx, dfimmx,  &
                  derivm, dfold, tol(intol), oldrt1,  &
                  onto6, smooth, callrt, strctr, oscchk, ddouble , reposs)
!
      IF (pdebug) THEN
         IF (smooth) WRITE(6,904)
         IF (callrt) WRITE(6,905)
         IF (oscchk) WRITE(6,906)
         IF (strctr) WRITE(6,907)
         IF (reposs) WRITE(6,908)
         IF (savedu) WRITE(6,909)
         IF (onto6) WRITE(6,910)
         WRITE(6,911) rat1, oldrt1
      END IF
      oldrt1 = rat1
      dfold = dfexmx
!
      IF (callrt) THEN
!  
!        ratcor calculates a more expensive implicit deferred correction.
!        The matrix bhold must be unchanged from the last call to jaccal.
!        If callrt is true, onto6 is always true also.
!
         CALL RATCOR(ncomp, nmsh, xx, defimp, bhold, def)
! 
      ELSE IF (linear) THEN
!
         IF (oscchk) THEN
!
            CALL OSC(ncomp, nmsh, dfexmx, incmp,  &
                     defexp, ratdc, ddouble , inmsh, onto6, trst6, smooth)
!
         ELSE IF (reposs)  THEN
!    
!           If reposs (`Richardson extrapolation possible') is true
!           for two successive 4th order solutions, a special termination
!           test may be used based on Richardson extrapolation.
!
!           If reposs is true for the first time, savedu will
!           be false; savedu can be true only when reposs is true for a
!           second consecutive mesh (a doubled version of the first).
!
            IF (savedu) THEN
! 
!              The flag savedu is .true. when the immediately preceding
!              converged 4th order solution was saved in the array usave,
!              and the mesh was then doubled.
!              In this case, the routine rerrvl is called to compute a
!              Richardson extrapolation (RE) error estimate remax.
!              The rerr array does not need to be adjusted, so adjrer is false.
!
               adjrer = .false.
               CALL RERRVL(ncomp, nmsh, nudim, u, usave, ntol,  &
                           ltol, rerr, remax, itlmx, adjrer )
! 
               IF (remax < rerfct*tol(itlmx)) THEN
                  succes = .true.
                  RETURN
               END IF
!              end of logic for savedu = .true.
!
            END IF
!    
!           Here, reposs is .true., but either we hadn't saved the previous
!           solution, or else the RE error estimate is not small.
!           Save u in usave, and set savedu to .true.
!           Set ddouble  to .true. to indicate that the mesh should be doubled.
!           Set .onto6. to .false. to return to the top of the 4th order
!           iteration loop.
! 
            CALL MATCOP(nudim, ncomp, ncomp, nmsh, u, usave)
            ddouble= .true.
            savedu = .true.
            onto6 = .false.
!
!           end of logic for reposs = .true.
         END IF
!
!     end of logic for linear
      END IF
!
      IF (pdebug .AND. reposs .AND. .NOT. onto6) WRITE(6,912)
!
!     if (stiff_cond .and. stab_cond) onto6 = .true.
!
      IF (.NOT. onto6) THEN
!
!        NB: onto6 can be false only for linear problems
!
         IF (ddouble) THEN
!
            IF ((use_c)) THEN
               IF (stiff_cond .AND. .NOT. stab_cond) THEN
                  CALL SELCONDMSH(ncomp, nmsh,  &
                                  nfxpnt, fixpnt, nmax, xx, irefin,  &
                                  nmold, xxold, ddouble, maxmsh, r4, amg)
                  ddouble= .false.
               ELSE
                  CALL DBLMSH(nmsh, nmax, xx, nmold, xxold, maxmsh)
               END IF
            ELSE
               CALL DBLMSH(nmsh, nmax, xx, nmold, xxold, maxmsh)
            END IF
!
         ELSE
! 
!           Refine the mesh near interval inmsh; the number of points
!           to be added depends on the relative size of dfexmx compared
!           to u and tol.
!    
!           If the conditining parameter are not stabilized and the problem
!           is considered stiff  the points are added and removed using
!           the monitor function given by the conditioning parameters
!           otherwise the old technique and the conditioning parameters
!           are used together.
!           new functions: selcondmsh, smpselcondmsh
! 
            IF (use_c) THEN
               IF ((stiff_cond)) THEN
                  drat = dfexmx/ (MAX(one, ABS(u(incmp,inmsh)))*tol(intol))
                  IF (pdebug) WRITE(6,913) drat, u(incmp,inmsh), tol(intol)
                  numadd = drat**power
                  CALL SMPSELCONDMSH(ncomp, nmsh,  &
                                     nfxpnt, fixpnt, nmax, xx, irefin, inmsh, numadd,  &
                                     nmold, xxold, ddouble , maxmsh,r4,amg)
               ELSE
                  drat = dfexmx/ (MAX(one, ABS(u(incmp,inmsh)))*tol(intol))
                  IF (pdebug) WRITE(6,913) drat, u(incmp,inmsh), tol(intol)
                  numadd = drat**power
                  CALL SMPMSH(nmsh, nmax, xx, inmsh, numadd, nmold, xxold, maxmsh)
               END IF
            ELSE
               drat = dfexmx/ (MAX(one, ABS(u(incmp,inmsh)))*tol(intol))
               IF (pdebug) WRITE(6,913) drat, u(incmp,inmsh), tol(intol)
               numadd = drat**power
               CALL SMPMSH(nmsh, nmax, xx, inmsh, numadd, nmold, xxold, maxmsh)
            END IF
!
         END IF
!
         IF (.NOT. maxmsh) CALL INITU(ncomp, nmsh, xx, nudim, u,rpar,ipar)
! 
!     end of logic for .not. onto6
      END IF
!
901   FORMAT(1H ,'conv4')
902   FORMAT(1H ,'smooth')
903   FORMAT(1H ,'reaft6')
904   FORMAT(1H ,'smooth')
905   FORMAT(1H ,'callrt')
906   FORMAT(1H ,'oscchk')
907   FORMAT(1H ,'strctr')
908   FORMAT(1H ,'reposs')
909   FORMAT(1H ,'savedu')
910   FORMAT(1H ,'onto6 after decid4')
911   FORMAT(1H ,'rat1,oldrt1',2(1PE11.3))
912   FORMAT(1H ,'reposs and not onto6')
913   FORMAT(1H ,'drat,u,tol',3(1PE11.3))
!
      RETURN
   END SUBROUTINE CONV4
!!!
!!!
   SUBROUTINE FAIL4(ncomp, nmsh, nlbc, ntol, ltol,  &
                    xx, nudim, u, rhs, linear, nmax, nmold, xxold, uold, tmwork,  &
                    iorder, iflnwt, itnwt, ddouble , maxmsh,  &
                    numbig, nummed, r4, amg, stab_cond, stiff_cond,  &
                    nfxpnt, fixpnt, irefin, itcond, itcondmax, rpar, ipar)
!
      USE BVPShared, ONLY: pdebug, use_c, comp_c, uval0, nminit, iprint, idum
!
      IMPLICIT NONE
!
      INTEGER(i4), INTENT(IN)    :: ncomp
      INTEGER(i4), INTENT(INOUT) :: nmsh
      INTEGER(i4), INTENT(IN)    :: nlbc
      INTEGER(i4), INTENT(IN)    :: ntol
      INTEGER(i4), INTENT(IN)    :: ltol(ntol)
      REAL(r8),    INTENT(INOUT) :: xx(*)
      INTEGER(i4), INTENT(IN)    :: nudim
      REAL(r8),    INTENT(INOUT) :: u(nudim,*)
      REAL(r8),    INTENT(INOUT) :: rhs(*)
      LOGICAL,     INTENT(IN)    :: linear
      INTEGER(i4), INTENT(INOUT) :: nmax
      INTEGER(i4), INTENT(INOUT) :: nmold
      REAL(r8),    INTENT(INOUT) :: xxold(*)
      REAL(r8),    INTENT(INOUT) :: uold(ncomp,*)
      REAL(r8),    INTENT(INOUT) :: tmwork(*)
      INTEGER(i4), INTENT(INOUT) :: iorder
      INTEGER(i4), INTENT(INOUT) :: iflnwt
      INTEGER(i4), INTENT(INOUT) :: itnwt
      LOGICAL,     INTENT(INOUT) :: ddouble
      LOGICAL,     INTENT(INOUT) :: maxmsh
      INTEGER(i4), INTENT(INOUT) :: numbig
      INTEGER(i4), INTENT(INOUT) :: nummed
      REAL(r8),    INTENT(INOUT) :: r4(*)
      REAL(r8),    INTENT(INOUT) :: amg(*)
      LOGICAL,     INTENT(INOUT) :: stab_cond
      LOGICAL,     INTENT(INOUT) :: stiff_cond
      INTEGER(i4), INTENT(IN)    :: nfxpnt
      REAL(r8),    INTENT(IN)    :: fixpnt(*)
      INTEGER(i4), INTENT(INOUT) :: irefin(*)
      INTEGER(i4), INTENT(INOUT) :: itcond
      INTEGER(i4), INTENT(IN)    :: itcondmax
      REAL(r8),    INTENT(INOUT) :: rpar(*)
      INTEGER(i4), INTENT(INOUT) :: ipar(*)
!!!
!
!     The Newton procedure failed to obtain a 4th order solution.
!
      IF (iprint == 1) WRITE(6,901)
!
      maxmsh = .false.
!
      IF (iflnwt == -1) THEN
! 
!        iflnwt = -1 means that the Jacobian was considered singular.
!        (This is the only possible failure for a linear problem.)
!
         CALL DBLMSH(nmsh, nmax, xx, nmold, xxold, maxmsh)
         CALL INITU(ncomp, nmsh, xx, nudim, u,rpar,ipar)
!
      ELSE
!
!        The routine mshref decides how to refine the mesh and then
!        performs the refinement, either by doubling or based on
!        the rhs vector at the best point.
!
         CALL MSHREF(ncomp, nmsh, nlbc, ntol, ltol, iorder, rhs, tmwork,  &
                     nmax, xx, nmold, xxold, ddouble, maxmsh,  &
                     numbig, nummed, amg, stab_cond, stiff_cond,  &
                     r4, nfxpnt,fixpnt, irefin, itcond, itcondmax)
! 
         IF (.NOT. maxmsh) THEN
            IF (linear  .OR. itnwt == 0) THEN
! .or.
!     *                            (itnwt.le.2 .and. iflnwt.ne.0)) then
               CALL INITU(ncomp, nmsh, xx, nudim, u,rpar,ipar)
            ELSE
!              Interpolate the partially converged solution.
               CALL MATCOP(nudim, ncomp, ncomp, nmold, u, uold)
               CALL INTERP(ncomp, nmsh, xx, nudim, u, nmold, xxold, uold)
            END IF
         END IF
! 
!        End of logic for failure because of some reason other than
!        singularity.
      END IF
!
901   FORMAT(1H ,'fail4')
!
      RETURN
   END SUBROUTINE FAIL4
!!!
!!!
   SUBROUTINE CONV6(ncomp, nmsh, ntol, ltol, tol,  &
                    nudim, u, uold, etest6, err6, trst6, onto8, reaft6, succes)
!
      USE BVPShared, ONLY: pdebug, use_c, comp_c,  &
                           uval0, nminit, iprint, idum
!
      IMPLICIT NONE
!
      INTEGER(i4), INTENT(IN)    :: ncomp
      INTEGER(i4), INTENT(INOUT) :: nmsh
      INTEGER(i4), INTENT(IN)    :: ntol
      INTEGER(i4), INTENT(IN)    :: ltol(ntol)
      REAL(r8),    INTENT(IN)    :: tol(ntol)
      INTEGER(i4), INTENT(IN)    :: nudim
      REAL(r8),    INTENT(INOUT) :: u(nudim,*)
      REAL(r8),    INTENT(INOUT) :: uold(ncomp,*)
      REAL(r8),    INTENT(INOUT) :: etest6(*)
      REAL(r8),    INTENT(INOUT) :: err6
      LOGICAL,     INTENT(INOUT) :: trst6
      LOGICAL,     INTENT(INOUT) :: onto8
      LOGICAL,     INTENT(INOUT) :: reaft6
      LOGICAL,     INTENT(INOUT) :: succes
!
      LOGICAL :: errok
!!!
!
!     The Newton iteration converged for a 6th-order solution.
!
      IF (iprint == 1) WRITE(6,901)
!
      succes = .false.
!
!     The logical flag reaft6 is true only when the 6th order solution
!     failed.  Since the 6th order solution converged, reaft6 is false.
!
      reaft6 = .false.
      onto8 = .true.
!
!     Calculate the error estimates for the 4th order solution.
!
      CALL ERREST(ncomp, nmsh, ntol, ltol, tol,  &
                  nudim, u, uold, etest6, err6, errok)
!
      IF (trst6 .AND. errok) THEN
         succes = .true.
         RETURN
      END IF
!
901 FORMAT(1H ,'conv6')
!
      RETURN
   END SUBROUTINE CONV6
!!!
!!!
   SUBROUTINE FAIL6(ncomp, nmsh, nlbc, ntol, ltol, tol,  &
                    nfxpnt, fixpnt, iorder, nmax, xx, nudim, u, rhs, usave, xxold, uold, nmold,  &
                    ihcomp, irefin, rerr, ermx, tmwork, reaft6, ddouble, succes, maxmsh,  &
                    numbig, nummed, r4, amg, stab_cond, ckappa1, gamma1, ckappa, stiff_cond,  &
                    itcond, itcondmax)
!
      USE BVPShared, ONLY: pdebug, use_c, comp_c, uval0,  &
                           nminit, iprint, idum
!
      IMPLICIT NONE
! 
      INTEGER(i4), INTENT(IN)    :: ncomp
      INTEGER(i4), INTENT(INOUT) :: nmsh
      INTEGER(i4), INTENT(IN)    :: nlbc
      INTEGER(i4), INTENT(IN)    :: ntol
      INTEGER(i4), INTENT(IN)    :: ltol(ntol)
      REAL(r8),    INTENT(IN)    :: tol(*)
      INTEGER(i4), INTENT(IN)    :: nfxpnt
      REAL(r8),    INTENT(IN)    :: fixpnt(*)
      INTEGER(i4), INTENT(INOUT) :: iorder
      INTEGER(i4), INTENT(INOUT) :: nmax
      REAL(r8),    INTENT(INOUT) :: xx(*)
      INTEGER(i4), INTENT(IN)    :: nudim
      REAL(r8),    INTENT(INOUT) :: u(nudim,*)
      REAL(r8),    INTENT(INOUT) :: rhs(*)
      REAL(r8),    INTENT(INOUT) :: usave(ncomp,*)
      REAL(r8),    INTENT(INOUT) :: xxold(*)
      REAL(r8),    INTENT(INOUT) :: uold(ncomp,*)
      INTEGER(i4), INTENT(INOUT) :: nmold
      INTEGER(i4), INTENT(INOUT) :: ihcomp(*)
      INTEGER(i4), INTENT(INOUT) :: irefin(*)
      REAL(r8),    INTENT(INOUT) :: rerr(ncomp,*)
      REAL(r8),    INTENT(INOUT) :: ermx(*)
      REAL(r8),    INTENT(INOUT) :: tmwork(*)
      LOGICAL,     INTENT(INOUT) :: reaft6
      LOGICAL,     INTENT(INOUT) :: ddouble
      LOGICAL,     INTENT(INOUT) :: succes
      LOGICAL,     INTENT(INOUT) :: maxmsh
      INTEGER(i4), INTENT(INOUT) :: numbig
      INTEGER(i4), INTENT(INOUT) :: nummed
      REAL(r8),    INTENT(INOUT) :: r4(*)
      REAL(r8),    INTENT(INOUT) :: amg(*)
      LOGICAL,     INTENT(INOUT) :: stab_cond
      REAL(r8),    INTENT(INOUT) :: ckappa1
      REAL(r8),    INTENT(INOUT) :: gamma1
      REAL(r8),    INTENT(INOUT) :: ckappa
      LOGICAL,     INTENT(INOUT) :: stiff_cond
      INTEGER(i4), INTENT(INOUT) :: itcond
      INTEGER(i4), INTENT(IN)    :: itcondmax
!
!  nmpt is the standard number of mesh points to be added to selected
!  intervals when the mesh is altered based on the distribution
!  in size of the rhs vector.
!
      INTEGER(i4), PARAMETER :: nmpt = 15
      REAL(r8),    PARAMETER :: eight = 8.0D+0
!
      REAL(r8) :: remax
      INTEGER(i4) :: itlmx, ipow
      LOGICAL :: adjrer
!!!
!
!     Non-convergence of 6th order.

      IF (iprint == 1) WRITE(6,901)
      succes = .false.
      maxmsh = .false.
!
!     NB: the problem must be nonlinear.  Linear problems will either
!     fail to converge for 4th order, or else, once they've converged
!     for 4th order, must converge for 6th order.
!
!     Restore the u array to the previous 4th order solution.
!
      CALL MATCOP(ncomp, nudim, ncomp, nmold, uold, u)
!
      IF (reaft6 .AND. iprint >= 0) WRITE(6,9999)
9999  FORMAT(1H ,'in fail6, reaft6is true')
      IF (.NOT.reaft6 .AND. iprint >= 0) WRITE(6,9998)
9998  FORMAT(1H ,'in fail6, not reaft6')
      IF (ddouble.AND. iprint >= 0) WRITE(6,9997)
9997  FORMAT(1H ,'in fail6, ddouble  is true')
      IF (.NOT.ddouble.AND. iprint >= 0 ) WRITE(6,9996)
9996  FORMAT(1H ,'in fail6, not double')
!
      IF (.NOT. reaft6 .OR. .NOT. ddouble ) THEN
! 
!        Here, either
!        (1) the mesh for which this 6th order solution failed
!        is not a doubled version of the immediately preceding mesh, or
!        (2) for the immediately preceding mesh, it is not true
!        that the 4th order converged and the 6th order failed.
! 
!        Setting reaft6 to .true. signals that Richardson extrapolation
!        may be possible if the next 6th order solution fails.  When
!        reaft6 is true, the routine conv4 immediately sets onto6 to true.
! 
         reaft6 = .true.
! 
!        Save the current 4th order solution in usave.
!
         CALL MATCOP(nudim, ncomp, ncomp, nmsh, u, usave)
! 
!        Use the distribution of sizes in the rhs vector to refine the mesh.
! 
         CALL MSHREF(ncomp, nmsh, nlbc, ntol, ltol, iorder, rhs, tmwork,  &
                     nmax, xx, nmold, xxold, ddouble, maxmsh,  &
                     numbig, nummed, amg, stab_cond, stiff_cond,  &
                     r4,nfxpnt, fixpnt, irefin, itcond, itcondmax)
!
         IF (.NOT. maxmsh) CALL INTERP(ncomp, nmsh, xx, nudim, u, nmold, xxold, uold)
!
      ELSE
!
!        Here, reaft6 and ddoubleare both true.  So for two consecutive
!        meshes, the 4th order converged and the 6th order failed,
!        and the second mesh is the ddoubleof the first.
!  
!        Calculate an error estimate from Richardson extrapolation
!        with the current and previous 4th order solutions.
!        (usave is the 4th order solution saved from the previous (halved)
!        mesh.)
!        Set addrer to .true. to signal that the rerr array should
!        be adjusted.
! 
         adjrer = .true.
! 
         CALL RERRVL(ncomp, nmsh, nudim, u, usave, ntol, ltol,  &
                     rerr, remax, itlmx, adjrer)
!
         WRITE(6,9994)
9994     FORMAT(1H ,'***in fail6')
         WRITE(6,9993) remax, eight*tol(itlmx)
9993     FORMAT(1H ,'remax',1PE14.4,5X,'8*tol',1PE14.4)
!
         IF (remax < eight*tol(itlmx)) THEN
!
            succes = .true.
!
         ELSE
!    
!           Richardson extrapolation did not give sufficient accuracy.
!           Perform selective mesh refinement on the OLD (NB: old!) mesh
!           and the old (saved) solution, using the error estimate from
!           Richardson extrapolation to guide where the mesh points are placed.
!
!f          controllare se posso mettere 2
            nmsh = 1 + (nmsh-1)/2
            CALL DCOPY(nmsh, xxold, 2, xx, 1)
            CALL DCOPY(nmsh, amg, 2, tmwork, 1)
            ipow = 4
! 
!           The rerr array is overwritten by selmsh.
            IF (use_c) THEN
               IF (.NOT. stiff_cond ) THEN
                  CALL SELMSH(ncomp, nmsh, ntol, ltol, tol,  &
                              nfxpnt, fixpnt, ipow, nmax, xx, ncomp, usave, rerr, irefin, ihcomp,  &
                              nmold, xxold, ermx, ddouble , maxmsh)
               ELSE
!                 The rerr array is overwritten by selconderrmsh.
                  CALL SELCONDERRMSH(ncomp, nmsh, ntol, ltol, tol,  &
                                     nfxpnt, fixpnt, ipow, nmax, xx, ncomp, usave, rerr, irefin, ihcomp,  &
                                     nmold, xxold, ermx, ddouble, maxmsh, r4, tmwork, stab_cond)
               END IF
            ELSE
               CALL SELMSH(ncomp, nmsh, ntol, ltol, tol, nfxpnt, fixpnt, ipow, nmax,  &
                           xx, ncomp, usave, rerr, irefin, ihcomp,  &
                           nmold, xxold, ermx, ddouble , maxmsh)
            END IF
! 
!           If ddoubleis false on exit from selmsh, the call to selmsh has
!           produced a different (non-doubled) mesh.   Interpolate the
!           saved solution (from the old mesh) onto the mesh newly created
!           by selmsh.
!           NB: Because ddouble is false, we won't try Richardson extrapolation
!           if the next 4th order converges and 6th order fails.
! 
            IF (.NOT. maxmsh) THEN
               IF (.NOT. ddouble ) THEN
                  CALL INTERP(ncomp, nmsh, xx, nudim, u, nmold, xxold, usave)
               ELSE
! 
!                 Selective mesh refinement based on the old mesh simply
!                 produced the same mesh that we started with.  So now refine
!                 starting with the doubled mesh (from before) and the solution.
! 
                  reaft6 = .true.
!
!                 Save the solution in usave in case we can carry out Richardson
!                 extrapolation in the same circumstances next time.
!
                  CALL MATCOP(nudim, ncomp, ncomp, nmsh, u, usave)
!
!                 Use the distribution of sizes in the rhs vector to refine the mesh.
! 
                  CALL MSHREF(ncomp, nmsh, nlbc, ntol, ltol, iorder, rhs, tmwork,  &
                              nmax, xx, nmold, xxold, ddouble, maxmsh,  &
                              numbig, nummed, amg, stab_cond, stiff_cond,  &
                              r4,nfxpnt, fixpnt, irefin, itcond, itcondmax)
! 
                  IF (.NOT. maxmsh) CALL INTERP(ncomp, nmsh, xx, nudim, u, nmold, xxold, usave)
! 
!                 end of logic for needing to refine (again) based on the
!                 current mesh
               END IF
! 
!           end of logic for not too many mesh points
            END IF
! 
!        end of logic for failure of Richardson extrapolation
!        to produce a converged solution
         END IF
!  
!     end of logic for both reaft6 and ddouble being true
      END IF
!
901   FORMAT(1H ,'fail6')
!
      RETURN
   END SUBROUTINE FAIL6
!!!
!!!
   SUBROUTINE CONV8(ncomp, nmsh, ntol, ltol, tol,  &
                    nfxpnt, fixpnt, linear, nmax, xx, nudim, u, def, def6, def8, uold,  &
                    ihcomp, irefin, ermx, err6, etest8, strctr,  &
                    ddouble, nmold, xxold, maxmsh, succes, first8,  &
                    r4, amg, stab_cond, ckappa1, gamma1, ckappa, stiff_cond, rpar, ipar)
!
      USE BVPShared, ONLY: pdebug, use_c, comp_c, uval0,  &
                           nminit, iprint, idum
!
      IMPLICIT NONE
!
      INTEGER(i4), INTENT(IN)    :: ncomp
      INTEGER(i4), INTENT(INOUT) :: nmsh
      INTEGER(i4), INTENT(IN)    :: ntol
      INTEGER(i4), INTENT(IN)    :: ltol(ntol)
      REAL(r8),    INTENT(IN)    :: tol(ntol)
      INTEGER(i4), INTENT(IN)    :: nfxpnt
      REAL(r8),    INTENT(IN)    :: fixpnt(*)
      LOGICAL,     INTENT(IN)    :: linear
      INTEGER(i4), INTENT(INOUT) :: nmax
      REAL(r8),    INTENT(INOUT) :: xx(*)
      INTEGER(i4), INTENT(IN)    :: nudim
      REAL(r8),    INTENT(INOUT) :: u(nudim,*)
      REAL(r8),    INTENT(INOUT) :: def(ncomp,*)
      REAL(r8),    INTENT(INOUT) :: def6(ncomp,*)
      REAL(r8),    INTENT(INOUT) :: def8(ncomp,*)
      REAL(r8),    INTENT(INOUT) :: uold(ncomp,*)
      INTEGER(i4), INTENT(INOUT) :: ihcomp(*)
      INTEGER(i4), INTENT(INOUT) :: irefin(*)
      REAL(r8),    INTENT(INOUT) :: ermx(*)
      REAL(r8),    INTENT(INOUT) :: err6
      REAL(r8),    INTENT(INOUT) :: etest8(ntol)
      LOGICAL,     INTENT(INOUT) :: strctr
      LOGICAL,     INTENT(INOUT) :: ddouble
      INTEGER(i4), INTENT(INOUT) :: nmold
      REAL(r8),    INTENT(INOUT) :: xxold(*)
      LOGICAL,     INTENT(INOUT) :: maxmsh
      LOGICAL,     INTENT(INOUT) :: succes
      LOGICAL,     INTENT(INOUT) :: first8
      REAL(r8),    INTENT(INOUT) :: r4(*)
      REAL(r8),    INTENT(INOUT) :: amg(*)
      LOGICAL,     INTENT(INOUT) :: stab_cond
      REAL(r8),    INTENT(INOUT) :: ckappa1
      REAL(r8),    INTENT(INOUT) :: gamma1
      REAL(r8),    INTENT(INOUT) :: ckappa
      LOGICAL,     INTENT(INOUT) :: stiff_cond
      REAL(r8),    INTENT(INOUT) :: rpar(*)
      INTEGER(i4), INTENT(INOUT) :: ipar(*)
!
      REAL(r8), PARAMETER :: one = 1.0D+0
      REAL(r8), PARAMETER :: fourth = 0.25D+0
      REAL(r8), PARAMETER :: quan8 = 0.025D+0
      REAL(r8), PARAMETER :: efact  = 100.0D+0
      REAL(r8), PARAMETER :: huge = 1.0D+30
!
      REAL(r8) :: err8
      REAL(r8), SAVE :: er6old, er8old
      LOGICAL :: errok
      INTEGER(i4) :: i, ipow
!!!
!
!     The Newton iteration converged for the 8th order solution.
!
      IF (iprint == 1) WRITE(6,901)
!
      IF (first8) THEN
         er6old = huge
         er8old = huge
         first8 = .false.
      END IF
!
      IF (.NOT. linear) THEN
         CALL DLOAD(ntol, one, etest8, 1)
      ELSE
         DO i = 1, ntol
            etest8(i) = one/MAX(quan8, tol(i)**fourth)
         END DO
      END IF
      succes = .false.
      maxmsh = .false.
!
!     Check estimated error.  For a nonlinear problem, all components
!     of etest8 (the ratios used in testing the error) are set to one.
!     For a linear problem, the components of etest8 are in general
!     larger than one.  But if strctr is .true. and the number of mesh
!     points decreased, we set the elements of etest8 to one (which
!     makes a stricter test).
!
      IF (linear .AND. strctr .AND. nmsh < nmold) CALL DLOAD(ntol, one, etest8, 1)
!
      CALL ERREST(ncomp, nmsh, ntol, ltol, tol,  &
                  nudim, u, uold, etest8, err8, errok)
      IF (errok)  THEN
         succes = .true.
         RETURN
      END IF
!
!     At this point, the 8th order solution converged, but did not
!     satisfy the test for termination.
!
      IF (pdebug) WRITE(6,902) err6, err8, er6old, er8old
      IF (nmsh < nmold .AND. err6 > efact*er6old .AND.  &
          err8 > efact*er8old) THEN
!
!        If the number of mesh points decreased and the errors in the
!        6th and 8th order solutions did not decrease sufficiently compared
!        to the previous mesh, double the mesh and go back to try to
!        calculate a 4th order solution.
!
         CALL DBLMSH(nmsh, nmax, xx, nmold, xxold, maxmsh)
         IF (.NOT. maxmsh) THEN
            er6old = err6
            er8old = err8
!           If the problem is not linear we use the old solution
!           instead the the first value
! old code
!              call INITU(ncomp, nmsh, xx, nudim, u,rpar,ipar)
! new code
            IF (linear) THEN
               CALL INITU(ncomp, nmsh, xx, nudim, u,rpar,ipar)
            ELSE
!              we do not use u but uold
               CALL MATCOP(nudim, ncomp, ncomp, nmold, u, uold)
               CALL INTERP(ncomp, nmsh, xx, nudim, u, nmold, xxold, uold)
            END IF
         END IF
         RETURN
      END IF
!
!     Here, we know that
!     (1) the number of mesh points exceeds that for the previous mesh; or
!     (2) the number of mesh points decreased, the 6th and 8th order
!         errors did not satisfy the termination tests, but they did not
!         increase so much that the mesh needed to be doubled.
!
      er6old = err6
      er8old = err8
!
      IF (err8 <= err6) THEN
!  
!        Perform selective mesh refinement based on the 8th order deferred
!        corrections.  The value of ipow indicates that the error estimate
!        is of order 6.  Then, for a nonlinear problem, interpolate the
!        latest solution onto the new mesh.
! 
         ipow = 6
!
!        NB: The array def8 will be overwritten by selmsh.
!
         IF (use_c) THEN
            IF (.NOT. stiff_cond ) THEN
               CALL SELMSH(ncomp, nmsh, ntol, ltol, tol, nfxpnt, fixpnt, ipow, nmax,  &
                           xx, ncomp, uold, def8, irefin, ihcomp,  &
                           nmold, xxold, ermx, ddouble , maxmsh)
            ELSE
!              NB: The array def8 will be overwritten by selconderrmsh.
               CALL SELCONDERRMSH(ncomp, nmsh, ntol, ltol, tol,  &
                                  nfxpnt, fixpnt, ipow, nmax, xx, ncomp, uold, def8, irefin, ihcomp,  &
                                  nmold, xxold, ermx, ddouble , maxmsh,r4,amg,stab_cond)
            END IF
         ELSE
            CALL SELMSH(ncomp, nmsh, ntol, ltol, tol, nfxpnt, fixpnt, ipow, nmax,  &
                        xx, ncomp, uold, def8, irefin, ihcomp,  &
                        nmold, xxold, ermx, ddouble , maxmsh)
         END IF
! 
         IF (.NOT. maxmsh) THEN
            IF (linear) THEN
               CALL INITU(ncomp, nmsh, xx, nudim, u,rpar,ipar)
            ELSE
               CALL MATCOP(nudim, ncomp, ncomp, nmold, u, uold)
               CALL INTERP(ncomp, nmsh, xx, nudim, u, nmold, xxold, uold)
            END IF
         END IF
!
      ELSE
! 
!        err8 is greater than err6
! 
!        For a linear problem, set all elements of etest8 to one,
!        which makes the error test stricter.  (The elements of etest8
!        may have already been set to one earlier in this routine.)
! 
         IF (linear) CALL DLOAD(ntol, one, etest8, 1)
! 
!        Selectively refine the mesh using the old solution and the
!        6th order deferred correction.  Then, for a nonlinear prpblem,
!        interpolate the old solution onto the new mesh.
! 
         ipow = 4
! 
!        The array def6 will be overwritten by selmsh.
         IF (use_c) THEN
            IF (.NOT. stiff_cond ) THEN
               CALL SELMSH(ncomp, nmsh, ntol, ltol, tol, nfxpnt, fixpnt, ipow, nmax,  &
                           xx, ncomp, uold, def6, irefin, ihcomp,  &
                           nmold, xxold, ermx, ddouble , maxmsh)
            ELSE
               CALL SELCONDERRMSH(ncomp, nmsh, ntol, ltol, tol,  &
                                  nfxpnt, fixpnt, ipow, nmax, xx, ncomp, uold, def6, irefin, ihcomp,  &
                                  nmold, xxold, ermx, ddouble , maxmsh,r4,amg,stab_cond)
            END IF
         ELSE
            CALL SELMSH(ncomp, nmsh, ntol, ltol, tol, nfxpnt, fixpnt, ipow, nmax,  &
                        xx, ncomp, uold, def6, irefin, ihcomp,  &
                        nmold, xxold, ermx, ddouble , maxmsh)
         END IF
         IF (.NOT. maxmsh) THEN
            IF (linear) THEN
               CALL INITU(ncomp, nmsh, xx, nudim, u, rpar,ipar)
            ELSE
               CALL INTERP(ncomp, nmsh, xx, nudim, u, nmold, xxold, uold)
            END IF
         END IF
!
!     end of logic for err8 greater than err6
      END IF
!
      IF (pdebug .AND. .NOT.succes) WRITE(6,903)
!
901   FORMAT(1H ,'conv8')
902   FORMAT(1H ,'err6, err8, er6old, er8old',4(1PE11.3))
903   FORMAT(1H ,'8th order fails error tests.')
!
      RETURN
   END SUBROUTINE CONV8
!!!
!!!
   SUBROUTINE FAIL8(ncomp, nmsh, nfxpnt, fixpnt, nmax, ntol, ltol, tol, nmold,  &
                    xx, nudim, u, def6, xxold, uold, ihcomp, irefin, ermx, ddouble , maxmsh,  &
                    r4,amg, stiff_cond, stab_cond)
!
      USE BVPShared, ONLY: pdebug, use_c, comp_c, uval0,  &
                           nminit, iprint, idum
!
      IMPLICIT NONE
!
      INTEGER(i4), INTENT(IN)    :: ncomp
      INTEGER(i4), INTENT(INOUT) :: nmsh
      INTEGER(i4), INTENT(IN)    :: nfxpnt
      REAL(r8),    INTENT(IN)    :: fixpnt(*)
      INTEGER(i4), INTENT(INOUT) :: nmax
      INTEGER(i4), INTENT(IN)    :: ntol
      INTEGER(i4), INTENT(IN)    :: ltol(ntol)
      REAL(r8),    INTENT(IN)    :: tol(ntol)
      INTEGER(i4), INTENT(INOUT) :: nmold
      REAL(r8),    INTENT(INOUT) :: xx(*)
      INTEGER(i4), INTENT(IN)    :: nudim
      REAL(r8),    INTENT(INOUT) :: u(nudim,*)
      REAL(r8),    INTENT(INOUT) :: def6(ncomp,*)
      REAL(r8),    INTENT(INOUT) :: xxold(*)
      REAL(r8),    INTENT(INOUT) :: uold(ncomp,*)
      INTEGER(i4), INTENT(INOUT) :: ihcomp(*)
      INTEGER(i4), INTENT(INOUT) :: irefin(*)
      REAL(r8),    INTENT(INOUT) :: ermx(*)
      LOGICAL,     INTENT(INOUT) :: ddouble
      LOGICAL,     INTENT(INOUT) :: maxmsh
      REAL(r8),    INTENT(INOUT) :: r4(*)
      REAL(r8),    INTENT(INOUT) :: amg(*)
      LOGICAL,     INTENT(INOUT) :: stiff_cond
      LOGICAL,     INTENT(INOUT) :: stab_cond
!
      INTEGER(i4) :: ipow
!!!
!
      IF (pdebug) WRITE(6,901)
!
!     8th order solution did not converge (the problem must be nonlinear)
!
      ipow = 4
!
!     Selectively refine the mesh based on the 6th order deferred
!     correction and the old solution.
!
!     The def6 array is overwritten by selmsh.
!
      IF (use_c .AND. stiff_cond) THEN
         CALL SELCONDERRMSH(ncomp, nmsh, ntol, ltol, tol,  &
                            nfxpnt, fixpnt, ipow, nmax, xx, ncomp, uold, def6, irefin, ihcomp,  &
                            nmold, xxold, ermx, ddouble , maxmsh, r4, amg, stab_cond)
      ELSE
         CALL SELMSH(ncomp, nmsh, ntol, ltol, tol, nfxpnt, fixpnt, ipow, nmax,  &
                     xx, ncomp, uold, def6, irefin, ihcomp, nmold, xxold, ermx, ddouble, maxmsh)
      END IF
!
!     Interpolate to obtain the new initial solution.
!
      IF (.NOT. maxmsh) THEN
         CALL INTERP(ncomp, nmsh, xx, nudim, u, nmold, xxold, uold)
      END IF
!
901 FORMAT(1H ,'fail8')
!
      RETURN
   END SUBROUTINE FAIL8
!!!
!!!
   SUBROUTINE DCCAL(ncomp, nmsh, ntol, ltol, defexp, defimp, dfctol, fval,  &
                    ratdc, dfexmx, incmp, inmsh, intol, derivm, dfimmx,  &
                    rat1, rat2)
!
      USE BVPShared, ONLY: pdebug, use_c, comp_c, uval0,  &
                           nminit, iprint, idum
!
      IMPLICIT NONE
!
      INTEGER(i4), INTENT(IN)    :: ncomp
      INTEGER(i4), INTENT(INOUT) :: nmsh
      INTEGER(i4), INTENT(IN)    :: ntol
      INTEGER(i4), INTENT(IN)    :: ltol(ntol)
      REAL(r8),    INTENT(IN)    :: defexp(ncomp,nmsh-1)
      REAL(r8),    INTENT(IN)    :: defimp(ncomp,nmsh-1)
      REAL(r8),    INTENT(IN)    :: dfctol
      REAL(r8),    INTENT(IN)    :: fval(ncomp,nmsh)
      REAL(r8),    INTENT(INOUT) :: ratdc(nmsh-1)
      REAL(r8),    INTENT(INOUT) :: dfexmx
      INTEGER(i4), INTENT(INOUT) :: incmp
      INTEGER(i4), INTENT(INOUT) :: inmsh
      INTEGER(i4), INTENT(INOUT) :: intol
      REAL(r8),    INTENT(INOUT) :: derivm
      REAL(r8),    INTENT(INOUT) :: dfimmx
      REAL(r8),    INTENT(INOUT) :: rat1
      REAL(r8),    INTENT(INOUT) :: rat2
!
      REAL(r8), PARAMETER :: zero = 0.0D+0
      REAL(r8), PARAMETER :: one = 1.0D+0
      REAL(r8), PARAMETER :: rtst = 50.0D+0
      REAL(r8), PARAMETER :: tstrat = 0.1D+0
!
      REAL(r8) :: dval, smtest, texp, timp, abtexp, abrat
      INTEGER(i4) :: it, icmp, idmx, im
!
!     Find dfexmx, the maximum-magnitude element of defexp
!     in components for which a tolerance is specified.
!     The component index of dfexmx (incmp), its mesh
!     interval index (inmsh), and its tolerance index (intol),
!     are output parameters.
!
      dfexmx = zero
      DO it = 1, ntol
         icmp = ltol(it)
         idmx = IDAMAX(nmsh-1, defexp(icmp, 1), ncomp)
         dval = ABS(defexp(icmp, idmx))
         IF (dval >= dfexmx) THEN
            dfexmx = dval
            incmp = icmp
            inmsh = idmx
            intol = it
         END IF
      END DO
!
      IF (pdebug) THEN
         WRITE(6,901)
         WRITE(6,902) dfexmx, incmp, inmsh, intol
      END IF
!
!     Find derivm (maximum-magnitude element of fval(incmp,*))
!     for all mesh points.
!
      idmx = IDAMAX(nmsh, fval(incmp, 1), ncomp)
      derivm = ABS(fval(incmp, idmx))
      IF (pdebug) WRITE(6,903) derivm
!
!     For component incmp, go through the mesh intervals to calculate
!     (1) dfimmx, the maximum implicit deferred correction;
!     (2) two crucial ratios, rat1 and rat2, used in deciding whether
!         to refine the mesh;
!     (3) the array ratdc of deferred-correction ratios (explicit to
!         implicit).
!
!     In defining rat1 and rat2, we consider only intervals for
!     which the explicit deferred correction (defexp) exceeds the
!     tolerance dfctol in magnitude.  If it does not, the associated
!     interval does not affect rat1 or rat2, and the value of ratdc
!     is taken as 1.
!
!     rat2 is the maximum-magnitude ratio of sufficiently large
!     defexp to the larger of defimp and dfctol.
!
!     rat1 is the maximum-magnitude ratio of sufficiently large
!     defexp to the larger in magnitude of defimp and dfctol, but only
!     for those values of defexp greater than tstrat*dfexmx.
!     Thus by construction rat1 is less than or equal to rat2.
!
      rat1 = zero
      rat2 = zero
      dfimmx = zero
      smtest = tstrat*dfexmx
!
      DO im = 1, nmsh-1
         texp = defexp(incmp, im)
         timp = defimp(incmp, im)
         dfimmx = MAX(dfimmx, ABS(timp))
         abtexp = ABS(texp)
         IF (abtexp <= dfctol) THEN
            ratdc(im) = one
         ELSE
            IF (ABS(timp) < dfctol) timp = dfctol
            ratdc(im) = texp/timp
            abrat = ABS(ratdc(im))
            rat2 = MAX(rat2, abrat)
            IF (ABS(texp) >= smtest .AND. abrat >= rat1)  rat1 = abrat
         END IF
      END DO
!
      IF (pdebug) WRITE(6,905) rat1, rat2, dfimmx
!
901   FORMAT(1H ,'dccal')
902   FORMAT(1H ,'dfexmx, incmp, inmsh, intol',1PE11.3,3I5)
903   FORMAT(1H ,'derivm',1PE11.3)
904   FORMAT(1H ,'im, texp, timp, ratdc, dfctol',i5,4(1PE11.3))
905   FORMAT(1H ,'rat1, rat2, dfimmx', 3(1PE11.3))
!
      RETURN
   END SUBROUTINE DCCAL
!!!
!!!
   SUBROUTINE DECID4(linear, rat1, rat2, dfexmx, dfimmx,  &
                     derivm, dfold, tolval, oldrt1,  &
                     onto6, smooth, callrt, strctr, oscchk, ddouble , reposs)
!
      USE BVPShared, ONLY: pdebug, use_c, comp_c, uval0,  &
                           nminit, iprint, idum
!
      IMPLICIT NONE
!
      LOGICAL,  INTENT(IN)    :: linear
      REAL(r8), INTENT(INOUT) :: rat1
      REAL(r8), INTENT(INOUT) :: rat2
      REAL(r8), INTENT(INOUT) :: dfexmx
      REAL(r8), INTENT(INOUT) :: dfimmx
      REAL(r8), INTENT(INOUT) :: derivm
      REAL(r8), INTENT(INOUT) :: dfold
      REAL(r8), INTENT(IN)    :: tolval
      REAL(r8), INTENT(INOUT) :: oldrt1
      LOGICAL,  INTENT(INOUT) :: onto6
      LOGICAL,  INTENT(INOUT) :: smooth
      LOGICAL,  INTENT(INOUT) :: callrt
      LOGICAL,  INTENT(INOUT) :: strctr
      LOGICAL,  INTENT(INOUT) :: oscchk
      LOGICAL,  INTENT(INOUT) :: ddouble
      LOGICAL,  INTENT(INOUT) :: reposs
!
      REAL(r8), PARAMETER :: tenth = 0.1D+0
      REAL(r8), PARAMETER :: one = 1.0D+0
      REAL(r8), PARAMETER :: two = 2.0D+0
      REAL(r8), PARAMETER :: thrtwo = 32.0D+0
      REAL(r8), PARAMETER :: rtst = 50.0D+0
      REAL(r8), PARAMETER :: derval = 50.0D+0
!
      REAL(r8) :: thttol
      LOGICAL :: stest
!!!
!
!     decid4 evaluates information about the deferred corrections
!     and the nature of the problem, and sets various logical
!     flags that control the subsequent logic of the algorithm.
!
!     This code has been written for clarity, and NOT for efficiency.
!
      onto6 = .true.
      callrt = .false.
      smooth = .false.
      oscchk = .false.
      strctr = .false.
      reposs = .false.
      ddouble = .false.
!
!     rat2 is always greater than or equal to rat1.
!
      IF (pdebug) WRITE(6,901)
      IF (pdebug) WRITE(6,902) tolval, rtst
!
      stest = .true.
      IF (linear) stest = dfexmx < tenth*dfold
!
      IF (rat2 < rtst) THEN
         IF (stest) THEN
            smooth = .true.
         ELSE
            oscchk = .true.
         END IF
         RETURN
      END IF
!
!     We know now that rat2 .ge. rtst.
!
      thttol = thrtwo*tolval
      IF (pdebug) WRITE(6,903) thttol
!
      IF (rat1 < rtst .AND. dfexmx < thttol) THEN
         IF (stest) THEN
            smooth = .true.
         ELSE
            oscchk = .true.
         END IF
         RETURN
      END IF
!
      IF (rat1 < rtst .AND. dfexmx >= thttol) THEN
         callrt = .true.
         RETURN
      END IF
!
!     We know now that rat1 .ge. rtst (and that rat2 .ge. rtst).
!
      IF (derivm > derval .AND. dfexmx < thttol) THEN
         IF (stest) THEN
            smooth = .true.
         ELSE
            oscchk = .true.
         END IF
         RETURN
      END IF
!
      IF (derivm > derval .AND. dfexmx > thttol) THEN
         IF (dfimmx < one) THEN
            callrt = .true.
         ELSE
            strctr = .true.
            IF (linear) THEN
               onto6 = .false.
               IF (two*rat1 >= oldrt1) ddouble = .true.
!              end of logic for linear
            END IF
!           end of logic for dfimmx .ge. one
         END IF
         RETURN
!        end of logic for derivm .gt. derval .and dfexmx .gt. thttol
      END IF
!
!     To reach this point in the code, both of the following must hold:
!     rat1 .ge. rtst (which means that rat2 .ge. rtst)
!     derivm .le. derval
!
!     On linear problems, a special termination rule is tried if two
!     conditions hold:
!       (1) the 4th order solution has been computed on two consecutive
!           meshes, one of which is the double of the other
!           (this holds when ddouble  is .true.), and
!       (2) on both meshes, rat1 .ge. rtst, and derivm is small.  When
!           the conditions in (2) hold for a particular mesh, decid4
!           sets reposs to .true. to indicate that Richardson
!           extrapolation may be possible.
!     This set of tests is to take care of transients kept out by
!     initial conditions.
!
      IF (linear) reposs = .true.
!
901   FORMAT(1H ,'decid4')
902   FORMAT(1H ,'tolval, rtst',2(1PE11.3))
903   FORMAT(1H ,'thttol',1PE11.3)
!
      RETURN
   END SUBROUTINE DECID4
!!!
!!!
   SUBROUTINE DFEXCL(ncomp, nmsh, xx, nudim, u, defexp, fval,  &
                     tmp, rpar,ipar)
!
      USE BVPShared, ONLY: pdebug, use_c, comp_c, uval0,  &
                           nminit, iprint, idum, &
                           alp1, alp2, alp3, bet0, bet2, bet3, bet4,  &
                           a1, b1, c1, d1, e2, f2, c2, d2, b2, a2,  &
                           p3, q3, e3, f3, c3, d3, a3, b3,  &
                           a4, p4, x4, e4, c4,  &
                           a5, b5, c5, d5, e5, f5, a6, b6, c6
      USE BVPExtern, ONLY: FSUB
!
      IMPLICIT NONE
!
      INTEGER(i4), INTENT(IN)    :: ncomp
      INTEGER(i4), INTENT(INOUT) :: nmsh
      REAL(r8),    INTENT(IN)    :: xx(nmsh)
      INTEGER(i4), INTENT(IN)    :: nudim
      REAL(r8),    INTENT(INOUT) :: u(nudim,nmsh)
      REAL(r8),    INTENT(INOUT) :: defexp(ncomp,nmsh-1)
      REAL(r8),    INTENT(INOUT) :: fval(ncomp, nmsh)
      REAL(r8),    INTENT(INOUT) :: tmp(ncomp,4)
      REAL(r8),    INTENT(INOUT) :: rpar(*)
      INTEGER(i4), INTENT(INOUT) :: ipar(*)
!
      REAL(r8), PARAMETER :: half = 0.5D+0
      REAL(r8), PARAMETER :: fourth = 0.25D+0
      REAL(r8), PARAMETER :: thfrth= 0.75D+0
!
      REAL(r8) :: hmsh
      INTEGER(i4) :: im, ic
!!!
!
!     Given the nmsh mesh points xx, the estimated solution
!     u and the array fval of function values at (xx(im), u(*,im)),
!     im = 1,...,nmsh, dfexcl calculates sixth-order explicit
!     deferred correction, stored in the array defexp, indexed
!     over the components and mesh intervals.
!
!     The array tmp is workspace for 4 intermediate vectors of
!     dimension ncomp.
!
      DO im = 1, nmsh-1
!
         hmsh = xx(im+1) - xx(im)
         DO ic = 1, ncomp
            tmp(ic,1) = (a5*u(ic, im+1) + b5*u(ic, im))  &
               + hmsh*(c5*fval(ic,im)- d5*fval(ic,im+1))
            tmp(ic,2) = (b5*u(ic,im+1) + a5*u(ic,im))  &
               + hmsh*(-c5*fval(ic,im+1) + d5*fval(ic,im))
         END DO
!        print *, tmp
! 
         CALL FSUB(ncomp, xx(im)+fourth*hmsh, tmp(1,1), tmp(1,3), rpar, ipar)
         CALL FSUB(ncomp, xx(im)+thfrth*hmsh, tmp(1,2), tmp(1,4), rpar, ipar)
! 
!        print *, tmp
         DO ic = 1, ncomp
            tmp(ic,1) = half*(u(ic,im+1) + u(ic,im))  &
               + e5*hmsh*(fval(ic,im+1) - fval(ic,im)) - f5*hmsh*(tmp(ic,4) - tmp(ic,3))
!           print *, half*(u(ic,im+1) + u(ic,im))
!           print *,  e5*hmsh*(fval(ic,im+1) - fval(ic,im))
!           print *, - f5*hmsh*(tmp(ic,4) - tmp(ic,3))
         END DO
!        print *, tmp
! 
         CALL FSUB(ncomp, half*(xx(im)+xx(im+1)), tmp(1,1), tmp(1,2), rpar, ipar)
!        print *, tmp
         DO ic = 1, ncomp
            defexp(ic,im) = hmsh*(a6*(fval(ic,im+1) + fval(ic,im))  &
               + b6*(tmp(ic,3) + tmp(ic,4)) + c6*tmp(ic,2)) - u(ic,im+1) + u(ic,im)
         END DO
!        print *, (defexp(ic,im),ic=1,ncomp)
!        read(5,*)
!
      END DO
!
!
      RETURN
   END SUBROUTINE DFEXCL
!!!
!!!
   SUBROUTINE DF8CAL(ncomp, nmsh, xx, nudim, u, fval, def8, tmp, rpar,ipar)
!
!   Given the mesh points xx, the solution u, and the function
!   values fval, df8cal computes eighth-order deferred corrections,
!   which are stored in def8.
!   The array tmp is workspace for 8 intermediate vectors.
!
      USE BVPShared, ONLY: pdebug, use_c, comp_c, uval0,  &
                           nminit, iprint, idum, &
                           alp1, alp2, alp3, bet0, bet2, bet3, bet4,  &
                           a1, b1, c1, d1, e2, f2, c2, d2, b2, a2,  &
                           p3, q3, e3, f3, c3, d3, a3, b3,   &
                           a4, p4, x4, e4, c4,  &
                           a5, b5, c5, d5, e5, f5, a6, b6, c6
      USE BVPExtern, ONLY: FSUB
!
      IMPLICIT NONE
!
      INTEGER(i4), INTENT(IN)    :: ncomp
      INTEGER(i4), INTENT(INOUT) :: nmsh
      REAL(r8),    INTENT(INOUT) :: xx(nmsh)
      INTEGER(i4), INTENT(IN)    :: nudim
      REAL(r8),    INTENT(INOUT) :: u(nudim,nmsh)
      REAL(r8),    INTENT(INOUT) :: fval(ncomp,nmsh)
      REAL(r8),    INTENT(INOUT) :: def8(ncomp,nmsh-1)
      REAL(r8),    INTENT(INOUT) :: tmp(ncomp,8)
      REAL(r8),    INTENT(INOUT) :: rpar(*)
      INTEGER(i4), INTENT(INOUT) :: ipar(*)
!
      REAL(r8), PARAMETER :: half = 0.5D+0
      REAL(r8), PARAMETER :: two = 2.0D+0
      REAL(r8), PARAMETER :: fc1 = 0.625D+0
      REAL(r8), PARAMETER :: fc2= 0.375D+0
!
      REAL(r8) :: hmsh
      INTEGER(i4) :: im, ic
!!!
!
      DO im = 1, nmsh-1
! 
         hmsh = xx(im+1) - xx(im)
! 
         DO ic = 1, ncomp
            tmp(ic,1) = a1*u(ic,im+1) + b1*u(ic,im)  &
               + hmsh*(c1*fval(ic,im+1) + d1*fval(ic,im))
            tmp(ic,2) = b1*u(ic,im+1) + a1*u(ic,im)  &
               - hmsh*(c1*fval(ic,im) + d1*fval(ic,im+1))
         END DO
! 
         CALL FSUB(ncomp, xx(im)+fc1*hmsh, tmp(1,1), tmp(1,3), rpar, ipar)
         CALL FSUB(ncomp, xx(im)+fc2*hmsh, tmp(1,2), tmp(1,4), rpar, ipar)
         DO ic = 1, ncomp
            tmp(ic,1) = a2*u(ic,im+1) + b2*u(ic,im)  &
               + hmsh*(c2*fval(ic,im+1) + d2*fval(ic,im) + e2*tmp(ic,3) + f2*tmp(ic,4))
            tmp(ic,2) = b2*u(ic,im+1) + a2*u(ic,im)  &
               - hmsh*(d2*fval(ic,im+1) + c2*fval(ic,im) + f2*tmp(ic,3) + e2*tmp(ic,4))
         END DO
! 
         CALL FSUB(ncomp, xx(im)+(half+alp2)*hmsh, tmp(1,1), tmp(1,5), rpar, ipar)
         CALL FSUB(ncomp, xx(im)+(half-alp2)*hmsh, tmp(1,2), tmp(1,6), rpar, ipar)
         DO ic = 1, ncomp
            tmp(ic,1) = a3*u(ic,im+1) + b3*u(ic,im)  &
               + hmsh*(c3*fval(ic,im+1) + d3*fval(ic,im) + e3*tmp(ic,3) + f3*tmp(ic,4)  &
               + p3*tmp(ic,5) + q3*tmp(ic,6))
            tmp(ic,2) = b3*u(ic,im+1) + a3*u(ic,im)  &
               - hmsh*(d3*fval(ic,im+1) + c3*fval(ic,im) + f3*tmp(ic,3) + e3*tmp(ic,4)  &
               + q3*tmp(ic,5) + p3*tmp(ic,6))
         END DO
! 
         CALL FSUB(ncomp, xx(im)+(half+alp3)*hmsh, tmp(1,1), tmp(1,7), rpar, ipar)
         CALL FSUB(ncomp, xx(im)+(half-alp3)*hmsh, tmp(1,2), tmp(1,8), rpar, ipar)
         DO ic = 1, ncomp
            tmp(ic,1) = a4*(u(ic,im+1) + u(ic,im))  &
               + hmsh*(c4*(fval(ic,im+1) - fval(ic,im)) + e4*(tmp(ic,3) - tmp(ic,4))  &
               + x4*(tmp(ic,7) - tmp(ic,8)))
         END DO
! 
         CALL FSUB(ncomp, xx(im)+half*hmsh, tmp(1,1), tmp(1,2), rpar, ipar)
         DO ic = 1, ncomp
            def8(ic,im) = hmsh*(bet0*(fval(ic,im) + fval(ic,im+1))  &
               + bet2*(tmp(ic,5) + tmp(ic,6)) + bet3*(tmp(ic,7) + tmp(ic,8))  &
               + two*bet4*tmp(ic,2)) - u(ic,im+1) + u(ic,im)
         END DO
! 
      END DO
!
      RETURN
   END SUBROUTINE df8cal
!!!
!!!
   SUBROUTINE DFIMCL( ncomp, nmsh, defexp, chold, dsq, dexr, ipivot, defimp)
!
      USE BVPShared, ONLY: pdebug, use_c, comp_c, uval0, nminit,  &
                           iprint, idum
!
      IMPLICIT NONE
!
      INTEGER(i4), INTENT(IN)    :: ncomp
      INTEGER(i4), INTENT(INOUT) :: nmsh
      REAL(r8),    INTENT(INOUT) :: defexp(ncomp,nmsh-1)
      REAL(r8),    INTENT(INOUT) :: chold(ncomp,ncomp,nmsh-1)
      REAL(r8),    INTENT(INOUT) :: dsq(ncomp,ncomp)
      REAL(r8),    INTENT(INOUT) :: dexr(ncomp)
      INTEGER(i4), INTENT(INOUT) :: ipivot(ncomp)
      REAL(r8),    INTENT(INOUT) :: defimp(ncomp,nmsh-1)
!
      REAL(r8), PARAMETER :: zero = 0.0D+0
      INTEGER(i4) :: im, ic, ierlu
!
!     dfimcl calculates the rational deferred correction array,
!     which is indexed over the components and mesh intervals.
!
      CALL MTLOAD(ncomp, nmsh-1, zero, ncomp, defimp)
!
      DO im = 1, nmsh-1
         CALL DCOPY(ncomp, defexp(1,im), 1, dexr(1), 1)
         DO  ic = 1, ncomp
            CALL DCOPY(ncomp, chold(1,ic,im), 1, dsq(1,ic), 1)
         END DO
         CALL LUFAC(ncomp, ncomp, dsq, ipivot, ierlu)
         IF (ierlu == 0) THEN
            CALL LUSOL(ncomp, ncomp, dsq, ipivot, dexr, defimp(1, im))
         END IF
      END DO
!
      RETURN
   END SUBROUTINE DFIMCL
!!!
!!!
   SUBROUTINE OSC(ncomp, nmsh, dfexmx, incmp,  &
                  defcor, ratdc, ddouble , inmsh, onto6, trst6, smooth)
!
      USE BVPShared, ONLY: pdebug, use_c, comp_c, uval0,  &
                           nminit, iprint, idum
!
      IMPLICIT NONE
! 
      INTEGER(i4), INTENT(IN)    :: ncomp
      INTEGER(i4), INTENT(INOUT) :: nmsh
      REAL(r8),    INTENT(INOUT) :: dfexmx
      INTEGER(i4), INTENT(INOUT) :: incmp
      REAL(r8),    INTENT(INOUT) :: defcor(ncomp,*)
      REAL(r8),    INTENT(INOUT) :: ratdc(*)
      LOGICAL,     INTENT(INOUT) :: ddouble
      INTEGER(i4), INTENT(INOUT) :: inmsh
      LOGICAL,     INTENT(INOUT) :: onto6
      LOGICAL,     INTENT(INOUT) :: trst6
      LOGICAL,     INTENT(INOUT) :: smooth
!
      REAL(r8), PARAMETER :: zero = 0.0D+0
      REAL(r8), PARAMETER :: half = 0.5D+0
      REAL(r8), PARAMETER :: frac1 = 0.1D+0
      REAL(r8), PARAMETER :: frac2 = 1.0D-2
!
      REAL(r8) :: rmax, allsum, smlsum, bigsum, abdef
      REAL(r8) :: avsm, avbg, ave
      INTEGER(i4) :: jsndif, ninter, ibig, ism, im
!!!
!
!     For linear problems, subroutine osc performs heuristic tests
!     to detect an oscillating solution.  the tests check whether
!     (significant) explicit and implicit deferred corrections have
!     different (componentwise) signs.
!
!     dfexmx is the maximum-magnitude explicit deferred correction,
!     and is known to occur in component incmp.
!
!     The array defcor contains the explicit deferred corrections.
!     The array ratdc contains the ratios of explicit to implicit
!     deferred corrections.
!
!     jsndif counts the number of differences in sign.
!
      jsndif = 0
      rmax = zero
!
!     allsum is the sum of the magnitudes of all deferred corrections,
!     smlsum is the sum of the magnitudes of the small deferred
!     corrections, and bigsum is the sum of the magnitudes of the
!     large deferred corrections.  Here, small is defined as less
!     than half of the maximum.
!
      ninter = nmsh - 1
!
      IF (pdebug) WRITE(6,901)
!
      allsum = zero
      smlsum = zero
      bigsum = zero
      ibig = 0
      ism = 0
!
      DO im = 1, ninter
!
         abdef = ABS(defcor(incmp,im))
         allsum = allsum + abdef
         IF (abdef < half*dfexmx) THEN
            ism = ism + 1
            smlsum = smlsum + abdef
         ELSE
            ibig = ibig + 1
            bigsum = bigsum + abdef
         END IF
!
!        The counter of sign differences is incremented if (1) ratdc is negative
!        (which means that the two deferred corrections have opposite
!        sign) and (2) the explicit deferred correction is not too small
!        relative to the maximum.
! 
         IF (pdebug) WRITE(6,902) im, ratdc(im), abdef, frac2*dfexmx
! 
         IF (ratdc(im) < zero .AND. abdef >= frac2*dfexmx) THEN
!
            jsndif = jsndif + 1
!    
!           If more than 4 sign differences have occurred, exit after setting
!           ddouble to .true., which signals that the mesh
!           should be doubled (i.e., twice as many intervals).
!
            IF (jsndif > 4) THEN
               onto6 = .false.
               ddouble = .true.
               RETURN
            END IF
            IF (ABS(ratdc(im)) >= rmax) THEN
               rmax = ABS(ratdc(im))
               inmsh = im
            END IF
!
         END IF
!
      END DO
!
      IF (pdebug) WRITE(6,903) rmax, jsndif
!
      avsm = zero
      IF (ism > 0) avsm = smlsum/ism
      avbg = zero
      IF (ibig > 0) avbg = bigsum/ibig
      ave = allsum/ninter
!
      IF (pdebug) WRITE(6,904) ave, avsm, avbg
!
      IF (avsm > frac1*avbg .OR. ave > half*avbg) THEN
! 
!        The error appears to be uniformly large.
!        Signal that the 6th order solution should be calculated.
!
         onto6 = .true.
!
      ELSE IF (jsndif == 0) THEN
!
!        If there were no sign changes, the problem appears to be smooth.
!
         smooth = .true.
         onto6 = .true.
!
      ELSE
!
!        If the sign changed at between 1 and 4 points, don't go on to
!        6th order, and don't ever accept a 6th order solution even if the
!        error estimate at a later stage indicates that it is OK to do so.
!        Set ddouble to .false., to signal that the mesh will not necessarily
!        be doubled.
!
         ddouble = .false.
         onto6 = .false.
         trst6 = .false.
!
      END IF
!
901   FORMAT(1H ,'osc')
902   FORMAT(1H ,'im, ratdc, abdef, val',i5,3(1PE11.3))
903   FORMAT(1H ,'rmax, jsndif', 1PE11.3,i5)
904   FORMAT(1H ,'ave, avsm, avbg', 3(1PE11.3))
!
      RETURN
   END SUBROUTINE OSC
!!!
!!!
   SUBROUTINE RATCOR(ncomp, nmsh, xx, defimp, bhold, dfrat)
!
      IMPLICIT NONE
!
      INTEGER(i4), INTENT(IN)    :: ncomp
      INTEGER(i4), INTENT(INOUT) :: nmsh
      REAL(r8),    INTENT(INOUT) :: xx(nmsh)
      REAL(r8),    INTENT(INOUT) :: defimp(ncomp,nmsh-1)
      REAL(r8),    INTENT(INOUT) :: bhold(ncomp,ncomp,nmsh-1)
      REAL(r8),    INTENT(INOUT) :: dfrat(ncomp,nmsh-1)
!
      REAL(r8), PARAMETER :: zero = 0.0D+0
      REAL(r8), PARAMETER :: half = 0.5D+0
      REAL(r8), PARAMETER :: one = 1.0D+0
!
      REAL(r8) :: hmsh
      INTEGER(i4) :: ninter, im, ic
!!!
!
      ninter = nmsh - 1
      DO im = 1, ninter
        hmsh = xx(im+1) - xx(im)
        CALL DSCAL(ncomp*ncomp, (-half*hmsh), bhold(1,1,im), 1)
      END DO
!
      DO im = 1, ninter
        DO ic = 1, ncomp
          bhold(ic,ic,im) = bhold(ic,ic,im) + one
        END DO
      END DO
!
      DO im = 1, ninter
        DO  ic = 1, ncomp
          dfrat(ic,im) = ddot(ncomp, bhold(ic,1,im), ncomp, defimp(1,im), 1)
        END DO
      END DO
!
      RETURN
   END SUBROUTINE RATCOR
!!!
!!!
   SUBROUTINE STCONS
!
!  stcons computes constants needed in integration formulae
!  and stores them in a labeled common area.
!
      USE BVPShared, ONLY: alp1, alp2, alp3, bet0, bet2, bet3, bet4,  &
                           a1, b1, c1, d1, e2, f2, c2, d2, b2, a2,  &
                           p3, q3, e3, f3, c3, d3, a3, b3, a4, p4, x4, e4, c4,  &
                           a5, b5, c5, d5, e5, f5, a6, b6, c6
!
      IMPLICIT NONE
!
      REAL(r8), PARAMETER :: one = 1.0D+0, four = 4.0D+0, two = 2.0D+0
      REAL(r8), PARAMETER :: five = 5.0D+0, three = 3.0D+0
      REAL(r8), PARAMETER :: half = 0.5D+0, fourth = 0.25D+0
!
      REAL(r8) :: alp1sq, alp2sq, alp3sq
      REAL(r8) :: uu, vv, rr, ss, ww, z1, z2, u1, v1, r1, s1, w1
!!!
!
      alp1 = one/8.0D+0
      alp1sq = one/64.0D+0
      alp2sq = five/28.0D+0
      alp3sq = five/84.0D+0
      alp2 = SQRT(alp2sq)
      alp3 = SQRT(alp3sq)
      bet0 = one/6.0D+0 - (10.0D+0 - 28.0D+0*(alp2sq + alp3sq))/  &
          (105.0D+0*(one - four*alp2sq)* (one - four*alp3sq))
      bet2 = -(28.0D+0*alp3sq - three)/ (1680.0D+0*alp2sq*  &
          (one - four*alp2sq)*(alp2sq - alp3sq))
      bet3 = (28.0D+0*alp2sq - three)/ (1680.0D+0*alp3sq*  &
          (one - four*alp3sq)*(alp2sq - alp3sq))
      bet4 = half - bet0 - bet2-bet3
      a1 = half*(one - alp1)*(two*alp1 + one)* (two*alp1 + one)
      b1 = half*(one + alp1)* (two*alp1 - one)*(two*alp1 - one)
      c1 = half*(alp1sq - fourth)*(two*alp1 + one)
      d1 = half*(alp1sq - fourth)*(two*alp1 - one)
      uu = alp2*((four*alp2sq - one)**2)/  &
          ((four*alp1sq - one)*(20.0D+0*alp1*alp1 - one))
      vv = ((four*alp2sq - one)**2)/ (16.0D+0*alp1*(four*alp1sq - one))
      e2 = half*(uu + vv)
      f2 = half*(uu - vv)
      rr = half*(alp2*(four*alp2sq - one) + (one - 12.0D+0*alp1sq)*(e2 + f2))
      ss = fourth*(four*alp2sq - one) - two*alp1*(e2 - f2)
      c2 = half*(rr + ss)
      d2 = half*(rr - ss)
      ww = two*(alp2 - (c2 + d2 + e2 + f2))
      b2 = half*(one - ww)
      a2 = half*(one + ww)
      z1 = (three - 28.0D+0*alp3sq)/ (1680.0D+0*alp2*(four*alp2sq - one)*  &
          (alp2*alp2 - alp3sq)*bet3)
      z2 = one/(105.0D+0*alp3*bet3* (20.0D+0*alp2sq - one)*(four*alp2sq - one))
      p3 = half*(z1 + z2)
      q3 = half*(z2 - z1)
      u1 = (alp3*((four*alp3sq - one)**2)- (p3 + q3)*(20.0D+0*alp2sq - one)  &
          *(four*alp2sq - one))/ ((four*alp1sq - one)*(20.0D+0*alp1sq - one))
      v1 = (alp3sq*(one - two*alp3sq)- two*alp2*(one - four*alp2sq)*(p3 - q3)  &
          -one/8.0D+0)/(two*alp1*(one - four*alp1sq))
      e3 = half*(u1 + v1)
      f3 = half*(u1 - v1)
      r1 = half*(alp3*(four*alp3sq - one) + (e3 + f3)*(one - 12.0D+0*alp1sq) +  &
          (p3 + q3)*(one - 12.0D+0*alp2sq))
      s1 = alp3sq - fourth - two*alp1*(e3 - f3) - two*alp2*(p3 - q3)
      c3 = half*(r1 + s1)
      d3 = half*(r1 - s1)
      w1 = two*(alp3 - (c3 + d3 + e3 + f3 + p3 + q3))
      a3 = half*(one + w1)
      b3 = half*(one - w1)
      a4 = half
      p4 = 0.0D+0
      x4 = (three - 28.0D+0*alp2sq)/ (3360.0D+0*alp3*bet4*(four*alp3sq - one)  &
          *(alp3sq - alp2sq))
      e4 = (0.125D+0 + four*alp2*p4* (one - four*alp2sq) +  &
          four*alp3*x4*(one - four*alp3sq))/ (four*alp1*(four*alp1sq - one))
      c4 = -(0.125D+0 + two*alp1*e4 + two*alp2*p4 + two*alp3*x4)
      a5 = five/32.0D+0
      b5 = 27.0D+0/32.0D+0
      c5 = 9.0D+0/64.0D+0
      d5 = three/64.0D+0
      e5 = five/24.0D+0
      f5 = two/three
      a6 = 7.0D+0/90.0D+0
      b6 = 16.0D+0/45.0D+0
      c6 = two/15.0D+0
!
      RETURN
   END SUBROUTINE STCONS
!!!
!!!
   SUBROUTINE FIXJAC(ncomp, nmsh, nlbc, iorder, ntol, ltol, tol,  &
                     xx, nudim, u, defcor, defnew, delu, rhs, fval, utrial, rhstri,  &
                     rnsq, uint, ftmp, tmprhs, ajac, topblk, botblk, ipivot,  &
                     iflag,rpar,ipar)
!
      USE BVPShared, ONLY: pdebug, use_c, comp_c, uval0,  &
                           nminit, iprint, idum,  &
                           flmin, flmax, epsmch,  &
                           DCOPY, DSSQ, CRSLVE, MAXPY, MATCOP
      USE BVPExtern, ONLY: FSUB, GSUB
!
!     Fixed Jacobian iterations.
!
      IMPLICIT NONE
!
      INTEGER(i4), INTENT(IN)    :: ncomp
      INTEGER(i4), INTENT(INOUT) :: nmsh
      INTEGER(i4), INTENT(IN)    :: nlbc
      INTEGER(i4), INTENT(IN)    :: iorder
      INTEGER(i4), INTENT(IN)    :: ntol
      INTEGER(i4), INTENT(IN)    :: ltol(ntol)
      REAL(r8),    INTENT(IN)    :: tol(ntol)
      REAL(r8),    INTENT(INOUT) :: xx(nmsh)
      INTEGER(i4), INTENT(IN)    :: nudim
      REAL(r8),    INTENT(INOUT) :: u(nudim,nmsh)
      REAL(r8),    INTENT(INOUT) :: defcor(ncomp,nmsh-1)
      REAL(r8),    INTENT(INOUT) :: defnew(ncomp,nmsh-1)
      REAL(r8),    INTENT(INOUT) :: delu(ncomp,nmsh)
      REAL(r8),    INTENT(INOUT) :: rhs(ncomp*nmsh)
      REAL(r8),    INTENT(INOUT) :: fval(ncomp,nmsh)
      REAL(r8),    INTENT(INOUT) :: utrial(ncomp,nmsh)
      REAL(r8),    INTENT(INOUT) :: rhstri(ncomp*nmsh)
      REAL(r8),    INTENT(INOUT) :: rnsq
      REAL(r8),    INTENT(INOUT) :: uint(ncomp)
      REAL(r8),    INTENT(INOUT) :: ftmp(ncomp)
      REAL(r8),    INTENT(INOUT) :: tmprhs(ncomp*nmsh)
      REAL(r8),    INTENT(INOUT) :: ajac(ncomp,2*ncomp,nmsh-1)
      REAL(r8),    INTENT(INOUT) :: topblk(nlbc,ncomp)
      REAL(r8),    INTENT(INOUT) :: botblk(ncomp-nlbc,ncomp)
      INTEGER(i4), INTENT(INOUT) :: ipivot(ncomp*nmsh)
      INTEGER(i4), INTENT(INOUT) :: iflag
      REAL(r8),    INTENT(INOUT) :: rpar(*)
      INTEGER(i4), INTENT(INOUT) :: ipar(*)
!
      REAL(r8), PARAMETER :: one    = 1.0D+0
      REAL(r8), PARAMETER :: xlarge = 1.0D+6
      REAL(r8), PARAMETER :: huge = 1.0D+30
      REAL(r8), PARAMETER :: rngrow = 16.0D+0
      REAL(r8), PARAMETER :: rfact = 100.0D+0
      REAL(r8), PARAMETER :: tolfct = 0.1D+0
      INTEGER(i4), PARAMETER :: lmtfrz = 8
!
      REAL(r8) :: rnold, scale, sumsq, er
      INTEGER(i4) :: ninter, isize, ind, im, ic, iter, job, itol, it
      LOGICAL :: better
!
!     The iteration scheme uses a fixed Jacobian matrix to solve for
!     correction vectors, once there has been convergence of the Newton
!     iterations on this mesh.   It is assumed that the LU
!     factors of the Jacobian have been left unaltered since
!     their calculation.
!
      IF (iprint == 1) WRITE(6,901)
      ninter = nmsh - 1
      rnold = flmax
      isize = nmsh*ncomp
!
!     Evaluate the right-hand side rhstri at the initial solution u by
!     adding the new deferred corrections to the already-calculated
!     rhs vector.
!
      CALL DCOPY(nlbc, rhs, 1, rhstri, 1)
      ind = nlbc
      DO im = 1, ninter
         DO ic = 1, ncomp
            ind = ind + 1
            rhstri(ind) = rhs(ind) + defnew(ic, im)
         END DO
      END DO
      ind = ninter*nmsh + nlbc + 1
      CALL DCOPY(ncomp-nlbc, rhs, 1, rhstri, 1)
!
      CALL DSSQ(nmsh*ncomp, rhstri, 1, scale, sumsq)
      rnsq = (scale**2)*sumsq
!
      iter = 0
!
!     If the initial right-hand side is too large, do not even attempt to
!     solve the nonlinear equations.
!
      IF (rnsq > huge .OR. (iorder == 8 .AND. rnsq > xlarge)) THEN
         IF (iprint == 1) WRITE (6,902) rnsq
         iflag = -2
         RETURN
      END IF
      CALL DCOPY(ncomp*nmsh, rhstri, 1, rhs, 1)
!
!     Statement 100 is the top of the iteration loop.
!
100   CONTINUE
!
!     If rnsq is sufficiently small, terminate immediately.
!
      IF (iprint == 1) WRITE(6,903) iter, rnsq
      IF (rnsq <= epsmch) THEN
         iflag = 0
         RETURN
      END IF
!
      iter = iter + 1
!
!     Solve for the step delu by solving a system involving the fixed
!     Jacobian (whose LU factors are saved).  Copy the rhs array into
!     tmprhs, which is then overwritten by blkslv.
!
      CALL DCOPY(ncomp*nmsh, rhs, 1, tmprhs, 1)
      CALL DCOPY(ncomp*nmsh,tmprhs,1,delu,1)
      job = 0
      CALL CRSLVE(topblk, nlbc, ncomp, ajac, ncomp, 2*ncomp,  &
                  ninter, botblk, ncomp-nlbc, ipivot, delu, job)
!
!     Compute the trial point utrial by adding delu to u.
! 
      CALL MATCOP(nudim, ncomp, ncomp, nmsh, u, utrial)
      CALL MAXPY(ncomp, nmsh, one, delu, ncomp, utrial)
!
!     Compute the right-hand side vector rhstri and its squared
!     two-norm at the trial point.
!
      rnold = rnsq
      CALL FNEVAL(ncomp, nmsh, xx, ncomp, utrial, fval, rpar, ipar)
      CALL RHSCAL(ncomp, nmsh, nlbc, xx, ncomp, utrial, defcor,  &
                  rhstri, rnsq, fval, ftmp, uint,rpar,ipar)
!
!     If rnsq strictly decreased, update the solution vector u
!     and the right-hand side rhs.
!
      better = .false.
      IF (rnsq < rnold) THEN
         better = .true.
         CALL MATCOP(ncomp, nudim, ncomp, nmsh, utrial, u)
         CALL DCOPY(ncomp*nmsh, rhstri, 1, rhs, 1)
      END IF
!
!     Stop the fixed Jacobian iterations if there have been too
!     many iterations, or if rnsq has not decreased by a factor
!     of at least rngrow.
!
      IF (iter >= lmtfrz .OR. rnsq > (rnold/rngrow)) THEN
         IF (better) THEN
!    
!           Setting iflag to -3 signals that, although the fixed Jacobian
!           iterations did not succeed, the current point was an improvement
!           on the previous one.  Hence, if we switch to a Newton procedure,
!           the right-hand side does not need to be recalculated.
! 
            iflag = -3
         ELSE
            iflag = -2
         END IF
         IF (iprint == 1) WRITE(6,904) iflag
         RETURN
      END IF
!
!     Test for convergence using the ratio abs((change in u)/max(u,1)).
!
      DO im = 1, nmsh
         DO it = 1, ntol
            itol = ltol(it)
            er = ABS(delu(itol,im))/MAX(ABS(u(itol,im)), one)
            IF (er > tolfct*tol(it)) GO TO 100
         END DO
      END DO
!
!     To exit from the loop here, the convergence tests have
!     been passed.
!
      IF (iprint >= 0) WRITE(6,905) iter, rnsq
!
      iflag = 0
!
901   FORMAT(1H ,'fixed Jacobian iterations')
902   FORMAT(1H ,'Large residual, rnsq =',1PE12.4)
903   FORMAT(1H ,'iter, rnsq',i5,1PE11.3)
904   FORMAT(1H ,'failure of fixed Jacobian, iflag =',i5)
905   FORMAT(1H ,'fixed Jacobian convergence',i5,1PE11.3)
!
      RETURN
   END SUBROUTINE FIXJAC
!!!
!!!
   SUBROUTINE LINEQ(ncomp, nmsh, nlbc, ludone, xx, nudim, u, defcor,  &
                    delu, rhs, fval, uint, ftmp, dftmp1, dftmp2, dgtm, tmprhs,  &
                    ajac, topblk, botblk, bhold, chold, ipivot,  &
                    iflag, rpar, ipar)
!
      USE BVPShared, ONLY: pdebug, use_c, comp_c, uval0,  &
                           nminit, iprint, idum,  &
                           DCOPY, COLROW, DLOAD, CRSLVE, MAXPY
      USE BVPExtern, ONLY: FSUB, DFSUB, GSUB, DGSUB
!
      IMPLICIT NONE
!
      INTEGER(i4), INTENT(IN)    :: ncomp
      INTEGER(i4), INTENT(INOUT) :: nmsh
      INTEGER(i4), INTENT(IN)    :: nlbc
      LOGICAL,     INTENT(INOUT) :: ludone
      REAL(r8),    INTENT(INOUT) :: xx(nmsh)
      INTEGER(i4), INTENT(IN)    :: nudim
      REAL(r8),    INTENT(INOUT) :: u(nudim,nmsh)
      REAL(r8),    INTENT(INOUT) :: defcor(ncomp,nmsh-1)
      REAL(r8),    INTENT(INOUT) :: delu(ncomp,nmsh)
      REAL(r8),    INTENT(INOUT) :: rhs(ncomp*nmsh)
      REAL(r8),    INTENT(INOUT) :: fval(ncomp,nmsh)
      REAL(r8),    INTENT(INOUT) :: uint(ncomp)
      REAL(r8),    INTENT(INOUT) :: ftmp(ncomp)
      REAL(r8),    INTENT(INOUT) :: dftmp1(ncomp,ncomp)
      REAL(r8),    INTENT(INOUT) :: dftmp2(ncomp,ncomp)
      REAL(r8),    INTENT(INOUT) :: dgtm(ncomp)
      REAL(r8),    INTENT(INOUT) :: tmprhs(ncomp*nmsh)
      REAL(r8),    INTENT(INOUT) :: ajac(ncomp,2*ncomp,nmsh-1)
      REAL(r8),    INTENT(INOUT) :: topblk(nlbc,*)
      REAL(r8),    INTENT(INOUT) :: botblk(ncomp-nlbc,*)
      REAL(r8),    INTENT(INOUT) :: bhold(ncomp,ncomp,nmsh-1)
      REAL(r8),    INTENT(INOUT) :: chold(ncomp,ncomp,nmsh-1)
      INTEGER(i4), INTENT(INOUT) :: ipivot(ncomp*nmsh)
      INTEGER(i4), INTENT(INOUT) :: iflag
      REAL(r8),    INTENT(INOUT) :: rpar(*)
      INTEGER(i4), INTENT(INOUT) :: ipar(*)
!
      REAL(r8), PARAMETER :: one = 1.0D+0
      REAL(r8), PARAMETER :: zero = 0.0D+0
!
      REAL(r8) :: rnsq
      INTEGER(i4) :: isize, ninter, im, loc, nrhs
      INTEGER(i4) :: job
!!!
!
!     The routine lineq calculates the Newton step for a linear
!     problem.  The Newton step is exact unless the Jacobian
!     matrix is singular.
!
      isize = nmsh*ncomp
      ninter = nmsh - 1
      iflag = 0
!
      IF (.NOT. ludone) THEN
! 
!        Compute the right-hand side vector rhs.
!
         CALL LNRHS(ncomp, nmsh, nlbc, xx, nudim, u,  &
                    rhs, rnsq, fval, ftmp, uint, rpar, ipar)
! 
!        If the Jacobian for this mesh has not previously been
!        calulated and factorized successfully, call jaccal.
!        The block-structured Jacobian matrix is stored in three
!        matrices (topblk, ajac, and botblk).
!        The matrices bhold and chold are also calculated in jaccal,
!        and are saved for later use in outer routines.
! 
         CALL JACCAL(ncomp, nmsh, nlbc,  &
                     xx, nudim, u, fval, dgtm, dftmp1, dftmp2, uint,  &
                     ajac, topblk, botblk, bhold, chold, rpar, ipar)
! 
!        Call blkdcm to calculate the LU factors of the Jacobian.
!        The factors are overwritten on the matrices topblk, ajac and botblk.
!        Interchanges are represented in the integer array ipivot.
! 
         CALL DCOPY(ncomp*nmsh,rhs,1,tmprhs,1)
         CALL DCOPY(ncomp*nmsh,tmprhs,1,delu,1)
! 
         job = 0
         CALL COLROW(isize, topblk, nlbc, ncomp, ajac, ncomp, 2*ncomp,  &
                     ninter, botblk, ncomp-nlbc, ipivot, delu, iflag, job)
!
         ludone = .true.
! 
!        Copy the rhs into the temporary vector tmprhs, which will be
!        overwritten by blkslv.
!
      ELSE
! 
!        The right-hand side is the deferred correction array,
!        padded with zeros at the boundary conditions.
! 
         CALL DLOAD(nlbc, zero, tmprhs(1), 1)
         DO im = 1, ninter
            loc = (im-1)*ncomp + nlbc + 1
            CALL DCOPY(ncomp, defcor(1,im), 1, tmprhs(loc), 1)
         END DO
         nrhs = ninter*ncomp + nlbc + 1
         CALL DLOAD(ncomp-nlbc, zero, tmprhs(nrhs), 1)
         CALL DCOPY(ncomp*nmsh,tmprhs,1,delu,1)
         job = 0
! 
         CALL CRSLVE(topblk,nlbc,ncomp,ajac,ncomp,2*ncomp,ninter,  &
                     botblk,ncomp-nlbc,ipivot,delu,job)
! 
      END IF
!
!     Since the problem is linear, the Newton step  is exact.  The
!     new u array is obtained by adding delu to u.
!
      CALL MAXPY(ncomp, nmsh, one, delu, nudim, u)
!
!     iflag = 0
!
901   FORMAT(1H ,'Singular matrix')
!
      RETURN
   END SUBROUTINE LINEQ
!!!
!!!
   SUBROUTINE NEWTEQ(ncomp, nmsh, nlbc, rhsgiv, ntol, ltol, tol,  &
                     xx, nudim, u, defcor, delu, rhs, fval,  &
                     utrial, rhstri, uint,  &
                     ftmp, dftmp1, dftmp2, dgtm, tmprhs, xmerit,  &
                     ajac, topblk, botblk, bhold, chold, ipivot,  &
                     iter, iflag, rpar, ipar, frscal)
!
      USE BVPShared, ONLY: pdebug, use_c, comp_c, uval0,  &
                           nminit, iprint, idum,  &
                           flmin, flmax, epsmch,  &
                           DCOPY, COLROW, MSSQ, MATCOP, MAXPY
      USE BVPExtern, ONLY: FSUB, DFSUB, GSUB, DGSUB
!
      IMPLICIT NONE
!
      INTEGER(i4), INTENT(IN)    :: ncomp
      INTEGER(i4), INTENT(INOUT) :: nmsh
      INTEGER(i4), INTENT(IN)    :: nlbc
      LOGICAL,     INTENT(INOUT) :: rhsgiv
      INTEGER(i4), INTENT(IN)    :: ntol
      INTEGER(i4), INTENT(IN)    :: ltol(*)
      REAL(r8),    INTENT(IN)    :: tol(*)
      REAL(r8),    INTENT(INOUT) :: xx(*)
      INTEGER(i4), INTENT(IN)    :: nudim
      REAL(r8),    INTENT(INOUT) :: u(nudim,*)
      REAL(r8),    INTENT(INOUT) :: defcor(ncomp,*)
      REAL(r8),    INTENT(INOUT) :: delu(ncomp,*)
      REAL(r8),    INTENT(INOUT) :: rhs(ncomp*nmsh)
      REAL(r8),    INTENT(INOUT) :: fval(ncomp,*)
      REAL(r8),    INTENT(INOUT) :: utrial(ncomp,nmsh)
      REAL(r8),    INTENT(INOUT) :: rhstri(ncomp*nmsh)
      REAL(r8),    INTENT(INOUT) :: uint(*)
      REAL(r8),    INTENT(INOUT) :: ftmp(*)
      REAL(r8),    INTENT(INOUT) :: dftmp1(ncomp, ncomp)
      REAL(r8),    INTENT(INOUT) :: dftmp2(ncomp, ncomp)
      REAL(r8),    INTENT(INOUT) :: dgtm(ncomp)
      REAL(r8),    INTENT(INOUT) :: tmprhs(ncomp*nmsh)
      REAL(r8),    INTENT(INOUT) :: xmerit(ncomp,nmsh)
      REAL(r8),    INTENT(INOUT) :: ajac(ncomp,2*ncomp,*)
      REAL(r8),    INTENT(INOUT) :: topblk(nlbc,*)
      REAL(r8),    INTENT(INOUT) :: botblk(ncomp-nlbc,*)
      REAL(r8),    INTENT(INOUT) :: bhold(ncomp,ncomp,*)
      REAL(r8),    INTENT(INOUT) :: chold(ncomp,ncomp,*)
      INTEGER(i4), INTENT(INOUT) :: ipivot(*)
      INTEGER(i4), INTENT(INOUT) :: iter
      INTEGER(i4), INTENT(INOUT) :: iflag
      REAL(r8),    INTENT(INOUT) :: rpar(*)
      INTEGER(i4), INTENT(INOUT) :: ipar(*)
      LOGICAL,     INTENT(INOUT) :: frscal
!
      REAL(r8), PARAMETER :: zero   = 0.0D+0
      REAL(r8), PARAMETER :: one    = 1.0D+0
      REAL(r8), PARAMETER :: two = 2.0D+0
      REAL(r8), PARAMETER :: half = 0.5D+0
      REAL(r8), PARAMETER :: fourth = 0.25D+0
      REAL(r8), PARAMETER :: tenth = 0.1D+0
      REAL(r8), PARAMETER :: ten = 10.0D+0
      REAL(r8), PARAMETER :: hund = 100.0D+0
!
      REAL(r8)    :: alfsml = 1.0D-4
      REAL(r8)    :: alfmax = 1.1D+0
      INTEGER(i4) :: imerit = 1
      INTEGER(i4) :: lmtnwt = 39
      REAL(r8)    :: shrfct = 100.0D+0
      REAL(r8)    :: stpfct = 2.0D+0
      LOGICAL     :: gtpdeb = .false.
      INTEGER(i4) :: mfsrch = 5
      REAL(r8)    :: eta = .999999D+0
      REAL(r8)    :: rmu = 1.0D-6
!
      REAL(r8) :: epsaf, epsag, tolabs, xmscal, xmsq
      REAL(r8) :: tolrel, toltny, alfold, alfa, rnsq, rnprev, rnbest
      REAL(r8) :: fmtry, fa, oldg, fmold, alfuzz, xscale, xsolsq
      REAL(r8) :: alfbst, fbest, alin, blin, factor, fv, fw, xtry, xv, xw, rnsqtr, er
!
      INTEGER(i4) :: ninter, itwtch, iflwat, job, inform, iwr
      INTEGER(i4) :: nfsrch, nsamea, nsameb
      INTEGER(i4) :: im, it, icmp
!
      LOGICAL :: imprvd, braktd, crampd, extrap, vset, wset
      SAVE  gtpdeb, mfsrch, epsaf, epsag, eta, rmu, tolabs, alfmax
      SAVE  tolrel, toltny
!
      REAL(r8), PARAMETER :: cnvfct = 0.1D+0
!
!     The routine newteq performs Newton iterations with a line
!     search, to solve the nonlinear equations.
!
!     Set up constants if this is the first call to newteq.
!
      IF (frscal) THEN
         frscal = .false.
         epsaf = epsmch
         epsag = epsmch
         tolabs = epsmch
         tolrel = epsmch
         toltny = epsmch
      END IF
      ninter = nmsh - 1
!
      IF (iprint == 1) WRITE(6,901)
!
!     A Newton method with line search and watchdog safeguarding
!     is used to try to solve the nonlinear equations.
!
!     Initialize iter (the counter of Newton iterations) and alfold
!     (the step taken at the previous iteration).
!
      iter = -1
      alfold = one
      alfa = zero
!
      IF (.NOT. rhsgiv) THEN
!        If necessary, evaluate the right-hand side at the initial u.
         CALL RHSCAL(ncomp, nmsh, nlbc, xx, nudim, u, defcor,  &
                     rhs, rnsq, fval, ftmp, uint, rpar, ipar)
      END IF
!
!     At any given Newton iteration, rnprev is the value of rnsq at
!     the immediately preceding Newton iteration.
!
      rnprev = flmax
      rnbest = flmax
      IF (.NOT. pdebug .AND. iprint >= 0) WRITE (6,902)
!
!     Initialize counter of watchdog iterations.
!
      itwtch = 0
!
!     Statement 100 is the top of the Newton iteration loop.
!
100   CONTINUE
!
      iter = iter + 1
!
      IF (iprint == 1) WRITE(6,910) iter
!
!     If there have been too many Newton iterations, terminate.
!
      IF (iter >= lmtnwt) THEN
         IF (iprint >= 0) WRITE(6,903)
         iflag = -2
         RETURN
      END IF
!
!     The vector rhs is the right-hand side at the current iterate,
!     and rnsq is its squared two-norm.
!     Perform watchdog tests, using the unscaled merit function (rnsq)
!     as the watchdog function.  The routine wtchdg updates rnbest
!     and itwtch.  If iflwat is not zero, this sequence of Newton
!     iterations is terminated.
!
      iflwat = 0
!
      CALL WTCHDG(iter, rnsq, rnbest, rnprev, itwtch, alfold, iflwat)
!
      IF (iflwat /= 0) THEN
         IF (iprint >= 0) WRITE(6,904) iter
         iflag = -3
         RETURN
      END IF
!
!     Watchdog tests are passed.  Proceed with the Newton iteration.
!     Call jaccal to evaluate the block-structured Jacobian matrix,
!     which is stored in three matrices (topblk, ajac, and botblk).
!     The matrices bhold and chold are saved for use in later
!     calculations in the outer routine.
!
!     If rnsq is sufficiently small, terminate immediately.
!     Note that the stored Jacobian does not correspond exactly
!     to the final point.
!
      IF (rnsq <= epsmch) THEN
         IF (iprint >= 0)  WRITE(6,906) iter, rnsq
         iflag = 0
         RETURN
      END IF
!
      CALL JACCAL(ncomp, nmsh, nlbc,  &
                  xx, nudim, u, fval, dgtm, dftmp1, dftmp2, uint,  &
                  ajac, topblk, botblk, bhold, chold, rpar, ipar)
!
!     blkdcm is called to calculate the LU factors of the Jacobian,
!     which are overwritten on topblk, ajac and botblk.
!     Interchanges are represented in the integer array ipivot.
!
!     Solve for the Newton step delu.  Copy the rhs array into tmprhs,
!     which is then overwritten by blkslv.
!
      CALL DCOPY(ncomp*nmsh, rhs, 1, tmprhs, 1)
      CALL DCOPY(ncomp*nmsh,tmprhs,1,delu,1)
      job = 0
      CALL COLROW(nmsh*ncomp,topblk,nlbc,ncomp,ajac,ncomp,2*ncomp,  &
                  ninter,botblk,ncomp-nlbc,ipivot,delu,iflag,job)
!
      IF (iprint >= 0 .AND. iflag /= 0) WRITE(6,905) iter
      IF (iflag /= 0) RETURN
!
!     the jacobian is singular
!
!     If imerit = 1, the line search is based on the scaled merit function,
!     the squared two-norm of the solution xmerit of the linear system
!      (Jacobian)*xmerit = rhs,
!     where (Jacobian) is the Jacobian at the current Newton iterate.
!     Thus the initial value of the scaled merit function is simply
!     the squared two-norm of the Newton step delu itself.
!
      IF (imerit == 1) THEN
         CALL MSSQ(ncomp, nmsh, delu, xmscal, xmsq)
         fmtry = (xmscal**2)*xmsq
      ELSE
!        The unscaled merit function is simply the squared two-norm of rhs.
         fmtry = rnsq
      END IF
!
!     fa and oldg represent the merit function and its gradient
!     at the initial point of the line search.
!
      fa = fmtry
      oldg = -two*fa
      alfa = zero
      IF (iprint == 1) WRITE (6,908) alfa, fmtry, rnsq
!
!     On the first Newton iteration, the initial trial step is unity.
!     On subsequent iterations, the initial step is not allowed to
!     be more than the factor stpfct larger than the final step at
!     the immediately preceding iteration.
!
      alfa = one
      IF (stpfct*alfold < one) alfa = stpfct*alfold
!
      IF (alfa < alfsml) alfa = alfsml
!
      fmold = fa
      inform = -1
!
!     Statement 150 is the top of the inner line search iteration.
!     The line search routine getptq has been altered so that it
!     terminates with an indication of success as soon as a
!     strictly lower value of the merit function is found.  Note that
!     this is a much less strict requirement than the usual sufficient
!     decrease conditions.
!
150   CONTINUE
!
      iwr = 6
      CALL GETPTQ(gtpdeb, mfsrch, iwr, alfmax, alfsml, alfuzz, epsaf, epsag,  &
                  eta, fmtry, fmold, oldg, rmu, tolabs, tolrel, toltny,  &
                  imprvd, inform, nfsrch, alfa, alfbst, fbest,  &
                  braktd, crampd, extrap, vset, wset, nsamea, nsameb,  &
                  alin, blin, fa, factor, fv, fw, xtry, xv, xw)
!
!     inform = 1, 2 or 3 indicates success in finding an acceptable point.
!     inform = 4 means alfmax is too small (this should never happen here,
!     since alfmax is set always to 1.1).
!     inform = 5 means that a decrease was not achieved for any step
!     greater than alfsml.
!     inform = 6 means a better point could not be found (the minimum
!     probably lies too close to alfa=0).
!     inform = 7 means that the gradient at alfa=0 (oldg) is positive
!     (this cannot happen here, since oldg=-two*fa, and fa is a non-negative
!     number)
!
      IF (pdebug) WRITE(6,907) inform, alfa
!
      IF (inform == 5) THEN
         iflag = -5
         RETURN
      ELSE IF (inform == 4 .OR. inform == 7) THEN
         iflag = -4
         RETURN
      ELSE IF (inform == 0) THEN
! 
!        inform = 0 means that a new function value should be obtained
!        with the step alfa.
!        We may override alfa from getptq by requiring that the step is not
!        allowed to decrease by more than a factor of shrfct during
!        a line search iteration.
! 
         IF (alfa < alfold/shrfct) alfa = alfold/shrfct
         alfold = alfa
! 
!        Define the next iterate utrial = u + alfa*delu.
!        Call fneval and rhscal to evaluate the right-hand side
!        rhstri at utrial.
!        The vector rhstri is stored separately, and rhs is overwritten
!        only when an improved point is found.
! 
         CALL MATCOP(nudim, ncomp, ncomp, nmsh, u, utrial)
         CALL MAXPY (ncomp, nmsh, alfa, delu, ncomp, utrial )
         CALL FNEVAL(ncomp, nmsh, xx, ncomp, utrial, fval, rpar,ipar)
         CALL RHSCAL(ncomp, nmsh, nlbc, xx, ncomp, utrial, defcor,  &
                     rhstri, rnsqtr, fval, ftmp, uint,rpar,ipar)
! 
         fmold = fmtry
!
         IF (imerit == 1) THEN
!
!           Solve a linear system to obtain the 2-d array xmerit whose squared
!           norm is the scaled merit function.   The LU factors of the Jacobian
!           have already been calculated by blkdcm.
!           Copy rhstri into tmprhs, which is overwritten by blkslv.
!
            CALL DCOPY(ncomp*nmsh, rhstri, 1, tmprhs, 1)
            CALL DCOPY(ncomp*nmsh,tmprhs,1,xmerit,1)
            job = 0
            CALL CRSLVE(topblk, nlbc, ncomp, ajac, ncomp, 2*ncomp, ninter,  &
                        botblk, ncomp-nlbc, ipivot, xmerit,job)
            CALL MSSQ(ncomp, nmsh, xmerit, xscale, xsolsq )
            fmtry = (xscale**2)*xsolsq
!
         ELSE
!  
!           The unscaled merit function is the squared two-norm of the right-hand
!           side.
!
            fmtry = rnsqtr
!
         END IF
         IF (iprint == 1) WRITE (6,908) alfa, fmtry, rnsqtr
         GO TO 150
!
      END IF
!
!     To reach here, inform must be 1, 2, 3, or 6, and the line search
!     has found a strictly lower value of the merit function.
!     Store the new Newton iterate in u, and the corresponding rhs
!     vector in rhs.
!
      rnprev = rnsq
      rnsq = rnsqtr
      CALL MATCOP(ncomp, nudim, ncomp, nmsh, utrial, u)
      CALL DCOPY(ncomp*nmsh, rhstri, 1, rhs, 1)
      IF (iprint >= 0) WRITE(6,909) iter, alfa, fmtry, rnsq
!
!     Now test for convergence using the ratio of the Newton step
!     for each component with max(1, abs(current solution estimate)).
!     If the test fails for any element of u, branch back to the
!     top of the Newton iteration.
!
      DO im = 1, nmsh
         DO it = 1, ntol
            icmp = ltol(it)
            er = ABS(delu(icmp,im))/MAX(ABS(u(icmp,im)), one)
            IF (er > cnvfct*tol(it)) GO TO 100
         END DO
      END DO
!
      IF (iprint >= 0) WRITE(6, 906) iter+1, rnsq
      iflag = 0
!
!     To fall through the above loop, the termination test for a
!     sufficiently small delu is satisfied.
!     Note that the stored Jacobian and its factorization do not
!     correspond to the final solution.
!
901   FORMAT(1H ,'start Newton iterations')
902   FORMAT(1H ,' iter', 7X,'alfa',6X,'merit',7X,'rnsq')
903   FORMAT(1H ,'Too many Newton iterations')
904   FORMAT(1H ,'Watchdog tests fail, iter =', i5)
905   FORMAT(1H ,'Singular Jacobian, iter=',i5)
906   FORMAT(1H ,'Convergence, iter =',i5,4X,'rnsq =',1PE12.3)
907   FORMAT(1H ,'inform, alfa after getptq',i5,3X, 1PE11.3)
908   FORMAT(1H ,'alfa, merit, rnsq',3(1PE11.3))
909   FORMAT(1H ,i5,3(1PE11.3))
910   FORMAT(1H ,'Newton iteration',i5)
!
      RETURN
   END SUBROUTINE NEWTEQ
!!!
!!!
   SUBROUTINE WTCHDG(iter, wmerit, wmbest, wmprev, itwtch, alfold, iflag)
!
!     Logic for watchdog tests.
!
      IMPLICIT NONE
!
      INTEGER(i4), INTENT(INOUT) :: iter
      REAL(r8),    INTENT(INOUT) :: wmerit
      REAL(r8),    INTENT(INOUT) :: wmbest
      REAL(r8),    INTENT(INOUT) :: wmprev
      INTEGER(i4), INTENT(INOUT) :: itwtch
      REAL(r8),    INTENT(INOUT) :: alfold
      INTEGER(i4), INTENT(INOUT) :: iflag
!
      INTEGER(i4), PARAMETER :: itonew = 5
      INTEGER(i4), PARAMETER :: itwtmx = 8
      REAL(r8),    PARAMETER :: grfct = 80.0D+0
      REAL(r8),    PARAMETER :: half = 0.5D+0
!
!     Perform watchdog tests in two forms:
!     (1) to determine whether a sufficient decrease in the
!     watchdog merit function has occurred within the most recent
!     sequence of itwtmx iterations;
!     (2) to determine whether the watchdog merit function has increased
!     too much in a single iteration after itonew Newton iterations
!     have been performed.  This allows the merit function to increase
!     wildly only during the first itonew iterations.
!
!     wmbest is the smallest watchdog merit function achieved in this
!     sequence of Newton iterations.
!     wmprev is the watchdog merit function from the immediately
!     preceding Newton iteration.
!
!     itwtch counts the number of iterations without an improvement
!     in the unscaled merit function.
!
!     write(6,99) iter, wmerit, wmbest, wmprev
!     write(6,98) itwtch, alfold
!  99 format(1h ,'iter,wmer,wbest,wprev',i5,3(1pe15.5))
!  98 format(1h ,'itwtch,alfold',i5,1pe15.5)
!!!
!
      iflag = 0
!
      IF (wmerit <= wmbest) THEN
!
!        The current watchdog merit function is the best.
!  
         wmbest = wmerit
         itwtch = 0
         RETURN
      END IF
!
!     The current merit function is not the best.
!
      itwtch = itwtch + 1
!
!     Do not apply watchdog tests if (1) the previous step alfold
!     exceeds 1/2, or (2) the watchdog merit function decreased in
!     the immediately preceding iteration and itwtch does not
!     exceed twice its maximum.
!
      IF (alfold >= half) RETURN
      IF (wmerit <= wmprev .AND. itwtch <= 2*itwtmx) RETURN
!
!     If more than itwtmx iterations have occurred without
!     an overall improvement in the watchdog merit function,
!     signal for termination.
!
      IF (itwtch >= itwtmx) THEN
         iflag = -1
! 
!        If a too-large increase in the watchdog merit function
!        compared to the best value occurred, and iter .ge. itonew,
!        signal for termination.
! 
      ELSE IF (iter >= itonew .AND.  &
               wmerit > grfct*wmbest) THEN
         iflag = -1
      END IF
!
      RETURN
   END SUBROUTINE WTCHDG
!!!
!!!
   SUBROUTINE FNEVAL(ncomp, nmsh, xx, nudim, u, fval, rpar, ipar)
!
      USE BVPExtern, ONLY: FSUB
!
      IMPLICIT NONE
!
      INTEGER(i4), INTENT(IN)    :: ncomp
      INTEGER(i4), INTENT(INOUT) :: nmsh
      REAL(r8),    INTENT(INOUT) :: xx(nmsh)
      INTEGER(i4), INTENT(IN)    :: nudim
      REAL(r8),    INTENT(INOUT) :: u(nudim,nmsh)
      REAL(r8),    INTENT(INOUT) :: fval(ncomp,nmsh)
      REAL(r8),    INTENT(INOUT) :: rpar(*)
      INTEGER(i4), INTENT(INOUT) :: ipar(*)
!
      REAL(r8) :: hmsh
      INTEGER(i4) :: im
!
!     fneval evaluates the function values (from fsub) for
!     a given mesh xx and array u, and stores the values
!     in the array fval.
!
      CALL FSUB(ncomp, xx(1), u(1,1), fval(1,1), rpar, ipar)
!
      DO im = 1, nmsh-1
         hmsh = xx(im+1) - xx(im)
         CALL FSUB(ncomp, xx(im+1), u(1,im+1), fval(1,im+1), rpar, ipar)
      END DO
!
      RETURN
   END SUBROUTINE FNEVAL
!!!
!!!
   SUBROUTINE JACCAL(ncomp, nmsh, nlbc, xx, nudim, u, fval,  &
                     dgtm, dftm1, dftm2, uint, ajac, topblk, botblk, bhold, chold,  &
                     rpar, ipar)
!
      USE BVPShared, ONLY: pdebug, use_c, comp_c, uval0,  &
                           nminit, iprint, idum,  &
                           DCOPY, DDOT
      USE BVPExtern, ONLY: DFSUB, DGSUB
!
      IMPLICIT NONE
!
      INTEGER(i4), INTENT(IN)    :: ncomp
      INTEGER(i4), INTENT(INOUT) :: nmsh
      INTEGER(i4), INTENT(IN)    :: nlbc
      REAL(r8),    INTENT(INOUT) :: xx(nmsh)
      INTEGER(i4), INTENT(IN)    :: nudim
      REAL(r8),    INTENT(INOUT) :: u(nudim,nmsh)
      REAL(r8),    INTENT(INOUT) :: fval(ncomp,nmsh)
      REAL(r8),    INTENT(INOUT) :: dgtm(ncomp)
      REAL(r8),    INTENT(INOUT) :: dftm1(ncomp,ncomp)
      REAL(r8),    INTENT(INOUT) :: dftm2(ncomp,ncomp)
      REAL(r8),    INTENT(INOUT) :: uint(ncomp)
      REAL(r8),    INTENT(INOUT) :: ajac(ncomp,2*ncomp,nmsh-1)
      REAL(r8),    INTENT(INOUT) :: topblk(nlbc, ncomp)
      REAL(r8),    INTENT(INOUT) :: botblk(ncomp-nlbc,ncomp)
      REAL(r8),    INTENT(INOUT) :: bhold(ncomp, ncomp, nmsh-1)
      REAL(r8),    INTENT(INOUT) :: chold(ncomp, ncomp, nmsh-1)
      REAL(r8),    INTENT(INOUT) :: rpar(*)
      INTEGER(i4), INTENT(INOUT) :: ipar(*)
!
      REAL(r8), PARAMETER :: zero = 0.0D+0
      REAL(r8), PARAMETER :: half = 0.5D+0
      REAL(r8), PARAMETER :: eighth = 0.125D+0
      REAL(r8), PARAMETER :: four = 4.0D+0
      REAL(r8), PARAMETER :: six = 6.0D+0
      REAL(r8), PARAMETER :: one = 1.0D+0
      REAL(r8), PARAMETER :: three = 3.0D+0
      REAL(r8), PARAMETER :: twelve = 12.0D+0
!
      REAL(r8) :: hmsh, xhalf, dsq
      INTEGER(i4) :: ninter, i, im, ic, jc
!!!
!
!     if (pdebug) write(6,901)
!
      ninter = nmsh - 1
!
!     if (pdebug) write(6,902)
!
      DO i = 1, nlbc
         CALL DGSUB(i, ncomp, u(1,1), dgtm,rpar,ipar)
         CALL DCOPY(ncomp, dgtm(1), 1, topblk(i,1), nlbc)
!        if (pdebug) write(6,903) i, (topblk(i,j),j=1,ncomp)
      END DO
!
      CALL DFSUB(ncomp, xx(1), u(1,1), dftm1(1,1),rpar,ipar)
!
!     on entry to jaccal, the array fval contains the function values
!     at (xx(im), u(ic,im)), ic=1,...,ncomp and im = 1,...,nmsh,
!     calculated by a preceding call of rhscal with the same xx and u
!     arrays.
!
!     if (pdebug) write(6,904)
!
      DO im = 1, ninter
! 
         hmsh = xx(im+1) - xx(im)
! 
         DO ic = 1, ncomp
            uint(ic) = half*(u(ic,im) + u(ic,im+1))  &
               - eighth*hmsh*(fval(ic,im+1) - fval(ic,im))
         END DO
         xhalf = half*(xx(im+1) + xx(im))
         CALL DFSUB(ncomp, xhalf, uint, dftm2(1,1), rpar, ipar)
         DO ic = 1, ncomp
            DO jc = 1, ncomp
               dsq = DDOT(ncomp, dftm2(ic,1), ncomp, dftm1(1,jc), 1)
               ajac(ic,jc,im) = -hmsh*(dftm1(ic,jc)/six  &
                  + dftm2(ic,jc)/three + hmsh*dsq/twelve)
            END DO
            ajac(ic,ic,im) = ajac(ic,ic,im) - one
!            if (pdebug) write(6,905) im, ic,
!     *            (ajac(ic,jc,im), jc=1,ncomp)
         END DO
! 
         CALL DFSUB(ncomp, xx(im+1), u(1,im+1), dftm1(1,1), rpar, ipar)
         DO ic = 1, ncomp
            DO jc = 1, ncomp
               dsq = DDOT(ncomp, dftm2(ic,1), ncomp, dftm1(1,jc), 1)
               ajac(ic,jc+ncomp,im) = -hmsh*(dftm1(ic,jc)/six  &
                  + dftm2(ic,jc)/three - hmsh*dsq/twelve)
            END DO
            CALL DCOPY(ncomp, ajac(ic,ncomp+1,im), ncomp, chold(ic,1,im), ncomp)
            CALL DCOPY(ncomp, dftm1(ic,1), ncomp, bhold(ic,1,im), ncomp)
            ajac(ic,ic+ncomp,im) = ajac(ic,ic+ncomp,im) + one
            chold(ic,ic,im) = ajac(ic,ic+ncomp,im)
!            if (pdebug) write(6,905) im, ic,
!     *                    (ajac(ic,jc+ncomp,im),jc=1,ncomp)
         END DO
!
      END DO
!     if (pdebug) write(6,906)
!
      DO i = nlbc+1, ncomp
         CALL DGSUB(i, ncomp, u(1, nmsh), dgtm,rpar,ipar)
         CALL DCOPY(ncomp, dgtm(1), 1, botblk(i-nlbc,1), ncomp-nlbc)
!        if (pdebug) write(6,903) i,(botblk(i-nlbc,j), j=1,ncomp)
      END DO
!
!     write(6,991)
!991  format(1x,'topblk')
!     write(6,992) topblk(1,1),topblk(1,2)
!992  format(1x,2g22.10)
!     write(6,993)
!993  format(1x,'main jacobian')
!     do 994 iu=1,ninter
!     do 994 iv=1,ncomp
!     write(6,996) ajac(iv,1,iu),ajac(iv,2,iu),ajac(iv,3,iu),
!    +ajac(iv,4,iu)
!996  format(1x,4f12.7)
!994  continue
!     write(6,995)
!995  format(1x,'botblk')
!     write(6,992) botblk(1,1),botblk(1,2)
!
901   FORMAT(1H ,'jaccal')
902   FORMAT(1H ,'topblk')
903   FORMAT(1H ,i5,6(1PE11.3))
904   FORMAT(1H ,'ajac')
905   FORMAT(1H ,2I5,5(1PE11.3))
906   FORMAT(1H ,'botblk')
!
      RETURN
   END SUBROUTINE JACCAL
!!!
!!!
   SUBROUTINE LNRHS(ncomp, nmsh, nlbc, xx, nudim, u, &
                    rhs, rnsq, fval, ftmp, uint, rpar, ipar)
!
      USE BVPShared, ONLY: flmin, flmax, epsmch
      USE BVPExtern, ONLY: FSUB, GSUB
!
      IMPLICIT NONE
!
      INTEGER(i4), INTENT(IN)    :: ncomp
      INTEGER(i4), INTENT(INOUT) :: nmsh
      INTEGER(i4), INTENT(IN)    :: nlbc
      REAL(r8),    INTENT(INOUT) :: xx(*)
      INTEGER(i4), INTENT(IN)    :: nudim
      REAL(r8),    INTENT(INOUT) :: u(nudim,*)
      REAL(r8),    INTENT(INOUT) :: rhs(*)
      REAL(r8),    INTENT(INOUT) :: rnsq
      REAL(r8),    INTENT(INOUT) :: fval(ncomp,*)
      REAL(r8),    INTENT(INOUT) :: ftmp(*)
      REAL(r8),    INTENT(INOUT) :: uint(*)
      REAL(r8),    INTENT(INOUT) :: rpar(*)
      INTEGER(i4), INTENT(INOUT) :: ipar(*)
!
!     This subroutine is designed to calculate the right-hand
!     side for linear problems.
!
      REAL(r8), PARAMETER :: zero = 0.0D+0
      REAL(r8), PARAMETER :: half = 0.5D+0
      REAL(r8), PARAMETER :: eighth = 0.125D+0
      REAL(r8), PARAMETER :: one = 1.0D+0
      REAL(r8), PARAMETER :: four = 4.0D+0
      REAL(r8), PARAMETER :: six = 6.0D+0
!
      REAL(r8) :: wg, hmsh, xhalf, scale, sumsq
      INTEGER(i4) :: nrhs, loc
      INTEGER(i4) :: i, im, ic, ii, ninter
!!!
!
      ninter = nmsh - 1
      rnsq = zero
!
!     first, process the left-hand boundary conditions.
!
      DO i = 1, nlbc
         CALL GSUB(i, ncomp, u(1,1), wg,rpar,ipar)
         rhs(i) = -wg
      END DO
!
!     Next, process the interior mesh points.  The fval array
!     contains the function values from fsub at xx and u.
!
      DO im = 1, ninter
         hmsh = xx(im+1) - xx(im)
         DO ic = 1, ncomp
            uint(ic) = half*(u(ic,im) + u(ic,im+1))  &
               - eighth*hmsh*(fval(ic,im+1) - fval(ic,im))
         END DO
         xhalf = half*(xx(im) + xx(im+1))
         CALL FSUB(ncomp, xhalf, uint, ftmp,rpar,ipar)
         loc = (im-1)*ncomp + nlbc
         DO ic = 1, ncomp
            rhs(loc+ic) = -u(ic,im+1) + u(ic,im) + hmsh*  &
               (fval(ic,im) + fval(ic,im+1) + four*ftmp(ic))/six
         END DO
      END DO
!
      nrhs = ninter*ncomp
      DO ii = nlbc+1, ncomp
         CALL GSUB(ii, ncomp, u(1,nmsh), wg,rpar,ipar)
         rhs(nrhs+ii) = -wg
      END DO
!
      CALL DSSQ(nmsh*ncomp, rhs, 1, scale, sumsq)
      rnsq = (scale**2)*sumsq
!
      RETURN
   END SUBROUTINE LNRHS
!!!
!!!
   SUBROUTINE RHSCAL(ncomp, nmsh, nlbc, xx, nudim, u, defcor,  &
                     rhs, rnsq, fval, ftmp, uint, rpar, ipar)
!
      USE BVPShared, ONLY: pdebug, use_c, comp_c, uval0,  &
                           nminit, iprint, idum,  &
                           flmin, flmax, epsmch
      USE BVPExtern, ONLY: FSUB, GSUB
!
      IMPLICIT NONE
!
      INTEGER(i4), INTENT(IN)    :: ncomp
      INTEGER(i4), INTENT(INOUT) :: nmsh
      INTEGER(i4), INTENT(IN)    :: nlbc
      REAL(r8),    INTENT(INOUT) :: xx(nmsh)
      INTEGER(i4), INTENT(IN)    :: nudim
      REAL(r8),    INTENT(INOUT) :: u(nudim,nmsh)
      REAL(r8),    INTENT(INOUT) :: defcor(ncomp,nmsh-1)
      REAL(r8),    INTENT(INOUT) :: rhs(ncomp*nmsh)
      REAL(r8),    INTENT(INOUT) :: rnsq
      REAL(r8),    INTENT(INOUT) :: fval(ncomp,nmsh)
      REAL(r8),    INTENT(INOUT) :: ftmp(ncomp)
      REAL(r8),    INTENT(INOUT) :: uint(ncomp)
      REAL(r8),    INTENT(INOUT) :: rpar(*)
      INTEGER(i4), INTENT(INOUT) :: ipar(*)
!
!     This subroutine constructs the (ncomp*nmsh)-dimensional
!     vector rhs, which is the right-hand side of the Newton equations.
!     The ncomp by nmsh array fval is assumed to have been calculated
!     elsewhere by routine fneval.
!
      REAL(r8), PARAMETER :: zero = 0.0D+0
      REAL(r8), PARAMETER :: half = 0.5D+0
      REAL(r8), PARAMETER :: eighth = 0.125D+0
      REAL(r8), PARAMETER :: one = 1.0D+0
      REAL(r8), PARAMETER :: four = 4.0D+0
      REAL(r8), PARAMETER :: six = 6.0D+0
!
      REAL(r8) :: wg, scale, sumsq, hmsh, xhalf
      INTEGER(i4) :: ninter, im, ic, i, ii, nrhs, loc
!!!
!
!     if (pdebug) write(6,901)
!
!     ninter is the number of intervals in the mesh (one less than the
!     number of mesh points)
!
      ninter = nmsh - 1
      rnsq = zero
!
!     First, process the left-hand boundary conditions.
!   
      DO i = 1, nlbc
         CALL GSUB(i, ncomp, u(1,1), wg,rpar,ipar)
         rhs(i) = -wg
      END DO
!
!     Next, process the interior mesh points.  The fval array
!     contains the function values from fsub at xx and u.
!
      DO im = 1, ninter
         hmsh = xx(im+1) - xx(im)
         DO ic = 1, ncomp
            uint(ic) = half*(u(ic,im) + u(ic,im+1))  &
               - eighth*hmsh*(fval(ic,im+1) - fval(ic,im))
         END DO
         xhalf = half*(xx(im) + xx(im+1))
         CALL FSUB(ncomp, xhalf, uint, ftmp,rpar,ipar)
         loc = (im-1)*ncomp + nlbc
         DO ic = 1, ncomp
            rhs(loc+ic) = -u(ic,im+1) + u(ic,im) + defcor(ic,im) + hmsh*  &
               (fval(ic,im) + fval(ic,im+1) + four*ftmp(ic))/six
         END DO
      END DO
!
      nrhs = ninter*ncomp
      DO ii = nlbc+1, ncomp
         CALL GSUB(ii, ncomp, u(1,nmsh), wg,rpar,ipar)
         rhs(nrhs+ii) = -wg
      END DO
!
      CALL DSSQ(nmsh*ncomp, rhs, 1, scale, sumsq)
      rnsq = (scale**2)*sumsq
!
      IF (pdebug) THEN
         WRITE (6,902) rnsq
         WRITE(6,903)
         WRITE(6,904) (rhs(i), i=1,ncomp*nmsh)
      END IF
!
901   FORMAT(1H ,'rhscal')
902   FORMAT(1H ,'rnsq',1PE11.3)
903   FORMAT(1H ,'rhs vector')
904   FORMAT(1H ,(7(1PE11.3)))
!
      RETURN
   END SUBROUTINE RHSCAL
!!!
!!!
   SUBROUTINE DBLMSH(nmsh, nmax, xx, nmold, xxold, maxmsh)
!
      USE BVPShared, ONLY: pdebug, use_c, comp_c, uval0,  &
                           nminit, iprint, idum
!
      IMPLICIT NONE
!
      INTEGER(i4), INTENT(INOUT) :: nmsh
      INTEGER(i4), INTENT(INOUT) :: nmax
      REAL(r8),    INTENT(INOUT) :: xx(*)
      INTEGER(i4), INTENT(INOUT) :: nmold
      REAL(r8),    INTENT(INOUT) :: xxold(*)
      LOGICAL,     INTENT(INOUT) :: maxmsh
!
      REAL(r8), PARAMETER :: half = 0.5D+0
      INTEGER(i4) :: ninnew, nmnew, id2, i
!
!     This routine is used to double the mesh, i.e., produce a mesh
!     with twice as many intervals in which each new interval is
!     half the corresponding old interval.
!
!     On entry to dblmsh, the integer nmsh and the array xx
!     specify a set of mesh points xx(1),..., xx(nmsh) (assumed
!     to be in ascending order).
!
!     If the number of mesh points in the doubled mesh would
!     exceed the maximum allowed number nmax, the flag maxmsh is
!     set to true, and we exit without changing any other parameters.
!
!     Otherwise, nmold is set to the old number of mesh points,
!     xxold is set to the old mesh, nmsh is the new number of mesh
!     points, and xx contains the new mesh points.
!
      nmold = nmsh
      CALL DCOPY(nmold, xx, 1, xxold, 1)
!
      ninnew = 2*(nmsh-1)
      nmnew = ninnew + 1
      IF (nmnew >= nmax) THEN
         IF (iprint >= 0) WRITE(6,901) nmnew
         nmsh = nmold
         CALL DCOPY(nmold, xxold, 1, xx, 1)
         maxmsh = .true.
         RETURN
      END IF
      maxmsh = .false.
!
!     Loop backwards through the old mesh points to create the new ones.
!
      xx(nmnew) = xx(nmsh)
      DO i = ninnew, 4, -2
         id2 = i/2
         xx(i) = half*(xx(i+1) + xx(id2))
         xx(i-1) = xx(id2)
      END DO
!
!     Calculate the new xx(2). xx(1) remains unchanged.
!
      xx(2) = half*(xx(3) + xx(1))
      nmsh = nmnew
      IF (iprint >= 0) WRITE(6,902) nmsh
!
901   FORMAT (1H , ' dblmsh.  maximum mesh exceeded, nmnew =', i8)
902   FORMAT (1H , ' dblmsh.  the doubled mesh has ', i8,' points.')
!
      RETURN
   END SUBROUTINE DBLMSH
!!!
!!!
   SUBROUTINE SELMSH(ncomp, nmsh, ntol, ltol, tol, nfxpnt, fixpnt, ipow, nmax,  &
                     xx, nudim, u, ermeas, irefin, ihcomp, nmold, xxold, ermx, ddouble , maxmsh)
!
      USE BVPShared, ONLY: pdebug, use_c, comp_c, uval0,  &
                           nminit, iprint, idum,  &
                           DCOPY
!
      IMPLICIT NONE
!
      INTEGER(i4), INTENT(IN)    :: ncomp
      INTEGER(i4), INTENT(INOUT) :: nmsh
      INTEGER(i4), INTENT(IN)    :: ntol
      INTEGER(i4), INTENT(IN)    :: ltol(ntol)
      REAL(r8),    INTENT(IN)    :: tol(ntol)
      INTEGER(i4), INTENT(IN)    :: nfxpnt
      REAL(r8),    INTENT(IN)    :: fixpnt(*)
      INTEGER(i4), INTENT(INOUT) :: ipow
      INTEGER(i4), INTENT(INOUT) :: nmax
      REAL(r8),    INTENT(INOUT) :: xx(*)
      INTEGER(i4), INTENT(IN)    :: nudim
      REAL(r8),    INTENT(INOUT) :: u(nudim,*)
      REAL(r8),    INTENT(INOUT) :: ermeas(ncomp,*)
      INTEGER(i4), INTENT(INOUT) :: irefin(nmsh-1)
      INTEGER(i4), INTENT(INOUT) :: ihcomp(nmsh-1)
      INTEGER(i4), INTENT(INOUT) :: nmold
      REAL(r8),    INTENT(INOUT) :: xxold(*)
      REAL(r8),    INTENT(INOUT) :: ermx(*)
      LOGICAL,     INTENT(INOUT) :: ddouble
      LOGICAL,     INTENT(INOUT) :: maxmsh
!
      REAL(r8), PARAMETER :: zero = 0.0D+0
      REAL(r8), PARAMETER :: one = 1.0D+0
      REAL(r8), PARAMETER :: onep1 = 1.1D+0
      REAL(r8), PARAMETER :: erdcid = 5.0D+0
      REAL(r8), PARAMETER :: phitst = 0.1D+0
!
      LOGICAL, SAVE :: first = .true.
      REAL(r8), SAVE :: rlndec
      REAL(r8) :: frcpow, thres, errmax, denom, ems, err, decii, rlen, slen, dx
      REAL(r8) :: fxnext, rold, rlold, phihat, val1
      INTEGER(i4) :: ninter, ithres, im, it, jcomp, ii, ilg, nmest, new
      INTEGER(i4) :: j, ifxcnt, jtkout, ind1, nmnew
!!!
!
!     The routine selmsh performs selective mesh refinement, depending
!     on the error measure ermeas.
!
      IF (first) THEN
         first = .false.
         rlndec = DLOG(erdcid)
      END IF
!
      maxmsh = .false.
!
      IF (pdebug) WRITE(6,901) nmsh, ipow
!
      frcpow = one/ipow
      ddouble = .false.
      nmold = nmsh
      ninter = nmsh - 1
!
!     Copy the current mesh into the xxold array.
!
      CALL DCOPY(nmold, xx, 1, xxold, 1)
!
      ithres = 0
      thres = one
!
!     On input, the array ermeas represents some error measure defined
!     over the components and mesh intervals (not mesh points).
!     It is normalized in the following loop with respect to the
!     tolerance array and the current solution.
!     The value errmax gives the maximum normalized error.
!
      errmax = zero
      DO im = 1, ninter
         ermx(im) = zero
         DO it = 1, ntol
            jcomp = ltol(it)
            denom = tol(it)*MAX(one, ABS(u(jcomp,im)))
            ems = ermeas(jcomp,im)
            ermeas(jcomp,im) = ABS(ems)/denom
!            if (pdebug .and. ermeas(jcomp,im) .ge. thres)
!     *             write(6,902) im,jcomp,ems,ermeas(jcomp,im)
            err = ermeas(jcomp, im)
            IF (err >= ermx(im)) THEN
               ermx(im) = err
               ihcomp(im) = jcomp
            END IF
         END DO
         errmax = MAX(ermx(im), errmax)
      END DO
!
      IF (pdebug) WRITE(6,903) errmax
!
      IF (errmax > zero .AND. errmax <= erdcid) THEN
! 
!        If errmax > 0 and .le. erdcid, find the smallest integer exponent ii
!        such that (erdcid**ii)*errmax > erdcid.
! 
         IF (errmax > one) THEN
            ii = 1
            decii = erdcid
         ELSE
            ilg = -DLOG(errmax)/rlndec
            ii = 2 + ilg
            decii = erdcid**ii
         END IF
!
!        Multiply error measures by erdcid**ii.
!
         errmax = decii*errmax
         DO im = 1, ninter
            ermx(im) = decii*ermx(im)
            DO it = 1, ntol
               jcomp = ltol(it)
               ermeas(jcomp,im) = decii*ermeas(jcomp,im)
            END DO
         END DO
!
      END IF
!
200   CONTINUE
!
!     For each interval im,  the integer irefin(im) is calculated
!     based on an equidistrbution procedure involving the
!     threshold value thres.  If irefin(im) > 1, we add
!     points to interval im.  If irefin(im) = 1, we check whether
!     point im can be removed from the mesh.
!
!     nmest is a lower bound on the number of points in the new mesh.
!     We do not know in advance how many points will be removed,
!     so nmest is computed by assuming all eligible points are removed.
!
      nmest = nmsh
      DO im = 1, ninter
         IF (ermx(im) >= thres) THEN
            irefin(im) = INT(ermx(im)**frcpow) + 1
            nmest = nmest + irefin(im) - 1
         ELSE
            irefin(im) = 1
            nmest = nmest - 1
         END IF
      END DO
!     if (pdebug) write(6,904) nmest, (irefin(i), i=1,ninter)
!
      IF (nmest > nmax) THEN
! 
         GO TO 360
!  
      ELSE IF (nmest-1 > 3*ninter) THEN
! 
         CALL DBLMSH(nmsh, nmax, xx, nmold, xxold, maxmsh)
         ddouble = .true.
         RETURN
!
      END IF
!
!     It appears that we can perform the desired selective mesh
!     refinement.
!
!     Now begin running through the mesh, adding and possibly deleting
!     points as indicated by the irefin array.
!
!     The integer new is a count of the number of intervals in
!     the tentative mesh being generated by the refinement strategy.
!
      new = 1
!
!     The first interval is treated as a special case, since xx(1)
!     always remains in the mesh, and cannot be a fixed point.
!
      rlen = xxold(2) - xx(1)
      slen = rlen
      IF (irefin(1) > 1) THEN
         dx = rlen/irefin(1)
         DO j = 2, irefin(1)
            new = new + 1
            xx(new) = xx(1) + (j-1)*dx
         END DO
      END IF
!
!     The fixed points specified by the fixpnt array cannot be
!     removed from the mesh.  The value fxnext indicates the 'next'
!     fixed point. When no further fixed points remain to be processed
!     (or if nfxpnt = 0), fxnext is set to a value strictly larger than
!     the last mesh point, so that no mesh point can equal fxnext.
!     This way we need to compare only the values of xxold(i)
!     and fxnext.
!
      ifxcnt = 1
      IF (nfxpnt == 0) THEN
         fxnext = onep1*ABS(xxold(nmsh))
      ELSE
         fxnext = fixpnt(ifxcnt)
      END IF
!
!     jtkout is a counter of the number of consecutive points that
!     have been removed from the mesh.
!
      jtkout = 0
!
      DO im = 2, ninter
!
         rlold = rlen
         rlen = xxold(im+1) - xxold(im)
!
!        If xxold(im) is the next fixed point, it cannot be removed
!        and so we don't test its error estimates.
! 
         IF (xxold(im) == fxnext) THEN
            ifxcnt = ifxcnt + 1
            IF (ifxcnt > nfxpnt) THEN
               fxnext = onep1*ABS(xxold(nmsh))
            ELSE
               fxnext = fixpnt(ifxcnt)
            END IF
         ELSE IF (irefin(im) == 1) THEN
!
!           If xxold(im) is not a fixed point and irefin(im) = 1, we may wish
!           to remove point im from the mesh.
! 
!           If we are considering removing points and jtkout = 0, this
!           is the first point in a possible consecutive set to be removed,
!           and we initialize phihat, which represents a maximum of
!           certain estimates.
!           If jtkout is not zero, previous points contiguous to this
!           point have been removed, and phihat retains its previous value.
! 
            slen = slen + rlen
! 
            IF (jtkout == 0) THEN
               ind1 = ihcomp(im-1)
               phihat = ermeas(ind1,im-1)/(rlold**ipow)
            END IF
            phihat = MAX(phihat, ermeas(ihcomp(im),im)/(rlen**ipow))
            val1 = phihat*(slen**ipow)
            IF (val1 <= phitst .AND. jtkout < 4) THEN
! 
!              Increment the counter of removed points.
!              'Remove' the mesh point xxold(im) by not including it.
! 
               jtkout = jtkout+1
               CYCLE
            END IF
!        end of logic for irefin(im) = 1.
         END IF
!
         jtkout = 0
         new = new + 1
         xx(new) = xxold(im)
! 
         IF (new + irefin(im)  > nmax) GO TO 360
! 
         IF (irefin(im) > 1) THEN
            dx = rlen/irefin(im)
            DO j = 2, irefin(im)
               new = new + 1
               xx(new) = xxold(im) + (j-1)*dx
            END DO
         END IF
         slen = rlen
! 
         IF (new  > nmax) THEN
!    
!           If the new mesh contains too many points, branch out of the
!           loop to try alternative strategies.
!    
            GO TO 360
!
         ELSE IF (new > 3*ninter) THEN
! 
!           Here, the new mesh does not exceed the specified maximum,
!           but has more than 3 times as many intervals as the old mesh.
!           Try doubling the mesh if possible.
! 
            CALL DCOPY(nmsh, xxold, 1, xx, 1)
            CALL DBLMSH(nmsh, nmax, xx, nmold, xxold, maxmsh)
            ddouble = .true.
            RETURN
! 
         END IF
!
      END DO
!
!     To end up here, we have processed the entire interval,
!     and have neither exceeded the specified maximum nor
!     exceeded three times the number of intervals in the old
!     mesh. The last mesh point remains unchanged.
!
      new = new + 1
      xx(new) = xxold(nmsh)
      nmsh = new
      maxmsh = .false.
      IF (iprint >= 0) WRITE(6,905) nmsh
      RETURN
!
360   CONTINUE
!
!     To reach here, the number of mesh points created at some stage
!     of the refinement process was larger than the maximum permitted
!     value nmax.
!
!     Check whether the mesh can safely be doubled.
!
      IF ((2*nmsh-1) < nmax) THEN
! 
!        Double the mesh.
         CALL DBLMSH(nmsh, nmax, xx, nmold, xxold, maxmsh)
         ddouble = .true.
! 
!        If the number of intervals is too large and the mesh cannot be
!        doubled, increase the threshold thres by a factor of erdcid and
!        try the selective refinement again.
!        If this happens three times without success or if thres exceeds
!        or is equal to errmax, stop.  (In this case, we know already
!        that doubling the mesh produces too many points.)
! 
      ELSE IF (thres < errmax .AND. ithres < 3) THEN
         ithres = ithres + 1
         thres = erdcid*thres
         IF(thres > errmax) thres = errmax
         CALL DCOPY(nmsh, xxold, 1, xx, 1)
         GO TO 200
      ELSE
!        nmsh = 2*nmsh - 1
         nmsh = nmold
         CALL DCOPY(nmold, xxold, 1, xx, 1)
         maxmsh = .true.
      END IF
!
901   FORMAT(1H ,'selmsh.  nmsh, ipow =',2I5)
902   FORMAT(1H ,'im, jcomp, ermeas, normalized er',2I5,2(1PE11.3))
903   FORMAT(1H ,'errmax',1PE11.3)
904   FORMAT(1H ,'nmest, irefin',(10I5))
905   FORMAT(1H ,'selmsh.  new nmsh =',i8)
910   FORMAT(1H ,'ihcomp',(10I5))
!
      RETURN
   END SUBROUTINE selmsh
!!!
!!!
   SUBROUTINE SMPMSH(nmsh, nmax, xx, intref, numadd, nmold, xxold, maxmsh)
!
      USE BVPShared, ONLY: pdebug, use_c, comp_c, uval0,  &
                           nminit, iprint, idum
!
      IMPLICIT NONE
!
      INTEGER(i4), INTENT(INOUT) :: nmsh
      INTEGER(i4), INTENT(INOUT) :: nmax
      REAL(r8),    INTENT(INOUT) :: xx(*)
      INTEGER(i4), INTENT(INOUT) :: intref
      INTEGER(i4), INTENT(INOUT) :: numadd
      INTEGER(i4), INTENT(INOUT) :: nmold
      REAL(r8),    INTENT(INOUT) :: xxold(*)
      LOGICAL,     INTENT(INOUT) :: maxmsh
!
      REAL(r8) :: dx
      INTEGER(i4) :: nmnew, nint, noalt, nochsm
      INTEGER(i4) :: i, j, innew, ninter
!
!     The routine smpmsh performs simple mesh refinement by adding
!     points to one or three interval(s) in the region indicated
!     by the integer intref.
!     numadd gives the trial number of points to be added in each
!     interval.
!!!
!
      IF (pdebug) WRITE(6,901)
!
      nmold = nmsh
      CALL DCOPY(nmold, xx, 1, xxold, 1)
!
!     numadd is altered if necessary so that it lies between 4 and 49
!
      IF (numadd > 49) THEN
         numadd = 49
      ELSE IF (numadd < 4) THEN
         numadd = 4
      END IF
      IF (pdebug) WRITE (6,902) nmsh, intref, numadd
!
      maxmsh = .false.
!
      IF (intref == 1) THEN
!
!        Add numadd points to the first interval if intref = 1.
!
         nmnew = nmsh + numadd
         IF (nmnew > nmax) THEN
            IF (iprint >= 0)  WRITE(6,903) nmnew
            maxmsh = .true.
            RETURN
         END IF
! 
!        Renumber the later points in reverse order.
! 
         nint = numadd + 1
         DO i = nmnew, numadd+2, -1
            xx(i) = xx(i-numadd)
         END DO
         dx = (xx(2) - xx(1))/nint
         DO i = 2, nint
            xx(i) = xx(1) + (i-1)*dx
         END DO
!
      ELSE IF (intref == nmsh-1) THEN
!  
!        Add numadd points to the last interval if intref = nmsh-1.
!
         nmnew = nmsh + numadd
         IF (nmnew > nmax) THEN
            IF (iprint >= 0)  WRITE(6,903) nmnew
            maxmsh = .true.
            RETURN
         END IF
         nint = numadd + 1
         dx = (xx(nmsh) - xx(nmsh-1))/nint
         xx(nmnew) = xx(nmsh)
         DO i = nmsh, nmnew-1
            xx(i) = xx(nmsh-1) + (i-nmsh+1)*dx
         END DO
! 
      ELSE
!
         IF (numadd > 9) numadd = 9
!
!        Here, intref lies between 2 and nmsh-2.  Add numadd points to
!        each of the three intervals intref-1, intref and intref+1.
! 
         nmnew = nmsh + 3*numadd
         IF (nmnew > nmax) THEN
            IF (iprint >= 0) WRITE(6,903) nmnew
            maxmsh = .true.
            RETURN
         END IF
!
!        noalt is the number of points at the right end of the interval
!        whose numerical values remain the same, but whose indices change.
!        nochsm is the smallest index in the new ordering of one of these
!        points.
! 
         noalt = nmsh - intref - 1
         nochsm = nmnew - noalt + 1
! 
!        Renumber the noalt unchanged points at the right end of the interval
!        (in reverse order).
! 
         j = 0
         DO i = nmnew, nochsm, -1
            xx(i) = xx(nmsh-j)
            j = j + 1
         END DO
! 
!        Add numadd points to the three contiguous intervals.
!        The remaining points at the left end of the interval retain
!        their original indices, and are left unchanged.
! 
         nint = numadd + 1
         innew = nochsm - nint
         DO i = intref+1, intref-1, -1
            xx(innew) = xx(i)
            dx = (xx(innew + nint) - xx(innew))/nint
            DO j = 1, numadd
               xx(innew + j) = xx(innew) + j*dx
            END DO
            innew = innew - nint
         END DO
!  
      END IF
!
      nmsh = nmnew
!
      IF(iprint >= 0)  WRITE(6,904) nmsh
!
901   FORMAT(1H , ' smpmsh')
902   FORMAT(1H ,'nmsh, intref, numadd',3I6)
903   FORMAT(1H , ' smpmsh.  maximum points exceeded, nmnew =',i6)
904   FORMAT(1H ,'smpmsh, new nmsh =',i7)
!
      RETURN
   END SUBROUTINE SMPMSH
!!!
!!!
   SUBROUTINE UNIMSH(nmsh, aleft, aright, nfxpnt, fixpnt, xx)
!
      USE BVPShared, ONLY: pdebug, use_c, comp_c, uval0,  &
                           nminit, iprint, idum
!
      IMPLICIT NONE
!
      INTEGER(i4), INTENT(INOUT) :: nmsh
      REAL(r8),    INTENT(IN)    :: aleft
      REAL(r8),    INTENT(IN)    :: aright
      INTEGER(i4), INTENT(IN)    :: nfxpnt
      REAL(r8),    INTENT(IN)    :: fixpnt(*)
      REAL(r8),    INTENT(INOUT) :: xx(nmsh)
!
      REAL(r8) :: dx, xleft, totint, xright
      INTEGER(i4) :: i, j, innew, ninter, ileft, ndif, nmin, iright, npt
!
!     Given a left endpoint aleft, a right endpoint aright,
!     a set of nfxpnt fixed points fixpnt(i), i = 1,...,nfxpnt,
!     (where fixpnt(i) is different from aleft and aright for all i),
!     and an initial target number nmsh of mesh points,
!     the subroutine unimsh generates a piecewise uniform mesh
!     beginning at aleft, ending at aright, and with equally
!     spaced points between aleft and fixpnt(1), then between
!     fixpnt(1) and fixpnt(2), ..., and finally between
!     fixpnt(nfxpnt) and aright.  The final number of intervals
!     is the maximum of nfxpnt+2 and the initial value of nmsh.
!
!     In the simplest case when nfxpnt = 0, unimsh generates a
!     uniform mesh with nmsh intervals in the closed interval
!     (aleft, aright).
!
!     On exit, the integer nmsh contains the number of mesh points
!     (which is the maximum of the initial nmsh and nfxpnt).
!     The array xx (of dimension nmsh) contains the mesh points.
!!!
!
      IF (iprint >= 0) WRITE(6,901) nmsh
!
      IF (nfxpnt == 0) THEN
! 
!        If there are no interior fixed points, the spacing is uniform
!        throughout the interval.  Calculate the spacing dx
!        and set up the xx array.
! 
         ninter = nmsh - 1
!
         dx = (aright - aleft)/ninter
         DO i = 1, ninter
            xx(i) = aleft + (i-1)*dx
         END DO
         xx(nmsh) = aright
         RETURN
      END IF
!
!     We know that there is at least one fixed point strictly between
!     the endpoints.
!
      IF (nmsh < nfxpnt+2)  nmsh = nfxpnt + 2
      ninter = nmsh - 1
      xx(1) = aleft
      ileft = 1
      xleft = aleft
      totint = aright - aleft
      ndif = ninter - nfxpnt
!
      DO j = 1, nfxpnt + 1
!
!        Deal in turn with the subintervals defined by the interval
!        boundaries and the fixed  points.
! 
         IF (j < nfxpnt+1) THEN
! 
!           The j-th fixed point is xright.  Calculate where it should
!           fall in the mesh.
! 
            xright = fixpnt(j)
            nmin = ninter*(xright-aleft)/totint + 1.5D+0
            IF (nmin > ndif+j) nmin = ndif + j
            iright = MAX(ileft+1, nmin)
         ELSE
            xright = aright
            iright = nmsh
         END IF
! 
!        npt is the number of equally spaced points that should
!        lie strictly between the (j-1)-th and j-th fixed points.
! 
         xx(iright) = xright
         npt = iright - ileft - 1
         dx = (xright - xleft)/(npt + 1)
         DO i = 1, npt
            xx(ileft+i) = xleft + i*dx
         END DO
         ileft = iright
         xleft = xright
!
      END DO
!
901   FORMAT (1H ,'unimsh.  nmsh =',i5)
!
      RETURN
   END SUBROUTINE UNIMSH
!!!
!!!
   SUBROUTINE STATS(len, elem, ebigst, esecnd, summod, index)
!
      IMPLICIT NONE
!
      INTEGER(i4), INTENT(IN)    :: len
      REAL(r8),    INTENT(INOUT) :: elem(len)
      REAL(r8),    INTENT(INOUT) :: ebigst
      REAL(r8),    INTENT(INOUT) :: esecnd
      REAL(r8),    INTENT(INOUT) :: summod
      INTEGER(i4), INTENT(INOUT) :: index
!
      REAL(r8), PARAMETER :: zero = 0.0D+0
      REAL(r8) :: elmod
      INTEGER(i4) :: i
!
!     Given the real array elem of length len, stats calculates
!     the following:
!      - summod, the sum of the magnitudes of the elements of elem;
!      - ebigst (the largest element in magnitude);
!      - index (the index in elem of ebigst); and
!      - esecnd (the second largest element in magnitude, strictly
!          less than ebigst unless both are zero).
!!!
!
      index = 1
      ebigst = zero
      esecnd = zero
      summod = zero
!
      DO i = 1, len
         elmod = ABS(elem(i))
         summod = summod + elmod
         IF (elmod > ebigst) THEN
            esecnd = ebigst
            ebigst = elmod
            index = i
         ELSE IF (elmod > esecnd) THEN
            esecnd = elmod
         END IF
      END DO
!
      RETURN
   END SUBROUTINE STATS
!!!
!!!
   SUBROUTINE MSHREF(ncomp, nmsh, nlbc, ntol, ltol, iorder, rhs, tmwork,  &
                     nmax, xx, nmold, xxold, ddouble , maxmsh,  &
                     numbig, nummed, amg, stab_cond, stiff_cond, r4,  &
                     nfxpnt, fixpnt, irefin, itcond, itcondmax)
!
      USE BVPShared, ONLY: pdebug, use_c, comp_c, uval0,  &
                           nminit, iprint, idum
!
      IMPLICIT NONE
!
      INTEGER(i4), INTENT(IN)    :: ncomp
      INTEGER(i4), INTENT(INOUT) :: nmsh
      INTEGER(i4), INTENT(IN)    :: nlbc
      INTEGER(i4), INTENT(IN)    :: ntol
      INTEGER(i4), INTENT(IN)    :: ltol(ntol)
      INTEGER(i4), INTENT(INOUT) :: iorder
      REAL(r8),    INTENT(INOUT) :: rhs(ncomp*nmsh)
      REAL(r8),    INTENT(INOUT) :: tmwork(nmsh-1)
      INTEGER(i4), INTENT(INOUT) :: nmax
      REAL(r8),    INTENT(INOUT) :: xx(nmsh)
      INTEGER(i4), INTENT(INOUT) :: nmold
      REAL(r8),    INTENT(INOUT) :: xxold(nmold)
      LOGICAL,     INTENT(INOUT) :: ddouble
      LOGICAL,     INTENT(INOUT) :: maxmsh
      INTEGER(i4), INTENT(INOUT) :: numbig
      INTEGER(i4), INTENT(INOUT) :: nummed
      REAL(r8),    INTENT(INOUT) :: amg(*)
      LOGICAL,     INTENT(INOUT) :: stab_cond
      LOGICAL,     INTENT(INOUT) :: stiff_cond
      REAL(r8),    INTENT(INOUT) :: r4(*)
      INTEGER(i4), INTENT(IN)    :: nfxpnt
      REAL(r8),    INTENT(IN)    :: fixpnt(*)
      INTEGER(i4), INTENT(INOUT) :: irefin(*)
      INTEGER(i4), INTENT(INOUT) :: itcond
      INTEGER(i4), INTENT(IN)    :: itcondmax
!
      REAL(r8),    PARAMETER :: two = 2.0D+0
      REAL(r8),    PARAMETER :: bigfac = 10.0D+0
      REAL(r8),    PARAMETER :: small = 1.0D-2
      INTEGER(i4), PARAMETER :: numpt = 14
!
      REAL(r8) :: rbigst, rsecnd, sumrhs, tstval, er, denom, errel, tol, deltaf
      REAL(r8) :: ctry, cbest, atrue, btrue, alfaw, gap
      LOGICAL :: forcedouble, nodouble
      INTEGER(i4) :: ninter, nup, ic, icmp, indrhs, intref, numadd
      INTEGER(i4) :: im, it
!
!     This routine performs calculations leading to a decision
!     about how the mesh will be refined, and then refines the mesh.
!
!     The choices for mesh refinement in this routine are either to
!     double the mesh (i.e., divide each existing interval in half),
!     or to add points to just a few intervals.
!
!     The decision is made based on two criteria:  (1) the distribution
!     of magnitudes of components of rhs (broadly speaking, if
!     the maximum component of rhs is much larger than the
!     average, points are added near the corresponding interval),
!     and (2) the history of previous mesh refinements (if points
!     have been added to only a few intervals already, this strategy
!     is abandoned for the moment and the mesh is doubled).
!
!     The decision is indicated by setting the logical flag double
!     and (if ddouble is .false.) the integer intref.
!     Setting ddouble to .true. means that the mesh should be doubled
!     (i.e., the new mesh should contain twice as many intervals).
!     The integer intref indicates the region of the mesh where the
!     points should be added (see the routine smpmsh), and numadd
!     indicates how many points are to be added.
!     The integers nummed and numbig represent running totals,
!     used in deciding on the mesh refinement strategy.
!
!     If iorder = 4, meaning that we were just performing a Newton
!     iteration for a 4th order solution, check all elements of rhs.
!     If iorder .gt. 4, signalling that we were trying Newton
!     iterations for order 6 or 8, check only elements of rhs
!     corresponding to components for which a tolerance is specified.
!
      IF (pdebug) WRITE(6,901) nummed, numbig
!
      nodouble = ((iorder == 4) .AND.  &
         (stiff_cond .AND. .NOT. stab_cond) .AND. (use_c))
!
!     nodouble = nodouble
!    *  .or.((iorder.gt.4) .and. ( stiff_cond .and. .not. stab_cond)
!    *   .and. (use_c))
!
      forcedouble = .false.
      IF (use_c) THEN
         IF (itcond == itcondmax) THEN
            itcond = 0
            forcedouble = .true.
         ELSE
            forcedouble = .false.
         END IF
      END IF
!
      ninter = nmsh-1
      nup = ncomp
      IF (iorder > 4) nup = ntol
!
!     Check the vector rhs for a non-negligible component
!     whose magnitude is significantly larger than the average.
!     (small defines negligible, and bigfac defines significantly larger.)
!
      DO ic = 1, nup
!
         icmp = ic
         IF (iorder > 4) icmp = ltol(ic)
! 
!        For component icmp, examine the ninter elements of rhs not
!        corresponding to boundary conditions.
!
!        subroutine stats calculates rbigst and rsecnd (the first- and
!        second-largest elements of rhs in magnitude), and the index
!        intref of the interval in which the largest value occurs.
!        The value sumrhs is the sum of the magnitudes of the components
!        of rhs.
!
         indrhs = nlbc + ic
!
!        Copy the elements of rhs corresponding to interior mesh
!        points for component icmp into a single vector tmwork.
! 
         CALL DCOPY(ninter, rhs(indrhs), ncomp, tmwork, 1)
         CALL STATS(ninter, tmwork, rbigst, rsecnd, sumrhs, intref)
!
         tstval = bigfac*(sumrhs-rbigst)/ninter
         IF (pdebug) WRITE(6,902) ic, tstval, rbigst, rsecnd
         IF (rbigst >= small .AND. rbigst >= tstval) GO TO 100
!
      END DO
!
!     If we reach this point, no interval has a significantly larger
!     element than average.  Set counters and double the mesh.
!
      numbig = 0
      nummed = 0
      ddouble = .true.
      IF (pdebug) WRITE(6,903)
!     f the mesh if not doubled if the problem is stiff and the order is 4
      IF (nodouble .AND. .NOT. forcedouble) THEN
         CALL SELCONDMSH(ncomp, nmsh, nfxpnt, fixpnt,  nmax, xx,  irefin,  &
                         nmold, xxold, ddouble , maxmsh,r4,amg)
         ddouble = .false.
         itcond = itcond + 1
      ELSE
         CALL DBLMSH(nmsh, nmax, xx, nmold, xxold, maxmsh)
         itcond = 0
      END IF
!     call dblmsh(nmsh, nmax, xx, nmold, xxold, maxmsh)
      RETURN
!
!     To reach statement 100, it must be true that for some component
!     icmp,  rbigst is non-negligible and large relative to the
!     average.  intref indicates the region to which points may be added.
!
!     If too many specialized refinements (adding a few points to
!     a small number of intervals) have been made, signal that
!     the mesh should be doubled.
!     Otherwise, increment counters and add numadd points as indicated.
!
100   CONTINUE
      IF (rbigst < two*rsecnd) nummed = nummed + 1
      numadd = numpt
      numbig = numbig + 1
      ddouble = .false.
!
      IF (rbigst <= bigfac*rsecnd .OR. numbig > 8) THEN
         numbig = 0
         nummed = nummed + 1
         IF (nummed >= 4 .AND. iorder == 4) THEN
            ddouble = .true.
            nummed = 0
         ELSE IF (nummed >= 8 .AND. iorder > 4) THEN
            ddouble = .true.
            nummed = 0
         END IF
      END IF
!
!     Refine the mesh.
!
      IF (pdebug) WRITE(6,904) numbig, nummed
!
      IF (ddouble) THEN
!
!f       the mesh if not doubled if the problem is stiff
!
         IF  (nodouble .AND. .NOT. forcedouble)  THEN
            numadd = numpt
            CALL SMPSELCONDMSH(ncomp, nmsh, nfxpnt, fixpnt, nmax, xx, irefin, intref, numadd,  &
                               nmold, xxold, ddouble, maxmsh, r4, amg)
!               itcond = itcond + 1
!           call selcondmsh(ncomp, nmsh,
!     *        nfxpnt, fixpnt,  nmax, xx,  irefin,
!     *        nmold, xxold, ddouble , maxmsh,r4,amg)
            ddouble = .false.
            itcond = itcond + 1
         ELSE
            CALL DBLMSH(nmsh, nmax, xx, nmold, xxold, maxmsh)
            itcond  = 0
         END IF
!
      ELSE
!
!f       if the problem is stiff  we use both the old technique
!f       and the conditioning
!f
         IF (nodouble .AND. .NOT. forcedouble)  THEN
            numadd = numpt
            CALL SMPSELCONDMSH(ncomp, nmsh, nfxpnt, fixpnt, nmax, xx, irefin, intref, numadd,  &
                               nmold, xxold, ddouble, maxmsh, r4, amg)
            itcond = itcond + 1
         ELSE IF (forcedouble .AND. use_c) THEN
            ddouble = .true.
            itcond = 0
            CALL DBLMSH(nmsh, nmax, xx, nmold, xxold, maxmsh)
         ELSE
            numadd = numpt
            CALL SMPMSH(nmsh, nmax, xx, intref, numadd, nmold, xxold, maxmsh)
         END IF
!
      END IF
!
      IF (ddouble  .AND. pdebug) WRITE(6,905)
!
901   FORMAT(1H ,'mshref. nummed, numbig =',2I5)
902   FORMAT(1H ,'ic, tst, bigst, second',i5, 3(1PE11.3))
903   FORMAT(1H ,'No significantly large value')
904   FORMAT(1H ,'numbig, nummed =',2I5)
905   FORMAT(1H ,'double the mesh')
!
      RETURN
   END SUBROUTINE MSHREF
!!!
!!!
   SUBROUTINE ERREST(ncomp, nmsh, ntol, ltol, tol,  &
                     nudim, u, uold, etest, errsum, errok)
!
      USE BVPShared, ONLY: pdebug, use_c, comp_c, uval0,  &
                           nminit, iprint, idum
!
      IMPLICIT NONE
!
      INTEGER(i4), INTENT(IN)    :: ncomp
      INTEGER(i4), INTENT(INOUT) :: nmsh
      INTEGER(i4), INTENT(IN)    :: ntol
      INTEGER(i4), INTENT(IN)    :: ltol(ntol)
      REAL(r8),    INTENT(IN)    :: tol(ntol)
      INTEGER(i4), INTENT(IN)    :: nudim
      REAL(r8),    INTENT(INOUT) :: u(nudim,nmsh)
      REAL(r8),    INTENT(INOUT) :: uold(ncomp,nmsh)
      REAL(r8),    INTENT(INOUT) :: etest(ntol)
      REAL(r8),    INTENT(INOUT) :: errsum
      LOGICAL,     INTENT(INOUT) :: errok
!
      REAL(r8), PARAMETER :: one = 1.0D+0
      REAL(r8), PARAMETER :: zero = 0.0D+0
!
      REAL(r8) :: er, denom, errel
      INTEGER(i4) :: im, it, icmp
!
!     Given current and previous solutions u and uold on the same
!     mesh, errest calculates an error measure for each
!     component for which a tolerance is specified.
!     The error measure is the usual relative error normalized
!     by dividing by the tolerance.  On exit, errsum is the
!     sum of these error measures.
!
!     The array etest specifies the error test to be applied to each
!     error measure.
!
!     On exit, the logical flag errok
!      -- is .false. if any of the error measures exceeds the
!         corresponding value in the array etest
!      -- is .true. if all error measures are less than the
!         corresponding values of etest.
!!!
!
      IF (pdebug) WRITE(6,900)
      errsum = zero
      errok = .true.
!
      DO im = 1, nmsh
         DO it = 1, ntol
            icmp = ltol(it)
            er = u(icmp,im) - uold(icmp,im)
            denom = MAX(one, ABS(uold(icmp,im)))
            errel = ABS(er/(tol(it)*denom))
            IF (pdebug)  WRITE(6,901) im, it, errel, etest(it)
            errsum = errsum + errel
            IF (errel > etest(it)) errok = .false.
         END DO
      END DO
!
      IF (pdebug) WRITE(6,902) errsum
!
900   FORMAT(1H ,'errest')
901   FORMAT(1H ,2I5,2(1PE11.3))
902   FORMAT(1H ,'errsum',1PE11.3)
!
      RETURN
   END SUBROUTINE ERREST
!!!
!!!
   SUBROUTINE GETPTQ(debug, mfsrch, nout, alfmax, alfsml, alfuzz,  &
                     epsaf, epsag, eta, ftry, oldf, oldg, rmu, tolabs, tolrel, toltny, imprvd,  &
                     inform, nfsrch, alfa, alfbst, fbest,  &
                     braktd, crampd, extrap,vset,wset,nsamea,nsameb,  &
                     a, b, fa, factor, fv, fw, xtry, xv, xw )
!
      IMPLICIT NONE
!
      LOGICAL,     INTENT(INOUT) :: debug
      INTEGER(i4), INTENT(INOUT) :: mfsrch
      INTEGER(i4), INTENT(INOUT) :: nout
      REAL(r8),    INTENT(INOUT) :: alfmax
      REAL(r8),    INTENT(INOUT) :: alfsml
      REAL(r8),    INTENT(INOUT) :: alfuzz
      REAL(r8),    INTENT(INOUT) :: epsaf
      REAL(r8),    INTENT(INOUT) :: epsag
      REAL(r8),    INTENT(INOUT) :: eta
      REAL(r8),    INTENT(INOUT) :: ftry
      REAL(r8),    INTENT(INOUT) :: oldf
      REAL(r8),    INTENT(INOUT) :: oldg
      REAL(r8),    INTENT(INOUT) :: rmu
      REAL(r8),    INTENT(INOUT) :: tolabs
      REAL(r8),    INTENT(INOUT) :: tolrel
      REAL(r8),    INTENT(INOUT) :: toltny
      LOGICAL,     INTENT(INOUT) :: imprvd
      INTEGER(i4), INTENT(INOUT) :: inform
      INTEGER(i4), INTENT(INOUT) :: nfsrch
      REAL(r8),    INTENT(INOUT) :: alfa
      REAL(r8),    INTENT(INOUT) :: alfbst
      REAL(r8),    INTENT(INOUT) :: fbest
      LOGICAL,     INTENT(INOUT) :: braktd
      LOGICAL,     INTENT(INOUT) :: crampd
      LOGICAL,     INTENT(INOUT) :: extrap
      LOGICAL,     INTENT(INOUT) :: vset
      LOGICAL,     INTENT(INOUT) :: wset
      INTEGER(i4), INTENT(INOUT) :: nsamea
      INTEGER(i4), INTENT(INOUT) :: nsameb
      REAL(r8),    INTENT(INOUT) :: a
      REAL(r8),    INTENT(INOUT) :: b
      REAL(r8),    INTENT(INOUT) :: fa
      REAL(r8),    INTENT(INOUT) :: factor
      REAL(r8),    INTENT(INOUT) :: fv
      REAL(r8),    INTENT(INOUT) :: fw
      REAL(r8),    INTENT(INOUT) :: xtry
      REAL(r8),    INTENT(INOUT) :: xv
      REAL(r8),    INTENT(INOUT) :: xw
!
!     getptq  is a step-length algorithm for minimizing a function of one
!     variable.  it will be called repeatedly by a search routine whose
!     purpose is to estimate a point  alfa = alfbst  that minimizes some
!     function  f(alfa)  over the closed interval (0, alfmax).
!
!     getptq  requires the function  f(alfa)  (but not its gradient)
!     to be evaluated at various points within the interval.  new
!     step-length estimates are computed using quadratic interpolation with
!     safeguards.
!
!     reverse communication is used to allow the calling program to
!     evaluate  f.  some of the parameters must be set or tested
!     by the calling program.  the remainder would ordinarily be local
!     variables.
!
!
!     input parameters (relevant to the calling program)
!     --------------------------------------------------
!
!     debug         specifies whether detailed output is wanted.
!
!     inform        must be nonzero on the first entry (e.g., -1).
!                   it will be altered by  getptq  for later entries.
!
!     mfsrch        is an upper limit on the number of times  getptq  is
!                   to be entered consecutively with  inform = 0
!                   (following an initial entry with  inform lt 0).
!
!     nout          is the file number to be used for printed output
!                   if debug is true.
!
!     alfa          is the first estimate of the step length.  alfa  is
!                   subsequently altered by  getptq  (see below).
!
!     alfmax        is the upper limit of the interval to be searched.
!
!     alfsml        is intended to prevent inefficiency when the optimum
!                   step is very small, for cases where the calling
!                   program would prefer to re-define  f(alfa).  alfsml is
!                   allowed to be zero. early termination will occur if
!                   getptq  determines that the optimum step lies
!                   somewhere in the interval  (0, alfsml)  (but not if
!                   alfmax .le. alfsml).
!
!     epsaf         is an estimate of the absolute precision in the
!                   computed values of  f.
!
!     eta           controls the accuracy of the search.  it must lie
!                   in the range   0.0  le  eta  lt  1.0.  decreasing
!                   eta  tends to increase the accuracy of the search.
!
!     oldf          is the value of  f(0).
!
!     oldg          is an estimate of the gradient of  f  at  alfa = 0.
!                   it should be non-positive.
!
!     rmu           controls what is meant by a significant decrease in  f.
!                   the final  f(alfbst)  should lie on or below the line
!                         l(alfa)  =  oldf + alfa*rmu*oldg.
!                   rmu  should be in the open interval (0, 0.5).
!                   the value  rmu = 1.0d-4  is good for most purposes.
!
!     tolabs,tolrel define a function  tol(alfa) = tolrel*alfa + tolabs
!                   such that if  f  has already been evaluated at step
!                   alfa,  then it will not be evaluated at any point
!                   closer than  tol(alfa).
!                   these values may be reduced by  getptq  if they seem
!                   to be too large.
!
!     toltny        is the smallest value that  tolabs  is allowed to be
!                   reduced to.
!
!
!     output parameters (relevant to the calling program)
!     ---------------------------------------------------
!
!     imprvd        is true if the previous step  alfa  was the best
!                   point so far.  any related quantities (e.g., arrays)
!                   should be saved by the calling program before paying
!                   attention to  inform.
!
!     inform = 0    means the calling program should evaluate
!                              ftry = f(alfa)
!                   for the new trial step  alfa,  and then re-enter
!                   getptq.
!
!     inform = 1    means the search has terminated successfully
!                   with a step  alfbst  that is less than the
!                   upper bound  alfmax.
!
!     inform = 2    means the search has terminated successfully
!                   with a step  alfbst  that is equal to the
!                   upper bound  alfmax.
!
!     inform = 3    means that the search failed to find a point of
!                   sufficient decrease in  mfsrch  functions, but an
!                   improved point was found.
!
!     inform = 4    means  alfmax  is so small that a search should
!                   not have been done.
!
!     inform = 5    means that the search was terminated prematurely
!                   because of the value of  alfsml  (see above).
!
!     inform = 6    means the search has failed to find a useful step.  if
!                   the subroutine for the function and gradient has been
!                   programmed correctly, this will usually occur if the
!                   minimum lies very close to  alfa = 0  or the gradient
!                   is not sufficiently accurate.
!
!     inform = 7    means that the value of  g(0) was positive on entry.
!
!     alfa          is the step at which the next function value must be
!                   computed.
!
!     alfbst        should be accepted by the calling program as the
!                   required step-length estimate, whenever  getptq
!                   returns  inform = 1,  2  or  3.
!
!     fbest         will be the corresponding value of  f.
!
!
!     the following parameters retain information between entries
!     -----------------------------------------------------------
!
!     alfuzz        is such that, if the final  alfa  lies in the interval
!                   (0,alfuzz)  and  abs( f(alfa)-oldf ) le epsaf,  alfa
!                   cannot be guaranteed to be a point of sufficient
!                   decrease.
!
!     braktd        is false if  f  has not been evaluated at the far end
!                   of the interval of uncertainty.  in this case, the
!                   point  b  will be at  alfmax + tol(alfmax).
!
!     crampd        is true if  alfmax  is very small (le tolabs).
!                   if the search fails, this indicates that a zero
!                   step should be taken.
!
!     extrap        is true if alfbst has moved at least once and  xv
!                   lies outside the interval of uncertainty.  in this
!                   case, extra safeguards are applied to allow for
!                   instability in the polynomial fit.
!
!     vset          records whether a third-best point has been
!                   determined.
!
!     wset          records whether a second-best point has been
!                   determined.  it will always be true by the
!                   time the convergence test is applied (label 300).
!
!     nsamea        is the number of consecutive times that the left-hand
!                   end of the interval of uncertainty has remained the
!                   same.
!
!     nsameb        similarly for the right-hand end.
!
!     a, b, alfbst  define the current interval of uncertainty.
!                   the required minimum lies somewhere within the
!                   closed interval  (alfbst + a, alfbst + b).
!
!     alfbst        is the best point so far.  it is strictly within the
!                   the interval of uncertainty except when it lies at the
!                   left-hand end when  alfbst  has not been moved.
!                   hence we have    a le 0,   b gt 0.
!
!     fbest         is the value of  f  at the point  alfbst.
!
!     fa            is the value of  f  at the point  alfbst + a.
!
!     factor        controls the rate at which extrapolated estimates of
!                   alfa  may expand into the interval of uncertainty.
!                   factor is not used if the minimum has been bracketed
!                   (i.e., when the variable  braktd  is true).
!
!     fv, fw        are the values of  f  at the points  alfbst + xv,
!                   alfbst + xw.  they are not defined until  vset
!                   or  wset  (respectively) is true.
!
!     ftry          is the value of  f  at the new point  alfbst + xtry.
!
!     xtry          is the trial point within the shifted interval (a, b).
!                   the new trial function value must be computed at the
!                   point  alfa  =  alfbst + xtry.
!
!     xv            is such that  alfbst + xv  is the third-best point.
!                   it is not defined until  vset  is true.
!
!     xw            is such that  alfbst + xw  is the second-best point.
!                   it is not defined until  wset  is true.
!                   in some cases,  xw  will replace a previous  xw  that
!                   has a lower function but has just been excluded from
!                   the interval of uncertainty.
!
!
!  systems optimization laboratory, stanford university, california.
!  original version february 1982.  rev. may 1983.
!
!!!
!
      REAL(r8), PARAMETER :: zero = 0.0D+0, point1 = 0.1D+0, half = 0.5D+0
      REAL(r8), PARAMETER :: one = 1.0D+0, two = 2.0D+0, five = 5.0D+0
      REAL(r8), PARAMETER :: ten = 10.0D+0, eleven = 11.0D+0
!
      REAL(r8) :: tol, deltaf, ctry, cbest, atrue, btrue, alfaw, gap
      REAL(r8) :: alfav, xmidpt, q, s, gw, gv, artifa, artifb, endpnt, dtry, daux
      REAL(r8) :: xdif, xint, xrat
      LOGICAL :: closef, conv1, conv2, conv3, convrg
      LOGICAL :: moved, sigdec, xinxw
      INTEGER(i4) :: im, k
!
!     local variables
!     ---------------
!
!     closef        is true if the worst function  fv  is within  epsaf
!                   of  fbest  (up or down).
!
!     convrg        will be set to true if at least one of the convergence
!                   conditions holds at  alfbst.
!
!     moved         is true if a better point has been found (alfbst gt 0).
!
!     sigdec        says whether  fbest  represents a significant decrease
!                   in the function, compared to the initial value  oldf.
!
!     xinxw         is true if  xtry  is in  (xw,0)  or  (0,xw).
!     ---------------------------------------------------------------------
!!!!
!
      imprvd = .false.
      IF (inform /= -1) GO TO 100
!
!     first entry.  initialize various quantities, check input data and
!     prepare to evaluate the function at the initial step  alfa.
!
      nfsrch = 0
      alfbst = zero
      fbest  = oldf
      IF (oldg   >      zero) GO TO 970
      IF (oldg   >= (- epsag)) GO TO 960
      IF (alfmax <=    toltny) GO TO 940
!
      braktd = .false.
      crampd = alfmax <= tolabs
      extrap = .false.
      vset   = .false.
      wset   = .false.
      nsamea = 0
      nsameb = 0
      alfuzz = two*epsaf/(rmu*ABS( oldg ))
      a      = zero
      b      = alfmax + (tolrel*alfmax + tolabs)
      fa     = oldf
      factor = five
      tol    = tolabs
      xtry   = alfa
      IF (debug) WRITE (nout, 1000) alfmax, oldf, oldg, tolabs,  &
         alfuzz, epsaf, epsag, tolrel, crampd
      GO TO 800
!
!  ---------------------------------------------------------------------
!  subsequent entries.
!  the function has just been evaluated at  alfa = alfbst + xtry,
!  giving  ftry.
!  ---------------------------------------------------------------------
!
100   nsamea = nsamea + 1
      nsameb = nsameb + 1
      xtry   = alfa - alfbst
      moved  = alfbst > zero
!
!     check if  xtry  is in the interval  (xw,0)  or  (0,xw).
!
      xinxw  = .false.
      IF (wset) xinxw =       zero < xtry  .AND.  xtry <= xw  &
          .OR.    xw <= xtry  .AND.  xtry < zero
!
!     see if the new step is better.
!
      deltaf = ftry   - oldf
      ctry   = deltaf - alfa*rmu*oldg
      IF (alfa <= alfuzz) sigdec = deltaf <= (- epsaf)
      IF (alfa > alfuzz) sigdec = ctry   <=    epsaf
      imprvd = sigdec  .AND.  ( ftry - fbest ) <= (- epsaf)
!
      IF (debug) WRITE (nout, 1100) alfa, ftry, ctry
      IF (.NOT. imprvd) GO TO 130
!
!     we seem to have an improvement.  the new point becomes the
!     origin and other points are shifted accordingly.
!
      IF (.NOT. wset) GO TO 110
      xv     = xw - xtry
      fv     = fw
      vset   = .true.
110   xw     = zero - xtry
      fw     = fbest
      wset   = .true.
      fbest  = ftry
      alfbst = alfa
      a      =    a - xtry
      b      =    b - xtry
      moved  = .true.
      extrap = .NOT. xinxw
!
!     decrease the length of the interval of uncertainty.
!
      IF (xtry < zero) GO TO 120
      a      = xw
      fa     = fw
      nsamea = 0
      GO TO 300
120   b      = xw
      nsameb = 0
      braktd = .true.
      GO TO 300
!
!     the new function value is no better than the best point found so far.
!     the point  xtry  must be a new end point of the interval of
!     uncertainty.
!
130   IF (xtry >= zero) GO TO 140
      a      = xtry
      fa     = ftry
      nsamea = 0
      GO TO 150
140   b      = xtry
      nsameb = 0
      braktd = .true.
!
!     the origin remains unchanged but  xtry  may qualify as  xw.
!
150   IF (.NOT. wset)   GO TO 160
      IF ((ftry - fw) > epsaf) GO TO 170
      xv     = xw
      fv     = fw
      vset   = .true.
160   xw     = xtry
      fw     = ftry
      wset   = .true.
      IF (moved) extrap = xinxw
      GO TO 300
!
!     ftry  is no better than  fbest  or  fw.  if the best point has not
!     been moved, there must be more than one minimum.
!
170   IF (moved) GO TO 175
      xw     = xtry
      fw     = ftry
      GO TO 300
!
!     ftry  is no better than  fbest  or  fw,  but  xtry  may become  xv.
!     extrap  has the value set in the previous entry.
!
175   IF (.NOT. vset) GO TO 180
      IF ((ftry - fv) > epsaf  .AND.  extrap) GO TO 300
180   IF (xinxw) GO TO 190
      xv     = xtry
      fv     = ftry
      vset   = .true.
      GO TO 300
190   IF (vset) xw = xv
      IF (vset) fw = fv
      xv     = xtry
      fv     = ftry
      vset   = .true.
!
!     ---------------------------------------------------------------------
!     check the termination criteria.
!     ---------------------------------------------------------------------
300   tol    = tolrel*alfbst + tolabs
      deltaf = fbest  - oldf
!
      cbest  = deltaf - alfbst*rmu*oldg
      IF (alfbst <= alfuzz) sigdec = deltaf <= (- epsaf)
      IF (alfbst > alfuzz) sigdec = cbest  <=    epsaf
      closef = .false.
      IF (vset) closef = ABS( fbest - fv ) <= epsaf
!
      conv1  = MAX( ABS( a ), b )  <=  (tol + tol)
!
!     conv2 changed by mhw, 20 sept 1992, to allow it to be
!     satified for any significant decrease in f
!      conv2  =  moved  .and.  sigdec
!     *                 .and.  abs( fa - fbest )  .le.  a*eta*oldg
!
      conv2 = moved .AND. sigdec
      conv3  = closef  .AND.  (sigdec  .OR. (.NOT. moved)  .AND.  (b <= alfuzz))
      convrg = conv1  .OR.  conv2  .OR.  conv3
!
      atrue  = alfbst + a
      btrue  = alfbst + b
      alfaw  = alfbst + xw
      gap    = b - a
      IF (debug) WRITE (nout, 1200) atrue, btrue, gap, tol,  &
          nsamea, nsameb, braktd, closef, imprvd, conv1, conv2, conv3,  &
          extrap, alfbst, fbest, cbest, alfaw, fw
      IF (vset) alfav  = alfbst + xv
      IF (debug  .AND.  vset) WRITE (nout, 1300) alfav, fv
      IF (convrg  .AND.  moved) GO TO 910
!
!     exit if the step is too small.
!
      IF (btrue   <  alfsml) GO TO 950
!
      IF (nfsrch  >=  mfsrch) GO TO 930
      IF (.NOT. convrg) GO TO 400
!
!     a better point has not yet been found (the step  xw  is no better
!     than step  zero).  check that the change in  f  is consistent with a
!     perturbation in  x  of  tol, the estimate of the minimum spacing
!     constant.  if the change in  f  is larger than  epsaf,  the value
!     of  tol  is reduced.
!
      tol    = xw/ten
      tolabs = tol
      IF (ABS(fw - oldf) > epsaf  .AND.  tol > toltny) GO TO 400
      IF (crampd) GO TO 940
      GO TO 960
!
!     ---------------------------------------------------------------------
!     proceed with the computation of a trial step length.
!     the choices are...
!     1. parabolic fit using function values only.
!     2. damped parabolic fit if the regular fit appears to be
!        consistently over-estimating the distance to the minimum.
!     3. bisection, geometric bisection, or a step of  tol  if the
!        parabolic fit is unsatisfactory.
!     ---------------------------------------------------------------------
400   xmidpt = half*(a + b)
      q      = zero
      s      = zero
!
!     ---------------------------------------------------------------------
!     fit a parabola.
!     ---------------------------------------------------------------------
!
!     check if there are two or three points for the parabolic fit.
!
      gw = (fw - fbest)/xw
      IF (vset  .AND.  moved) GO TO 450
!
!     only two points available.  use  fbest,  fw  and the derivative
!     oldg.
!
      IF (.NOT. moved) s = oldg
      IF (      moved) s = oldg - two*gw
      q = two*(oldg - gw)
      IF (debug) WRITE (nout, 2100)
      GO TO 600
!
!     three points available.  use  fbest,  fw  and  fv.
!
450   gv = (fv - fbest)/xv
      s  = gv - (xv/xw)*gw
      q  = two*(gv - gw)
      IF (debug) WRITE (nout, 2200)
!
!     ---------------------------------------------------------------------
!     construct an artificial interval  (artifa, artifb)  in which the
!     new estimate of the step length must lie.  set a default value of
!     xtry  that will be used if the polynomial fit is rejected.  in the
!     following, the interval  (a,b)  is considered the sum of two
!     intervals of lengths  dtry  and  daux, with common end point at the
!     best point (zero).  dtry  is the length of the interval into which
!     the default  xtry  will be placed and  endpnt  denotes its non-zero
!     end point.  the magnitude of  xtry  is computed so that the exponents
!     of  dtry  and  daux  are approximately bisected.
!     ---------------------------------------------------------------------
600   artifa = a
      artifb = b
      IF (braktd) GO TO 610
!
!     the minimum has not been bracketed.  set an artificial upper bound
!     by expanding the interval  xw  by a suitable factor.
!
      xtry   = - factor*xw
      artifb =   xtry
      IF (alfbst + xtry < alfmax) factor = five*factor
      GO TO 700
!
!     the minimum has been bracketed.
!     if the gradient at the origin is being used for the
!     polynomial fit, the default  xtry  is one tenth of  xw.
!
610   IF (vset  .AND.  moved) GO TO 620
      xtry   = xw/ten
      IF (debug) WRITE (nout, 2400) xtry
      GO TO 700
!
!     three points exist in the interval of uncertainty.  check whether
!     the points are configured for an extrapolation or interpolation.
!
620   IF (extrap) GO TO 660
!
!     if the interpolation appears to be consistently over-estimating the
!     distance to the minimum,  damp the interpolation step.
!
      IF (nsamea < 3  .AND.  nsameb < 3) GO TO 630
      factor = factor / five
      s      = factor * s
      GO TO 640
630   factor = one
!
!     the points are configured for an interpolation.  the artificial
!     interval will be just  (a,b).   set  endpnt  so that  xtry
!     lies in the larger of the intervals  (a,0)  and  (0,b).
!
640   IF (xmidpt < zero) endpnt = a
      IF (xmidpt > zero) endpnt = b
!
!     if a bound has remained the same for three iterations, set  endpnt
!     so that  xtry  is likely to replace the offending bound.
!
      IF (nsamea >= 3) endpnt = a
      IF (nsameb >= 3) endpnt = b
      GO TO 680
!
!     the points are configured for an extrapolation.
!
660   IF (xw < zero) endpnt = b
      IF (xw > zero) endpnt = a
!
!     compute the default value of  xtry.
!
680   dtry = ABS( endpnt )
      daux = gap - dtry
      IF (daux >= dtry)   xtry = five*dtry*(point1 + dtry/daux)/eleven
      IF (daux < dtry)   xtry = half*SQRT( daux )*SQRT( dtry )
      IF (endpnt < zero) xtry = - xtry
      IF (debug) WRITE (nout, 2500) xtry, daux, dtry
!
!     if the points are configured for an extrapolation set the artificial
!     bounds so that the artificial interval lies strictly within  (a,b).
!     if the polynomial fit is rejected,  xtry  will remain at the relevant
!     artificial bound.
!
      IF (extrap  .AND.  xtry <= zero) artifa = xtry
      IF (extrap  .AND.  xtry > zero) artifb = xtry
!
!     ---------------------------------------------------------------------
!     the polynomial fits give  (s/q)*xw  as the new step.
!     reject this step if it lies outside  (artifa, artifb).
!     ---------------------------------------------------------------------
700   IF (q == zero) GO TO 800
      IF (q < zero) s = - s
      IF (q < zero) q = - q
      IF (s*xw < q*artifa   .OR.   s*xw > q*artifb) GO TO 800
!
!     accept the polynomial fit.
!
      xtry = zero
      IF (ABS( s*xw ) >= q*tol) xtry = (s/q)*xw
      IF (debug) WRITE (nout, 2600) xtry
!
!     ---------------------------------------------------------------------
!     test for  xtry  being larger than  alfmax  or too close to  a  or  b.
!     ---------------------------------------------------------------------
800   IF (braktd) GO TO 810
!
!     if the step is close to or larger than  alfmax,  replace it by
!     alfmax  (to force evaluation of the function at the boundary).
!
      alfa   = alfbst + xtry
      IF (alfmax - alfa > (tolrel*alfmax + tolabs)) GO TO 810
      braktd = .true.
      xtry   = alfmax - alfbst
      alfa   = alfmax
      GO TO 900
!
!     otherwise, the function must not be evaluated at a point too close
!     to  a  or  b.  (it has already been evaluated at both those points.)
!
810   xmidpt = half*(a + b)
      IF (xtry > a + tol  .AND.  xtry < b - tol) GO TO 820
      IF (xmidpt > zero) xtry =   tol
      IF (xmidpt <= zero) xtry = - tol
!
!     f  must not be calculated too close to  alfbst.
!
820   IF (ABS( xtry ) < tol  .AND.  xmidpt < zero) xtry = - tol
      IF (ABS( xtry ) < tol  .AND.  xmidpt >= zero) xtry =   tol
      alfa   = alfbst + xtry
!
!     ---------------------------------------------------------------------
!     exit.
!     ---------------------------------------------------------------------
!
!     new function value required.
!
900   inform = 0
      GO TO 990
!
!     convergence test satisfied.
!
910   inform = 1
      IF (alfa == alfmax) inform = 2
      GO TO 990
!
!     mfsrch  function evaluations without sufficient decrease, but an
!     improved point was found.
!
930   IF (.NOT. moved) GO TO 960
      inform = 3
      GO TO 990
!
!     zero step (alfmax too small).
!
940   inform = 4
      GO TO 990
!
!     premature termination.  the step is smaller than  alfsml.
!
950   inform = 5
      GO TO 990
!
!     zero step (a sufficiently better point could not be found).
!
960   inform = 6
      GO TO 990
!
!     zero step (positive gradient at the starting point).
!
970   inform = 7
!
!     exit.
!
990   IF (debug) WRITE (nout, 3000)

1000  FORMAT(/ 31H alfmax  oldf    oldg    tolabs, 1P2E22.14, 1P2E16.8  &
    / 31H alfuzz  epsaf   epsag   tolrel, 1P2E22.14, 1P2E16.8  &
    / 31H crampd                        ,  l6)
1100  FORMAT(/ 31H alfa    ftry    ctry          , 1P2E22.14, 1PE16.8)
1200  FORMAT(/ 31H a       b       b - a   tol   , 1P2E22.14, 1P2E16.8  &
    / 31H nsamea  nsameb  braktd  closef, 2I3, 2L6  &
    / 31H imprvd  convrg  extrap        ,  l6, 3X, 3L1, l6  &
    / 31H alfbst  fbest   cbest         , 1P2E22.14, 1PE16.8  &
    / 31H alfaw   fw                    , 1P2E22.14)
1300  FORMAT(  31H alfav   fv                    , 1P2E22.14 /)

2100  FORMAT(30H parabolic fit,    two points.)
2200  FORMAT(30H parabolic fit,  three points.)
2400  FORMAT(31H exponent reduced.  trial point, 1P1E22.14)
2500  FORMAT(31H geo. bisection. xtry,daux,dtry, 1P3E22.14)
2600  FORMAT(31H polynomial fit accepted.  xtry, 1P1E22.14)
3000  FORMAT(53H ---------------------------------------------------- /)
!
      RETURN
!  end of getptq
   END SUBROUTINE GETPTQ
!!!
!!!
   SUBROUTINE INTERP(ncomp, nmsh, xx, nudim, u, nmold, xxold, uold)
!
      USE BVPShared, ONLY: pdebug, use_c, comp_c, uval0,  &
                           nminit, iprint, idum,  &
                           DCOPY
!
      IMPLICIT NONE
!
      INTEGER(i4), INTENT(IN)    :: ncomp
      INTEGER(i4), INTENT(INOUT) :: nmsh
      REAL(r8),    INTENT(INOUT) :: xx(*)
      INTEGER(i4), INTENT(IN)    :: nudim
      REAL(r8),    INTENT(INOUT) :: u(nudim,*)
      INTEGER(i4), INTENT(INOUT) :: nmold
      REAL(r8),    INTENT(INOUT) :: xxold(*)
      REAL(r8),    INTENT(INOUT) :: uold(ncomp,*)
!
      REAL(r8), PARAMETER :: zero = 0.0D+0
!
      REAL(r8) :: xdif, xint, xrat
      INTEGER(i4) :: i, im, k, it, icmp, imnew
!
!     interp performs piecewise linear interpolation of the old
!     solution uold at the nmold old mesh points xxold onto the nmsh
!     new mesh points xx, producing an interpolated solution u.
!     Note that no assumption is made that the new mesh has
!     more points than the old, nor that the new and old mesh
!     points are related in a specific way (except that their first
!     and last points are identical).
!
      IF (pdebug) WRITE(6,900)
!
!     By construction, xx(1) = xxold(1).  Copy the first ncomp
!     components of uold into those of u.
!
      CALL DCOPY(ncomp, uold(1,1), 1, u(1,1), 1)
!
      i = 2
      DO im = 2, nmsh-1
!  
50       CONTINUE
         IF (i > nmold) RETURN
! 
!        Check whether the im-th point in the new mesh lies strictly
!        to the right of, or to the left of (or exactly on) the
!        i-th point in the old mesh.
! 
         IF (xx(im) > xxold(i)) THEN
            i = i + 1
            GO TO 50
         ELSE
            xdif = xxold(i) - xx(im)
            IF (xdif == zero) THEN
! 
!              xx(im) and xxold(i) are identical.
! 
               CALL DCOPY(ncomp, uold(1,i), 1, u(1,im), 1)
               i = i + 1
            ELSE
               xint = xxold(i) - xxold(i-1)
               xrat = xdif/xint
               DO k = 1, ncomp
                  u(k,im) = uold(k,i) + xrat*(uold(k,i-1)-uold(k,i))
               END DO
            END IF
         END IF
!
      END DO
      CALL DCOPY(ncomp, uold(1,nmold), 1, u(1,nmsh), 1)
!
900   FORMAT(1H ,'interp')
!
      RETURN
   END SUBROUTINE INTERP
!!!
!!!
   SUBROUTINE RERRVL(ncomp, nmsh, nudim, u, usvrex, ntol, ltol,  &
                     rerr, remax, itlmx, adjrer)
!
      IMPLICIT NONE
!
      INTEGER(i4), INTENT(IN)    :: ncomp
      INTEGER(i4), INTENT(INOUT) :: nmsh
      INTEGER(i4), INTENT(IN)    :: nudim
      REAL(r8),    INTENT(INOUT) :: u(nudim, *)
      REAL(r8),    INTENT(INOUT) :: usvrex(ncomp, *)
      INTEGER(i4), INTENT(IN)    :: ntol
      INTEGER(i4), INTENT(IN)    :: ltol(*)
      REAL(r8),    INTENT(INOUT) :: rerr(ncomp, *)
      REAL(r8),    INTENT(INOUT) :: remax
      INTEGER(i4), INTENT(INOUT) :: itlmx
      LOGICAL,     INTENT(INOUT) :: adjrer
!
      REAL(r8), PARAMETER :: one = 1.0D+0
      REAL(r8), PARAMETER :: zero = 0.0D+0
!
      REAL(r8) :: denom, rerel, rlndec, frcpow
      INTEGER(i4) :: nmold, it, icmp, imnew, im, ninter, ithres
!
!     rerrvl is used in considering Richardson extrapolation.
!     The two solutions u and usvrex have a special relationship:
!     u corresponds to a doubled mesh, with twice as many
!     intervals, where each interval is half the size of that for
!     usvrex's mesh.   nmsh is the number of mesh points in the
!     mesh for u.
!
!     remax is the maximum relative error, and itlmx is the
!     index of the tolerance for which the maximum occurs.
!
!     The array rerr contains the absolute value of the difference
!     between u and usvrex at the common mesh points, but is defined
!     only for components for which an error tolerance is specified.
!
!     The logical variable adjrer is true on entry if the values in
!     rerr are to be adjusted for later use in selective mesh refinement.
!
!!!
!
      itlmx = 1
      remax = zero
      nmold = 1 + (nmsh-1)/2
      DO it = 1, ntol
         icmp = ltol(it)
         imnew = 1
         DO im = 1, nmold
            rerr(icmp, im) = ABS(usvrex(icmp,im) - u(icmp,imnew))
            denom = MAX(one, ABS(usvrex(icmp,im)))
            rerel = rerr(icmp,im)/denom
            IF (rerel > remax) THEN
               remax = rerel
               itlmx = it
            END IF
            imnew = imnew + 2
         END DO
      END DO
!
      IF (adjrer) THEN
!
!        Adjust the rerr array if it may be used later in selective
!        mesh refinement.
!
         DO it = 1, ntol
            icmp = ltol(it)
            DO im = 1, nmold - 1
               rerr(icmp,im) = MAX(rerr(icmp,im), rerr(icmp, im+1))
            END DO
         END DO
!
      END IF
!
      RETURN
   END SUBROUTINE RERRVL
!!!
!!!
   SUBROUTINE SELCONDERRMSH(ncomp, nmsh, ntol, ltol, tol,  &
                            nfxpnt, fixpnt, ipow, nmax,  &
                            xx, nudim, u, ermeas, irefin, ihcomp, nmold, xxold,  &
                            ermx, ddouble, maxmsh, r4, amg, stab_cond)
!
      USE BVPShared, ONLY: pdebug, use_c, comp_c, uval0,  &
                           nminit, iprint, idum
!
      IMPLICIT NONE
!
      INTEGER(i4), INTENT(IN)    :: ncomp
      INTEGER(i4), INTENT(INOUT) :: nmsh
      INTEGER(i4), INTENT(IN)    :: ntol
      INTEGER(i4), INTENT(IN)    :: ltol(ntol)
      REAL(r8),    INTENT(IN)    :: tol(ntol)
      INTEGER(i4), INTENT(IN)    :: nfxpnt
      REAL(r8),    INTENT(IN)    :: fixpnt(*)
      INTEGER(i4), INTENT(INOUT) :: ipow
      INTEGER(i4), INTENT(INOUT) :: nmax
      REAL(r8),    INTENT(INOUT) :: xx(*)
      INTEGER(i4), INTENT(IN)    :: nudim
      REAL(r8),    INTENT(INOUT) :: u(nudim,*)
      REAL(r8),    INTENT(INOUT) :: ermeas(ncomp,*)
      INTEGER(i4), INTENT(INOUT) :: irefin(nmsh-1)
      INTEGER(i4), INTENT(INOUT) :: ihcomp(nmsh-1)
      INTEGER(i4), INTENT(INOUT) :: nmold
      REAL(r8),    INTENT(INOUT) :: xxold(*)
      REAL(r8),    INTENT(INOUT) :: ermx(*)
      LOGICAL,     INTENT(INOUT) :: ddouble
      LOGICAL,     INTENT(INOUT) :: maxmsh
      REAL(r8),    INTENT(INOUT) :: r4(nmsh)
      REAL(r8),    INTENT(INOUT) :: amg(*)
      LOGICAL,     INTENT(INOUT) :: stab_cond
!
      REAL(r8), PARAMETER :: zero = 0.0D+0
      REAL(r8), PARAMETER :: one = 1.0D+0
      REAL(r8), PARAMETER :: onep1 = 1.1D+0
      REAL(r8), PARAMETER :: erdcid = 5.0D+0
      REAL(r8), PARAMETER :: phitst = 0.1D+0
!
      LOGICAL, SAVE :: first = .true.
      REAL(r8), SAVE :: rlndec
      REAL(r8) :: frcpow, thres, errmax, denom, ems, err
      REAL(r8) :: r1, r2, r3, fatt_r1r3, fatt_r3, decii, rlen, slen, dx
      REAL(r8) :: fxnext, rlold, phihat, val1
      INTEGER(i4) :: ninter, ithres, im, it, jcomp, nptcond, ii, ilg, nmest
      INTEGER(i4) :: new, i, j, ifxcnt, jtkout, ind1
      LOGICAL :: add
!
!     The routine selconderrmsh performs selective mesh refinement, depending
!     on the error measure ermeas and the monitor function based on the
!     conditioning.
!
      IF (first) THEN
         first = .false.
         rlndec = DLOG(erdcid)
      END IF
!
      maxmsh = .false.
!
      IF (pdebug) WRITE(6,901) nmsh, ipow
!
      frcpow = one/ipow
      ddouble = .false.
      nmold = nmsh
      ninter = nmsh - 1
!
!     Copy the current mesh into the xxold array.
!
      CALL DCOPY(nmold, xx, 1, xxold, 1)
      ithres = 0
      thres = one
!
!     On input, the array ermeas represents some error measure defined
!     over the components and mesh intervals (not mesh points).
!     It is normalized in the following loop with respect to the
!     tolerance array and the current solution.
!     The value errmax gives the maximum normalized error.
!
      errmax = zero
      DO im = 1, ninter
         ermx(im) = zero
         DO it = 1, ntol
            jcomp = ltol(it)
            denom = tol(it)*MAX(one, ABS(u(jcomp,im)))
            ems = ermeas(jcomp,im)
            ermeas(jcomp,im) = ABS(ems)/denom
!             if (pdebug .and. ermeas(jcomp,im) .ge. thres)
!     *              write(6,902) im,jcomp,ems,ermeas(jcomp,im)
            err = ermeas(jcomp, im)
            IF (err >= ermx(im)) THEN
               ermx(im) = ERR
               ihcomp(im) = jcomp
            END IF
         END DO
         errmax = MAX(ermx(im), errmax)
      END DO

!       if ((.not. stab_cond) .and. errmax .ge. 1.0d6  ) then
!f  only the conditioning
!            call  selcondmsh(ncomp, nmsh,
!     *         nfxpnt, fixpnt,  nmax, xx,  irefin,
!     *         nmold, xxold, ddouble , maxmsh,r4,amg)

!      else
!f   the conditining and the error
!
      CALL MONCONDMSH(nmsh,xx,r1,r2,r3,fatt_r1r3,fatt_r3,nptcond,r4,amg)
!
      IF (pdebug) WRITE(6,903) errmax
!
      IF (errmax > zero .AND. errmax <= erdcid) THEN
! 
!        If errmax > 0 and .le. erdcid, find the smallest integer exponent ii
!        such that (erdcid**ii)*errmax > erdcid.
! 
         IF (errmax > one) THEN
            ii = 1
            decii = erdcid
         ELSE
            ilg = -DLOG(errmax)/rlndec
            ii = 2 + ilg
            decii = erdcid**ii
         END IF
! 
!        Multiply error measures by erdcid**ii.
! 
         errmax = decii*errmax
         DO im = 1, ninter
            ermx(im) = decii*ermx(im)
            DO it = 1, ntol
               jcomp = ltol(it)
               ermeas(jcomp,im) = decii*ermeas(jcomp,im)
            END DO
         END DO
!
      END IF
!
200   CONTINUE
!
!     For each interval im,  the integer irefin(im) is calculated
!     based on an equidistrbution procedure involving the
!     threshold value thres.  If irefin(im) > 1, we add
!     points to interval im.  If irefin(im) = 1, we check whether
!     point im can be removed from the mesh.
!
!     nmest is a lower bound on the number of points in the new mesh.
!     We do not know in advance how many points will be removed,
!     so nmest is computed by assuming all eligible points are removed.
!
      nmest = nmsh
      DO im = 1, ninter
         IF (ermx(im) >= thres) THEN
            irefin(im) = INT(ermx(im)**frcpow) + 1
            nmest = nmest + irefin(im) - 1
         ELSE
            irefin(im) = 1
            nmest = nmest - 1
         END IF
      END DO
!
      IF (nptcond >= 4 ) THEN
         add = .false.
         nptcond = nptcond/2
         DO im = 1, ninter-1
            IF (MAX(r4(im), r4(im+1)) > fatt_r1r3) THEN
               IF (.NOT. add) THEN
                  irefin(im) = MAX(nptcond, irefin(im))
!                 nmest = nmest + nptcond - 1
               END IF
               irefin(im+1) = MAX(nptcond, irefin(im+1))
!              nmest = nmest + nptcond - 1
               add = .true.
            ELSE
               irefin(im) = MAX(1, irefin(im))
               irefin(im+1) = MAX(1, irefin(im+1))
!              nmest = nmest - 1
               add = .false.
            END IF
         END DO
! 
!        do 221 im = 1, ninter
!         if ( r4(im) .gt. fatt_r1r3) then
!              irefin(im) = max(nptcond, irefin(im))
!c              nmest = nmest + nptcond - 1
!         else
!              irefin(im) = max(1, irefin(im))
!c              nmest = nmest - 1
!         endif
! 221     continue
      END IF
!
      IF (pdebug) WRITE(6,904) nmest, (irefin(i), i=1,ninter)
!
      IF (nmest > nmax) THEN
         GO TO 360
      ELSE IF (nmest-1 > 3*ninter) THEN
         CALL DBLMSH(nmsh, nmax, xx, nmold, xxold, maxmsh)
         ddouble = .true.
         RETURN
      END IF
!
!     It appears that we can perform the desired selective mesh
!     refinement.
!
!     Now begin running through the mesh, adding and possibly deleting
!     points as indicated by the irefin array.
!
!     The integer new is a count of the number of intervals in
!     the tentative mesh being generated by the refinement strategy.
!
      new = 1
!
!     The first interval is treated as a special case, since xx(1)
!     always remains in the mesh, and cannot be a fixed point.
!
      rlen = xxold(2) - xx(1)
      slen = rlen
      IF (irefin(1) > 1) THEN
         dx = rlen/irefin(1)
         DO  j = 2, irefin(1)
            new = new + 1
            xx(new) = xx(1) + (j-1)*dx
         END DO
      END IF
!
!     The fixed points specified by the fixpnt array cannot be
!     removed from the mesh.  The value fxnext indicates the 'next'
!     fixed point.  When no further fixed points remain to be processed
!     (or if nfxpnt = 0), fxnext is set to a value strictly larger than
!     the last mesh point, so that no mesh point can equal fxnext.
!     This way we need to compare only the values of xxold(i)
!     and fxnext.
!
      ifxcnt = 1
      IF (nfxpnt == 0) THEN
         fxnext = onep1*ABS(xxold(nmsh))
      ELSE
         fxnext = fixpnt(ifxcnt)
      END IF
!
!     jtkout is a counter of the number of consecutive points that
!     have been removed from the mesh.
!
      jtkout = 0
!
      DO im = 2, ninter
!
         rlold = rlen
         rlen = xxold(im+1) - xxold(im)
! 
!        If xxold(im) is the next fixed point, it cannot be removed
!        and so we don't test its error estimates.
! 
         IF (xxold(im) == fxnext) THEN
!
            ifxcnt = ifxcnt + 1
            IF (ifxcnt > nfxpnt) THEN
               fxnext = onep1*ABS(xxold(nmsh))
            ELSE
               fxnext = fixpnt(ifxcnt)
            END IF
!
         ELSE IF (irefin(im) == 1) THEN
! 
!           If xxold(im) is not a fixed point and irefin(im) = 1, we may wish
!           to remove point im from the mesh.
! 
!           If we are considering removing points and jtkout = 0, this
!           is the first point in a possible consecutive set to be removed,
!           and we initialize phihat, which represents a maximum of
!           certain estimates.
!           If jtkout is not zero, previous points contiguous to this
!           point have been removed, and phihat retains its previous value.
! 
            slen = slen + rlen
! 
            IF (jtkout == 0) THEN
               ind1 = ihcomp(im-1)
               phihat = ermeas(ind1,im-1)/(rlold**ipow)
            END IF
            phihat = MAX(phihat, ermeas(ihcomp(im),im)/(rlen**ipow))
            val1 = phihat*(slen**ipow)
            IF (val1 <= phitst .AND. jtkout < 4  &
               .AND. r4(im) <= r3) THEN
! 
!              Increment the counter of removed points.
!              'Remove' the mesh point xxold(im) by not including it.
! 
               jtkout = jtkout+1
               CYCLE
            END IF
!        end of logic for irefin(im) = 1.
         END IF
! 
         jtkout = 0
         new = new + 1
         xx(new) = xxold(im)
         IF (irefin(im) > 1) THEN
            dx = rlen/irefin(im)
            DO j = 2, irefin(im)
               new = new + 1
               xx(new) = xxold(im) + (j-1)*dx
            END DO
         END IF
         slen = rlen
! 
         IF (new > nmax) THEN
! 
!           If the new mesh contains too many points, branch out of the
!           loop to try alternative strategies.
! 
            GO TO 360
! 
         ELSE IF (new > 3*ninter)  THEN
! 
!           Here, the new mesh does not exceed the specified maximum,
!           but has more than 3 times as many intervals as the old mesh.
!           Try doubling the mesh if possible.
!
            nmsh = nmold
            CALL DCOPY(nmsh, xxold, 1, xx, 1)
            CALL DBLMSH(nmsh, nmax, xx, nmold, xxold, maxmsh)
            ddouble = .true.
            RETURN
!
         END IF
! 
      END DO
!
!     To end up here, we have processed the entire interval,
!     and have neither exceeded the specified maximum nor
!     exceeded three times the number of intervals in the old
!     mesh.  The last mesh point remains unchanged.
!
      new = new + 1
      xx(new) = xxold(nmsh)
      nmsh = new
      maxmsh = .false.
      IF (iprint >= 0) WRITE(6,905) nmsh
      RETURN
!
360   CONTINUE
!
!     To reach here, the number of mesh points created at some stage
!     of the refinement process was larger than the maximum permitted
!     value nmax.
!
!     Check whether the mesh can safely be doubled.
!
      IF ((2*nmsh-1) < nmax) THEN
! 
!        Double the mesh.
         CALL DBLMSH(nmsh, nmax, xx, nmold, xxold, maxmsh)
         ddouble = .true.
! 
!        If the number of intervals is too large and the mesh cannot be
!        doubled, increase the threshold thres by a factor of erdcid and
!        try the selective refinement again.
!        If this happens three times without success or if thres exceeds
!        or is equal to errmax, stop.  (In this case, we know already
!        that doubling the mesh produces too many points.)
! 
      ELSE IF (thres < errmax .AND. ithres < 3) THEN
         ithres = ithres + 1
         thres = erdcid*thres
         IF (thres > errmax) thres = errmax
         CALL DCOPY(nmsh, xxold, 1, xx, 1)
         GO TO 200
      ELSE
!        nmsh = 2*nmsh - 1
         nmsh = nmold
         CALL DCOPY(nmsh, xxold, 1, xx, 1)
         maxmsh = .true.
      END IF
!
!     endif use the conditioning and the error
901   FORMAT(1H ,'selconderrmsh.  nmsh, ipow =',2I5)
902   FORMAT(1H ,'im, jcomp, ermeas, normalized er',2I5,2(1PE11.3))
903   FORMAT(1H ,'errmax',1PE11.3)
904   FORMAT(1H ,'nmest, irefin',(10I5))
905   FORMAT(1H ,'selconderrmsh.  new nmsh =',i8)
910   FORMAT(1H ,'ihcomp',(10I5))
!
      RETURN
   END SUBROUTINE SELCONDERRMSH
!!!
!!!
   SUBROUTINE SELCONDMSH(ncomp, nmsh, nfxpnt, fixpnt,  nmax, xx,  irefin,  &
                         nmold, xxold, ddouble , maxmsh, r4, amg)
!
      USE BVPShared, ONLY: pdebug, use_c, comp_c, uval0,  &
                           nminit, iprint, idum
!
      IMPLICIT NONE
!
      INTEGER(i4), INTENT(IN)    :: ncomp
      INTEGER(i4), INTENT(INOUT) :: nmsh
      INTEGER(i4), INTENT(IN)    :: nfxpnt
      REAL(r8),    INTENT(IN)    :: fixpnt(*)
      INTEGER(i4), INTENT(INOUT) :: nmax
      REAL(r8),    INTENT(INOUT) :: xx(*)
      INTEGER(i4), INTENT(INOUT) :: irefin(*)
      INTEGER(i4), INTENT(INOUT) :: nmold
      REAL(r8),    INTENT(INOUT) :: xxold(*)
      LOGICAL,     INTENT(INOUT) :: ddouble
      LOGICAL,     INTENT(INOUT) :: maxmsh
      REAL(r8),    INTENT(INOUT) :: r4(nmsh)
      REAL(r8),    INTENT(INOUT) :: amg(*)
!
      REAL(r8), PARAMETER :: zero = 0.0D+0
      REAL(r8), PARAMETER :: one = 1.0D+0
      REAL(r8), PARAMETER :: onep1 = 1.1D+0
      REAL(r8), PARAMETER :: erdcid = 5.0D+0
!
      REAL(r8) :: thres, r1, r2, r3, fatt_r1r3, fatt_r3, rlen, slen, dx
      REAL(r8) :: fxnext, rlold
      LOGICAL :: add
      INTEGER(i4) :: ninter, ithres, nptcond, nmest, im, new, j
      INTEGER(i4) :: ifxcnt, jtkout
!!!
!
!     The routine selcondmsh performs selective mesh refinement, depending
!     on the monitor function based on the conditioning parameters.
!
      maxmsh = .false.
!
      IF (pdebug) WRITE(6,901) nmsh
!
      ddouble = .false.
      nmold = nmsh
      ninter = nmsh - 1
!
!     Copy the current mesh into the xxold array.
!
      CALL DCOPY(nmold, xx, 1, xxold, 1)
!
      ithres = 0
      thres = one
!
!     On input, the array amg represents the conditioning vector defined
!     over the components and mesh intervals (not mesh points).
!     We compute the monitor function and the related parameters in the
!     followinf function
!
      CALL MONCONDMSH(nmsh, xx, r1, r2, r3, fatt_r1r3, fatt_r3, nptcond, r4, amg)
!
!     For each interval im,  the integer irefin(im) is calculated
!     based on an equidistrbution procedure involving the
!     threshold value thres.  If irefin(im) > 1, we add
!     points to interval im.  If irefin(im) = 1, we check whether
!     point im can be removed from the mesh.
!
!     nmest is a lower bound on the number of points in the new mesh.
!     We do not know in advance how many points will be removed,
!     so nmest is computed by assuming all eligible points are removed.
!
      add = .false.
      nmest = nmsh
      DO im = 1, ninter-1
         IF (MAX(r4(im), r4(im+1)) > fatt_r1r3) THEN
            IF (.NOT. add) THEN
               irefin(im) = nptcond
               nmest = nmest + nptcond - 1
            END IF
            irefin(im+1) = nptcond
            nmest = nmest + nptcond - 1
            add = .true.
         ELSE
            irefin(im) = 1
            irefin(im+1) = 1
            nmest = nmest - 1
            add = .false.
         END IF
      END DO
!
      IF (nmest > nmax) THEN
         GO TO 360
      END IF
!
!     It appears that we can perform the desired selective mesh
!     refinement.
!
!     Now begin running through the mesh, adding and possibly deleting
!     points as indicated by the irefin array.
!
!     The integer new is a count of the number of intervals in
!     the tentative mesh being generated by the refinement strategy.
!
      new = 1
!
!     The first interval is treated as a special case, since xx(1)
!     always remains in the mesh, and cannot be a fixed point.
!
      rlen = xxold(2) - xx(1)
      slen = rlen
      IF (irefin(1) > 1) THEN
         dx = rlen/irefin(1)
         DO j = 2, irefin(1)
            new = new + 1
            xx(new) = xx(1) + (j-1)*dx
         END DO
      END IF
!
!     The fixed points specified by the fixpnt array cannot be
!     removed from the mesh.  The value fxnext indicates the 'next'
!     fixed point.  When no further fixed points remain to be processed
!     (or if nfxpnt = 0), fxnext is set to a value strictly larger than
!     the last mesh point, so that no mesh point can equal fxnext.
!     This way we need to compare only the values of xxold(i)
!     and fxnext.
!
      ifxcnt = 1
      IF (nfxpnt == 0) THEN
         fxnext = onep1*ABS(xxold(nmsh))
      ELSE
         fxnext = fixpnt(ifxcnt)
      END IF
!
!     jtkout is a counter of the number of consecutive points that
!     have been removed from the mesh.
!
      jtkout = 0
!
      DO im = 2, ninter
! 
         rlold = rlen
         rlen = xxold(im+1) - xxold(im)
! 
!        If xxold(im) is the next fixed point, it cannot be removed
!        and so we don't test its error estimates.
!
         IF (xxold(im) == fxnext) THEN
! 
            ifxcnt = ifxcnt + 1
            IF (ifxcnt > nfxpnt) THEN
               fxnext = onep1*ABS(xxold(nmsh))
            ELSE
               fxnext = fixpnt(ifxcnt)
            END IF
!
         ELSE IF (irefin(im) == 1) THEN
!
!           If xxold(im) is not a fixed point and irefin(im) = 1, we may wish
!           to remove point im from the mesh.
!
!           If we are considering removing points and jtkout = 0, this
!           is the first point in a possible consecutive set to be removed,
!           and we initialize phihat, which represents a maximum of
!           certain estimates.
!           If jtkout is not zero, previous points contiguous to this
!           point have been removed, and phihat retains its previous value.
! 
            slen = slen + rlen
!
            IF (jtkout < 1 .AND. r4(im) <= fatt_r3) THEN
! 
!              Increment the counter of removed points.
!              'Remove' the mesh point xxold(im) by not including it.
!
               jtkout = jtkout+1
               CYCLE
!
            END IF
!
!        end of logic for irefin(im) = 1.
         END IF
!
         IF (new+irefin(im) > nmax) THEN
! 
!           If the new mesh contains too many points, branch out of the
!           loop to try alternative strategies.
!
            GO TO 360
!
         END IF
!
         jtkout = 0
         new = new + 1
         xx(new) = xxold(im)
         IF (irefin(im) > 1) THEN
            dx = rlen/irefin(im)
            DO j = 2, irefin(im)
               new = new + 1
               xx(new) = xxold(im) + (j-1)*dx
            END DO
         END IF
         slen = rlen
!
         IF (new > nmax) THEN
!
!           If the new mesh contains too many points, branch out of the
!           loop to try alternative strategies.
!
            GO TO 360
! 
         ELSE IF (new > 3*ninter)  THEN
! 
!           Here, the new mesh does not exceed the specified maximum,
!           but has more than 3 times as many intervals as the old mesh.
!           Try doubling the mesh if possible.
!
            nmsh = nmold
            CALL DCOPY(nmsh, xxold, 1, xx, 1)
            CALL DBLMSH(nmsh, nmax, xx, nmold, xxold, maxmsh)
            ddouble = .true.
            RETURN
! 
         END IF
!
      END DO
!
!     To end up here, we have processed the entire interval,
!     and have neither exceeded the specified maximum nor
!     exceeded three times the number of intervals in the old
!     mesh.  The last mesh point remains unchanged.
!
      new = new + 1
      xx(new) = xxold(nmsh)
      nmsh = new
      maxmsh = .false.
      IF (iprint >= 0) WRITE(6,905) nmsh
      RETURN
!
360   CONTINUE
!
!     To reach here, the number of mesh points created at some stage
!     of the refinement process was larger than the maximum permitted
!     value nmax.
!
      nmsh = nmold
      CALL DCOPY(nmsh, xxold, 1, xx, 1)
      maxmsh = .true.
!
901   FORMAT(1H ,'selcondmsh.  nmsh, ipow =',2I5)
902   FORMAT(1H ,'im, jcomp, ermeas, normalized er',2I5,2(1PE11.3))
903   FORMAT(1H ,'errmax',1PE11.3)
904   FORMAT(1H ,'nmest, irefin',(10I5))
905   FORMAT(1H ,'selcondmsh.  new nmsh =',i8)
910   FORMAT(1H ,'ihcomp',(10I5))
!
      RETURN
   END SUBROUTINE SELCONDMSH
!!!
!!!
   SUBROUTINE SMPSELCONDMSH(ncomp, nmsh,  &
                            nfxpnt, fixpnt,  nmax, xx,  irefin, intref, numadd,  &
                            nmold, xxold, ddouble, maxmsh, r4, amg)
!
      USE BVPShared, ONLY: pdebug, use_c, comp_c, uval0,  &
                           nminit, iprint, idum,  &
                           DCOPY
!
      IMPLICIT NONE
!
      INTEGER(i4), INTENT(IN)    :: ncomp
      INTEGER(i4), INTENT(INOUT) :: nmsh
      INTEGER(i4), INTENT(IN)    :: nfxpnt
      REAL(r8),    INTENT(IN)    :: fixpnt(*)
      INTEGER(i4), INTENT(INOUT) :: nmax
      REAL(r8),    INTENT(INOUT) :: xx(*)
      INTEGER(i4), INTENT(INOUT) :: irefin(*)
      INTEGER(i4), INTENT(INOUT) :: intref
      INTEGER(i4), INTENT(INOUT) :: numadd
      INTEGER(i4), INTENT(INOUT) :: nmold
      REAL(r8),    INTENT(INOUT) :: xxold(*)
      LOGICAL,     INTENT(INOUT) :: ddouble
      LOGICAL,     INTENT(INOUT) :: maxmsh
      REAL(r8),    INTENT(INOUT) :: r4(*)
      REAL(r8),    INTENT(INOUT) :: amg(*)
!
      REAL(r8), PARAMETER :: zero = 0.0D+0
      REAL(r8), PARAMETER :: one = 1.0D+0
      REAL(r8), PARAMETER :: onep1 = 1.1D+0
      REAL(r8), PARAMETER :: erdcid = 5.0D+0
!
      REAL(r8) :: thres, r1, r2, r3, fatt_r1r3, fatt_r3, rlen, slen, dx
      REAL(r8) :: fxnext, rlold
      LOGICAL :: add
      INTEGER(i4) :: ninter, ithres, nptcond, nmest, im, new, j
      INTEGER(i4) :: ifxcnt, jtkout
!
!     The routine smpselcondmsh performs selective mesh refinement, by adding
!     points to one or three interval(s) in the region indicated
!     by the integer intref (numadd gives the trial number of points to
!     added in each  interval)  and by adding and removing point
!     using  the monitor function based on the conditioning parameters.
!
      maxmsh = .false.
!
      IF (pdebug) WRITE(6,901) nmsh
!
      ddouble = .false.
      nmold = nmsh
      ninter = nmsh - 1
!
!     Copy the current mesh into the xxold array.
!
      CALL DCOPY(nmold, xx, 1, xxold, 1)
!
      ithres = 0
      thres = one
!
!     On input, the array amg represents the conditioning vector defined
!     over the components and mesh intervals (not mesh points).
!
      CALL MONCONDMSH(nmsh,xx,r1,r2,r3,fatt_r1r3,fatt_r3,nptcond,r4,amg)
!
!     For each interval im,  the integer irefin(im) is calculated
!     based on an equidistrbution procedure involving the
!     threshold value thres.  If irefin(im) > 1, we add
!     points to interval im.  If irefin(im) = 1, we check whether
!     point im can be removed from the mesh.
!
!     nmest is a lower bound on the number of points in the new mesh.
!     We do not know in advance how many points will be removed,
!     so nmest is computed by assuming all eligible points are removed.
!
      nmest = nmsh
!
      DO im = 1, ninter
         IF (r4(im) > fatt_r1r3) THEN
            irefin(im) = nptcond
            nmest = nmest + nptcond - 1
         ELSE
            irefin(im) = 1
            nmest = nmest - 1
         END IF
      END DO
!
      IF (numadd > 49) THEN
         numadd = 49
      ELSE IF (numadd < 4) THEN
         numadd = 4
      END IF
!
      IF (intref == 1) THEN
         irefin(1) = MAX(numadd+1,irefin(1))
         nmest = nmest + numadd - 1
      ELSE IF (intref == ninter) THEN
         irefin(ninter) = MAX(numadd+1,irefin(ninter))
         nmest = nmest + numadd - 1
      ELSE
         IF (numadd > 9) THEN
            numadd = 9
         ELSE IF (numadd < 4) THEN
            numadd = 4
         END IF
         irefin(intref-1) = MAX(numadd+1 , irefin(intref-1))
         irefin(intref)   = MAX(numadd+1 , irefin(intref))
         irefin(intref+1) = MAX(numadd+1 , irefin(intref+1))
         nmest = nmest + 3*numadd - 1
      END IF
!
      IF (nmest > nmax) THEN
         GO TO 360
      END IF
!
!     It appears that we can perform the desired selective mesh
!     refinement.
!
!     Now begin running through the mesh, adding and possibly deleting
!     points as indicated by the irefin array.
!
!     The integer new is a count of the number of intervals in
!     the tentative mesh being generated by the refinement strategy.
!
      new = 1
!
!     The first interval is treated as a special case, since xx(1)
!     always remains in the mesh, and cannot be a fixed point.
!
      rlen = xxold(2) - xx(1)
      slen = rlen
      IF (irefin(1) > 1) THEN
         dx = rlen / irefin(1)
         DO  j = 2, irefin(1)
            new = new + 1
            xx(NEW) = xx(1) + (j-1)*dx
         END DO
      END IF
!
!     The fixed points specified by the fixpnt array cannot be
!     removed from the mesh.  The value fxnext indicates the 'next'
!     fixed point.  When no further fixed points remain to be processed
!     (or if nfxpnt = 0), fxnext is set to a value strictly larger than
!     the last mesh point, so that no mesh point can equal fxnext.
!     This way we need to compare only the values of xxold(i)
!     and fxnext.
!
      ifxcnt = 1
      IF (nfxpnt == 0) THEN
         fxnext = onep1*ABS(xxold(nmsh))
      ELSE
         fxnext = fixpnt(ifxcnt)
      END IF
!
!     jtkout is a counter of the number of consecutive points that
!     have been removed from the mesh.
!
      jtkout = 0
!
      DO im = 2, ninter
!  
         rlold = rlen
         rlen = xxold(im+1) - xxold(im)
! 
!        If xxold(im) is the next fixed point, it cannot be removed
!        and so we don't test its error estimates.
! 
         IF (xxold(im) == fxnext) THEN
!    
            ifxcnt = ifxcnt + 1
            IF (ifxcnt > nfxpnt) THEN
               fxnext = onep1*ABS(xxold(nmsh))
            ELSE
               fxnext = fixpnt(ifxcnt)
            END IF
! 
         ELSE IF (irefin(im) == 1) THEN
! 
!           If xxold(im) is not a fixed point and irefin(im) = 1, we may wish
!           to remove point im from the mesh.
! 
!           If we are considering removing points and jtkout = 0, this
!           is the first point in a possible consecutive set to be removed,
!           and we initialize phihat, which represents a maximum of
!           certain estimates.
!           If jtkout is not zero, previous points contiguous to this
!           point have been removed, and phihat retains its previous value.
! 
            slen = slen + rlen
! 
            IF (jtkout < 1 .AND. r4(im) <= fatt_r3) THEN
! 
!              Increment the counter of removed points.
!              'Remove' the mesh point xxold(im) by not including it.
! 
               jtkout = jtkout+1
               CYCLE
            END IF
!        end of logic for irefin(im) = 1
         END IF
!
         IF (new + irefin(im) > nmax) THEN
!    
!           If the new mesh contains too many points, branch out of the
!           loop to try alternative strategies.
! 
            GO TO 360
         END IF
!
         jtkout = 0
         new = new + 1
         xx(new) = xxold(im)
         IF (irefin(im) > 1) THEN
            dx = rlen/irefin(im)
            DO j = 2, irefin(im)
               new = new + 1
               xx(new) = xxold(im) + (j-1)*dx
            END DO
         END IF
         slen = rlen
! 
         IF (new > nmax) THEN
! 
!           If the new mesh contains too many points, branch out of the
!           loop to try alternative strategies.
! 
            GO TO 360
! 
         ELSE IF (NEW > 3*ninter)  THEN
! 
!           Here, the new mesh does not exceed the specified maximum,
!           but has more than 3 times as many intervals as the old mesh.
!           Try doubling the mesh if possible.
!
            IF (iprint == 1) WRITE(6,*) 'smpselcondmsh'
            nmsh = nmold
            CALL DCOPY(nmsh, xxold, 1, xx, 1)
            CALL DBLMSH(nmsh, nmax, xx, nmold, xxold, maxmsh)
            ddouble = .true.
            RETURN
! 
         END IF
!
      END DO
!
!     To end up here, we have processed the entire interval,
!     and have neither exceeded the specified maximum nor
!     exceeded three times the number of intervals in the old
!     mesh.  The last mesh point remains unchanged.
!
      new = new + 1
      xx(new) = xxold(nmsh)
      nmsh = new
      maxmsh = .false.
      IF (iprint >= 0) WRITE(6,905) nmsh
      RETURN
!
360   CONTINUE
!
!     To reach here, the number of mesh points created at some stage
!     of the refinement process was larger than the maximum permitted
!     value nmax.
!
      nmsh = nmold
      CALL DCOPY(nmsh, xxold, 1, xx, 1)
      maxmsh = .true.
!
901   FORMAT(1H ,'smpselcondmsh.  nmsh, ipow =',2I5)
902   FORMAT(1H ,'im, jcomp, ermeas, normalized er',2I5,2(1PE11.3))
903   FORMAT(1H ,'errmax',1PE11.3)
904   FORMAT(1H ,'nmest, irefin',(10I5))
905   FORMAT(1H ,'smpselcondmsh.  new nmsh =',i8)
910   FORMAT(1H ,'ihcomp',(10I5))
!
      RETURN
   END SUBROUTINE SMPSELCONDMSH
!!!
!!!
   SUBROUTINE MONCONDMSH(nmsh,xx,r1,r2,r3,fatt_r1r3,fatt_r3, nptcond,r4,amg)
!
      USE BVPShared, ONLY: pdebug, use_c, comp_c, uval0,  &
                           nminit, iprint, idum
!
      IMPLICIT NONE
!
      INTEGER(i4), INTENT(INOUT) :: nmsh
      REAL(r8),    INTENT(INOUT) :: xx(*)
      REAL(r8),    INTENT(INOUT) :: r1
      REAL(r8),    INTENT(INOUT) :: r2
      REAL(r8),    INTENT(INOUT) :: r3
      REAL(r8),    INTENT(INOUT) :: fatt_r1r3
      REAL(r8),    INTENT(INOUT) :: fatt_r3
      INTEGER(i4), INTENT(INOUT) :: nptcond
      REAL(r8),    INTENT(INOUT) :: r4(*)
      REAL(r8),    INTENT(INOUT) :: amg(*)
!
      REAL(r8), PARAMETER :: zero = 0.0D+0
      REAL(r8), PARAMETER :: one = 1.0D+0
!
      REAL(r8) :: r1m
      INTEGER(i4) :: i, nptm, nptr
!
!     the function moncond compute the monitor function based on the
!     conditioning parameters and the factor used to perform the mesh selection
!
      DO i = 1, nmsh-1
         r4(i) = (xx(i+1)-xx(i))*DABS(amg(i+1)- amg(i))
      END DO
!
      r2 = r4(1)
      DO i = 2, nmsh-1
         r2 = r2 + r4(i)
      END DO
      DO i = 1, nmsh-1
         r4(i) = r4(i)+(xx(i+1)-xx(i))*(r2/(xx(nmsh)-xx(1)))*0.08D0
      END DO
      r1 = r4(1)
      DO i = 2, nmsh-1
         r1 = MAX(r1,r4(i))
      END DO
      DO i = 1, nmsh-1
         r4(i) = r4(i)/r1
      END DO
      r1 = one
!
      r2 = r4(1)
      r1m = r4(1)
      DO i = 2, nmsh-1
         r1m = MIN(r1m,r4(i))
         r2 = r2 + r4(i)
      END DO
!
      r3 = r2/(nmsh-1)
      fatt_r3  = r3*0.00001D0
      fatt_r1r3= MAX(r3,r1*0.5D0)
!
      nptm = 0
      nptr = 0
      DO i = 1, nmsh-1
         IF (r4(i) >= fatt_r1r3)  nptm = nptm + 1
         IF (r4(i) <= fatt_r3)  nptr = nptr + 1
      END DO
!
      IF (nptm <= 1) THEN
         nptcond =  14
      ELSE IF (nptm <= 2) THEN
         nptcond =  10
      ELSE IF (nptm <= 4) THEN
         nptcond =  8
      ELSE IF (nptm <= 8) THEN
         nptcond = 6
      ELSE IF (nptm <= nmsh/20 ) THEN
         nptcond = 4
      ELSE
         nptcond = 2
      END IF
!
      IF (iprint == 1) WRITE(6,901)r1,r3,fatt_r1r3,nptcond,nptm,nptr
!
901   FORMAT(1H ,'moncondmsh.', (1PE11.3), 2(1PE11.3), 3I10)
!
      RETURN
   END SUBROUTINE MONCONDMSH
!
END
