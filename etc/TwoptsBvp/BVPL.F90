MODULE BVPL
!
   USE Base,      ONLY: i4, r4, r8
   USE BVPShared, ONLY: flmin, flmax, epsmch,  &
                        pdebug, use_c, comp_c,  &
                        uval0, nminit, iprint, idum,  &
                        a21a, a22a, a23a, a24a, a31a, a32a, a33a, a34a,  &
                        c1a, c2a, c16a, c26a, c123a, c223a, c14a, c24a,  &
                        a21b, a22b, a23b, a24b, a25b, a31b, a32b, a34b,  &
                        a35b, a41b, a42b, a43b, a44b, a45b,  &
                        b1b, b2b, b3b, c1b, c2b, c3b, c16b, c26b, c36b,  &
                        c123b, c223b, c323b,c14b,c24b,c34b,  &
                        ifinal, iback, iprec,  &
                        ABDNRM, DASUM, DONEST,  &
                        SPRT, MPRT, COLROW1, INVERSE, CRDCMP1, CRSLVE,  &
                        LUFAC, LUSOL, DCOPY, DAXPY, DDOT,  &
                        DSCAL, DSWAP, IDAMAX, DLOAD,  &
                        MAXPY, MATCOP, MTLOAD, MSSQ, DSSQ
   USE BVPExtern, ONLY: INITU, FSUB, DFSUB, GSUB, DGSUB
!
   IMPLICIT NONE
!
   PRIVATE
   PUBLIC :: TWPBVPL_INIT, TWPBVPL
!
CONTAINS
!
   SUBROUTINE TWPBVPL_INIT
!
      USE BVPShared, ONLY: pdebug, use_c, comp_c,  &
                           uval0, nminit, iprint, idum
!
      IMPLICIT NONE
!
!     BLOCK DATA
!
      nminit = 10
      pdebug = .true.
      iprint = 1
      uval0 = 0.0_r8
      use_c = .false.
      comp_c = .false.
!
      RETURN
   END SUBROUTINE TWPBVPL_INIT
!
!-------------------------------------------------------------------------------
!
!  SUBROUTINE TWPBVPL
!
!  Code converted using TO_F90 by Alan Miller
!  Date: 2017-02-06  Time: 18:03:10
!
!  The subroutine TWPBVPLC is intended to solve two-point boundary value problems
!
!  References:
!     Cash, J. R., and Mazzia, F., 2006:
!        Hybrid mesh selection algorithms based on conditiong for two-point
!        boundary value problems. J. Numer. Analy. Industr. Appl. Math.,
!        1, 1, 81--90.
!
!  Revision History:
!
!  10JUL2006: Added rpar and ipar in the functions
!             DFSUB, DGSUB, FSUB, GSUB
!             changed the name of the variable double in ddouble
!
!  31AUG2004: This is a modified version of twpbvp that uses the conditioning
!             in the mesh selection.
!
!  New subroutines not included in the old version:
!     CONDESTIM,             MONCOND, SELCONDMSH, SELCONDERRMSH, SMPSELCONDMSH
!
!  Updates subroutines:
!     BVPSOL, CONV4, FAIL4, FAIL6, CONV8, FAIL8, NEWTEQ, MSHREF, WTCHDG
!     BVPSOL,        FAIL4, FAIL6, DECID4,       NEWTEQ, MSHREF
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
   SUBROUTINE TWPBVPL(ncomp, nlbc, aleft, aright,  &
                      nfxpnt, fixpnt, ntol, ltol, tol,  &
                      linear, givmsh, giveu, nmsh,  &
                      nxxdim, xx, nudim, u, nmax,  &
                      lwrkfl, wrk, lwrkin, iwrk,  &
                      ckappa1, gamma1, ckappa, rpar, ipar, iflbvp,  &
                      liseries, iseries, indnms)
!
      IMPLICIT NONE
!
      INTEGER(i4), INTENT(IN)    :: ncomp
      INTEGER(i4), INTENT(IN)    :: nlbc
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
      INTEGER(i4), INTENT(IN)    :: liseries
      INTEGER(i4), INTENT(INOUT) :: iseries(*)
      INTEGER(i4), INTENT(INOUT) :: indnms
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
      INTEGER(i4) :: iipvlu, lipvlu, iisign, lisign, ir4, lr4
      INTEGER(i4) :: i
!
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
      isp = lwrkfl - 3 - 2*ntol - 22*ncomp - 6*ncomp*ncomp
      iden = 6*ncomp*ncomp + 22*ncomp + 3
      nmax1 = isp/iden
!
      isp = lwrkin - 2*ncomp-3
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
!     ldefex = ncomp*(nmax-1)
!
!     def6 uses the same space as defexp
!
      idef6 = idefex
      ldef6 = ncomp*(nmax-1)
!
      idefim = idef6 + ldef6
!     ldefim = ncomp*(nmax-1)
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
!     lxmer = ncomp*nmax
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
      ir4 = iwrkrhs + lwrkrhs
      lr4 = nmax
!
      ibhold = ir4 + lr4
      lbhold = 9*ncomp*ncomp
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
!     lisign = ncomp*nmax
!
      CALL BVPSOL(ncomp, nmsh, nlbc, aleft, aright,  &
                  nfxpnt, fixpnt, ntol, ltol, tol, nmax,   &
                  linear, giveu, givmsh, xx, nudim, u,  &
                  wrk(idefex), wrk(idefim), wrk(idef),  &
                  wrk(idelu), wrk(irhs), wrk(ifval),  &
                  wrk(itpblk), wrk(ibtblk), wrk(iajac), wrk(ibhold),  &
                  wrk(ichold), wrk(ibhold), iwrk(iipvbk), iwrk(iipvlu),  &
                  iwrk(iisign), wrk(iuint), wrk(iftmp), wrk(itmrhs),  &
                  wrk(idftm1), wrk(idftm2), wrk(idgtm),  &
                  wrk(iutri), wrk(irhtri), wrk(ixmer),  &
                  wrk(ixxold), wrk(iuold), wrk(iusve),  &
                  wrk(itmp), wrk(idsq), wrk(idexr), wrk(irtdc),  &
                  wrk(irerr), wrk(ietst6), wrk(ietst8), wrk(iermx),  &
                  iwrk(iihcom), iwrk(iiref), wrk(idef6), wrk(idef8),  &
                  iflbvp,  &
                  wrk(iamg), wrk(ic1), wrk(iwrkrhs),  &
                  ckappa1, gamma1, ckappa, wrk(ir4), rpar, ipar,  &
                  liseries, iseries, indnms)
!
      RETURN
   END SUBROUTINE TWPBVPL
!!!
!!!
   SUBROUTINE BVPSOL(ncomp, nmsh, nlbc, aleft, aright,  &
                     nfxpnt, fixpnt, ntol, ltol, tol, nmax,  &
                     linear, giveu, givmsh, xx, nudim, u,  &
                     defexp, defimp, def,  &
                     delu, rhs, fval,  &
                     topblk, botblk, ajac, bhold,  &
                     chold, dhold, ipvblk, ipivlu, isign,  &
                     uint, ftmp, tmprhs,  &
                     dftmp1, dftmp2, dgtm,  &
                     utrial, rhstri, xmerit,  &
                     xxold, uold, usave,  &
                     tmp, dsq, dexr, ratdc,  &
                     rerr, etest6, etest8, ermx,  &
                     ihcomp, irefin, def6, def8,  &
                     iflbvp,  &
                     amg, c1, wrkrhs,  &
                     ckappa1, gamma1, ckappa, r4, rpar, ipar,  &
                     liseries, iseries, indnms)
!
      USE BVPShared, ONLY: flmin, flmax, epsmch,  &
                           pdebug, use_c, comp_c, uval0, nminit, iprint, idum
      USE BVPExtern, ONLY: FSUB, DFSUB, GSUB, DGSUB
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
      REAL(r8),    INTENT(INOUT) :: dhold(3*ncomp,*)
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
      REAL(r8),    INTENT(INOUT) :: r4(*)
      REAL(r8),    INTENT(INOUT) :: rpar(*)
      INTEGER(i4), INTENT(INOUT) :: ipar(*)
      INTEGER(i4), INTENT(IN)    :: liseries
      INTEGER(i4), INTENT(INOUT) :: iseries(*)
      INTEGER(i4), INTENT(INOUT) :: indnms
!
      REAL(r8),    PARAMETER :: zero = 0.0D+0
      REAL(r8),    PARAMETER :: one = 1.0D+0
      REAL(r8),    PARAMETER :: third = 0.33D+0
      REAL(r8),    PARAMETER :: fourth = 0.25D+0
      REAL(r8),    PARAMETER :: quan6 = 0.5D+0
      REAL(r8),    PARAMETER :: power = 1.0d+0/6.0d+0
      INTEGER(i4), PARAMETER :: itcondmax = 1
!
      REAL(r8)    :: df(ncomp,ncomp), smaldef, smalldef, bigdef, holdef
      REAL(r8)    :: siz, rat, drat
      REAL(r8)    :: gamma1old, ckappa1old, tolmin, rnsq, err6
      REAL(r8)    :: fxfct = 10.0D+0
      LOGICAL     :: ddouble, chstif, nodouble, forcedouble, reposs
      LOGICAL     :: smooth, succes, strctr, trst6, reaft6
      LOGICAL     :: onto6, onto8, ludone, rhsgiv
      LOGICAL     :: first4, first8
      LOGICAL     :: frscal = .TRUE.
      LOGICAL     :: stab_kappa, stab_gamma, stab_cond, stiff_cond
      LOGICAL, SAVE :: mchset = .TRUE.
      LOGICAL     :: maxmsh = .FALSE.
      INTEGER(i4) :: itcond = 0
      INTEGER(i4) :: nmold, numbig, nummed, iorder, iflnwt, itnwt
      INTEGER(i4) :: ninter, n, iter
      INTEGER(i4) :: i, j, jflag, indnmsold, icmph, ix, iv, iu, ipoint
      INTEGER(i4) :: ixx, intol, numadd, jc
!!!
!
!     frscal = .true.
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
      CALL STCON1
      CALL STCON2
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
      DO i = 1, ncomp
         DO j = 1, ncomp
            df(i,j) = 0.0d+0
         END DO
      END DO
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
      maxmsh = .false.
      Chstif = .true.
      ddouble = .false.
!
!     Initialize parameter for the conditioning estimation
!
      IF (comp_c) THEN
      gamma1old  = flmax
      gamma1     = flmax
      ckappa1old  = flmax
      ckappa1     = flmax
      ckappa      = flmax
      stiff_cond = .false.
      stab_cond  = .false.
      END IF
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
      IF (.NOT. giveu) CALL INITU(ncomp, nmsh, xx, nudim, u, rpar, ipar)
      indnms = 0
      indnmsold = 0
!
!     top of logic for 4th order solution ****
!
400   CONTINUE
!
      IF (indnmsold /= nmsh) THEN
         indnms = indnms + 1
         iseries(indnms) = nmsh
         indnmsold = nmsh
      END IF
!
      IF (indnms >= liseries) THEN
         WRITE(6,1008) nmsh
         GO TO 1900
      END IF
      IF (maxmsh) GO TO 900
!
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
                     itnwt, iflnwt, isign, rpar, ipar, frscal)
!
      END IF   ! if (linear) then
!
!     these flags are used in the mesh selection strategy
!
      IF (iflnwt == 0) THEN
!
!        COMPUTE ESTIMATIONS OF CONDITIONING NUMBERS: norms of inverse
!        jacobian matrix BY BRUGNANO & TRIGIANTE, AND HIGHAM
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
!              WRITE(6,*) 'amg', (amg(i), i = 1, nmsh)
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
            stiff_cond = (( (ckappa1/gamma1 >= 1d2  )))
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
!        The subroutine DFEXCL substitute CONV4
!
         CALL DFEXCL(ncomp, nmsh, xx, nudim, u, def8, def, linear, fval,  &
                     tmp, df, ipivlu, dhold,  &
                     ntol, ltol, tol, jflag, rpar, ipar)
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
            stiff_cond = (( (ckappa1/gamma1 >= 1.0d1  )))
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
!  note: ratdc in subroutines fail4, fail6 does not related to 
!        ratdc=dfexmx/defimp , they only use the storage
         GO TO 400
!
      END IF  !   if (iflnwt == 0) then
!
!     IF (succes) THEN
!        IF (comp_c) THEN
!           CALL ESTIMKAPPA(nmsh, ncomp, n, xx, topblk, nlbc, ncomp, ajac,  &
!                           ncomp, 2*ncomp, ninter, botblk, ncomp-nlbc, ipvblk,  &
!                           isign, c1, wrkrhs, ckappa)
!        END IF
!        IF (iprint /= -1 .AND. comp_c) THEN
!           IF (ckappa >= tolmin/epsmch) WRITE(6,1005)
!        END IF
!        iflbvp = 0
!        RETURN
!     ELSE IF (maxmsh) THEN
!
      IF (jflag == 1) THEN
!
         nodouble = ( (stiff_cond .AND. .NOT. stab_cond) .AND. (use_c))
!        nodouble = ((iorder.eq.4) .and.   &
!                   (stiff_cond .and. .not. stab_cond) .and. (use_c))
!
         forcedouble = .false.
!
         WRITE(6,*) 'nodouble', itcond
         IF (use_c) THEN
            IF (itcond == itcondmax) THEN
               itcond = 0
               forcedouble = .true.
            ELSE
               itcond = itcond + 1
               forcedouble = .false.
            END IF
         END IF
         IF (nodouble .AND. .NOT. forcedouble) THEN
            CALL SELCONDMSH(ncomp, nmsh, nfxpnt, fixpnt, nmax, xx, irefin,  &
                            nmold, xxold, ddouble, maxmsh, r4, amg)
            ddouble = .false.
         ELSE
            CALL DBLMSH(nmsh, nmax, xx, nmold, xxold, maxmsh)
            itcond = 0
         END IF
         CALL MATCOP(nudim, ncomp, ncomp, nmold, u, uold)
         Call INTERP(ncomp, nmsh, xx, nudim, u, nmold, xxold, uold)
         IF (maxmsh) GO TO 900
         GO TO 400
!
      ELSE
!
!        find where biggest deferred correction is
!        If(Chstif) Then
!
         smaldef=1.0D+40
         bigdef=0.0D+0
         icmph = 1
         ix = 1
         DO iv = 1, nmsh-1
            DO iu = 1, ntol
               ipoint = ltol(iu)
               holdef = ABS(def8(ipoint,iv))
               IF (smaldef > holdef) smalldef = holdef
               IF (holdef > bigdef) THEN
                  bigdef = holdef
                  icmph = ipoint
                  ixx = iv
                  intol = iu
               END IF
            END DO
         END DO
!
!        Biggest deferred correction is in component icmph and
!        at the mesh interval ix.
!        Now compute an explicit deferred correction for this.
!
         CALL EXPL(ncomp, nmsh, xx, nudim, u, dgtm, fval, ixx, rpar, ipar)
!
         ix = ixx
         siz = ABS(dgtm(icmph))
!
!        write(6,*) 'siz=',siz
         rat = siz / bigdef
!
!        write(6,*) 'rat=',rat
!        write(6,*) ' bigdef=',bigdef,' > ',dsqrt(tol(Icmph))
!        write(6,*) ' dsq=', dsqrt(tol(Icmph)) 
!
         IF (rat > 50.0D+0 .AND. bigdef > DSQRT(tol(Icmph)) .AND.  &
             siz > 1.0D+0) THEN
!
            IF (use_c) THEN
               IF ((stiff_cond)) THEN
                  drat = bigdef / (MAX(one, ABS(u(icmph,ix)))*tol(intol))
                  IF (pdebug) WRITE(6,913) drat, u(icmph,ix), tol(intol)
!                 numadd = drat**power
                  numadd = 15
                  CALL SMPSELCONDMSH(ncomp, nmsh, nfxpnt, fixpnt,  nmax, xx, irefin, ix, numadd,  &
                                     nmold, xxold, ddouble, maxmsh, r4,amg)
               ELSE
                  drat = bigdef / (MAX(one, ABS(u(icmph,Ix)))*tol(intol))
                  IF (pdebug) WRITE(6,913) drat, u(icmph,ix), tol(intol)
!                 numadd = drat**power
                  numadd = 15
                  CALL SMPMSH(nmsh, nmax, xx, ix, numadd, nmold, xxold, maxmsh)
               END IF
            ELSE
               drat = bigdef / (MAX(one, ABS(u(icmph,ix)))*tol(intol))
               IF (pdebug) WRITE(6,913) drat, u(icmph,ix), tol(intol)
               numadd = 15
               CALL SMPMSH(nmsh, nmax, xx, ix, numadd, nmold, xxold, maxmsh)
            END IF 
913         FORMAT(1H ,'drat,u,tol',3(1pe11.3))
!
            GO TO 400
!
         ELSE
!
!           Chstif = .false.
            onto6 = .true.
            IF (linear .AND. ddouble) reposs = .true.
!
         END IF
!
!        Endif
! endif of Chstif=true

!       Iorder = 6
!       if (iflnwt.ne.0) then
!         call mshref(ncomp, nmsh, nlbc, ntol, ltol,
!     *             iorder, rhs, ratdc,
!     *             nmax, xx, nmold, xxold, ddouble, maxmsh,
!     *             numbig, nummed,
!     *             amg,stab_cond,stiff_cond,
!     *             r4, nfxpnt,fixpnt, irefin,itcond,itcondmax)
!       elseif((.not.onto6)) then
!         Call Matcop(Nudim, Ncomp, Ncomp, Nmold, U, Uold)
!         Call Interp(Ncomp, Nmsh, Xx, Nudim, U, Nmold, Xxold, Uold)
!         goto 400
!       endif
!
      END IF
408   CONTINUE
!
      CALL MATCOP(ncomp, ncomp, ncomp, nmsh-1, def, def6)
!
!        if (succes) then
!            if (iprint .ne. -1 .and. comp_c ) then
!              if (ckappa .ge. tolmin/epsmch) write(6,1005)  
!            end if
!            iflbvp = 0
!            return
!
      IF (maxmsh) THEN
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
                    iflnwt, rpar, ipar)
!
      ELSE
!
         CALL FIXJAC(ncomp, nmsh, nlbc, iorder, ntol, ltol, tol, xx,  &
                     nudim, u, def, def, delu, rhs, fval, utrial, rhstri,  &
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
                        iter, iflnwt, isign, rpar, ipar, frscal)
! 
         END IF
!
      END IF
!
      IF (iflnwt == 0) THEN
!
         itcond = 0
         CALL CONV6(ncomp, nmsh, ntol, ltol, tol, nudim, u, uold,  &
                    etest6, err6, trst6, onto8, reaft6, linear, succes)
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
      IF (maxmsh) GO TO 900
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
      IF (linear) CALL FNEVAL(ncomp, nmsh, xx, nudim, u, fval, rpar, ipar)
!
!     Calculate 8th order deferred corrections (the array def8).
!
      CALL DF8CAL(ncomp, nmsh, xx, nudim, u, fval, def8, linear,  &
                  tmp, df, ipivlu, dhold, ntol, ltol, tol, jc, rpar, ipar)
!
      IF (jc == 1) THEN
!
         nodouble = ((stiff_cond) .AND. .NOT. stab_cond  .AND. (use_c))
!      nodouble = nodouble
!     *  .or.((iorder.gt.4) .and. .not.stab_cond .and. (stiff_cond)
!     *   .and. (use_c))
         forcedouble = .false.
!
         IF (use_c) THEN
            IF (itcond == itcondmax) THEN
               itcond = 0
               forcedouble = .true.
            ELSE
               itcond = itcond + 1
               forcedouble = .false.
            END IF
         END IF
         IF (nodouble .AND. .NOT. forcedouble) THEN
            CALL SELCONDMSH(ncomp, nmsh, nfxpnt, fixpnt, nmax, xx,  irefin,  &
                            nmold, xxold, ddouble, maxmsh, r4, amg)
            ddouble = .false.
         ELSE
            CALL DBLMSH(nmsh, nmax, xx, nmold, xxold, maxmsh)
            itcond = 0
         END IF
         IF (maxmsh) GO TO 900
         CALL MATCOP(nudim, ncomp, ncomp, nmold, uold, u)
         CALL INTERP(ncomp, nmsh, xx, nudim, u, nmold, xxold, uold)
         GO TO 400
      END IF
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
                    iflnwt, rpar, ipar)
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
                        iter, iflnwt, isign, rpar, ipar, frscal)
!
         END IF
!
      END IF
!
      IF (iflnwt == 0) THEN
!
         itcond = 0
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
      IF (comp_c) THEN
         CALL ESTIMKAPPA(nmsh, ncomp, n, xx, topblk, nlbc,ncomp, ajac, ncomp, 2*ncomp,  &
                         ninter, botblk, ncomp-nlbc, ipvblk, isign, c1, wrkrhs, ckappa)
         IF (linear .AND. ckappa >= tolmin/epsmch) WRITE(6,1006)
         IF (.NOT.linear .AND. ckappa >= tolmin/epsmch) WRITE(6,1007)
      END IF
!
      RETURN
1900  CONTINUE
!
!     Error exit---too many meshes  .return
!
      iflbvp = 1
      WRITE(6,*) 'Terminated, too many meshes'
      IF (comp_c) THEN
         CALL ESTIMKAPPA(nmsh, ncomp, n, xx, topblk, nlbc, ncomp, ajac, ncomp, 2*ncomp,  &
                         ninter, botblk, ncomp-nlbc, ipvblk, isign, c1, wrkrhs, ckappa)
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
1008  FORMAT(1H ,'Terminated too many meshes, nmsh',i5)
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
                    nudim, u, uold, etest6, err6, trst6, onto8, reaft6, linear, succes)
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
      LOGICAL,     INTENT(IN)    :: linear
      LOGICAL,     INTENT(INOUT) :: succes
!
      REAL(r8), PARAMETER :: quan6 = 0.5d+0, third = 0.33d+0, one = 1.d+0
      LOGICAL :: errok
      INTEGER(i4) :: i
!!!
!
!     The Newton iteration converged for a 6th-order solution.
!
      IF (.NOT. linear) THEN
         CALL DLOAD(ntol, one, etest6, 1)
      ELSE
         DO i = 1, ntol
            etest6(i) = one/MAX(quan6, tol(i)**third)
         END DO
      END IF
!
      CALL DLOAD(ntol, 10D0, etest6, 1)
      IF (iprint == 1) WRITE(6,901)
!
      succes = .false.
!
!     The logical flag reaft6 is true only when the 6th order solution
!     failed.  Since the 6th order solution converged, reaft6 is false.
!
!     reaft6 = .false.
      onto8 = .true.
!
!     Calculate the error estimates for the 4th order solution.
!
      CALL ERREST(ncomp, nmsh, ntol, ltol, tol,  &
                  nudim, u, uold, etest6, err6, errok)
!
      IF (trst6 .AND. errok) THEN
!      IF (errok) THEN
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
!     no possibility of richardson extrapolation error test
!     reaft6 = .false.
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
         IF (iprint == 1) WRITE(6,9994)
9994     FORMAT(1H ,'***in fail6')
         IF (iprint == 1) WRITE(6,9993) remax, eight*tol(itlmx)
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
!
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
                                     nmold, xxold, ermx, ddouble, maxmsh, r4, amg, stab_cond)
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
      REAL(r8), PARAMETER :: quan8 = 0.5D+0
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
!     IF (first8) THEN
!        er6old = huge
!        er8old = huge
!        first8 = .false.
!     END IF
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
      CALL DLOAD(ntol, 10d0, etest8, 1)
!
      CALL ERREST(ncomp, nmsh, ntol, ltol, tol,  &
                  nudim, u, uold, etest8, err8, errok)
      IF (errok)  THEN
         succes = .true.
         RETURN
      END IF
!
!      write(6,*) ' err8', err8, nmsh, stiff_cond, stab_cond
!      if ( (use_c .and. err8 .le. 5e2 )) then
!        do im=1,ncomp
!         do jm = 1,nmsh
!              def8(im,jm)=abs(u(im,jm)-uold(im,jm))
!            def8(im,jm)=max(abs(def8(im,jm)),
!     *       abs(u(im,jm)-uold(im,jm)))
!           enddo
!        enddo
!
!      elseif ( .not. use_c .and. err8 .le. 5e2 ) then
!       do im=1,ncomp
!         do jm = 1,nmsh
!           def8(im,jm)=abs(u(im,jm)-uold(im,jm))
!           def8(im,jm)=max(abs(def8(im,jm)),abs(u(im,jm)-uold(im,jm)))
!           enddo
!         enddo
!      endif
!
!     At this point, the 8th order solution converged, but did not
!     satisfy the test for termination.
!
!     IF (pdebug) WRITE(6,902) err6, err8, er6old, er8old
!     IF (nmsh < nmold .AND. err6 > efact*er6old .AND.  &
!         err8 > efact*er8old) THEN
!
!        If the number of mesh points decreased and the errors in the
!        6th and 8th order solutions did not decrease sufficiently compared
!        to the previous mesh, double the mesh and go back to try to
!        calculate a 4th order solution.
!
!        CALL DBLMSH(nmsh, nmax, xx, nmold, xxold, maxmsh)
!        IF (.NOT. maxmsh) THEN
!           er6old = err6
!           er8old = err8
!           If the problem is not linear we use the old solution
!           instead the the first value
! old code
!              call INITU(ncomp, nmsh, xx, nudim, u,rpar,ipar)
! new code
!           IF (linear) THEN
!              CALL INITU(ncomp, nmsh, xx, nudim, u,rpar,ipar)
!           ELSE
!              we do not use u but uold
!              CALL MATCOP(nudim, ncomp, ncomp, nmold, u, uold)
!              CALL INTERP(ncomp, nmsh, xx, nudim, u, nmold, xxold, uold)
!           END IF
!        END IF
!        RETURN
!     END IF
!
!     Here, we know that
!     (1) the number of mesh points exceeds that for the previous mesh; or
!     (2) the number of mesh points decreased, the 6th and 8th order
!         errors did not satisfy the termination tests, but they did not
!         increase so much that the mesh needed to be doubled.
!
!     er6old = err6
!     er8old = err8
!
!     IF (err8 <= err6) THEN
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
         RETURN
!
!     ELSE
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
!     END IF
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
   SUBROUTINE DCCAL(ncomp, nmsh, ntol, ltol, defexp, dfctol, fval,  &
                    dfexmx, incmp, inmsh, intol, derivm)
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
      REAL(r8),    INTENT(IN)    :: dfctol
      REAL(r8),    INTENT(IN)    :: fval(ncomp,nmsh)
      REAL(r8),    INTENT(INOUT) :: dfexmx
      INTEGER(i4), INTENT(INOUT) :: incmp
      INTEGER(i4), INTENT(INOUT) :: inmsh
      INTEGER(i4), INTENT(INOUT) :: intol
      REAL(r8),    INTENT(INOUT) :: derivm
!
      REAL(r8), PARAMETER :: zero = 0.0D+0
      REAL(r8), PARAMETER :: one = 1.0D+0
      REAL(r8), PARAMETER :: rtst = 50.0D+0
      REAL(r8), PARAMETER :: tstrat = 0.1D+0
!
      REAL(r8) :: dval, smtest, texp, timp, abtexp, abrat
      REAL(r8) :: defimp(ncomp,nmsh-1), ratdc(nmsh-1)
      REAL(r8) :: rat1, rat2, dfimmx
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
      RETURN
!     do not want to calculate variables corresponding to implicit dc
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
   SUBROUTINE DECID4(linear, dfexmx,  &
                     derivm, dfold, tolval,  &
                     onto6, smooth, callrt, strctr, oscchk, ddouble , reposs)
!
      USE BVPShared, ONLY: pdebug, use_c, comp_c, uval0,  &
                           nminit, iprint, idum
!
      IMPLICIT NONE
!
      LOGICAL,  INTENT(IN)    :: linear
      REAL(r8), INTENT(INOUT) :: dfexmx
      REAL(r8), INTENT(INOUT) :: derivm
      REAL(r8), INTENT(INOUT) :: dfold
      REAL(r8), INTENT(IN)    :: tolval
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
!     IF (rat2 < rtst) THEN
!        IF (stest) THEN
!           smooth = .true.
!        ELSE
!           oscchk = .true.
!        END IF
!        RETURN
!     END IF
!
!     We know now that rat2 .ge. rtst.
!
      thttol = thrtwo*tolval
      IF (pdebug) WRITE(6,903) thttol
!
!     IF (rat1 < rtst .AND. dfexmx < thttol) THEN
!        IF (stest) THEN
!           smooth = .true.
!        ELSE
!           oscchk = .true.
!        END IF
!        RETURN
!     END IF
!
!     IF (rat1 < rtst .AND. dfexmx >= thttol) THEN
!        callrt = .true.
!        RETURN
!     END IF
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
      RETURN
!
!     IF (derivm > derval .AND. dfexmx > thttol) THEN
!        IF (dfimmx < one) THEN
!           callrt = .true.
!        ELSE
!           strctr = .true.
!           IF (linear) THEN
!              onto6 = .false.
!              IF (two*rat1 >= oldrt1) ddouble = .true.
!              end of logic for linear
!           END IF
!           end of logic for dfimmx .ge. one
!        END IF
!        RETURN
!        end of logic for derivm .gt. derval .and dfexmx .gt. thttol
!     END IF
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
!     IF (linear) reposs = .true.
!
901   FORMAT(1H ,'decid4')
902   FORMAT(1H ,'tolval, rtst',2(1PE11.3))
903   FORMAT(1H ,'thttol',1PE11.3)
!
      RETURN
   END SUBROUTINE DECID4
!!!
!!!
   SUBROUTINE DFEXCL(ncomp, nmsh, xx, nudim, u, def8, def6, linear, fval,  &
                     tmp, df, ip, dhold, ntol, ltol, tol, jc, rpar, ipar)
!
      USE BVPShared, ONLY: pdebug, use_c, comp_c, uval0,  &
                           nminit, iprint, idum, &
                           a21a, a22a, a23a, a24a, a31a, a32a, a33a, a34a, c1a, c2a,  &
                           c16a, c26a, c123a, c223a, c14a, c24a,  &
                           ifinal, iback, iprec,  &
                           LUFAC, LUSOL
!
      IMPLICIT NONE
!
      INTEGER(i4), INTENT(IN)    :: ncomp
      INTEGER(i4), INTENT(INOUT) :: nmsh
      REAL(r8),    INTENT(IN)    :: xx(nmsh)
      INTEGER(i4), INTENT(IN)    :: nudim
      REAL(r8),    INTENT(INOUT) :: u(nudim,nmsh)
      REAL(r8),    INTENT(INOUT) :: def8(ncomp,nmsh-1)
      REAL(r8),    INTENT(INOUT) :: def6(ncomp,nmsh-1)
      LOGICAL,     INTENT(IN)    :: linear
      REAL(r8),    INTENT(INOUT) :: fval(ncomp,nmsh)
      REAL(r8),    INTENT(INOUT) :: tmp(ncomp,8)
      REAL(r8),    INTENT(INOUT) :: df(ncomp,ncomp)
      INTEGER(i4), INTENT(INOUT) :: ip(2*ncomp)
      REAL(r8),    INTENT(INOUT) :: dhold(2*ncomp,2*ncomp)
      INTEGER(i4), INTENT(IN)    :: ntol
      INTEGER(i4), INTENT(IN)    :: ltol(ntol)
      REAL(r8),    INTENT(IN)    :: tol(ntol)
      INTEGER(i4), INTENT(INOUT) :: jc
      REAL(r8),    INTENT(INOUT) :: rpar(*)
      INTEGER(i4), INTENT(INOUT) :: ipar(*)
!
      REAL(r8), PARAMETER :: one = 1.0D+0
      REAL(r8), PARAMETER :: two = 2.0D+0
!
      REAL(r8) :: st1(200), st2(200), st3(200)
!
      REAL(r8) :: hmsh, c16h, c26h
      REAL(r8) :: fvim, fvim1, uim, uim1, xxc1, xxc2, dfij
      REAL(r8) :: tmp3, tmp4, er
      INTEGER(i4) :: im, ic, nit, i, j, ier, ii
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
         c16h = c16a/hmsh
         c26h = c26a/hmsh
         DO ic = 1, ncomp
            fvim = fval(ic,im)
            fvim1 = fval(ic,im+1)
            uim = u(ic,im)
            uim1 = u(ic,im+1)
            tmp(ic,3) = c16h*(uim1-uim)+c123a*fvim1+c14a*fvim
            tmp(ic,4) = c26h*(uim1-uim)+c223a*fvim1+c24a*fvim
!
!           Put cubic Hermite approximations to the unknowns in
!           tmp(ic,3) and tmp(ic,4).
!
            st1(ic) = (uim+uim1)/two
            st2(ic) = a21a*fvim + a24a*fvim1
            st3(ic) = a31a*fvim + a34a*fvim1
         END DO
!
         xxc1 = xx(im)+c1a*hmsh
         xxc2 = xx(im)+c2a*hmsh
!
         DO nit = 1, 10
!
            DO ic = 1, ncomp
               tmp3 = tmp(ic,3)
               tmp4 = tmp(ic,4)
               tmp(ic,1)  = st1(ic) + hmsh*(st2(ic) + a22a*tmp3 + a23a*tmp4)
               tmp(ic,2)  = st1(ic) + hmsh*(st3(ic) + a32a*tmp3 + a33a*tmp4)
            END DO
!
            CALL FSUB(ncomp,xxc1,tmp(1,1),tmp(1,5),rpar,ipar)
            CALL FSUB(ncomp,xxc2,tmp(1,2),tmp(1,6),rpar,ipar)
!
            CALL DFSUB(ncomp,xxc1,tmp(1,1),df(1,1),rpar,ipar)
            DO i = 1, ncomp
               tmp(i,5) = tmp(i,5)-tmp(i,3)
               tmp(i,6) = tmp(i,6)-tmp(i,4)
               DO j = 1,ncomp
                  dfij = hmsh*df(i,j)
                  dhold(i,j) = -a22a*dfij
                  dhold(i,j+ncomp) = -a23a*dfij
               END DO
            END DO
!
            CALL DFSUB(ncomp,xxc2,tmp(1,2),df,rpar,ipar)
            DO i = 1, ncomp
               DO j = 1, ncomp
                  dfij = hmsh*df(i,j)
                  dhold(i+ncomp,j) = -a32a*dfij
                  dhold(i+ncomp,j+ncomp) = -a33a*dfij
               END DO
            END DO
!
            DO i = 1, ncomp
               dhold(i,i) = dhold(i,i) + one
               dhold(i+ncomp,i+ncomp) = dhold(i+ncomp,i+ncomp) + one
            END DO
!
            CALL LUFAC(2*ncomp,2*ncomp,dhold,ip,ier)
            CALL LUSOL(2*ncomp,2*ncomp,dhold,ip,tmp(1,5),tmp(1,7))
!
            DO i = 1, ncomp
               tmp(i,3) = tmp(i,3) + tmp(i,7)
               tmp(i,4) = tmp(i,4) + tmp(i,8)
            END DO
!
            jc = 0
            IF (linear) GO TO 70
!
            DO i = 1, ntol
               ii = ltol(i)
               er = tol(i)/hmsh
               IF (ABS(tmp(ii,7)) > er*MAX(one,ABS(tmp(ii,3)))  .OR.  &
                   ABS(tmp(ii,8)) > er*MAX(one,ABS(tmp(ii,4)))) jc = 1
            END DO
!
            IF (jc == 0) GO TO 70
!
         END DO    ! nit
!
         IF (iprint == 1) WRITE(6,75)
75       FORMAT(1X,'no convergence of corrections')
!
         RETURN
!
70       CONTINUE
!
         DO ic = 1, ncomp
            def6(ic,im) = (hmsh/12.d+0)*(fval(ic,im)+  &
               5.d+0*(tmp(ic,3)+tmp(ic,4))+fval(ic,im+1))- u(ic,im+1)+u(ic,im)
!           write(6,*) 'Def6(',Ic,',',Im,') =', Def6(Ic,Im)
         END DO
!
         DO ic = 1, ncomp
            tmp(ic,5) = def6(ic,im)
            tmp(ic,6) = def6(ic,im)
         END DO
         CALL LUSOL(2*ncomp,2*ncomp,dhold,ip,tmp(1,5),tmp(1,7))
         def8(1,im)=tmp(1,7)
         def8(2,im)=tmp(2,7)
!
      END DO
!
      RETURN
   END SUBROUTINE DFEXCL
!!!
!!!
   SUBROUTINE EXPL(ncomp, nmsh, xx, nudim, u, defexp, fval, im, rpar, ipar)
!
      USE BVPShared, ONLY: pdebug, use_c, comp_c, uval0,  &
                           nminit, iprint, idum
!
      IMPLICIT NONE
!
      INTEGER(i4), INTENT(IN)    :: ncomp
      INTEGER(i4), INTENT(INOUT) :: nmsh
      REAL(r8),    INTENT(INOUT) :: xx(nmsh)
      INTEGER(i4), INTENT(IN)    :: nudim
      REAL(r8),    INTENT(INOUT) :: u(nudim,nmsh)
      REAL(r8),    INTENT(INOUT) :: defexp(ncomp)
      REAL(r8),    INTENT(INOUT) :: fval(ncomp, nmsh)
      INTEGER(i4), INTENT(IN)    :: im
      REAL(r8),    INTENT(INOUT) :: rpar(*)
      INTEGER(i4), INTENT(INOUT) :: ipar(*)
!
!     SCMODIFIED: increased the number of dimensions
!      dimension t1(2),t2(2),t3(2),t4(2)
!
      REAL(r8) :: t1(ncomp),t2(ncomp),t3(ncomp),t4(ncomp)
!
      REAL(r8), PARAMETER :: half = 0.5D+0
      REAL(r8), PARAMETER :: fourth = 0.25D+0
      REAL(r8), PARAMETER :: thfrth= 0.75D+0
      REAL(r8) :: a5, b5, c5, d5, e5, f5, a6, b6, c6
      REAL(r8) :: hmsh, au
      INTEGER(i4) :: ic
!!!
!
      a5=5.0D+0/32.0D+0
      b5=27.0D+0/32.0D+0
      c5=9.0D+0/64.0D+0
      d5=3.0D+0/64.0D+0
      e5=5.0D+0/24.0D+0
      f5=2.0D+0/3.0D+0
      a6=7.0D+0/90.0D+0
      b6=16.0D+0/45.0D+0
      c6=2.0D+0/15.0D+0
      hmsh = xx(im+1) - xx(im)
!
      DO ic = 1, ncomp
         t1(ic) = (a5*u(ic, im+1) + b5*u(ic, im))  &
            + hmsh*(c5*fval(ic,im)- d5*fval(ic,im+1))
         t2(ic) = (b5*u(ic,im+1) + a5*u(ic,im))  &
            + hmsh*(-c5*fval(ic,im+1) + d5*fval(ic,im))
      END DO
!
      CALL FSUB(ncomp, xx(im)+fourth*hmsh, t1, t3,rpar,ipar)
      CALL FSUB(ncomp, xx(im)+thfrth*hmsh, t2, t4,rpar,ipar)
!
      DO ic=1,ncomp
         t1(ic) = half*(u(ic,im+1) + u(ic,im))  &
            + e5*hmsh*(fval(ic,im+1) - fval(ic,im)) - f5*hmsh*(t4(ic) - t3(ic))
      END DO
!
      CALL FSUB(ncomp, half*(xx(im)+xx(im+1)), t1, t2,rpar,ipar)
!
      DO ic = 1, ncomp
         defexp(ic) = hmsh*(a6*(fval(ic,im+1) + fval(ic,im))  &
            + b6*(t3(ic) + t4(ic)) + c6*t2(ic)) - u(ic,im+1) + u(ic,im)
         au = MAX(ABS(u(ic,im)),ABS(u(ic,im+1)))
         au = 0.0D+0
         defexp(ic)=defexp(ic)/MAX(1.0D+0,au)
      END DO
!
78    CONTINUE
9000  CONTINUE
!
      RETURN
   END SUBROUTINE EXPL
!!!
!!!
   SUBROUTINE DF8CAL(ncomp, nmsh, xx, nudim, u, fval, def8, linear,  &
                     tmp, df, ip, dhold, ntol, ltol, tol, jc, rpar,ipar)
!
!   Given the mesh points xx, the solution u, and the function
!   values fval, df8cal computes eighth-order deferred corrections,
!   which are stored in def8.
!   The array tmp is workspace for 9 intermediate vectors.
!
      USE BVPShared, ONLY: pdebug, use_c, comp_c, uval0,  &
                           nminit, iprint, idum, &
                           a21b, a22b, a23b, a24b, a25b, a31b, a32b, a34b,  &
                           a35b, a41b, a42b, a43b, a44b, a45b,  &
                           b1b, b2b, b3b, c1b, c2b, c3b, c16b, c26b, c36b, &
                           c123b, c223b, c323b, c14b, c24b, c34b,  &
                           ifinal, iback, iprec,  &
                           LUFAC, LUSOL
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
      LOGICAL,     INTENT(IN)    :: linear
      REAL(r8),    INTENT(INOUT) :: tmp(ncomp,8)
      REAL(r8),    INTENT(INOUT) :: df(ncomp,ncomp)
      INTEGER(i4), INTENT(INOUT) :: ip(3*ncomp)
      REAL(r8),    INTENT(INOUT) :: dhold(3*ncomp,3*ncomp)
      INTEGER(i4), INTENT(IN)    :: ntol
      INTEGER(i4), INTENT(IN)    :: ltol(ntol)
      REAL(r8),    INTENT(IN)    :: tol(ntol)
      INTEGER(i4), INTENT(INOUT) :: jc
      REAL(r8),    INTENT(INOUT) :: rpar(*)
      INTEGER(i4), INTENT(INOUT) :: ipar(*)
!
      REAL(r8), PARAMETER :: one = 1.0D+0
      REAL(r8), PARAMETER :: two = 2.0D+0
!
      REAL(r8) :: st1(200), st2(200), st3(200), st4(200)
      REAL(r8) :: hmsh, c16h, c26h, c36h, fvim, fvim1, uim, uim1
      REAL(r8) :: xxc1, xxc2, xxc3, dfij, tmp4, tmp5, tmp6, er
      INTEGER(i4) :: im, ic, nit, i, j, ier, ii
!!!
!
      DO im = 1, nmsh-1
! 
         hmsh = xx(im+1) - xx(im)
! 
         c16h = c16b/hmsh
         c26h = c26b/hmsh
         c36h = c36b/hmsh
!
         DO ic = 1, ncomp
            fvim = fval(ic,im)
            fvim1 = fval(ic,im+1)
            uim = u(ic,im)
            uim1 = u(ic,im+1)
            tmp(ic,4) = c16h*(uim1-uim)+c123b*fvim1+c14b*fvim
            tmp(ic,5) = c26h*(uim1-uim)+c223b*fvim1+c24b*fvim
            tmp(ic,6) = c36h*(uim1-uim)+c323b*fvim1+c34b*fvim
            st1(ic) = (uim+uim1)/two
            st2(ic) = a21b*fvim + a25b*fvim1
            st3(ic) = a31b*fvim + a35b*fvim1
            st4(ic) = a41b*fvim + a45b*fvim1
         END DO
!
         xxc1 = xx(im)+c1b*hmsh
         xxc2 = xx(im)+c2b*hmsh
         xxc3 = xx(im)+c3b*hmsh
!
         DO nit = 1, 10
!
            DO ic = 1, ncomp
!
               tmp4 = tmp(ic,4)
               tmp5 = tmp(ic,5)
               tmp6 = tmp(ic,6)
               tmp(ic,1) = st1(ic) + hmsh*(st2(ic) + a22b*tmp4 + a23b*tmp5 + a24b*tmp6)
               tmp(ic,2) = st1(ic) + hmsh*(st3(ic) + a32b*tmp4 + a34b*tmp6)
               tmp(ic,3) = st1(ic) + hmsh*(st4(ic) + a42b*tmp4 + a43b*tmp5 + a44b*tmp6)
!
            END DO
!
            CALL FSUB(ncomp,xxc1,tmp(1,1),tmp(1,7),rpar,ipar)
            CALL FSUB(ncomp,xxc2,tmp(1,2),tmp(1,8),rpar,ipar)
            CALL FSUB(ncomp,xxc3,tmp(1,3),tmp(1,9),rpar,ipar)
!
            CALL DFSUB(ncomp,xxc1,tmp(1,1),df,rpar,ipar)
            DO i = 1, ncomp
               tmp(i,7) = tmp(i,7)-tmp(i,4)
               tmp(i,8) = tmp(i,8)-tmp(i,5)
               tmp(i,9) = tmp(i,9)-tmp(i,6)
               DO j = 1,ncomp
                  dfij = hmsh*df(i,j)
                  dhold(i,j) = -a22b*dfij
                  dhold(i,j+ncomp) = -a23b*dfij
                  dhold(i,j+2*ncomp) = -a24b*dfij
               END DO
            END DO
!
            CALL DFSUB(ncomp,xxc2,tmp(1,2),df,rpar,ipar)
            DO i = 1, ncomp
               DO j = 1, ncomp
                  dfij = hmsh*df(i,j)
                  dhold(i+ncomp,j) = -a32b*dfij
                  dhold(i+ncomp,j+ncomp) = 0.d+0
                  dhold(i+ncomp,j+2*ncomp) = -a34b*dfij
               END DO
            END DO
!
            CALL DFSUB(ncomp,xxc3,tmp(1,3),df,rpar,ipar)
            DO i = 1, ncomp
               DO j =  1, ncomp
                  dfij = hmsh*df(i,j)
                  dhold(i+2*ncomp,j) = -a42b*dfij
                  dhold(i+2*ncomp,j+ncomp) = -a43b*dfij
                  dhold(i+2*ncomp,j+2*ncomp) = -a44b*dfij
               END DO
            END DO
!
            DO i = 1, ncomp
               dhold(i,i) = dhold(i,i) + one
               dhold(i+ncomp,i+ncomp) = dhold(i+ncomp,i+ncomp) + one
               dhold(i+2*ncomp,i+2*ncomp) = dhold(i+2*ncomp,i+2*ncomp) + one
            END DO
!
            CALL LUFAC(3*ncomp,3*ncomp,dhold,ip,ier)
            CALL LUSOL(3*ncomp,3*ncomp,dhold,ip,tmp(1,7),tmp(1,10))
!
            DO i = 1,ncomp
               tmp(i,4) = tmp(i,4) + tmp(i,10)
               tmp(i,5) = tmp(i,5) + tmp(i,11)
               tmp(i,6) = tmp(i,6) + tmp(i,12)
            END DO
!
            jc = 0
            IF (linear) GO TO 90
            DO i = 1, ntol
               ii = ltol(i)
               er = tol(i)/hmsh
               IF (ABS(tmp(ii,10)) > er*MAX(one,ABS(tmp(ii,4))) .OR.  &
                   ABS(tmp(ii,11)) > er*MAX(one,ABS(tmp(ii,5))) .OR.  &
                   ABS(tmp(ii,12)) > er*MAX(one,ABS(tmp(ii,6)))) jc = 1
            END DO
!
            IF (jc == 0) GO TO 90
!
         END DO
!
         IF (iprint == 1) WRITE(6,8930)
8930     FORMAT(1X,'no convergence of 8th order defcors')
!
         RETURN
90       CONTINUE
!
         DO ic = 1, ncomp
            def8(ic,im) = hmsh*(b1b*(fval(ic,im)+fval(ic,im+1))+  &
                          b2b*(tmp(ic,4)+tmp(ic,6))+b3b*tmp(ic,5))- u(ic,im+1)+u(ic,im)
         END DO
!
      END DO
!
900   FORMAT(/,'** Warning - Possibly Approaching Machine Precision ',  &
      'Beyond Epsilon  = ',d10.3,/)
!
      RETURN
   END SUBROUTINE DF8CAL
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
   SUBROUTINE STCON1
!
!  stcons computes constants needed in integration formulae
!  and stores them in a labeled common area.
!
      USE BVPShared, ONLY: a21a, a22a, a23a, a24a, a31a, a32a, a33a, a34a,  &
                           c1a, c2a, c16a, c26a, c123a, c223a, c14a, c24a
!
      IMPLICIT NONE
!
      REAL(r8), PARAMETER :: one = 1.0D+0, two = 2.0D+0, three = 3.0D+0
      REAL(r8), PARAMETER :: four = 4.0D+0, five = 5.0D+0, six = 6.0D+0
!
      REAL(r8) :: rt5, c12a, c22a
!!!
!
      rt5 = SQRT(5.0D+0)
!
      a21a = (six + rt5)/120.0D0
      a22a = -rt5/120.0D0
      a23a = (-13.d0*rt5)/120.0D0
      a24a = (-six + rt5)/120.d0
!
      a31a = (six-rt5)/120.0D0
      a32a = (13.0D0*rt5)/120.0D0
      a33a = rt5 / 120.0D0
      a34a = (-six - rt5)/120.d0
!
      c1a = (five - rt5)/10.0D0
      c2a = (five + rt5)/10.0D0
!
      c12a = c1a*c1a
      c22a = c2a*c2a
!
      c16a = six*(c1a - c12a)
      c26a = six*(c2a - c22a)
!
      c123a = three*c12a - two*c1a
      c223a = three*c22a - two*c2a
!
      c14a = one - four*c1a + three*c12a
      c24a = one - four*c2a + three*c22a
!
      RETURN
   END SUBROUTINE STCON1
!!!
!!!
   SUBROUTINE STCON2
!
!  stcons computes constants needed in integration formulae
!  and stores them in a labeled common area.
!
      USE BVPShared, ONLY: a21b, a22b, a23b, a24b, a25b, a31b, a32b, a34b,  &
                           a35b, a41b, a42b, a43b, a44b, a45b,  &
                           b1b, b2b, b3b, c1b, c2b, c3b, c16b, c26b, c36b,  &
                           c123b, c223b, c323b, c14b, c24b, c34b
!
      IMPLICIT NONE
!
      REAL(r8), PARAMETER :: one = 1.0D+0, two = 2.0D+0, three = 3.0D+0
      REAL(r8), PARAMETER :: four = 4.0D+0, six = 6.0D+0
!
      REAL(r8) :: rt21, c12b, c22b, c32b
!!!
!
      rt21 = SQRT(21.0D+0)
!
      a21b = one/28.d0 + three*rt21/1960.d0
      a22b = -rt21/280.d0
      a23b = -32.d0*rt21/735.d0
      a24b = -23.d0*rt21/840.d0
      a25b = -one/28.d0 + three*rt21/1960.d0
!
      a31b = one/64.d0
      a32b = 7.d0*rt21/192.d0
      a34b = -7.d0*rt21/192.d0
      a35b = -one/64.d0
!
      a41b = one/28.d0 - three*rt21/1960.d0
      a42b = 23.d0*rt21/840.d0
      a43b = 32.d0*rt21/735.d0
      a44b = rt21/280.d0
      a45b = -(one/28.d0) - three*rt21/1960.d0
!
      b1b = one/20.0D0
      b2b = 49.0D0/180.d0
      b3b = 16.0D0/45.d0
!
      c1b = one/two - rt21/14.d0
      c2b = one/two
      c3b = one/two + rt21/14.d0
!
      c12b = c1b*c1b
      c22b = c2b*c2b
      c32b = c3b*c3b
!
      c16b = six*(c1b - c12b)
      c26b = six*(c2b- c22b)
      c36b = six*(c3b - c32b)
!
      c123b = three*c12b - two*c1b
      c223b = three*c22b - two*c2b
      c323b = three*c32b - two*c3b
!
      c14b = one - four*c1b + three*c12b
      c24b = one - four*c2b + three*c22b
      c34b = one - four*c3b + three*c32b
!
      RETURN
   END SUBROUTINE STCON2
!!!
!!!
   SUBROUTINE FIXJAC(ncomp, nmsh, nlbc, iorder, ntol, ltol, tol,  &
                     xx, nudim, u, defcor, defnew, delu, rhs, fval, utrial, rhstri,  &
                     rnsq, uint, ftmp, tmprhs, ajac, topblk, botblk, ipivot,  &
                     iflag, rpar, ipar)
!
      USE BVPShared, ONLY: pdebug, use_c, comp_c, uval0,  &
                           nminit, iprint, idum,  &
                           flmin, flmax, epsmch,  &
                           DCOPY, DSSQ, CRSLVE, MAXPY, MATCOP
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
      REAL(r8), PARAMETER :: huge = 1.0D+20
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
      CALL FNEVAL(ncomp, nmsh, xx, ncomp, utrial, fval, rpar,ipar)
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
                           DCOPY, COLROW1, DLOAD, CRSLVE, MAXPY
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
!
      IF (.NOT. ludone) THEN
! 
!        Compute the right-hand side vector rhs.
!
         iflag = 0
         CALL LNRHS(ncomp, nmsh, nlbc, xx, nudim, u,  &
                    rhs, rnsq, fval, ftmp, uint,rpar,ipar)
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
         CALL COLROW1(isize, topblk,nlbc, ncomp, ajac, ncomp, 2*ncomp,  &
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
         iflag = 0
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
                     iter, iflag, isign, rpar, ipar, frscal)
!
      USE BVPShared, ONLY: pdebug, use_c, comp_c, uval0,  &
                           nminit, iprint, idum,  &
                           flmin, flmax, epsmch,  &
                           DCOPY, COLROW1, MSSQ, MATCOP, MAXPY
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
      INTEGER(i4), INTENT(INOUT) :: isign(*)
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
      INTEGER(i4), PARAMETER :: itcondmax = 5
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
                     rhs, rnsq, fval, ftmp, uint,rpar,ipar)
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
      CALL COLROW1(nmsh*ncomp,topblk,nlbc,ncomp,ajac,ncomp,2*ncomp,  &
                   ninter,botblk,ncomp-nlbc,ipivot,delu,iflag,job)
!
      IF (iprint >= 0 .AND. iflag /= 0) WRITE(6,905) iter
      IF (iflag /= 0) THEN
         iflag = -1
         RETURN
      END IF
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
         CALL FNEVAL(ncomp, nmsh, xx, ncomp, utrial, fval, rpar, ipar)
         CALL RHSCAL(ncomp, nmsh, nlbc, xx, ncomp, utrial, defcor,  &
                     rhstri, rnsqtr, fval, ftmp, uint, rpar, ipar)
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
      REAL(r8),    PARAMETER :: grfct = 100.0D+0
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
   SUBROUTINE LNRHS(ncomp, nmsh, nlbc, xx, nudim, u,  &
                    rhs, rnsq, fval, ftmp, uint,rpar,ipar)
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
                           ifinal, iback, iprec,  &
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
!     nref = .true.
!
      IF (pdebug) WRITE(6,901) nmsh, ipow
!
      frcpow = one/ipow
      ddouble = .false.
      nmold = nmsh
      ninter = nmsh - 1
!     iprec = min(iprec,1)
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
         nmsh = 2*nmsh - 1
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
      IF (use_c) THEN
!        nodouble = ((iorder.eq.4) .and. (stiff_cond) .and. (use_c))
         nodouble = ((iorder == 4) .AND.  &
            (stiff_cond .AND. .NOT. stab_cond) .AND. (use_c))
!     nodouble = nodouble
!    *  .or.((iorder.gt.4) .and. ( stiff_cond .and. .not. stab_cond)
!    *   .and. (use_c))
         forcedouble = .false.
      END IF
      IF (use_c) THEN
         IF (itcond == itcondmax) THEN
            itcond = 0
            forcedouble = .true.
         ELSE
            itcond = itcond + 1
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
!     If(Iflnwt.eq.0) then
!       onto6=.true.
!     return
!     endif
      ddouble = .true.
      IF (pdebug) WRITE(6,903)
!     f the mesh if not doubled if the problem is stiff and the order is 4
!
      IF (.NOT. use_c .AND. .NOT. comp_c) nodouble = .false.
      IF (nodouble .AND. .NOT. forcedouble) THEN
         CALL SELCONDMSH(ncomp, nmsh, nfxpnt, fixpnt,  nmax, xx,  irefin,  &
                         nmold, xxold, ddouble , maxmsh,r4,amg)
         ddouble = .false.
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
      IF (.NOT. use_c .AND. .NOT. comp_c) nodouble = .false.
      IF (ddouble) THEN
!
!f       the mesh if not doubled if the problem is stiff
!
         IF  (nodouble .AND. .NOT. forcedouble)  THEN
            numadd = numpt
            CALL SMPSELCONDMSH(ncomp, nmsh, nfxpnt, fixpnt, nmax, xx, irefin, intref, numadd,  &
                               nmold, xxold, ddouble, maxmsh, r4, amg)
            ddouble = .false.
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
         ELSE
            numadd = numpt
            CALL SMPMSH(nmsh, nmax, xx, intref, numadd, nmold, xxold, maxmsh)
            itcond = 0
         END IF
!
      END IF
!
      IF (ddouble  .AND. pdebug) WRITE(6,905)
!     onto6 = .false.
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
!           IF (pdebug)  WRITE(6,901) im, it, errel, etest(it)
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
!     The routine selmsh performs selective mesh refinement, depending
!     on the error measure ermeas.
!
      IF (first) THEN
         first = .false.
         rlndec = DLOG(erdcid)
      END IF
!
      maxmsh = .false.
!     nref = .true.
!
      IF (pdebug) WRITE(6,901) nmsh, ipow
!
      frcpow = one/ipow
      ddouble = .false.
      nmold = nmsh
      ninter = nmsh - 1
!     iprec = min(iprec,1)
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
!
      IF (pdebug) WRITE(6,903) errmax
!
      CALL MONCONDMSH(nmsh,xx,r1,r2,r3,fatt_r1r3,fatt_r3,nptcond,r4,amg)
!
      IF (.NOT. stab_cond .AND.  errmax >= 1.0E20 .AND. r1 >= 1.0D0 ) THEN
!
!f       only the conditioning
         CALL SELCONDMSH(ncomp, nmsh, nfxpnt, fixpnt,  nmax, xx,  irefin,  &
                         nmold, xxold, ddouble, maxmsh,r4,amg)
!
      ELSE
!
         IF (errmax > zero .AND. errmax <= erdcid) THEN
! 
!           If errmax > 0 and .le. erdcid, find the smallest integer exponent ii
!           such that (erdcid**ii)*errmax > erdcid.
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
!           Multiply error measures by erdcid**ii.
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
200      CONTINUE
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
         IF (nptcond >= 4 .AND. .NOT. stab_cond) THEN
            add = .false.
            DO im = 1, ninter-1
               IF (MAX(r4(im), r4(im+1)) > fatt_r1r3) THEN
                  IF (.NOT. add) THEN
                     irefin(im) = MAX(nptcond, irefin(im))
!                    nmest = nmest + nptcond - 1
                  END IF
                  irefin(im+1) = MAX(nptcond, irefin(im+1))
!                 nmest = nmest + nptcond - 1
                  add = .true.
               ELSE
                  irefin(im) = MAX(1, irefin(im))
                  irefin(im+1) = MAX(1, irefin(im+1))
!                 nmest = nmest - 1
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
!        IF (pdebug) WRITE(6,904) nmest, (irefin(i), i=1,ninter)
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
               IF (jtkout == 0) THEN
                  ind1 = ihcomp(im-1)
                  phihat = ermeas(ind1,im-1)/(rlold**ipow)
               END IF
               phihat = MAX(phihat, ermeas(ihcomp(im),im)/(rlen**ipow))
               val1 = phihat*(slen**ipow)
               IF (val1 <= phitst .AND. jtkout < 4  &
                  .AND. r4(im) < fatt_r1r3) THEN
! 
!              Increment the counter of removed points.
!              'Remove' the mesh point xxold(im) by not including it.
! 
                  jtkout = jtkout+1
                  CYCLE
               END IF
!           end of logic for irefin(im) = 1.
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
360      CONTINUE
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
            nmsh = 2*nmsh - 1
            maxmsh = .true.
         END IF
!
!f only the conditioning
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
!     nptcond = nptcond/2
      DO im = 1, ninter-1
         IF (MAX(r4(im), r4(im+1)) >= fatt_r1r3) THEN
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
      add = .false.
      nmest = nmsh
!     nptcond = nptcond/2
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
            IF (jtkout < 1 .AND. (r4(im) <= 5d-1*fatt_r3 .AND. r1 >= 1.0d0)) THEN
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
      REAL(r8) :: r1m, r2m
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
         r4(i) = r4(i)+(xx(i+1)-xx(i))*(r2/(xx(nmsh)-xx(1)))*1.d-5
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
      r2m = r4(1)
      r1m = r4(1)
      DO i = 2, nmsh-1
         r1m = MIN(r1m,r4(i))
         r2m = r2m + r4(i)
      END DO
!
      r3 = r2m/(nmsh-1)
      fatt_r3  = r3*1.0d-4
      fatt_r1r3= MAX(r3,0.65D0)
!
      IF (r1 > 1.0d0) THEN
         r1m = one
         nptm = 0
         nptr = 0
         DO i = 1, nmsh-1
            IF (r4(i) >= fatt_r1r3)  nptm = nptm + 1
            IF (r4(i) <= fatt_r3)  nptr = nptr + 1
         END DO
         IF (nptm <= 1) THEN
            nptcond =  14
         ELSE IF (nptm <= 2) THEN
            nptcond =  10
         ELSE IF (nptm <= 4) THEN
            nptcond =  8
         ELSE IF (nptm <= 8) THEN
            nptcond = 6
         ELSE IF (nptm <= INT(nmsh/20) ) THEN
            nptcond = 4
         ELSE
            nptcond = 2
         END IF
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
