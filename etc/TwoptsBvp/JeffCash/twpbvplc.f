      subroutine twpbvplc(ncomp, nlbc, aleft, aright,
     *       nfxpnt, fixpnt, ntol, ltol, tol,
     *       linear, givmsh, giveu, nmsh,
     *       nxxdim, xx, nudim, u, nmax,
     *       lwrkfl, wrk, lwrkin, iwrk, 
     *       fsub, dfsub, gsub, dgsub, 
     *       ckappa1,gamma1,ckappa,rpar,ipar,iflbvp,
     *       liseries,iseries,indnms)

*  The subroutine twpbvplc is intended to solve two-point boundary
*  value problems. 
*
* References:  
*  Cash, J. R.; Mazzia, F. 
*      Hybrid Mesh Selection Algorithms Based on Conditioning for 
*      Two-Point Boundary Value Problems, Journal of Numerical Analysis,
*     Industrial and Applied Mathematics, in press
*
*   Revision History
*
* revision  July 10,  2006 
*   added rpar and ipar in the functions
*   DFSUB, DGSUB, FSUB, GSUB
*   changed the name of the variable double in ddouble 
*   
*
* revision 31 August 2004
* This is a modified version of twpbvp that uses the conditioning
* in the mesh selection and lobatto methods.
*
*     New subroutines not included in the old version:
*           condestim
*           moncond
*           selcondmsh
*           selconderrmsh
*           smpselcondmsh
* 
*     Updates subroutines:
*           bvpsol
*           fail4
*           decid4
*           fail6
*           newteq
*           mshref
*
*     Auxiliary function not used by the old version
*           donest
*           abdnrm
*    
*     The common block algprs contains two more variable
*     logical pdebug, use_c, comp_c
*     common/algprs/ nminit, pdebug, iprint, idum, uval0, use_c, comp_c
*


      implicit double precision (a-h,o-z)
      dimension rpar(*),ipar(*)
      dimension fixpnt(*), ltol(*), tol(*)
      dimension xx(*), u(nudim,*)
      dimension wrk(lwrkfl), iwrk(lwrkin)
      dimension iseries(*)
      logical linear, givmsh, giveu
      external fsub, dfsub, gsub, dgsub

      logical pdebug, use_c, comp_c
      common/algprs/ nminit, pdebug, iprint, idum, uval0, use_c, comp_c

      intrinsic abs, min

      parameter ( zero = 0.0d+0 )

*  Check for invalid input parameters.  If any parameters are
*  invalid, exit with the flag iflbvp set to -1.

      iflbvp = -1
      if (ncomp .le. 0)  return
      if (nlbc .lt. 0 .or. nlbc .gt. ncomp) return
      if (aleft .ge. aright) return

      if (nfxpnt .lt. 0)  return
      if (givmsh .and. nmsh .lt. nfxpnt+2) return
      if (givmsh .and. xx(1) .ne. aleft) return
C     SCMODIFIED add an extra condition to avoid accessing xx(0)
      if (nmsh .gt. 0) then
        if (givmsh .and. xx(nmsh) .ne. aright) return
      end if
      if (nfxpnt .gt. 0) then
         if (fixpnt(1) .le. aleft) return
         if (fixpnt(nfxpnt) .ge. aright) return
         do 50 i = 1, nfxpnt-1
            if (fixpnt(i+1) .le. fixpnt(i)) return
 50             continue
      endif

      if (ntol .lt. 1) return
      do 60 i = 1, ntol
         if (ltol(i) .lt. 0 .or. ltol(i) .gt. ncomp) return
         if (tol(i) .le. zero) return
 60       continue

      if (giveu .and. .not. givmsh) return

      if (use_c .and. .not. comp_c) return 

      if (nudim .le. 0) return
      if (lwrkfl .le. 0 .or. lwrkin .le. 0) return


*  Calculate maximum number of mesh points possible with the
*  given floating-point and integer workspace.

      isp = lwrkfl - 3 - 2*ntol - 22*ncomp - 6*ncomp*ncomp
      iden = 6*ncomp*ncomp + 22*ncomp + 3 
      nmax1 = isp/iden

      isp = lwrkin - 2*ncomp-3
      nmax2 = isp/(2*ncomp+3)

      nmax = min(nmax1,nmax2)
* nmax from workspace
      nmax = min(nmax, nxxdim)
* nmax from size of u and xx 

      if (iprint .ge. 0) write(6,901) nmax
 901   format(1h ,'nmax from workspace =',i8)

      if (nmax .le. 1) return


*  Partition floating point workspace.

      irhs = 1
      lrhs = ncomp*nmax

      itpblk = irhs + lrhs
      ltpblk = ncomp*nlbc

      ibtblk = itpblk + ltpblk
      lbtblk = ncomp*(ncomp - Nlbc)

      iajac = ibtblk + lbtblk
      lajac = 2*ncomp*ncomp*nmax

      ibhold = iajac + lajac
      lbhold = ncomp*ncomp*nmax

      ichold = ibhold + lbhold
      lchold = ncomp*ncomp*nmax

      ifval = ichold + lchold
      lfval = ncomp*nmax

      idef = ifval + lfval
      ldef = ncomp*(nmax-1)

      idefex = idef + ldef
*      ldefex = ncomp*(nmax-1)

*  def6 uses the same space as defexp

      idef6 = idefex
      ldef6 = ncomp*(nmax-1)

      idefim = idef6 + ldef6
*      ldefim = ncomp*(nmax-1)

*  def8 uses the same space as defimp

      idef8 = idefim
      ldef8 = ncomp*(nmax-1)

      iusve = idef8 + ldef8
      lusve = ncomp*nmax

      iuold = iusve + lusve
      luold = ncomp*nmax

      itmrhs = iuold + luold
      ltmrhs = ncomp*nmax

      irhtri = itmrhs + ltmrhs
      lrhtri = ncomp*nmax

      idelu = irhtri + lrhtri
      ldelu = ncomp*nmax

      ixmer = idelu + ldelu
*      lxmer = ncomp*nmax

*  rerr occupies the same space as xmerit 
      irerr = ixmer
      lrerr = ncomp*nmax

      iutri = irerr + lrerr
      lutri = ncomp*nmax

      iermx = iutri + lutri
      lermx = nmax

      irtdc = iermx + lermx
      lrtdc = nmax

      ixxold = irtdc + lrtdc
      lxxold = nmax
 
      iuint = ixxold + lxxold
      luint = ncomp

      iftmp = iuint + luint
      lftmp = ncomp

      idgtm = iftmp + lftmp
      ldgtm = ncomp

      idftm1 = idgtm + ldgtm
      ldftm1 = ncomp*ncomp

      idftm2 = idftm1 + ldftm1
      ldftm2 = ncomp*ncomp

      itmp = idftm2 + ldftm2
      ltmp = ncomp*8

      idsq = itmp + ltmp
      ldsq = ncomp*ncomp

      idexr = idsq + ldsq
      ldexr = ncomp

      ietst6 = idexr + ldexr
      letst6 = ntol

      ietst8 = ietst6 + letst6
      letst8 = ntol

      iamg = ietst8 + letst8
      lamg = ncomp*nmax

      ic1 = iamg + lamg
      lc1 = ncomp*ncomp*nmax

    
      iwrkrhs = ic1+lc1
      lwrkrhs = ncomp*nmax
       
      ir4 = iwrkrhs + lwrkrhs
      lr4 = nmax

      ibhold = ir4 + lr4
      lbhold = 9*ncomp*ncomp

      ilast = iwrkrhs +  lwrkrhs


      if (iprint .eq. 1) write(6,903) ilast
 903     format(1h ,'ilast',i10)


*  Partition integer workspace.
      
      iiref = 1
      liref = nmax

      iihcom = iiref + liref
      lihcom = nmax

      iipvbk = iihcom + lihcom
      lipvbk = ncomp*nmax

      iipvlu = iipvbk + lipvbk
      lipvlu = ncomp

      iisign = iipvlu + lipvlu
*      lisign = ncomp*nmax
   
      call bvpsol(ncomp, nmsh, nlbc, aleft, aright, 
     *   nfxpnt, fixpnt, ntol, ltol, tol, nmax, linear, 
     *   giveu, givmsh, xx, nudim, u, 
     *   wrk(idefex), wrk(idefim), wrk(idef), wrk(idelu), 
     *   wrk(irhs), wrk(ifval),
     *   wrk(itpblk), wrk(ibtblk), wrk(iajac), wrk(ibhold), 
     *   wrk(ichold), wrk(ibhold), iwrk(iipvbk), iwrk(iipvlu),
     *   iwrk(iisign), wrk(iuint), wrk(iftmp), wrk(itmrhs), 
     *   wrk(idftm1), wrk(idftm2), wrk(idgtm), 
     *   wrk(iutri), wrk(irhtri), wrk(ixmer), 
     *   wrk(ixxold), wrk(iuold), wrk(iusve),
     *   wrk(itmp), wrk(idsq), wrk(idexr), wrk(irtdc), 
     *   wrk(irerr), wrk(ietst6), wrk(ietst8), wrk(iermx), 
     *   iwrk(iihcom), iwrk(iiref), wrk(idef6), wrk(idef8),
     *   fsub,dfsub,gsub,dgsub,iflbvp,
     *   wrk(iamg),wrk(ic1),wrk(iwrkrhs),
     *   ckappa1,gamma1,ckappa,wrk(ir4),rpar,ipar,
     *    liseries,iseries,indnms)
      return
      end



      subroutine bvpsol(ncomp, nmsh, nlbc, aleft, aright, 
     *   nfxpnt, fixpnt, 
     *   ntol, ltol, tol, nmax, linear, giveu, givmsh, 
     *   xx, nudim, u, defexp, defimp, def, delu, rhs, fval,
     *   topblk, botblk, ajac, bhold, chold, dhold, 
     *   ipvblk, ipivlu,isign,
     *   uint, ftmp, tmprhs, dftmp1, dftmp2, dgtm, 
     *   utrial, rhstri, xmerit, xxold, uold, usave,
     *   tmp, dsq, dexr, ratdc, rerr, 
     *   etest6, etest8, ermx, ihcomp, irefin, 
     *   def6, def8, fsub, dfsub, gsub, dgsub, iflbvp,
     *   amg, c1, wrkrhs,ckappa1,gamma1,ckappa,r4,rpar,ipar,
     *    liseries,iseries,indnms)

      implicit double precision (a-h,o-z)

      dimension rpar(*), ipar(*)
      dimension  fixpnt(*), ltol(ntol), tol(ntol)
      dimension  xx(*), u(nudim, *)
      dimension  defexp(ncomp,*), defimp(ncomp,*), def(ncomp,*)
      dimension  delu(ncomp, *), rhs(*), fval(ncomp,*)
      dimension  topblk(nlbc,*), botblk(ncomp-nlbc,*)
      dimension  ajac(ncomp, 2*ncomp, *)
      dimension  bhold(ncomp, ncomp, *), chold(ncomp, ncomp, *)
      dimension  ipivlu(*), ipvblk(*), isign(*)
      dimension  uint(ncomp), ftmp(ncomp)
      dimension  dgtm(ncomp), tmprhs(*)
      dimension  dftmp1(ncomp, ncomp), dftmp2(ncomp, ncomp)
      dimension  utrial(ncomp,*), rhstri(*)
      dimension  xmerit(ncomp, *)
      dimension  xxold(*), uold(ncomp,*), usave(ncomp,*)
      dimension  tmp(ncomp,8)
      dimension  dsq(ncomp,ncomp), dexr(ncomp)
      dimension  ratdc(*), rerr(ncomp,*)
      dimension  etest6(*), etest8(*), ermx(*)
      dimension  ihcomp(*), irefin(*)
      dimension  def6(ncomp,*), def8(ncomp,*)
      dimension  amg(*), c1(ncomp,*), wrkrhs(*)
      Dimension  Dhold(3*ncomp,*),df(ncomp,ncomp)
      Dimension  r4(*), iseries(*)

      logical linear, giveu, givmsh, ddouble

      external fsub, dfsub, gsub, dgsub

      common/mchprs/flmin, flmax, epsmch

      logical pdebug, use_c, comp_c, Chstif
      common/algprs/ nminit, pdebug, iprint, idum, uval0, use_c, comp_c

      intrinsic max

      logical smooth, succes, strctr, trst6, reaft6
      logical onto6, onto8, ludone, rhsgiv, maxmsh
      logical first4, first8
      logical nodouble, forcedouble
      logical reposs

      logical mchset
      save mchset
      logical frscal
c      save frscal
    

      logical stab_kappa, stab_gamma, stab_cond, stiff_cond
     
      parameter (zero = 0.0d+0, one = 1.0d+0)
      parameter (third = 0.33d+0, fourth = 0.25d+0)
      parameter (quan6 = 0.5d+0 )
      parameter (itcondmax = 1)      
      parameter (power = 1.0d+0/6.0d+0)
*  blas: dload
*  double precision d1mach
      
      data mchset/.true./
      data fxfct/10.0d+0/
      data maxmsh/.false./
      data itcond/0/
      data frscal/.true./
c      frscal = .true. 
      if (mchset) then
         flmin = d1mach(1)
         flmax = d1mach(2)
         epsmch = d1mach(3)
         if (pdebug) write(6,901) epsmch
         mchset = .false.
      endif

*  The routine stcons calculates integration constants stored in
*  labeled common consts.

      Call Stcon1
      Call Stcon2

*  Set up arrays for the error tests.


      if (.not. linear) then
         call dload(ntol, one, etest6, 1)
      else
         do 10 i = 1, ntol
            etest6(i) = one/max(quan6, tol(i)**third)
 10             continue
      endif          

      do i=1,ncomp
         do j=1,ncomp
            df(i,j)=0.0d+0
         enddo
      enddo   
      nmold = 1
      smooth = .false.
      strctr = .false.
      trst6 = .true.
      reaft6 = .false.
      numbig = 0
      nummed = 0
      first4 = .true.
      first8 = .true.
      onto6  = .false.
      maxmsh = .false.
      Chstif = .true.
      ddouble = .false.

      if (comp_c) then
*     initialize parameter for the conditioning estimation  
      gamma1old  = flmax
      gamma1     = flmax
      ckappa1old = flmax
      ckappa1    = flmax
      ckappa     = flmax
      stiff_cond = .false.
      stab_cond  = .false.
      endif

      tolmin = flmax
      do i=1,ntol
         tolmin = min(tol(i),tolmin)
      end do
*  If givmsh is .true., the initial number of mesh points must be 
*  provided by the user in nmsh, and the mesh points must be
*  contained in the array xx (of dimension nmsh).
*  Otherwise, nmsh is set to its default value, and a
*  uniform initial mesh is created.
      
      if (.not. giveu .and. .not. givmsh) then
         nmsh = nminit
         if (nmsh .lt. nfxpnt+2) nmsh = nfxpnt + 2
         call unimsh(nmsh, aleft, aright, nfxpnt, fixpnt, xx)
      endif
      if (pdebug) then
         write(6,902)
         call sprt(nmsh, xx)
      endif
      
      if (.not. giveu) call initu(ncomp, nmsh, xx, nudim, u,rpar,ipar)
      indnms = 0
      indnmsold = 0
***** top of logic for 4th order solution ****
     
 400   continue
       if (indnmsold.ne.nmsh) then
       indnms = indnms + 1
       iseries(indnms) = nmsh
       indnmsold = nmsh
       endif
      
       if (indnms .ge. liseries) then
          write(6,1008) nmsh
          goto 1900
      end if
      If(Maxmsh) goto 900
      if (iprint .eq. 1) write(6,903) nmsh

*  Set the def (deferred correction) array to zero.
         
      call mtload(ncomp, nmsh-1, zero, ncomp, def)
      iorder = 4

*  The routine fneval calls fsub at the mesh xx and the
*  solution u, and saves the values in the array fval.
     
      call fneval(ncomp, nmsh, xx, nudim, u, fval, fsub,rpar,ipar)
      
     
*  Try to compute a 4th order solution by solving a system of nonlinear
*  equations.

      if (linear) then
         ludone = .false.
        
          call lineq( ncomp, nmsh, nlbc, 
     *    ludone, xx, nudim, u, def, 
     *    delu, rhs, fval, 
     *    uint, ftmp, dftmp1, dftmp2, dgtm, tmprhs,
     *    ajac, topblk, botblk, bhold, chold, ipvblk,
     *    fsub, dfsub, gsub, dgsub, iflnwt,rpar,ipar)

*  Call fneval to evaluate the fval array at the new solution u.
*  (Such a call is not necessary for the nonlinear case because
*  fval is called within newteq for the new u.)

         call fneval(ncomp, nmsh, xx, nudim, u, fval, fsub,rpar,ipar)
        
      else
      
         rhsgiv = .false.
         call newteq(ncomp, nmsh, nlbc, 
     *    rhsgiv, ntol, ltol, tol, 
     *    xx, nudim, u, def, 
     *    delu, rhs, fval,
     *    utrial, rhstri, 
     *    uint, ftmp, dftmp1, dftmp2, dgtm, tmprhs, xmerit, 
     *    ajac, topblk, botblk, bhold, chold, ipvblk,
     *    fsub, dfsub, gsub, dgsub, itnwt, iflnwt,isign,rpar,ipar,
     *    frscal)

      endif
*
*  these flags are used in the mesh selection strategy
* 
        


      if (iflnwt .eq. 0) then

c       COMPUTE ESTIMATIONS OF CONDITIONING NUMBERS: norms of inverse
c       jacobian matrix
c       BY BRUGNANO & TRIGIANTE, AND HIGHAM

        N =nmsh*ncomp
        ninter=nmsh-1        
        if (comp_c) then 
          gamma1old = gamma1
          ckappa1old = ckappa1        
            call CONDESTIM(aleft,aright,nmsh,ncomp,N,xx,topblk,nlbc,
     *       ncomp, ajac, ncomp,2*ncomp,ninter,botblk,ncomp-nlbc,
     *       ipvblk,isign,amg,c1,wrkrhs,ckappa1,gamma1)
           if (iprint .ge.0) then
            call ESTIMKAPPA(nmsh,ncomp,N,xx,topblk,
     *    	nlbc,ncomp,ajac, ncomp,2*ncomp,
     *          ninter,botblk,ncomp-nlbc,ipvblk,isign,c1,wrkrhs,ckappa)
           end if
          if (iprint .ge. 0) then
            write(6,1001) ckappa1/gamma1
            write(6,1002) gamma1 
            write(6,1003) ckappa1
            write(6,1004) ckappa
c            write(6,*) 'amg',( amg(i),i=1,nmsh)
          end if
          stab_kappa = abs(ckappa1old-ckappa1)/(1d0+ckappa1).lt.5d-2
     *      .and. ckappa1 .lt. flmax .and. gamma1.lt.flmax
 
          stab_gamma = abs(gamma1old-gamma1)/(1d0+gamma1).lt.5d-2
     *     .and. gamma1.lt.flmax .and. ckappa1 .lt. flmax

          stab_cond = stab_kappa .and. stab_gamma
         
     
          stiff_cond = (( (ckappa1/gamma1 .ge. 1d2 )))


          if (iprint .eq. 1) then
           write(6,*) 'stab_kappa = ',stab_kappa
           write(6,*) 'stab_gamma = ', stab_gamma
           write(6,*) 'stiff_cond = ', stiff_cond
          end if
        end if
c endif if (comp_c)
c
c  The subroutine Dfexcl substitute conv4
c
         Call Dfexcl(Ncomp, Nmsh, Xx, Nudim, U,def8,Def, Linear, Fval,
     *        Tmp, Fsub, Dfsub, Df, Ipivlu, Dhold,
     *        Ntol, Ltol, Tol, Jflag,rpar,ipar)

         if (reaft6) then
            onto6=.true.
            goto 408
         endif   
      else

       if (comp_c) then
         if (iflnwt .ne. -1) then
           gamma1old = gamma1
           ckappa1old = ckappa1
              N =nmsh*ncomp
           ninter=nmsh-1
          call CONDESTIM(aleft,aright,nmsh,ncomp,N,xx,topblk,nlbc,
     *       ncomp, ajac, ncomp,2*ncomp,ninter,botblk,ncomp-nlbc,
     *       ipvblk,isign,amg,c1,wrkrhs,ckappa1,gamma1)
           if (iprint .ge.0) then
            call ESTIMKAPPA(nmsh,ncomp,N,xx,topblk,
     *    	nlbc,ncomp,ajac, ncomp,2*ncomp,
     *          ninter,botblk,ncomp-nlbc,ipvblk,isign,c1,wrkrhs,ckappa)
           end if
         end if

         if (iprint .ge. 0) then
            write(6,1001) ckappa1/gamma1
            write(6,1002) gamma1 
            write(6,1003) ckappa1
            write(6,1004) ckappa
         end if
        
         stab_kappa = abs(ckappa1old-ckappa1)/(1d0+ckappa1).lt.5d-2
     *      .and. ckappa1 .lt. flmax .and. gamma1.lt.flmax
 
         stab_gamma = abs(gamma1old-gamma1)/(1d0+gamma1).lt.5d-2
     *     .and. gamma1.lt.flmax .and. ckappa1 .lt. flmax

         stab_cond = stab_kappa .and. stab_gamma
       
         stiff_cond = (( (ckappa1/gamma1 .ge. 1.0d1  )))
         if (iprint .eq. 1) then
           write(6,*) 'stab_kappa = ',stab_kappa
           write(6,*) 'stab_gamma = ', stab_gamma
           write(6,*) 'stiff_cond = ', stiff_cond
         end if
       end if
c end if if(comp_c)

         succes = .false.
         onto6 = .false.
         reaft6 = .false.

         call fail4( ncomp, nmsh, nlbc, ntol, ltol,
     *           xx, nudim, u, rhs, linear, nmax,
     *           nmold, xxold, uold, ratdc, 
     *           iorder, iflnwt, itnwt, ddouble, maxmsh,
     *           numbig, nummed,wrkrhs,amg,stab_cond,stiff_cond,
     *           nfxpnt, fixpnt, irefin,itcond,itcondmax,rpar,ipar)

c note: ratdc in subroutines fail4, fail6 does not related to 
c       ratdc=dfexmx/defimp , they only use the storage
         goto 400

      endif

      If(Jflag.eq.1) then

      nodouble = ( (stiff_cond .and. .not. stab_cond)
     * .and. (use_c))
c      nodouble = ((iorder.eq.4) .and. 
c     *       (stiff_cond .and. .not. stab_cond) .and. (use_c))

       forcedouble = .false.
 
       write(6,*) 'nodouble', itcond
       if (use_c) then
          if ( itcond .eq. itcondmax) then
             itcond = 0
             forcedouble = .true.
          else
             itcond = itcond + 1
             forcedouble = .false.
          endif
       endif
       if (nodouble .and. .not. forcedouble) then
          call selcondmsh(ncomp, nmsh,
     *     nfxpnt, fixpnt,  nmax, xx,  irefin,
     *     nmold, xxold, ddouble, maxmsh,r4,amg)
           ddouble = .false.
        else
         Call Dblmsh(Nmsh, Nmax, Xx, Nmold, Xxold, Maxmsh)
         itcond = 0
        endif 
         Call Matcop(Nudim, Ncomp, Ncomp, Nmold, U, Uold)
         Call Interp(Ncomp, Nmsh, Xx, Nudim, U, Nmold, Xxold, Uold)
         If(Maxmsh) goto 900
         goto 400
      Else
c    find where biggest deferred correction is
c        If(Chstif) Then
           smaldef=1.0D+40
           Bigdef=0.0D+0
           Icmph = 1
           Ix = 1
           Do 405 Iv=1,Nmsh-1
           Do 405 Iu = 1,Ntol
             Ipoint = Ltol(Iu)
             Holdef=Abs(Def8(Ipoint,Iv))
             if(smaldef.gt.holdef) smalldef=holdef
             If(Holdef.gt.Bigdef) Then
              Bigdef = Holdef
              Icmph=Ipoint
              Ixx = Iv
              intol = Iu
             Endif
 405       Continue
c   Biggest deferred correction is in component Icmph and
c   at the mesh interval Ix.
c   Now compute an explicit deferred correction for this.
         Call expl(Ncomp,Nmsh,xx,Nudim,U,dgtm,fval,fsub,Ixx,rpar,ipar)

           Ix=Ixx
           siz=abs(dgtm(Icmph))
c           write(6,*) 'siz=',siz
           rat=siz/bigdef
c           write(6,*) 'rat=',rat
c           write(6,*) ' bigdef=',bigdef,' > ',dsqrt(tol(Icmph))
c           write(6,*) ' dsq=', dsqrt(tol(Icmph)) 
          If(Rat.gt.50.0D+0.and.bigdef.gt.dsqrt(tol(Icmph)).and.
     +      siz.gt.1.0D+0) Then

        if (use_c) then
         if ((stiff_cond)) then
             drat = bigdef/
     *               (max(one, abs(u(icmph,Ix)))*tol(intol))

            if (pdebug) write(6,913) drat, u(icmph,Ix), tol(intol)
c            numadd = drat**power
             numadd = 15
            call smpselcondmsh(ncomp, nmsh,
     *        nfxpnt, fixpnt,  nmax, xx,  irefin,Ix,numadd,
     *        nmold, xxold, ddouble, maxmsh,r4,amg)
         else
             drat = bigdef/
     *               (max(one, abs(u(icmph,Ix)))*tol(intol))
            if (pdebug) write(6,913) drat, u(icmph,Ix), tol(intol)
c            numadd = drat**power
             numadd = 15
                call smpmsh (nmsh, nmax, xx, ix, numadd,
     *             nmold, xxold, maxmsh)
         endif
         else
             drat = bigdef/
     *               (max(one, abs(u(icmph,Ix)))*tol(intol))
            if (pdebug) write(6,913) drat, u(icmph,Ix), tol(intol)
            numadd = 15
c            numadd = drat**power
            call smpmsh (nmsh, nmax, xx, Ix, numadd,
     *             nmold, xxold, maxmsh)
         end if
 913     format(1h ,'drat,u,tol',3(1pe11.3))

             goto 400
           else
c             Chstif = .false.
             onto6=.true.
             if (linear.and.ddouble) reposs=.true.
           endif
c        Endif
c endif of Chstif=true

c       Iorder = 6
c       if (iflnwt.ne.0) then
c         call mshref(ncomp, nmsh, nlbc, ntol, ltol,
c     *             iorder, rhs, ratdc,
c     *             nmax, xx, nmold, xxold, ddouble, maxmsh,
c     *             numbig, nummed,
c     *             amg,stab_cond,stiff_cond,
c     *             r4, nfxpnt,fixpnt, irefin,itcond,itcondmax)
c       elseif((.not.onto6)) then
c         Call Matcop(Nudim, Ncomp, Ncomp, Nmold, U, Uold)
c         Call Interp(Ncomp, Nmsh, Xx, Nudim, U, Nmold, Xxold, Uold)
c         goto 400
c       endif
       Endif
 408     continue
      Call Matcop(Ncomp, Ncomp, Ncomp, Nmsh-1, Def, Def6)

c      if (succes) then
c          if (iprint .ne. -1 .and. comp_c ) then
c            if (ckappa .ge. tolmin/epsmch) write(6,1005)  
c          end if
c          iflbvp = 0
c          return
      if (maxmsh) then
          go to 900
      elseif (.not. onto6)  then
          go to 400
       endif

*  To reach here, onto6 must be .true.

**** logic for 6th order ****

      if (iprint .eq. 1) write(6,905)

*  Save the 4th order solution on this mesh in uold. 

      call matcop(nudim, ncomp, ncomp, nmsh, u, uold)
*  Copy the current mesh into the xxold array.
      nmold = nmsh
      call dcopy(nmold, xx, 1, xxold, 1)

      iorder = 6

      if (linear) then
         call lineq( ncomp, nmsh, nlbc, 
     *    ludone, xx, nudim, u, def, 
     *    delu, rhs, fval, 
     *    uint, ftmp, dftmp1, dftmp2, dgtm, tmprhs,
     *    ajac, topblk, botblk, bhold, chold, ipvblk,
     *    fsub, dfsub, gsub, dgsub, iflnwt,rpar,ipar)

      else

         call fixjac(ncomp, nmsh, nlbc, 
     *     iorder, ntol, ltol, tol, 
     *     xx, nudim, u, def, def, delu, 
     *     rhs, fval, utrial, rhstri, 
     *     rnsq, uint, ftmp, tmprhs,
     *     ajac, topblk, botblk, ipvblk,
     *     fsub, gsub, iflnwt,rpar,ipar)
          
*  If the fixed Jacobian iterations fail but rnsq is small,
*  try a Newton procedure.  Set rhsgiv to indicate that
*  the right-hand side and fval have already been evaluated
*  at the given u.
        
         if (iflnwt .eq. -3 .and. rnsq .lt. fxfct*epsmch) then
            rhsgiv = .true.
             
            call newteq(ncomp, nmsh, nlbc,
     *           rhsgiv, ntol, ltol, tol, 
     *           xx, nudim, u, def, 
     *           delu, rhs, fval,
     *           utrial, rhstri, 
     *           uint, ftmp, dftmp1, dftmp2, dgtm, tmprhs, xmerit, 
     *           ajac, topblk, botblk, bhold, chold, ipvblk,
     *           fsub, dfsub, gsub, dgsub, iter, iflnwt,isign,rpar,ipar,
     *              frscal)
   
         endif
      endif
   
      if (iflnwt .eq. 0) then
c         write(6,*) 'before etest=',etest6(1)
         itcond=0
         call conv6(ncomp, nmsh, ntol, ltol, tol,
     *             nudim, u, uold, etest6, err6,
     *             trst6, onto8, reaft6, linear, succes)
c         write(6,*) 'etest6=',etest6(1)

      else 

         onto8 = .false.
         
         call fail6( ncomp, nmsh, nlbc, ntol, ltol, tol,
     *              nfxpnt, fixpnt, iorder, nmax,
     *              xx, nudim, u, rhs, usave, xxold, uold, nmold, 
     *              ihcomp, irefin,
     *              rerr, ermx, ratdc,
     *              reaft6, ddouble, succes, maxmsh,
     *              numbig, nummed,
     *           wrkrhs,amg, stab_cond,ckappa1,gamma1,ckappa,stiff_cond,
     *              itcond, itcondmax)

   

      endif
      If(Maxmsh) goto 900
      if (succes) then
       if (comp_c) then
           call ESTIMKAPPA(nmsh,ncomp,N,xx,topblk,
     *    	nlbc,ncomp,ajac, ncomp,2*ncomp,
     *          ninter,botblk,ncomp-nlbc,ipvblk,isign,c1,wrkrhs,ckappa)
         end if
         if (iprint .ne. -1 .and. comp_c ) then
           if ( ckappa .ge. tolmin/epsmch) write(6,1005)  
         end if
         iflbvp = 0
         return
       elseif (.not. onto8) then
         go to 400
       endif

***** logic for trying to calculate 8th order solution *****

      if (iprint .eq. 1) write(6,906)

      call matcop(nudim, ncomp, ncomp, nmsh, u, uold)
*  Copy the current mesh into the xxold array.
      nmold = nmsh
      call dcopy(nmold, xx, 1, xxold, 1)


*  Save the old deferred correction vector def in def6.

      call matcop(ncomp, ncomp, ncomp, nmsh-1, def, def6)

*  For linear problems, calculate the fval array for the
*  new solution u.

      if (linear) call fneval(ncomp, nmsh, xx, nudim, u, fval, fsub,
     *       rpar,ipar)

*  Calculate 8th order deferred corrections (the array def8).

      Call Df8cal (Ncomp, Nmsh, Xx, Nudim, U, Fval, Def8, Linear,
     *             Tmp, Fsub, Dfsub, Df, Ipivlu, Dhold,
     *             Ntol, Ltol, Tol,JC,rpar,ipar)

      If(Jc.eq.1) then
      nodouble = ((stiff_cond) .and. .not.stab_cond  .and. (use_c))
c      nodouble = nodouble
c     *  .or.((iorder.gt.4) .and. .not.stab_cond .and. (stiff_cond)
c     *   .and. (use_c))
      forcedouble = .false.

      if (use_c) then
         if ( itcond .eq. itcondmax) then
            itcond = 0
            forcedouble = .true.
         else
            itcond = itcond + 1
            forcedouble = .false.
         endif
      endif
       if (nodouble .and. .not. forcedouble) then
          call selcondmsh(ncomp, nmsh,
     *     nfxpnt, fixpnt,  nmax, xx,  irefin,
     *     nmold, xxold, ddouble, maxmsh,r4,amg)
           ddouble = .false.
        else
         Call Dblmsh(Nmsh, Nmax, Xx, Nmold, Xxold, Maxmsh)
          itcond = 0
        endif 
         If(Maxmsh) goto 900
         Call Matcop(Nudim, Ncomp, Ncomp, Nmold, Uold, U)
         Call Interp(Ncomp, Nmsh, Xx, Nudim, U, Nmold, Xxold, Uold)
      goto 400
      endif


*  For linear problems, the def array is the def8 array.
*  For nonlinear problems, add the def8 array to the 
*  already-calculated def array.

      if (linear) then
         call matcop(ncomp, ncomp, ncomp, nmsh-1, def8, def)
      else
         call maxpy(ncomp, nmsh-1, one, def8, ncomp, def)
      endif

      iorder = 8

      if (linear) then
         call lineq( ncomp, nmsh, nlbc,  
     *    ludone, xx, nudim, u, def, 
     *    delu, rhs, fval, 
     *    uint, ftmp, dftmp1, dftmp2, dgtm, tmprhs,
     *    ajac, topblk, botblk, bhold, chold, ipvblk,
     *    fsub, dfsub, gsub, dgsub, iflnwt,rpar,ipar)

      else

         call fixjac(ncomp, nmsh, nlbc, 
     *     iorder, ntol, ltol, tol, 
     *     xx, nudim, u, def, def8, delu, 
     *     rhs, fval, utrial, rhstri, 
     *     rnsq, uint, ftmp, tmprhs,
     *     ajac, topblk, botblk, ipvblk,
     *     fsub, gsub, iflnwt,rpar,ipar)
        
*  If the fixed Jacobian iterations fail but rnsq is small,
*  try a Newton procedure.  Set rhsgiv to indicate that
*  the right-hand side and fval have already been evaluated
*  at the given u.

         if (iflnwt .eq. -3 .and. rnsq .lt. fxfct*epsmch) then
            rhsgiv = .true.
             
            call newteq(ncomp, nmsh, nlbc, 
     *           rhsgiv, ntol, ltol, tol, 
     *           xx, nudim, u, def, 
     *           delu, rhs, fval,
     *           utrial, rhstri, 
     *           uint, ftmp, dftmp1, dftmp2, dgtm, tmprhs, xmerit, 
     *           ajac, topblk, botblk, bhold, chold, ipvblk,
     *           fsub, dfsub, gsub, dgsub, iter, iflnwt,isign,rpar,ipar,
     *            frscal)
             
         endif
      endif
      
      if (iflnwt .eq. 0) then
         itcond = 0
         call conv8( ncomp, nmsh, ntol, ltol, tol,
     *              nfxpnt, fixpnt, linear, nmax,
     *              xx, nudim, u, def, def6, def8, uold, 
     *              ihcomp, irefin, ermx, err6,
     *              etest8, strctr,
     *              ddouble, nmold, xxold, maxmsh, succes, first8,
     *              wrkrhs,amg, stab_cond,ckappa1,gamma1,ckappa,
     *                                   stiff_cond,rpar,ipar)

      else 

         succes = .false.
         call  fail8(ncomp, nmsh, nfxpnt, fixpnt, nmax,
     *      ntol, ltol, tol, nmold,
     *      xx, nudim, u, def6, xxold, uold,
     *      ihcomp, irefin, ermx, ddouble, maxmsh,
     *       wrkrhs,amg, stiff_cond, stab_cond)

      endif

      if (maxmsh) then
         go to 900
      elseif (.not. succes) then
         go to 400
      endif

*  Successful termination.

       if (comp_c) then
           call ESTIMKAPPA(nmsh,ncomp,N,xx,topblk,
     *    	nlbc,ncomp,ajac, ncomp,2*ncomp,
     *          ninter,botblk,ncomp-nlbc,ipvblk,isign,c1,wrkrhs,ckappa)
         end if
      if (iprint .ne. -1 .and. comp_c ) then
        if ( ckappa .ge. tolmin/epsmch) write(6,1005)  
      end if
      iflbvp = 0

      return

 900   continue

* Error exit---too many mesh points.

      iflbvp = 1
      write(6,*) 'Terminated, too many mesh points'
      if (comp_c) then
          call ESTIMKAPPA(nmsh,ncomp,N,xx,topblk,
     *    	nlbc,ncomp,ajac, ncomp,2*ncomp,
     *          ninter,botblk,ncomp-nlbc,ipvblk,isign,c1,wrkrhs,ckappa)
        if (linear .and. ckappa .ge. tolmin/epsmch) write(6,1006) 
        if (.not.linear .and. ckappa .ge. tolmin/epsmch) write(6,1007) 
      end if
      return
 1900   continue
 
* Error exit---too many meshes  .return

      iflbvp = 1
      write(6,*) 'Terminated, too many meshes'
      if (comp_c) then
          call ESTIMKAPPA(nmsh,ncomp,N,xx,topblk,
     *    	nlbc,ncomp,ajac, ncomp,2*ncomp,
     *          ninter,botblk,ncomp-nlbc,ipvblk,isign,c1,wrkrhs,ckappa)
        if (linear .and. ckappa .ge. tolmin/epsmch) write(6,1006) 
        if (.not.linear .and. ckappa .ge. tolmin/epsmch) write(6,1007) 
      end if
 901  format(1h ,'epsmch',1pe10.3)  
 902  format(1h ,'initial mesh')
 903  format(1h ,'start 4th order, nmsh',i5)
 904  format(1h ,'do not go on to 6th')
 905  format(1h ,'start 6th order')
 906  format(1h ,'start 8th order')
 1001 format(1h ,'stiffness = ',1pe11.3)
 1002 format(1h ,'gamma1    = ',1pe11.3)
 1003 format(1h ,'kappa1    = ',1pe11.3)
 1004 format(1h ,'kappa     = ',1pe11.3)
 1005 format(1h ,'The problem is ill-conditioned',
     *     'the solution could be inaccurate')
 1006 format(1h ,'The problem is ill-conditioned',
     *     'try with a less stringent tolerance')
 1007 format(1h ,'The problem is ill-conditioned',
     *     'try with a less stringent tolerance',
     *     ' or with a different initial guess' )
 1008 format(1h ,'Terminated too many meshes, nmsh',i5)
      end



      block data

*  This block data routine initializes nminit (the initial number 
*  of mesh points), pdebug (a logical variable indicating whether
*  debug printout is desired), uval0 (the initial value for the trial
*  solution) to their default values, use_c (if the conditioning is used 
*  in the mesh selection), and comp_c (if the conditining parameters
*  are computed)

      double precision uval0
      logical pdebug, use_c, comp_c
      common/algprs/nminit, pdebug, iprint, idum, uval0, use_c, comp_c

      data nminit/10/
      data pdebug/.true./
      data iprint/1/
      data uval0/0.0d+0/
      data use_c/.false./
      data comp_c/.false./
      end

        SUBROUTINE CONDESTIM(ALEFT,ARIGHT,NMSH,NCOMP,N,XX,TOPBLK,
     *    	NRWTOP,NOVRLP,ARRAY,
     *          NRWBLK,NCLBLK,NBLOKS,BOTBLK,NRWBOT,IPVCD,ISIGN,OMG,
     *          C1,WORK,KPPA,GAMMA)

C     **************************************************************
C
C     COMPUTES THE FIRST AND LAST BLOCK COLUMN OF INVERSE MATRIX AND
C     ESTIMATE THE CONDITION NUMBERS KPPA, GAMMA 
C
c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        IMPLICIT NONE
        DOUBLE PRECISION TOPBLK,ARRAY,BOTBLK,WORK,C1,
     *          OMG,GAMMA,gamma1,KPPA,MINMG
        DOUBLE PRECISION ALEFT,ARIGHT,XX,  CSUM
        INTEGER ISIGN(*)
        INTEGER N,NRWTOP,NOVRLP,NRWBLK,NCLBLK,NBLOKS,NRWBOT,IPVCD
        INTEGER NCOMP,NMSH,idmx,idmn,idamax,idomg,job
        DIMENSION TOPBLK(NRWTOP,*),ARRAY(NRWBLK,NCLBLK,*),OMG(*),
     *          BOTBLK(NRWBOT,*),WORK(*),C1(ncomp,*),XX(*),
     *          IPVCD(*)


        INTEGER k,i,j,l
       
        MINMG=0.0d+0
        KPPA=0.0d+0
        do 300 k=1,ncomp
           do 330 l=1,N
              WORK(l)=0.0d0
              if (k.le.NRWTOP) then
                   WORK(k)=1.0d0
              else
                   WORK(N-ncomp+k)=1.0d0
              endif
 330    continue
  
        CALL CRSLVE(TOPBLK,NRWTOP,NOVRLP,ARRAY,NRWBLK,
     *                 NCLBLK,NBLOKS,BOTBLK,NRWBOT,IPVCD,
     *                 WORK,0)

       
           do 340 l=1,N
             C1(k,l)=WORK(l)
 340       continue
      
 300            continue

C
C     ESTIMATION OF CONDITION NUMBER BY BRUGNANO,TRIGIANTE & MAZZIA
C
c     infinity-norm
c 
        l=1
        do 400 j=1,N-ncomp+1,ncomp
           OMG(l) = 0d0
           do  i = j,j+ncomp-1
               CSUM = 0d0
               do k=1,ncomp
                  CSUM = CSUM + DABS(C1(k,i))
               end do
               OMG(l) = DMAX1( OMG(l),CSUM)
           end do
           

           l=l+1
 400                  continue
c
c     1-norm
c
c          l=1
c        do 400 j=1,N-ncomp+1,ncomp
c           OMG(l) = 0d0
c           do  k = 1,ncomp
c               CSUM = 0d0
c               do i=j,j+ncomp-1
c                  CSUM = CSUM + DABS(C1(k,i))
c               end do
c               OMG(l) = DMAX1( OMG(l),CSUM)
c           end do
c           
c
c           l=l+1
c 400                  continue     
  
      
        GAMMA=0.0D0
        do 450 i=2,NBLOKS+1
           if (OMG(I).GT.OMG(I-1)) then
               gamma1=OMG(I)* (XX(I)-XX(I-1))
           else
               gamma1=OMG(I-1)* (XX(I)-XX(I-1))
           end if
           GAMMA=GAMMA + gamma1

 450                  continue
        IDOMG=IDAMAX(NBLOKS+1,OMG,1)

        GAMMA=GAMMA/(ARIGHT-ALEFT)

        idmx=1
        idmn=1
        MINMG = OMG(1)
        KPPA = OMG(1)

        do  500 j=2,NBLOKS+1
           if (OMG(J).GT.KPPA) then
               KPPA=OMG(J)
               idmx=j
           endif
           if (OMG(J).LT.MINMG) then
               MINMG=OMG(J)
               idmn=j
           endif
 500    continue
         


         RETURN
          END

        SUBROUTINE ESTIMKAPPA(NMSH,NCOMP,N,XX,TOPBLK,
     *    	NRWTOP,NOVRLP,ARRAY, NRWBLK,
     *          NCLBLK,NBLOKS,BOTBLK,NRWBOT,IPVCD,ISIGN,C1,WORK,ckappa)

C     **************************************************************
C
C  
C     ESTIMATE THE CONDITION NUMBER  ckappa
C
c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        IMPLICIT NONE
        DOUBLE PRECISION TOPBLK,ARRAY,BOTBLK,WORK,C1
        DOUBLE PRECISION XX
        INTEGER ISIGN(*)
        INTEGER N,NRWTOP,NOVRLP,NRWBLK,NCLBLK,NBLOKS,NRWBOT,IPVCD
        INTEGER NCOMP,NMSH,job
        DIMENSION TOPBLK(NRWTOP,*),ARRAY(NRWBLK,NCLBLK,*),
     *          BOTBLK(NRWBOT,*),WORK(*),XX(*),C1(NCOMP,*),
     *          IPVCD(*)


        double precision ckappa
        INTEGER KASE,ISOLVE
        INTEGER i,j

    
c
c       DO THE CONDITION NUMBER ESTIMATION BY HIGHAM:
C       (DONE IN COLROW) infinity-norm
c        ckappa = 0.0D0

        ISOLVE = 0
        KASE = 0
 55               CALL DONEST(N,C1,WORK,ISIGN,ckappa,KASE)

        IF (KASE .NE. 0) THEN
           ISOLVE = ISOLVE+1

           IF (KASE .eq. 1) THEN
               JOB = 1
           ELSE
               JOB = 0
           END IF
         
           IF (JOB .eq. 0) THEN
             DO i=1,NBLOKS
               DO j=1,NRWBLK
                  WORK( (i-1)*NRWBLK+NRWTOP+j ) = 
     *             (XX(i+1)-XX(i))*WORK( (i-1)*NRWBLK+NRWTOP+j) 
               ENDDO
             END DO 
           END IF
           CALL CRSLVE(TOPBLK,NRWTOP,NOVRLP,ARRAY,NRWBLK,
     *                 NCLBLK,NBLOKS,BOTBLK,NRWBOT,IPVCD,
     *                 WORK,JOB)
           IF (JOB .eq. 1) THEN
             DO i=1,NBLOKS
               DO j=1,NRWBLK
                  WORK( (i-1)*NRWBLK+NRWTOP+j ) = 
     *             (XX(i+1)-XX(i))*WORK( (i-1)*NRWBLK+NRWTOP+j) 
               ENDDO
             END DO 
           END IF
           GOTO 55
        END IF
         

        RETURN
        END






      subroutine fail4( ncomp, nmsh, nlbc, ntol, ltol,
     *             xx, nudim, u, rhs, linear, nmax,
     *             nmold, xxold, uold, tmwork,
     *             iorder, iflnwt, itnwt, ddouble, maxmsh,
     *             numbig, nummed,r4,amg,stab_cond,stiff_cond,
     *              nfxpnt, fixpnt, irefin,itcond,itcondmax,rpar,ipar)


      implicit double precision (a-h,o-z)
      dimension rpar(*),ipar(*)
      dimension ltol(ntol)
      dimension xx(*), u(nudim, *), rhs(*)
      dimension xxold(*), uold(ncomp, *), tmwork(*)
      logical linear, ddouble, maxmsh
      logical stab_cond,stiff_cond
      dimension amg(*),r4(*), fixpnt(*), irefin(*)
      logical pdebug, use_c, comp_c
      common/algprs/ nminit, pdebug, iprint, idum, uval0, use_c, comp_c

*  The Newton procedure failed to obtain a 4th order solution.

      if (iprint .eq. 1) write(6,901)

      maxmsh = .false.
c      write(6,*) 'iflnwt',iflnwt, itnwt

      if (iflnwt .eq. -1) then

*  iflnwt = -1 means that the Jacobian was considered singular.
*  (This is the only possible failure for a linear problem.) 
        
         call dblmsh (nmsh, nmax, xx, nmold, xxold, maxmsh) 
         call initu(ncomp, nmsh, xx, nudim, u,rpar,ipar)
      else
*  The routine mshref decides how to refine the mesh and then
*  performs the refinement, either by doubling or based on 
*  the rhs vector at the best point.

         
         call mshref(ncomp, nmsh, nlbc, ntol, ltol,
     *             iorder, rhs, tmwork,
     *             nmax, xx, nmold, xxold, ddouble, maxmsh,
     *             numbig, nummed,
     *             amg,stab_cond,stiff_cond,
     *             r4, nfxpnt,fixpnt, irefin,itcond,itcondmax)
    

         if (.not. maxmsh) then
              if (linear  .or. itnwt .eq.0)  then 
c     *		   .or. (iflnwt .le.-4 .and. comp_c) ) then
                call initu(ncomp, nmsh, xx, nudim, u,rpar,ipar)
              else
*  Interpolate the partially converged solution.
               call matcop(nudim, ncomp, ncomp, nmold, u, uold)
               call interp(ncomp, nmsh, xx, nudim, u, nmold,
     *                 xxold, uold) 
              endif
         endif

*     End of logic for failure because of some reason other than
*     singularity.
      endif

      return
 901   format(1h ,'fail4')
      end





      subroutine conv6(ncomp, nmsh, ntol, ltol, tol,
     *             nudim, u, uold, etest6, err6,
     *             trst6, onto8, reaft6, linear,succes)

      implicit double precision (a-h,o-z)
      dimension ltol(ntol)
      dimension u(nudim,*), tol(ntol) 
      dimension uold(ncomp,*), etest6(*)
      logical trst6, onto8, reaft6, succes 
      parameter (quan6=0.5d+0,Third = 0.33d+0,one=1.d+0)
      logical linear
      logical pdebug, use_c, comp_c
      common/algprs/ nminit, pdebug, iprint, idum, uval0, use_c, comp_c

      logical errok

*  The Newton iteration converged for a 6th-order solution.

      if (.not. linear) then
         call dload(ntol, one, etest6, 1)
      else
         do 10 i = 1, ntol
            etest6(i) = one/max(quan6, tol(i)**third)
 10                      continue
      endif

      call dload(ntol, 10d0, etest6, 1)
      if (iprint .eq. 1) write(6,901)

      succes = .false.

*  The logical flag reaft6 is true only when the 6th order solution
*  failed.  Since the 6th order solution converged, reaft6 is false.

c      reaft6 = .false.
      onto8 = .true.

* Calculate the error estimates for the 4th order solution.

      call errest (ncomp, nmsh, ntol, ltol, tol,
     *          nudim, u, uold, etest6, err6, errok)

      if (trst6 .and. errok) then
c      if (errok) then
         succes = .true.
         return
      endif

      return
 901   format(1h ,'conv6')
      end



      subroutine fail6( ncomp, nmsh, nlbc, ntol, ltol, tol,
     *              nfxpnt, fixpnt, iorder, nmax,
     *              xx, nudim, u, rhs, usave, xxold, uold, nmold, 
     *              ihcomp, irefin, rerr, ermx, tmwork,
     *              reaft6, ddouble, succes, maxmsh, 
     *              numbig, nummed,
     *          r4,amg, stab_cond,ckappa1,gamma1,ckappa,stiff_cond,
     *          itcond,itcondmax)
      implicit double precision (a-h,o-z)
      dimension ltol(ntol)
      dimension fixpnt(*), tol(*), rhs(*)
      dimension xx(*), u(nudim,*), xxold(*)
      dimension uold(ncomp,*), usave(ncomp,*)
      dimension ihcomp(*), irefin(*)
      dimension rerr(ncomp,*), ermx(*), tmwork(*)
      logical reaft6, ddouble, succes, maxmsh
      logical adjrer
      logical stab_cond, stiff_cond
      dimension amg(*), r4(*)
     
      logical pdebug, use_c, comp_c
      common/algprs/ nminit, pdebug, iprint, idum, uval0, use_c, comp_c

*  blas: dcopy

      parameter (eight = 8.0d+0)

*  nmpt is the standard number of mesh points to be added to selected
*  intervals when the mesh is altered based on the distribution
*  in size of the rhs vector.
      parameter ( nmpt = 15)

*  Non-convergence of 6th order.

      if (iprint .eq. 1) write(6,901)
      succes = .false.
      maxmsh = .false.
  
*  NB: the problem must be nonlinear.  Linear problems will either
*  fail to converge for 4th order, or else, once they've converged 
*  for 4th order, must converge for 6th order.

*  Restore the u array to the previous 4th order solution.

      call matcop(ncomp, nudim, ncomp, nmold, uold, u)

      if (reaft6 .and. iprint .ge. 0) write(6,9999)
 9999  format(1h ,'in fail6, reaft6is true')
      if (.not.reaft6 .and. iprint .ge. 0) write(6,9998)
 9998  format(1h ,'in fail6, not reaft6')
      if (ddouble .and. iprint .ge. 0) write(6,9997)
 9997  format(1h ,'in fail6, ddouble is true')
      if (.not.ddouble .and. iprint.ge.0 ) write(6,9996)
 9996  format(1h ,'in fail6, not ddouble')
c no possibility of richardson extrapolation error test 
c       reaft6=.false.
      if (.not. reaft6 .or. .not. ddouble) then

*  Here, either 
*  (1) the mesh for which this 6th order solution failed
*  is not a doubled version of the immediately preceding mesh, or
*  (2) for the immediately preceding mesh, it is not true
*  that the 4th order converged and the 6th order failed.

*  Setting reaft6 to .true. signals that Richardson extrapolation 
*  may be possible if the next 6th order solution fails.  When 
*  reaft6 is true, the routine conv4 immediately sets onto6 to true.

         reaft6 = .true.

*  Save the current 4th order solution in usave.

         call matcop(nudim, ncomp, ncomp, nmsh, u, usave)

*  Use the distribution of sizes in the rhs vector to refine the mesh.
         

  
         call mshref(ncomp, nmsh, nlbc, ntol, ltol,
     *              iorder, rhs, tmwork,
     *              nmax, xx, nmold, xxold, ddouble, maxmsh,
     *              numbig, nummed,amg,stab_cond,stiff_cond,
     *              r4,nfxpnt, fixpnt,irefin,itcond,itcondmax)
         if (.not. maxmsh) call interp(ncomp, nmsh, xx, nudim, u, 
     *                        nmold, xxold, uold)
      
      else

*  Here, reaft6 and ddouble are both true.  So for two consecutive
*  meshes, the 4th order converged and the 6th order failed,
*  and the second mesh is the double of the first.

*  Calculate an error estimate from Richardson extrapolation
*  with the current and previous 4th order solutions.
*  (usave is the 4th order solution saved from the previous (halved) 
*  mesh.)
*  Set addrer to .true. to signal that the rerr array should
*  be adjusted.

         adjrer = .true.

         call rerrvl( ncomp, nmsh, nudim, u, usave, ntol, ltol,
     *          rerr, remax, itlmx, adjrer )
         if (iprint.eq.1) write(6,9994)
 9994        format(1h ,'***in fail6')
         if (iprint.eq.1) write(6,9993) remax, eight*tol(itlmx)
 9993        format(1h ,'remax',1pe14.4,5x,'8*tol',1pe14.4)
         if (remax .lt. eight*tol(itlmx)) then
            succes = .true.
         else

*  Richardson extrapolation did not give sufficient accuracy.
*  Perform selective mesh refinement on the OLD (NB: old!) mesh
*  and the old (saved) solution, using the error estimate from 
*  Richardson extrapolation to guide where the mesh points are placed.
cf controllare se posso mettere 2 
            nmsh = 1 + (nmsh-1)/2
            call dcopy(nmsh, xxold, 2, xx, 1)
            ipow = 4 

*  The rerr array is overwritten by selmsh.
        if (use_c) then
          if (.not. stiff_cond ) then
             call selmsh(ncomp, nmsh, ntol, ltol, tol,
     *             nfxpnt, fixpnt, ipow, nmax,
     *             xx, ncomp, usave, rerr, irefin, ihcomp,
     *             nmold, xxold, ermx, ddouble, maxmsh)

          else
*  The rerr array is overwritten by selconderrmsh.
           call selconderrmsh(ncomp, nmsh, ntol, ltol, tol,
     *             nfxpnt, fixpnt, ipow, nmax,
     *             xx, ncomp, usave, rerr, irefin, ihcomp,
     *             nmold, xxold, ermx, ddouble, maxmsh,
     *             r4,amg,stab_cond)

          end if 
        else
                 call selmsh(ncomp, nmsh, ntol, ltol, tol,
     *             nfxpnt, fixpnt, ipow, nmax,
     *             xx, ncomp, usave, rerr, irefin, ihcomp,
     *             nmold, xxold, ermx, ddouble, maxmsh)
        end if   
 
*  If ddouble is false on exit from selmsh, the call to selmsh has
*  produced a different (non-doubled) mesh.   Interpolate the
*  saved solution (from the old mesh) onto the mesh newly created
*  by selmsh.
*  NB: Because ddouble is false, we won't try Richardson extrapolation
*  if the next 4th order converges and 6th order fails.
            
            if (.not. maxmsh) then
               if (.not. ddouble) then
                  call interp(ncomp, nmsh, xx, nudim, u, nmold, 
     *                    xxold, usave)
                else

*  Selective mesh refinement based on the old mesh simply 
*  produced the same mesh that we started with.  So now refine
*  starting with the doubled mesh (from before) and the solution.

                  reaft6 = .true.
*  Save the solution in usave in case we can carry out Richardson
*  extrapolation in the same circumstances next time.
                  call matcop(nudim, ncomp, ncomp, nmsh, u, usave)

*  Use the distribution of sizes in the rhs vector to refine the mesh.

                  call mshref(ncomp, nmsh, nlbc, ntol, ltol,
     *                   iorder, rhs, tmwork,
     *                   nmax, xx, nmold, xxold, ddouble, maxmsh,
     *                   numbig, nummed,amg,stab_cond,stiff_cond,
     *                   r4,nfxpnt, fixpnt,irefin,itcond,itcondmax)
               
                  if (.not. maxmsh)
     *                call interp(ncomp, nmsh, xx, nudim, u, 
     *                        nmold, xxold, usave)


*              end of logic for needing to refine (again) based on the 
*              current mesh
               endif

*           end of logic for not too many mesh points
            endif

*        end of logic for failure of Richardson extrapolation
*        to produce a converged solution
         endif
      
*     end of logic for both reaft6 and ddouble being true
      endif

      return
 901   format(1h ,'fail6')
      end



  
      subroutine conv8( ncomp, nmsh, ntol, ltol, tol,
     *              nfxpnt, fixpnt, linear, nmax,
     *              xx, nudim, u, def, def6, def8, uold, 
     *              ihcomp, irefin, ermx, err6,
     *              etest8, strctr, 
     *              ddouble, nmold, xxold, maxmsh, succes, first8,
     *              r4, amg,stab_cond,ckappa1,gamma1,ckappa,stiff_cond,
     *                       rpar,ipar)

      implicit double precision (a-h,o-z)
      dimension rpar(*),ipar(*)
      dimension ltol(ntol), tol(ntol)
      dimension fixpnt(*)
      dimension etest8(ntol)
      dimension xx(*), u(nudim,*), def(ncomp,*)
      dimension def6(ncomp,*), def8(ncomp,*), uold(ncomp,*)
      dimension ihcomp(*), irefin(*)
      dimension ermx(*), xxold(*), amg(*), r4(*)
      logical linear, strctr, ddouble, maxmsh, succes, first8
      logical stab_cond, stiff_cond

      logical pdebug, use_c, comp_c
      common/algprs/ nminit, pdebug, iprint, idum, uval0, use_c, comp_c

      intrinsic max

      logical errok

      parameter (one = 1.0d+0, fourth = 0.25d+0, quan8 = 0.5d+0)
      parameter ( efact  = 100.0d+0, huge = 1.0d+30 )

*  blas: dload

      save er6old, er8old

*  The Newton iteration converged for the 8th order solution.
       
    
      if (iprint .eq. 1) write(6,901) 

c      if (first8) then
c         er6old = huge
c         er8old = huge
c         first8 = .false.
c      endif

      if (.not. linear) then
         call dload(ntol, one, etest8, 1)
      else
         do 10 i = 1, ntol
            etest8(i) = one/max(quan8, tol(i)**fourth)
 10             continue
      endif
      
      
      succes = .false.
      maxmsh = .false.
      
*  Check estimated error.  For a nonlinear problem, all components
*  of etest8 (the ratios used in testing the error) are set to one.  
*  For a linear problem, the components of etest8 are in general
*  larger than one.  But if strctr is .true. and the number of mesh
*  points decreased, we set the elements of etest8 to one (which 
*  makes a stricter test).

      if (linear .and. strctr .and. nmsh .lt. nmold)
     *   call dload(ntol, one, etest8, 1)

      call dload(ntol, 10d0, etest8, 1)

      call errest (ncomp, nmsh, ntol, ltol, tol,
     *          nudim, u, uold, etest8, err8, errok)

      if (errok)  then
         succes = .true.
         return
      endif

c      write(6,*) ' err8', err8, nmsh, stiff_cond, stab_cond
c      if ( (use_c .and. err8 .le. 5e2 )) then
c        do im=1,ncomp
c         do jm = 1,nmsh
c              def8(im,jm)=abs(u(im,jm)-uold(im,jm))
c            def8(im,jm)=max(abs(def8(im,jm)),
c     *       abs(u(im,jm)-uold(im,jm)))           
c           enddo
c        enddo
      
c      elseif ( .not. use_c .and. err8 .le. 5e2 ) then
c       do im=1,ncomp
c         do jm = 1,nmsh
c           def8(im,jm)=abs(u(im,jm)-uold(im,jm))
c           def8(im,jm)=max(abs(def8(im,jm)),abs(u(im,jm)-uold(im,jm)))
c           enddo
c         enddo
c      endif
*  At this point, the 8th order solution converged, but did not
*  satisfy the test for termination.
c      if (pdebug) write(6,902) err6, err8, er6old, er8old
c      if (nmsh .lt. nmold. and.
c     *         err6 .gt. efact*er6old .and.
c     *         err8 .gt. efact*er8old) then 

*  If the number of mesh points decreased and the errors in the
*  6th and 8th order solutions did not decrease sufficiently compared
*  to the previous mesh, double the mesh and go back to try to
*  calculate a 4th order solution.

         
c         call dblmsh (nmsh, nmax, xx, nmold, xxold, maxmsh) 
c         if (.not. maxmsh) then
c            er6old = err6
c            er8old = err8
c If the problem is not linear we use the old solution 
c instead the the first value
c old code
c            call initu(ncomp, nmsh, xx, nudim, u,rpar,ipar)
c new code
c            if (linear) then 
c               call initu(ncomp, nmsh, xx, nudim, u,rpar,ipar)
c            else
c we do not use u but uold
c               call matcop(nudim, ncomp, ncomp, nmold, u, uold)
c               call interp(ncomp, nmsh, xx, nudim, u,
c     *                       nmold, xxold, uold)
c            endif
c         endif
c         return
c      endif

*   Here, we know that
*   (1) the number of mesh points exceeds that for the previous mesh; or
*   (2) the number of mesh points decreased, the 6th and 8th order 
*       errors did not satisfy the termination tests, but they did not
*       increase so much that the mesh needed to be doubled.
         
c      er6old = err6
c      er8old = err8
c      if (err8 .le. err6) then
 
*  Perform selective mesh refinement based on the 8th order deferred
*  corrections.  The value of ipow indicates that the error estimate
*  is of order 6.  Then, for a nonlinear problem, interpolate the
*  latest solution onto the new mesh.

         ipow = 6

*  NB: The array def8 will be overwritten by selmsh.
       
      if (use_c) then 
          if (.not. stiff_cond ) then
            call selmsh(ncomp, nmsh, ntol, ltol, tol,
     *            nfxpnt, fixpnt, ipow, nmax,
     *            xx, ncomp, uold, def8, irefin, ihcomp,
     *            nmold, xxold, ermx, ddouble, maxmsh)
          else
*  NB: The array def8 will be overwritten by selconderrmsh.
 
            call selconderrmsh(ncomp, nmsh, ntol, ltol, tol,
     *            nfxpnt, fixpnt, ipow, nmax,
     *            xx, ncomp, uold, def8, irefin, ihcomp,
     *            nmold, xxold, ermx, ddouble, maxmsh,r4,amg,stab_cond)
          end if 
      else
            call selmsh(ncomp, nmsh, ntol, ltol, tol,
     *            nfxpnt, fixpnt, ipow, nmax,
     *            xx, ncomp, uold, def8, irefin, ihcomp,
     *            nmold, xxold, ermx, ddouble, maxmsh)
      end if
        
      if (.not. maxmsh) then
            if (linear) then 
               call initu(ncomp, nmsh, xx, nudim, u,rpar,ipar)
            else
               call matcop(nudim, ncomp, ncomp, nmold, u, uold)
               call interp(ncomp, nmsh, xx, nudim, u,
     *                       nmold, xxold, uold)
            endif
      endif

      return


c      else 
         
*  err8 is greater than err6

*  For a linear problem, set all elements of etest8 to one, 
*  which makes the error test stricter.  (The elements of etest8
*  may have already been set to one earlier in this routine.)

         if (linear) call dload(ntol, one, etest8, 1)

*  Selectively refine the mesh using the old solution and the
*  6th order deferred correction.  Then, for a nonlinear prpblem,
*  interpolate the old solution onto the new mesh.

         ipow = 4

*  The array def6 will be overwritten by selmsh.
      if (use_c) then 
         if (.not. stiff_cond ) then
          call selmsh(ncomp, nmsh, ntol, ltol, tol,
     *            nfxpnt, fixpnt, ipow, nmax,
     *            xx, ncomp, uold, def6, irefin, ihcomp,
     *            nmold, xxold, ermx, ddouble, maxmsh)
         else
            call selconderrmsh(ncomp, nmsh, ntol, ltol, tol,
     *            nfxpnt, fixpnt, ipow, nmax,
     *            xx, ncomp, uold, def6, irefin, ihcomp,
     *        nmold, xxold, ermx, ddouble, maxmsh,r4,amg,stab_cond)

         end if
      else
          call selmsh(ncomp, nmsh, ntol, ltol, tol,
     *            nfxpnt, fixpnt, ipow, nmax,
     *            xx, ncomp, uold, def6, irefin, ihcomp,
     *            nmold, xxold, ermx, ddouble, maxmsh)
      end if
      if (.not. maxmsh) then
            if (linear) then
               call initu(ncomp, nmsh, xx, nudim, u,rpar,ipar)
            else
               call interp(ncomp, nmsh, xx, nudim, u,
     *                       nmold, xxold, uold)
            endif
      endif
*     end of logic for err8 greater than err6
c      endif

      if (pdebug .and. .not.succes) write(6,903)
      return

 901  format(1h ,'conv8')
 902  format(1h ,'err6, err8, er6old, er8old',4(1pe11.3))
 903  format(1h ,'8th order fails error tests.')
      end




      subroutine fail8(ncomp, nmsh, nfxpnt, fixpnt, nmax,
     *      ntol, ltol, tol, nmold,
     *      xx, nudim, u, def6, xxold, uold,
     *      ihcomp, irefin, ermx, ddouble, maxmsh,
     *       r4,amg, stiff_cond, stab_cond)

      implicit double precision(a-h,o-z)
      dimension fixpnt(*), ltol(ntol), tol(ntol)
      dimension xx(*), u(nudim,*), def6(ncomp,*)
      dimension xxold(*), uold(ncomp,*)
      dimension ihcomp(*), irefin(*)
      dimension ermx(*)
      dimension amg(*), r4(*)
      logical ddouble, maxmsh
      logical stiff_cond, stab_cond

      logical pdebug, use_c, comp_c
      common/algprs/ nminit, pdebug, iprint, idum, uval0, use_c, comp_c

      if (pdebug) write(6,901)

*  8th order solution did not converge (the problem must be nonlinear)

      ipow = 4

*  Selectively refine the mesh based on the 6th order deferred
*  correction and the old solution.

*  The def6 array is overwritten by selmsh.
    
        if ( use_c .and. stiff_cond) then 
         call selconderrmsh(ncomp, nmsh, ntol, ltol, tol,
     *            nfxpnt, fixpnt, ipow, nmax,
     *            xx, ncomp, uold, def6, irefin, ihcomp,
     *      nmold, xxold, ermx, ddouble, maxmsh,r4,amg,stab_cond)
        else
         call selmsh(ncomp, nmsh, ntol, ltol, tol,
     *        nfxpnt, fixpnt, ipow, nmax,
     *        xx, ncomp, uold, def6, irefin, ihcomp,
     *        nmold, xxold, ermx, ddouble, maxmsh)
        end if

*  Interpolate to obtain the new initial solution.

      if (.not. maxmsh) then
         call interp(ncomp, nmsh, xx, nudim, u, nmold, xxold, uold)
      endif

      return
 901   format(1h ,'fail8')
      end

      subroutine dccal( ncomp, nmsh, ntol, ltol,
     *                     defexp, dfctol, fval,
     *                     dfexmx, incmp, inmsh, intol,
     *                     derivm)

      implicit double precision  (a-h,o-z)

      dimension  ltol(ntol)
      dimension  defexp(ncomp,nmsh-1), defimp(ncomp,nmsh-1)
      dimension  fval(ncomp,nmsh), ratdc(nmsh-1)

      logical pdebug, use_c, comp_c
      common/algprs/ nminit, pdebug, iprint, idum, uval0, use_c, comp_c

      parameter ( zero = 0.0d+0, one = 1.0d+0 )
      parameter ( rtst = 50.0d+0, tstrat = 0.1d+0 )
      
      intrinsic abs, max

*  blas: idamax

*  Find dfexmx, the maximum-magnitude element of defexp
*  in components for which a tolerance is specified.
*  The component index of dfexmx (incmp), its mesh
*  interval index (inmsh), and its tolerance index (intol),
*  are output parameters.

      dfexmx = zero
      do 10 it = 1, ntol
         icmp = ltol(it)
         idmx = idamax(nmsh-1, defexp(icmp, 1), ncomp)
         dval = abs(defexp(icmp, idmx))
         if (dval .ge. dfexmx) then
            dfexmx = dval
            incmp = icmp
            inmsh = idmx
            intol = it
         endif
 10       continue

      if (pdebug) then
         write(6,901)
         write(6,902) dfexmx, incmp, inmsh, intol
      endif

*  Find derivm (maximum-magnitude element of fval(incmp,*)) 
*  for all mesh points.

      idmx = idamax(nmsh, fval(incmp, 1), ncomp)
      derivm = abs(fval(incmp, idmx))
      if (pdebug) write(6,903) derivm

      return
c     do not want to calculate variables corresponding to implicit dc

*  For component incmp, go through the mesh intervals to calculate
*  (1) dfimmx, the maximum implicit deferred correction;
*  (2) two crucial ratios, rat1 and rat2, used in deciding whether 
*      to refine the mesh; 
*  (3) the array ratdc of deferred-correction ratios (explicit to
*      implicit).
 
*  In defining rat1 and rat2, we consider only intervals for 
*  which the explicit deferred correction (defexp) exceeds the 
*  tolerance dfctol in magnitude.  If it does not, the associated 
*  interval does not affect rat1 or rat2, and the value of ratdc 
*  is taken as 1.

*  rat2 is the maximum-magnitude ratio of sufficiently large
*  defexp to the larger of defimp and dfctol.

*  rat1 is the maximum-magnitude ratio of sufficiently large
*  defexp to the larger in magnitude of defimp and dfctol, but only 
*  for those values of defexp greater than tstrat*dfexmx. 
*  Thus by construction rat1 is less than or equal to rat2.
 
      rat1 = zero
      rat2 = zero
      dfimmx = zero
      smtest = tstrat*dfexmx

      do 100 im = 1, nmsh-1
         texp = defexp(incmp, im) 
         timp = defimp(incmp, im)
         dfimmx = max(dfimmx, abs(timp))
         abtexp = abs(texp)
         if (abtexp .le. dfctol) then
            ratdc(im) = one
         else
            if (abs(timp) .lt. dfctol) timp = dfctol 
            ratdc(im) = texp/timp
            abrat = abs(ratdc(im))
            rat2 = max(rat2, abrat)
            if (abs(texp) .ge. smtest 
     *                .and. abrat .ge. rat1)  rat1 = abrat
         endif
*         if (pdebug) write(6,904) im, texp, timp, ratdc(im), dfctol
 100      continue

      if (pdebug) write(6,905) rat1, rat2, dfimmx
     
      return

 901  format(1h ,'dccal')
 902  format(1h ,'dfexmx, incmp, inmsh, intol',1pe11.3,3i5)
 903  format(1h ,'derivm',1pe11.3)
 904  format(1h ,'im, texp, timp, ratdc, dfctol',i5,4(1pe11.3))
 905  format(1h ,'rat1, rat2, dfimmx', 3(1pe11.3))
      end


      subroutine decid4(linear, dfexmx, 
     *     derivm, dfold, tolval, 
     *     onto6, smooth, callrt, strctr, oscchk, ddouble, reposs)

      
      implicit double precision (a-h,o-z)

      logical linear, onto6, smooth, callrt, strctr, 
     *     oscchk, ddouble, reposs

      logical pdebug, use_c, comp_c
      common/algprs/ nminit, pdebug, iprint, idum, uval0, use_c, comp_c

      logical stest

      parameter ( tenth = 0.1d+0, one = 1.0d+0, two = 2.0d+0 )
      parameter ( thrtwo = 32.0d+0 )
      parameter ( rtst = 50.0d+0, derval = 50.0d+0 )

*  decid4 evaluates information about the deferred corrections
*  and the nature of the problem, and sets various logical
*  flags that control the subsequent logic of the algorithm.

*  This code has been written for clarity, and NOT for efficiency.

      onto6 = .true.
      callrt = .false.
      smooth = .false.
      oscchk = .false.
      strctr = .false.
      reposs = .false.
      ddouble = .false.

*  rat2 is always greater than or equal to rat1.

      if (pdebug) write(6,901)
      if (pdebug) write(6,902) tolval, rtst

      stest = .true.
      if (linear) stest = dfexmx .lt. tenth*dfold

    
      
c      if (rat2 .lt. rtst) then
c         if (stest) then
c            smooth = .true.
c         else
c            oscchk = .true.
c         endif
c         return
c      endif

*  We know now that rat2 .ge. rtst.

      thttol = thrtwo*tolval
      if (pdebug) write(6,903) thttol
      
c      if (rat1 .lt. rtst .and. dfexmx .lt. thttol) then
c          if (stest) then
c             smooth = .true.
c          else
c             oscchk = .true.
c          endif
c          return
c      endif

c      if (rat1 .lt. rtst .and. dfexmx .ge. thttol) then
c         callrt = .true.
c         return
c      endif

*  We know now that rat1 .ge. rtst (and that rat2 .ge. rtst).

      if (derivm .gt. derval .and. dfexmx .lt. thttol) then
         if (stest) then
            smooth = .true.
         else
            oscchk = .true.
         endif
         return
      endif

      return
c      if (derivm .gt. derval .and. dfexmx .gt. thttol) then
c         if (dfimmx .lt. one) then
c            callrt = .true.
            
c         else
c            strctr = .true.
c            if (linear) then
c               onto6 = .false.
c               if (two*rat1 .ge. oldrt1) ddouble = .true.
*           end of logic for linear
c            endif
*        end of logic for dfimmx .ge. one
c         endif
        
c         return
*     end of logic for derivm .gt. derval .and dfexmx .gt. thttol
c      endif

*  To reach this point in the code, both of the following must hold:
*    rat1 .ge. rtst (which means that rat2 .ge. rtst)
*    derivm .le. derval

*  On linear problems, a special termination rule is tried if two 
*  conditions hold:
*    (1) the 4th order solution has been computed on two consecutive 
*        meshes, one of which is the double of the other 
*        (this holds when ddouble is .true.), and 
*    (2) on both meshes, rat1 .ge. rtst, and derivm is small.  When 
*        the conditions in (2) hold for a particular mesh, decid4 
*        sets reposs to .true. to indicate that Richardson 
*        extrapolation may be possible.
*  This set of tests is to take care of transients kept out by 
*  initial conditions.
      
c      if (linear) reposs = .true.
c      return

 901   format(1h ,'decid4')
 902    format(1h ,'tolval, rtst',2(1pe11.3))
 903     format(1h ,'thttol',1pe11.3)
      end


      Subroutine Dfexcl (Ncomp, Nmsh, Xx, Nudim, U, Def8,Def6, Linear,
     *                 Fval, Tmp, Fsub, Dfsub, Df, Ip, Dhold,
     *                 Ntol, Ltol, Tol,JC,rpar,ipar)

      Implicit Double Precision (A-H, O-Z)
      Integer Ncomp, Nmsh
      dimension rpar(*),ipar(*)
      Dimension Xx(Nmsh), U(Nudim,Nmsh), Fval(Ncomp, Nmsh)
      Dimension Def6(Ncomp,Nmsh-1), Def8(Ncomp,Nmsh-1),Tmp(Ncomp,8)
      Dimension Df(Ncomp,Ncomp), Ltol(Ntol), Tol(Ntol)
      Dimension Ip(2*Ncomp), Dhold(2*Ncomp,2*Ncomp)
      Dimension St1(200), St2(200), St3(200)
      External Fsub, Dfsub

      Parameter ( One = 1.0d+0, Two = 2.0d+0 )

      
      Common /Cons1/ A21,A22,A23,A24,A31,A32,A33,A34,C1,C2,
     +       C16,C26,C123,C223,C14,C24
      Common /Algprs/nminit,pdebug,Iprint,idum,Uval0,use_c,comp_c
      Common /Flags/ Ifinal,Iback,Iprec
      
      Logical Linear,pdebug,use_c,comp_c

*  Given The Nmsh Mesh Points Xx, The Estimated Solution
*  U And The Array Fval Of Function Values At (Xx(Im), U(*,Im)),
*  Im = 1,...,Nmsh, Dfexcl Calculates Sixth-Order Implicit
*  Deferred Correction, Stored In The Array Def6, Indexed
*  Over The Components And Mesh Intervals.

*  The Array Tmp Is Workspace For 4 Intermediate Vectors Of
*  Dimension Ncomp.

      Do 90 Im = 1, Nmsh-1
         Hmsh = Xx(Im+1) - Xx(Im)                                              
         C16h = C16/Hmsh 
         C26h = C26/Hmsh
         Do 10 Ic = 1, Ncomp
            Fvim = Fval(Ic,Im)
            Fvim1 = Fval(Ic,Im+1)
            Uim = U(Ic,Im)
            Uim1 = U(Ic,Im+1)
            Tmp(Ic,3) = C16h*(Uim1-Uim)+C123*Fvim1+C14*Fvim
            Tmp(Ic,4) = C26h*(Uim1-Uim)+C223*Fvim1+C24*Fvim
c
c     Put cubic Hermite approximations to the unknowns in
c     Tmp(Ic,3) and Tmp(Ic,4).
c
            St1(Ic) = (Uim+Uim1)/Two
            St2(Ic) = A21*Fvim + A24*Fvim1
            St3(Ic) = A31*Fvim + A34*Fvim1
 10                        Continue

         Xxc1 = Xx(Im)+C1*Hmsh
         Xxc2 = Xx(Im)+C2*Hmsh

         Do 60 Nit = 1, 10

         Do 20 Ic = 1,Ncomp
            Tmp3 = Tmp(Ic,3)
            Tmp4 = Tmp(Ic,4)
            Tmp(Ic,1)  = St1(Ic) + Hmsh*(St2(Ic) + A22*Tmp3 + A23*Tmp4)
            Tmp(Ic,2)  = St1(Ic) + Hmsh*(St3(Ic) + A32*Tmp3 + A33*Tmp4)
 20                        Continue

         Call Fsub (Ncomp,Xxc1,Tmp(1,1),Tmp(1,5),rpar,ipar)
         Call Fsub (Ncomp,Xxc2,Tmp(1,2),Tmp(1,6),rpar,ipar)

         Call Dfsub (Ncomp,Xxc1,Tmp(1,1),Df(1,1),rpar,ipar)

         Do 30 I = 1, Ncomp
            Tmp(I,5) = Tmp(I,5)-Tmp(I,3)
            Tmp(I,6) = Tmp(I,6)-Tmp(I,4)
            Do 25 J = 1,Ncomp
               Dfij = Hmsh*Df(I,J)
               Dhold(I,J) = -A22*Dfij
               Dhold(I,J+Ncomp) = -A23*Dfij
 25         Continue
 30      Continue

         Call Dfsub(Ncomp,Xxc2,Tmp(1,2),Df,rpar,ipar)
         Do 35 I = 1, Ncomp
            Do 32 J = 1, Ncomp
               Dfij = Hmsh*Df(I,J)
               Dhold(I+Ncomp,J) = -A32*Dfij
               Dhold(I+Ncomp,J+Ncomp) = -A33*Dfij
 32         Continue
 35      Continue

         Do 40 I = 1,Ncomp
           Dhold(I,I) = Dhold(I,I) + One
           Dhold(I+Ncomp,I+Ncomp) = Dhold(I+Ncomp,I+Ncomp) + One
 40      Continue

         Call Lufac(2*Ncomp,2*Ncomp,Dhold,Ip,Ier)
         Call Lusol(2*Ncomp,2*Ncomp,Dhold,Ip,Tmp(1,5),Tmp(1,7))

         Do 45 I = 1,Ncomp
            Tmp(I,3) = Tmp(I,3) + Tmp(I,7)
            Tmp(I,4) = Tmp(I,4) + Tmp(I,8)
 45      Continue
           
         Jc = 0
         If (Linear) Goto 70

         
         Do 50 I = 1, Ntol
           Ii = Ltol(I)
           Er = Tol(I)/Hmsh
           If (Abs(Tmp(Ii,7)) .gt. Er*Max(One,Abs(Tmp(Ii,3)))  .or.
     *         Abs(Tmp(Ii,8)) .gt. Er*Max(One,Abs(Tmp(Ii,4)))) Jc = 1
 50      Continue

         If (Jc .eq. 0) Goto 70

 60      Continue

         if (iprint.eq.1) write(6,75)
 75      format(1x,'no convergence of corrections')

         Return

 70      Continue

         Do 80 Ic = 1, Ncomp
            Def6(Ic,Im) = (Hmsh/12.d+0)*(Fval(Ic,Im)+
     *              5.d+0*(Tmp(Ic,3)+Tmp(Ic,4))+Fval(Ic,Im+1))-
     *              U(Ic,Im+1)+U(Ic,Im)
c           write(6,*) 'Def6(',Ic,',',Im,') =', Def6(Ic,Im)
 80      Continue
      do 85 ic=1,ncomp
      tmp(ic,5)=def6(ic,im)
      tmp(ic,6)=def6(ic,im)
 85         continue
      call lusol(2*ncomp,2*ncomp,dhold,Ip,tmp(1,5),tmp(1,7))
      def8(1,Im)=tmp(1,7)
      def8(2,im)=tmp(2,7)
 90   Continue

      Return

 900  Format(/,'** Warning - Possibly Approaching Machine Precision ',
     *         'Beyond Epsilon  = ',D10.3,/)

      End


      subroutine expl (ncomp, nmsh, xx, nudim, u, defexp, fval,
     *                fsub,Im,rpar,ipar)

      implicit double precision (a-h, o-z)
      integer ncomp, nmsh
      dimension rpar(*),ipar(*)
      dimension xx(nmsh), u(nudim,nmsh), fval(ncomp, nmsh)
C     SCMODIFIED: increased the number of dimensions
C      dimension t1(2),t2(2),t3(2),t4(2)
      dimension t1(ncomp),t2(ncomp),t3(ncomp),t4(ncomp)
      dimension defexp(Ncomp)
      external fsub

      logical pdebug,use_c,comp_c
      common/algprs/ nminit, pdebug, iprint, idum, uval0,use_c,comp_c
      parameter ( half = 0.5d+0, fourth = 0.25d+0, thfrth= 0.75d+0 )
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
        Do 10 Ic=1,Ncomp
            t1(ic) = (a5*u(ic, im+1) + b5*u(ic, im))
     *         + hmsh*(c5*fval(ic,im)- d5*fval(ic,im+1))
            t2(ic) = (b5*u(ic,im+1) + a5*u(ic,im))
     *         + hmsh*(-c5*fval(ic,im+1) + d5*fval(ic,im))

 10     continue
        call fsub (ncomp, xx(im) + fourth*hmsh, t1,
     *          t3,rpar,ipar)
        call fsub (ncomp, xx(im) + thfrth*hmsh, t2,
     *          t4,rpar,ipar)

      Do 20 Ic=1,Ncomp
            t1(ic) = half*(u(ic,im+1) + u(ic,im))
     *          + e5*hmsh*(fval(ic,im+1) - fval(ic,im))
     *          - f5*hmsh*(t4(ic) - t3(ic))
 20   continue

         call fsub(ncomp, half*(xx(im) + xx(im+1)), t1,
     *          t2,rpar,ipar)
      do 30 Ic=1,Ncomp
            defexp(ic) = hmsh*(a6*(fval(ic,im+1) + fval(ic,im))
     *           + b6*(t3(ic) + t4(ic)) + c6*t2(ic))
     *           - u(ic,im+1) + u(ic,im)
      Au=Max(abs(u(Ic,Im)),abs(u(ic,Im+1)))
      au=0.0D+0
      Defexp(ic)=defexp(ic)/max(1.0D+0,Au)

 30   continue
 78   continue
 9000 continue
      return
      end



      Subroutine Df8cal (Ncomp, Nmsh, Xx, Nudim, U, Fval, Def8, Linear,
     *                   Tmp, Fsub, Dfsub, Df, Ip, Dhold, Ntol,
     *                   Ltol, Tol,JC,rpar,ipar)

*   Given The Mesh Points Xx, The Solution U, And The Function
*   Values Fval, Df8cal Computes Eighth-Order Deferred Corrections,
*   Which Are Stored In Def8.
*   The Array Tmp Is Workspace For 9 Intermediate Vectors.

      Implicit Double Precision (A-H, O-Z)
      Integer Ncomp, Nmsh
      dimension rpar(*),ipar(*)
      Dimension Xx(Nmsh), U(Nudim,Nmsh), Fval(Ncomp,Nmsh)
      Dimension Def8(Ncomp, Nmsh-1), Tmp(Ncomp,12)
      Dimension Df(Ncomp,Ncomp), Ltol(Ntol), Tol(Ntol)
      Dimension Ip(3*Ncomp), Dhold(3*Ncomp,3*Ncomp)
      Dimension St1(200), St2(200), St3(200), St4(200)
      External Fsub, Dfsub

      
      Parameter ( One = 1.0d+0, Two = 2.0d+0 )

      Common/Cons2/ A21,A22,A23,A24,A25,A31,A32,A34,A35,A41,A42,
     +      A43,A44,A45,B1,B2,B3,C1,C2,C3,C16,C26,C36,C123,C223,
     +      C323,C14,C24,C34
      logical pdebug,use_c,comp_c
      Common /Algprs/nminit,pdebug,Iprint,idum,Uval0,use_c,comp_c
      Common /Flags/ Ifinal,Iback,Iprec

      Logical Linear

      Do 110 Im = 1, Nmsh-1
         Hmsh = Xx(Im+1) - Xx(Im)

         C16h = C16/Hmsh
         C26h = C26/Hmsh
         C36h = C36/Hmsh

         Do 10 Ic = 1, Ncomp
             Fvim = Fval(Ic,Im)
             Fvim1 = Fval(Ic,Im+1)
             Uim = U(Ic,Im)
             Uim1 = U(Ic,Im+1)
             Tmp(Ic,4) = C16h*(Uim1-Uim)+C123*Fvim1+C14*Fvim
             Tmp(Ic,5) = C26h*(Uim1-Uim)+C223*Fvim1+C24*Fvim
             Tmp(Ic,6) = C36h*(Uim1-Uim)+C323*Fvim1+C34*Fvim
             St1(Ic) = (Uim+Uim1)/Two
             St2(Ic) = A21*Fvim + A25*Fvim1
             St3(Ic) = A31*Fvim + A35*Fvim1
             St4(Ic) = A41*Fvim + A45*Fvim1
 10                          Continue

         Xxc1 = Xx(Im)+C1*Hmsh
         Xxc2 = Xx(Im)+C2*Hmsh
         Xxc3 = Xx(Im)+C3*Hmsh

         Do 80 Nit = 1, 10

            Do 20 Ic = 1, Ncomp
              Tmp4 = Tmp(Ic,4)
              Tmp5 = Tmp(Ic,5)
              Tmp6 = Tmp(Ic,6)
              Tmp(Ic,1) = St1(Ic) + Hmsh*(St2(Ic) + A22*Tmp4 + A23*Tmp5
     +             + A24*Tmp6)
              Tmp(Ic,2) = St1(Ic) + Hmsh*(St3(Ic) + A32*Tmp4 + A34*Tmp6)
              Tmp(Ic,3) = St1(Ic) + Hmsh*(St4(Ic) + A42*Tmp4 + A43*Tmp5
     +             + A44*Tmp6)
 20                               Continue

            Call Fsub (Ncomp,Xxc1,Tmp(1,1),Tmp(1,7),rpar,ipar)
            Call Fsub (Ncomp,Xxc2,Tmp(1,2),Tmp(1,8),rpar,ipar)
            Call Fsub (Ncomp,Xxc3,Tmp(1,3),Tmp(1,9),rpar,ipar)

            Call Dfsub(Ncomp,Xxc1,Tmp(1,1),Df,rpar,ipar)
            Do 30 I = 1, Ncomp
               Tmp(I,7) = Tmp(I,7)-Tmp(I,4)
               Tmp(I,8) = Tmp(I,8)-Tmp(I,5)
               Tmp(I,9) = Tmp(I,9)-Tmp(I,6)
               Do 25 J = 1,Ncomp
                  Dfij = Hmsh*Df(I,J)
                  Dhold(I,J) = -A22*Dfij
                  Dhold(I,J+Ncomp) = -A23*Dfij
                  Dhold(I,J+2*Ncomp) = -A24*Dfij
 25            Continue
 30         Continue

            Call Dfsub(Ncomp,Xxc2,Tmp(1,2),Df,rpar,ipar)
            Do 40 I = 1, Ncomp
                Do 35 J = 1, Ncomp
                   Dfij = Hmsh*Df(I,J)
                   Dhold(I+Ncomp,J) = -A32*Dfij
                   Dhold(I+Ncomp,J+Ncomp) = 0.d+0
                   Dhold(I+Ncomp,J+2*Ncomp) = -A34*Dfij
 35             Continue
 40         Continue

            Call Dfsub(Ncomp,Xxc3,Tmp(1,3),Df,rpar,ipar)
            Do 50 I = 1, Ncomp
                 Do 45 J =  1, Ncomp
                    Dfij = Hmsh*Df(I,J)
                    Dhold(I+2*Ncomp,J) = -A42*Dfij
                    Dhold(I+2*Ncomp,J+Ncomp) = -A43*Dfij
                    Dhold(I+2*Ncomp,J+2*Ncomp) = -A44*Dfij
 45              Continue
 50         Continue

            Do 60 I = 1, Ncomp
               Dhold(I,I) = Dhold(I,I) + One
               Dhold(I+Ncomp,I+Ncomp) = Dhold(I+Ncomp,I+Ncomp) + One
               Dhold(I+2*Ncomp,I+2*Ncomp) =
     *              Dhold(I+2*Ncomp,I+2*Ncomp) + One
 60         Continue

            Call Lufac(3*Ncomp,3*Ncomp,Dhold,Ip,Ier)
            Call Lusol(3*Ncomp,3*Ncomp,Dhold,Ip,Tmp(1,7),Tmp(1,10))

            Do 65 I = 1,Ncomp
               Tmp(I,4) = Tmp(I,4) + Tmp(I,10)
               Tmp(I,5) = Tmp(I,5) + Tmp(I,11)
               Tmp(I,6) = Tmp(I,6) + Tmp(I,12)
 65         Continue

            Jc = 0
            If (Linear) Goto 90

            Do 70 I = 1, Ntol
               Ii = Ltol(I)
               Er = Tol(I)/Hmsh
               If (Abs(Tmp(Ii,10)) .gt. Er*Max(One,Abs(Tmp(Ii,4))) .or.
     *           Abs(Tmp(Ii,11)) .gt. Er*Max(One,Abs(Tmp(Ii,5))) .or.
     *           Abs(Tmp(Ii,12)) .gt. Er*Max(One,Abs(Tmp(Ii,6)))) Jc = 1
 70                                 Continue

            If (Jc .eq. 0) Goto 90

 80         Continue

         if (iprint.eq.1) write(6,8930)
 8930    format(1x,'no convergence of 8th order defcors')

      Return
 90   Continue

         Do 100 Ic = 1, Ncomp
           Def8(Ic,Im) = Hmsh*(B1*(Fval(Ic,Im)+Fval(Ic,Im+1))+
     *                   B2*(Tmp(Ic,4)+Tmp(Ic,6))+B3*Tmp(Ic,5))-
     *                   U(Ic,Im+1)+U(Ic,Im)
 100     Continue

 110     Continue

      Return

 900  Format(/,'** Warning - Possibly Approaching Machine Precision ',
     *         'Beyond Epsilon  = ',D10.3,/)

      End


      Subroutine Stcon1
      Implicit Double Precision(A-H,O-Z)
      Common/Cons1/A21,A22,A23,A24,A31,A32,A33,A34,C1,C2,
     +      C16,C26,C123,C223,C14,C24

      Parameter ( One = 1.0d+0, Two = 2.0d+0,  Three = 3.0d+0 )
      Parameter ( Four = 4.0d+0, Five = 5.0d+0,  Six = 6.0d+0 )

      Rt5 = Sqrt(5.0d+0)

      A21 = (Six + Rt5)/120.0d0
      A22 = -Rt5/120.0d0
      A23 = (-13.d0*Rt5)/120.0d0
      A24 = (-Six + Rt5)/120.d0

      A31 = (Six-Rt5)/120.0d0
      A32 = (13.0d0*Rt5)/120.0d0
      A33 = Rt5 / 120.0d0
      A34 = (-Six - Rt5)/120.d0

      C1 = (Five - Rt5)/10.0d0
      C2 = (Five + Rt5)/10.0d0

      C12 = C1*C1
      C22 = C2*C2

      C16 = Six*(C1 - C12)
      C26 = Six*(C2 - C22)

      C123 = Three*C12 - Two*C1
      C223 = Three*C22 - Two*C2

      C14 = One - Four*C1 + Three*C12
      C24 = One - Four*C2 + Three*C22

      Return
      End

      Subroutine Stcon2
      Implicit Double Precision(A-H,O-Z)
      Common/Cons2/ A21,A22,A23,A24,A25,A31,A32,A34,A35,A41,A42,
     +      A43,A44,A45,B1,B2,B3,C1,C2,C3,C16,C26,C36,C123,C223,
     +      C323,C14,C24,C34

      Parameter ( One = 1.0d+0, Two = 2.0d+0, Three = 3.0d+0 )
      Parameter ( Four = 4.0d+0, Six = 6.0d+0 )

      Rt21 = Sqrt(21.0d+0)

      A21 = One/28.d0 + Three*Rt21/1960.d0
      A22 = -Rt21/280.d0
      A23 = -32.d0*Rt21/735.d0
      A24 = -23.d0*Rt21/840.d0
      A25 = -One/28.d0 + Three*Rt21/1960.d0

      A31 = One/64.d0
      A32 = 7.d0*Rt21/192.d0
      A34 = -7.d0*Rt21/192.d0
      A35 = -One/64.d0

      A41 = One/28.d0 - Three*Rt21/1960.d0
      A42 = 23.d0*Rt21/840.d0
      A43 = 32.d0*Rt21/735.d0
      A44 = Rt21/280.d0
      A45 = -(One/28.d0) - Three*Rt21/1960.d0

      B1 = One/20.0d0
      B2 = 49.0d0/180.d0
      B3 = 16.0d0/45.d0

      C1 = One/Two - Rt21/14.d0
      C2 = One/Two
      C3 = One/Two + Rt21/14.d0

      C12 = C1*C1
      C22 = C2*C2
      C32 = C3*C3

      C16 = Six*(C1 - C12)
      C26 = Six*(C2- C22)
      C36 = Six*(C3 - C32)

      C123 = Three*C12 - Two*C1
      C223 = Three*C22 - Two*C2
      C323 = Three*C32 - Two*C3

      C14 = One - Four*C1 + Three*C12
      C24 = One - Four*C2 + Three*C22
      C34 = One - Four*C3 + Three*C32

      Return
      End

      subroutine dfimcl( ncomp, nmsh, defexp, chold, dsq, dexr,
     *                 ipivot, defimp)
      implicit double precision (a-h,o-z)
      dimension defexp(ncomp, nmsh-1), chold(ncomp, ncomp, nmsh-1)
      dimension dsq(ncomp, ncomp), dexr(ncomp)
      dimension ipivot(ncomp), defimp(ncomp, nmsh-1)

      logical pdebug, use_c, comp_c
      common/algprs/ nminit, pdebug, iprint, idum, uval0, use_c, comp_c

*  blas: dcopy

      parameter (zero = 0.0d+0)

*  dfimcl calculates the rational deferred correction array,
*  which is indexed over the components and mesh intervals.

      call mtload(ncomp, nmsh-1, zero, ncomp, defimp)

      do 100 im = 1, nmsh-1

         call dcopy(ncomp, defexp(1,im), 1, dexr(1), 1 )
         do 50 ic = 1, ncomp
            call dcopy(ncomp, chold(1,ic,im), 1, dsq(1,ic), 1)
 50             continue
         call lufac (ncomp, ncomp, dsq, ipivot, ierlu)
         if (ierlu .eq. 0) then
            call lusol (ncomp, ncomp, dsq, ipivot, dexr, 
     *          defimp(1, im))
         endif 
 100      continue
      return
      end


      subroutine osc (ncomp, nmsh, dfexmx, incmp, 
     *     defcor, ratdc, ddouble, inmsh, onto6, trst6, smooth) 

      implicit double precision (a-h, o-z)
      dimension defcor(ncomp, *), ratdc(*)

      logical ddouble, onto6, trst6, smooth

      logical pdebug, use_c, comp_c
      common/algprs/ nminit, pdebug, iprint, idum, uval0, use_c, comp_c

      intrinsic abs

      parameter ( zero = 0.0d+0, half = 0.5d+0 )
      parameter ( frac1 = 0.1d+0, frac2 = 1.0d-2 )
  
*  For linear problems, subroutine osc performs heuristic tests  
*  to detect an oscillating solution.  the tests check whether 
*  (significant) explicit and implicit deferred corrections have 
*  different (componentwise) signs.

*  dfexmx is the maximum-magnitude explicit deferred correction,
*  and is known to occur in component incmp.

*  The array defcor contains the explicit deferred corrections.
*  The array ratdc contains the ratios of explicit to implicit 
*  deferred corrections.

*  jsndif counts the number of differences in sign.
            
      jsndif = 0
      rmax = zero
      
*  allsum is the sum of the magnitudes of all deferred corrections,
*  smlsum is the sum of the magnitudes of the small deferred
*  corrections, and bigsum is the sum of the magnitudes of the
*  large deferred corrections.  Here, small is defined as less 
*  than half of the maximum.

      ninter = nmsh - 1

      if (pdebug) write(6,901)

      allsum = zero
      smlsum = zero 
      bigsum = zero
      ibig = 0
      ism = 0

      do 30 im = 1, ninter
         abdef = abs(defcor(incmp,im))
         allsum = allsum + abdef
         if (abdef .lt. half*dfexmx) then
              ism = ism + 1 
              smlsum = smlsum + abdef
         else 
              ibig = ibig + 1 
              bigsum = bigsum + abdef
         endif

*  The counter of sign differences is incremented if (1) ratdc is negative 
*  (which means that the two deferred corrections have opposite 
*  sign) and (2) the explicit deferred correction is not too small
*  relative to the maximum.

         if (pdebug) write(6,902) im, ratdc(im), abdef, frac2*dfexmx

         if (ratdc(im).lt.zero .and. abdef.ge.frac2*dfexmx) then
            jsndif = jsndif + 1

*  If more than 4 sign differences have occurred, exit after setting
*  ddouble to .true., which signals that the mesh
*  should be doubled (i.e., twice as many intervals).

            if (jsndif.gt.4) then
               onto6 = .false.
               ddouble = .true. 
               return
            endif
            if (abs(ratdc(im)).ge.rmax) then
               rmax = abs(ratdc(im))
               inmsh = im
            endif
         endif
 30       continue

      if (pdebug) write(6,903) rmax, jsndif
      
      avsm = zero  
      if (ism.gt.0) avsm = smlsum/ism 
      avbg = zero
      if (ibig.gt.0) avbg = bigsum/ibig
      ave = allsum/ninter

      if (pdebug) write(6,904) ave, avsm, avbg

      if (avsm.gt.frac1*avbg .or. ave.gt.half*avbg) then

*  The error appears to be uniformly large.
*  Signal that the 6th order solution should be calculated.
         onto6 = .true.

      elseif (jsndif.eq.0) then
*  If there were no sign changes, the problem appears to be smooth.
         smooth = .true.
         onto6 = .true.
      else

*  If the sign changed at between 1 and 4 points, don't go on to
*  6th order, and don't ever accept a 6th order solution even if the
*  error estimate at a later stage indicates that it is OK to do so.
*  Set ddouble to .false., to signal that the mesh will not necessarily
*  be doubled.

          ddouble = .false.
          onto6 = .false.
          trst6 = .false.
      endif 
      return
 901   format(1h ,'osc')
 902    format(1h ,'im, ratdc, abdef, val',i5,3(1pe11.3))
 903     format(1h ,'rmax, jsndif', 1pe11.3,i5)
 904      format(1h ,'ave, avsm, avbg', 3(1pe11.3))
      end 


      subroutine ratcor ( ncomp, nmsh, xx, defimp, bhold, dfrat)
      implicit double precision (a-h, o-z)
      dimension xx(nmsh), defimp(ncomp,nmsh-1) 
      dimension dfrat(ncomp,nmsh-1), bhold(ncomp,ncomp,nmsh-1)

*  blas: ddot, dscal

      parameter (zero = 0.0d+0, half = 0.5d+0, one = 1.0d+0)

      ninter = nmsh - 1
      do 10 im = 1, ninter
         hmsh = xx(im+1) - xx(im)
         call dscal(ncomp*ncomp, (-half*hmsh), bhold(1,1,im), 1)     
 10       continue

      do 20 im = 1, ninter
      do 20 ic = 1, ncomp
         bhold(ic,ic,im) = bhold(ic,ic,im) + one
 20       continue

      do 30 im = 1, ninter
      do 30 ic = 1, ncomp
         dfrat(ic,im) = ddot(ncomp, bhold(ic,1,im), ncomp,
     *                        defimp(1,im), 1) 
 30       continue
      return
      end 


      subroutine fixjac(ncomp, nmsh, nlbc,  
     *    iorder, ntol, ltol, tol, 
     *    xx, nudim, u, defcor, defnew, delu, 
     *    rhs, fval, utrial, rhstri, 
     *    rnsq, uint, ftmp, tmprhs,
     *    ajac, topblk, botblk, ipivot,
     *    fsub, gsub, iflag,rpar,ipar)

* Fixed Jacobian iterations.

      implicit double precision (a-h,o-z)
      dimension rpar(*),ipar(*)
      dimension  ltol(ntol), tol(ntol)
      dimension  xx(nmsh), u(nudim,nmsh), defcor(ncomp,nmsh-1) 
      dimension  defnew(ncomp,nmsh-1), delu(ncomp,nmsh)
      dimension  rhs(ncomp*nmsh), fval(ncomp,nmsh)
      dimension  utrial(ncomp,nmsh), rhstri(ncomp*nmsh)
      dimension  uint(ncomp), ftmp(ncomp), tmprhs(ncomp*nmsh)
      dimension  ajac(ncomp, 2*ncomp, nmsh-1)
      dimension  topblk(nlbc,ncomp), botblk(ncomp-nlbc,ncomp)
      dimension  ipivot(ncomp*nmsh)
      logical    better

      external   fsub, gsub

      logical pdebug, use_c, comp_c
      common/algprs/ nminit, pdebug, iprint, idum, uval0, use_c, comp_c
      common/mchprs/flmin, flmax, epsmch

      intrinsic  abs, max 

*  blas: dcopy, dssq

      parameter ( one    = 1.0d+0 )
      parameter ( xlarge = 1.0d+6, huge = 1.0d+20, lmtfrz = 8)
      parameter ( rngrow = 16.0d+0, rfact = 100.0d+0 )
      parameter ( tolfct = 0.1d+0 )

*  The iteration scheme uses a fixed Jacobian matrix to solve for 
*  correction vectors, once there has been convergence of the Newton
*  iterations on this mesh.   It is assumed that the LU
*  factors of the Jacobian have been left unaltered since
*  their calculation.

      if (iprint .eq. 1) write(6,901)
      ninter = nmsh - 1
      rnold = flmax
      isize=nmsh*ncomp

*  Evaluate the right-hand side rhstri at the initial solution u by
*  adding the new deferred corrections to the already-calculated
*  rhs vector.

      call dcopy(nlbc, rhs, 1, rhstri, 1)
      ind = nlbc
      do 10 im = 1, ninter
      do 10 ic = 1, ncomp
         ind = ind + 1
         rhstri(ind) = rhs(ind) + defnew(ic, im)
 10       continue
      ind = ninter*nmsh + nlbc + 1
      call dcopy(ncomp-nlbc, rhs, 1, rhstri, 1)

      call dssq  ( nmsh*ncomp, rhstri, 1, scale, sumsq )
      rnsq = (scale**2)*sumsq

      iter = 0
      
*  If the initial right-hand side is too large, do not even attempt to 
*  solve the nonlinear equations.

      if (rnsq.gt.huge .or. 
     *      (iorder.eq. 8 .and. rnsq.gt.xlarge)) then
         if (iprint .eq. 1) write (6,902) rnsq
         iflag = -2
         return
      end if
      call dcopy(ncomp*nmsh, rhstri, 1, rhs, 1)

*  Statement 100 is the top of the iteration loop.

 100   continue

*  If rnsq is sufficiently small, terminate immediately.

      if (iprint .eq. 1) write(6,903) iter, rnsq
      if (rnsq .le. epsmch) then
         iflag = 0
         return
      endif

      iter = iter + 1

*  Solve for the step delu by solving a system involving the fixed 
*  Jacobian (whose LU factors are saved).  Copy the rhs array into 
*  tmprhs, which is then overwritten by blkslv.

      call dcopy(ncomp*nmsh, rhs, 1, tmprhs, 1) 
      call dcopy(ncomp*nmsh,tmprhs,1,delu,1)
      job = 0
      call crslve(Topblk, Nlbc, Ncomp, Ajac, Ncomp, 2*Ncomp,
     + Ninter, Botblk, Ncomp-Nlbc, Ipivot, Delu, job)

*  Compute the trial point utrial by adding delu to u.

      call matcop( nudim, ncomp, ncomp, nmsh, u, utrial )
      call maxpy( ncomp, nmsh, one, delu, ncomp, utrial )

*  compute the right-hand side vector rhstri and its squared
*  two-norm at the trial point.

      rnold = rnsq
      call fneval(ncomp, nmsh, xx, ncomp, utrial, fval, fsub,rpar,ipar)
      call rhscal (ncomp, nmsh, nlbc, xx, ncomp, utrial, defcor,
     *   fsub, gsub, rhstri, rnsq, fval, ftmp, uint,rpar,ipar) 

*  If rnsq strictly decreased, update the solution vector u 
*  and the right-hand side rhs.

      better = .false.
      if (rnsq .lt. rnold) then
         better = .true.
         call matcop( ncomp, nudim, ncomp, nmsh, utrial, u )
         call dcopy( ncomp*nmsh, rhstri, 1, rhs, 1 )
      endif

*  Stop the fixed Jacobian iterations if there have been too
*  many iterations, or if rnsq has not decreased by a factor
*  of at least rngrow.

      if (iter .ge. lmtfrz .or. rnsq .gt. (rnold/rngrow)) then
         if (better) then

*  Setting iflag to -3 signals that, although the fixed Jacobian
*  iterations did not succeed, the current point was an improvement
*  on the previous one.  Hence, if we switch to a Newton procedure,
*  the right-hand side does not need to be recalculated.

            iflag = -3
         else
            iflag = -2
         endif
         if (iprint .eq. 1) write(6,904) iflag
         return
      endif

*  Test for convergence using the ratio abs((change in u)/max(u,1)).

      do 150 im = 1, nmsh
      do 150 it = 1, ntol
         itol = ltol(it)
         er = abs(delu(itol,im))/max(abs(u(itol,im)), one)
         if (er .gt. tolfct*tol(it)) go to 100
 150      continue

*  To exit from the loop here, the convergence tests have
*  been passed.

      if (iprint .ge. 0) write(6,905) iter, rnsq

      iflag = 0
      return
 901  format(1h ,'fixed Jacobian iterations')
 902  format(1h ,'Large residual, rnsq =',1pe12.4)
 903  format(1h ,'iter, rnsq',i5,1pe11.3)
 904  format(1h ,'failure of fixed Jacobian, iflag =',i5)
 905  format(1h ,'fixed Jacobian convergence',i5,1pe11.3)
      end 

************************************************************************

      subroutine lineq( ncomp, nmsh, nlbc, 
     *    ludone, xx, nudim, u, defcor, 
     *    delu, rhs, fval, 
     *    uint, ftmp, dftmp1, dftmp2, dgtm, tmprhs,
     *    ajac, topblk, botblk, bhold, chold, ipivot,
     *    fsub, dfsub, gsub, dgsub, iflag,rpar,ipar)
***********************************************************************
*     call by: bvpsol
*     calls to: lnrhs, jaccal, dcopy, colrow, dload, dload, clrslve
*           maxpy
**********************************************************************
      implicit double precision (a-h,o-z)
      dimension rpar(*),ipar(*)
      dimension  xx(nmsh), u(nudim, nmsh), defcor(ncomp, nmsh-1)
      dimension  delu(ncomp, nmsh), rhs(ncomp*nmsh)
      dimension  fval(ncomp,nmsh), uint(ncomp), ftmp(ncomp)
      dimension  topblk(nlbc,*), botblk(ncomp-nlbc,*)
      dimension  dftmp1(ncomp, ncomp), dftmp2(ncomp, ncomp)
      dimension  dgtm(ncomp), tmprhs(ncomp*nmsh)
      dimension  ajac(ncomp, 2*ncomp, nmsh-1),
     *               bhold(ncomp, ncomp, nmsh-1),
     *               chold(ncomp, ncomp, nmsh-1)
      dimension  ipivot(ncomp*nmsh)
      integer    job
      logical    ludone
      external   fsub, dfsub, gsub, dgsub

      logical pdebug, use_c, comp_c
      common/algprs/ nminit, pdebug, iprint, idum, uval0, use_c, comp_c

*  blas: dcopy, dload

      parameter  ( one = 1.0d+0, zero = 0.0d+0 )

*  The routine lineq calculates the Newton step for a linear
*  problem.  The Newton step is exact unless the Jacobian
*  matrix is singular.

      isize=nmsh*ncomp
      ninter = nmsh - 1

      if (.not. ludone) then

*  Compute the right-hand side vector rhs.
         iflag = 0 
         call lnrhs (ncomp, nmsh, nlbc, xx, nudim, u, 
     *          fsub, gsub, rhs, rnsq, fval, ftmp, uint,rpar,ipar) 

*  If the Jacobian for this mesh has not previously been
*  calulated and factorized successfully, call jaccal.
*  The block-structured Jacobian matrix is stored in three 
*  matrices (topblk, ajac, and botblk).
*  The matrices bhold and chold are also calculated in jaccal,
*  and are saved for later use in outer routines.

         call jaccal (ncomp, nmsh, nlbc, 
     *      xx, nudim, u, fval, dgtm, dftmp1, dftmp2, uint,
     *      ajac, topblk, botblk, bhold, chold,
     *      dfsub, dgsub,rpar,ipar)

*  Call blkdcm to calculate the LU factors of the Jacobian.
*  The factors are overwritten on the matrices topblk, ajac and botblk.
*  Interchanges are represented in the integer array ipivot.

      call dcopy(ncomp*nmsh,rhs,1,tmprhs,1)
      call dcopy(ncomp*nmsh,tmprhs,1,delu,1)
     
      job = 0
      call colrow(isize, topblk,nlbc, ncomp, ajac, ncomp, 2*ncomp,
     +   ninter, botblk, ncomp-nlbc, ipivot, delu, iflag, job)

         ludone = .true.

*  Copy the rhs into the temporary vector tmprhs, which will be
*  overwritten by blkslv.


      else

*  The right-hand side is the deferred correction array,
*  padded with zeros at the boundary conditions.
        iflag = 0
         call dload(nlbc, zero, tmprhs(1), 1)
         do 100 im = 1, ninter
            loc = (im-1)*ncomp + nlbc + 1
            call dcopy(ncomp, defcor(1,im), 1, tmprhs(loc), 1)
 100            continue
         nrhs = ninter*ncomp + nlbc + 1
         call dload(ncomp-nlbc, zero, tmprhs(nrhs), 1)
         call dcopy(ncomp*nmsh,tmprhs,1,delu,1)
      job = 0
      
      call crslve(topblk,nlbc,ncomp,ajac,ncomp,2*ncomp,ninter,
     +   botblk,ncomp-nlbc,ipivot,delu,job)

      endif

*  Since the problem is linear, the Newton step  is exact.  The 
*  new u array is obtained by adding delu to u.

      call maxpy ( ncomp, nmsh, one, delu, nudim, u )
     
c      iflag = 0
      return

 901   format(1h ,'Singular matrix')
      end

******************************************************************
      subroutine newteq(ncomp, nmsh, nlbc,  
     *    rhsgiv, ntol, ltol, tol, 
     *    xx, nudim, u, defcor, 
     *    delu, rhs, fval,
     *    utrial, rhstri, 
     *    uint, ftmp, dftmp1, dftmp2, dgtm, tmprhs, xmerit, 
     *    ajac, topblk, botblk, bhold, chold, ipivot,
     *    fsub, dfsub, gsub, dgsub, iter, iflag,
     *    isign,rpar,ipar,frscal)
******************************************************************
*     call by: bvpsol
*     calls to:
******************************************************************
      implicit double precision (a-h,o-z)
      dimension rpar(*),ipar(*)
      dimension  ltol(*), tol(*), xx(*)
      dimension  fval(ncomp,*)
      dimension  u(nudim, *), delu(ncomp, *), utrial(ncomp,nmsh)
      dimension  rhs(ncomp*nmsh),  defcor(ncomp,*)
      dimension  ajac(ncomp, 2*ncomp, *)
      dimension  topblk(nlbc,*), botblk(ncomp-nlbc,*)
      dimension  ftmp(*), uint(*), dgtm(ncomp)
      dimension  dftmp1(ncomp, ncomp), dftmp2(ncomp, ncomp)
      dimension  bhold(ncomp, ncomp, *), chold(ncomp, ncomp, *)
      dimension  ipivot(*), isign(*)
      dimension  rhstri(ncomp*nmsh)
      dimension  tmprhs(ncomp*nmsh), xmerit(ncomp, nmsh)
      logical rhsgiv
    

      external   fsub, dfsub, gsub, dgsub

      parameter  ( zero   = 0.0d+0, one    = 1.0d+0 )
      parameter ( two = 2.0d+0, half = 0.5d+0, fourth = 0.25d+0 )
      parameter ( tenth = 0.1d+0, ten = 10.0d+0, hund = 100.0d+0 )
      parameter (itcondmax = 5)      

      logical pdebug, use_c, comp_c
      common/algprs/ nminit, pdebug, iprint, idum, uval0, use_c, comp_c
      common/mchprs/flmin, flmax, epsmch

      logical gtpdeb, imprvd, braktd, crampd, extrap, vset, wset
      save  gtpdeb, mfsrch, epsaf, epsag, eta, rmu, tolabs, alfmax
      save  tolrel, toltny
      
      intrinsic  abs, max

*  blas: dcopy
  
      logical frscal
c      save frscal
      parameter (cnvfct = 0.1d+0 )
c      data  frscal/.true./
c frscal has been updated in bvpsol
c
      data  alfsml/1.0d-4/,  alfmax/1.1d+0/
      data  imerit/1/, lmtnwt/39/
      data  shrfct/100.0d+0/, stpfct/2.0d+0/
      data  gtpdeb/.false./, mfsrch/5/
      data  eta/.999999d+0/, rmu/1.0d-6/

    
*  The routine newteq performs Newton iterations with a line
*  search, to solve the nonlinear equations.


*  Set up constants if this is the first call to newteq.
 
      if (frscal) then
         frscal = .false.
         epsaf = epsmch
         epsag = epsmch
         tolabs = epsmch
         tolrel = epsmch
         toltny = epsmch
      endif
      ninter = nmsh - 1
    
      if (iprint .eq. 1) write(6,901) 
  
*  A Newton method with line search and watchdog safeguarding
*  is used to try to solve the nonlinear equations.

*  Initialize iter (the counter of Newton iterations) and alfold
*  (the step taken at the previous iteration).

      iter = -1
      alfold = one
      alfa = zero

      

      if (.not. rhsgiv) then

*  If necessary, evaluate the right-hand side at the initial u. 
         
         call rhscal (ncomp, nmsh, nlbc, xx, nudim, u, defcor,
     *      fsub, gsub, rhs, rnsq, fval, ftmp, uint,rpar,ipar) 
      endif

*  At any given Newton iteration, rnprev is the value of rnsq at 
*  the immediately preceding Newton iteration.

      rnprev = flmax
      rnbest = flmax
      if (.not. pdebug .and. iprint .ge. 0) write (6,902)

*  Initialize counter of watchdog iterations.

      itwtch = 0

*  Statement 100 is the top of the Newton iteration loop.
   
 100   continue

      iter = iter + 1

      if (iprint .eq. 1) write(6,910) iter

*  If there have been too many Newton iterations, terminate.

      if (iter .ge. lmtnwt) then
         if (iprint .ge. 0) write(6,903)
         iflag = -2
         return
      endif

*  The vector rhs is the right-hand side at the current iterate,
*  and rnsq is its squared two-norm.
*  Perform watchdog tests, using the unscaled merit function (rnsq)
*  as the watchdog function.  The routine wtchdg updates rnbest
*  and itwtch.  If iflwat is not zero, this sequence of Newton
*  iterations is terminated.

      iflwat = 0 
    
      call wtchdg ( iter, rnsq, rnbest, rnprev, itwtch, 
     *                alfold, iflwat )
 
      if (iflwat .ne. 0) then
         if (iprint .ge. 0) write(6,904) iter
         iflag = -3
         return
      endif

*  Watchdog tests are passed.  Proceed with the Newton iteration.
*  Call jaccal to evaluate the block-structured Jacobian matrix, 
*  which is stored in three matrices (topblk, ajac, and botblk).
*  The matrices bhold and chold are saved for use in later 
*  calculations in the outer routine.    
      

*  If rnsq is sufficiently small, terminate immediately.
*  Note that the stored Jacobian does not correspond exactly
*  to the final point.

      if (rnsq .le. epsmch) then
         if (iprint .ge. 0)  write(6,906) iter, rnsq
         iflag = 0
         return
      endif

      call jaccal (ncomp, nmsh, nlbc, 
     *    xx, nudim, u, fval, dgtm, dftmp1, dftmp2, uint,
     *    ajac, topblk, botblk, bhold, chold,
     *    dfsub, dgsub,rpar,ipar)

*  blkdcm is called to calculate the LU factors of the Jacobian, 
*  which are overwritten on topblk, ajac and botblk.
*  Interchanges are represented in the integer array ipivot.  

*   Solve for the Newton step delu.  Copy the rhs array into tmprhs,
*   which is then overwritten by blkslv.

      call dcopy(ncomp*nmsh, rhs, 1, tmprhs, 1)
      call dcopy(ncomp*nmsh,tmprhs,1,delu,1)
      job = 0  
      call colrow(Nmsh*ncomp,topblk,nlbc,ncomp,ajac,ncomp,2*ncomp,
     +   ninter,botblk,ncomp-nlbc,ipivot,delu,iflag,job)
  
      if (iprint .ge. 0 .and. iflag.ne.0) write(6,905) iter
      if (iflag .ne. 0) then
         iflag = -1
         return
      end if
*     the jacobian is singular   


*  If imerit = 1, the line search is based on the scaled merit function,
*  the squared two-norm of the solution xmerit of the linear system
*      (Jacobian)*xmerit = rhs, 
*  where (Jacobian) is the Jacobian at the current Newton iterate. 
*  Thus the initial value of the scaled merit function is simply
*  the squared two-norm of the Newton step delu itself. 

      if (imerit.eq.1) then
         call mssq( ncomp, nmsh, delu, xmscal, xmsq )
         fmtry = (xmscal**2)*xmsq
      else
*  The unscaled merit function is simply the squared two-norm of rhs. 
         fmtry = rnsq
      end if

c  fa and oldg represent the merit function and its gradient
c  at the initial point of the line search.

      fa = fmtry
      oldg = -two*fa
      alfa = zero
      if (iprint .eq. 1) write (6,908) alfa, fmtry, rnsq

*  On the first Newton iteration, the initial trial step is unity.   
*  On subsequent iterations, the initial step is not allowed to 
*  be more than the factor stpfct larger than the final step at 
*  the immediately preceding iteration.


      alfa = one
      if(stpfct*alfold .lt. one) alfa = stpfct*alfold

      if (alfa .lt. alfsml) alfa = alfsml

      fmold = fa
      inform = -1

*  Statement 150 is the top of the inner line search iteration.
*  The line search routine getptq has been altered so that it
*  terminates with an indication of success as soon as a
*  strictly lower value of the merit function is found.  Note that
*  this is a much less strict requirement than the usual sufficient
*  decrease conditions.

      
 150   continue
      iwr = 6
      call getptq (gtpdeb, mfsrch, iwr, alfmax, alfsml, alfuzz,
     *      epsaf, epsag, 
     *      eta, fmtry, fmold, oldg, rmu, tolabs, tolrel, toltny,
     *      imprvd, inform, nfsrch, alfa, alfbst, fbest, 
     *      braktd, crampd, extrap, vset, wset, nsamea, nsameb,
     *      alin, blin, fa, factor, fv, fw, xtry, xv, xw)

*  inform = 1, 2 or 3 indicates success in finding an acceptable point.
*  inform = 4 means alfmax is too small (this should never happen here,
*  since alfmax is set always to 1.1).
*  inform = 5 means that a decrease was not achieved for any step
*  greater than alfsml.
*  inform = 6 means a better point could not be found (the minimum
*  probably lies too close to alfa=0).
*  inform = 7 means that the gradient at alfa=0 (oldg) is positive 
*  (this cannot happen here, since oldg=-two*fa, and fa is a non-negative
*  number)

      if (pdebug) write(6,907) inform, alfa

      if (inform .eq. 5) then
          iflag = -5
          return
      elseif (inform .eq. 4 .or. inform .eq. 7) then
         iflag = -4
         return
      elseif (inform .eq. 0) then

*  inform = 0 means that a new function value should be obtained
*  with the step alfa.
*  We may override alfa from getptq by requiring that the step is not 
*  allowed to decrease by more than a factor of shrfct during 
*  a line search iteration.
*
         if (alfa .lt. alfold/shrfct) alfa = alfold/shrfct
         alfold = alfa

*  Define the next iterate utrial = u + alfa*delu. 
*  Call fneval and rhscal to evaluate the right-hand side 
*  rhstri at utrial. 
*  The vector rhstri is stored separately, and rhs is overwritten
*  only when an improved point is found.

         call matcop ( nudim, ncomp, ncomp, nmsh, u, utrial)
         call maxpy ( ncomp, nmsh, alfa, delu, ncomp, utrial )
         call fneval(ncomp, nmsh, xx, ncomp, utrial, fval, fsub,
     *        rpar,ipar)
         call rhscal (ncomp, nmsh, nlbc, xx, ncomp, utrial, defcor,
     *      fsub, gsub, rhstri, rnsqtr, fval, ftmp, uint,rpar,ipar) 

         fmold = fmtry
         if (imerit .eq. 1) then

*  Solve a linear system to obtain the 2-d array xmerit whose squared 
*  norm is the scaled merit function.   The LU factors of the Jacobian 
*  have already been calculated by blkdcm.  
*  Copy rhstri into tmprhs, which is overwritten by blkslv. 

      call dcopy(ncomp*nmsh, rhstri, 1, tmprhs, 1)
      call dcopy(ncomp*nmsh,tmprhs,1,xmerit,1)
      job = 0 
      call crslve(topblk, nlbc, ncomp, ajac, ncomp, 2*ncomp, ninter,
     +   botblk, ncomp-nlbc, ipivot, xmerit,job)
      call mssq( ncomp, nmsh, xmerit, xscale, xsolsq )
      fmtry = (xscale**2)*xsolsq
         else

*  The unscaled merit function is the squared two-norm of the right-hand
*  side.
            fmtry = rnsqtr
         end if
         if (iprint .eq. 1) write (6,908) alfa, fmtry, rnsqtr
         go to 150
      endif

*  To reach here, inform must be 1, 2, 3, or 6, and the line search 
*  has found a strictly lower value of the merit function.
*  Store the new Newton iterate in u, and the corresponding rhs
*  vector in rhs.
      
      rnprev = rnsq
      rnsq = rnsqtr
      call matcop (ncomp, nudim, ncomp, nmsh, utrial, u)
      call dcopy(ncomp*nmsh, rhstri, 1, rhs, 1)
      if (iprint .ge. 0) write(6,909) iter, alfa, fmtry, rnsq

*  Now test for convergence using the ratio of the Newton step 
*  for each component with max(1, abs(current solution estimate)).
*  If the test fails for any element of u, branch back to the
*  top of the Newton iteration.

      do 160 im = 1, nmsh
      do 160 it = 1, ntol
         icmp = ltol(it)
         er = abs(delu(icmp,im))/max(abs(u(icmp,im)), one)
         if (er .gt. cnvfct*tol(it)) go to 100
 160      continue
       
  
      if (iprint .ge. 0) write(6, 906) iter+1, rnsq
      iflag = 0

*  To fall through the above loop, the termination test for a 
*  sufficiently small delu is satisfied. 
*  Note that the stored Jacobian and its factorization do not
*  correspond to the final solution. 
 

      return

 901  format(1h ,'start Newton iterations')
 902  format(1h ,' iter',
     *              7x,'alfa',6x,'merit',7x,'rnsq')
 903  format(1h ,'Too many Newton iterations')
 904  format(1h ,'Watchdog tests fail, iter =', i5)
 905  format(1h ,'Singular Jacobian, iter=',i5)
 906  format(1h ,'Convergence, iter =',i5,4x,'rnsq =',1pe12.3)
 907  format(1h ,'inform, alfa after getptq',i5,3x, 1pe11.3)
 908  format(1h ,'alfa, merit, rnsq',3(1pe11.3))
 909  format(1h ,i5,3(1pe11.3))
 910  format(1h ,'Newton iteration',i5)
      end 


      subroutine wtchdg ( iter, wmerit, wmbest, wmprev, 
     *      itwtch, alfold, iflag )

*  Logic for watchdog tests.

      implicit double precision (a-h,o-z)
      parameter ( itonew = 5, itwtmx = 8, grfct = 100.0d+0 )
      parameter ( half = 0.5d+0 )

*  Perform watchdog tests in two forms: 
*  (1) to determine whether a sufficient decrease in the 
*  watchdog merit function has occurred within the most recent
*  sequence of itwtmx iterations;
*  (2) to determine whether the watchdog merit function has increased 
*  too much in a single iteration after itonew Newton iterations 
*  have been performed.  This allows the merit function to increase
*  wildly only during the first itonew iterations.

*  wmbest is the smallest watchdog merit function achieved in this 
*  sequence of Newton iterations.
*  wmprev is the watchdog merit function from the immediately 
*  preceding Newton iteration.

*  itwtch counts the number of iterations without an improvement
*  in the unscaled merit function.

*      write(6,99) iter, wmerit, wmbest, wmprev
*      write(6,98) itwtch, alfold
*   99 format(1h ,'iter,wmer,wbest,wprev',i5,3(1pe15.5))
*   98 format(1h ,'itwtch,alfold',i5,1pe15.5)
      iflag = 0      
      if (wmerit .le. wmbest) then 

*  The current watchdog merit function is the best.

         wmbest = wmerit
         itwtch = 0 
         return
      endif

*  The current merit function is not the best.

      itwtch = itwtch + 1

*  Do not apply watchdog tests if (1) the previous step alfold 
*  exceeds 1/2, or (2) the watchdog merit function decreased in 
*  the immediately preceding iteration and itwtch does not
*  exceed twice its maximum.

      if (alfold .ge. half) return
      if (wmerit .le. wmprev .and. itwtch .le. 2*itwtmx) return


*  If more than itwtmx iterations have occurred without 
*  an overall improvement in the watchdog merit function,
*  signal for termination.

      if (itwtch .ge. itwtmx) then
         iflag = -1

*  If a too-large increase in the watchdog merit function 
*  compared to the best value occurred, and iter .ge. itonew,
*  signal for termination.

      elseif (iter .ge. itonew .and. 
     *          wmerit .gt. grfct*wmbest) then
          iflag = -1
      endif
      return
      end

      subroutine fneval(ncomp, nmsh, xx, nudim, u, fval, fsub,rpar,ipar)
      implicit double precision (a-h,o-z)
      dimension rpar(*),ipar(*)
      dimension xx(nmsh), u(nudim,nmsh), fval(ncomp,nmsh)
      external fsub

*  fneval evaluates the function values (from fsub) for
*  a given mesh xx and array u, and stores the values
*  in the array fval.

      call fsub (ncomp, xx(1), u(1,1), fval(1,1),rpar,ipar)
      do 50 im = 1, nmsh-1
         hmsh = xx(im+1) - xx(im)
         call fsub (ncomp, xx(im+1), u(1,im+1), fval(1,im+1),rpar,ipar)
 50       continue
      return
      end


      subroutine jaccal (ncomp, nmsh, nlbc, 
     *   xx, nudim, u, fval,
     *   dgtm, dftm1, dftm2, uint,
     *   ajac, topblk, botblk, bhold, chold,
     *   dfsub, dgsub,rpar,ipar)

      implicit double precision (a-h,o-z)
      dimension rpar(*), ipar(*)                                                            
      dimension xx(nmsh), u(nudim,nmsh), fval(ncomp,nmsh)
      dimension dgtm(ncomp)
      dimension dftm1(ncomp, ncomp), dftm2(ncomp, ncomp), 
     *             uint(ncomp)
      dimension ajac(ncomp, 2*ncomp, nmsh-1)
      dimension topblk(nlbc, ncomp), botblk(ncomp-nlbc,ncomp)
      dimension bhold(ncomp, ncomp, nmsh-1), 
     *             chold(ncomp, ncomp, nmsh-1)

      external  dfsub, dgsub 

      logical pdebug, use_c, comp_c
      common/algprs/ nminit, pdebug, iprint, idum, uval0, use_c, comp_c

*  blas: dcopy, ddot

      parameter ( zero = 0.0d+0, half = 0.5d+0, eighth = 0.125d+0 )
      parameter ( four = 4.0d+0, six = 6.0d+0 )
      parameter ( one = 1.0d+0, three = 3.0d+0, twelve = 12.0d+0 )


*      if (pdebug) write(6,901)

      ninter = nmsh - 1
      
*      if (pdebug) write(6,902)
      do 110 i = 1, nlbc
         call dgsub (i, ncomp, u(1,1), dgtm,rpar,ipar)
         call dcopy(ncomp, dgtm(1), 1, topblk(i,1), nlbc)
*         if (pdebug) write(6,903) i, (topblk(i,j),j=1,ncomp)
 110      continue

      call dfsub (ncomp, xx(1), u(1,1), dftm1(1,1),rpar,ipar)

*  on entry to jaccal, the array fval contains the function values
*  at (xx(im), u(ic,im)), ic=1,...,ncomp and im = 1,...,nmsh,
*  calculated by a preceding call of rhscal with the same xx and u
*  arrays.

*      if (pdebug) write(6,904)

      do 200 im = 1, ninter

         hmsh = xx(im+1) - xx(im)

         do 120 ic = 1, ncomp
            uint(ic) = half*(u(ic,im) + u(ic,im+1))  
     *         - eighth*hmsh*(fval(ic,im+1) - fval(ic,im))
 120            continue
         xhalf = half*(xx(im+1) + xx(im))
         call dfsub (ncomp, xhalf, uint, dftm2(1,1),rpar,ipar)
         do 140 ic = 1, ncomp
            do 130 jc = 1, ncomp
               dsq = ddot(ncomp, dftm2(ic,1), ncomp,
     *                       dftm1(1,jc), 1)
               ajac(ic,jc,im) = -hmsh*(dftm1(ic,jc)/six
     *             + dftm2(ic,jc)/three + hmsh*dsq/twelve)
 130                  continue
            ajac(ic,ic,im) = ajac(ic,ic,im) - one
*            if (pdebug) write(6,905) im, ic, 
*     *            (ajac(ic,jc,im), jc=1,ncomp)
 140            continue

         call dfsub (ncomp, xx(im+1), u(1,im+1), dftm1(1,1),rpar,ipar)
         do 170 ic = 1, ncomp
            do 160 jc = 1, ncomp
               dsq = ddot(ncomp, dftm2(ic,1), ncomp,
     *                         dftm1(1,jc), 1)
               ajac(ic,jc+ncomp,im) = -hmsh*(dftm1(ic,jc)/six
     *               + dftm2(ic,jc)/three - hmsh*dsq/twelve)
 160                  continue
            call dcopy(ncomp, ajac(ic,ncomp+1,im), ncomp,
     *                   chold(ic,1,im), ncomp)
            call dcopy(ncomp, dftm1(ic,1), ncomp,
     *                   bhold(ic,1,im), ncomp)
            ajac(ic,ic+ncomp,im) = ajac(ic,ic+ncomp,im) + one
            chold(ic,ic,im) = ajac(ic,ic+ncomp,im)
*            if (pdebug) write(6,905) im, ic,
*     *                    (ajac(ic,jc+ncomp,im),jc=1,ncomp)
 170            continue


 200             continue
*      if (pdebug) write(6,906)
      do 220 i = nlbc+1, ncomp
         call dgsub (i, ncomp, u(1, nmsh), dgtm,rpar,ipar) 
         call dcopy(ncomp, dgtm(1), 1, botblk(i-nlbc,1), ncomp-nlbc)
*         if (pdebug) write(6,903) i,(botblk(i-nlbc,j), j=1,ncomp)
 220      continue
c      write(6,991)
 991         format(1x,'topblk')
c      write(6,992) topblk(1,1),topblk(1,2)
 992            format(1x,2g22.10)
c      write(6,993)
 993               format(1x,'main jacobian')
      do 994 iu=1,ninter
      do 994 iv=1,ncomp
c      write(6,996) ajac(iv,1,iu),ajac(iv,2,iu),ajac(iv,3,iu),
c     +ajac(iv,4,iu)
 996        format(1x,4f12.7)
 994           continue
c      write(6,995)
 995              format(1x,'botblk')
c      write(6,992) botblk(1,1),botblk(1,2)


      return

 901  format(1h ,'jaccal')
 902  format(1h ,'topblk')
 903  format(1h ,i5,6(1pe11.3))
 904  format(1h ,'ajac')
 905  format(1h ,2i5,5(1pe11.3))
 906  format(1h ,'botblk')
      end 



      subroutine lnrhs (ncomp, nmsh, nlbc,
     *   xx, nudim, u, fsub, gsub, 
     *   rhs, rnsq, fval, ftmp, uint,rpar,ipar) 

       implicit double precision(a-h,o-z)

*  This subroutine is designed to calculate the right-hand
*  side for linear problems.
      dimension rpar(*), ipar(*)    
      dimension xx(*), u(nudim,*)
      dimension rhs(*), fval(ncomp,*), ftmp(*), uint(*)
      external fsub, gsub 

      parameter ( zero = 0.0d+0, half = 0.5d+0, eighth = 0.125d+0 )
      parameter ( one = 1.0d+0, four = 4.0d+0, six = 6.0d+0 )

      common/mchprs/flmin, flmax, epsmch
      intrinsic abs

*  blas: dssq

      ninter = nmsh - 1
      rnsq = zero

*  first, process the left-hand boundary conditions.

      do 20 i = 1, nlbc
         call gsub (i, ncomp, u(1,1), wg,rpar,ipar)
         rhs(i) = -wg
 20       continue

*  Next, process the interior mesh points.  The fval array
*  contains the function values from fsub at xx and u.

      do 50 im = 1, ninter
         hmsh = xx(im+1) - xx(im)
         do 30 ic = 1, ncomp
            uint(ic) = half*(u(ic,im) + u(ic,im+1))  
     *         - eighth*hmsh*(fval(ic,im+1) - fval(ic,im))
 30             continue
         xhalf = half*(xx(im) + xx(im+1)) 
         call fsub (ncomp, xhalf, uint, ftmp,rpar,ipar)
         loc = (im-1)*ncomp + nlbc
         do 40 ic = 1, ncomp
            rhs(loc+ic) = -u(ic,im+1) + u(ic,im)
     *        + hmsh*
     *            (fval(ic,im) + fval(ic,im+1) + four*ftmp(ic))/six
 40             continue
 50              continue

      nrhs = ninter*ncomp
      do 60 ii = nlbc+1, ncomp
         call gsub (ii, ncomp, u(1,nmsh), wg,rpar,ipar) 
         rhs(nrhs+ii) = -wg
 60       continue

      call dssq  ( nmsh*ncomp, rhs, 1, scale, sumsq )
      rnsq = (scale**2)*sumsq

      return
      end


      subroutine rhscal (ncomp, nmsh, nlbc,
     *   xx, nudim, u, defcor,
     *   fsub, gsub, 
     *   rhs, rnsq, fval, ftmp, uint,rpar,ipar) 

       implicit double precision(a-h,o-z)

*  This subroutine constructs the (ncomp*nmsh)-dimensional
*  vector rhs, which is the right-hand side of the Newton equations.
*  The ncomp by nmsh array fval is assumed to have been calculated
*  elsewhere by routine fneval.
          dimension rpar(*),ipar(*)
      dimension  xx(nmsh), u(nudim,nmsh), defcor(ncomp,nmsh-1) 
      dimension  rhs(ncomp*nmsh), fval(ncomp,nmsh)
      dimension  ftmp(ncomp), uint(ncomp)
      external   fsub, gsub 

      logical pdebug, use_c, comp_c
      common/algprs/ nminit, pdebug, iprint, idum, uval0, use_c, comp_c
      common/mchprs/flmin, flmax, epsmch
      intrinsic abs

*  blas: dssq

      parameter ( zero = 0.0d+0, half = 0.5d+0, eighth = 0.125d+0 )
      parameter ( one = 1.0d+0, four = 4.0d+0, six = 6.0d+0 )

*      if (pdebug) write(6,901)

*  ninter is the number of intervals in the mesh (one less than the
*  number of mesh points)

      ninter = nmsh - 1
      rnsq = zero

*  First, process the left-hand boundary conditions.

      do 20 i = 1, nlbc
         call gsub (i, ncomp, u(1,1), wg,rpar,ipar)
         rhs(i) = -wg
 20       continue

*  Next, process the interior mesh points.  The fval array
*  contains the function values from fsub at xx and u.

      do 50 im = 1, ninter
         hmsh = xx(im+1) - xx(im)
         do 30 ic = 1, ncomp
            uint(ic) = half*(u(ic,im) + u(ic,im+1))  
     *         - eighth*hmsh*(fval(ic,im+1) - fval(ic,im))
 30             continue
         xhalf = half*(xx(im) + xx(im+1)) 
         call fsub (ncomp, xhalf, uint, ftmp,rpar,ipar)
         loc = (im-1)*ncomp + nlbc
         do 40 ic = 1, ncomp
            rhs(loc+ic) = -u(ic,im+1) + u(ic,im) + defcor(ic,im)
     *        + hmsh*
     *            (fval(ic,im) + fval(ic,im+1) + four*ftmp(ic))/six

 40             continue
 50              continue

      nrhs = ninter*ncomp
      do 60 ii = nlbc+1, ncomp
         call gsub (ii, ncomp, u(1,nmsh), wg,rpar,ipar) 
         rhs(nrhs+ii) = -wg
 60       continue

     
      call dssq  ( nmsh*ncomp, rhs, 1, scale, sumsq )
      rnsq = (scale**2)*sumsq
*f
     
      if (pdebug) then
         write (6,902) rnsq
         write(6,903)
         write(6,904) (rhs(i), i=1,ncomp*nmsh)
      endif
      return

 901   format(1h ,'rhscal')
 902    format(1h ,'rnsq',1pe11.3)
 903     format(1h ,'rhs vector')
 904      format(1h ,(7(1pe11.3)))
      end

      subroutine dblmsh (nmsh, nmax, xx, nmold, xxold, maxmsh) 
      implicit double precision (a-h,o-z)
      dimension xx(*), xxold(*)
      logical maxmsh

      logical pdebug, use_c, comp_c
      common/algprs/ nminit, pdebug, iprint, idum, uval0, use_c, comp_c

*  blas: dcopy

      parameter (half = 0.5d+0)


*  This routine is used to double the mesh, i.e., produce a mesh
*  with twice as many intervals in which each new interval is
*  half the corresponding old interval.

*  On entry to dblmsh, the integer nmsh and the array xx
*  specify a set of mesh points xx(1),..., xx(nmsh) (assumed 
*  to be in ascending order).

*  If the number of mesh points in the doubled mesh would
*  exceed the maximum allowed number nmax, the flag maxmsh is
*  set to true, and we exit without changing any other parameters.

*  Otherwise, nmold is set to the old number of mesh points,
*  xxold is set to the old mesh, nmsh is the new number of mesh
*  points, and xx contains the new mesh points.

      nmold = nmsh
      call dcopy(nmold, xx, 1, xxold, 1)

      ninnew = 2*(nmsh-1)
      nmnew = ninnew + 1
      if(nmnew .ge. nmax) then
         if (iprint .ge. 0)  write(6,901) nmnew
         maxmsh = .true.
         return
      endif
      maxmsh = .false.

*  Loop backwards through the old mesh points to create the new ones.

      xx(nmnew) = xx(nmsh)
      do 100 i = ninnew, 4, -2
         id2 = i/2
         xx(i) = half*(xx(i+1) + xx(id2))
         xx(i-1) = xx(id2)
 100      continue

*  Calculate the new xx(2). xx(1) remains unchanged.
      xx(2) = half*(xx(3) + xx(1))
      nmsh = nmnew
      if(iprint .ge. 0)  write(6,902) nmsh
      return
 901  format (1h , ' dblmsh.  maximum mesh exceeded, nmnew =', i8)
 902  format (1h , ' dblmsh.  the doubled mesh has ', i8,' points.') 
      end 


      subroutine selmsh(ncomp, nmsh, ntol, ltol, tol,
     *     nfxpnt, fixpnt, ipow, nmax,
     *     xx, nudim, u, ermeas, irefin, ihcomp,
     *     nmold, xxold, ermx, ddouble, maxmsh)

      implicit double precision (a-h,o-z)

      dimension  ltol(ntol), tol(ntol), fixpnt(*)
      dimension  xx(*), u(nudim, *), ermeas(ncomp,*)
      dimension  irefin(nmsh-1), ihcomp(nmsh-1)
      dimension  xxold(*), ermx(*)
      logical    ddouble, maxmsh

      intrinsic abs, max, int

*  blas: dcopy
*  double precision dlog

      logical pdebug, use_c, comp_c
      common/algprs/ nminit, pdebug, iprint, idum, uval0, use_c, comp_c
      Common /Flags/ Ifinal,Iback,Iprec
      

      parameter  ( zero = 0.0d+0, one = 1.0d+0, onep1 = 1.1d+0 ) 
      parameter  ( erdcid = 5.0d+0 )
      parameter  ( phitst = 0.1d+0 )
  
      logical first
      save    first, rlndec
      data    first / .true. /

*  The routine selmsh performs selective mesh refinement, depending
*  on the error measure ermeas.
   
      if (first) then
         first = .false.
         rlndec = dlog(erdcid)
      endif

      maxmsh = .false.
*      Nref = .true.

      if (pdebug) write(6,901) nmsh, ipow

      frcpow = one/ipow
      ddouble = .false.
      nmold = nmsh
      ninter = nmsh - 1
c      Iprec = Min(Iprec,1)

*  Copy the current mesh into the xxold array.


      call dcopy(nmold, xx, 1, xxold, 1)
      ithres = 0
      thres = one

*  On input, the array ermeas represents some error measure defined
*  over the components and mesh intervals (not mesh points).  
*  It is normalized in the following loop with respect to the
*  tolerance array and the current solution.
*  The value errmax gives the maximum normalized error.

      errmax = zero
      do 120 im = 1, ninter
         ermx(im) = zero
         do 110 it = 1, ntol
            jcomp = ltol(it)
            denom = tol(it)*max(one, abs(u(jcomp,im)))
            ems = ermeas(jcomp,im)
            ermeas(jcomp,im) = abs(ems)/denom
*            if (pdebug .and. ermeas(jcomp,im) .ge. thres)
*     *             write(6,902) im,jcomp,ems,ermeas(jcomp,im)
            err = ermeas(jcomp, im)
            if (err .ge. ermx(im)) then
                ermx(im) = err
                ihcomp(im) = jcomp
            endif
 110            continue
         errmax = max(ermx(im), errmax)
 120      continue


      if (pdebug) write(6,903) errmax
c       write(6,903) errmax
      if (errmax .gt. zero .and. errmax .le. erdcid) then

*  If errmax > 0 and .le. erdcid, find the smallest integer exponent ii 
*  such that (erdcid**ii)*errmax > erdcid.

         if(errmax .gt. one) then
            ii = 1
            decii = erdcid
         else
            ilg = -dlog(errmax)/rlndec
            ii = 2 + ilg
            decii = erdcid**ii
         endif
         
*  Multiply error measures by erdcid**ii.

         errmax = decii*errmax
         do 140 im = 1, ninter
            ermx(im) = decii*ermx(im)
            do 140 it = 1, ntol
               jcomp = ltol(it)
               ermeas(jcomp,im) = decii*ermeas(jcomp,im)
 140               continue
      endif

 200   continue

*  For each interval im,  the integer irefin(im) is calculated 
*  based on an equidistrbution procedure involving the
*  threshold value thres.  If irefin(im) > 1, we add 
*  points to interval im.  If irefin(im) = 1, we check whether
*  point im can be removed from the mesh.

*  nmest is a lower bound on the number of points in the new mesh.
*  We do not know in advance how many points will be removed,
*  so nmest is computed by assuming all eligible points are removed.

      nmest = nmsh
      do 220 im = 1, ninter
         if (ermx(im).ge.thres) then
            irefin(im) = int(ermx(im)**frcpow) + 1
            nmest = nmest + irefin(im) - 1
         else
            irefin(im) = 1
            nmest = nmest - 1
         endif
 220      continue
*      if (pdebug) write(6,904) nmest, (irefin(i), i=1,ninter)

      if (nmest .gt. nmax) then

         go to 360

      elseif (nmest-1 .gt. 3*ninter) then

         call dblmsh (nmsh, nmax, xx, nmold, xxold, maxmsh) 
         ddouble = .true.
         return
      endif

*  It appears that we can perform the desired selective mesh
*  refinement.

*  Now begin running through the mesh, adding and possibly deleting
*  points as indicated by the irefin array.

*  The integer new is a count of the number of intervals in
*  the tentative mesh being generated by the refinement strategy.

      new = 1

*  The first interval is treated as a special case, since xx(1)
*  always remains in the mesh, and cannot be a fixed point.

      rlen = xxold(2) - xx(1)
      slen = rlen
      if (irefin(1).gt.1) then
         dx = rlen/irefin(1)
         do 230 j = 2, irefin(1) 
            new = new + 1
            xx(new) = xx(1) + (j-1)*dx
 230            continue
      endif

*  The fixed points specified by the fixpnt array cannot be
*  removed from the mesh.  The value fxnext indicates the 'next'
*  fixed point.  When no further fixed points remain to be processed
*  (or if nfxpnt = 0), fxnext is set to a value strictly larger than
*  the last mesh point, so that no mesh point can equal fxnext.  
*  This way we need to compare only the values of xxold(i)
*  and fxnext.

      ifxcnt = 1
      if (nfxpnt .eq. 0) then
         fxnext = onep1*abs(xxold(nmsh))
      else
         fxnext = fixpnt(ifxcnt)
      endif

*  jtkout is a counter of the number of consecutive points that 
*  have been removed from the mesh.

      jtkout = 0
      do 330 im = 2, ninter
         rlold = rlen
         rlen = xxold(im+1) - xxold(im)

*  If xxold(im) is the next fixed point, it cannot be removed
*  and so we don't test its error estimates.

         if(xxold(im) .eq. fxnext)  then

            ifxcnt = ifxcnt + 1
            if(ifxcnt .gt. nfxpnt) then
               fxnext = onep1*abs(xxold(nmsh))
            else
               fxnext = fixpnt(ifxcnt)
            endif

         elseif (irefin(im).eq.1) then

*  If xxold(im) is not a fixed point and irefin(im) = 1, we may wish 
*  to remove point im from the mesh.

*  If we are considering removing points and jtkout = 0, this
*  is the first point in a possible consecutive set to be removed,
*  and we initialize phihat, which represents a maximum of
*  certain estimates.
*  If jtkout is not zero, previous points contiguous to this
*  point have been removed, and phihat retains its previous value.

            slen = slen + rlen

            if (jtkout .eq. 0) then
                ind1 = ihcomp(im-1)
                phihat = ermeas(ind1,im-1)/(rlold**ipow)
            endif
            phihat = max(phihat,
     *                 ermeas(ihcomp(im),im)/(rlen**ipow))
            val1 = phihat*(slen**ipow)
            if (val1 .le. phitst 
     *         .and. jtkout .lt. 4 ) then

*  Increment the counter of removed points.
*  'Remove' the mesh point xxold(im) by not including it.

               jtkout = jtkout+1
               go to 330
            endif
*        end of logic for irefin(im) = 1.
         endif

         jtkout = 0 
         new = new + 1
         xx(new) = xxold(im)
         if (irefin(im) .gt. 1) then
            dx = rlen/irefin(im)
            do 300 j = 2, irefin(im)
              new = new + 1
              xx(new) = xxold(im) + (j-1)*dx
 300                 continue
         endif
         slen = rlen
         
         if (new .gt. nmax) then

*  If the new mesh contains too many points, branch out of the
*  loop to try alternative strategies.

            go to 360

         elseif (new .gt. 3*ninter)  then

*  Here, the new mesh does not exceed the specified maximum,
*  but has more than 3 times as many intervals as the old mesh.
*  Try doubling the mesh if possible.

            call dcopy(nmsh, xxold, 1, xx, 1)
            call dblmsh (nmsh, nmax, xx, nmold, xxold, maxmsh) 
            ddouble = .true.
            return
   
          end if
 330       continue

*  To end up here, we have processed the entire interval,
*  and have neither exceeded the specified maximum nor
*  exceeded three times the number of intervals in the old
*  mesh.  The last mesh point remains unchanged.
     
      new = new + 1
      xx(new) = xxold(nmsh)
      nmsh = new
      maxmsh = .false.
      if (iprint .ge. 0) write(6,905) nmsh
      return

 360   continue
*  To reach here, the number of mesh points created at some stage
*  of the refinement process was larger than the maximum permitted 
*  value nmax.  

*  Check whether the mesh can safely be doubled.
  
      if ((2*nmsh-1) .lt. nmax) then

*  Double the mesh.
         call  dblmsh (nmsh, nmax, xx, nmold, xxold, maxmsh) 
         ddouble = .true.

*  If the number of intervals is too large and the mesh cannot be 
*  doubled, increase the threshold thres by a factor of erdcid and
*  try the selective refinement again.  
*  If this happens three times without success or if thres exceeds
*  or is equal to errmax, stop.  (In this case, we know already 
*  that doubling the mesh produces too many points.)

      elseif (thres .lt. errmax .and. ithres .lt. 3) then
         ithres = ithres + 1
         thres = erdcid*thres
         if(thres .gt. errmax) thres = errmax
         call dcopy(nmsh, xxold, 1, xx, 1)
         go to 200
      else
         nmsh = 2*nmsh - 1
         maxmsh = .true.
      endif
      return

 901  format(1h ,'selmsh.  nmsh, ipow =',2i5)
 902  format(1h ,'im, jcomp, ermeas, normalized er',2i5,2(1pe11.3))
 903  format(1h ,'errmax',1pe11.3)
 904  format(1h ,'nmest, irefin',(10i5))
 905  format(1h ,'selmsh.  new nmsh =',i8)
 910  format(1h ,'ihcomp',(10i5))
      end 


       subroutine smpmsh (nmsh, nmax, xx, intref, numadd,
     *      nmold, xxold, maxmsh) 

      implicit double precision (a-h,o-z)
      logical maxmsh
      dimension xx(*), xxold(*)

      logical pdebug, use_c, comp_c
      common/algprs/ nminit, pdebug, iprint, idum, uval0, use_c, comp_c

*  blas: dcopy

*  The routine smpmsh performs simple mesh refinement by adding 
*  points to one or three interval(s) in the region indicated
*  by the integer intref.
*  numadd gives the trial number of points to be added in each
*  interval. 

      if (pdebug) write(6,901)

      nmold = nmsh
      call dcopy(nmold, xx, 1, xxold, 1)

*  numadd is altered if necessary so that it lies between 4 and 49

      if(numadd .gt. 49) then
         numadd = 49
      elseif (numadd .lt. 4) then
         numadd = 4
      endif
      if (pdebug) write (6,902) nmsh, intref, numadd

      maxmsh = .false.
      if (intref .eq. 1) then

*  Add numadd points to the first interval if intref = 1.

         nmnew = nmsh + numadd
         if (nmnew .gt. nmax) then
            if (iprint .ge. 0)  write(6,903) nmnew
            maxmsh = .true.
            return
         endif

*  Renumber the later points in reverse order.

         nint = numadd + 1
         do 100 i = nmnew, numadd+2, -1
            xx(i) = xx(i-numadd)
 100            continue
         dx = (xx(2) - xx(1))/nint
         do 110 i = 2, nint
            xx(i) = xx(1) + (i-1)*dx
 110            continue

      elseif (intref .eq. nmsh-1) then

*  Add numadd points to the last interval if intref = nmsh-1.

         nmnew = nmsh + numadd
         if (nmnew .gt. nmax) then
            if (iprint .ge. 0)  write(6,903) nmnew
            maxmsh = .true.
            return
         endif
         nint = numadd + 1
         dx = (xx(nmsh) - xx(nmsh-1))/nint
         xx(nmnew) = xx(nmsh)
         do 200 i = nmsh, nmnew-1 
            xx(i) = xx(nmsh-1) + (i-nmsh+1)*dx
 200            continue
      
      else

         if (numadd .gt. 9) numadd = 9

*  Here, intref lies between 2 and nmsh-2.  Add numadd points to
*  each of the three intervals intref-1, intref and intref+1.

         nmnew = nmsh + 3*numadd 
         if (nmnew .gt. nmax) then
            if (iprint .ge. 0)  write(6,903) nmnew
            maxmsh = .true.
            return
         endif

*  noalt is the number of points at the right end of the interval 
*  whose numerical values remain the same, but whose indices change.  
*  nochsm is the smallest index in the new ordering of one of these
*  points.

         noalt = nmsh - intref - 1
         nochsm = nmnew - noalt + 1

*  Renumber the noalt unchanged points at the right end of the interval
*  (in reverse order).

         j = 0
         do 300 i = nmnew, nochsm, -1
            xx(i) = xx(nmsh-j)
            j = j + 1
 300            continue

*  Add numadd points to the three contiguous intervals.
*  The remaining points at the left end of the interval retain 
*  their original indices, and are left unchanged.
         
         nint = numadd + 1
         innew = nochsm - nint
         do 320 i = intref+1, intref-1, -1
            xx(innew) = xx(i)
            dx = (xx(innew + nint) - xx(innew))/nint
            do 310 j = 1, numadd
               xx(innew + j) = xx(innew) + j*dx
 310                  continue
            innew = innew - nint
 320            continue

      endif
      nmsh = nmnew

      if(iprint .ge. 0)  write(6,904) nmsh
      return
 901  format(1h , ' smpmsh') 
 902  format(1h ,'nmsh, intref, numadd',3i6)
 903  format(1h , ' smpmsh.  maximum points exceeded, nmnew =',i6)
 904  format(1h ,'smpmsh, new nmsh =',i7)
      end

   

      subroutine unimsh(nmsh, aleft, aright, nfxpnt, fixpnt, xx)
      implicit double precision (a-h,o-z)
      integer  nmsh, nfxpnt
      dimension fixpnt(*), xx(nmsh)

      logical pdebug, use_c, comp_c
      common/algprs/ nminit, pdebug, iprint, idum, uval0, use_c, comp_c

      intrinsic max


*  Given a left endpoint aleft, a right endpoint aright,
*  a set of nfxpnt fixed points fixpnt(i), i = 1,...,nfxpnt,
*  (where fixpnt(i) is different from aleft and aright for all i),
*  and an initial target number nmsh of mesh points,
*  the subroutine unimsh generates a piecewise uniform mesh 
*  beginning at aleft, ending at aright, and with equally 
*  spaced points between aleft and fixpnt(1), then between 
*  fixpnt(1) and fixpnt(2), ..., and finally between 
*  fixpnt(nfxpnt) and aright.  The final number of intervals
*  is the maximum of nfxpnt+2 and the initial value of nmsh.

*  In the simplest case when nfxpnt = 0, unimsh generates a
*  uniform mesh with nmsh intervals in the closed interval
*  (aleft, aright).

*  On exit, the integer nmsh contains the number of mesh points
*  (which is the maximum of the initial nmsh and nfxpnt).
*  The array xx (of dimension nmsh) contains the mesh points.

      if (iprint .ge. 0) write(6,901) nmsh

      if (nfxpnt .eq. 0) then

*  If there are no interior fixed points, the spacing is uniform 
*  throughout the interval.  Calculate the spacing dx
*  and set up the xx array.

        ninter = nmsh - 1

         dx = (aright - aleft)/ninter
         do 10 i = 1, ninter
            xx(i) = aleft + (i-1)*dx
 10             continue
         xx(nmsh) = aright
         return
      endif

*  We know that there is at least one fixed point strictly between
*  the endpoints.

      if (nmsh .lt. nfxpnt+2)  nmsh = nfxpnt + 2
      ninter = nmsh - 1
      xx(1) = aleft
      ileft = 1
      xleft = aleft
      totint = aright - aleft
      ndif = ninter - nfxpnt
      do 50 j = 1, nfxpnt + 1

*  Deal in turn with the subintervals defined by the interval
*  boundaries and the fixed  points.

         if (j .lt. nfxpnt+1) then

*  The j-th fixed point is xright.  Calculate where it should
*  fall in the mesh.

            xright = fixpnt(j)
            nmin = ninter*(xright-aleft)/totint + 1.5d+0
            if (nmin .gt. ndif+j) nmin = ndif + j
            iright = max(ileft+1, nmin)
         else
            xright = aright
            iright = nmsh
         endif

*  npt is the number of equally spaced points that should
*  lie strictly between the (j-1)-th and j-th fixed points.

         xx(iright) = xright
         npt = iright - ileft - 1
         dx = (xright - xleft)/(npt + 1)
         do 30 i = 1, npt
            xx(ileft+i) = xleft + i*dx
 30             continue
         ileft = iright
         xleft = xright
 50       continue

      return
c
 901   format (1h ,'unimsh.  nmsh =',i5) 
      end 


      subroutine stats(len, elem, ebigst, esecnd, summod, index)
      implicit double precision (a-h, o-z)
      integer index
      dimension elem(len)

      intrinsic abs

      parameter  ( zero = 0.0d+0 )

*  Given the real array elem of length len, stats calculates
*  the following:
*      - summod, the sum of the magnitudes of the elements of elem;
*      - ebigst (the largest element in magnitude);
*      - index (the index in elem of ebigst); and 
*      - esecnd (the second largest element in magnitude, strictly
*          less than ebigst unless both are zero). 

      index = 1
      ebigst = zero
      esecnd = zero
      summod = zero

      do 100 i = 1, len
         elmod = abs(elem(i))
         summod = summod + elmod
         if (elmod .gt. ebigst) then
            esecnd = ebigst
            ebigst = elmod
            index = i
         elseif (elmod .gt. esecnd) then
            esecnd = elmod
         endif
 100      continue
      return
      end
      

      subroutine mshref(ncomp, nmsh, nlbc, ntol, ltol,
     *              iorder, rhs, tmwork,
     *              nmax, xx, nmold, xxold, ddouble, maxmsh,
     *              numbig, nummed,
     *              amg,stab_cond,stiff_cond,r4,
     *              nfxpnt, fixpnt,irefin,itcond,itcondmax)

      implicit double precision (a-h,o-z)

      dimension  ltol(ntol)
      dimension  rhs(ncomp*nmsh), tmwork(nmsh)
      dimension  xx(nmsh), xxold(nmold)
      dimension  fixpnt(*), irefin(*),r4(*)
      logical    ddouble, maxmsh
      dimension  amg(*)
      logical stab_cond, stiff_cond, forcedouble
      logical pdebug, use_c, comp_c
      common/algprs/ nminit, pdebug, iprint, idum, uval0, use_c, comp_c

*  blas: dcopy
      logical nodouble 
      parameter (two = 2.0d+0)
      parameter (bigfac = 10.0d+0, small = 1.0d-2, numpt = 14)


*  This routine performs calculations leading to a decision
*  about how the mesh will be refined, and then refines the mesh.

*  The choices for mesh refinement in this routine are either to 
*  double the mesh (i.e., divide each existing interval in half),
*  or to add points to just a few intervals.

*  The decision is made based on two criteria:  (1) the distribution 
*  of magnitudes of components of rhs (broadly speaking, if
*  the maximum component of rhs is much larger than the
*  average, points are added near the corresponding interval),
*  and (2) the history of previous mesh refinements (if points
*  have been added to only a few intervals already, this strategy
*  is abandoned for the moment and the mesh is doubled).

*  The decision is indicated by setting the logical flag ddouble
*  and (if ddouble is .false.) the integer intref.  
*  Setting ddouble to .true. means that the mesh should be doubled 
*  (i.e., the new mesh should contain twice as many intervals).
*  The integer intref indicates the region of the mesh where the 
*  points should be added (see the routine smpmsh), and numadd 
*  indicates how many points are to be added.
*  The integers nummed and numbig represent running totals,
*  used in deciding on the mesh refinement strategy.

*  If iorder = 4, meaning that we were just performing a Newton 
*  iteration for a 4th order solution, check all elements of rhs.  
*  If iorder .gt. 4, signalling that we were trying Newton
*  iterations for order 6 or 8, check only elements of rhs 
*  corresponding to components for which a tolerance is specified.

      if (pdebug) write(6,901) nummed, numbig
      if (use_c) then
c      nodouble = ((iorder.eq.4) .and. (stiff_cond) .and. (use_c))
       nodouble = ((iorder.eq.4) .and. 
     *       (stiff_cond .and. .not. stab_cond) .and. (use_c))
c       nodouble = (stiff_cond .and. .not. stab_cond) .and. (use_c)
c      nodouble = nodouble 
c    *  .or.((iorder.gt.4) .and. .not.stab_cond .and. (stiff_cond)
c    *   .and. (use_c))
      forcedouble = .false.
      endif
c 	write(6,*) 'mshref', itcond
      if (use_c) then
         if ( itcond .eq. itcondmax) then
            itcond = 0
            forcedouble = .true.
         else
            itcond = itcond + 1
            forcedouble = .false.
         endif
      endif 
      
      ninter = nmsh-1
      nup = ncomp
      if (iorder .gt. 4) nup = ntol

*  Check the vector rhs for a non-negligible component 
*  whose magnitude is significantly larger than the average.  
*  (small defines negligible, and bigfac defines significantly larger.)
            
      do 50 ic = 1, nup
         icmp = ic
         if (iorder .gt. 4) icmp = ltol(ic)

*  For component icmp, examine the ninter elements of rhs not 
*  corresponding to boundary conditions.
 
*  subroutine stats calculates rbigst and rsecnd (the first- and 
*  second-largest elements of rhs in magnitude), and the index 
*  intref of the interval in which the largest value occurs.
*  The value sumrhs is the sum of the magnitudes of the components
*  of rhs.

         indrhs = nlbc + ic
*  Copy the elements of rhs corresponding to interior mesh
*  points for component icmp into a single vector tmwork.

         call dcopy(ninter, rhs(indrhs), ncomp, tmwork, 1)

         call stats(ninter, tmwork, rbigst, rsecnd, 
     *                 sumrhs, intref)
         tstval = bigfac*(sumrhs-rbigst)/ninter
         if (pdebug) write(6,902) ic, tstval, rbigst, rsecnd
         if (rbigst .ge. small .and. rbigst .ge. tstval) go to 100
 50       continue

*  If we reach this point, no interval has a significantly larger 
*  element than average.  Set counters and double the mesh.

      numbig = 0
      nummed = 0
c      If(Iflnwt.eq.0) then
c        onto6=.true.
c      return
c      endif
      ddouble = .true.
      if (pdebug) write(6,903)
cf the mesh if not doubled if the problem is stiff and the order is 4

      if (.not.use_c.and..not.comp_c) nodouble=.false.
       if (nodouble .and. .not. forcedouble) then
          call selcondmsh(ncomp, nmsh, 
     *     nfxpnt, fixpnt,  nmax, xx,  irefin,
     *     nmold, xxold, ddouble, maxmsh,r4,amg)
           ddouble = .false.
        else
          call dblmsh (nmsh, nmax, xx, nmold, xxold, maxmsh)
          itcond = 0
        end if 
c      call dblmsh(nmsh, nmax, xx, nmold, xxold, maxmsh)
      return

*  To reach statement 100, it must be true that for some component 
*  icmp,  rbigst is non-negligible and large relative to the 
*  average.  intref indicates the region to which points may be added.
 
*  If too many specialized refinements (adding a few points to
*  a small number of intervals) have been made, signal that 
*  the mesh should be doubled.  
*  Otherwise, increment counters and add numadd points as indicated.

 100   continue
      if (rbigst .lt. two*rsecnd) nummed = nummed + 1
      numadd = numpt 
      numbig = numbig + 1
      ddouble = .false.
     
      if (rbigst .le. bigfac*rsecnd .or. numbig .gt. 8) then
         numbig = 0
         nummed = nummed + 1
         if (nummed .ge. 4 .and. iorder .eq. 4) then
            ddouble = .true.
            nummed = 0
         elseif (nummed .ge. 8 .and. iorder .gt. 4) then
            ddouble = .true.
            nummed = 0
         endif
      endif

*  Refine the mesh.
      if (pdebug) write(6,904) numbig, nummed
 
      if (.not.use_c.and..not.comp_c) nodouble=.false.
      if (ddouble) then
cf the mesh if not doubled if the problem is stiff and the order is 4
         if  (nodouble .and. .not. forcedouble)  then
           call selcondmsh(ncomp, nmsh, 
     *        nfxpnt, fixpnt,  nmax, xx,  irefin,
     *        nmold, xxold, ddouble, maxmsh,r4,amg)
           ddouble = .false.
         else
           call dblmsh (nmsh, nmax, xx, nmold, xxold, maxmsh)
           itcond = 0
         end if 
      else
cf if the problem is stiff and the order is 4 we use both the old technique
cf and the conditioning  
c         forcedouble = .false.
         if (nodouble .and. .not. forcedouble)  then
           numadd = numpt
           call smpselcondmsh(ncomp, nmsh, 
     *       nfxpnt, fixpnt,  nmax, xx,  irefin,intref,numadd,
     *       nmold, xxold, ddouble, maxmsh,r4,amg)
c            call selcondmsh(ncomp, nmsh, 
c     *        nfxpnt, fixpnt,  nmax, xx,  irefin,
c     *        nmold, xxold, ddouble, maxmsh,r4,amg)
c        elseif (forcedouble .and. use_c) then
c            ddouble = .true.
c            call dblmsh (nmsh, nmax, xx, nmold, xxold, maxmsh)
c            itcond = 0
        else
                numadd = numpt
                call smpmsh (nmsh, nmax, xx, intref, numadd,
     *          nmold, xxold, maxmsh)
                itcond =0
        endif
 
      endif
      if (ddouble .and. pdebug) write(6,905)
c      onto6=.false.
      return

 901  format(1h ,'mshref. nummed, numbig =',2i5)
 902  format(1h ,'ic, tst, bigst, second',i5, 3(1pe11.3))
 903  format(1h ,'No significantly large value')
 904  format(1h ,'numbig, nummed =',2i5)
 905  format(1h ,'double the mesh')
      end 

      subroutine errest (ncomp, nmsh, ntol, ltol, tol,
     *   nudim, u, uold, etest, errsum, errok)

      implicit double precision (a-h,o-z)
      dimension ltol(ntol), tol(ntol), u(nudim,nmsh), 
     *              uold(ncomp,nmsh), etest(ntol)
      logical errok

      logical pdebug, use_c, comp_c
      common/algprs/ nminit, pdebug, iprint, idum, uval0, use_c, comp_c
      intrinsic abs, max

      parameter ( one = 1.0d+0, zero = 0.0d+0 )

 
*  Given current and previous solutions u and uold on the same
*  mesh, errest calculates an error measure for each 
*  component for which a tolerance is specified.
*  The error measure is the usual relative error normalized 
*  by dividing by the tolerance.  On exit, errsum is the
*  sum of these error measures.

*  The array etest specifies the error test to be applied to each
*  error measure.

*  On exit, the logical flag errok
*   -- is .false. if any of the error measures exceeds the 
*      corresponding value in the array etest
*   -- is .true. if all error measures are less than the 
*      corresponding values of etest.

      if (pdebug) write(6,900)
      errsum = zero
      errok = .true.

      do 10 im = 1, nmsh
      do 10 it = 1, ntol
         icmp = ltol(it)
         er = u(icmp,im) - uold(icmp,im)
         denom = max(one, abs(uold(icmp,im)))
         errel = abs(er/(tol(it)*denom)) 
c         if (pdebug)  write(6,901) im, it, errel, etest(it)
         errsum = errsum + errel
c         if (errel .gt. etest(it)) errok = .false.
         If (Errel .gt. Etest(It)) Errok = .false.
 10       continue

      if (pdebug) write(6,902) errsum
      return
 900   format(1h ,'errest')
 901    format(1h ,2i5,2(1pe11.3))
 902     format(1h ,'errsum',1pe11.3)
      end 


      subroutine getptq( debug, mfsrch, nout, alfmax, alfsml, alfuzz, 
     *                   epsaf, epsag, eta, ftry, oldf, oldg,
     *                   rmu, tolabs, tolrel, toltny, imprvd,
     *                   inform, nfsrch, alfa, alfbst, fbest,
     *                   braktd, crampd, extrap,vset,wset,nsamea,nsameb,
     *                   a, b, fa, factor, fv, fw, xtry, xv, xw )

      implicit double precision (a-h,o-z)
      logical            debug, imprvd
      logical            braktd, crampd, extrap, vset, wset 
      integer            mfsrch, nout, inform, nfsrch, nsamea, nsameb 
c
c  *********************************************************************
c  getptq  is a step-length algorithm for minimizing a function of one
c  variable.  it will be called repeatedly by a search routine whose
c  purpose is to estimate a point  alfa = alfbst  that minimizes some 
c  function  f(alfa)  over the closed interval (0, alfmax). 
c
c  getptq  requires the function  f(alfa)  (but not its gradient)
c  to be evaluated at various points within the interval.  new
c  step-length estimates are computed using quadratic interpolation with
c  safeguards.
c
c  reverse communication is used to allow the calling program to
c  evaluate  f.  some of the parameters must be set or tested
c  by the calling program.  the remainder would ordinarily be local
c  variables.
c
c
c  input parameters (relevant to the calling program)
c  --------------------------------------------------
c
c  debug         specifies whether detailed output is wanted.
c
c  inform        must be nonzero on the first entry (e.g., -1).
c                it will be altered by  getptq  for later entries.
c
c  mfsrch        is an upper limit on the number of times  getptq  is 
c                to be entered consecutively with  inform = 0
c                (following an initial entry with  inform lt 0).
c
c  nout          is the file number to be used for printed output
c                if debug is true.
c
c  alfa          is the first estimate of the step length.  alfa  is
c                subsequently altered by  getptq  (see below).
c
c  alfmax        is the upper limit of the interval to be searched.
c
c  alfsml        is intended to prevent inefficiency when the optimum 
c                step is very small, for cases where the calling
c                program would prefer to re-define  f(alfa).  alfsml is
c                allowed to be zero. early termination will occur if
c                getptq  determines that the optimum step lies
c                somewhere in the interval  (0, alfsml)  (but not if
c                alfmax .le. alfsml).
c
c  epsaf         is an estimate of the absolute precision in the
c                computed values of  f. 
c
c  eta           controls the accuracy of the search.  it must lie
c                in the range   0.0  le  eta  lt  1.0.  decreasing
c                eta  tends to increase the accuracy of the search.
c
c  oldf          is the value of  f(0). 
c
c  oldg          is an estimate of the gradient of  f  at  alfa = 0.
c                it should be non-positive.
c
c  rmu           controls what is meant by a significant decrease in  f.
c                the final  f(alfbst)  should lie on or below the line
c                      l(alfa)  =  oldf + alfa*rmu*oldg.
c                rmu  should be in the open interval (0, 0.5).
c                the value  rmu = 1.0d-4  is good for most purposes.
c
c  tolabs,tolrel define a function  tol(alfa) = tolrel*alfa + tolabs
c                such that if  f  has already been evaluated at step
c                alfa,  then it will not be evaluated at any point
c                closer than  tol(alfa).
c                these values may be reduced by  getptq  if they seem 
c                to be too large.
c
c  toltny        is the smallest value that  tolabs  is allowed to be 
c                reduced to.
c
c
c  output parameters (relevant to the calling program)
c  ---------------------------------------------------
c
c  imprvd        is true if the previous step  alfa  was the best
c                point so far.  any related quantities (e.g., arrays) 
c                should be saved by the calling program before paying 
c                attention to  inform.
c
c  inform = 0    means the calling program should evaluate
c                           ftry = f(alfa)
c                for the new trial step  alfa,  and then re-enter
c                getptq.
c
c  inform = 1    means the search has terminated successfully
c                with a step  alfbst  that is less than the 
c                upper bound  alfmax.
c
c  inform = 2    means the search has terminated successfully
c                with a step  alfbst  that is equal to the
c                upper bound  alfmax.
c
c  inform = 3    means that the search failed to find a point of
c                sufficient decrease in  mfsrch  functions, but an
c                improved point was found.
c
c  inform = 4    means  alfmax  is so small that a search should
c                not have been done.
c
c  inform = 5    means that the search was terminated prematurely
c                because of the value of  alfsml  (see above).
c
c  inform = 6    means the search has failed to find a useful step.  if
c                the subroutine for the function and gradient has been
c                programmed correctly, this will usually occur if the 
c                minimum lies very close to  alfa = 0  or the gradient
c                is not sufficiently accurate.
c
c  inform = 7    means that the value of  g(0) was positive on entry. 
c
c  alfa          is the step at which the next function value must be 
c                computed.
c
c  alfbst        should be accepted by the calling program as the
c                required step-length estimate, whenever  getptq
c                returns  inform = 1,  2  or  3.
c
c  fbest         will be the corresponding value of  f.
c
c
c  the following parameters retain information between entries
c  -----------------------------------------------------------
c
c  alfuzz        is such that, if the final  alfa  lies in the interval
c                (0,alfuzz)  and  abs( f(alfa)-oldf ) le epsaf,  alfa 
c                cannot be guaranteed to be a point of sufficient
c                decrease.
c
c  braktd        is false if  f  has not been evaluated at the far end
c                of the interval of uncertainty.  in this case, the
c                point  b  will be at  alfmax + tol(alfmax).
c
c  crampd        is true if  alfmax  is very small (le tolabs).
c                if the search fails, this indicates that a zero
c                step should be taken.
c
c  extrap        is true if alfbst has moved at least once and  xv
c                lies outside the interval of uncertainty.  in this
c                case, extra safeguards are applied to allow for
c                instability in the polynomial fit.
c
c  vset          records whether a third-best point has been
c                determined.
c
c  wset          records whether a second-best point has been
c                determined.  it will always be true by the 
c                time the convergence test is applied (label 300).
c
c  nsamea        is the number of consecutive times that the left-hand
c                end of the interval of uncertainty has remained the
c                same.
c
c  nsameb        similarly for the right-hand end.
c
c  a, b, alfbst  define the current interval of uncertainty.
c                the required minimum lies somewhere within the
c                closed interval  (alfbst + a, alfbst + b). 
c
c  alfbst        is the best point so far.  it is strictly within the 
c                the interval of uncertainty except when it lies at the
c                left-hand end when  alfbst  has not been moved.
c                hence we have    a le 0,   b gt 0.
c
c  fbest         is the value of  f  at the point  alfbst.
c
c  fa            is the value of  f  at the point  alfbst + a.
c
c  factor        controls the rate at which extrapolated estimates of 
c                alfa  may expand into the interval of uncertainty.
c                factor is not used if the minimum has been bracketed 
c                (i.e., when the variable  braktd  is true).
c
c  fv, fw        are the values of  f  at the points  alfbst + xv,
c                alfbst + xw.  they are not defined until  vset
c                or  wset  (respectively) is true.
c
c  ftry          is the value of  f  at the new point  alfbst + xtry. 
c
c  xtry          is the trial point within the shifted interval (a, b).
c                the new trial function value must be computed at the 
c                point  alfa  =  alfbst + xtry.
c
c  xv            is such that  alfbst + xv  is the third-best point.
c                it is not defined until  vset  is true.
c
c  xw            is such that  alfbst + xw  is the second-best point. 
c                it is not defined until  wset  is true.
c                in some cases,  xw  will replace a previous  xw  that
c                has a lower function but has just been excluded from 
c                the interval of uncertainty.
c
c
c  systems optimization laboratory, stanford university, california.
c  original version february 1982.  rev. may 1983.
c  *********************************************************************
c
      logical            closef, conv1, conv2, conv3, convrg
      logical            moved, sigdec, xinxw
      data               zero, point1, half/ 0.0d+0,  0.1d+0, 0.5d+0/ 
      data                one,  two,  five   / 1.0d+0, 2.0d+0,  5.0d+0/
      data                ten, eleven      /10.0d+0, 11.0d+0        / 
c
c
c  local variables
c  ---------------
c
c  closef        is true if the worst function  fv  is within  epsaf
c                of  fbest  (up or down).
c
c  convrg        will be set to true if at least one of the convergence
c                conditions holds at  alfbst.
c
c  moved         is true if a better point has been found (alfbst gt 0).
c
c  sigdec        says whether  fbest  represents a significant decrease
c                in the function, compared to the initial value  oldf.
c
c  xinxw         is true if  xtry  is in  (xw,0)  or  (0,xw).
c  ---------------------------------------------------------------------
c
      imprvd = .false.
      if (inform .ne. -1) go to 100
c
c  ---------------------------------------------------------------------
c  first entry.  initialize various quantities, check input data and
c  prepare to evaluate the function at the initial step  alfa.
c  ---------------------------------------------------------------------
      nfsrch = 0
      alfbst = zero 
      fbest  = oldf 
      if (oldg   .gt.      zero) go to 970
      if (oldg   .ge. (- epsag)) go to 960
      if (alfmax .le.    toltny) go to 940
c
      braktd = .false.
      crampd = alfmax .le. tolabs
      extrap = .false.
      vset   = .false.
      wset   = .false.
      nsamea = 0
      nsameb = 0
      alfuzz = two*epsaf/(rmu*abs( oldg ))
      a      = zero 
      b      = alfmax + (tolrel*alfmax + tolabs)
      fa     = oldf 
      factor = five 
      tol    = tolabs
      xtry   = alfa 
      if (debug) write (nout, 1000) alfmax, oldf, oldg, tolabs,
     *   alfuzz, epsaf, epsag, tolrel, crampd
      go to 800
c
c  ---------------------------------------------------------------------
c  subsequent entries.
c  the function has just been evaluated at  alfa = alfbst + xtry,
c  giving  ftry.
c  ---------------------------------------------------------------------
 100   nsamea = nsamea + 1
      nsameb = nsameb + 1
      xtry   = alfa - alfbst
      moved  = alfbst .gt. zero
c
c  check if  xtry  is in the interval  (xw,0)  or  (0,xw).
c
      xinxw  = .false.
      if (wset) xinxw =       zero .lt. xtry  .and.  xtry .le. xw
     *                  .or.    xw .le. xtry  .and.  xtry .lt. zero
c
c  see if the new step is better.
c
      deltaf = ftry   - oldf
      ctry   = deltaf - alfa*rmu*oldg
      if (alfa .le. alfuzz) sigdec = deltaf .le. (- epsaf)
      if (alfa .gt. alfuzz) sigdec = ctry   .le.    epsaf
      imprvd = sigdec  .and.  ( ftry - fbest ) .le. (- epsaf)
c
      if (debug) write (nout, 1100) alfa, ftry, ctry
      if (.not. imprvd) go to 130
c
c  we seem to have an improvement.  the new point becomes the
c  origin and other points are shifted accordingly.
c
      if (.not. wset) go to 110
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
      extrap = .not. xinxw
c
c  decrease the length of the interval of uncertainty.
c
      if (xtry .lt. zero) go to 120
      a      = xw
      fa     = fw
      nsamea = 0
      go to 300
 120   b      = xw
      nsameb = 0
      braktd = .true.
      go to 300
c
c  the new function value is no better than the best point found so far.
c  the point  xtry  must be a new end point of the interval of
c  uncertainty.
c
 130   if (xtry .ge. zero) go to 140
      a      = xtry 
      fa     = ftry 
      nsamea = 0
      go to 150
 140   b      = xtry 
      nsameb = 0
      braktd = .true.
c
c  the origin remains unchanged but  xtry  may qualify as  xw.
c
 150   if (.not. wset)   go to 160
      if ((ftry - fw) .gt. epsaf) go to 170
      xv     = xw
      fv     = fw
      vset   = .true.
 160   xw     = xtry 
      fw     = ftry 
      wset   = .true.
      if (moved) extrap = xinxw
      go to 300
c
c  ftry  is no better than  fbest  or  fw.  if the best point has not 
c  been moved, there must be more than one minimum.
c
 170   if (moved) go to 175
      xw     = xtry 
      fw     = ftry 
      go to 300
c
c  ftry  is no better than  fbest  or  fw,  but  xtry  may become  xv.
c  extrap  has the value set in the previous entry.
c
 175   if (.not. vset) go to 180
      if ((ftry - fv) .gt. epsaf  .and.  extrap) go to 300
 180   if (xinxw) go to 190
      xv     = xtry 
      fv     = ftry 
      vset   = .true.
      go to 300
 190   if (vset) xw = xv
      if (vset) fw = fv
      xv     = xtry 
      fv     = ftry 
      vset   = .true.
c
c  ---------------------------------------------------------------------
c  check the termination criteria.
c  ---------------------------------------------------------------------
 300   tol    = tolrel*alfbst + tolabs
      deltaf = fbest  - oldf

      cbest  = deltaf - alfbst*rmu*oldg 
      if (alfbst .le. alfuzz) sigdec = deltaf .le. (- epsaf)
      if (alfbst .gt. alfuzz) sigdec = cbest  .le.    epsaf 
      closef = .false.
      if (vset) closef = abs( fbest - fv ) .le. epsaf
c
      conv1  = max( abs( a ), b )  .le.  (tol + tol)

*  conv2 changed by mhw, 20 sept 1992, to allow it to be
*  satified for any significant decrease in f
*      conv2  =  moved  .and.  sigdec
*     *                 .and.  abs( fa - fbest )  .le.  a*eta*oldg
      conv2 = moved .and. sigdec
      conv3  = closef  .and.  (sigdec  .or.
     *                        (.not. moved)  .and.  (b .le. alfuzz))
      convrg = conv1  .or.  conv2  .or.  conv3
c
      atrue  = alfbst + a
      btrue  = alfbst + b
      alfaw  = alfbst + xw
      gap    = b - a
      if (debug) write (nout, 1200) atrue, btrue, gap, tol, 
     *   nsamea, nsameb, braktd, closef, imprvd, conv1, conv2, conv3, 
     *   extrap, alfbst, fbest, cbest, alfaw, fw
      if (vset) alfav  = alfbst + xv
      if (debug  .and.  vset) write (nout, 1300) alfav, fv
      if (convrg  .and.  moved) go to 910
c
c  exit if the step is too small.
c
      if (btrue   .lt.  alfsml) go to 950
      
      if (nfsrch  .ge.  mfsrch) go to 930
      if (.not. convrg) go to 400
c
c  a better point has not yet been found (the step  xw  is no better
c  than step  zero).  check that the change in  f  is consistent with a
c  perturbation in  x  of  tol, the estimate of the minimum spacing
c  constant.  if the change in  f  is larger than  epsaf,  the value
c  of  tol  is reduced.
c
      tol    = xw/ten
      tolabs = tol
      if (abs(fw - oldf) .gt. epsaf  .and.  tol .gt. toltny) go to 400
      if (crampd) go to 940
      go to 960
c
c  ---------------------------------------------------------------------
c  proceed with the computation of a trial step length.
c  the choices are...
c  1. parabolic fit using function values only.
c  2. damped parabolic fit if the regular fit appears to be 
c     consistently over-estimating the distance to the minimum.
c  3. bisection, geometric bisection, or a step of  tol  if the
c     parabolic fit is unsatisfactory.
c  ---------------------------------------------------------------------
 400   xmidpt = half*(a + b)
      q      = zero 
      s      = zero 
c
c  ---------------------------------------------------------------------
c  fit a parabola.
c  ---------------------------------------------------------------------
c
c  check if there are two or three points for the parabolic fit.
c
      gw = (fw - fbest)/xw
      if (vset  .and.  moved) go to 450 
c
c  only two points available.  use  fbest,  fw  and the derivative
c  oldg.
c
      if (.not. moved) s = oldg
      if (      moved) s = oldg - two*gw
      q = two*(oldg - gw)
      if (debug) write (nout, 2100)
      go to 600
c
c  three points available.  use  fbest,  fw  and  fv.
c
 450   gv = (fv - fbest)/xv
      s  = gv - (xv/xw)*gw
      q  = two*(gv - gw)
      if (debug) write (nout, 2200)
c
c  ---------------------------------------------------------------------
c  construct an artificial interval  (artifa, artifb)  in which the
c  new estimate of the step length must lie.  set a default value of
c  xtry  that will be used if the polynomial fit is rejected.  in the 
c  following, the interval  (a,b)  is considered the sum of two
c  intervals of lengths  dtry  and  daux, with common end point at the
c  best point (zero).  dtry  is the length of the interval into which 
c  the default  xtry  will be placed and  endpnt  denotes its non-zero
c  end point.  the magnitude of  xtry  is computed so that the exponents
c  of  dtry  and  daux  are approximately bisected.
c  ---------------------------------------------------------------------
 600   artifa = a
      artifb = b
      if (braktd) go to 610
c
c  the minimum has not been bracketed.  set an artificial upper bound 
c  by expanding the interval  xw  by a suitable factor.
c
      xtry   = - factor*xw
      artifb =   xtry
      if (alfbst + xtry .lt. alfmax) factor = five*factor
      go to 700
c
c  the minimum has been bracketed.
c  if the gradient at the origin is being used for the
c  polynomial fit, the default  xtry  is one tenth of  xw.
c
 610   if (vset  .and.  moved) go to 620 
      xtry   = xw/ten
      if (debug) write (nout, 2400) xtry
      go to 700
c
c  three points exist in the interval of uncertainty.  check whether
c  the points are configured for an extrapolation or interpolation.
c
 620   if (extrap) go to 660
c
c  if the interpolation appears to be consistently over-estimating the
c  distance to the minimum,  damp the interpolation step.
c
      if (nsamea .lt. 3  .and.  nsameb .lt. 3) go to 630
      factor = factor / five
      s      = factor * s
      go to 640
 630   factor = one
c
c  the points are configured for an interpolation.  the artificial
c  interval will be just  (a,b).   set  endpnt  so that  xtry
c  lies in the larger of the intervals  (a,0)  and  (0,b).
c
 640    if (xmidpt .lt. zero) endpnt = a
      if (xmidpt .gt. zero) endpnt = b
c
c  if a bound has remained the same for three iterations, set  endpnt 
c  so that  xtry  is likely to replace the offending bound. 
c
      if (nsamea .ge. 3) endpnt = a
      if (nsameb .ge. 3) endpnt = b
      go to 680
c
c  the points are configured for an extrapolation.
c
 660   if (xw .lt. zero) endpnt = b
      if (xw .gt. zero) endpnt = a
c
c  compute the default value of  xtry.
c
 680   dtry = abs( endpnt )
      daux = gap - dtry
      if (daux .ge. dtry)   xtry = five*dtry*(point1 + dtry/daux)/eleven
      if (daux .lt. dtry)   xtry = half*sqrt( daux )*sqrt( dtry )
      if (endpnt .lt. zero) xtry = - xtry
      if (debug) write (nout, 2500) xtry, daux, dtry
c
c  if the points are configured for an extrapolation set the artificial
c  bounds so that the artificial interval lies strictly within  (a,b).
c  if the polynomial fit is rejected,  xtry  will remain at the relevant
c  artificial bound.
c
      if (extrap  .and.  xtry .le. zero) artifa = xtry
      if (extrap  .and.  xtry .gt. zero) artifb = xtry
c
c  ---------------------------------------------------------------------
c  the polynomial fits give  (s/q)*xw  as the new step.
c  reject this step if it lies outside  (artifa, artifb).
c  ---------------------------------------------------------------------
 700   if (q .eq. zero) go to 800
      if (q .lt. zero) s = - s
      if (q .lt. zero) q = - q
      if (s*xw .lt. q*artifa   .or.   s*xw .gt. q*artifb) go to 800
c
c  accept the polynomial fit. 
c
      xtry = zero
      if (abs( s*xw ) .ge. q*tol) xtry = (s/q)*xw 
      if (debug) write (nout, 2600) xtry
c
c  ---------------------------------------------------------------------
c  test for  xtry  being larger than  alfmax  or too close to  a  or  b.
c  ---------------------------------------------------------------------
 800   if (braktd) go to 810
c
c  if the step is close to or larger than  alfmax,  replace it by
c  alfmax  (to force evaluation of the function at the boundary).
c
      alfa   = alfbst + xtry
      if (alfmax - alfa .gt. (tolrel*alfmax + tolabs)) go to 810
      braktd = .true.
      xtry   = alfmax - alfbst
      alfa   = alfmax
      go to 900
c
c  otherwise, the function must not be evaluated at a point too close 
c  to  a  or  b.  (it has already been evaluated at both those points.)
c
 810   xmidpt = half*(a + b)
      if (xtry .gt. a + tol  .and.  xtry .lt. b - tol) go to 820
      if (xmidpt .gt. zero) xtry =   tol
      if (xmidpt .le. zero) xtry = - tol
c
c
c  f  must not be calculated too close to  alfbst.
c
 820   if (abs( xtry ) .lt. tol  .and.  xmidpt .lt. zero) xtry = - tol 
      if (abs( xtry ) .lt. tol  .and.  xmidpt .ge. zero) xtry =   tol 
      alfa   = alfbst + xtry
c
c  ---------------------------------------------------------------------
c  exit.
c  ---------------------------------------------------------------------
c
c  new function value required.
c
 900   inform = 0
      go to 990
c
c  convergence test satisfied.
c
 910   inform = 1
      if (alfa .eq. alfmax) inform = 2
      go to 990
c
c  mfsrch  function evaluations without sufficient decrease, but an
c  improved point was found.
c
 930   if (.not. moved) go to 960
      inform = 3
      go to 990
c
c  zero step (alfmax too small).
c
 940   inform = 4
      go to 990
c
c  premature termination.  the step is smaller than  alfsml.
c
 950   inform = 5
      go to 990
c
c  zero step (a sufficiently better point could not be found).
c
 960   inform = 6
      go to 990
c
c  zero step (positive gradient at the starting point).
c
 970   inform = 7
c
c  exit.
c
 990    if (debug) write (nout, 3000)
      return
c
 1000 format(/ 31h alfmax  oldf    oldg    tolabs, 1p2e22.14, 1p2e16.8
     *       / 31h alfuzz  epsaf   epsag   tolrel, 1p2e22.14, 1p2e16.8
     *       / 31h crampd                        ,  l6)
 1100 format(/ 31h alfa    ftry    ctry          , 1p2e22.14, 1pe16.8)
 1200 format(/ 31h a       b       b - a   tol   , 1p2e22.14, 1p2e16.8
     *       / 31h nsamea  nsameb  braktd  closef, 2i3, 2l6 
     *       / 31h imprvd  convrg  extrap        ,  l6, 3x, 3l1, l6
     *       / 31h alfbst  fbest   cbest         , 1p2e22.14, 1pe16.8 
     *       / 31h alfaw   fw                    , 1p2e22.14)
 1300 format(  31h alfav   fv                    , 1p2e22.14 /)

 2100 format(30h parabolic fit,    two points.)
 2200 format(30h parabolic fit,  three points.)
 2400 format(31h exponent reduced.  trial point, 1p1e22.14) 
 2500 format(31h geo. bisection. xtry,daux,dtry, 1p3e22.14) 
 2600 format(31h polynomial fit accepted.  xtry, 1p1e22.14) 
 3000 format(53h ---------------------------------------------------- /)
c
c  end of getptq
      end 


      subroutine interp(ncomp, nmsh, xx, nudim, u, 
     *                    nmold, xxold, uold) 

      implicit double precision (a-h, o-z)
      dimension xx(*), u(nudim,*), xxold(*), uold(ncomp,*) 
      logical pdebug, use_c, comp_c
      common/algprs/ nminit, pdebug, iprint, idum, uval0, use_c, comp_c

* blas: dcopy

      parameter (zero = 0.0d+0)

*  interp performs piecewise linear interpolation of the old
*  solution uold at the nmold old mesh points xxold onto the nmsh
*  new mesh points xx, producing an interpolated solution u.
*  Note that no assumption is made that the new mesh has
*  more points than the old, nor that the new and old mesh
*  points are related in a specific way (except that their first
*  and last points are identical).

      if (pdebug) write(6,900)

*  By construction, xx(1) = xxold(1).  Copy the first ncomp
*  components of uold into those of u.

      call dcopy(ncomp, uold(1,1), 1, u(1,1), 1)  

      i = 2 
      do 100 im = 2, nmsh-1

 50          continue
         if (i .gt. nmold) return

*  Check whether the im-th point in the new mesh lies strictly
*  to the right of, or to the left of (or exactly on) the
*  i-th point in the old mesh.


         if (xx(im) .gt. xxold(i)) then
            i = i + 1
            go to 50
         else
            xdif = xxold(i) - xx(im)
            if (xdif .eq. zero) then

*  xx(im) and xxold(i) are identical.
  
               call dcopy(ncomp, uold(1,i), 1, u(1,im), 1)
               i = i + 1
            else
               xint = xxold(i) - xxold(i-1)
               xrat = xdif/xint
               do 70 k = 1, ncomp
                  u(k,im) = uold(k,i) + xrat*(uold(k,i-1)-uold(k,i))
 70                         continue
            endif
         endif

 100      continue
      call dcopy(ncomp, uold(1,nmold), 1, u(1,nmsh), 1)
      return
 900   format(1h ,'interp')
      end 


      subroutine rerrvl( ncomp, nmsh, nudim, u, usvrex, ntol, ltol,
     *          rerr, remax, itlmx, adjrer )
      implicit double precision (a-h,o-z)
      dimension ltol(*)
      dimension u(nudim, *), usvrex(ncomp, *)
      dimension rerr(ncomp, *)
      logical adjrer

      intrinsic abs, max

      parameter ( one = 1.0d+0, zero = 0.0d+0 )
  
*  rerrvl is used in considering Richardson extrapolation.
*  The two solutions u and usvrex have a special relationship:
*  u corresponds to a doubled mesh, with twice as many 
*  intervals, where each interval is half the size of that for 
*  usvrex's mesh.   nmsh is the number of mesh points in the
*  mesh for u.

*  remax is the maximum relative error, and itlmx is the
*  index of the tolerance for which the maximum occurs.

*  The array rerr contains the absolute value of the difference
*  between u and usvrex at the common mesh points, but is defined 
*  only for components for which an error tolerance is specified.

*  The logical variable adjrer is true on entry if the values in
*  rerr are to be adjusted for later use in selective mesh refinement.


      itlmx = 1
      remax = zero
      nmold = 1 + (nmsh-1)/2
      do 100 it = 1, ntol
         icmp = ltol(it)
         imnew = 1
         do 50 im = 1, nmold
            rerr(icmp, im) = abs(usvrex(icmp,im) - u(icmp,imnew))
            denom = max(one, abs(usvrex(icmp,im)))
            rerel = rerr(icmp,im)/denom
            if (rerel .gt. remax) then
               remax = rerel
               itlmx = it
            endif
            imnew = imnew + 2
 50             continue
 100             continue

      if (adjrer) then

*  Adjust the rerr array if it may be used later in selective
*  mesh refinement.

         do 150 it = 1, ntol
            icmp = ltol(it)
            do 150 im = 1, nmold - 1
               rerr(icmp,im) = max(rerr(icmp,im), 
     *                              rerr(icmp, im+1))
 150               continue
      endif

      return
      end

      double precision function dasum(n,dx,incx)
c#
c#     takes the sum of the absolute values.
c#     jack dongarra, linpack, 3/11/78.
c#     modified 3/93 to return if incx .le. 0.
c#     modified 12/3/93, array(1) declarations changed to array(*)
c#
      double precision dx(*),dtemp
      integer i,incx,m,mp1,n,nincx
c#
      dasum = 0.0d0
      dtemp = 0.0d0
      if( n .le. 0 .or. incx .le. 0 )return
      if(incx.eq.1)go to 20
c#
c#        code for increment not equal to 1
c#
      nincx = n*incx
      do 10 i = 1,nincx,incx
        dtemp = dtemp + dabs(dx(i))
 10      continue
      dasum = dtemp
      return
c#
c#        code for increment equal to 1
c#
c#
c#        clean-up loop
c#
 20    m = mod(n,6)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dtemp = dtemp + dabs(dx(i))
 30      continue
      if( n .lt. 6 ) go to 60
 40    mp1 = m + 1
      do 50 i = mp1,n,6
        dtemp = dtemp + dabs(dx(i)) + dabs(dx(i + 1)) + dabs(dx(i + 2))
     *  + dabs(dx(i + 3)) + dabs(dx(i + 4)) + dabs(dx(i + 5))
 50      continue
 60       dasum = dtemp
      return
      end


      subroutine selconderrmsh(ncomp, nmsh, ntol, ltol, tol,
     *     nfxpnt, fixpnt, ipow, nmax,
     *     xx, nudim, u, ermeas, irefin, ihcomp,
     *     nmold, xxold, ermx, ddouble, maxmsh,r4,amg,stab_cond)

      implicit double precision (a-h,o-z)

      dimension  ltol(ntol), tol(ntol), fixpnt(*)
      dimension  xx(*), u(nudim, *), ermeas(ncomp,*)
      dimension  irefin(nmsh-1), ihcomp(nmsh-1)
      dimension  xxold(*), ermx(*), r4(*), amg(*)
      logical    ddouble, maxmsh, stab_cond

      intrinsic abs, max, int

*  blas: dcopy
*  double precision dlog

      logical pdebug, use_c, comp_c
      common/algprs/ nminit, pdebug, iprint, idum, uval0, use_c, comp_c

      parameter  ( zero = 0.0d+0, one = 1.0d+0, onep1 = 1.1d+0 ) 
      parameter  ( erdcid = 5.0d+0 )
      parameter  ( phitst = 0.1d+0 )
  
      logical first, add
      save    first, rlndec
      data    first / .true. /

*  The routine selmsh performs selective mesh refinement, depending
*  on the error measure ermeas.
   
      if (first) then
         first = .false.
         rlndec = dlog(erdcid)
      endif

      maxmsh = .false.
*      Nref = .true.

      if (pdebug) write(6,901) nmsh, ipow

      frcpow = one/ipow
      ddouble = .false.
      nmold = nmsh
      ninter = nmsh - 1
c      Iprec = Min(Iprec,1)

*  Copy the current mesh into the xxold array.


      call dcopy(nmold, xx, 1, xxold, 1)
      ithres = 0
      thres = one

*  On input, the array ermeas represents some error measure defined
*  over the components and mesh intervals (not mesh points).  
*  It is normalized in the following loop with respect to the
*  tolerance array and the current solution.
*  The value errmax gives the maximum normalized error.

      errmax = zero
      do 120 im = 1, ninter
         ermx(im) = zero
         do 110 it = 1, ntol
            jcomp = ltol(it)
            denom = tol(it)*max(one, abs(u(jcomp,im)))
            ems = ermeas(jcomp,im)
            ermeas(jcomp,im) = abs(ems)/denom
*            if (pdebug .and. ermeas(jcomp,im) .ge. thres)
*     *             write(6,902) im,jcomp,ems,ermeas(jcomp,im)
            err = ermeas(jcomp, im)
            if (err .ge. ermx(im)) then
                ermx(im) = err
                ihcomp(im) = jcomp
            endif
 110            continue
         errmax = max(ermx(im), errmax)
 120      continue
      if (pdebug) write(6,903) errmax
c      write(6,903) errmax
      call moncondmsh(nmsh,xx,r1,r2,r3,fatt_r1r3,fatt_r3,nptcond,r4,amg)

      if (.not. stab_cond .and.  errmax .ge. 1.0e20
     *           .and. r1 .ge. 1.0d0 ) then
cf  only the conditioning
            call  selcondmsh(ncomp, nmsh, 
     *         nfxpnt, fixpnt,  nmax, xx,  irefin,
     *         nmold, xxold, ddouble, maxmsh,r4,amg)
      else
cf   the conditioning and the error
c      call moncondmsh(nmsh,xx,r1,r2,r3,fatt_r1r3,fatt_r3,nptcond,r4,amg)

      
       if (errmax .gt. zero .and. errmax .le. erdcid) then
c      if (errmax .gt. zero .and. errmax .le. 500d0) then
*  If errmax > 0 and .le. erdcid, find the smallest integer exponent ii 
*  such that (erdcid**ii)*errmax > erdcid.

         if(errmax .gt. one) then
            ii = 1
            decii = erdcid
         else
            ilg = -dlog(errmax)/rlndec
            ii = 2 + ilg
            decii = erdcid**ii
         endif
         write(6,*) 'decii',decii, errmax
*  Multiply error measures by erdcid**ii.

         errmax = decii*errmax
         do 140 im = 1, ninter
            ermx(im) = decii*ermx(im)
            do 140 it = 1, ntol
               jcomp = ltol(it)
               ermeas(jcomp,im) = decii*ermeas(jcomp,im)
 140               continue
      endif

 200   continue

*  For each interval im,  the integer irefin(im) is calculated 
*  based on an equidistrbution procedure involving the
*  threshold value thres.  If irefin(im) > 1, we add 
*  points to interval im.  If irefin(im) = 1, we check whether
*  point im can be removed from the mesh.

*  nmest is a lower bound on the number of points in the new mesh.
*  We do not know in advance how many points will be removed,
*  so nmest is computed by assuming all eligible points are removed.

      nmest = nmsh
      do 220 im = 1, ninter
         if (ermx(im).ge.thres) then
            irefin(im) = int((ermx(im))**frcpow) + 1
            nmest = nmest + irefin(im) - 1
         else
            irefin(im) = 1
            nmest = nmest - 1
         endif
 220      continue
c        write(6,*) 'irefin',( irefin(imj),imj=1,ninter )

       if (nptcond .ge. 4  .and. .not. stab_cond) then
c           nptcond = nptcond/2
       add = .false.
       do 221 im = 1, ninter-1
         if  (max(r4(im), r4(im+1)) .ge. fatt_r1r3) then
            if (.not. add) then
              irefin(im) = max(nptcond, irefin(im)) 
c              nmest = nmest + nptcond - 1
            endif 
              irefin(im+1) = max(nptcond, irefin(im+1))
c	        nmest = nmest + nptcond - 1
              add = .true.
         else
            irefin(im) = max(1, irefin(im))
            irefin(im+1) = max(1, irefin(im+1))
c            nmest = nmest - 1
            add = .false. 
         endif
 221     continue

c        do 221 im = 1, ninter
c         if ( r4(im) .gt. fatt_r1r3) then
c              irefin(im) =max( nptcond, irefin(im)) 
cc              nmest = nmest + nptcond - 1
c         else
c              irefin(im) =max( 1, irefin(im)) 
cc              nmest = nmest - 1
c         endif
c        
c 221  continue

          endif
        
*      if (pdebug) write(6,904) nmest, (irefin(i), i=1,ninter)

      if (nmest .gt. nmax) then

         go to 360

      elseif (nmest-1 .gt. 3*ninter ) then

         call dblmsh (nmsh, nmax, xx, nmold, xxold, maxmsh) 
         ddouble = .true.
         return
      endif

*  It appears that we can perform the desired selective mesh
*  refinement.

*  Now begin running through the mesh, adding and possibly deleting
*  points as indicated by the irefin array.

*  The integer new is a count of the number of intervals in
*  the tentative mesh being generated by the refinement strategy.

      new = 1

*  The first interval is treated as a special case, since xx(1)
*  always remains in the mesh, and cannot be a fixed point.

      rlen = xxold(2) - xx(1)
      slen = rlen
      if (irefin(1).gt.1) then
         dx = rlen/irefin(1)
         do 230 j = 2, irefin(1) 
            new = new + 1
            xx(new) = xx(1) + (j-1)*dx
 230            continue
      endif

*  The fixed points specified by the fixpnt array cannot be
*  removed from the mesh.  The value fxnext indicates the 'next'
*  fixed point.  When no further fixed points remain to be processed
*  (or if nfxpnt = 0), fxnext is set to a value strictly larger than
*  the last mesh point, so that no mesh point can equal fxnext.  
*  This way we need to compare only the values of xxold(i)
*  and fxnext.

      ifxcnt = 1
      if (nfxpnt .eq. 0) then
         fxnext = onep1*abs(xxold(nmsh))
      else
         fxnext = fixpnt(ifxcnt)
      endif

*  jtkout is a counter of the number of consecutive points that 
*  have been removed from the mesh.

      jtkout = 0
      do 330 im = 2, ninter
         rlold = rlen
         rlen = xxold(im+1) - xxold(im)

*  If xxold(im) is the next fixed point, it cannot be removed
*  and so we don't test its error estimates.

         if(xxold(im) .eq. fxnext)  then

            ifxcnt = ifxcnt + 1
            if(ifxcnt .gt. nfxpnt) then
               fxnext = onep1*abs(xxold(nmsh))
            else
               fxnext = fixpnt(ifxcnt)
            endif

         elseif (irefin(im).eq.1) then

*  If xxold(im) is not a fixed point and irefin(im) = 1, we may wish 
*  to remove point im from the mesh.

*  If we are considering removing points and jtkout = 0, this
*  is the first point in a possible consecutive set to be removed,
*  and we initialize phihat, which represents a maximum of
*  certain estimates.
*  If jtkout is not zero, previous points contiguous to this
*  point have been removed, and phihat retains its previous value.

            slen = slen + rlen

            if (jtkout .eq. 0) then
                ind1 = ihcomp(im-1)
                phihat = ermeas(ind1,im-1)/(rlold**ipow)
            endif
            phihat = max(phihat,
     *                 ermeas(ihcomp(im),im)/(rlen**ipow))
            val1 = phihat*(slen**ipow)
            if (val1 .le. phitst 
     *             .and. jtkout .lt. 4 ) then
c     *             .and. r4(im) .lt. fatt_r1r3) then
*  Increment the counter of removed points.
*  'Remove' the mesh point xxold(im) by not including it.

               jtkout = jtkout+1
               go to 330
            endif
*        end of logic for irefin(im) = 1.
         endif

         jtkout = 0 
         new = new + 1
         xx(new) = xxold(im)
         if (irefin(im) .gt. 1) then
            dx = rlen/irefin(im)
            do 300 j = 2, irefin(im)
              new = new + 1
              xx(new) = xxold(im) + (j-1)*dx
 300                 continue
         endif
         slen = rlen
         
         if (new .gt. nmax) then

*  If the new mesh contains too many points, branch out of the
*  loop to try alternative strategies.

            go to 360

         elseif (new .gt. 3*ninter)  then

*  Here, the new mesh does not exceed the specified maximum,
*  but has more than 3 times as many intervals as the old mesh.
*  Try doubling the mesh if possible.

            call dcopy(nmsh, xxold, 1, xx, 1)
            call dblmsh (nmsh, nmax, xx, nmold, xxold, maxmsh) 
            ddouble = .true.
            return
   
          end if
 330       continue

*  To end up here, we have processed the entire interval,
*  and have neither exceeded the specified maximum nor
*  exceeded three times the number of intervals in the old
*  mesh.  The last mesh point remains unchanged.
     
      new = new + 1
      xx(new) = xxold(nmsh)
      nmsh = new
      maxmsh = .false.
      if (iprint .ge. 0) write(6,905) nmsh
      return

 360   continue
*  To reach here, the number of mesh points created at some stage
*  of the refinement process was larger than the maximum permitted 
*  value nmax.  

*  Check whether the mesh can safely be doubled.
  
      if ((2*nmsh-1) .lt. nmax) then

*  Double the mesh.
         call  dblmsh (nmsh, nmax, xx, nmold, xxold, maxmsh) 
         ddouble = .true.

*  If the number of intervals is too large and the mesh cannot be 
*  doubled, increase the threshold thres by a factor of erdcid and
*  try the selective refinement again.  
*  If this happens three times without success or if thres exceeds
*  or is equal to errmax, stop.  (In this case, we know already 
*  that doubling the mesh produces too many points.)

      elseif (thres .lt. errmax .and. ithres .lt. 3) then
         ithres = ithres + 1
         thres = erdcid*thres
         if(thres .gt. errmax) thres = errmax
         call dcopy(nmsh, xxold, 1, xx, 1)
         go to 200
      else
         nmsh = 2*nmsh - 1
         maxmsh = .true.
      endif
cf only the conditioning
      end if
      return

 901  format(1h ,'selconderrmsh.  nmsh, ipow =',2i5)
 902  format(1h ,'im, jcomp, ermeas, normalized er',2i5,2(1pe11.3))
 903  format(1h ,'errmax',1pe11.3)
 904  format(1h ,'nmest, irefin',(10i5))
 905  format(1h ,'selconderrmsh.  new nmsh =',i8)
 910  format(1h ,'ihcomp',(10i5))
      end 



   
  
      subroutine selcondmsh(ncomp, nmsh, 
     *     nfxpnt, fixpnt,  nmax, xx,  irefin,
     *     nmold, xxold, ddouble, maxmsh, r4, amg)

      implicit double precision (a-h,o-z)

      dimension  fixpnt(*)
      dimension  xx(*)
      dimension  irefin(*)
      dimension  xxold(*),  amg(*), r4(nmsh)
      logical    ddouble, maxmsh, add

      intrinsic abs, max, int

*  blas: dcopy
*  double precision dlog

      logical pdebug, use_c, comp_c
      common/algprs/ nminit, pdebug, iprint, idum, uval0, use_c, comp_c

      parameter  ( zero = 0.0d+0, one = 1.0d+0, onep1 = 1.1d+0 ) 
      parameter  ( erdcid = 5.0d+0 )
  
  
   

*  The routine selcondmsh performs selective mesh refinement, depending
*  on the monitor function based on the conditioning parameters.
   
   
      maxmsh = .false.
 

      if (pdebug) write(6,901) nmsh

      ddouble = .false.
      nmold = nmsh
      ninter = nmsh - 1

*  Copy the current mesh into the xxold array.


      call dcopy(nmold, xx, 1, xxold, 1)

      ithres = 0
      thres = one

*  On input, the array amg represents the conditioning vector defined
*  over the components and mesh intervals (not mesh points).  
*  We compute the monitor function and the related parameters in the
*  followinf function

      call moncondmsh(nmsh,xx,r1,r2,r3,fatt_r1r3,fatt_r3,nptcond,r4,amg)
 


*  For each interval im,  the integer irefin(im) is calculated 
*  based on an equidistrbution procedure involving the
*  threshold value thres.  If irefin(im) > 1, we add 
*  points to interval im.  If irefin(im) = 1, we check whether
*  point im can be removed from the mesh.

*  nmest is a lower bound on the number of points in the new mesh.
*  We do not know in advance how many points will be removed,
*  so nmest is computed by assuming all eligible points are removed.

      add = .false.
      nmest = nmsh
c      nptcond = nptcond/2
      do 220 im = 1, ninter-1
         if  (max(r4(im), r4(im+1)) .ge. fatt_r1r3) then
            if (.not. add) then
              irefin(im) = nptcond 
              nmest = nmest + nptcond - 1
            endif 
              irefin(im+1) = nptcond
	        nmest = nmest + nptcond - 1
              add = .true.
         else
            irefin(im) = 1
            irefin(im+1) = 1
            nmest = nmest - 1
            add = .false. 
         endif
        
  220 continue


      if (nmest .gt. nmax) then

         go to 360

      endif

*  It appears that we can perform the desired selective mesh
*  refinement.

*  Now begin running through the mesh, adding and possibly deleting
*  points as indicated by the irefin array.

*  The integer new is a count of the number of intervals in
*  the tentative mesh being generated by the refinement strategy.

      new = 1

*  The first interval is treated as a special case, since xx(1)
*  always remains in the mesh, and cannot be a fixed point.

      rlen = xxold(2) - xx(1)
      slen = rlen
      if (irefin(1).gt.1) then
         dx = rlen/irefin(1)
         do 230 j = 2, irefin(1) 
            new = new + 1
            xx(new) = xx(1) + (j-1)*dx
  230    continue
      endif

*  The fixed points specified by the fixpnt array cannot be
*  removed from the mesh.  The value fxnext indicates the 'next'
*  fixed point.  When no further fixed points remain to be processed
*  (or if nfxpnt = 0), fxnext is set to a value strictly larger than
*  the last mesh point, so that no mesh point can equal fxnext.  
*  This way we need to compare only the values of xxold(i)
*  and fxnext.

      ifxcnt = 1
      if (nfxpnt .eq. 0) then
         fxnext = onep1*abs(xxold(nmsh))
      else
         fxnext = fixpnt(ifxcnt)
      endif

*  jtkout is a counter of the number of consecutive points that 
*  have been removed from the mesh.

      jtkout = 0
      do 330 im = 2, ninter
         
         rlold = rlen
         rlen = xxold(im+1) - xxold(im)

*  If xxold(im) is the next fixed point, it cannot be removed
*  and so we don't test its error estimates.

         if(xxold(im) .eq. fxnext)  then

            ifxcnt = ifxcnt + 1
            if(ifxcnt .gt. nfxpnt) then
               fxnext = onep1*abs(xxold(nmsh))
            else
               fxnext = fixpnt(ifxcnt)
            endif

         elseif (irefin(im).eq.1) then

*  If xxold(im) is not a fixed point and irefin(im) = 1, we may wish 
*  to remove point im from the mesh.

*  If we are considering removing points and jtkout = 0, this
*  is the first point in a possible consecutive set to be removed,
*  and we initialize phihat, which represents a maximum of
*  certain estimates.
*  If jtkout is not zero, previous points contiguous to this
*  point have been removed, and phihat retains its previous value.

            slen = slen + rlen
        
            if ( jtkout .lt. 1
     *             .and. r4(im) .le. fatt_r3 ) then

*  Increment the counter of removed points.
*  'Remove' the mesh point xxold(im) by not including it.

               jtkout = jtkout+1
               go to 330
            endif
*        end of logic for irefin(im) = 1.
         endif

         jtkout = 0 
         new = new + 1
         xx(new) = xxold(im)
         if (irefin(im) .gt. 1) then
            dx = rlen/irefin(im)
            do 300 j = 2, irefin(im)
              new = new + 1
              xx(new) = xxold(im) + (j-1)*dx
  300       continue
         endif
         slen = rlen
                  
         if (new .gt. nmax ) then

*  If the new mesh contains too many points, branch out of the
*  loop to try alternative strategies.

c            go to 360
c
c         elseif (new .gt. 3*ninter)  then

*  Here, the new mesh does not exceed the specified maximum,
*  but has more than 3 times as many intervals as the old mesh.
*  Try doubling the mesh if possible.
            nmsh = nmold
            call dcopy(nmsh, xxold, 1, xx, 1)
            call dblmsh (nmsh, nmax, xx, nmold, xxold, maxmsh) 
            ddouble = .true.
            return

 
         endif

 

  330 continue

*  To end up here, we have processed the entire interval,
*  and have neither exceeded the specified maximum nor
*  exceeded three times the number of intervals in the old
*  mesh.  The last mesh point remains unchanged.
     
      new = new + 1
      xx(new) = xxold(nmsh)
      nmsh = new
      maxmsh = .false.
      if (iprint .ge. 0) write(6,905) nmsh
      return

      
  360 continue
*  To reach here, the number of mesh points created at some stage
*  of the refinement process was larger than the maximum permitted 
*  value nmax.  

      nmsh = nmold
      call dcopy(nmsh, xxold, 1, xx, 1)
      maxmsh = .true.
   

      return



  901 format(1h ,'selcondmsh.  nmsh, ipow =',2i5)
  902 format(1h ,'im, jcomp, ermeas, normalized er',2i5,2(1pe11.3))
  903 format(1h ,'errmax',1pe11.3)
  904 format(1h ,'nmest, irefin',(10i5))
  905 format(1h ,'selcondmsh.  new nmsh =',i8)
  910 format(1h ,'ihcomp',(10i5))
      end 



      subroutine smpselcondmsh(ncomp, nmsh, 
     *     nfxpnt, fixpnt,  nmax, xx,  irefin, intref, numadd, 
     *     nmold, xxold, ddouble, maxmsh, r4, amg)

      implicit double precision (a-h,o-z)

      dimension  fixpnt(*)
      dimension  xx(*)
      dimension  irefin(*)
      dimension  xxold(*),  amg(*), r4(*)
      logical    ddouble, maxmsh, add

      intrinsic abs, max, int

*  blas: dcopy
*  double precision dlog

      logical pdebug, use_c, comp_c
      common/algprs/ nminit, pdebug, iprint, idum, uval0, use_c, comp_c

      parameter  ( zero = 0.0d+0, one = 1.0d+0,onep1 = 1.1d+0) 
      parameter  ( erdcid = 5.0d+0 )
  
  
  

*  The routine smpselcondmsh performs selective mesh refinement, by adding 
*  points to one or three interval(s) in the region indicated
*  by the integer intref (numadd gives the trial number of points to 
*  added in each  interval)  and by adding and removing point 
*  using  the monitor function based on the conditioning parameters.
  
  
  
      maxmsh = .false.
 

      if (pdebug) write(6,901) nmsh

      ddouble = .false.
      nmold = nmsh
      ninter = nmsh - 1

*  Copy the current mesh into the xxold array.


         call dcopy(nmold, xx, 1, xxold, 1)

      ithres = 0
      thres = one

*  On input, the array amg represents the conditioning vector defined
*  over the components and mesh intervals (not mesh points).  

      call moncondmsh(nmsh,xx,r1,r2,r3,fatt_r1r3,fatt_r3,nptcond,r4,amg)

           

*  For each interval im,  the integer irefin(im) is calculated 
*  based on an equidistrbution procedure involving the
*  threshold value thres.  If irefin(im) > 1, we add 
*  points to interval im.  If irefin(im) = 1, we check whether
*  point im can be removed from the mesh.

*  nmest is a lower bound on the number of points in the new mesh.
*  We do not know in advance how many points will be removed,
*  so nmest is computed by assuming all eligible points are removed.

      add = .false.
      nmest = nmsh
c      nptcond = nptcond/2
      do 220 im = 1, ninter-1
         if  (max(r4(im), r4(im+1)) .gt. fatt_r1r3) then
            if (.not. add) then
              irefin(im) = nptcond 
              nmest = nmest + nptcond - 1
            endif 
              irefin(im+1) = nptcond
	        nmest = nmest + nptcond - 1
              add = .true.
         else
            irefin(im) = 1
            irefin(im+1) = 1
            nmest = nmest - 1
            add = .false. 
         endif
        
  220 continue

c      do 220 im = 1, ninter
c         if ( r4(im) .gt. fatt_r1r3) then
c              irefin(im) = nptcond 
c              nmest = nmest + nptcond - 1
c         else
c              irefin(im) = 1 
c              nmest = nmest - 1
c         endif
c        
c  220 continue

      if(numadd .gt. 49) then
         numadd = 49
      elseif (numadd .lt. 4) then
         numadd = 4
      endif

        if ( intref .eq. 1) then
            irefin(1) = max(numadd+1,irefin(1))
            nmest = nmest + numadd - 1
        elseif (intref .eq. ninter) then
            irefin(ninter) = max(numadd+1,irefin(ninter))
            nmest = nmest + numadd - 1 
        else
	       if(numadd .gt. 9) then
             numadd = 9
           elseif (numadd .lt. 4) then
              numadd = 4
            endif
            irefin(intref-1) = max(numadd+1 , irefin(intref-1))
            irefin(intref)   = max(numadd+1 , irefin(intref))
            irefin(intref+1) = max(numadd+1 , irefin(intref+1))
            nmest = nmest + 3*numadd - 1 
        end if



      if (nmest .gt. nmax) then

         go to 360

      endif

*  It appears that we can perform the desired selective mesh
*  refinement.

*  Now begin running through the mesh, adding and possibly deleting
*  points as indicated by the irefin array.

*  The integer new is a count of the number of intervals in
*  the tentative mesh being generated by the refinement strategy.

      new = 1

*  The first interval is treated as a special case, since xx(1)
*  always remains in the mesh, and cannot be a fixed point.

      rlen = xxold(2) - xx(1)
      slen = rlen
      if (irefin(1).gt.1) then
         dx = rlen/irefin(1)
         do 230 j = 2, irefin(1) 
            new = new + 1
            xx(new) = xx(1) + (j-1)*dx
  230    continue
      endif

*  The fixed points specified by the fixpnt array cannot be
*  removed from the mesh.  The value fxnext indicates the 'next'
*  fixed point.  When no further fixed points remain to be processed
*  (or if nfxpnt = 0), fxnext is set to a value strictly larger than
*  the last mesh point, so that no mesh point can equal fxnext.  
*  This way we need to compare only the values of xxold(i)
*  and fxnext.

      ifxcnt = 1
      if (nfxpnt .eq. 0) then
         fxnext = onep1*abs(xxold(nmsh))
      else
         fxnext = fixpnt(ifxcnt)
      endif

*  jtkout is a counter of the number of consecutive points that 
*  have been removed from the mesh.

      jtkout = 0
      do 330 im = 2, ninter
         
         rlold = rlen
         rlen = xxold(im+1) - xxold(im)

*  If xxold(im) is the next fixed point, it cannot be removed
*  and so we don't test its error estimates.

         if(xxold(im) .eq. fxnext)  then

            ifxcnt = ifxcnt + 1
            if(ifxcnt .gt. nfxpnt) then
               fxnext = onep1*abs(xxold(nmsh))
            else
               fxnext = fixpnt(ifxcnt)
            endif

         elseif (irefin(im).eq.1) then

*  If xxold(im) is not a fixed point and irefin(im) = 1, we may wish 
*  to remove point im from the mesh.

*  If we are considering removing points and jtkout = 0, this
*  is the first point in a possible consecutive set to be removed,
*  and we initialize phihat, which represents a maximum of
*  certain estimates.
*  If jtkout is not zero, previous points contiguous to this
*  point have been removed, and phihat retains its previous value.

            slen = slen + rlen
        
            if ( jtkout .lt. 1
     *          .and. (r4(im) .le. 5d-1*fatt_r3 
     *          .and. r1 .ge. 1.0d0)) then

*  Increment the counter of removed points.
*  'Remove' the mesh point xxold(im) by not including it.

               jtkout = jtkout+1
               go to 330
            endif
*        end of logic for irefin(im) = 1
         endif

         jtkout = 0 
         new = new + 1
         xx(new) = xxold(im)
         if (irefin(im) .gt. 1) then
            dx = rlen/irefin(im)
            do 300 j = 2, irefin(im)
              new = new + 1
              xx(new) = xxold(im) + (j-1)*dx
  300       continue
         endif
         slen = rlen
                  
         if (new .gt. nmax) then

*  If the new mesh contains too many points, branch out of the
*  loop to try alternative strategies.

            go to 360

         elseif (new .gt. 3*ninter)  then

*  Here, the new mesh does not exceed the specified maximum,
*  but has more than 3 times as many intervals as the old mesh.
*  Try doubling the mesh if possible.
            if (iprint .eq. 1) write(6,*) 'smpselcondmsh'
            nmsh = nmold
            call dcopy(nmsh, xxold, 1, xx, 1)
            call dblmsh (nmsh, nmax, xx, nmold, xxold, maxmsh) 
            ddouble = .true.
            return

         endif

 

  330 continue

*  To end up here, we have processed the entire interval,
*  and have neither exceeded the specified maximum nor
*  exceeded three times the number of intervals in the old
*  mesh.  The last mesh point remains unchanged.
     
      new = new + 1
      xx(new) = xxold(nmsh)
      nmsh = new
      maxmsh = .false.
      if (iprint .ge. 0) write(6,905) nmsh
      return

  360 continue
*  To reach here, the number of mesh points created at some stage
*  of the refinement process was larger than the maximum permitted 
*  value nmax.  

      nmsh = nmold
      call dcopy(nmsh, xxold, 1, xx, 1)
      maxmsh = .true.
   

      return



  901 format(1h ,'smpselcondmsh.  nmsh, ipow =',2i5)
  902 format(1h ,'im, jcomp, ermeas, normalized er',2i5,2(1pe11.3))
  903 format(1h ,'errmax',1pe11.3)
  904 format(1h ,'nmest, irefin',(10i5))
  905 format(1h ,'smpselcondmsh.  new nmsh =',i8)
  910 format(1h ,'ihcomp',(10i5))
      end 





      subroutine moncondmsh(nmsh,xx,r1,r2,r3,fatt_r1r3,fatt_r3,
     *                        nptcond,r4,amg)

      implicit double precision (a-h,o-z)

      dimension  xx(*)
      dimension  amg(*), r4(*)
 
      intrinsic abs, max, int

*  blas: dcopy
*  double precision dlog

      logical pdebug, use_c, comp_c
      common/algprs/ nminit, pdebug, iprint, idum, uval0, use_c, comp_c

      parameter  ( zero = 0.0d+0, one = 1.0d+0 ) 
      
* the function moncond compute the monitor function based on the 
* conditioning parameters and the factor used to perform the mesh selection
    
       do i=1,nmsh-1
         r4(i) = (xx(i+1)-xx(i))*dabs(amg(i+1)- amg(i))
 
       end do
     
       r2 = r4(1)
       do i=2,nmsh-1
          r2 = r2 + r4(i)
       end do
  
       do i=1,nmsh-1
          r4(i) = r4(i)+(xx(i+1)-xx(i))*(r2/(xx(nmsh)-xx(1)))*1e-5
       end do
       r1 = r4(1)
       do i=2,nmsh-1
          r1 = max(r1,r4(i))
       end do 
       
      
       do i=1,nmsh-1
         r4(i) = r4(i)/r1
       end do
      
     
        r2m = r4(1)
        r1m = r4(1)
       do i=2,nmsh-1
          r1m = min(r1m,r4(i))
          r2m = r2m + r4(i)
       end do 
     
       r3 = r2m/(nmsh-1)
       fatt_r3  = r3*1.0d-4
       fatt_r1r3= max(r3,0.65d0)
       if (r1 .gt. 1.0d0) then
         r1m = one
         nptm = 0
	    nptr = 0 
        do i=1,nmsh-1
	    if (r4(i) .ge. fatt_r1r3)  nptm = nptm + 1
	    if (r4(i) .le. fatt_r3)  nptr = nptr + 1
	    enddo
	
	   if (nptm .le. 1) then
	   nptcond =  14
        elseif (nptm .le. 2) then
	   nptcond =  10
	    elseif (nptm .le. 4) then
	   nptcond =  8
    	elseif (nptm .le. 8) then
	    nptcond = 6
	    elseif (nptm .le. int(nmsh/20) ) then
	     nptcond = 4
	    else
        nptcond = 2
	    endif   
       else
         nptcond=2
       endif 
 
       if (iprint .eq. 1) write(6,901)r1,r3,fatt_r1r3,nptcond,nptm,nptr

       

  901 format(1h ,'moncondmsh.', (1pe11.3), 2(1pe11.3), 3i10)
 
      end 


C==================================================================
C
        DOUBLE PRECISION FUNCTION ABDNRM(NBLOKS,NTOP,NBOT,NOVRLP,
     *                            NRWBLK,NCLBLK,TOP,A,BOT)

C******************************************************************
C
C       ABDNRM IS USED IN CONJUNCTION WITH DONEST TO COMPUTE THE
C       CONDITION NUMBER OF AN ALMOST BLOCK DIAGONAL MATRIX LIKE
C       THE ONES HANDLED BY COLROW. [SEE COMMENTS IN COLROW, CRDCMP]

C               *****  AUXILIARY PROGRAMS  *****

C       THE BLAS ARE REQUIRED BY ABDNRM. 
C       SPECIFICALLY, DASUM IS USED.
C             

C******************************************************************

        INTEGER NBLOKS,NTOP,NBOT,NOVRLP,NRWBLK,NCLBLK
        DOUBLE PRECISION TOP(NTOP,*),A(NRWBLK,NCLBLK,*),BOT(NBOT,*)
        INTEGER J,K
        DOUBLE PRECISION DASUM,TEMP
        TEMP = 0.0D0
C
C       FIRST, GO OVER THE COLUMNS OF TOP AND THE FIRST BLOCK:
C
        DO 10 J = 1,NOVRLP
           TEMP = MAX(TEMP, DASUM(NTOP,TOP(1,J),1) +
     *                      DASUM(NRWBLK,A(1,J,1),1))
 10                   CONTINUE
        DO 40 K = 1,NBLOKS-1
C
C          IN EACH BLOCK:
C
C          FIRST, THE COLUMNS FROM THE KTH BLOCK ALONE:
C
           DO 20 J = NOVRLP+1,NRWBLK
              TEMP = MAX(TEMP, DASUM(NRWBLK,A(1,J,K),1))
 20                            CONTINUE
C
C          NOW, TH COLUMNS WHICH INTERSECT BOTH THE KTH AND
C          (K+1)ST BLOCKS.
C
           DO 30 J = NRWBLK+1,NCLBLK
              TEMP = MAX(TEMP, DASUM(NRWBLK,A(1,J,K),1) +
     *                         DASUM(NRWBLK,A(1,J-NRWBLK,K+1),1))
 30                            CONTINUE
 40                                                CONTINUE
C
C       FINALLY, THE COLUMNS OF THE LAST BLOCK WHICH DO NOT OVERLAP
C       WITH ANYTHING.
C
        DO 50 J = NOVRLP+1,NRWBLK
           TEMP = MAX(TEMP, DASUM(NRWBLK,A(1,J,NBLOKS),1))
 50                   CONTINUE
C
C       AND THOSE COLUMNS OVERLAPPING WITH BOTH THE LAST BLOCK AND THE
C       BOTTOM BLOCK.
C
        DO 60 J = NRWBLK+1,NCLBLK
           TEMP = MAX(TEMP, DASUM(NRWBLK,A(1,J,NBLOKS),1) +
     *                      DASUM(NBOT,BOT(1,J-NRWBLK),1))
 60                   CONTINUE
        ABDNRM = TEMP
        RETURN
        END

      SUBROUTINE DONEST (N, V, X, ISGN, EST, KASE)
      INTEGER N, ISGN(N), KASE
      DOUBLE PRECISION V(N), X(N), EST

C
C     DONEST ESTIMATES THE 1-NORM OF A SQUARE, DOUBLE PRECISION MATRIX  A.
C     REVERSE COMMUNICATION IS USED FOR EVALUATING
C     MATRIX-VECTOR PRODUCTS. 
C
C     ON ENTRY
C
C        N       INTEGER
C                THE ORDER OF THE MATRIX.  N .GE. 1.
C
C        ISGN    INTEGER(N)
C                USED AS WORKSPACE.
C
C        KASE    INTEGER
C                = 0.
C
C     ON INTERMEDIATE RETURNS 
C
C        KASE    = 1 OR 2.
C
C        X       DOUBLE PRECISION(N)
C                MUST BE OVERWRITTEN BY 
C
C                     A*X,             IF KASE=1, 
C                     TRANSPOSE(A)*X,  IF KASE=2, 
C
C                AND DONEST MUST BE RE-CALLED, WITH ALL THE OTHER
C                PARAMETERS UNCHANGED.
C
C     ON FINAL RETURN
C
C        KASE    = 0.
C
C        EST     DOUBLE PRECISION
C                CONTAINS AN ESTIMATE (A LOWER BOUND) FOR NORM(A).
C
C        V       DOUBLE PRECISION(N)
C                = A*W,   WHERE  EST = NORM(V)/NORM(W)
C                         (W  IS NOT RETURNED).
C
C     THIS VERSION DATED MARCH 16, 1988.
C     NICK HIGHAM, UNIVERSITY OF MANCHESTER.
C
C     MODIFIED FOR DOUBLE PRECISION ON JUNE 11, 1996.
C
C     REFERENCE
C     N.J. HIGHAM (1987) FORTRAN CODES FOR ESTIMATING
C     THE ONE-NORM OF A REAL OR COMPLEX MATRIX, WITH APPLICATIONS
C     TO CONDITION  ESTIMATION, NUMERICAL ANALYSIS REPORT NO. 135,
C     UNIVERSITY OF MANCHESTER, MANCHESTER M13 9PL, ENGLAND.
C
C     SUBROUTINES AND FUNCTIONS
C     BLAS     IDAMAX, DASUM, DCOPY
C     GENERIC  ABS, NINT, DFLOAT, SIGN
C
        INTRINSIC DFLOAT
        DOUBLE PRECISION DFLOAT

        INTRINSIC ABS
        DOUBLE PRECISION ABS

        INTRINSIC SIGN
        DOUBLE PRECISION SIGN


      DOUBLE PRECISION DASUM
      INTEGER IDAMAX

      INTEGER ITMAX
      PARAMETER (ITMAX = 5)
      DOUBLE PRECISION ZERO,ONE,TWO
      PARAMETER (ZERO = 0.0D0)
      PARAMETER (ONE = 1.0D0)
      PARAMETER (TWO = 2.0D0)
C
C     INTERNAL VARIABLES
      INTEGER I, ITER, J, JLAST, JUMP
      DOUBLE PRECISION ALTSGN, ESTOLD, TEMP
C
      SAVE
C
      IF (KASE .EQ. 0) THEN
         DO 10,I = 1,N
            X(I) = ONE/DFLOAT(N)
 10                      CONTINUE
         KASE = 1
         JUMP = 1
         RETURN
      ENDIF
C
      GOTO (100, 200, 300, 400, 500) JUMP
C
C     ................ ENTRY   (JUMP = 1)
C     FIRST ITERATION.  X HAS BEEN OVERWRITTEN BY A*X.
C
 100     CONTINUE
      IF (N .EQ. 1) THEN
         V(1) = X(1)
         EST = ABS(V(1))
C        ... QUIT
         GOTO 510
      ENDIF
      EST = DASUM(N,X,1)
C
      DO 110,I = 1,N
         X(I) = SIGN(ONE,X(I))
         ISGN(I) = NINT(X(I)) 
 110           CONTINUE
      KASE = 2
      JUMP = 2
      RETURN
C
C     ................ ENTRY   (JUMP = 2)
C     FIRST ITERATION.  X HAS BEEN OVERWRITTEN BY TRANSPOSE(A)*X.
C
 200     CONTINUE
      J = IDAMAX(N,X,1)
      ITER = 2
C
C     MAIN LOOP - ITERATIONS 2,3,...,ITMAX.
C
 220     CONTINUE
      DO 230,I = 1,N
         X(I) = ZERO 
 230           CONTINUE
      X(J) = ONE
      KASE = 1
      JUMP = 3
      RETURN
C
C     ................ ENTRY   (JUMP = 3)
C     X HAS BEEN OVERWRITTEN BY A*X.
C
 300     CONTINUE
      CALL DCOPY(N,X,1,V,1)
      ESTOLD = EST
      EST = DASUM(N,V,1)
      DO 310,I = 1,N
         IF ( NINT( SIGN(ONE,X(I)) ) .NE. ISGN(I) ) GOTO 320
 310           CONTINUE
C     REPEATED SIGN VECTOR DETECTED, HENCE ALGORITHM HAS CONVERGED.
      GOTO 410
C
 320     CONTINUE
C     TEST FOR CYCLING.
      IF (EST .LE. ESTOLD) GOTO 410
C      
      DO 330,I = 1,N
         X(I) = SIGN(ONE,X(I))
         ISGN(I) = NINT(X(I)) 
 330           CONTINUE
      KASE = 2
      JUMP = 4
      RETURN
C
C     ................ ENTRY   (JUMP = 4)
C     X HAS BEEN OVERWRITTEN BY TRANSPOSE(A)*X.
C
 400     CONTINUE
      JLAST = J
      J = IDAMAX(N,X,1)
      IF (   (  X(JLAST) .NE. ABS(X(J))  ) .AND.
     +       (ITER .LT. ITMAX)   ) THEN
         ITER = ITER + 1
         GOTO 220
      ENDIF
C
C     ITERATION COMPLETE.  FINAL STAGE. 
C
 410     CONTINUE
      ALTSGN = ONE
      DO 420,I = 1,N
         X(I) = ALTSGN * (ONE + DFLOAT(I-1)/DFLOAT(N-1))
         ALTSGN = -ALTSGN
 420           CONTINUE
      KASE = 1
      JUMP = 5
      RETURN
C
C     ................ ENTRY   (JUMP = 5)
C     X HAS BEEN OVERWRITTEN BY A*X.
C
 500     CONTINUE
      TEMP = TWO*DASUM(N,X,1)/DFLOAT(3*N) 
      IF (TEMP. GT. EST) THEN 
         CALL DCOPY(N,X,1,V,1)
         EST = TEMP 
      ENDIF
C
 510     KASE = 0
      RETURN
C
      END

c This is the code that has gone into netlib as a replacement.

      DOUBLE PRECISION FUNCTION D1MACH(I)
      INTEGER I
C
C  DOUBLE-PRECISION MACHINE CONSTANTS
C
C  D1MACH( 1) = B**(EMIN-1), THE SMALLEST POSITIVE MAGNITUDE.
C
C  D1MACH( 2) = B**EMAX*(1 - B**(-T)), THE LARGEST MAGNITUDE.
C
C  D1MACH( 3) = B**(-T), THE SMALLEST RELATIVE SPACING.
C
C  D1MACH( 4) = B**(1-T), THE LARGEST RELATIVE SPACING.
C
C  D1MACH( 5) = LOG10(B)
C
C  TO ALTER THIS FUNCTION FOR A PARTICULAR ENVIRONMENT,
C  THE DESIRED SET OF DATA STATEMENTS SHOULD BE ACTIVATED BY
C  REMOVING THE C FROM COLUMN 1.
C  ON RARE MACHINES A STATIC STATEMENT MAY NEED TO BE ADDED.
C  (BUT PROBABLY MORE SYSTEMS PROHIBIT IT THAN REQUIRE IT.)
C
C  FOR IEEE-ARITHMETIC MACHINES (BINARY STANDARD), ONE OF THE FIRST
C  TWO SETS OF CONSTANTS BELOW SHOULD BE APPROPRIATE.  IF YOU DO NOT
C  KNOW WHICH SET TO USE, TRY BOTH AND SEE WHICH GIVES PLAUSIBLE
C  VALUES.
C
C  WHERE POSSIBLE, DECIMAL, OCTAL OR HEXADECIMAL CONSTANTS ARE USED
C  TO SPECIFY THE CONSTANTS EXACTLY.  SOMETIMES THIS REQUIRES USING
C  EQUIVALENT INTEGER ARRAYS.  IF YOUR COMPILER USES HALF-WORD
C  INTEGERS BY DEFAULT (SOMETIMES CALLED INTEGER*2), YOU MAY NEED TO
C  CHANGE INTEGER TO INTEGER*4 OR OTHERWISE INSTRUCT YOUR COMPILER
C  TO USE FULL-WORD INTEGERS IN THE NEXT 5 DECLARATIONS.
C
C  COMMENTS JUST BEFORE THE END STATEMENT (LINES STARTING WITH *)
C  GIVE C SOURCE FOR D1MACH.
C
      INTEGER SMALL(2)
      INTEGER LARGE(2)
      INTEGER RIGHT(2)
      INTEGER DIVER(2)
      INTEGER LOG10(2)
      INTEGER SC, CRAY1(38), J
      COMMON /D9MACH/ CRAY1
C/6S
C/7S
      SAVE SMALL, LARGE, RIGHT, DIVER, LOG10, SC
C/
      DOUBLE PRECISION DMACH(5)
C
      EQUIVALENCE (DMACH(1),SMALL(1))
      EQUIVALENCE (DMACH(2),LARGE(1))
      EQUIVALENCE (DMACH(3),RIGHT(1))
      EQUIVALENCE (DMACH(4),DIVER(1))
      EQUIVALENCE (DMACH(5),LOG10(1))
C
C     MACHINE CONSTANTS FOR BIG-ENDIAN IEEE ARITHMETIC (BINARY FORMAT)
C     MACHINES IN WHICH THE MOST SIGNIFICANT BYTE IS STORED FIRST,
C     SUCH AS THE AT&T 3B SERIES, MOTOROLA 68000 BASED MACHINES (E.G.
C     SUN 3), AND MACHINES THAT USE SPARC, HP, OR IBM RISC CHIPS.
C
C      DATA SMALL(1),SMALL(2) /    1048576,          0 /
C      DATA LARGE(1),LARGE(2) / 2146435071,         -1 /
C      DATA RIGHT(1),RIGHT(2) / 1017118720,          0 /
C      DATA DIVER(1),DIVER(2) / 1018167296,          0 /
C      DATA LOG10(1),LOG10(2) / 1070810131, 1352628735 /, SC/987/
C
C     MACHINE CONSTANTS FOR LITTLE-ENDIAN (BINARY) IEEE ARITHMETIC
C     MACHINES IN WHICH THE LEAST SIGNIFICANT BYTE IS STORED FIRST,
C     E.G. IBM PCS AND OTHER MACHINES THAT USE INTEL 80X87 OR DEC
C     ALPHA CHIPS.
C
      DATA SMALL(1),SMALL(2) /          0,    1048576 /
      DATA LARGE(1),LARGE(2) /         -1, 2146435071 /
      DATA RIGHT(1),RIGHT(2) /          0, 1017118720 /
      DATA DIVER(1),DIVER(2) /          0, 1018167296 /
      DATA LOG10(1),LOG10(2) / 1352628735, 1070810131 /, SC/987/
C
C     MACHINE CONSTANTS FOR AMDAHL MACHINES.
C
C      DATA SMALL(1),SMALL(2) /    1048576,          0 /
C      DATA LARGE(1),LARGE(2) / 2147483647,         -1 /
C      DATA RIGHT(1),RIGHT(2) /  856686592,          0 /
C      DATA DIVER(1),DIVER(2) /  873463808,          0 /
C      DATA LOG10(1),LOG10(2) / 1091781651, 1352628735 /, SC/987/
C
C     MACHINE CONSTANTS FOR THE BURROUGHS 1700 SYSTEM.
C
C      DATA SMALL(1) / ZC00800000 /
C      DATA SMALL(2) / Z000000000 /
C
C      DATA LARGE(1) / ZDFFFFFFFF /
C      DATA LARGE(2) / ZFFFFFFFFF /
C
C      DATA RIGHT(1) / ZCC5800000 /
C      DATA RIGHT(2) / Z000000000 /
C
C      DATA DIVER(1) / ZCC6800000 /
C      DATA DIVER(2) / Z000000000 /
C
C      DATA LOG10(1) / ZD00E730E7 /
C      DATA LOG10(2) / ZC77800DC0 /, SC/987/
C
C     MACHINE CONSTANTS FOR THE BURROUGHS 5700 SYSTEM.
C
C      DATA SMALL(1) / O1771000000000000 /
C      DATA SMALL(2) / O0000000000000000 /
C
C      DATA LARGE(1) / O0777777777777777 /
C      DATA LARGE(2) / O0007777777777777 /
C
C      DATA RIGHT(1) / O1461000000000000 /
C      DATA RIGHT(2) / O0000000000000000 /
C
C      DATA DIVER(1) / O1451000000000000 /
C      DATA DIVER(2) / O0000000000000000 /
C
C      DATA LOG10(1) / O1157163034761674 /
C      DATA LOG10(2) / O0006677466732724 /, SC/987/
C
C     MACHINE CONSTANTS FOR THE BURROUGHS 6700/7700 SYSTEMS.
C
C      DATA SMALL(1) / O1771000000000000 /
C      DATA SMALL(2) / O7770000000000000 /
C
C      DATA LARGE(1) / O0777777777777777 /
C      DATA LARGE(2) / O7777777777777777 /
C
C      DATA RIGHT(1) / O1461000000000000 /
C      DATA RIGHT(2) / O0000000000000000 /
C
C      DATA DIVER(1) / O1451000000000000 /
C      DATA DIVER(2) / O0000000000000000 /
C
C      DATA LOG10(1) / O1157163034761674 /
C      DATA LOG10(2) / O0006677466732724 /, SC/987/
C
C     MACHINE CONSTANTS FOR FTN4 ON THE CDC 6000/7000 SERIES.
C
C      DATA SMALL(1) / 00564000000000000000B /
C      DATA SMALL(2) / 00000000000000000000B /
C
C      DATA LARGE(1) / 37757777777777777777B /
C      DATA LARGE(2) / 37157777777777777774B /
C
C      DATA RIGHT(1) / 15624000000000000000B /
C      DATA RIGHT(2) / 00000000000000000000B /
C
C      DATA DIVER(1) / 15634000000000000000B /
C      DATA DIVER(2) / 00000000000000000000B /
C
C      DATA LOG10(1) / 17164642023241175717B /
C      DATA LOG10(2) / 16367571421742254654B /, SC/987/
C
C     MACHINE CONSTANTS FOR FTN5 ON THE CDC 6000/7000 SERIES.
C
C      DATA SMALL(1) / O"00564000000000000000" /
C      DATA SMALL(2) / O"00000000000000000000" /
C
C      DATA LARGE(1) / O"37757777777777777777" /
C      DATA LARGE(2) / O"37157777777777777774" /
C
C      DATA RIGHT(1) / O"15624000000000000000" /
C      DATA RIGHT(2) / O"00000000000000000000" /
C
C      DATA DIVER(1) / O"15634000000000000000" /
C      DATA DIVER(2) / O"00000000000000000000" /
C
C      DATA LOG10(1) / O"17164642023241175717" /
C      DATA LOG10(2) / O"16367571421742254654" /, SC/987/
C
C     MACHINE CONSTANTS FOR CONVEX C-1
C
C      DATA SMALL(1),SMALL(2) / '00100000'X, '00000000'X /
C      DATA LARGE(1),LARGE(2) / '7FFFFFFF'X, 'FFFFFFFF'X /
C      DATA RIGHT(1),RIGHT(2) / '3CC00000'X, '00000000'X /
C      DATA DIVER(1),DIVER(2) / '3CD00000'X, '00000000'X /
C      DATA LOG10(1),LOG10(2) / '3FF34413'X, '509F79FF'X /, SC/987/
C
C     MACHINE CONSTANTS FOR THE CRAY 1, XMP, 2, AND 3.
C
C      DATA SMALL(1) / 201354000000000000000B /
C      DATA SMALL(2) / 000000000000000000000B /
C
C      DATA LARGE(1) / 577767777777777777777B /
C      DATA LARGE(2) / 000007777777777777776B /
C
C      DATA RIGHT(1) / 376434000000000000000B /
C      DATA RIGHT(2) / 000000000000000000000B /
C
C      DATA DIVER(1) / 376444000000000000000B /
C      DATA DIVER(2) / 000000000000000000000B /
C
C      DATA LOG10(1) / 377774642023241175717B /
C      DATA LOG10(2) / 000007571421742254654B /, SC/987/
C
C     MACHINE CONSTANTS FOR THE DATA GENERAL ECLIPSE S/200
C
C     SMALL, LARGE, RIGHT, DIVER, LOG10 SHOULD BE DECLARED
C     INTEGER SMALL(4), LARGE(4), RIGHT(4), DIVER(4), LOG10(4)
C
C     NOTE - IT MAY BE APPROPRIATE TO INCLUDE THE FOLLOWING LINE -
C     STATIC DMACH(5)
C
C      DATA SMALL/20K,3*0/,LARGE/77777K,3*177777K/
C      DATA RIGHT/31420K,3*0/,DIVER/32020K,3*0/
C      DATA LOG10/40423K,42023K,50237K,74776K/, SC/987/
C
C     MACHINE CONSTANTS FOR THE HARRIS SLASH 6 AND SLASH 7
C
C      DATA SMALL(1),SMALL(2) / '20000000, '00000201 /
C      DATA LARGE(1),LARGE(2) / '37777777, '37777577 /
C      DATA RIGHT(1),RIGHT(2) / '20000000, '00000333 /
C      DATA DIVER(1),DIVER(2) / '20000000, '00000334 /
C      DATA LOG10(1),LOG10(2) / '23210115, '10237777 /, SC/987/
C
C     MACHINE CONSTANTS FOR THE HONEYWELL DPS 8/70 SERIES.
C
C      DATA SMALL(1),SMALL(2) / O402400000000, O000000000000 /
C      DATA LARGE(1),LARGE(2) / O376777777777, O777777777777 /
C      DATA RIGHT(1),RIGHT(2) / O604400000000, O000000000000 /
C      DATA DIVER(1),DIVER(2) / O606400000000, O000000000000 /
C      DATA LOG10(1),LOG10(2) / O776464202324, O117571775714 /, SC/987/
C
C     MACHINE CONSTANTS FOR THE IBM 360/370 SERIES,
C     THE XEROX SIGMA 5/7/9 AND THE SEL SYSTEMS 85/86.
C
C      DATA SMALL(1),SMALL(2) / Z00100000, Z00000000 /
C      DATA LARGE(1),LARGE(2) / Z7FFFFFFF, ZFFFFFFFF /
C      DATA RIGHT(1),RIGHT(2) / Z33100000, Z00000000 /
C      DATA DIVER(1),DIVER(2) / Z34100000, Z00000000 /
C      DATA LOG10(1),LOG10(2) / Z41134413, Z509F79FF /, SC/987/
C
C     MACHINE CONSTANTS FOR THE INTERDATA 8/32
C     WITH THE UNIX SYSTEM FORTRAN 77 COMPILER.
C
C     FOR THE INTERDATA FORTRAN VII COMPILER REPLACE
C     THE Z'S SPECIFYING HEX CONSTANTS WITH Y'S.
C
C      DATA SMALL(1),SMALL(2) / Z'00100000', Z'00000000' /
C      DATA LARGE(1),LARGE(2) / Z'7EFFFFFF', Z'FFFFFFFF' /
C      DATA RIGHT(1),RIGHT(2) / Z'33100000', Z'00000000' /
C      DATA DIVER(1),DIVER(2) / Z'34100000', Z'00000000' /
C      DATA LOG10(1),LOG10(2) / Z'41134413', Z'509F79FF' /, SC/987/
C
C     MACHINE CONSTANTS FOR THE PDP-10 (KA PROCESSOR).
C
C      DATA SMALL(1),SMALL(2) / "033400000000, "000000000000 /
C      DATA LARGE(1),LARGE(2) / "377777777777, "344777777777 /
C      DATA RIGHT(1),RIGHT(2) / "113400000000, "000000000000 /
C      DATA DIVER(1),DIVER(2) / "114400000000, "000000000000 /
C      DATA LOG10(1),LOG10(2) / "177464202324, "144117571776 /, SC/987/
C
C     MACHINE CONSTANTS FOR THE PDP-10 (KI PROCESSOR).
C
C      DATA SMALL(1),SMALL(2) / "000400000000, "000000000000 /
C      DATA LARGE(1),LARGE(2) / "377777777777, "377777777777 /
C      DATA RIGHT(1),RIGHT(2) / "103400000000, "000000000000 /
C      DATA DIVER(1),DIVER(2) / "104400000000, "000000000000 /
C      DATA LOG10(1),LOG10(2) / "177464202324, "047674776746 /, SC/987/
C
C     MACHINE CONSTANTS FOR PDP-11 FORTRANS SUPPORTING
C     32-BIT INTEGERS (EXPRESSED IN INTEGER AND OCTAL).
C
C      DATA SMALL(1),SMALL(2) /    8388608,           0 /
C      DATA LARGE(1),LARGE(2) / 2147483647,          -1 /
C      DATA RIGHT(1),RIGHT(2) /  612368384,           0 /
C      DATA DIVER(1),DIVER(2) /  620756992,           0 /
C      DATA LOG10(1),LOG10(2) / 1067065498, -2063872008 /, SC/987/
C
C      DATA SMALL(1),SMALL(2) / O00040000000, O00000000000 /
C      DATA LARGE(1),LARGE(2) / O17777777777, O37777777777 /
C      DATA RIGHT(1),RIGHT(2) / O04440000000, O00000000000 /
C      DATA DIVER(1),DIVER(2) / O04500000000, O00000000000 /
C      DATA LOG10(1),LOG10(2) / O07746420232, O20476747770 /, SC/987/
C
C     MACHINE CONSTANTS FOR PDP-11 FORTRANS SUPPORTING
C     16-BIT INTEGERS (EXPRESSED IN INTEGER AND OCTAL).
C
C     SMALL, LARGE, RIGHT, DIVER, LOG10 SHOULD BE DECLARED
C     INTEGER SMALL(4), LARGE(4), RIGHT(4), DIVER(4), LOG10(4)
C
C      DATA SMALL(1),SMALL(2) /    128,      0 /
C      DATA SMALL(3),SMALL(4) /      0,      0 /
C
C      DATA LARGE(1),LARGE(2) /  32767,     -1 /
C      DATA LARGE(3),LARGE(4) /     -1,     -1 /
C
C      DATA RIGHT(1),RIGHT(2) /   9344,      0 /
C      DATA RIGHT(3),RIGHT(4) /      0,      0 /
C
C      DATA DIVER(1),DIVER(2) /   9472,      0 /
C      DATA DIVER(3),DIVER(4) /      0,      0 /
C
C      DATA LOG10(1),LOG10(2) /  16282,   8346 /
C      DATA LOG10(3),LOG10(4) / -31493, -12296 /, SC/987/
C
C      DATA SMALL(1),SMALL(2) / O000200, O000000 /
C      DATA SMALL(3),SMALL(4) / O000000, O000000 /
C
C      DATA LARGE(1),LARGE(2) / O077777, O177777 /
C      DATA LARGE(3),LARGE(4) / O177777, O177777 /
C
C      DATA RIGHT(1),RIGHT(2) / O022200, O000000 /
C      DATA RIGHT(3),RIGHT(4) / O000000, O000000 /
C
C      DATA DIVER(1),DIVER(2) / O022400, O000000 /
C      DATA DIVER(3),DIVER(4) / O000000, O000000 /
C
C      DATA LOG10(1),LOG10(2) / O037632, O020232 /
C      DATA LOG10(3),LOG10(4) / O102373, O147770 /, SC/987/
C
C     MACHINE CONSTANTS FOR THE PRIME 50 SERIES SYSTEMS
C     WITH 32-BIT INTEGERS AND 64V MODE INSTRUCTIONS,
C     SUPPLIED BY IGOR BRAY.
C
C      DATA SMALL(1),SMALL(2) / :10000000000, :00000100001 /
C      DATA LARGE(1),LARGE(2) / :17777777777, :37777677775 /
C      DATA RIGHT(1),RIGHT(2) / :10000000000, :00000000122 /
C      DATA DIVER(1),DIVER(2) / :10000000000, :00000000123 /
C      DATA LOG10(1),LOG10(2) / :11504046501, :07674600177 /, SC/987/
C
C     MACHINE CONSTANTS FOR THE SEQUENT BALANCE 8000
C
C      DATA SMALL(1),SMALL(2) / $00000000,  $00100000 /
C      DATA LARGE(1),LARGE(2) / $FFFFFFFF,  $7FEFFFFF /
C      DATA RIGHT(1),RIGHT(2) / $00000000,  $3CA00000 /
C      DATA DIVER(1),DIVER(2) / $00000000,  $3CB00000 /
C      DATA LOG10(1),LOG10(2) / $509F79FF,  $3FD34413 /, SC/987/
C
C     MACHINE CONSTANTS FOR THE UNIVAC 1100 SERIES.
C
C      DATA SMALL(1),SMALL(2) / O000040000000, O000000000000 /
C      DATA LARGE(1),LARGE(2) / O377777777777, O777777777777 /
C      DATA RIGHT(1),RIGHT(2) / O170540000000, O000000000000 /
C      DATA DIVER(1),DIVER(2) / O170640000000, O000000000000 /
C      DATA LOG10(1),LOG10(2) / O177746420232, O411757177572 /, SC/987/
C
C     MACHINE CONSTANTS FOR THE VAX UNIX F77 COMPILER
C
C      DATA SMALL(1),SMALL(2) /        128,           0 /
C      DATA LARGE(1),LARGE(2) /     -32769,          -1 /
C      DATA RIGHT(1),RIGHT(2) /       9344,           0 /
C      DATA DIVER(1),DIVER(2) /       9472,           0 /
C      DATA LOG10(1),LOG10(2) /  546979738,  -805796613 /, SC/987/
C
C     MACHINE CONSTANTS FOR THE VAX-11 WITH
C     FORTRAN IV-PLUS COMPILER
C
C      DATA SMALL(1),SMALL(2) / Z00000080, Z00000000 /
C      DATA LARGE(1),LARGE(2) / ZFFFF7FFF, ZFFFFFFFF /
C      DATA RIGHT(1),RIGHT(2) / Z00002480, Z00000000 /
C      DATA DIVER(1),DIVER(2) / Z00002500, Z00000000 /
C      DATA LOG10(1),LOG10(2) / Z209A3F9A, ZCFF884FB /, SC/987/
C
C     MACHINE CONSTANTS FOR VAX/VMS VERSION 2.2
C
C      DATA SMALL(1),SMALL(2) /       '80'X,        '0'X /
C      DATA LARGE(1),LARGE(2) / 'FFFF7FFF'X, 'FFFFFFFF'X /
C      DATA RIGHT(1),RIGHT(2) /     '2480'X,        '0'X /
C      DATA DIVER(1),DIVER(2) /     '2500'X,        '0'X /
C      DATA LOG10(1),LOG10(2) / '209A3F9A'X, 'CFF884FB'X /, SC/987/
C
C  ***  ISSUE STOP 779 IF ALL DATA STATEMENTS ARE COMMENTED...
      IF (SC .NE. 987) THEN
         DMACH(1) = 1.D13
         IF (      SMALL(1) .EQ. 1117925532
     *       .AND. SMALL(2) .EQ. -448790528) THEN
*           *** IEEE BIG ENDIAN ***
            SMALL(1) = 1048576
            SMALL(2) = 0
            LARGE(1) = 2146435071
            LARGE(2) = -1
            RIGHT(1) = 1017118720
            RIGHT(2) = 0
            DIVER(1) = 1018167296
            DIVER(2) = 0
            LOG10(1) = 1070810131
            LOG10(2) = 1352628735
         ELSE IF ( SMALL(2) .EQ. 1117925532
     *       .AND. SMALL(1) .EQ. -448790528) THEN
*           *** IEEE LITTLE ENDIAN ***
            SMALL(2) = 1048576
            SMALL(1) = 0
            LARGE(2) = 2146435071
            LARGE(1) = -1
            RIGHT(2) = 1017118720
            RIGHT(1) = 0
            DIVER(2) = 1018167296
            DIVER(1) = 0
            LOG10(2) = 1070810131
            LOG10(1) = 1352628735
         ELSE IF ( SMALL(1) .EQ. -2065213935
     *       .AND. SMALL(2) .EQ. 10752) THEN
*               *** VAX WITH D_FLOATING ***
            SMALL(1) = 128
            SMALL(2) = 0
            LARGE(1) = -32769
            LARGE(2) = -1
            RIGHT(1) = 9344
            RIGHT(2) = 0
            DIVER(1) = 9472
            DIVER(2) = 0
            LOG10(1) = 546979738
            LOG10(2) = -805796613
         ELSE IF ( SMALL(1) .EQ. 1267827943
     *       .AND. SMALL(2) .EQ. 704643072) THEN
*               *** IBM MAINFRAME ***
            SMALL(1) = 1048576
            SMALL(2) = 0
            LARGE(1) = 2147483647
            LARGE(2) = -1
            RIGHT(1) = 856686592
            RIGHT(2) = 0
            DIVER(1) = 873463808
            DIVER(2) = 0
            LOG10(1) = 1091781651
            LOG10(2) = 1352628735
         ELSE IF ( SMALL(1) .EQ. 1120022684
     *       .AND. SMALL(2) .EQ. -448790528) THEN
*           *** CONVEX C-1 ***
            SMALL(1) = 1048576
            SMALL(2) = 0
            LARGE(1) = 2147483647
            LARGE(2) = -1
            RIGHT(1) = 1019215872
            RIGHT(2) = 0
            DIVER(1) = 1020264448
            DIVER(2) = 0
            LOG10(1) = 1072907283
            LOG10(2) = 1352628735
         ELSE IF ( SMALL(1) .EQ. 815547074
     *       .AND. SMALL(2) .EQ. 58688) THEN
*           *** VAX G-FLOATING ***
            SMALL(1) = 16
            SMALL(2) = 0
            LARGE(1) = -32769
            LARGE(2) = -1
            RIGHT(1) = 15552
            RIGHT(2) = 0
            DIVER(1) = 15568
            DIVER(2) = 0
            LOG10(1) = 1142112243
            LOG10(2) = 2046775455
         ELSE
            DMACH(2) = 1.D27 + 1
            DMACH(3) = 1.D27
            LARGE(2) = LARGE(2) - RIGHT(2)
            IF (LARGE(2) .EQ. 64 .AND. SMALL(2) .EQ. 0) THEN
               CRAY1(1) = 67291416
               DO 10 J = 1, 20
 10                              CRAY1(J+1) = CRAY1(J) + CRAY1(J)
               CRAY1(22) = CRAY1(21) + 321322
               DO 20 J = 22, 37
 20                              CRAY1(J+1) = CRAY1(J) + CRAY1(J)
               IF (CRAY1(38) .EQ. SMALL(1)) THEN
*                  *** CRAY ***
*                 SMALL(1) = 2332160919536140288
                  SMALL(1) = 2332160
                  SMALL(1) = 1000000*SMALL(1) + 919536
                  SMALL(1) = 1000000*SMALL(1) + 140288
                  SMALL(2) = 0
*                 LARGE(1) = 6917247552664371199
                  LARGE(1) = 6917247
                  LARGE(1) = 1000000*LARGE(1) + 552664
                  LARGE(1) = 1000000*LARGE(1) + 371199
*                 LARGE(2) = 281474976710654
                  LARGE(2) = 28147497
                  LARGE(2) = 10000000*LARGE(2) + 6710654
*                 RIGHT(1) = 4585649583081652224
                  RIGHT(1) = 4585649
                  RIGHT(1) = 1000000*RIGHT(1) + 583081
                  RIGHT(1) = 1000000*RIGHT(1) + 652224
                  RIGHT(2) = 0
*                 DIVER(1) = 4585931058058362880
                  DIVER(1) = 4585931
                  DIVER(1) = 1000000*DIVER(1) + 058058
                  DIVER(1) = 1000000*DIVER(1) + 362880
                  DIVER(2) = 0
*                 LOG10(1) = 4611574008272714703
                  LOG10(1) = 4611574
                  LOG10(1) = 1000000*LOG10(1) +   8272
                  LOG10(1) = 1000000*LOG10(1) + 714703
*                 LOG10(2) = 272234615232940
                  LOG10(2) = 27223461
                  LOG10(2) = 10000000*LOG10(2) + 5232940
               ELSE
                  WRITE(*,9000)
                  STOP 779
                  END IF
            ELSE
               WRITE(*,9000)
               STOP 779
               END IF
            END IF
         SC = 987
         END IF
C
C  ***  ISSUE STOP 778 IF ALL DATA STATEMENTS ARE OBVIOUSLY WRONG...
      IF (DMACH(4) .GE. 1.0D0) STOP 778
*C/6S
*C     IF (I .LT. 1  .OR.  I .GT. 5)
*C    1   CALL SETERR(24HD1MACH - I OUT OF BOUNDS,24,1,2)
*C/7S
*      IF (I .LT. 1  .OR.  I .GT. 5)
*     1   CALL SETERR('D1MACH - I OUT OF BOUNDS',24,1,2)
*C/
      IF (I .LT. 1 .OR. I .GT. 5) THEN
         WRITE(*,*) 'D1MACH(I): I =',I,' is out of bounds.'
         STOP
         END IF
      D1MACH = DMACH(I)
      RETURN
C/6S
C9000 FORMAT(/46H Adjust D1MACH by uncommenting data statements/
C    *30H appropriate for your machine.)
C/7S
 9000  FORMAT(/' Adjust D1MACH by uncommenting data statements'/
     *' appropriate for your machine.')
C/
C
* /* C source for D1MACH -- remove the * in column 1 */
*#include <stdio.h>
*#include <float.h>
*#include <math.h>
*
*double d1mach_(long *i)
*{
*     switch(*i){
*       case 1: return DBL_MIN;
*       case 2: return DBL_MAX;
*       case 3: return DBL_EPSILON/FLT_RADIX;
*       case 4: return DBL_EPSILON;
*       case 5: return log10(FLT_RADIX);
*       }
*
*     fprintf(stderr, "invalid argument: d1mach(%ld)\n", *i);
*     exit(1);
*     return 0; /* for compilers that complain of missing return values */
*     }
      END


      subroutine sprt(n, array)
* sprt prints an array.
      implicit double precision (a-h,o-z)
      dimension array(n)
      
      write(6,900) (array(i), i=1,n)
 900   format(1h ,(7(1pe11.3)))
      return
      end

      subroutine mprt(nrowd, nrow, ncol, array)
* mprt prints a matrix.
      implicit double precision (a-h,o-z)
      dimension array(nrowd, ncol)

      do 400 i = 1, nrow
         write(6,900) i,(array(i,j), j=1,ncol)
 400      continue
 900       format(1h ,i5,(6(1pe11.3)))
      return
      end


C
C
C-----------------------------------------------------------------------
C
C     THE AUGUST 27 1992 VERSION OF COLROW IN WHICH X IS NO LONGER
C     REQUIRED, WITH THE SOLUTION BEING RETURNED IN B, THE RIGHT
C     HAND SIDE.  IN ADDITION, ALL VARIABLES ARE EXPLICITLY DECLARED.
C     A PARAMETER "JOB" IS INCLUDED, TO SPECIFY WHICH OF A.X = B OR
C     TRANSPOSE(A).X = B IS TO BE SOLVED.

      SUBROUTINE COLROW(N,TOPBLK,NRWTOP,NOVRLP,ARRAY,NRWBLK,
     *             NCLBLK,NBLOKS,BOTBLK,NRWBOT,PIVOT,B,IFLAG,JOB)
C
C***************************************************************
C
C  THIS PROGRAM SOLVES ONE OF THE LINEAR SYSTEMS  A*X = B OR
C  TRANSPOSE(A)*X = B, WHERE  A IS AN ALMOST BLOCK DIAGONAL
C  MATRIX OF THE FORM
C
C               TOPBLK
C               ARRAY(1)
C                     ARRAY(2)
C                          .
C                             .
C                                .
C                                   .
C                                    ARRAY(NBLOKS)
C                                           BOTBLK
C
C  WHERE
C           TOPBLK IS  NRWTOP  BY NOVRLP
C           ARRAY(K), K=1,NBLOKS, ARE NRWBLK BY NRWBLK+NOVRLP
C           BOTBLK IS NRWBOT BY NOVRLP,
C  AND
C           NOVRLP = NRWTOP + NRWBOT
C  WITH
C           NOVRLP.LE.NRWBLK .
C
C  THE LINEAR SYSTEM IS OF ORDER  N = NBLOKS*NRWBLK + NOVRLP.
C
C  THE METHOD IMPLEMENTED IS BASED ON GAUSS ELIMINATION WITH
C  ALTERNATE ROW AND COLUMN ELIMINATION WITH PARTIAL PIVOTING,
C  WHICH PRODUCES A STABLE DECOMPOSITION OF THE MATRIX  A
C  WITHOUT INTRODUCING FILL-IN.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  TO OBTAIN A SINGLE PRECISION VERSION OF THIS PACKAGE, REMOVE
C  ALL DOUBLE PRECISION STATEMENTS.  THERE IS ONE SUCH STATEMENT
C  IN C O L R O W, THREE IN C R D C M P, AND TWO IN C R S O L V.
C  IN ADDITION, REFERENCES TO BUILT-IN FUNCTIONS DABS AND DMAX1
C  MUST BE REPLACED BY ABS AND AMAX1, RESPECTIVELY.  DABS OCCURS
C  NINE TIMES, IN C R D C M P.  DMAX1 OCCURS FOUR TIMES, IN
C  C R D C M P.  FINALLY, ZERO IS INITIALISED TO 0.D0 IN A
C  DATA STATEMENT IN C R D C M P.  THIS MUST BE REPLACED BY:
C               DATA ZERO/0.0/
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               *****  PARAMETERS  *****
C
C       *** ON ENTRY ...
C
C               N      - INTEGER
C                         THE ORDER OF THE LINEAR SYSTEM,
C                         GIVEN BY NBLOKS*NRWBLK + NOVRLP
C
C               TOPBLK - DOUBLE PRECISION(NRWTOP,NOVRLP)
C                         THE FIRST BLOCK OF THE ALMOST BLOCK
C                         DIAGONAL MATRIX A
C
C               NRWTOP - INTEGER
C                         NUMBER OF ROWS IN THE BLOCK TOPBLK
C
C               NOVRLP - INTEGER
C                         THE NUMBER OF COLUMNS IN WHICH SUCC-
C                         ESSIVE BLOCKS OVERLAP, WHERE
C                                NOVRLP = NRWTOP + NRWBOT
C
C               ARRAY  - DOUBLE PRECISION(NRWBLK,NCLBLK,NBLOKS)
C                         ARRAY(,,K) CONTAINS THE K-TH NRWBLK
C                         BY NCLBLK BLOCK OF THE MATRIX A
C
C               NRWBLK - INTEGER
C                         NUMBER OF ROWS IN K-TH BLOCK
C
C               NCLBLK - INTEGER
C                         NUMBER OF COLUMNS IN K-TH BLOCK
C
C               NBLOKS - INTEGER
C                         NUMBER OF NRWBLK BY NCLBLK BLOCKS IN
C                         THE MATRIX A
C
C               BOTBLK - DOUBLE PRECISION(NRWBOT,NOVRLP)
C                         THE LAST BLOCK OF THE MATRIX A
C
C               NRWBOT - INTEGER
C                         NUMBER OF ROWS IN THE BLOCK BOTBLK
C
C                PIVOT - INTEGER(N)
C                         WORK SPACE
C
C                    B - DOUBLE PRECISION(N)
C                         THE RIGHT HAND SIDE VECTOR
C
C               JOB    - INTEGER, INDICATING:
C                      = 0: SOLVE A*X = B;
C                      NON-ZERO: SOLVE TRANSPOSE(A)*X = B.
C
C       *** ON RETURN  ...
C
C               TOPBLK,ARRAY,BOTBLK - ARRAYS CONTAINING THE
C                        DESIRED DECOMPOSITION OF THE MATRIX A
C                        (IF IFLAG = 0)
C
C                PIVOT - INTEGER(N)
C                         RECORDS THE PIVOTING INDICES DETER-
C                         MINED IN THE DECOMPOSITION
C
C                    B - DOUBLE PRECISION(N)
C                         THE SOLUTION VECTOR (IF IFLAG = 0)
C
C               IFLAG  - INTEGER
C                         =  1, IF INPUT PARAMETERS ARE INVALID
C                         = -1, IF MATRIX IS SINGULAR
C                         =  0, OTHERWISE
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               *****  AUXILIARY PROGRAMS  *****
C
C       CRDCMP(N,TOPBLK,NRWTOP,NOVRLP,ARRAY,NRWBLK,NCLBLK,NBLOKS,
C    *     BOTBLK,NRWBOT,PIVOT,IFLAG)
C            - DECOMPOSES THE MATRIX  A  USING MODIFIED
C              ALTERNATE ROW AND COLUMN ELIMINATON WITH
C              PARTIAL PIVOTING, AND IS USED FOR THIS
C              PURPOSE IN C O L R O W.
C              THE ARGUMENTS ARE AS IN C O L R O W.
C
C       CRSLVE(TOPBLK,NRWTOP,NOVRLP,ARRAY,NRWBLK,NCLBLK,NBLOKS,
C    *     BOTBLK,NRWBOT,PIVOT,B,JOB)
C            - SOLVES THE SYSTEM A*X = B ONCE A IS DECOMPOSED.
C              THE ARGUMENTS ARE ALL AS IN C O L R O W.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C       THE SUBROUTINE  C O L R O W  AUTOMATICALLY SOLVES THE
C  INPUT SYSTEM WHEN IFLAG=0.  C O L R O W  IS CALLED ONLY ONCE
C  FOR A GIVEN SYSTEM. THE SOLUTION FOR A SEQUENCE OF P RIGHT
C  HAND SIDES CAN BE OBTAINED BY ONE CALL TO  C O L R O W  AND
C  P-1 CALLS TO CRSLVE ONLY. SINCE THE ARRAYS TOPBLK,ARRAY,
C  BOTBLK AND PIVOT CONTAIN THE DECOMPOSITION OF THE GIVEN
C  COEFFICIENT MATRIX AND PIVOTING INFORMATION ON RETURN FROM
C  C O L R O W , THEY MUST NOT BE ALTERED BETWEEN SUCCESSIVE
C  CALLS TO CRSLVE WITH THE SAME LEFT HAND SIDES. FOR THE
C  SAME REASON, IF THE USER WISHES TO SAVE THE COEFFICIENT
C  MATRIX, THE ARRAYS TOPBLK,ARRAY,BOTBLK MUST BE COPIED
C  BEFORE A CALL TO  C O L R O W .
C
C*************************************************************************
        IMPLICIT NONE
        DOUBLE PRECISION TOPBLK,ARRAY,BOTBLK,B
        INTEGER N,NRWTOP,NOVRLP,NRWBLK,NCLBLK,NBLOKS,NRWBOT,PIVOT(*),
     *          IFLAG,JOB, IDAMAX,i,IFAIL,j,jj
        DIMENSION TOPBLK(NRWTOP,*),ARRAY(NRWBLK,NCLBLK,*),
     *          BOTBLK(NRWBOT,*),B(*)

c        NOUT = 6

C       DO THE FACTORIZATION USING CRDCMP:
C
        CALL CRDCMP(N,TOPBLK,NRWTOP,NOVRLP,ARRAY,NRWBLK,NCLBLK,NBLOKS,
     *          BOTBLK,NRWBOT,PIVOT,IFLAG)
        IF(IFLAG.NE.0) RETURN
C
c     *****************solving the linear system********************
        job=0

        CALL CRSLVE(TOPBLK,NRWTOP,NOVRLP,ARRAY,NRWBLK,NCLBLK,NBLOKS,
     *          BOTBLK,NRWBOT,PIVOT,B,JOB)
        RETURN
        END


      SUBROUTINE INVERSE(N,TOPBLK,NRWTOP,NOVRLP,ARRAY,NRWBLK,NCLBLK,
     *          NBLOKS,BOTBLK,NRWBOT,PIVOT,INMAT)
C******************************************************************
C     INVERSE COMPUTES THE INVERSE OF A MATRIX 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        IMPLICIT NONE
        DOUBLE PRECISION TOPBLK,ARRAY,BOTBLK,WORK, 
     *          INMAT
        INTEGER N,NRWTOP,NOVRLP,NRWBLK,NCLBLK,NBLOKS,NRWBOT,PIVOT(*)
        DIMENSION TOPBLK(NRWTOP,*),ARRAY(NRWBLK,NCLBLK,*),
     *          BOTBLK(NRWBOT,*),WORK(N),
     *          INMAT(N,N) 
        INTEGER K,L,J

        do 300 k=1,N
           do 330 l=1,N
              WORK(l)=0.0d0
              if (l.eq.k) WORK(l)=1.0d0
 330                 continue

          CALL CRSLVE(TOPBLK,NRWTOP,NOVRLP,ARRAY,NRWBLK,
     *                 NCLBLK,NBLOKS,BOTBLK,NRWBOT,PIVOT,
     *                 WORK,0)

          do 440 l=1,N
             INMAT(l,k)=WORK(l)
 440                  continue
 300                        continue

          RETURN
          END

C******************************************************************



      SUBROUTINE CRDCMP(N,TOPBLK,NRWTOP,NOVRLP,ARRAY,NRWBLK,NCLBLK,
     *   NBLOKS,BOTBLK,NRWBOT,PIVOT,IFLAG)
C
C***************************************************************
C
C  C R D C M P DECOMPOSES THE ALMOST BLOCK DIAGONAL MATRIX A
C  USING MODIFIED ALTERNATE ROW AND COLUMN ELIMINATION WITH 
C  PARTIAL PIVOTING.  THE MATRIX  A  IS STORED IN THE ARRAYS
C  TOPBLK, ARRAY, AND BOTBLK. 
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               *****  PARAMETERS  *****
C
C       *** ON ENTRY ...
C
C               N      - INTEGER
C                         THE ORDER OF THE LINEAR SYSTEM,
C                         GIVEN BY NBLOKS*NRWBLK + NOVRLP
C
C               TOPBLK - DOUBLE PRECISION(NRWTOP,NOVRLP)
C                         THE FIRST BLOCK OF THE ALMOST BLOCK
C                         DIAGONAL MATRIX A TO BE DECOMPOSED
C
C               NRWTOP - INTEGER
C                         NUMBER OF ROWS IN THE BLOCK TOPBLK
C
C               NOVRLP - INTEGER
C                         THE NUMBER OF COLUMNS IN WHICH SUCC-
C                         ESSIVE BLOCKS OVERLAP, WHERE
C                                NOVRLP = NRWTOP + NRWBOT
C
C               ARRAY  - DOUBLE PRECISION(NRWBLK,NCLBLK,NBLOKS)
C                         ARRAY(,,K) CONTAINS THE K-TH NRWBLK
C                         BY NCLBLK BLOCK OF THE MATRIX A
C
C               NRWBLK - INTEGER
C                         NUMBER OF ROWS IN K-TH BLOCK
C
C               NCLBLK - INTEGER
C                         NUMBER OF COLUMNS IN K-TH BLOCK
C
C               NBLOKS - INTEGER
C                         NUMBER OF NRWBLK BY NCLBLK BLOCKS IN
C                         THE MATRIX A
C
C               BOTBLK - DOUBLE PRECISION(NRWBOT,NOVRLP)
C                         THE LAST BLOCK OF THE MATRIX A
C
C               NRWBOT - INTEGER
C                         NUMBER OF ROWS IN THE BLOCK BOTBLK
C
C                PIVOT - INTEGER(N)
C                         WORK SPACE
C
C       *** ON RETURN  ...
C
C               TOPBLK,ARRAY,BOTBLK - ARRAYS CONTAINING THE 
C                        DESIRED DECOMPOSITION OF THE MATRIX A
C                        (IF IFLAG = 0) 
C
C                PIVOT - INTEGER(N)
C                         RECORDS THE PIVOTING INDICES DETER-
C                         MINED IN THE DECOMPOSITION
C
C               IFLAG  - INTEGER
C                         =  1, IF INPUT PARAMETERS ARE INVALID
C                         = -1, IF MATRIX IS SINGULAR
C                         =  0, OTHERWISE
C
C***************************************************************
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER PIVOT(*)
      DIMENSION TOPBLK(NRWTOP,*),ARRAY(NRWBLK,NCLBLK,*),BOTBLK(NRWBOT,*)
      DATA ZERO / 0.0D+0 /
      Common/Mchprs/Flmin, Flmax, Epsmch
C
C***************************************************************
C
C          ****  DEFINE THE CONSDTANTS USED THROUGHOUT  **** 
C
C***************************************************************
C
      IFLAG = 0
      PIVMAX = ZERO 
      NRWTP1 = NRWTOP+1
      NROWEL = NRWBLK-NRWTOP
      NRWEL1 = NROWEL+1
      NVRLP0 = NOVRLP-1
      Pivtol = 10.0d+0*Epsmch
C
C***************************************************************
C
C          ****  CHECK VALIDITY OF THE INPUT PARAMETERS.... 
C
C               IF PARAMETERS ARE INVALID THEN TERMINATE AT 10;
C                                         ELSE CONTINUE AT 100.
C
C***************************************************************
C
      IF (N.NE.NBLOKS*NRWBLK+NOVRLP) GO TO 10
      IF (NOVRLP.NE.NRWTOP+NRWBOT) GO TO 10
      IF (NCLBLK.NE.NOVRLP+NRWBLK) GO TO 10
      IF (NOVRLP.GT.NRWBLK) GO TO 10
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C          PARAMETERS ARE ACCEPTABLE - CONTINUE AT 100.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      GO TO 20
 10    CONTINUE
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C          PARAMETERS ARE INVALID.  SET IFLAG = 1, AND TERMINATE
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      IFLAG = 1
      RETURN
 20    CONTINUE
C
C***************************************************************
C
C               ****  FIRST, IN TOPBLK....
C
C***************************************************************
C
C          ***  APPLY NRWTOP COLUMN ELIMINATIONS WITH COLUMN
C                 PIVOTING ....
C
C***************************************************************
C
      DO 110 I = 1, NRWTOP
         IPLUS1 = I+1
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               DETERMINE COLUMN PIVOT AND PIVOT INDEX
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
         IPVT = I
         COLMAX = DABS(TOPBLK(I,I))
         DO 30 J = IPLUS1, NOVRLP
            TEMPIV = DABS(TOPBLK(I,J))
            IF (TEMPIV.LE.COLMAX) GO TO 30
            IPVT = J
            COLMAX = TEMPIV
 30             CONTINUE
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               TEST FOR SINGULARITY: 
C
C                       IF SINGULAR THEN TERMINATE AT 1000; 
C                                   ELSE CONTINUE.
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C

         If (Colmax .le. Pivtol) Then
            Iflag = -1
            Return
         Endif

c         IF (PIVMAX+COLMAX.EQ.PIVMAX) GO TO 410
         PIVMAX = DMAX1(COLMAX,PIVMAX)
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               IF NECESSARY INTERCHANGE COLUMNS
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
         PIVOT(I) = IPVT
         IF (IPVT.EQ.I) GO TO 60
         DO 40 L = I, NRWTOP
            SWAP = TOPBLK(L,IPVT)
            TOPBLK(L,IPVT) = TOPBLK(L,I)
            TOPBLK(L,I) = SWAP
 40             CONTINUE
         DO 50 L = 1, NRWBLK
            SWAP = ARRAY(L,IPVT,1)
            ARRAY(L,IPVT,1) = ARRAY(L,I,1)
            ARRAY(L,I,1) = SWAP
 50             CONTINUE
 60                 CONTINUE
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               COMPUTE MULTIPLIERS AND PERFORM COLUMN
C                       ELIMINATION
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
         COLPIV = TOPBLK(I,I) 
         DO 100 J = IPLUS1, NOVRLP
            COLMLT = TOPBLK(I,J)/COLPIV 
            TOPBLK(I,J) = COLMLT
            IF (IPLUS1.GT.NRWTOP) GO TO 80
            DO 70 L = IPLUS1, NRWTOP
               TOPBLK(L,J) = TOPBLK(L,J)-COLMLT*TOPBLK(L,I) 
 70                   CONTINUE
 80                          CONTINUE
            DO 90 L = 1, NRWBLK
               ARRAY(L,J,1) = ARRAY(L,J,1)-COLMLT*ARRAY(L,I,1)
 90                   CONTINUE
 100                      CONTINUE
 110                       CONTINUE
C
C***************************************************************
C
C          ****  IN EACH BLOCK ARRAY(,,K)....
C
C***************************************************************
C
      INCR = 0
      DO 320 K = 1, NBLOKS
         KPLUS1 = K+1
C
C          *****************************************************
C
C          ***  FIRST APPLY NRWBLK-NRWTOP ROW ELIMINATIONS WITH
C                       ROW PIVOTING....
C
C          *****************************************************
C
         DO 180 J = NRWTP1, NRWBLK
            JPLUS1 = J+1
            JMINN = J-NRWTOP
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               DETERMINE ROW PIVOT AND PIVOT INDEX
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
            IPVT = JMINN
            ROWMAX = DABS(ARRAY(JMINN,J,K))
            LOOP = JMINN+1
            DO 120 I = LOOP, NRWBLK
               TEMPIV = DABS(ARRAY(I,J,K))
               IF (TEMPIV.LE.ROWMAX) GO TO 120
               IPVT = I
               ROWMAX = TEMPIV
 120                  CONTINUE
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               TEST FOR SINGULARITY: 
C
C                       IF SINGULAR THEN TERMINATE AT 1000; 
C                                   ELSE CONTINUE.
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
         If (Rowmax .le. Pivtol) Then
            Iflag = -1
            Return
         Endif

c            IF (PIVMAX+ROWMAX.EQ.PIVMAX) GO TO 410
            PIVMAX = DMAX1(ROWMAX,PIVMAX)
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               IF NECESSARY INTERCHANGE ROWS
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
            INCRJ = INCR+J
            PIVOT(INCRJ) = INCR+IPVT+NRWTOP
            IF (IPVT.EQ.JMINN) GO TO 140
            DO 130 L = J, NCLBLK
               SWAP = ARRAY(IPVT,L,K)
               ARRAY(IPVT,L,K) = ARRAY(JMINN,L,K) 
               ARRAY(JMINN,L,K) = SWAP
 130                  CONTINUE
 140                         CONTINUE
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               COMPUTE MULTIPLERS
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
            ROWPIV = ARRAY(JMINN,J,K)
            DO 150 I = LOOP, NRWBLK
               ARRAY(I,J,K) = ARRAY(I,J,K)/ROWPIV 
 150                  CONTINUE
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               PERFORM ROW ELIMINATION WITH COLUMN INDEXING
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
            DO 170 L = JPLUS1, NCLBLK
               ROWMLT = ARRAY(JMINN,L,K)
               DO 160 I = LOOP, NRWBLK
                  ARRAY(I,L,K) = ARRAY(I,L,K)-ROWMLT*ARRAY(I,J,K)
 160                        CONTINUE
 170                               CONTINUE
 180                                   CONTINUE
C
C          *****************************************************
C
C          ***  NOW APPLY NRWTOP COLUMN ELIMINATIONS WITH
C                      COLUMN PIVOTING....
C
C          *****************************************************
C
         DO 310 I = NRWEL1, NRWBLK
            IPLUSN = I+NRWTOP 
            IPLUS1 = I+1
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               DETERMINE COLUMN PIVOT AND PIVOT INDEX
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
            IPVT = IPLUSN
            COLMAX = DABS(ARRAY(I,IPVT,K))
            LOOP = IPLUSN+1
            DO 190 J = LOOP, NCLBLK
               TEMPIV = DABS(ARRAY(I,J,K))
               IF (TEMPIV.LE.COLMAX) GO TO 190
               IPVT = J
               COLMAX = TEMPIV
 190                  CONTINUE
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               TEST FOR SINGULARITY: 
C
C                       IF SINGULAR THEN TERMINATE AT 1000; 
C                                   ELSE CONTINUE.
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C

         If (Colmax .le. Pivtol) Then
            Iflag = -1
            Return
         Endif

c            IF (PIVMAX+COLMAX.EQ.PIVMAX) GO TO 410
            PIVMAX = DMAX1(COLMAX,PIVMAX)
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               IF NECESSARY INTERCHANGE COLUMNS
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
            INCRN = INCR+IPLUSN
            PIVOT(INCRN) = INCR+IPVT
            IRWBLK = IPLUSN-NRWBLK
            IF (IPVT.EQ.IPLUSN) GO TO 240
            DO 200 L = I, NRWBLK
               SWAP = ARRAY(L,IPVT,K)
               ARRAY(L,IPVT,K) = ARRAY(L,IPLUSN,K)
               ARRAY(L,IPLUSN,K) = SWAP 
 200                  CONTINUE
            IPVBLK = IPVT-NRWBLK
            IF (K.EQ.NBLOKS) GO TO 220
            DO 210 L = 1, NRWBLK
               SWAP = ARRAY(L,IPVBLK,KPLUS1)
               ARRAY(L,IPVBLK,KPLUS1) = ARRAY(L,IRWBLK,KPLUS1)
               ARRAY(L,IRWBLK,KPLUS1) = SWAP
 210                  CONTINUE
            GO TO 240
 220               CONTINUE
            DO 230 L = 1, NRWBOT
               SWAP = BOTBLK(L,IPVBLK)
               BOTBLK(L,IPVBLK) = BOTBLK(L,IRWBLK)
               BOTBLK(L,IRWBLK) = SWAP
 230                  CONTINUE
 240                         CONTINUE
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               COMPUTE MULTIPLIERS AND PERFORM COLUMN
C                       ELIMINATION
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
            COLPIV = ARRAY(I,IPLUSN,K)
            DO 300 J = LOOP, NCLBLK
               COLMLT = ARRAY(I,J,K)/COLPIV
               ARRAY(I,J,K) = COLMLT
               IF (I.EQ.NRWBLK) GO TO 260
               DO 250 L = IPLUS1, NRWBLK
                  ARRAY(L,J,K) = ARRAY(L,J,K)-COLMLT*ARRAY(L,IPLUSN,K)
 250                        CONTINUE
 260                                  CONTINUE
               JRWBLK = J-NRWBLK
               IF (K.EQ.NBLOKS) GO TO 280
               DO 270 L = 1, NRWBLK
                  ARRAY(L,JRWBLK,KPLUS1) = ARRAY(L,JRWBLK,KPLUS1)-COLMLT
     *               *ARRAY(L,IRWBLK,KPLUS1)
 270                        CONTINUE
               GO TO 300
 280                     CONTINUE
               DO 290 L = 1, NRWBOT
                  BOTBLK(L,JRWBLK) = BOTBLK(L,JRWBLK)-COLMLT*BOTBLK(L,
     *               IRWBLK)
 290                        CONTINUE
 300                               CONTINUE
 310                                   CONTINUE
         INCR = INCR+NRWBLK
 320      CONTINUE
C
C***************************************************************
C
C          ****  FINALLY, IN BOTBLK.... 
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C          ***  APPLY NRWBOT ROW ELIMINATIONS WITH ROW
C                  PIVOTING....
C
C               IF BOT HAS JUST ONE ROW GO TO 500 
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C

      IF (NRWBOT.EQ.1) GO TO 400
      DO 390 J = NRWTP1, NVRLP0
         JPLUS1 = J+1
         JMINN = J-NRWTOP
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               DETERMINE ROW PIVOT AND PIVOT INDEX
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
         IPVT = JMINN
         ROWMAX = DABS(BOTBLK(JMINN,J))
         LOOP = JMINN+1
         DO 330 I = LOOP, NRWBOT
            TEMPIV = DABS(BOTBLK(I,J))
            IF (TEMPIV.LE.ROWMAX) GO TO 330
            IPVT = I
            ROWMAX = TEMPIV
 330            CONTINUE
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               TEST FOR SINGULARITY: 
C
C                       IF SINGULAR THEN TERMINATE AT 1000; 
C                                   ELSE CONTINUE.
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
         If (Rowmax .le. Pivtol) Then
            Iflag = -1
            Return
         Endif

c         IF (PIVMAX+ROWMAX.EQ.PIVMAX) GO TO 410
         PIVMAX = DMAX1(ROWMAX,PIVMAX)
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               IF NECESSARY INTERCHANGE ROWS
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
         INCRJ = INCR+J
         PIVOT(INCRJ) = INCR+IPVT+NRWTOP
         IF (IPVT.EQ.JMINN) GO TO 350
         DO 340 L = J, NOVRLP 
            SWAP = BOTBLK(IPVT,L)
            BOTBLK(IPVT,L) = BOTBLK(JMINN,L)
            BOTBLK(JMINN,L) = SWAP
 340            CONTINUE
 350                CONTINUE
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               COMPUTE MULTIPLIERS
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
         ROWPIV = BOTBLK(JMINN,J)
         DO 360 I = LOOP, NRWBOT
            BOTBLK(I,J) = BOTBLK(I,J)/ROWPIV
 360            CONTINUE
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               PERFORM ROW ELIMINATION WITH COLUMN INDEXING
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
         DO 380 L = JPLUS1, NOVRLP
            ROWMLT = BOTBLK(JMINN,L)
            DO 370 I = LOOP, NRWBOT
               BOTBLK(I,L) = BOTBLK(I,L)-ROWMLT*BOTBLK(I,J)
 370                  CONTINUE
 380                      CONTINUE
 390                       CONTINUE
 400                        CONTINUE
C
C***************************************************************
C
C          DONE PROVIDED THE LAST ELEMENT IS NOT ZERO
C
C***************************************************************
C
       if (DABS(BOTBLK(NRWBOT,NOVRLP)).gt.pivtol) return                      
c      IF (PIVMAX+DABS(BOTBLK(NRWBOT,NOVRLP)).NE.PIVMAX) RETURN

C
C***************************************************************
C
C       ****  MATRIX IS SINGULAR - SET IFLAG = - 1.
C                                  TERMINATE AT 1000.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
 410   CONTINUE
      IFLAG = -1
      RETURN
      END 

        SUBROUTINE CRSLVE(TOPBLK,NRWTOP,NOVRLP,ARRAY,NRWBLK,
     *             NCLBLK,NBLOKS,BOTBLK,NRWBOT,PIVOT,B,JOB)
C
C***************************************************************
C
C  C R S L V E  SOLVES THE LINEAR SYSTEM
C                       A*X = B
C  USING THE DECOMPOSITION ALREADY GENERATED IN  C R D C M P.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               *****  PARAMETERS  *****
C
C       *** ON ENTRY  ...
C
C               TOPBLK - DOUBLE PRECISION(NRWTOP,NOVRLP)
C                         OUTPUT FROM  C R D C M P
C
C               NOVRLP - INTEGER
C                         THE NUMBER OF COLUMNS IN WHICH SUCC-
C                         ESSIVE BLOCKS OVERLAP, WHERE
C                                NOVRLP = NRWTOP + NRWBOT
C
C               NRWTOP - INTEGER
C                         NUMBER OF ROWS IN THE BLOCK TOPBLK
C
C               ARRAY  - DOUBLE PRECISION(NRWBLK,NCLBLK,NBLOKS)
C                         OUTPUT FROM  C R D C M P
C
C               NRWBLK - INTEGER
C                         NUMBER OF ROWS IN K-TH BLOCK
C
C               NCLBLK - INTEGER
C                         NUMBER OF COLUMNS IN K-TH BLOCK
C
C               NBLOKS - INTEGER
C                         NUMBER OF NRWBLK BY NCLBLK BLOCKS IN
C                         THE MATRIX A
C
C               BOTBLK - DOUBLE PRECISION(NRWBOT,NOVRLP)
C                         OUTPUT FROM  C R D C M P
C
C               NRWBOT - INTEGER
C                         NUMBER OF ROWS IN THE BLOCK BOTBLK
C
C                PIVOT - INTEGER(N)
C                         THE PIVOT VECTOR FROM  C R D C M P
C
C                    B - DOUBLE PRECISION(N)
C                         THE RIGHT HAND SIDE VECTOR
C
C               JOB    - INTEGER, INDICATING:
C                      = 0: SOLVE A*X = B;
C                      NON-ZERO: SOLVE TRANSPOSE(A)*X = B.
C
C       *** ON RETURN  ...
C
C                    B - DOUBLE PRECISION(N)
C                         THE SOLUTION VECTOR
C
C***************************************************************
C
        IMPLICIT NONE
        DOUBLE PRECISION TOPBLK,ARRAY,BOTBLK,B
        DOUBLE PRECISION DOTPRD,BJ,XINCRJ,BINCRJ,SWAP,BI
        INTEGER NRWTOP,NOVRLP,NRWBLK,NCLBLK,NBLOKS,NRWBOT,PIVOT(*),
     *          JOB
        DIMENSION TOPBLK(NRWTOP,*),ARRAY(NRWBLK,NCLBLK,*),
     *          BOTBLK(NRWBOT,*),B(*)
        INTEGER NRWTP1,NRWBK1,NVRLP1,NRWBT1,NROWEL,NVRLP0,NBLKS1,
     *          NBKTOP,J,I,LOOP,INCR,INCRJ,INCRI,JPIVOT,JRWTOP,
     *          LL,L1,IPLUSN,INCRN,NRWTP0,NRWEL1,K,INCRTP,NRWBTL,
     *          IPVTN,NRWELL,IPVTI,L
C
C***************************************************************
C
C          ****  DEFINE THE CONSTANTS USED THROUGHOUT  ****
C
C***************************************************************
C
        NRWTP1 = NRWTOP+1
        NRWBK1 = NRWBLK+1
        NVRLP1 = NOVRLP+1
        NRWTP0 = NRWTOP-1
        NRWBT1 = NRWBOT+1
        NROWEL = NRWBLK-NRWTOP
        NRWEL1 = NROWEL+1
        NVRLP0 = NOVRLP-1
        NBLKS1 = NBLOKS+1
        NBKTOP = NRWBLK+NRWTOP
C
C       IF JOB IS NON-ZERO, TRANSFER TO THE SECTION DEALING WITH
C       TRANSPOSE(A)*X = B.
C
        IF ( JOB .NE. 0 ) GO TO 530
C
C***************************************************************
C
C               ****  FORWARD RECURSION  ****
C
C***************************************************************
C
C          ***  FIRST, IN TOPBLK....
C
C       CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               FORWARD SOLUTION
C
C       CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
        DO 130 J = 1,NRWTOP
           B(J) = B(J)/TOPBLK(J,J)
           IF(J.EQ.NRWTOP)GO TO 120
              BJ = -B(J)
              LOOP = J+1
              DO 110 I = LOOP,NRWTOP
                 B(I) = B(I)+TOPBLK(I,J)*BJ
 110          CONTINUE
 120       CONTINUE
 130    CONTINUE
C
C       ********************************************************
C
C          ***  IN EACH BLOCK ARRAY(,,K)....
C
C       CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
        INCR = 0
        DO 280 K = 1,NBLOKS
           INCRTP = INCR+NRWTOP
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               FORWARD MODIFICATION
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
           DO 220 J = 1,NRWTOP
              INCRJ = INCR+J
              XINCRJ = -B(INCRJ)
              DO 210 I = 1,NRWBLK
                 INCRI = INCRTP+I
                 B(INCRI) = B(INCRI)+ARRAY(I,J,K)*XINCRJ
 210          CONTINUE
 220       CONTINUE
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               FORWARD ELIMINATION
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
           DO 240 J = NRWTP1,NRWBLK
              INCRJ = INCR+J
              JPIVOT = PIVOT(INCRJ)
              IF(JPIVOT.EQ.INCRJ)GO TO 225
                 SWAP = B(INCRJ)
                 B(INCRJ) = B(JPIVOT)
                 B(JPIVOT) = SWAP
 225          CONTINUE
              BINCRJ = -B(INCRJ)
              LOOP = J-NRWTP0
              DO 230 I = LOOP,NRWBLK
                 INCRI = INCRTP+I
                 B(INCRI) = B(INCRI)+ARRAY(I,J,K)*BINCRJ
 230          CONTINUE
 240       CONTINUE
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               FORWARD SOLUTION
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
           DO 270 J = NRWBK1,NBKTOP
              INCRJ = INCR+J
              JRWTOP = J -NRWTOP
              B(INCRJ) = B(INCRJ)/ARRAY(JRWTOP,J,K)
              IF(J.EQ.NBKTOP)GO TO 260
                 XINCRJ = -B(INCRJ)
                 LOOP = J-NRWTP0
                 DO 250 I = LOOP,NRWBLK
                    INCRI = INCRTP+I
                    B(INCRI) = B(INCRI)+ARRAY(I,J,K)*XINCRJ
 250             CONTINUE
 260          CONTINUE
 270       CONTINUE
           INCR = INCR+NRWBLK
 280       CONTINUE
C
C       ********************************************************
C
C          ***  FINALLY, IN BOTBLK....
C
C       CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               FORWARD MODIFICATION
C
C       CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
        INCRTP = INCR+NRWTOP
        DO 320 J = 1,NRWTOP
           INCRJ = INCR+J
           XINCRJ = -B(INCRJ)
           DO 310 I = 1,NRWBOT
              INCRI = INCRTP+I
              B(INCRI) = B(INCRI)+BOTBLK(I,J)*XINCRJ
 310       CONTINUE
 320     CONTINUE
C
C       CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               FORWARD ELIMINATION
C
C       CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
        IF(NRWBOT.EQ.1)GO TO 350
           DO 340 J = NRWTP1,NVRLP0
              INCRJ = INCR+J
              JPIVOT = PIVOT(INCRJ)
              IF(JPIVOT.EQ.INCRJ)GO TO 325
                 SWAP = B(INCRJ)
                 B(INCRJ) = B(JPIVOT)
                 B(JPIVOT) = SWAP
 325          CONTINUE
              BINCRJ = -B(INCRJ)
              LOOP = J-NRWTP0
              DO 330 I = LOOP,NRWBOT
                 INCRI = INCRTP+I
                 B(INCRI) = B(INCRI)+BOTBLK(I,J)*BINCRJ
 330          CONTINUE
 340          CONTINUE
 350       CONTINUE
C
C***************************************************************
C
C               ****  BACKWARD RECURSION  ****
C
C***************************************************************
C
C          ***  FIRST IN BOTBLK....
C
C       CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               BACKWARD SOLUTION
C
C       CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
        DO 430 LL = 1,NRWBOT
           J = NVRLP1-LL
           INCRJ = INCR+J
           NRWBTL = NRWBT1-LL
           B(INCRJ) = B(INCRJ)/BOTBLK(NRWBTL,J)
           IF(LL.EQ.NRWBOT)GO TO 420
              XINCRJ = -B(INCRJ)
              LOOP = NRWBOT-LL
              DO 410 I = 1,LOOP
                 INCRI = INCRTP+I
                 B(INCRI) = B(INCRI)+BOTBLK(I,J)*XINCRJ
 410          CONTINUE
 420       CONTINUE
 430    CONTINUE
C
C       ********************************************************
C
C          ***  THEN IN EACH BLOCK ARRAY(,,K)....
C
C       CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
        DO 490 L = 1,NBLOKS
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               BACKWARD ELIMINATION
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
           K = NBLKS1-L
           INCR = INCR-NRWBLK
           DO 450 L1 = NRWEL1,NRWBLK
              I = NRWBLK+NRWEL1-L1
              IPLUSN = I+NRWTOP
              LOOP = IPLUSN+1
              INCRN = INCR+IPLUSN
              DOTPRD = B(INCRN)
              DO 440 J = LOOP,NCLBLK
                 INCRJ = INCR+J
                 DOTPRD = DOTPRD-ARRAY(I,J,K)*B(INCRJ)
 440          CONTINUE
              B(INCRN) = DOTPRD
              IPVTN = PIVOT(INCRN)
              IF(INCRN.EQ.IPVTN)GO TO 445
                 SWAP = B(INCRN)
                 B(INCRN) = B(IPVTN)
                 B(IPVTN) = SWAP
 445          CONTINUE
 450      CONTINUE
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               BACKWARD MODIFICATION
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
           INCRTP = INCR+NRWTOP
           DO 460 J = NRWBK1,NCLBLK
              INCRJ = INCR+J
              XINCRJ = -B(INCRJ)
              DO 455 I = 1,NROWEL
                 INCRI = INCRTP+I
                 B(INCRI) = B(INCRI)+ARRAY(I,J,K)*XINCRJ
 455          CONTINUE
 460       CONTINUE
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               BACKWARD SOLUTION
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
           DO 480 LL = 1,NROWEL
              J = NRWBK1-LL
              INCRJ = INCR+J
              NRWELL = NRWEL1-LL
              B(INCRJ) = B(INCRJ)/ARRAY(NRWELL,J,K)
              IF(LL.EQ.NROWEL)GO TO 470
                 XINCRJ = -B(INCRJ)
                 LOOP = NROWEL-LL
                 DO 465 I = 1,LOOP
                    INCRI = INCRTP+I
                    B(INCRI) = B(INCRI)+ARRAY(I,J,K)*XINCRJ
 465             CONTINUE
 470          CONTINUE
 480       CONTINUE
 490      CONTINUE
C
C       ********************************************************
C
C          ***  IN TOPBLK FINISH WITH....
C
C       CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               BACKWARD ELIMINATION
C
C       CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
        DO 520 L = 1,NRWTOP
           I = NRWTP1-L
           LOOP = I+1
           DOTPRD = B(I)
           DO 510 J = LOOP,NOVRLP
              DOTPRD = DOTPRD-TOPBLK(I,J)*B(J)
 510       CONTINUE
           B(I) = DOTPRD
           IPVTI = PIVOT(I)
           IF(I.EQ.IPVTI)GO TO 515
                 SWAP = B(I)
                 B(I) = B(IPVTI)
                 B(IPVTI) = SWAP
 515       CONTINUE
 520    CONTINUE
C
C       RETURN FROM THE SOLUTION OF A.X = B.
        RETURN
C
C       IF JOB IS NON-ZERO, SOLVE TRANSPOSE(A)*X = B:
C
 530           CONTINUE

C       FIRST, FORWARD ELIMINATION OF RHS USING TRANSPOSE(U).

        DO 540 I = 1,NRWTOP
           IPVTI = PIVOT(I)
           IF ( I .NE. IPVTI ) THEN
              SWAP = B(I)
              B(I) = B(IPVTI)
              B(IPVTI) = SWAP
           ENDIF
           BI = -B(I)
           LOOP = I+1
           DO 535 J = LOOP,NOVRLP
              B(J) = B(J) + BI*TOPBLK(I,J)
 535       CONTINUE
 540    CONTINUE

C       IN EACH BLOCK, K = 1,..,NBLOKS:

        INCR = NRWTOP
        DO 590 K = 1,NBLOKS

C          FIRST, THE FORWARD SOLUTION.

           DO 550 J = 1,NROWEL
              INCRJ = INCR + J
              DO 545 I = 1,J-1
                 B(INCRJ) = B(INCRJ) - ARRAY(I,NRWTOP+J,K)*B(INCR+I)
 545          CONTINUE
              B(INCRJ) = B(INCRJ)/ARRAY(J,NRWTOP+J,K)
 550       CONTINUE

C           FORWARD MODIFICATION.

            DO 570 I = 1,NOVRLP
               INCRI = INCR + NROWEL + I
               LOOP = NRWBLK + I
               DO 560 J = 1,NROWEL
                  INCRJ = INCR + J
                  B(INCRI) = B(INCRI) - ARRAY(J,LOOP,K)*B(INCRJ)
 560           CONTINUE
 570        CONTINUE

C           NOW, FORWARD ELIMINATION OF RHS USING TRANSPOSE(U). THIS
C           CORRESPONDS TO THE LOOP 540 ABOVE.

            INCR = INCR + NROWEL
            DO 580 I = 1,NRWTOP
               INCRI = INCR + I
               IPVTI = PIVOT(INCRI)
               IF ( INCRI .NE. IPVTI ) THEN
                  SWAP = B(INCRI)
                  B(INCRI) = B(IPVTI)
                  B(IPVTI) = SWAP
               ENDIF
               LOOP = NROWEL + I
               BI = -B(INCRI)
               DO 575 J = I+1,NOVRLP
                  INCRJ = INCR+J
                  L = NRWBLK + J


                  B(INCRJ) = B(INCRJ) + BI*ARRAY(LOOP,L,K)
 575           CONTINUE
 580        CONTINUE
            INCR = INCR + NRWTOP
 590      CONTINUE

C       FINALLY, FINISH WITH NRWBOT SOLUTIONS:

        DO 600 J = 1,NRWBOT
           INCRJ = INCR + J
           DO 595 I = 1,J-1
              B(INCRJ) = B(INCRJ) - BOTBLK(I,J+NRWTOP)*B(INCR+I)
 595       CONTINUE
           B(INCRJ) = B(INCRJ)/BOTBLK(J,J+NRWTOP)
 600    CONTINUE


C       NOW, THE BACKWARD PASS:


C       FIRST, BACKWARD SOLUTION IN BOTBLK:

        INCRJ = INCR + NRWBOT
        DO 610 J = 1,NRWBOT-1
           INCRJ = INCRJ - 1
           DO 605 I = NRWBOT-J+1,NRWBOT
              INCRI = INCR + I
              B(INCRJ) = B(INCRJ) - BOTBLK(I,NOVRLP-J)*B(INCRI)
 605       CONTINUE

           IF ( INCRJ .NE. PIVOT(INCRJ) ) THEN
              SWAP = B(INCRJ)
              B(INCRJ) = B(PIVOT(INCRJ))
              B(PIVOT(INCRJ)) = SWAP
           ENDIF
 610    CONTINUE

C       NOW DO THE DEFERRED OPERATIONS IN BOTBLOK:

        DO 620 J = 1,NRWTOP
           INCRJ = INCR - J + 1
           DO 615 I = 1,NRWBOT
              INCRI = INCR + I
              B(INCRJ) = B(INCRJ) - BOTBLK(I,NRWTP1-J)*B(INCRI)
 615       CONTINUE
 620    CONTINUE


C       NOW, IN EACH BLOCK, K = NBLOKS,..,1:
        DO 800 K = NBLOKS,1,-1

C          FIRST, THE BACKSUBSTITUIONS:

           DO 630 J = 1,NRWTOP
              INCRJ = INCR - J + 1
              LOOP = NBKTOP - J + 1
              DO 625 I = 1,J-1
                 INCRI = INCR - I + 1
                 B(INCRJ) = B(INCRJ) - ARRAY(NRWBLK-I+1,LOOP,K)*B(INCRI)
 625          CONTINUE
              B(INCRJ) = B(INCRJ)/ARRAY(NRWBLK-J+1,LOOP,K)
 630       CONTINUE

C          THEN THE BACKWARD SOLUTION IN THE KTH BLOCK:

           DO 650 J = 1,NROWEL
              INCRJ = INCR - NRWTOP -J + 1
              DO 645 I = 1,J+NRWTOP-1
                 INCRI = INCRJ + I
                 B(INCRJ) = B(INCRJ) -
     *           ARRAY(NRWBLK-NRWTOP-J+1+I,NRWBLK-J+1,K)*B(INCRI)
 645          CONTINUE
              IF ( INCRJ .NE. PIVOT(INCRJ) ) THEN
                 SWAP = B(INCRJ)
                 B(INCRJ) = B(PIVOT(INCRJ))
                 B(PIVOT(INCRJ)) = SWAP
              ENDIF
 650       CONTINUE

C          NOW, THE DEFERRED OPERATIONS ON B:

           INCR = INCR - NRWBLK
           DO 660 J = 1,NRWTOP
              INCRJ = INCR + J - NRWTOP
              DO 655 I = 1,NRWBLK
                 INCRI = INCR + I
                 B(INCRJ) = B(INCRJ) - ARRAY(I,J,K)*B(INCRI)
 655          CONTINUE
 660        CONTINUE
 800    CONTINUE

C       FINALLY, THE LAST SET OF BACK-SUBSTITUTIONS IN TOPBLK:

        DO 900 J = 1,NRWTOP
           INCRJ = NRWTOP -J + 1
           DO 850 I = INCRJ+1,NRWTOP
              B(INCRJ) = B(INCRJ) - TOPBLK(I,INCRJ)*B(I)
 850       CONTINUE
           B(INCRJ) = B(INCRJ)/TOPBLK(INCRJ,INCRJ)
 900    CONTINUE
C
C       RETURN FROM THE SOLUTION OF A-TRANSPOSE.X = B

        RETURN
        END


      subroutine lufac(n, ndim, a, ip, ier)
      implicit double precision (a-h,o-z)
      dimension a(ndim,n), ip(n)
      intrinsic abs
*  blas: daxpy, dscal, dswap. idamax

      parameter ( zero = 0.0d+0, one = 1.0d+0 )

*  The subroutine lufac is a very simple code to compute the 
*  LU decomposition (with partial pivoting) of the n by n matrix a.

*  The LU factors are overwritten on a.  The integer array ip
*  reflects the pairwise interchanges performed.  note that ip(k)
*  therefore does not give the index in the original array of
*  the k-th pivot.

*  On exit, the error flag ier is zero when no zero pivots are
*  encountered.  Otherwise, ier is equal to the index of the 
*  step at which a zero pivot occurred.

      ier = 0
      ip(n) = 0

*  Begin loop over columns 1 through n-1.  k is the current
*  column index.

      do 100 k = 1, n-1

*  Find the row index ipiv of the element of largest magnitude in
*  column k.

         ipiv = k-1 + idamax(n-k+1, a(k,k), 1)
         piv = a(ipiv,k) 
         if (piv .eq. zero) then
            ier = k
            return
         endif
         ip(k) = ipiv

*  Perform interchanges if necessary.

         if (ipiv .ne. k) then
            call dswap(n-k+1, a(ipiv,k), ndim, a(k,k), ndim)
         endif

*  Save the (negative) multipliers in the subdiagonal elements of 
*  column k.

         call dscal(n-k, (-one/piv), a(k+1,k), 1)

*  Update the remaining matrix.  Note that a(i,k) now contains
*  the negative multipliers.

         do 50 j = k+1, n
            call daxpy(n-k, a(k,j), a(k+1,k), 1, a(k+1,j), 1)
 50             continue

*  End of loop over columns.

 100             continue
      if (a(n,n).eq.zero) ier = n
      return
      end 

      subroutine lusol (n, ndim, a, ip, b, x) 
      implicit double precision (a-h, o-z)
      dimension a(ndim,n), ip(n), b(n), x(n)

*  blas:  daxpy, dcopy

*  The subroutine lusol is a simple-minded routine to solve a
*  linear system whose LU factors have been computed by lufac.
*  On entry, the matrix a should contain the LU factors, and
*  ip should contain the interchange array constructed by lufac.


*  Copy the right-hand side b into x.

      call dcopy (n, b, 1, x, 1)

*  Forward solution with l (unit lower-triangular factor), which
*  is stored in the strict lower triangle of a.

      do 20 k = 1, n-1
         ipiv = ip(k)
         if (ipiv .ne. k) then
            tem = x(ipiv)
            x(ipiv) = x(k)
            x(k) = tem
         endif
         call daxpy ( n-k, x(k), a(k+1,k), 1, x(k+1), 1 )
 20       continue

*  Backward solution with u (upper-triangular factor), which is stored
*  in the upper triangle of a.

      do 40 kb = n, 1, -1
         x(kb) = x(kb)/a(kb,kb)
         call daxpy(kb-1, (-x(kb)), a(1,kb), 1, x(1), 1)
 40       continue

      return
      end 

      subroutine dcopy ( n, x, incx, y, incy )
      implicit double precision (a-h,o-z)
      integer            n, incx, incy
      dimension          x( * ), y( * )

c  dcopy  performs the operation
c
c     y := x
c
c  nag fortran 77 version of the blas routine dcopy .
c  nag fortran 77 o( n ) basic linear algebra routine.
c
c  -- written on 26-november-1982.
c     sven hammarling, nag central office.

      integer            i     , ix    , iy

      if( n.lt.1 )return

      if( ( incx.eq.incy ).and.( incy.gt.0 ) )then
         do 10, iy = 1, 1 + ( n - 1 )*incy, incy
            y( iy ) = x( iy )
 10             continue
      else
         if( incx.ge.0 )then
            ix = 1
         else
            ix = 1 - ( n - 1 )*incx
         end if
         if( incy.gt.0 )then
            do 20, iy = 1, 1 + ( n - 1 )*incy, incy
               y( iy ) = x( ix )
               ix      = ix + incx
 20                   continue
         else
            iy = 1 - ( n - 1 )*incy
            do 30, i = 1, n
               y( iy ) = x( ix )
               iy      = iy + incy
               ix      = ix + incx
 30                   continue
         end if
      end if

      return

*     end of dcopy .

      end

      subroutine daxpy ( n, alpha, x, incx, y, incy )
      implicit double precision (a-h,o-z)
      integer            n, incx, incy
      dimension   x( * ), y( * )

c  daxpy  performs the operation
c
c     y := alpha*x + y
c
c
c  modified nag fortran 77 version of the blas routine daxpy .
c
c  -- written on 3-september-1982.
c     sven hammarling, nag central office.

      integer            i     , ix    , iy
      parameter        ( zero  = 0.0+0 )

      if( n    .lt.1    )return
      if( alpha.eq.zero )return

      if( ( incx.eq.incy ).and.( incx.gt.0 ) )then
         do 10, ix = 1, 1 + ( n - 1 )*incx, incx
            y( ix ) = alpha*x( ix ) + y( ix )
 10             continue
      else
         if( incy.ge.0 )then
            iy = 1
         else
            iy = 1 - ( n - 1 )*incy
         end if
         if( incx.gt.0 )then
            do 20, ix = 1, 1 + ( n - 1 )*incx, incx
               y( iy ) = alpha*x( ix ) + y( iy )
               iy      = iy + incy
 20                   continue
         else
            ix = 1 - ( n - 1 )*incx
            do 30, i = 1, n
               y( iy ) = alpha*x( ix ) + y( iy )
               ix      = ix + incx
               iy      = iy + incy
 30                   continue
         end if
      end if
      return

*     end of daxpy .

      end

      double precision function ddot  ( n, x, incx, y, incy )
      implicit double precision (a-h,o-z)
      integer           n, incx, incy
      dimension         x( * ), y( * )

c  ddot   returns the value
c
c     ddot   = x'y
c
c
c  modified nag fortran 77 version of the blas routine ddot  .
c
c  -- written on 21-september-1982.
c     sven hammarling, nag central office.

      integer             i     , ix    , iy
      parameter         ( zero  = 0.0d+0 )

      sum = zero
      if( n.ge.1 )then
         if( ( incx.eq.incy ).and.( incx.gt.0 ) )then
            do 10, ix = 1, 1 + ( n - 1 )*incx, incx
               sum = sum + x( ix )*y( ix )
 10                   continue
         else
            if( incy.ge.0 )then
               iy = 1
            else
               iy = 1 - ( n - 1 )*incy
            end if
            if( incx.gt.0 )then
               do 20, ix = 1, 1 + ( n - 1 )*incx, incx
                  sum = sum + x( ix )*y( iy )
                  iy  = iy  + incy
 20                         continue
            else
               ix = 1 - ( n - 1 )*incx
               do 30, i = 1, n
                  sum = sum + x( ix )*y( iy )
                  ix  = ix  + incx
                  iy  = iy  + incy
 30                         continue
            end if
         end if
      end if

      ddot   = sum
      return

*     end of ddot  .

      end

      subroutine dscal ( n, alpha, x, incx )
      implicit double precision (a-h,o-z)
      integer          n, incx
      dimension        x( * )

c  dscal  performs the operation
c
c     x := alpha*x
c
c
c  modified nag fortran 77 version of the blas routine dscal .
c
c  -- written on 26-november-1982.
c     sven hammarling, nag central office.

      integer            ix
      parameter        ( one   = 1.0d+0, zero  = 0.0d+0 )

      if( n.ge.1 )then
         if( alpha.eq.zero )then
            do 10, ix = 1, 1 + ( n - 1 )*incx, incx
               x( ix ) = zero
 10                   continue
         else if( alpha.eq.( -one ) )then
            do 20, ix = 1, 1 + ( n - 1 )*incx, incx
               x( ix ) = -x( ix )
 20                   continue
         else if( alpha.ne.one )then
            do 30, ix = 1, 1 + ( n - 1 )*incx, incx
               x( ix ) = alpha*x( ix )
 30                   continue
         end if
      end if

      return

*     end of dscal .

      end

      subroutine dswap ( n, x, incx, y, incy )
      implicit double precision (a-h,o-z)
      integer      n, incx, incy
      dimension    x( * ), y( * )

c  dswap  performs the operations
c
c     temp := x,   x := y,   y := temp.
c
c
c  modified nag fortran 77 version of the blas routine dswap .
c
c  -- written on 26-november-1982.
c     sven hammarling, nag central office.

      integer            i     , ix    , iy

      if( n.lt.1 )return

      if( ( incx.eq.incy ).and.( incy.gt.0 ) )then
         do 10, iy = 1, 1 + ( n - 1 )*incy, incy
            temp    = x( iy )
            x( iy ) = y( iy )
            y( iy ) = temp
 10             continue
      else
         if( incx.ge.0 )then
            ix = 1
         else
            ix = 1 - ( n - 1 )*incx
         end if
         if( incy.gt.0 )then
            do 20, iy = 1, 1 + ( n - 1 )*incy, incy
               temp    = x( ix )
               x( ix ) = y( iy )
               y( iy ) = temp
               ix      = ix + incx
 20                   continue
         else
            iy = 1 - ( n - 1 )*incy
            do 30, i = 1, n
               temp    = x( ix )
               x( ix ) = y( iy )
               y( iy ) = temp
               iy      = iy + incy
               ix      = ix + incx
 30                   continue
         end if
      end if

      return

*     end of dswap .

      end

      integer function idamax( n, x, incx )
      implicit double precision (a-h,o-z)
      integer         n, incx
      dimension       x( * )

c  idamax returns the smallest value of i such that
c
c     abs( x( i ) ) = max( abs( x( j ) ) )
c                      j
c
c  nag fortran 77 version of the blas routine idamax.
c  nag fortran 77 o( n ) basic linear algebra routine.
c
c  -- written on 31-may-1983.
c     sven hammarling, nag central office.

      intrinsic           abs
      integer             i     , imax  , ix

      if( n.lt.1 )then
         idamax = 0
         return
      end if

      imax = 1
      if( n.gt.1 )then
         xmax = abs( x( 1 ) )
         ix   = 1
         do 10, i = 2, n
            ix = ix + incx
            if( xmax.lt.abs( x( ix ) ) )then
               xmax = abs( x( ix ) )
               imax = i
            end if
 10             continue
      end if

      idamax = imax
      return

*     end of idamax.

      end
      subroutine dload ( n, const, x, incx )
      implicit double precision (a-h,o-z)
      dimension  x(*)
c
c  dload  performs the operation
c
c     x = const*e,   e' = ( 1  1 ... 1 ).
c
c
c  nag fortran 77 o( n ) basic linear algebra routine.
c
c  -- written on 22-september-1983.
c     sven hammarling, nag central office.
c
      parameter        ( zero = 0.0d+0 )

      if( n.lt.1 )return

      if( const.ne.zero )then
         do 10, ix = 1, 1 + ( n - 1 )*incx, incx
            x( ix ) = const
 10             continue
      else
         do 20, ix = 1, 1 + ( n - 1 )*incx, incx
            x( ix ) = zero
 20             continue
      end if

      return

*     end of dload .
      end



      subroutine maxpy ( nrow, ncol, alpha, xmat, nrowy, ymat )
      implicit double precision (a-h,o-z)
      dimension xmat(nrow, ncol), ymat(nrowy, ncol)

*  Subroutine maxpy takes as input the scalar alpha and two matrices,
*  xmat and ymat.  xmat has declared row dimension nrow, and
*  ymat has declared row dimension nrowy, but both are
*  conceptually nrow by ncol.  
*  On output, (new ymat) is alpha*xmat+ (old ymat), by analogy 
*  with the vector blas routine saxpy.

      do 100 j = 1, ncol
      do 100 i = 1, nrow
         ymat(i,j) = ymat(i,j) + alpha*xmat(i,j)
 100      continue
      return
      end


      subroutine matcop( nrow1, nrow2, nrow, ncol, xmat1, xmat2 )
      implicit double precision (a-h,o-z)
      dimension xmat1(nrow1, ncol), xmat2(nrow2, ncol)

*  Given 2 matrices xmat1 and xmat2, where xmat1 has declared 
*  row dimension nrow1, xmat2 has declared row dimension nrow2,
*  and both have column dimension ncol, the routine matcop copies
*  rows 1 through nrow, and columns 1 through ncol from xmat1 into
*  xmat2.


      if (nrow .le. 0 .or. ncol .le. 0) return
      do 100 j = 1, ncol
      do 100 i = 1, nrow
           xmat2(i,j) = xmat1(i,j)
 100        continue
      return
      end

      subroutine mtload( nrow, ncol, const, nrowx, xmat )
      implicit double precision (a-h,o-z)
      dimension xmat(nrowx, ncol)

*  mtload sets elements 1 through nrow, 1 through ncol, of the
*  matrix xmat (whose declared row dimension is nrowx) to the 
*  scalar value const.  

      if (nrow .le. 0 .or. ncol .le. 0)  return
      do 100 j = 1, ncol
      do 100 i = 1, nrow
         xmat(i,j) = const
 100      continue
      return
      end

      subroutine mtload1( nrow, ncol, xx, nrowx, xmat )
      implicit double precision (a-h,o-z)
      dimension xmat(nrowx, ncol), xx(*)

*  mtload sets elements 1 through nrow, 1 through ncol, of the
*  matrix xmat (whose declared row dimension is nrowx) to the 
*  scalar value const.  

      if (nrow .le. 0 .or. ncol .le. 0)  return
      do 100 j = 1, ncol
      do 100 i = 1, nrow
      if (i.eq.1) then
        xmat(i,j) = 2.D0 * xx(i-1)-1.D0
      else if (i.eq.2) then
        xmat(i,j) = 2.D0
      else
        xmat(i,j) = 0.D0
      end if
 100   continue
      return
      end


      subroutine mssq  ( nrow, ncol, xmat, scale, sumsq )
      implicit double precision (a-h,o-z)
      dimension   xmat(nrow, *)

*  Given the nrow by ncol matrix xmat, mssq returns values
*  scale and sumsq such that
*     (scale**2) * sumsq = sum of squares of xmat(i,j),
*  where  scale = max  abs(xmat(i,j)).

*  mssq is a stripped-down matrix version of the blas routine sssq.

      intrinsic    abs
      parameter   ( one   = 1.0d+0, zero  = 0.0d+0 )

      scale = zero
      sumsq = one
      if( nrow.ge.1 .and. ncol .ge. 1) then
         do 10 i = 1, nrow
         do 10 j = 1, ncol 
            if( xmat(i,j) .ne. zero )then
               absxij = abs(xmat(i,j))
               if( scale .lt. absxij ) then
                  sumsq = one + sumsq* (scale/absxij)**2
                  scale = absxij
               else
                  sumsq = sumsq + (absxij/scale)**2
               end if
            end if
 10             continue
      end if
      return
      end

      subroutine dssq  ( n, x, incx, scale, sumsq )
      implicit double precision (a-h,o-z)
      integer            n, incx
      dimension   x( * )

*  Given the n-vector x, dssq returns values scale and sumsq such that
*     (scale**2) * sumsq = sum of squares of x(i),
*  where  scale = max  abs(x(i)).

      intrinsic          abs
      parameter        ( one   = 1.0d+0, zero  = 0.0d+0 )

      scale = zero
      sumsq = one
      if( n.ge.1 )then
         do 10, ix = 1, 1 + ( n - 1 )*incx, incx
            if( x( ix ).ne.zero )then
               absxi = abs( x( ix ) )
               if( scale.lt.absxi )then
                  sumsq = one   + sumsq*( scale/absxi )**2
                  scale = absxi
               else
                  sumsq = sumsq + ( absxi/scale )**2
               end if
            end if
 10             continue
      end if

      return

*     end of dssq  .

      end


