#define VERBOSE
c-----------------------------------------------------------------------
c
c     module mud2crf90
c
c     . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
c     .                                                             .
c     .                  copyright (c) 1999                         .
c     .                                                             .
c     .       UNIVERSITY CORPORATION for ATMOSPHERIC RESEARCH       .
c     .                                                             .
c     .                      all rights reserved                    .
c     .                                                             .
c     .                                                             .
c     .                      MUDPACK version 5.0                    .
c     .                                                             .
c     . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
c
c
c     Author for this fortran 90 module
c
c     In-Sun Song
c     Laboratory for Mesoscale Dynamics
c     Department of Atmospheric Sciences, Yonsei University,
c     Seoul, Korea.
c
c     January 19, 2004
c     This module is first written.
c     Obsolete features such as equivanlence, common blocks, and 
c     external function declaration that has been used in Fortran 77
c     are removed. To facilitate the implement of this routine into
c     users' own programs, the whole of routines are included within
c     a Fortran 90/95 module. In addition, each module procedure
c     is slightly modified as to have explicit interface to calling
c     environments.
c
c     Author for original code and specialist
c
c     John C. Adams (National Center for Atmospheric Research)
c     email: johnad@ucar.edu, phone: 303-497-1213
c
c     For MUDPACK 5.0 information, visit the website:
c     (http://www.scd.ucar.edu/css/software/mudpack)
c
c     Purpose (see mud2cr.d for details)
c
c     mud2cr attempts to produce a second order finite difference
c     approximation to the two dimensional nonseparable elliptic
c     partial differential equation with cross derivative
c
c       cxx(x,y)*pxx + cxy(x,y)*pxy + cyy(x,y)*pyy +
c
c       cx(x,y)*px + cy(x,y)*py + ce(x,y)*pe(x,y) = r(x,y)
c
c     Documentation and test files
c
c     see the documentation file "mud2cr.d" for a complete discussion
c     of how to use subroutine mud2cr.  file "tmud2cr.f" is a test/driver
c     sample program illustrating use of mud2cr
c
c     Required MUDPACK files
c     
c     mudcom.F
c
c     Coefficients and boundary conditions of elliptic equation
c     are determined in this module (at the end of this module).
c
c     Test system
c
c     OSF5.1 (HP Tru64)
c     Cygwin with g95 (PC)
c
c-----------------------------------------------------------------------
c
      module mud2crf90
c
      use mudcom    
c
      type :: iptr
        integer, pointer :: p
      end type
      type :: fptr
        real(r8), pointer :: p
      end type
c
c     In original tmud2cr.f, iprm is defined as an integer array with
c     16 cells. However, this definition is inconsistent with iparm 
c     (with 17 cells) in mud2cr.f. In former programs, each file is
c     compiled separately, and thus such inconsistency seems not to
c     appear. To remove this, here, iprm is defined as an integer
c     array with 17 cells. Note that this modification will not affect
c     results.
c
      type(iptr), dimension(17) :: iprm
      type(fptr), dimension(6)  :: fprm
c
c     Variables that were in common block itmud2cr in tmud2cr.f
c
      integer, target :: intl0 ,nxa0  ,nxb0  ,nyc0  ,nyd0  
      integer, target :: ixp0  ,jyq0  ,iex0  ,jey0
      integer, target :: nx0   ,ny0   ,iguess0,maxcy0,method0,nwork0
      integer, target :: lwrkqd0,itero0
c
c     Variables that were in common block ftmud2cr in tmud2cr.f
c
      real(r8),  target :: xa0   ,xb0   ,yc0   ,yd0   ,tolmax0,relmax0
c
c     Variables that were in common block mud2crc in mud2cr.f
c
      integer, dimension(50) :: kpbgn ,kcbgn ,ktxbgn,ktybgn
      integer, dimension(50) :: nxk   ,nyk
      integer :: isx,jsy
c
c     Variables for coefficients and boundary conditions
c
      real(r8), allocatable, dimension(:)   :: xusr  ,yusr
      real(r8), allocatable, dimension(:,:) :: cxxusr,cxyusr,cyyusr   
      real(r8), allocatable, dimension(:,:) :: cxusr ,cyusr ,ceusr   
c
      private :: iprm  ,fprm
c
      private :: intl0 ,nxa0  ,nxb0  ,nyc0  ,nyd0
      private :: ixp0  ,jyq0  ,iex0  ,jey0  
      private :: nx0   ,ny0   ,iguess0,maxcy0,method0,nwork0
      private :: lwrkqd0,itero0
c
      private :: xa0   ,xb0   ,yc0   ,yd0   ,tolmax0,relmax0
c
      private :: kpbgn ,kcbgn ,ktxbgn,ktybgn,nxk,nyk
      private :: isx   ,jsy
c
      public :: xusr  ,yusr
      public :: cxxusr,cxyusr,cyyusr
      public :: cxusr ,cyusr ,ceusr
c
      contains
c
c-----------------------------------------------------------------------
c
      subroutine sol2cr(iixp  ,jjyq  ,iiex  ,jjey  ,nnx   ,nny   ,
     +                  nnxa  ,nnxb  ,nnyc  ,nnyd  ,iigues,mmaxcy,
     +                  mmethd,mgopt ,xval  ,yval  ,frc   ,var   ,
     +                  niter )

      implicit none
c
      integer, intent(in)    :: iixp  ,jjyq  ,iiex  ,jjey  ,nnx   ,nny
      integer, intent(inout) :: nnxa  ,nnxb  ,nnyc  ,nnyd
      integer, intent(inout) :: iigues,mmaxcy,mmethd
      integer, dimension(4), intent(inout) :: mgopt
c
      real(r8), dimension(nnx),     intent(inout) :: xval
      real(r8), dimension(nny),     intent(inout) :: yval
      real(r8), dimension(nnx,nny), intent(inout) :: frc
      real(r8), dimension(nnx,nny), intent(inout) :: var
c
      integer, intent(in) :: niter
c
      real(r8), allocatable, dimension(:) :: work
c
      integer :: i     ,j     ,k     ,ierror
      logical :: check1,check2
c
c     declare coefficient and boundary condition input routines external
c
      iprm( 1)%p => intl0
      iprm( 2)%p => nxa0
      iprm( 3)%p => nxb0
      iprm( 4)%p => nyc0
      iprm( 5)%p => nyd0
      iprm( 6)%p => ixp0
      iprm( 7)%p => jyq0
      iprm( 8)%p => iex0
      iprm( 9)%p => jey0
      iprm(10)%p => nx0
      iprm(11)%p => ny0
      iprm(12)%p => iguess0
      iprm(13)%p => maxcy0
      iprm(14)%p => method0
      iprm(15)%p => nwork0
      iprm(16)%p => lwrkqd0
      iprm(17)%p => itero0
c
      fprm( 1)%p => xa0
      fprm( 2)%p => xb0
      fprm( 3)%p => yc0
      fprm( 4)%p => yd0
      fprm( 5)%p => tolmax0
      fprm( 6)%p => relmax0
c
c     set input integer arguments
c
      intl0 = 0
c
c     set grid sizes from parameter statements
c
      ixp0  = iixp
      jyq0  = jjyq
      iex0  = iiex
      jey0  = jjey
      nx0   = nnx
      ny0   = nny
c
c     check grid size consistency
c
      check1 = nx0 /= ixp0*(2**(iex0-1))+1
      check2 = ny0 /= jyq0*(2**(jey0-1))+1
      if ( check1 .or. check2 ) then
        write(6,'(a)') '(sol2cr): Check nnx and nny.'
        stop
      end if
c
c     set bourndary condition flags
c
      nxa0 = nnxa
      nxb0 = nnxb
      nyc0 = nnyc
      nyd0 = nnyd
c
c     set performance flags
c
      maxcy0 = mmaxcy
      iguess0 = iigues
      method0 = mmethd
c
      nwork0 = (7*(nnx+2)*(nny+2)+44*nnx*nny)/3
      allocate(work(1:nwork0))
c
      xa0 = xval(1)
      xb0 = xval(nx0)
      yc0 = yval(1)
      yd0 = yval(ny0)
c
      tolmax0 = 0.
      relmax0 = 0.
c
#if ( defined VERBOSE )
      write(6,'( )')
      write(6,100) '------------------------------------------------'
      write(6,100) '                    SOL2CR'
      write(6,100) '------------------------------------------------'
      write(6,101) ' INITIALIZATION ARGUMENT   : ',iprm(1)%p
      write(6,101) ' NUMBER OF X-GRID (XA-XB)  : ',iprm(10)%p
      write(6,101) ' NUMBER OF Y-GRID (YC-YD)  : ',iprm(11)%p
      write(6,101) ' BOUNDARY CONDITION (XA)   : ',iprm(2)%p
      write(6,101) ' BOUNDARY CONDITION (XB)   : ',iprm(3)%p
      write(6,101) ' BOUNDARY CONDITION (YC)   : ',iprm(4)%p
      write(6,101) ' BOUNDARY CONDITION (YD)   : ',iprm(5)%p
      write(6,101) ' INITIAL GUESS FIELD       : ',iprm(12)%p
      write(6,101) ' MULTIGRID MAXIMUM CYCLE   : ',iprm(13)%p
      write(6,101) ' RELAXATION METHOD         : ',iprm(14)%p
      write(6,101) ' WORK SPACE ESTIMATE       : ',iprm(15)%p
      write(6,102) ' XA                        : ',fprm(1)%p
      write(6,102) ' XB                        : ',fprm(2)%p
      write(6,102) ' YC                        : ',fprm(3)%p
      write(6,102) ' YD                        : ',fprm(4)%p
      write(6,102) ' TOLERANCE (ERROR CONTROL) : ',fprm(5)%p
      write(6,101) ' KCYCLE (MULTIGRID OPTION) : ',mgopt(1)
      write(6,101) ' IPRER  (MULTIGRID OPTION) : ',mgopt(2)
      write(6,101) ' IPOST  (MULTIGRID OPTION) : ',mgopt(3)
      write(6,101) ' INTPOL (MULTIGRID OPTION) : ',mgopt(4)
#endif
c
c     Initialization call
c
#if ( defined VERBOSE )
      write(6,*)
      write(6,101) ' INITIALIZATION ARGUMENT   : ',iprm(1)%p
#endif
      call mud2cr(iprm,fprm,work,frc,var,mgopt,ierror)
c
#if ( defined VERBOSE )
      write(6,101) ' ERROR CODE                : ',ierror
      write(6,101) ' MINIMUM WORK SPACE        : ',iprm(16)%p
#endif
      if (ierror > 0) call exit(0)
c
c     ATTEMPT SOLUTION
c
      intl0 = 1
#if ( defined VERBOSE )
      write(6,*)
      write(6,101) ' INITIALIZATION ARGUMENT   : ',iprm(1)%p
      write(6,101) ' RELAXATION METHOD         : ',iprm(14)%p
      write(6,101) ' INITIAL GUESS FIELD       : ',iprm(12)%p
#endif
      call mud2cr(iprm,fprm,work,frc,var,mgopt,ierror)
c
#if ( defined VERBOSE )
      write(6,101) ' ERROR CODE                : ',ierror
#endif
c
      iguess0 = 1
c
      do i=1,niter
c
#if ( defined VERBOSE )
      write(6,*)
      write(6,101) ' NUMBER OF ITERATION       : ',i
      write(6,101) ' INITIALIZATION ARGUMENT   : ',iprm(1)%p
      write(6,101) ' RELAXATION METHOD         : ',iprm(14)%p
      write(6,101) ' INITIAL GUESS FIELD       : ',iprm(12)%p
#endif
      call mud2cr(iprm,fprm,work,frc,var,mgopt,ierror)
c
#if ( defined VERBOSE )
      WRITE(6,101) ' ERROR CODE                : ',ierror
#endif
c
      end do
c
      deallocate(work)
c
  100 format(A)
  101 format(A,I10)
  102 format(A,E17.10)
c      
      return 
      end subroutine sol2cr
c      
c-----------------------------------------------------------------------
c
      subroutine mud2cr(iparm,fparm,work,rhs,phi,mgopt,ierror)
c
      implicit none
      type(iptr) iparm
      type(fptr) fparm
      integer mgopt,ierror
      integer int,iw,k,kb,nx,ny,ic,itx,ity
      real(r8) phmax
      dimension iparm(17),fparm(6),mgopt(4)
      real(r8) work(*),phi(*),rhs(*)
c
      data int / 0 /
      save int
c
      ierror = 1
      intl = iparm(1)%p    ! set and check intl on all calls
      if (intl*(intl-1).ne.0) return
      if (int.eq.0) then
      int = 1
      if (intl.ne.0) return  ! very first call is not intl=0
      end if
      ierror = 0
c
c     set  arguments internally
c     these will not be rechecked if intl=1!
c
      nxa = iparm(2)%p
      nxb = iparm(3)%p
      nyc = iparm(4)%p
      nyd = iparm(5)%p
      ixp = iparm(6)%p
      jyq = iparm(7)%p
      iex = iparm(8)%p
      jey = iparm(9)%p
      ngrid = max0(iex,jey)
      nfx = iparm(10)%p
      nfy = iparm(11)%p
      iguess = iparm(12)%p
      maxcy = iparm(13)%p
      method = iparm(14)%p
      nwork = iparm(15)%p
      kcycle = mgopt(1)
      if (kcycle .eq. 0) then
c       set defaults
        kcycle = 2
        iprer = 2
        ipost = 1
        intpol = 3
      else
        iprer = mgopt(2)
        ipost = mgopt(3)
        intpol = mgopt(4)
      end if
      xa = fparm(1)%p
      xb = fparm(2)%p
      yc = fparm(3)%p
      yd = fparm(4)%p
      tolmax = fparm(5)%p
      if (intl .eq. 0) then  ! intialization call
c
c     check input arguments
c
        ierror = 2   ! check boundary condition flags
        if (max0(nxa,nxb,nyc,nyd).gt.2) return
        if (min0(nxa,nxb,nyc,nyd).lt.0) return
        if (nxa.eq.0.and.nxb.ne.0) return
        if (nxa.ne.0.and.nxb.eq.0) return
        if (nyc.eq.0.and.nyd.ne.0) return
        if (nyc.ne.0.and.nyd.eq.0) return
        ierror = 3   ! check grid sizes
        if (ixp.lt.2) return
        if (jyq.lt.2) return
        ierror = 4
        ngrid = max0(iex,jey)
        if (iex.lt.1) return
        if (jey.lt.1) return
        if (ngrid.gt.50) return
        ierror = 5
        if (nfx.ne.ixp*2**(iex-1)+1) return
        if (nfy.ne.jyq*2**(jey-1)+1) return
        ierror = 6
        if (iguess*(iguess-1).ne.0) return
        ierror = 7
        if (maxcy.lt.1) return
        ierror = 8
        if (method.lt.0 .or. method.gt.3) return
        ierror = 9
c       compute and test minimum work space
        isx = 0
        if (method.eq.1 .or. method.eq.3) then
          if (nxa.ne.0) isx = 3
          if (nxa.eq.0) isx = 5
        end if
        jsy = 0
        if (method.eq.2 .or. method.eq.3) then
          if (nyc.ne.0) jsy = 3
          if (nyc.eq.0) jsy = 5
        end if
        kps = 1
        do k=1,ngrid
c       set subgrid sizes
        nxk(k) = ixp*2**(max0(k+iex-ngrid,1)-1)+1
        nyk(k) = jyq*2**(max0(k+jey-ngrid,1)-1)+1
        nx = nxk(k)
        ny = nyk(k)
        kps = kps+(nx+2)*(ny+2)+nx*ny*(10+isx+jsy)
      end do
      iparm(16)%p = kps+(nfx+2)*(nfy+2)   ! exact minimum work space
      lwork = iparm(16)%p
      if (lwork .gt. nwork) return
      ierror = 10   ! check solution region
      if (xb.le.xa .or. yd.le.yc) return
      ierror = 11
      if (tolmax .lt. 0.0) return
      ierror = 12   ! multigrid parameters
      if (kcycle.lt.0) return
      if (min0(iprer,ipost).lt.1) return
      if ((intpol-1)*(intpol-3).ne.0) return
      if (max0(kcycle,iprer,ipost).gt.2) then
        ierror = -5   ! inefficient multigrid cycling
      end if
      if (ierror .gt. 0) ierror = 0   ! no fatal errors
c
c     set work space pointers and discretize pde at each grid level
c
      iw = 1
      do kb=1,ngrid
        k = ngrid-kb+1
        nx = nxk(k)
        ny = nyk(k)
        kpbgn(k) = iw
        kcbgn(k) = kpbgn(k)+(nx+2)*(ny+2)
        ktxbgn(k) = kcbgn(k)+10*nx*ny
        ktybgn(k) = ktxbgn(k)+isx*nx*ny
        iw = ktybgn(k)+jsy*nx*ny
        ic = kcbgn(k)
        itx = ktxbgn(k)
        ity = ktybgn(k)
        klevel = k
        call dismd2cr(nx,ny,work(ic),work(itx),work(ity),
     +                work,ierror)
        end do
      return
      end if   ! end of intl=0 initialization call block
      nx = nfx
      ny = nfy
      call mud2cr1(nx,ny,rhs,phi,work)
      iparm(17)%p = itero
      if (tolmax.gt.0.0) then   ! check for convergence
        fparm(6)%p = relmax
        if (relmax.gt.tolmax) ierror = -1   ! flag convergenc failure
      end if
      return

      end subroutine mud2cr

c-----------------------------------------------------------------------

      subroutine mud2cr1(nx,ny,rhsf,phif,wk)
      implicit none
      integer nx,ny
      real(r8) phif(nx,ny),rhsf(nx,ny),wk(*)
      real(r8) phmax
      integer k,kb,ip,ic,ir,ipc,irc,icc
      integer ncx,ncy,jj,ij,i,j,iter
c
      nx = nxk(ngrid)
      ny = nyk(ngrid)
      ip = kpbgn(ngrid)
      ic = kcbgn(ngrid)
      ir = ic+9*nx*ny
c
c     set phif,rhsf in wk and adjust right hand side
c
      call swk2(nx,ny,phif,rhsf,wk(ip),wk(ir))
      if (iguess.eq.0) then
c
c     no initial guess at finest grid level!
c
        do kb=2,ngrid
          k = ngrid-kb+1
          nx = nxk(k+1)
          ny = nyk(k+1)
          ip = kpbgn(k+1)
          ir = kcbgn(k+1)+9*nx*ny
          ncx = nxk(k)
          ncy = nyk(k)
          ipc = kpbgn(k)
          icc = kcbgn(k)
          irc = icc+9*ncx*ncy
c
c     transfer down to all grid levels
c
          call trsfc2(nx,ny,wk(ip),wk(ir),ncx,ncy,
     +                wk(ipc),wk(irc))
        end do
c
c     adjust right hand side at all grid levels in case
c     rhs or specified b.c. in phi or gbdy changed
c
        do k=1,ngrid
          nx = nxk(k)
          ny = nyk(k)
          ip = kpbgn(k)
          ic = kcbgn(k)
          call adjmd2cr(nx,ny,wk(ip),wk(ic))
        end do
c
c     execute one full multigrid cycle
c
        do k=1,ngrid-1
          kcur = k
          call kcymd2cr(wk)
          nx = nxk(k+1)
          ny = nyk(k+1)
          ip = kpbgn(k+1)
          ipc = kpbgn(k)
          ncx = nxk(k)
          ncy = nyk(k)
c
c     lift or prolong approximation from k to k+1
c
          call prolon2(ncx,ncy,wk(ipc),nx,ny,wk(ip),nxa,nxb,
     +                 nyc,nyd,intpol)
        end do
      else
c
c     adjust rhs at finest grid level only
c
        nx = nxk(ngrid)
        ny = nyk(ngrid)
        ip = kpbgn(ngrid)
        ic = kcbgn(ngrid)
        call adjmd2cr(nx,ny,wk(ip),wk(ic))
      end if
c
c     execute maxcy more multigrid k cycles from finest level
c
      kcur = ngrid
      do iter=1,maxcy
        itero = iter
        call kcymd2cr(wk)
        if (tolmax.gt.0.0) then
c
c      error control
c
          relmax = 0.0
          phmax = 0.0
          do j=1,nfy
            jj = j*(nfx+2)
            do i=1,nfx
              ij = jj+i+1
              phmax = max(phmax,abs(wk(ij)))
              relmax = max(relmax,abs(wk(ij)-phif(i,j)))
              phif(i,j) = wk(ij)
            end do
          end do
c
c     set maximum relative difference and check for convergence
c
          if (phmax.gt.0.0) relmax = relmax/phmax
          if (relmax.le.tolmax) return
        end if
      end do
c
c     set final interate after maxcy cycles in phif
c
      do j=1,nfy
         jj = j*(nfx+2)
         do i=1,nfx
         ij = jj+i+1
         phif(i,j) = wk(ij)
        end do
      end do
      return
      end subroutine mud2cr1
c
c----------------------------------------------------------------------
c
      subroutine kcymd2cr(wk)
c
c     execute multigrid k cycle from kcur grid level
c     kcycle=1 for v cycles, kcycle=2 for w cycles
c
      implicit none
      real(r8) wk(*)
      integer nx,ny,ip,ic,ipc,irc,itx,ity,ncx,ncy,l,nrel
      integer kount(50)
      klevel = kcur
      nx = nxk(klevel)
      ny = nyk(klevel)
      ip = kpbgn(klevel)
      ic = kcbgn(klevel)
      itx = ktxbgn(klevel)
      ity = ktybgn(klevel)
c
c     prerelax at current finest grid level
c
      do l=1,iprer
        call relmd2cr(nx,ny,wk(ip),wk(ic),wk(itx),wk(ity),wk(kps))
      end do
      if (kcur .eq. 1) go to 5
c
c     restrict residual to kcur-1 level
c
      ipc = kpbgn(klevel-1)
      ncx = nxk(klevel-1)
      ncy = nyk(klevel-1)
      irc = kcbgn(klevel-1)+9*ncx*ncy
      call resmd2cr(nx,ny,wk(ip),ncx,ncy,wk(ipc),wk(irc),wk(ic),wk(kps))
c
c    set counter for grid levels to zero
c
      do l = 1,kcur
        kount(l) = 0
      end do
c
c    set new grid level and continue k-cycling
c
      klevel = kcur-1
      nrel = iprer
c
c   kcycle control point
c
   10 continue
c
c      post relax when kcur revisited
c
      if (klevel .eq. kcur) go to 5
c
c   count hit at current level
c
      kount(klevel) = kount(klevel)+1
c
c   relax at current level
c
      nx = nxk(klevel)
      ny = nyk(klevel)
      ip = kpbgn(klevel)
      ic = kcbgn(klevel)
      itx = ktxbgn(klevel)
      ity = ktybgn(klevel)
      do l=1,nrel
        call relmd2cr(nx,ny,wk(ip),wk(ic),wk(itx),wk(ity),wk(kps))
      end do
      if (kount(klevel) .eq. kcycle+1) then
c
c     kcycle complete at klevel
c
        ipc = ip
        ip = kpbgn(klevel+1)
        ncx = nxk(klevel)
        ncy = nyk(klevel)
        nx = nxk(klevel+1)
        ny = nyk(klevel+1)
c
c    inject correction to finer grid
c
        call cor2(nx,ny,wk(ip),ncx,ncy,wk(ipc),nxa,nxb,nyc,nyd,
     +            intpol,wk(kps))
c
c    reset counter to zero
c
        kount(klevel) = 0
c
c     ascend to next higher level and set to postrelax there
c
        klevel = klevel+1
        nrel = ipost
        go to 10
      else
        if (klevel .gt. 1) then
c
c    kcycle not complete so descend unless at coarsest grid
c
         ipc = kpbgn(klevel-1)
         ncx = nxk(klevel-1)
         ncy = nyk(klevel-1)
         irc = kcbgn(klevel-1)+9*ncx*ncy
         call resmd2cr(nx,ny,wk(ip),ncx,ncy,wk(ipc),wk(irc),wk(ic),
     +                wk(kps))
c
c     prerelax at next coarser level
c
          klevel = klevel-1
          nrel = iprer
          go to 10
        else
c
c    postrelax at coarsest level
c
        do l=1,ipost
          call relmd2cr(nx,ny,wk(ip),wk(ic),wk(itx),wk(ity),wk(kps))
        end do
        ipc = ip
        ip = kpbgn(2)
        ncx = nxk(1)
        ncy = nyk(1)
        nx = nxk(2)
        ny = nyk(2)
c
c    inject correction to level 2
c
        call cor2(nx,ny,wk(ip),ncx,ncy,wk(ipc),nxa,nxb,nyc,nyd,
     +            intpol,wk(kps))
c
c     set to postrelax at level 2
c
        nrel = ipost
        klevel = 2
        go to 10
      end if
      end if
    5 continue
c
c     post relax at current finest grid level
c
      nx = nxk(kcur)
      ny = nyk(kcur)
      ip = kpbgn(kcur)
      ic = kcbgn(kcur)
      itx = ktxbgn(kcur)
      ity = ktybgn(kcur)
      do l=1,ipost
        call relmd2cr(nx,ny,wk(ip),wk(ic),wk(itx),wk(ity),wk(kps))
      end do
      return
      end subroutine kcymd2cr
c
c-----------------------------------------------------------------------
c
      subroutine dismd2cr(nx,ny,cf,tx,ty,wk,ier)
c
c     discretize elliptic pde for mud2cr, set nonfatal errors
c
      implicit none
      integer nx,ny,i,j,kbdy,l,im1,jm1,ier,jc
      real(r8) cf(nx,ny,10),tx(nx,ny,*),ty(ny,nx,*)
      real(r8) wk(*),dlx,dlx2,dlxx,dly,dly2,dlyy,cmin,alfmax,cemax
      real(r8) x,y,cxx,cxy,cyy,cx,cy,ce,c1,c2,c3,c4,c5
      real(r8) c6,c7,c8,c9
      real(r8) alfaa,alfab,alfac,alfad,betaa,betab,betac,betad,det
      real(r8) gamaa,gamab,gamac,gamad,dxoy,dyox,dlxy,dlxy2,dlxy4
      real(r8) alfim1,alfi,alfip1,betim1,beti,betip1,gamim1,gami,gamip1
      real(r8) alfjm1,alfj,alfjp1,betjm1,betj,betjp1,gamjm1,gamj,gamjp1
      real(r8) gammax,gbdim1,gbdi,gbdip1,gbdj,gbdjm1,gbdjp1
      real(r8) gbdya,gbdyb,gbdyc,gbdyd
c
      dlx = (xb-xa)/(nx-1)
      dlx2 = dlx+dlx
      dlxx = dlx*dlx
      dly = (yd-yc)/(ny-1)
      dly2 = dly+dly
      dlyy = dly*dly
      dlxy = dlx*dly
      dlxy2 = dlxy+dlxy
      dlxy4 = dlxy2+dlxy2
      dyox = dly/dlx
      dxoy = dlx/dly
      cmin = 1.0
      alfmax = 0.0
      cemax = 0.0
c
c     compute discretization coefficients on interior and
c     nonspecified boundaries
c
      do j=1,ny
        y = yc+(j-1)*dly
        do i=1,nx
          x = xa+(i-1)*dlx
          call coef(x,y,cxx,cxy,cyy,cx,cy,ce)
          cmin = min(cmin,cxx*cyy)
          cemax = max(abs(ce),cemax)
c
c     flag hyperbolic pde
c
          if (klevel.eq.ngrid) then
            if (abs(cx)*dlx.gt.2*abs(cxx) .or.
     +          abs(cy)*dly.gt.2*abs(cyy)) then
              ier = -4
            end if
          end if
          cxx = max(cxx,abs(cx)*dlx*0.5)
          cyy = max(cyy,abs(cy)*dly*0.5)
          c1=cxx/dlxx+cx/dlx2
          c2=cxy/dlxy4
          c3=cyy/dlyy+cy/dly2
          c4=-c2
          c5=cxx/dlxx-cx/dlx2
          c6=c2
          c7=cyy/dlyy-cy/dly2
          c8=-c2
          c9=ce-(c1+c3+c5+c7)
          cf(i,j,1)=c1
          cf(i,j,2)=c2
          cf(i,j,3)=c3
          cf(i,j,4)=c4
          cf(i,j,5)=c5
          cf(i,j,6)=c6
          cf(i,j,7)=c7
          cf(i,j,8)=c8
          cf(i,j,9)=c9
        end do
      end do
c
c     adjust at cornors virtual-virtual points not used
c
      call coef(xa,yc,cxx,cxy,cyy,cx,cy,ce)
      c4=-cxy/dlxy2
      cf(1,1,1)=cf(1,1,1)-c4
      cf(1,1,2)=0.0
      cf(1,1,3)=cf(1,1,3)-c4
      cf(1,1,4)=c4
      cf(1,1,5)=cf(1,1,5)-c4
      cf(1,1,6)=0.0
      cf(1,1,7)=cf(1,1,7)-c4
      cf(1,1,8)=c4
      cf(1,1,9)=cf(1,1,9)+2.0*c4
      call coef(xb,yc,cxx,cxy,cyy,cx,cy,ce)
      c2=cxy/dlxy2
      cf(nx,1,1)=cf(nx,1,1)-c2
      cf(nx,1,2)=c2
      cf(nx,1,3)=cf(nx,1,3)-c2
      cf(nx,1,4)=0.0
      cf(nx,1,5)=cf(nx,1,5)-c2
      cf(nx,1,6)=c2
      cf(nx,1,7)=cf(nx,1,7)-c2
      cf(nx,1,8)=0.0
      cf(nx,1,9)=cf(nx,1,9)+2.0*c2
      call coef(xa,yd,cxx,cxy,cyy,cx,cy,ce)
      c2=cxy/dlxy2
      cf(1,ny,1)=cf(1,ny,1)-c2
      cf(1,ny,2)=c2
      cf(1,ny,3)=cf(1,ny,3)-c2
      cf(1,ny,4)=0.0
      cf(1,ny,5)=cf(1,ny,5)-c2
      cf(1,ny,6)=c2
      cf(1,ny,7)=cf(1,ny,7)-c2
      cf(1,ny,8)=0.0
      cf(1,ny,9)=cf(1,ny,9)+2.0*c2
      call coef(xb,yd,cxx,cxy,cyy,cx,cy,ce)
      i = nx
      c4=-cxy/dlxy2
      cf(nx,ny,1)=cf(nx,ny,1)-c4
      cf(nx,ny,2)=0.0
      cf(nx,ny,3)=cf(nx,ny,3)-c4
      cf(nx,ny,4)=c4
      cf(nx,ny,5)=cf(nx,ny,5)-c4
      cf(nx,ny,6)=0.0
      cf(nx,ny,7)=cf(nx,ny,7)-c4
      cf(nx,ny,8)=c4
      cf(nx,ny,9)=cf(nx,ny,9)+2.0*c4
c
c     adjust discretization for mixed derivative b.c.
c
      gammax = 0.0
      if (nxa.eq.2) then
        kbdy = 1
        x = xa
        i = 1
        call bndyc(kbdy,yc,alfjm1,betjm1,gamjm1,gbdjm1)
        call bndyc(kbdy,yc+dly,alfj,betj,gamj,gbdj)
        gammax = max(abs(gamjm1),abs(gamj),gammax)
        do j=2,ny-1
          jc = j
          y = yc+j*dly
          call bndyc(kbdy,y,alfjp1,betjp1,gamjp1,gbdjp1)
          gammax = max(abs(gamjp1),gammax)
c
c     check for illegal tangential derivative b.c.
c
          if (alfjm1*alfj*alfjp1.eq.0.0) then
            ier = 13
            return
          end if
          c4=cf(1,jc,4)
          c5=cf(1,jc,5)
          c6=cf(1,jc,6)
          cf(1,jc,1)=cf(1,jc,1)+c5
          cf(1,jc,2)=cf(1,jc,2)+c4
          cf(1,jc,3)=cf(1,jc,3)+c6*(-betjm1/alfjm1*dxoy)+
     +    c5*(betj/alfj*dxoy)+c4*(3.*betjp1/alfjp1*dxoy+dlx2*gamjp1/
     +    alfjp1)
          cf(1,jc,4)=0.0
          cf(1,jc,5)=0.0
          cf(1,jc,6)=0.0
          cf(1,jc,7)=cf(1,jc,7)+c6*(-3.*betjm1/alfjm1*dxoy+dlx2*gamjm1/
     +    alfjm1)+c5*(-betj/alfj*dxoy)+c4*(betjp1/alfjp1*dxoy)
          cf(1,jc,8)=cf(1,jc,8)+c6
          cf(1,jc,9)=cf(1,jc,9)+c6*(4.*betjm1/alfjm1*dxoy)+
     +    c5*(dlx2*gamj/alfj)+c4*(-4.*betjp1/alfjp1*dxoy)
          alfjm1 = alfj
          betjm1 = betj
          gamjm1 = gamj
          gbdjm1 = gbdj
          alfj = alfjp1
          betj = betjp1
          gamj = gamjp1
          gbdj = gbdjp1
        end do
      end if

      if (nxb.eq.2) then
        kbdy = 2
        x = xb
        i = nx
        call bndyc(kbdy,yc,alfjm1,betjm1,gamjm1,gbdjm1)
        call bndyc(kbdy,yc+dly,alfj,betj,gamj,gbdj)
        gammax = max(abs(gamjm1),abs(gamj),gammax)
        do j=2,ny-1
          jc = j
          y = yc+j*dly
          call bndyc(kbdy,y,alfjp1,betjp1,gamjp1,gbdjp1)
          gammax = max(abs(gamjp1),gammax)
c
c     check for illegal tangential derivative b.c.
c
          if (alfjm1*alfj*alfjp1.eq.0.0) then
            ier = 13
            return
          end if
          c1=cf(nx,jc,1)
          c2=cf(nx,jc,2)
          c8=cf(nx,jc,8)
          cf(nx,jc,1)=0.0
          cf(nx,jc,2)=0.0
          cf(nx,jc,3)=cf(nx,jc,3)-c8*(-betjm1/alfjm1*dxoy)-c1*(betj/
     +    alfj*dxoy)-c2*(3.0*betjp1/alfjp1*dxoy+dlx2*gamjp1/alfjp1)
          cf(nx,jc,4)=cf(nx,jc,4)+c2
          cf(nx,jc,5)=cf(nx,jc,5)+c1
          cf(nx,jc,6)=cf(nx,jc,6)+c8
          cf(nx,jc,7)=cf(nx,jc,7)+c8*(3.0*betjm1/alfjm1*dxoy-dlx2*gamjm1
     +    /alfjm1)+c1*(betj/alfj*dxoy)+c2*(-betjp1/alfjp1*dxoy)
          cf(nx,jc,8)=0.0
          cf(nx,jc,9)=cf(nx,jc,9)-c8*(4.0*betjm1/alfjm1*dxoy)-
     +    c1*(dlx2*gamj/alfj)-c2*(-4.0*betjp1/alfjp1*dxoy)
          alfjm1 = alfj
          betjm1 = betj
          gamjm1 = gamj
          gbdjm1 = gbdj
          alfj = alfjp1
          betj = betjp1
          gamj = gamjp1
          gbdj = gbdjp1
        end do
      end if
      if (nyc.eq.2) then
        kbdy = 3
        y = yc
        j = 1
        jc = 1
        call bndyc(kbdy,xa,alfim1,betim1,gamim1,gbdim1)
        call bndyc(kbdy,xa+dlx,alfi,beti,gami,gbdi)
        gammax = max(abs(gamim1),abs(gami),gammax)
        do i=2,nx-1
          x=xa+i*dlx
          call bndyc(kbdy,x,alfip1,betip1,gamip1,gbdip1)
          gammax = max(abs(gamip1),gammax)
c
c     check for illegal tangential derivative b.c.
c
          if (betim1*beti*betip1.eq.0.0) then
            ier = 13
            return
          end if
          c6=cf(i,jc,6)
          c7=cf(i,jc,7)
          c8=cf(i,jc,8)
          cf(i,jc,1)=cf(i,jc,1)+c6*(-alfim1/betim1*dyox)+c7*(alfi/beti*
     +    dyox)+c8*(3.0*alfip1/betip1*dyox+dly2*gamip1/betip1)
          cf(i,jc,2)=cf(i,jc,2)+c8
          cf(i,jc,3)=cf(i,jc,3)+c7
          cf(i,jc,4)=cf(i,jc,4)+c6
          cf(i,jc,5)=cf(i,jc,5)+c6*(-3.0*alfim1/betim1*dyox+dly2*gamim1/
     +    betim1)+c7*(-alfi/beti*dyox)+c8*(alfip1/betip1*dyox)
          cf(i,jc,6)=0.0
          cf(i,jc,7)=0.0
          cf(i,jc,8)=0.0
          cf(i,jc,9)=cf(i,jc,9)+c6*(4.0*alfim1/betim1*dyox)+c7*dly2*gami
     +    /beti+c8*(-4.0*alfip1/betip1*dyox)
c     advance scalars for next pass
          alfim1=alfi
          betim1=beti
          gamim1=gami
          gbdim1 = gbdi
          alfi=alfip1
          beti=betip1
          gami=gamip1
          gbdi = gbdip1
        end do
      end if

      if (nyd.eq.2) then
        kbdy = 4
        y = yd
        j = ny
        jc = ny
        kbdy=4
        call bndyc(kbdy,xa,alfim1,betim1,gamim1,gbdim1)
        call bndyc(kbdy,xa+dlx,alfi,beti,gami,gbdi)
        gammax = max(abs(gamim1),abs(gami),gammax)
        do i=2,nx-1
          x=xa+i*dlx
          call bndyc(kbdy,x,alfip1,betip1,gamip1,gbdip1)
          gammax = max(abs(gamip1),gammax)
c
c     check for illegal tangential derivative b.c.
c
          if (betim1*beti*betip1.eq.0.0) then
            ier = 13
            return
          end if
          c2=cf(i,jc,2)
          c3=cf(i,jc,3)
          c4=cf(i,jc,4)
          cf(i,jc,1)=cf(i,jc,1)+c4*(alfim1/betim1*dyox)+c3*(-alfi/beti*
     +    dyox)+c2*(-3.*alfip1/betip1*dyox-dly2*gamip1/betip1)
c     set virtual point coefficients to zero
          cf(i,jc,2)=0.0
          cf(i,jc,3)=0.0
          cf(i,jc,4)=0.0
          cf(i,jc,5)=cf(i,jc,5)+c4*(3.0*alfim1/betim1*dyox-dly2*gamim1/
     +    betim1)+c3*(alfi/beti*dyox)+c2*(-alfip1/betip1*dyox)
          cf(i,jc,6)=cf(i,jc,6)+c4
          cf(i,jc,7)=cf(i,jc,7)+c3
          cf(i,jc,8)=cf(i,jc,8)+c2
          cf(i,jc,9)=cf(i,jc,9)+c4*(-4.*alfim1/betim1*dyox)+
     +    c3*(-dly2*gami/beti)+c2*(4.*alfip1/betip1*dyox)
          alfim1=alfi
          betim1=beti
          gamim1=gami
          gbdim1 = gbdi
          alfi=alfip1
          beti=betip1
          gami=gamip1
          gbdi = gbdip1
        end do
      end if
c
c     flag singular pde
c
      if (cemax.eq.0.0.and.alfmax.eq.0.0) then
        if (nxa.eq.0.or.(nxa.eq.2.and.nxb.eq.2)) then
          if (nyc.eq.0.or.(nyc.eq.2.and.nyd.eq.2)) then
            ier = -3
          end if
        end if
      end if
c
c     flag non-ellipticity
c
      if (cmin.le.0.0) then
        ier = -2
      end if
c
c     make cornor adjustments if necessary
c
      if (nyc.eq.2) then
        j = 1
        jc = 1
        if (nxa.eq.0) then
c
c     periodic-mixed at (xa,yc)
c
          kbdy = 3
          call bndyc(kbdy,xa,alfi,beti,gami,gbdi)
          call bndyc(kbdy,xa+dlx,alfip1,betip1,gamip1,gbdip1)
          c7=cf(1,jc,7)
          c8=cf(1,jc,8)
          cf(1,jc,1)=cf(1,jc,1)+c7*(alfi/beti*dyox)+
     +    c8*(3.0*alfip1/betip1*dyox+dly2*gamip1/betip1)
          cf(1,jc,2)=cf(1,jc,2)+c8
          cf(1,jc,3)=cf(1,jc,3)+c7
          cf(1,jc,5)=cf(1,jc,5)+c7*(-alfi/beti*dyox)+
     +    c8*(alfip1/betip1*dyox)
          cf(1,jc,6)=0.0
          cf(1,jc,7)=0.0
          cf(1,jc,8)=0.0
          cf(1,jc,9)=cf(1,jc,9)+c8*(-4.0*alfip1/betip1*dyox)+c7*(
     +    dly2*gami/beti)
c
c     adjust periodic-mixed at (xb,yc)
c
          call bndyc(kbdy,xb-dlx,alfim1,betim1,gamim1,gbdim1)
          call bndyc(kbdy,xb,alfi,beti,gami,gbdi)
          c6=cf(nx,jc,6)
          c7=cf(nx,jc,7)
          cf(nx,jc,1)=cf(nx,jc,1)+c6*(-alfim1/betim1*dyox)+c7*(alfi/beti
     +    *dyox)
          cf(nx,jc,3)=cf(nx,jc,3)+c7
          cf(nx,jc,4)=cf(nx,jc,4)+c6
          cf(nx,jc,5)=cf(nx,jc,5)+c6*(-3.0*alfim1/betim1*dyox+dly2*
     +    gamim1/betim1)+c7*(-alfi/beti*dyox)
          cf(nx,jc,6)=0.0
          cf(nx,jc,7)=0.0
          cf(nx,jc,8)=0.0
          cf(nx,jc,9)=cf(nx,jc,9)+c6*(4.0*alfim1/betim1*dyox)+c7*
     +    (dly2*gami/beti)
        else if (nxa.eq.2) then
c
c     mixed-mixed at (xa,yc)
c
        kbdy = 3
c     phase 1
        call bndyc(kbdy,xa+dlx,alfac,betac,gamac,gbdyc )
        kbdy = 1
        call bndyc(kbdy,yc+dly,alfaa,betaa,gamaa,gbdya)
        c4=cf(1,jc,4)
        c8=cf(1,jc,8)
        cf(1,jc,1)=cf(1,jc,1)+c8*(3.0*alfac/betac*dyox+dly2*gamac/
     +    betac)
        cf(1,jc,2)=cf(1,jc,2)+c8 +c4
        cf(1,jc,3)=cf(1,jc,3)+c4*(3.0*betaa/alfaa*dxoy+dlx2*gamaa
     +    /alfaa)
        cf(1,jc,4)=0.0
        cf(1,jc,5)=cf(1,jc,5)+c8*(alfac/betac*dyox)
        cf(1,jc,6)=0.0
        cf(1,jc,7)=cf(1,jc,7)+c4*(betaa/alfaa*dxoy)
        cf(1,jc,8)=0.0
        cf(1,jc,9)=cf(1,jc,9)+c8*(-4.0*alfac/betac*dyox)+
     +    c4*(-4.0*betaa/alfaa*dxoy)
c     phase 2
        c5=cf(1,jc,5)
        c7=cf(1,jc,7)
        kbdy = 3
        call bndyc(kbdy,xa,alfac,betac,gamac,gbdyc)
        kbdy = 1
        call bndyc(kbdy,yc,alfaa,betaa,gamaa,gbdya)
        det=alfaa*betac-betaa*alfac
        if (det.eq.0.0) then
          ier = 14
          return
        end if
        cf(1,jc,1)=cf(1,jc,1)+c5
        cf(1,jc,3)=cf(1,jc,3)+c7
        cf(1,jc,5)=0.0
        cf(1,jc,7)=0.0
        cf(1,jc,9)=cf(1,jc,9)+c5*(dlx2*(gamaa*betac-betaa*gamac)/det)+
     +    c7*(dly2*(alfaa*gamac-gamaa*alfac)/det)
      end if

      if (nxb.eq.2) then
c
c     mixed-mixed at (xb,yc)
c     phase 1
c
        kbdy=2
        call bndyc(kbdy,yc+dly,alfab,betab,gamab,gbdyb)
        kbdy=3
        call bndyc(kbdy,xb-dlx,alfac,betac,gamac,gbdyc)
        c2=cf(nx,jc,2)
        c6=cf(nx,jc,6)
        cf(nx,jc,1)=cf(nx,jc,1)+c6*(-alfac/betac*dyox)
        cf(nx,jc,2)=0.0
        cf(nx,jc,3)=cf(nx,jc,3)+c2*(-3.0*betab/alfab*dxoy-dlx2*
     +    gamab/alfab)
        cf(nx,jc,4)=cf(nx,jc,4)+c6+c2
        cf(nx,jc,5)=cf(nx,jc,5)+c6*(-3.0*alfac/betac*dyox+dly2*
     +    gamac/betac)
        cf(nx,jc,6)=0.0
        cf(nx,jc,7)=cf(nx,jc,7)-c2*(betab/alfab*dxoy)
        cf(nx,jc,8)=0.0
        cf(nx,jc,9)=cf(nx,jc,9)+c6*(4.0*alfac/betac*dyox)+
     +    c2*( 4.0*betab/alfab*dxoy)
c     phase 2
        kbdy=3
        call bndyc(kbdy,xb,alfac,betac,gamac,gbdyc)
        kbdy=2
        call bndyc(kbdy,yc,alfab,betab,gamab,gbdyb)
        det=betac*alfab-alfac*betab
        if (det.eq.0.0) then
          ier = 14
          return
        end if
        c1=cf(nx,jc,1)
        c7=cf(nx,jc,7)
        cf(nx,jc,1)=0.0
        cf(nx,jc,3)=cf(nx,jc,3)+c7
        cf(nx,jc,5)=cf(nx,jc,5)+c1
        cf(nx,jc,7)=0.0
        cf(nx,jc,9)=cf(nx,jc,9)+c1*(dlx2*(betab*gamac-gamab*betac)/
     +    det)+c7*(dly2*(alfab*gamac-gamab*alfac)/det)
      end if
      end if

      if (nyd.eq.2) then
      j = ny
      jc = ny
      if (nxa.eq.0) then
c     periodic-mixed at (xa,yd) and (xb,yd)
        kbdy=4
        call bndyc(kbdy,xa,alfi,beti,gami,gbdi)
        call bndyc(kbdy,xa+dlx,alfip1,betip1,gamip1,gbdip1)
        c2=cf(1,jc,2)
        c3=cf(1,jc,3)
        cf(1,jc,1)=cf(1,jc,1)+c3*(-alfi/beti*dyox)+
     +    c2*(-3.0*alfip1/betip1*dyox-dly2*gamip1/betip1)
        cf(1,jc,2) = 0.0
        cf(1,jc,3) = 0.0
        cf(1,jc,4) = 0.0
        cf(1,jc,5)=cf(1,jc,5)+c3*(alfi/beti*dyox)+
     +    c2*(-alfip1/betip1*dyox)
        cf(1,jc,7)=cf(1,jc,7)+c3
        cf(1,jc,8)=cf(1,jc,8)+c2
        cf(1,jc,9) = cf(1,jc,9)+c3*(-dly2*gami/beti)+
     +    c2*(4.0*alfip1/betip1*dyox)
        call bndyc(kbdy,xb-dlx,alfim1,betim1,gamim1,gbdim1)
        call bndyc(kbdy,xb,alfi,beti,gami,gbdi)
        c3=cf(nx,jc,3)
        c4=cf(nx,jc,4)
        cf(nx,jc,1)=cf(nx,jc,1)+c4*(alfim1/betim1*dyox)+
     +    c3*(-alfi/beti*dyox)
        cf(nx,jc,2)=0.0
        cf(nx,jc,3)=0.0
        cf(nx,jc,4)=0.0
        cf(nx,jc,5)=cf(nx,jc,5)+c4*(3.0*alfim1/betim1*dyox-dly2*
     +    gamim1/betim1)+c3*(alfi/beti*dyox)
        cf(nx,jc,6)=cf(nx,jc,6)+c4
        cf(nx,jc,7)=cf(nx,jc,7)+c3
        cf(nx,jc,9)=cf(nx,jc,9)+c4*(-4.0*alfim1/betim1*dyox)+
     +    c3*(-dly2*gami/beti)
      else if (nxa.eq.2) then
c     mixed-mixed at (xa,yd)
c     phase 1
        kbdy=4
        call bndyc(kbdy,xa+dlx,alfad,betad,gamad,gbdyd)
        kbdy=1
        call bndyc(kbdy,yd-dly,alfaa,betaa,gamaa,gbdya)
        c2=cf(1,jc,2)
        c6=cf(1,jc,6)
        cf(1,jc,1)=cf(1,jc,1)+c2*(-3.0*alfad/betad*dyox-dly2*
     +    gamad/betad)
        cf(1,jc,2)=0.0
        cf(1,jc,3)=cf(1,jc,3)+c6*(-betaa/alfaa*dxoy)
        cf(1,jc,4)=0.0
        cf(1,jc,5)=cf(1,jc,5)+c2*(-alfad/betad*dyox)
        cf(1,jc,6)=0.0
        cf(1,jc,7)=cf(1,jc,7)+c6*(-3.0*betaa/alfaa*dxoy+dlx2*
     +    gamaa/alfaa)
        cf(1,jc,8)=cf(1,jc,8)+c6+c2
        cf(1,jc,9)=cf(1,jc,9)+c2*(4.0*alfad/betad*dyox)+
     +    c6*(4.0*betaa/alfaa*dxoy)
c     phase 2
        kbdy=1
        call bndyc(kbdy,yd,alfaa,betaa,gamaa,gbdya)
        kbdy=4
        call bndyc(kbdy,xa,alfad,betad,gamad,gbdyd)
        det=alfad*betaa-betad*alfaa
        if (det.eq.0.0) then
          ier = 14
          return
        end if
        c3=cf(1,jc,3)
        c5=cf(1,jc,5)
        cf(1,jc,1)=cf(1,jc,1)+c5
        cf(1,jc,3)=0.0
        cf(1,jc,5)=0.0
        cf(1,jc,7)=cf(1,jc,7)+c3
        cf(1,jc,9)=cf(1,jc,9)+c5*(dlx2*(betaa*gamad-gamaa*betad)/det)+
     +    c3*(dly2*(alfaa*gamad-gamaa*alfad)/det)
      end if

      if (nxb.eq.2) then
c     mixed-mixed at (xb,yd)
        kbdy=4
        call bndyc(kbdy,xb-dlx,alfad,betad,gamad,gbdyd)
        kbdy=2
        call bndyc(kbdy,yd-dly,alfab,betab,gamab,gbdyb)
        c4=cf(nx,jc,4)
        c8=cf(nx,jc,8)
        cf(nx,jc,1)=cf(nx,jc,1)+c4*(alfad/betad*dyox)
        cf(nx,jc,2) = 0.0
        cf(nx,jc,3) = cf(nx,jc,3)+c8*(betab/alfab*dxoy)
        cf(nx,jc,4) = 0.0
        cf(nx,jc,5) = cf(nx,jc,5)+c4*(3.0*alfad/betad*dyox-dly2*
     +    gamad/betad)
        cf(nx,jc,6) = cf(nx,jc,6)+c4+c8
        cf(nx,jc,7) = cf(nx,jc,7)+c8*(3.0*betab/alfab*dxoy-dlx2*
     +    gamab/alfab)
        cf(nx,jc,8) = 0.0
        cf(nx,jc,9) = cf(nx,jc,9)+c4*(-4.0*alfad/betad*dyox)+
     +    c8*(-4.0*betab/alfab*dxoy)
c     phase 2
        kbdy=4
        call bndyc(kbdy,xb,alfad,betad,gamad,gbdyd)
        kbdy=2
        call bndyc(kbdy,yd,alfab,betab,gamab,gbdyb)
        det=alfad*betab-betad*alfab
        if (det.eq.0.0) then
          ier = 14
          return
        end if
        c1=cf(nx,jc,1)
        c3=cf(nx,jc,3)
        cf(nx,jc,1)=0.0
        cf(nx,jc,3)=0.0
        cf(nx,jc,5)=cf(nx,jc,5)+c1
        cf(nx,jc,7)=cf(nx,jc,7)+c3
        cf(nx,jc,9)=cf(nx,jc,9)+c1*(dlx2*(betad*gamab-gamad*betab)/
     +    det)+c3*(dly2*(alfab*gamad-gamab*alfad)/det)
      end if
      end if

      if (nxa.eq.2.and.nyc.eq.0) then
c
c     mixed-periodic at (xa,yc)
c
      kbdy = 1
      call bndyc(kbdy,yc,alfj,betj,gamj,gbdj)
      call bndyc(kbdy,yc+dly,alfjp1,betjp1,gamjp1,gbdjp1)
      c5 = cf(1,1,5)
      c4 = cf(1,1,4)
      cf(1,1,3) = cf(1,1,3)+c5*(betj/alfj*dxoy)+
     +  c4*(3.*betjp1/alfjp1*dxoy+dlx2*gamjp1/alfjp1)
      cf(1,1,2) = cf(1,1,2)+c4
      cf(1,1,1) = cf(1,1,1)+c5
      cf(1,1,7) = cf(1,1,7)+c5*(-betj/alfj*dxoy)+
     +  c4*(betjp1/alfjp1*dxoy)
      cf(1,1,6) = 0.0
      cf(1,1,5) = 0.0
      cf(1,1,4) = 0.0
      cf(1,1,9) = cf(1,1,9)+c4*(-4.*betjp1/alfjp1*dxoy)+c5*(dlx2*
     +  gamj/alfj)
c
c     adjust mixed-periodic at (xa,yd)
c
      call bndyc(kbdy,yd-dly,alfjm1,betjm1,gamjm1,gbdjm1)
      call bndyc(kbdy,yd,alfj,betj,gamj,gbdj)
      c6=cf(1,ny,6)
      c5=cf(1,ny,5)
      cf(1,ny,3) = cf(1,ny,3)+c6*(-betjm1/alfjm1*dxoy)+c5*(betj/alfj*
     +  dxoy)
      cf(1,ny,1) = cf(1,ny,1)+c5
      cf(1,ny,8) = cf(1,ny,8)+c6
      cf(1,ny,7) = cf(1,ny,7)+c6*(-3.*betjm1/alfjm1*dxoy+dlx2*
     +  gamjm1/alfjm1)+c5*(-betj/alfj*dxoy)
      cf(1,ny,6) = 0.0
      cf(1,ny,5) = 0.0
      cf(1,ny,4) = 0.0
      cf(1,ny,9) = cf(1,ny,9)+c6*(4.*betjm1/alfjm1*dxoy)+c5*
     +  (dlx2*gamj/alfj)
      end if

      if (nxb.eq.2.and.nyc.eq.0) then
c
c     mixed-periodic at (xb,yc) and (xb,yd)
c
      kbdy = 2
      call bndyc(kbdy,yc,alfj,betj,gamj,gbdj)
      call bndyc(kbdy,yc+dly,alfjp1,betjp1,gamjp1,gbdjp1)
      c2 = cf(nx,1,2)
      c1 = cf(nx,1,1)
      cf(nx,1,3) = cf(nx,1,3)+c1*(-betj/alfj*dxoy)+
     +  c2*(-3.*betjp1/alfjp1*dxoy-dlx2*gamjp1/alfjp1)
      cf(nx,1,2) = 0.0
      cf(nx,1,1) = 0.0
      cf(nx,1,8) = 0.0
      cf(nx,1,7) = cf(nx,1,7)+c1*(betj/alfj*dxoy)+
     +  c2*(-betjp1/alfjp1*dxoy)
      cf(nx,1,5) = cf(nx,1,5)+c1
      cf(nx,1,4) = cf(nx,1,4)+c2
      cf(nx,1,9) = cf(nx,1,9)+c1*(-dlx2*gamj/alfj)+
     +  c2*(4.*betjp1/alfjp1*dxoy)
      call bndyc(kbdy,yd-dly,alfjm1,betjm1,gamjm1,gbdjm1)
      call bndyc(kbdy,yd,alfj,betj,gamj,gbdj)
      c1 = cf(nx,ny,1)
      c8 = cf(nx,ny,8)
      cf(nx,ny,3) = cf(nx,ny,3)+c8*(betjm1/alfjm1*dxoy)+
     +  c1*(-betj/alfj*dxoy)
      cf(nx,ny,2) = 0.0
      cf(nx,ny,1) = 0.0
      cf(nx,ny,8) = 0.0
      cf(nx,ny,7) = cf(nx,ny,7)+c8*(3.*betjm1/alfjm1*dxoy-dlx2*
     +  gamjm1/alfjm1)+c1*(betj/alfj*dxoy)
      cf(nx,ny,6) = cf(nx,ny,6)+c8
      cf(nx,ny,5) = cf(nx,ny,5)+c1
      cf(nx,ny,9) = cf(nx,ny,9)+c8*(-4.*betjm1/alfjm1*dxoy)+
     +  c1*(-dlx2*gamj/alfj)
      end if
c
c     set coefficient for specified boundaries
c
      if (nxa.eq.1) then
      i = 1
      do j=1,ny
        do l=1,9
          cf(i,j,l) = 0.0
        end do
        cf(i,j,9) = 1.0
      end do
      end if
      if (nxb.eq.1) then
      i = nx
      do j=1,ny
        do l=1,9
          cf(i,j,l) = 0.0
        end do
        cf(i,j,9) = 1.0
      end do
      end if
      if (nyc.eq.1) then
      j = 1
      do i=1,nx
        do l=1,9
          cf(i,j,l) = 0.0
        end do
        cf(i,j,9) = 1.0
      end do
      end if
      if (nyd.eq.1) then
      j = ny
      do i=1,nx
        do l=1,9
          cf(i,j,l) = 0.0
        end do
        cf(i,j,9) = 1.0
      end do
      end if
c
c     set and factor tridiagonal matrices for line relaxation(s) if flagged
c
      if (method.eq.1.or.method.eq.3) then
      if (nxa.ne.0) then
c
c    nonperiodic x line relaxation
c
        do i=1,nx
          im1 = max0(i-1,1)
          do j=1,ny
            tx(im1,j,1) = cf(i,j,5)
            tx(i,j,2) = cf(i,j,9)
            tx(i,j,3) = cf(i,j,1)
          end do
        end do
        call factri(ny,nx,tx(1,1,1),tx(1,1,2),tx(1,1,3))
      else
c
c     periodic x line relaxation
c
        if (nx .gt. 3) then
c
c     set and factor iff nx > 3
c
          do i=1,nx-1
            do j=1,ny
            tx(i,j,1) = cf(i,j,5)
            tx(i,j,2) = cf(i,j,9)
            tx(i,j,3) = cf(i,j,1)
            end do
          end do
          call factrp(ny,nx,tx,tx(1,1,2),tx(1,1,3),tx(1,1,4),
     +                  tx(1,1,5),wk(kps))
        end if
      end if
      end if

      if (method.eq.2.or.method.eq.3) then
      if (nyc.ne.0) then
c
c     nonperiodic y line relaxation
c
        do j=1,ny
          jm1 = max0(j-1,1)
          do i=1,nx
            ty(jm1,i,1) = cf(i,j,7)
            ty(j,i,2) = cf(i,j,9)
            ty(j,i,3) = cf(i,j,3)
          end do
        end do
        call factri(nx,ny,ty(1,1,1),ty(1,1,2),ty(1,1,3))
      else
c
c      periodic y line relaxation
c
        if (ny .gt. 3) then
c
c     set and factor iff ny > 3
c
          do j=1,ny-1
            do i=1,nx
            ty(j,i,1) = cf(i,j,7)
            ty(j,i,2) = cf(i,j,9)
            ty(j,i,3) = cf(i,j,3)
            end do
          end do
          call factrp(nx,ny,ty,ty(1,1,2),ty(1,1,3),ty(1,1,4),
     +                  ty(1,1,5),wk(kps))
        end if
      end if
      end if
      return
      end subroutine dismd2cr
c
c-----------------------------------------------------------------------
c
      subroutine adjmd2cr(nx,ny,phi,cf)
c
c     adjust righthand side in cf(i,j,10) for boundary conditions
c
      implicit none
      integer nx,ny,i,j,kbdy
      real(r8) cf(nx,ny,10),phi(0:nx+1,0:ny+1)
      real(r8) dlx,dlx2,dlxx,dly,dly2,dlyy,dlxy,dlxy2,dlxy4,dxoy,dyox
      real(r8) x,y,cxx,cxy,cyy,cx,cy,ce,c1,c2,c3,c4,c5
      real(r8) c6,c7,c8
      real(r8) alfaa,alfab,alfac,alfad,betaa,betab,betac,betad,det
      real(r8) gamaa,gamab,gamac,gamad
      real(r8) alfim1,alfi,alfip1,betim1,beti,betip1,gamim1,gami,gamip1
      real(r8) alfjm1,alfj,alfjp1,betjm1,betj,betjp1,gamjm1,gamj,gamjp1
      real(r8) gbdim1,gbdi,gbdip1,gbdj,gbdjm1,gbdjp1
      real(r8) gbdya,gbdyb,gbdyc,gbdyd
c
      dlx = (xb-xa)/(nx-1)
      dlx2 = dlx+dlx
      dlxx = dlx*dlx
      dly = (yd-yc)/(ny-1)
      dly2 = dly+dly
      dlyy = dly*dly
      dlxy = dlx*dly
      dlxy2 = dlxy+dlxy
      dlxy4 = dlxy2+dlxy2
      dyox = dly/dlx
      dxoy = dlx/dly
c
c     adjust at derivative boundaries
c
      if (nyc.eq.2) then
      kbdy=3
      call bndyc(kbdy,xa,alfim1,betim1,gamim1,gbdim1)
      call bndyc(kbdy,xa+dlx,alfi,beti,gami,gbdi)
      do i=2,nx-1
        x=xa+i*dlx
        call bndyc(kbdy,x,alfip1,betip1,gamip1,gbdip1)
        call coef(x-dlx,yc,cxx,cxy,cyy,cx,cy,ce)
        cxx = max(cxx,abs(cx)*dlx*0.5)
        cyy = max(cyy,abs(cy)*dly*0.5)
        c6=cxy/dlxy4
        c7=cyy/dlyy-cy/dly2
        c8=-c6
        cf(i,1,10)=cf(i,1,10)+dly2*(c6*gbdim1/betim1+c7*gbdi/beti+
     +    c8*gbdip1/betip1)
        betim1=beti
        gbdim1=gbdi
        beti=betip1
        gbdi=gbdip1
      end do
      if (nxa.eq.0) then
        kbdy=3
        call bndyc(kbdy,xa,alfi,beti,gami,gbdi)
        call bndyc(kbdy,xa+dlx,alfip1,betip1,gamip1,gbdip1)
        call coef(xa,yc,cxx,cxy,cyy,cx,cy,ce)
        cxx = max(cxx,abs(cx)*dlx*0.5)
        cyy = max(cyy,abs(cy)*dly*0.5)
        c7=cyy/dlyy-cy/dly2+cxy/dlxy2
        c8=-cxy/dlxy2
        cf(1,1,10)=cf(1,1,10)+dly2*(c8*gbdip1/betip1+c7*gbdi/beti)
        call bndyc(kbdy,xb-dlx,alfim1,betim1,gamim1,gbdim1)
        call bndyc(kbdy,xb,alfi,beti,gami,gbdi)
        call coef(xb,yc,cxx,cxy,cyy,cx,cy,ce)
        cxx = max(cxx,abs(cx)*dlx*0.5)
        cyy = max(cyy,abs(cy)*dly*0.5)
        c6=cxy/dlxy2
        c7=cyy/dlyy-cy/dly2-c6
        cf(nx,1,10)=cf(nx,1,10)+dly2*(c6*gbdim1/betim1+c7*gbdi/beti)
      else if (nxa.eq.2) then
        kbdy=3
        call bndyc(kbdy,xa+dlx,alfac,betac,gamac,gbdyc)
        kbdy=1
        call bndyc(kbdy,yc+dly,alfaa,betaa,gamaa,gbdya)
        call coef(xa,yc,cxx,cxy,cyy,cx,cy,ce)
        cxx = max(cxx,abs(cx)*dlx*0.5)
        cyy = max(cyy,abs(cy)*dly*0.5)
        c4=-cxy/dlxy2
        c8=c4
        c5=cxx/dlxx-cx/dlx2-c4
        c7=cyy/dlyy-cy/dly2-c4
        c5=c5+c8*(alfac/betac*dyox)
        c7=c7+c4*(betaa/alfaa*dxoy)
        cf(1,1,10)=cf(1,1,10)+2.0*(c8*dly*gbdyc/betac+c4*dlx*gbdya/
     +    alfaa)
        kbdy=3
        call bndyc(kbdy,xa,alfac,betac,gamac,gbdyc)
        kbdy=1
        call bndyc(kbdy,yc,alfaa,betaa,gamaa,gbdya)
        det=alfaa*betac-betaa*alfac
        cf(1,1,10)=cf(1,1,10)-c5*(dlx2*(betaa*gbdyc-gbdya*betac)/det)-
     +    c7*(dly2*(alfac*gbdya-gbdyc*alfaa)/det)
      end if
      if (nxb.eq.2) then
c     correct for mixed-mixed at (xb,yc)
        kbdy=3
        call bndyc(kbdy,xb-dlx,alfac,betac,gamac,gbdyc)
        kbdy=2
        call bndyc(kbdy,yc+dly,alfab,betab,gamab,gbdyb)
        call coef(xb,yc,cxx,cxy,cyy,cx,cy,ce)
        cxx = max(cxx,abs(cx)*dlx*0.5)
        cyy = max(cyy,abs(cy)*dly*0.5)
        c2=cxy/dlxy2
        c6=c2
        c1=cxx/dlxx+cx/dlx2-c2
        c7=cyy/dlyy-cy/dly2-c2
        c1=c1+c6*(-alfac/betac*dyox)
        c7=c7-c2*(betab/alfab*dxoy)
        cf(nx,1,10)=cf(nx,1,10)+c6*(dly2*gbdyc/betac)-
     +    c2*(dlx2*gbdyb/alfab)
c     phase 2
        kbdy=3
        call bndyc(kbdy,xb,alfac,betac,gamac,gbdyc)
        kbdy=2
        call bndyc(kbdy,yc,alfab,betab,gamab,gbdyb)
        det=betac*alfab-alfac*betab
        cf(nx,1,10)=cf(nx,1,10)+c1*(dlx2*(betab*gbdyc-gbdyb*betac)/
     +    det)-c7*(dly2*(alfac*gbdyb-gbdyc*alfab)/det)
      end if
      end if

      if (nxb.eq.2) then
c     mixed along x=xb
      kbdy=2
      call bndyc(kbdy,yc,alfjm1,betjm1,gamjm1,gbdjm1)
      call bndyc(kbdy,yc+dly,alfj,betj,gamj,gbdj)
      do j=2,ny-1
        y=yc+j*dly
        call bndyc(kbdy,y,alfjp1,betjp1,gamjp1,gbdjp1)
        call coef(xb,y-dly,cxx,cxy,cyy,cx,cy,ce)
        cxx = max(cxx,abs(cx)*dlx*0.5)
        cyy = max(cyy,abs(cy)*dly*0.5)
        c1=cxx/dlxx+cx/dlx2
        c2=cxy/dlxy4
        c8=-c2
        cf(nx,j,10)=cf(nx,j,10)-dlx2*(c8*gbdjm1/alfjm1+c1*gbdj/alfj+
     +    c2*gbdjp1/alfjp1)
        alfjm1=alfj
        gbdjm1=gbdj
        alfj=alfjp1
        gbdj=gbdjp1
      end do
      end if

      if (nyd.eq.2) then
      kbdy=4
      call bndyc(kbdy,xa,alfim1,betim1,gamim1,gbdim1)
      call bndyc(kbdy,xa+dlx,alfi,beti,gami,gbdi)
      do i=2,nx-1
        x=xa+i*dlx
        call bndyc(kbdy,x,alfip1,betip1,gamip1,gbdip1)
        call coef(x-dlx,yd,cxx,cxy,cyy,cx,cy,ce)
        cxx = max(cxx,abs(cx)*dlx*0.5)
        cyy = max(cyy,abs(cy)*dly*0.5)
        c2=cxy/dlxy4
        c3=cyy/dlyy+cy/dly2
        c4=-c2
        cf(i,ny,10)=cf(i,ny,10)-dly2*(c4*gbdim1/betim1+c3*gbdi/beti+
     +    c2*gbdip1/betip1)
        betim1=beti
        gbdim1=gbdi
        beti=betip1
        gbdi=gbdip1
      end do
      if (nxa.eq.0) then
c     correct for periodic-mixed at (xa,yd), (xb,yd)
        kbdy=4
        call bndyc(kbdy,xa,alfi,beti,gami,gbdi)
        call bndyc(kbdy,xa+dlx,alfip1,betip1,gamip1,gbdip1)
        call coef(xa,yd,cxx,cxy,cyy,cx,cy,ce)
        cxx = max(cxx,abs(cx)*dlx*0.5)
        cyy = max(cyy,abs(cy)*dly*0.5)
        c2=cxy/dlxy2
        c3=cyy/dlyy+cy/dly2-c2
        cf(1,ny,10)=cf(1,ny,10)-dly2*(c3*gbdi/beti+c2*gbdip1/betip1)
        call bndyc(kbdy,xb-dlx,alfim1,betim1,gamim1,gbdim1)
        call bndyc(kbdy,xb,alfi,beti,gami,gbdi)
        call coef(xb,yd,cxx,cxy,cyy,cx,cy,ce)
        cxx = max(cxx,abs(cx)*dlx*0.5)
        cyy = max(cyy,abs(cy)*dly*0.5)
        c4=-cxy/dlxy2
        c3=cyy/dlyy+cy/dly2-c4
        cf(nx,ny,10)=cf(nx,ny,10)-dly2*(c3*gbdi/beti+c4*gbdim1/betim1)
      else if (nxa.eq.2) then
c     correct for mixed-mixed at (xa,yd)
c     phase 1
        call coef(xa,yd,cxx,cxy,cyy,cx,cy,ce)
        cxx = max(cxx,abs(cx)*dlx*0.5)
        cyy = max(cyy,abs(cy)*dly*0.5)
        c2=cxy/dlxy2
        c6=c2
        c3=cyy/dlyy+cy/dly2-c2
        c5=cxx/dlxx-cx/dlx2-c2
        kbdy=4
        call bndyc(kbdy,xa+dlx,alfad,betad,gamad,gbdyd)
        kbdy=1
        call bndyc(kbdy,yd-dly,alfaa,betaa,gamaa,gbdya)
        c3=c3+c6*(-betaa/alfaa*dxoy)
        c5=c5+c2*(-alfad/betad*dyox)
        cf(1,ny,10)=cf(1,ny,10)+c2*(-dly2*gbdyd/betad)+
     +    c6*(dlx2*gbdya/alfaa)
c     phase 2
        kbdy=4
        call bndyc(kbdy,xa,alfad,betad,gamad,gbdyd)
        kbdy=1
        call bndyc(kbdy,yd,alfaa,betaa,gamaa,gbdya)
        det=alfad*betaa-betad*alfaa
        cf(1,ny,10)=cf(1,ny,10)-c5*(dlx2*(betad*gbdya-gbdyd*betaa)/
     +    det)-c3*(dly2*(alfad*gbdya-gbdyd*alfaa)/det)
      end if
      if (nxb.eq.2) then
c     correct ofr mixed-mixed at (xb,yd)
        call coef(xb,yd,cxx,cxy,cyy,cx,cy,ce)
        cxx = max(cxx,abs(cx)*dlx*0.5)
        cyy = max(cyy,abs(cy)*dly*0.5)
        c4=-cxy/dlxy2
        c8=c4
        c1=cxx/dlxx+cx/dlx2-c4
        c3=cyy/dlyy+cy/dly2-c4
        kbdy=4
        call bndyc(kbdy,xb-dlx,alfad,betad,gamad,gbdyd)
        kbdy=2
        call bndyc(kbdy, yd-dly,alfab,betab,gamab,gbdyb)
        c1=c1+c4*(alfad/betad*dyox)
        c3=c3+c8*(betab/alfab*dxoy)
        cf(nx,ny,10)=cf(nx,ny,10)+c4*(-dly2*gbdyd/betad)+
     +    c8*(-dlx2*gbdyb/alfab)
c     phase 2
        kbdy=4
        call bndyc(kbdy,xb,alfad,betad,gamad,gbdyd)
        kbdy=2
        call bndyc(kbdy,yd,alfab,betab,gamab,gbdyb)
        det=alfad*betab-betad*alfab
        cf(nx,ny,10)=cf(nx,ny,10)-c1*(dlx2*(betab*gbdyd-gbdyb*betad)
     +    /det)-c3*(dly2*(alfad*gbdyb-gbdyd*alfab)/det)
      end if
      end if

      if (nxa.eq.2) then
c     mixed along x=xa
      kbdy=1
      call bndyc(kbdy,yc,alfjm1,betjm1,gamjm1,gbdjm1)
      call bndyc(kbdy,yc+dly,alfj,betj,gamj,gbdj)
      do j=2,ny-1
        y=yc+j*dly
        call bndyc(kbdy,y,alfjp1,betjp1,gamjp1,gbdjp1)
        call coef(xa,y-dly,cxx,cxy,cyy,cx,cy,ce)
        cxx = max(cxx,abs(cx)*dlx*0.5)
        cyy = max(cyy,abs(cy)*dly*0.5)
        c4=-cxy/dlxy4
        c5=cxx/dlxx-cx/dlx2
        c6=-c4
        cf(1,j,10)=cf(1,j,10)+dlx2*(c6*gbdjm1/alfjm1+c5*gbdj/alfj+
     +    c4*gbdjp1/alfjp1)
        alfjm1=alfj
        gbdjm1=gbdj
        alfj=alfjp1
        gbdj=gbdjp1
      end do
      end if

      if (nxa.eq.2.and.nyc.eq.0) then
c
c     mixed-periodic at (xa,yc)
c
      kbdy = 1
      call bndyc(kbdy,yc,alfj,betj,gamj,gbdj)
      call bndyc(kbdy,yc+dly,alfjp1,betjp1,gamjp1,gbdjp1)
      call coef(xa,yc,cxx,cxy,cyy,cx,cy,ce)
      cxx = max(cxx,abs(cx)*dlx*0.5)
      cyy = max(cyy,abs(cy)*dly*0.5)
      c4 = -cxy/dlxy2
      c5 = cxx/dlxx-cx/dlx2-c4
      cf(1,1,10) = cf(1,1,10)+dlx2*(c4*gbdjp1/alfjp1+c5*gbdj/alfj)
c
c     mixed-periodic at (xa,yd)
c
      call bndyc(kbdy,yd-dly,alfjm1,betjm1,gamjm1,gbdjm1)
      call bndyc(kbdy,yd,alfj,betj,gamj,gbdj)
      call coef(xa,yd,cxx,cxy,cyy,cx,cy,ce)
      cxx = max(cxx,abs(cx)*dlx*0.5)
      cyy = max(cyy,abs(cy)*dly*0.5)
      c6=cxy/dlxy2
      c5 = cxx/dlxx-cx/dlx2-c6
      cf(1,ny,10)=cf(1,ny,10)+dlx2*(c6*gbdjm1/alfjm1+c5*gbdj/alfj)
      end if

      if (nxb.eq.2.and.nyc.eq.0) then
c
c     mixed-periodic at (xb,yc) and (xb,yd)
c
      kbdy = 2
      call bndyc(kbdy,yc,alfj,betj,gamj,gbdj)
      call bndyc(kbdy,yc+dly,alfjp1,betjp1,gamjp1,gbdjp1)
      call coef(xb,yc,cxx,cxy,cyy,cx,cy,ce)
      cxx = max(cxx,abs(cx)*dlx*0.5)
      cyy = max(cyy,abs(cy)*dly*0.5)
      c2 = cxy/dlxy2
      c1 = cxx/dlxx+cx/dlx2-c2
      cf(nx,1,10) = cf(nx,1,10)-dlx2*(c1*gbdj/alfj+c2*gbdjp1/alfjp1)
      call bndyc(kbdy,yd-dly,alfjm1,betjm1,gamjm1,gbdjm1)
      call bndyc(kbdy,yd,alfj,betj,gamj,gbdj)
      call coef(xb,yd,cxx,cxy,cyy,cx,cy,ce)
      cxx = max(cxx,abs(cx)*dlx*0.5)
      cyy = max(cyy,abs(cy)*dly*0.5)
      c8 = -cxy/dlxy2
      c1 = cxx/dlxx+cx/dlx2-c8
      cf(nx,ny,10) = cf(nx,ny,10)-dlx2*(c1*gbdj/alfj+c8*gbdjm1/alfjm1)
      end if
c
c     set specified boundaries in rhs from phi
c
      if (nxa.eq.1) then
      i = 1
      do j=1,ny
        cf(i,j,10) = phi(i,j)
      end do
      end if
      if (nxb.eq.1) then
      i = nx
      do j=1,ny
        cf(i,j,10) = phi(i,j)
      end do
      end if
      if (nyc.eq.1) then
      j = 1
      do i=1,nx
        cf(i,j,10) = phi(i,j)
      end do
      end if
      if (nyd.eq.1) then
      j = ny
      do i=1,nx
        cf(i,j,10) = phi(i,j)
      end do
      end if
      return
      end subroutine adjmd2cr
c
c-----------------------------------------------------------------------
c
      subroutine resmd2cr(nx,ny,phi,ncx,ncy,phic,rhsc,cof,resf)
c
c     restrict residual from fine to coarse mesh using fully weighted
c     residual restriction
c
      implicit none
      integer nx,ny,ncx,ncy,i,j,ic,jc
      real(r8) rhsc(ncx,ncy),resf(nx,ny)
      real(r8) phi(0:nx+1,0:ny+1),phic(0:ncx+1,0:ncy+1)
      real(r8) cof(nx,ny,10)
c
c     set phic zero
c
      do jc=0,ncy+1
      do ic=0,ncx+1
        phic(ic,jc) = 0.0
      end do
      end do
c
c     compute residual on fine mesh in resf
c
!$OMP PARALLEL DO SHARED(resf,cof,phi,nx,ny) PRIVATE(i,j)
      do j=1,ny
      do i=1,nx
        resf(i,j) = cof(i,j,10)-(
     +                cof(i,j,1)*phi(i+1,j)+
     +                cof(i,j,2)*phi(i+1,j+1)+
     +                cof(i,j,3)*phi(i,j+1)+
     +                cof(i,j,4)*phi(i-1,j+1)+
     +                cof(i,j,5)*phi(i-1,j)+
     +                cof(i,j,6)*phi(i-1,j-1)+
     +                cof(i,j,7)*phi(i,j-1)+
     +                cof(i,j,8)*phi(i+1,j-1)+
     +                cof(i,j,9)*phi(i,j))
      end do
      end do
c
c     restrict resf to coarse mesh in rhsc
c
      call res2(nx,ny,resf,ncx,ncy,rhsc,nxa,nxb,nyc,nyd)
      return
      end subroutine resmd2cr
c
c-----------------------------------------------------------------------
c
      subroutine relmd2cr(nx,ny,phi,cof,tx,ty,sum)
c
c     relaxation for mud2
c
      implicit none
      integer nx,ny
      real(r8) phi(*),cof(*),tx(*),ty(*),sum(*)
c
      if (method.eq.0) then                ! point relaxation
      call relmd2crp(nx,ny,phi,cof)
      else if (method.eq.1) then           ! line x relaxation
      call slxmd2cr(nx,ny,phi,cof,tx,sum)
      else if (method.eq.2) then           ! line y relaxation
      call slymd2cr(nx,ny,phi,cof,ty,sum)
      else if (method.eq.3) then           ! line x&y relaxation
      call slxmd2cr(nx,ny,phi,cof,tx,sum)
      call slymd2cr(nx,ny,phi,cof,ty,sum)
      end if
      return
      end subroutine relmd2cr
c
c-----------------------------------------------------------------------
c
      subroutine relmd2crp(nx,ny,phi,cof)
c
c     gauss-seidel four color point relaxation
c
      implicit none
      integer nx,ny,i,j,lcolor,i1,i2,i3,i4,it
c
      real(r8) phi(0:nx+1,0:ny+1),cof(nx,ny,10)
      i1 = 1
      i2 = 4
      i3 = 3
      i4 = 2
c
c     sweep four colored grid points
c
      do lcolor=1,4
!$OMP PARALLEL DO SHARED(i1,cof,phi,nx,ny) PRIVATE(i,j)
      do j=1,ny,4
        do i=i1,nx,4
            phi(i,j) = (cof(i,j,10) - (
     +                    cof(i,j,1)*phi(i+1,j)   +
     +                    cof(i,j,2)*phi(i+1,j+1) +
     +                    cof(i,j,3)*phi(i,j+1)   +
     +                    cof(i,j,4)*phi(i-1,j+1) +
     +                    cof(i,j,5)*phi(i-1,j)   +
     +                    cof(i,j,6)*phi(i-1,j-1) +
     +                    cof(i,j,7)*phi(i,j-1)   +
     +                    cof(i,j,8)*phi(i+1,j-1)))/cof(i,j,9)
        end do
      end do
!$OMP PARALLEL DO SHARED(i2,cof,phi,nx,ny) PRIVATE(i,j)
      do j=2,ny,4
        do i=i2,nx,4
            phi(i,j) = (cof(i,j,10) - (
     +                    cof(i,j,1)*phi(i+1,j)   +
     +                    cof(i,j,2)*phi(i+1,j+1) +
     +                    cof(i,j,3)*phi(i,j+1)   +
     +                    cof(i,j,4)*phi(i-1,j+1) +
     +                    cof(i,j,5)*phi(i-1,j)   +
     +                    cof(i,j,6)*phi(i-1,j-1) +
     +                    cof(i,j,7)*phi(i,j-1)   +
     +                    cof(i,j,8)*phi(i+1,j-1)))/cof(i,j,9)
        end do
      end do
!$OMP PARALLEL DO SHARED(i3,cof,phi,nx,ny) PRIVATE(i,j)
      do j=3,ny,4
        do i=i3,nx,4
            phi(i,j) = (cof(i,j,10) - (
     +                    cof(i,j,1)*phi(i+1,j)   +
     +                    cof(i,j,2)*phi(i+1,j+1) +
     +                    cof(i,j,3)*phi(i,j+1)   +
     +                    cof(i,j,4)*phi(i-1,j+1) +
     +                    cof(i,j,5)*phi(i-1,j)   +
     +                    cof(i,j,6)*phi(i-1,j-1) +
     +                    cof(i,j,7)*phi(i,j-1)   +
     +                    cof(i,j,8)*phi(i+1,j-1)))/cof(i,j,9)
        end do
      end do
!$OMP PARALLEL DO SHARED(i4,cof,phi,nx,ny) PRIVATE(i,j)
      do j=4,ny,4
        do i=i4,nx,4
            phi(i,j) = (cof(i,j,10) - (
     +                    cof(i,j,1)*phi(i+1,j)   +
     +                    cof(i,j,2)*phi(i+1,j+1) +
     +                    cof(i,j,3)*phi(i,j+1)   +
     +                    cof(i,j,4)*phi(i-1,j+1) +
     +                    cof(i,j,5)*phi(i-1,j)   +
     +                    cof(i,j,6)*phi(i-1,j-1) +
     +                    cof(i,j,7)*phi(i,j-1)   +
     +                    cof(i,j,8)*phi(i+1,j-1)))/cof(i,j,9)
        end do
      end do
c
c     set periodic virtual boundaries as necessary
c
      if (nxa.eq.0) then
        do j=1,ny
          phi(0,j) = phi(nx-1,j)
          phi(nx+1,j) = phi(2,j)
        end do
      end if
      if (nyc.eq.0) then
        do i=1,nx
          phi(i,0) = phi(i,ny-1)
          phi(i,ny+1) = phi(i,2)
        end do
      end if
c
c    permute (i1,i2,i3,i4) for next color
c
      it = i4
      i4 = i3
      i3 = i2
      i2 = i1
      i1 = it
      end do
      return
      end subroutine relmd2crp
c
c-----------------------------------------------------------------------
c
      subroutine slxmd2cr(nx,ny,phi,cof,tx,sum)
c
c     line relaxation in the x direction (periodic or nonperiodic)
c
      implicit none
      integer nx,ny,i,ib,j
c
      real(r8) phi(0:nx+1,0:ny+1),cof(nx,ny,10),tx(nx,ny,*),sum(ny)
c
c     replace line x with point gauss-seidel if
c     x direction is periodic and nx = 3 (coarsest)
c
      if (nxa .eq. 0 .and. nx .eq. 3) then
      call relmd2crp(nx,ny,phi,cof)
      return
      end if
c
c     set periodic y virtual boundary if necessary
c
      if (nyc.eq.0) then
      do i=1,nx
        phi(i,0) = phi(i,ny-1)
        phi(i,ny+1) = phi(i,2)
      end do
      end if

      if (nxa.ne.0) then
c
c     x direction not periodic, sweep odd j lines
c
!$OMP PARALLEL DO SHARED(cof,phi,tx,nx,ny) PRIVATE(i,ib,j)
      do j=1,ny,2
        do i=1,nx
          phi(i,j) = cof(i,j,10)-(cof(i,j,2)*phi(i+1,j+1)+
     +                              cof(i,j,3)*phi(i,j+1)+
     +                              cof(i,j,4)*phi(i-1,j+1)+
     +                              cof(i,j,6)*phi(i-1,j-1)+
     +                              cof(i,j,7)*phi(i,j-1)+
     +                              cof(i,j,8)*phi(i+1,j-1))
        end do
c
c     forward sweep
c
        do i=2,nx
          phi(i,j) = phi(i,j)-tx(i-1,j,1)*phi(i-1,j)
        end do
c
c     backward sweep
c
        phi(nx,j) = phi(nx,j)/tx(nx,j,2)
        do ib=2,nx
          i = nx-ib+1
          phi(i,j) = (phi(i,j)-tx(i,j,3)*phi(i+1,j))/tx(i,j,2)
        end do
      end do
c
c     sweep even j lines forward and back
c
!$OMP PARALLEL DO SHARED(cof,phi,tx,nx,ny) PRIVATE(i,ib,j)
      do j=2,ny,2
        do i=1,nx
          phi(i,j) = cof(i,j,10)-(cof(i,j,2)*phi(i+1,j+1)+
     +                              cof(i,j,3)*phi(i,j+1)+
     +                              cof(i,j,4)*phi(i-1,j+1)+
     +                              cof(i,j,6)*phi(i-1,j-1)+
     +                              cof(i,j,7)*phi(i,j-1)+
     +                              cof(i,j,8)*phi(i+1,j-1))
        end do
        do i=2,nx
          phi(i,j) = phi(i,j)-tx(i-1,j,1)*phi(i-1,j)
        end do
        phi(nx,j) = phi(nx,j)/tx(nx,j,2)
        do ib=2,nx
          i = nx-ib+1
          phi(i,j) = (phi(i,j)-tx(i,j,3)*phi(i+1,j))/tx(i,j,2)
        end do
      end do
      else
c
c     x direction periodic
c
      do j=1,ny
        sum(j) = 0.0
        phi(0,j) = phi(nx-1,j)
        phi(nx+1,j) = phi(2,j)
      end do
c
c      sweep odd lines forward and back
c
!$OMP PARALLEL DO SHARED(sum,cof,phi,tx,nx,ny) PRIVATE(i,j,ib)
      do j=1,ny,2
        do i=1,nx-1
          phi(i,j) = cof(i,j,10)-(cof(i,j,2)*phi(i+1,j+1)+
     +                              cof(i,j,3)*phi(i,j+1)+
     +                              cof(i,j,4)*phi(i-1,j+1)+
     +                              cof(i,j,6)*phi(i-1,j-1)+
     +                              cof(i,j,7)*phi(i,j-1)+
     +                              cof(i,j,8)*phi(i+1,j-1))
        end do
c
c     forward sweep
c
        do i=2,nx-2
          phi(i,j) = phi(i,j)-tx(i,j,1)*phi(i-1,j)
        end do
        do i=1,nx-2
          sum(j) = sum(j)+tx(i,j,5)*phi(i,j)
        end do
        phi(nx-1,j) = phi(nx-1,j)-sum(j)
c
c     backward sweep
c
        phi(nx-1,j) = phi(nx-1,j)/tx(nx-1,j,2)
        phi(nx-2,j) = (phi(nx-2,j)-tx(nx-2,j,4)*phi(nx-1,j))/
     +                   tx(nx-2,j,2)
        do ib=4,nx
          i = nx-ib+1
          phi(i,j) = (phi(i,j)-tx(i,j,3)*phi(i+1,j)-tx(i,j,4)*
     +                 phi(nx-1,j))/tx(i,j,2)
        end do
      end do
c
c     set periodic and virtual points for j odd
c
      do j=1,ny,2
        phi(nx,j) = phi(1,j)
        phi(0,j) = phi(nx-1,j)
        phi(nx+1,j) = phi(2,j)
      end do
c
c     sweep even j lines
c
!$OMP PARALLEL DO SHARED(sum,cof,phi,tx,nx,ny) PRIVATE(i,j,ib)
      do j=2,ny,2
        do i=1,nx-1
          phi(i,j) = cof(i,j,10)-(cof(i,j,2)*phi(i+1,j+1)+
     +                              cof(i,j,3)*phi(i,j+1)+
     +                              cof(i,j,4)*phi(i-1,j+1)+
     +                              cof(i,j,6)*phi(i-1,j-1)+
     +                              cof(i,j,7)*phi(i,j-1)+
     +                              cof(i,j,8)*phi(i+1,j-1))
        end do
c
c     forward sweep
c
        do i=2,nx-2
          phi(i,j) = phi(i,j)-tx(i,j,1)*phi(i-1,j)
        end do
        do i=1,nx-2
          sum(j) = sum(j)+tx(i,j,5)*phi(i,j)
        end do
        phi(nx-1,j) = phi(nx-1,j)-sum(j)
c
c     backward sweep
c
        phi(nx-1,j) = phi(nx-1,j)/tx(nx-1,j,2)
        phi(nx-2,j) = (phi(nx-2,j)-tx(nx-2,j,4)*phi(nx-1,j))/
     +                   tx(nx-2,j,2)
        do ib=4,nx
          i = nx-ib+1
          phi(i,j) = (phi(i,j)-tx(i,j,3)*phi(i+1,j)-tx(i,j,4)*
     +                 phi(nx-1,j))/tx(i,j,2)
        end do
      end do
c
c     set periodic and virtual points for j even
c
      do j=2,ny,2
        phi(nx,j) = phi(1,j)
        phi(0,j) = phi(nx-1,j)
        phi(nx+1,j) = phi(2,j)
      end do
      end if
c
c     set periodic y virtual boundaries if necessary
c
      if (nyc.eq.0) then
      do i=1,nx
        phi(i,0) = phi(i,ny-1)
        phi(i,ny+1) = phi(i,2)
      end do
      end if
      return
      end subroutine slxmd2cr
c
c-----------------------------------------------------------------------
c
      subroutine slymd2cr(nx,ny,phi,cof,ty,sum)
      implicit none
      integer nx,ny,i,j,jb
c
      real(r8) phi(0:nx+1,0:ny+1),cof(nx,ny,10),ty(ny,nx,*),sum(nx)
c
c     replace line y with point gauss-seidel if
c     y direction is periodic and ny = 3
c
      if (nyc .eq. 0 .and. ny .eq. 3) then
      call relmd2crp(nx,ny,phi,cof)
      return
      end if
c
c      set periodic and virtual x boundaries if necessary
c
      if (nxa.eq.0) then
      do j=1,ny
        phi(0,j) = phi(nx-1,j)
        phi(nx,j) = phi(1,j)
        phi(nx+1,j) = phi(2,j)
      end do
      end if

      if (nyc.ne.0) then
c
c     y direction not periodic
c
!$OMP PARALLEL DO SHARED(cof,phi,ty,nx,ny) PRIVATE(i,j,jb)
      do i=1,nx,2
        do j=1,ny
          phi(i,j) = cof(i,j,10)-(cof(i,j,1)*phi(i+1,j)+
     +                              cof(i,j,2)*phi(i+1,j+1)+
     +                              cof(i,j,4)*phi(i-1,j+1)+
     +                              cof(i,j,5)*phi(i-1,j)+
     +                              cof(i,j,6)*phi(i-1,j-1)+
     +                              cof(i,j,8)*phi(i+1,j-1))
        end do
c
c     forward sweep thru odd x lines
c
        do j=2,ny
          phi(i,j) = phi(i,j)-ty(j-1,i,1)*phi(i,j-1)
        end do
c
c      backward sweep
c
        phi(i,ny) = phi(i,ny)/ty(ny,i,2)
        do jb=2,ny
          j = ny-jb+1
          phi(i,j) = (phi(i,j)-ty(j,i,3)*phi(i,j+1))/ty(j,i,2)
        end do
      end do
c
c     forward sweep even x lines
c
!$OMP PARALLEL DO SHARED(cof,phi,ty,nx,ny) PRIVATE(i,j,jb)
      do i=2,nx,2
        do j=1,ny
          phi(i,j) = cof(i,j,10)-(cof(i,j,1)*phi(i+1,j)+
     +                              cof(i,j,2)*phi(i+1,j+1)+
     +                              cof(i,j,4)*phi(i-1,j+1)+
     +                              cof(i,j,5)*phi(i-1,j)+
     +                              cof(i,j,6)*phi(i-1,j-1)+
     +                              cof(i,j,8)*phi(i+1,j-1))
        end do
        do j=2,ny
          phi(i,j) = phi(i,j)-ty(j-1,i,1)*phi(i,j-1)
        end do
c
c      backward sweep
c
        phi(i,ny) = phi(i,ny)/ty(ny,i,2)
        do jb=2,ny
          j = ny-jb+1
          phi(i,j) = (phi(i,j)-ty(j,i,3)*phi(i,j+1))/ty(j,i,2)
        end do
      end do
      else
c
c     y direction periodic
c
      do i=1,nx
        sum(i) = 0.0
        phi(i,0) = phi(i,ny-1)
        phi(i,ny) = phi(i,1)
        phi(i,ny+1) = phi(i,2)
      end do
c
c     forward sweep odd x lines
c
!$OMP PARALLEL DO SHARED(sum,cof,phi,ty,nx,ny) PRIVATE(i,j,jb)
      do i=1,nx,2
        do j=1,ny-1
          phi(i,j) = cof(i,j,10)-(cof(i,j,1)*phi(i+1,j)+
     +                              cof(i,j,2)*phi(i+1,j+1)+
     +                              cof(i,j,4)*phi(i-1,j+1)+
     +                              cof(i,j,5)*phi(i-1,j)+
     +                              cof(i,j,6)*phi(i-1,j-1)+
     +                              cof(i,j,8)*phi(i+1,j-1))
        end do
        do j=2,ny-2
          phi(i,j) = phi(i,j)-ty(j,i,1)*phi(i,j-1)
        end do
        do j=1,ny-2
          sum(i) = sum(i)+ty(j,i,5)*phi(i,j)
        end do
        phi(i,ny-1) = phi(i,ny-1)-sum(i)
c
c     backward sweep
c
        phi(i,ny-1) = phi(i,ny-1)/ty(ny-1,i,2)
        phi(i,ny-2) = (phi(i,ny-2)-ty(ny-2,i,4)*phi(i,ny-1))/
     +                   ty(ny-2,i,2)
        do jb=4,ny
          j = ny-jb+1
          phi(i,j) = (phi(i,j)-ty(j,i,3)*phi(i,j+1)-ty(j,i,4)*
     +                  phi(i,ny-1))/ty(j,i,2)
        end do
      end do
c
c       set odd periodic and virtual y boundaries
c
      do i=1,nx,2
        phi(i,0) = phi(i,ny-1)
        phi(i,ny) = phi(i,1)
        phi(i,ny+1) = phi(i,2)
      end do
c
c     forward sweep even x lines
c
!$OMP PARALLEL DO SHARED(sum,cof,phi,ty,nx,ny) PRIVATE(i,j,jb)
      do i=2,nx,2
        do j=1,ny-1
          phi(i,j) = cof(i,j,10)-(cof(i,j,1)*phi(i+1,j)+
     +                              cof(i,j,2)*phi(i+1,j+1)+
     +                              cof(i,j,4)*phi(i-1,j+1)+
     +                              cof(i,j,5)*phi(i-1,j)+
     +                              cof(i,j,6)*phi(i-1,j-1)+
     +                              cof(i,j,8)*phi(i+1,j-1))
        end do
        do j=2,ny-2
          phi(i,j) = phi(i,j)-ty(j,i,1)*phi(i,j-1)
        end do
        do j=1,ny-2
          sum(i) = sum(i)+ty(j,i,5)*phi(i,j)
        end do
        phi(i,ny-1) = phi(i,ny-1)-sum(i)
c
c     backward sweep
c
        phi(i,ny-1) = phi(i,ny-1)/ty(ny-1,i,2)
        phi(i,ny-2) = (phi(i,ny-2)-ty(ny-2,i,4)*phi(i,ny-1))/
     +                   ty(ny-2,i,2)
        do jb=4,ny
          j = ny-jb+1
          phi(i,j) = (phi(i,j)-ty(j,i,3)*phi(i,j+1)-ty(j,i,4)*
     +                  phi(i,ny-1))/ty(j,i,2)
        end do
      end do
c
c       set even periodic and virtual y boundaries
c
      do i=2,nx,2
        phi(i,0) = phi(i,ny-1)
        phi(i,ny) = phi(i,1)
        phi(i,ny+1) = phi(i,2)
      end do
      end if
c
c      set periodic and virtual x boundaries if necessary
c
      if (nxa.eq.0) then
      do j=1,ny
        phi(0,j) = phi(nx-1,j)
        phi(nx+1,j) = phi(2,j)
      end do
      end if
      return
      end subroutine slymd2cr
c
c-----------------------------------------------------------------------
c
      subroutine coef(x,y,cxx,cxy,cyy,cx,cy,ce)
c
      implicit none
c
      real(r8), intent(in)    :: x,y
      real(r8), intent(inout) :: cxx,cxy,cyy,cx,cy,ce
c
      integer :: i,j,i1,i2,j1,j2
      real(r8) :: faci1,faci2,facj1,facj2
c
      do i=1,nx0-1
        if ( xusr(i) <= x .and. x <= xusr(i+1) ) then
          i1 = i
          i2 = i+1
          faci1 = (xusr(i+1)-x)/(xusr(i+1)-xusr(i))
          faci2 = (x-xusr(i))/(xusr(i+1)-xusr(i))
          exit
        end if
      end do
c
      do j=1,ny0-1
        if ( yusr(j) <= y .and. y <= yusr(j+1) ) then
          j1 = j
          j2 = j+1
          facj1 = (yusr(j+1)-y)/(yusr(j+1)-yusr(j))
          facj2 = (y-yusr(j))/(yusr(j+1)-yusr(j))
          exit
        end if
      end do
c
      cxx = (cxxusr(i1,j1)*faci1+cxxusr(i2,j1)*faci2)*facj1 +
     +      (cxxusr(i1,j2)*faci1+cxxusr(i2,j2)*faci2)*facj2
      cxy = (cxyusr(i1,j1)*faci1+cxyusr(i2,j1)*faci2)*facj1 +
     +      (cxyusr(i1,j2)*faci1+cxyusr(i2,j2)*faci2)*facj2
      cyy = (cyyusr(i1,j1)*faci1+cyyusr(i2,j1)*faci2)*facj1 +
     +      (cyyusr(i1,j2)*faci1+cyyusr(i2,j2)*faci2)*facj2
      cx  = (cxusr (i1,j1)*faci1+cxusr (i2,j1)*faci2)*facj1 +
     +      (cxusr (i1,j2)*faci1+cxusr (i2,j2)*faci2)*facj2
      cy  = (cyusr (i1,j1)*faci1+cyusr (i2,j1)*faci2)*facj1 +
     +      (cyusr (i1,j2)*faci1+cyusr (i2,j2)*faci2)*facj2
      ce  = 0._r8
c
      return
      end subroutine coef
c
c-----------------------------------------------------------------------
c
      subroutine bndyc(kbdy,xory,alfa,beta,gama,gbdy) 
c
      implicit none
c
      integer, intent(in)    :: kbdy
      real(r8),    intent(in)    :: xory
      real(r8),    intent(inout) :: alfa,beta,gama,gbdy
c
c     kbdy == 1 : x=xa
c     kbdy == 2 : x=xb
c     kbdy == 3 : y=yc
c     kbdy == 4 : y=yd
c
c     For example, when kbdy == 4 (i.e., if y == yd),
c     b.c. has the form alfyd(x)*px+betyd(x)*py+gamyd(x)*pe = gbdyd(x)
c     where x = yorx.   alfa,beta,gama,gbdy corresponding to alfyd(x),
c     betyd(x),gamyd(x),gbdyd(y) must be output.
c
      IF (kbdy == 1) THEN
        alfa = 0.
        beta = 0.
        gama = 0.
        gbdy = 0.
      ELSE IF (kbdy == 2) THEN
        alfa = 0.
        beta = 0.
        gama = 0.
        gbdy = 0.
      ELSE IF (kbdy == 3) THEN
        alfa = 0.
        beta = 0.
        gama = 0.
        gbdy = 0.
      ELSE IF (kbdy == 4) THEN
        alfa = 0.
        beta = 1.
        gama = 0.
        gbdy = 0.
      END IF
c
      return
      end subroutine bndyc
c
      end module mud2crf90
