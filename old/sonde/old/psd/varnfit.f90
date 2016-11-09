!2345
    module varnfit 

    integer                                       :: nma   ,info    ,lwa
    integer                                       :: nregst,nregtr
    integer, allocatable, dimension(:)            :: iwa   ,ipvt
    double precision                              :: tol
    double precision, allocatable, dimension(:)   :: na
    double precision, allocatable, dimension(:)   :: wa

    double precision, allocatable, dimension(:)   :: uxst   ,uyst
    double precision, allocatable, dimension(:)   :: uxtr   ,uytr
    double precision, allocatable, dimension(:)   :: dfitst ,dfittr
    double precision, allocatable, dimension(:)   :: djacst ,djactr

    real, allocatable, dimension(:)   :: mfitst ,mfittr
    real, allocatable, dimension(:,:) :: fitst  ,fittr
    real, allocatable, dimension(:)   :: f0st   ,mstarst ,tslopest
    real, allocatable, dimension(:)   :: f0tr   ,mstartr ,tslopetr
    real, allocatable, dimension(:)   :: f0stm  ,mstarstm,tslopestm
    real, allocatable, dimension(:)   :: f0trm  ,mstartrm,tslopetrm

    contains

!---------------------------------------------------------------------------------
!
!   INITVARNFIT
!
!---------------------------------------------------------------------------------

    subroutine initvarnfit(ndat  ,nrst  ,nrtr  )

!---------------------------------------------------------------------------------
!   nma    : The number of parameters in a nonlinear fitting function
!   nnca   : 2*nma
!   nrst   : The number of points in the vertical in the stratosphere
!   nrtr   : The number of points in the vertical in the troposphere
!   nreg   : The number of points to be used in the nonlinear fitting process
!---------------------------------------------------------------------------------

    integer, intent(in) :: ndat  ,nrst  ,nrtr

    allocate(na (1:nma))             ; na (1:nma) = 0.d0
    allocate(ipvt(1:nma))            ; ipvt(1:nma) = 0
    allocate(iwa(1:nma))             ; iwa(1:nma) = 0
    allocate(wa (1:lwa))             ; wa(1:lwa)  = 0.d0

    allocate(uxst(1:nregst))         ; uxst(1:nregst) = 0.d+00
    allocate(uyst(1:nregst))         ; uyst(1:nregst) = 0.d+00
    allocate(uxtr(1:nregtr))         ; uxtr(1:nregtr) = 0.d+00
    allocate(uytr(1:nregtr))         ; uytr(1:nregtr) = 0.d+00
    allocate(dfitst(1:nregst))       ; dfitst(1:nregst) = 0.d+00
    allocate(dfittr(1:nregtr))       ; dfittr(1:nregtr) = 0.d+00
    allocate(djacst(1:nregst))       ; djacst(1:nregst) = 0.d+00
    allocate(djactr(1:nregtr))       ; djactr(1:nregtr) = 0.d+00

    allocate(mfitst(1:nrst/2+1))         ; mfitst(1:nrst/2+1) = 0.0
    allocate(mfittr(1:nrtr/2+1))         ; mfittr(1:nrtr/2+1) = 0.0
    allocate(fitst(1:nrst/2+1,1:ndat))   ; fitst(1:nrst/2+1,1:ndat) = 0.0
    allocate(fittr(1:nrtr/2+1,1:ndat))   ; fittr(1:nrtr/2+1,1:ndat) = 0.0

    allocate(f0st(1:ndat))           ; f0st(1:ndat) = 0.0
    allocate(f0tr(1:ndat))           ; f0tr(1:ndat) = 0.0
    allocate(mstarst(1:ndat))        ; mstarst(1:ndat) = 0.0
    allocate(mstartr(1:ndat))        ; mstartr(1:ndat) = 0.0
    allocate(tslopest(1:ndat))       ; tslopest(1:ndat) = 0.0
    allocate(tslopetr(1:ndat))       ; tslopetr(1:ndat) = 0.0
    allocate(f0stm(1:1))             ; f0stm(1) = 0.0
    allocate(mstarstm(1:1))          ; mstarstm(1) = 0.0
    allocate(tslopestm(1:1))         ; tslopestm(1) = 0.0
    allocate(f0trm(1:1))             ; f0trm(1) = 0.0
    allocate(mstartrm(1:1))          ; mstartrm(1) = 0.0
    allocate(tslopetrm(1:1))         ; tslopetrm(1) = 0.0

    return
    end subroutine initvarnfit

    end module varnfit
 
