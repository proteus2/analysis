!2345
    module varspec

    real, allocatable, dimension(:)   :: normt
    real, allocatable, dimension(:)   :: zst

    real, allocatable, dimension(:)   :: avgstet,avgstn

    real, allocatable, dimension(:)   :: ustprt1,vstprt1
    real, allocatable, dimension(:)   :: tnstprt,tnsthil
    real, allocatable, dimension(:)   :: ustbas1,vstbas1
    real, allocatable, dimension(:)   :: rstbas1

    real, allocatable, dimension(:)   :: ubst   ,ubdgst

    real, allocatable, dimension(:,:) :: ustprt ,vstprt ,tstprt 
    real, allocatable, dimension(:,:) :: ustbas ,vstbas ,tstbas
    real, allocatable, dimension(:,:) :: rstbas ,nstbas

    real, allocatable, dimension(:)   :: wavmst
    real, allocatable, dimension(:,:) :: rstpsd

    real, allocatable, dimension(:)   :: rstpsd1

    real, allocatable, dimension(:)   :: phist  ,degst  ,degnst
    real, allocatable, dimension(:)   :: w_fst  ,dfst
    real, allocatable, dimension(:)   :: mbst   ,kbst
    real, allocatable, dimension(:)   :: mbast  ,kbast
    real, allocatable, dimension(:)   :: czst   ,cist   ,cixst  ,ciyst
    real, allocatable, dimension(:)   :: cgxst  ,cgyst  ,wgrst  ,cgrst
    real, allocatable, dimension(:)   :: cxst   ,cyst   ,cpxst  ,cpyst

    real, allocatable, dimension(:)   :: angs
    real, allocatable, dimension(:)   :: angst   

    real, allocatable, dimension(:)   :: cclst  ,clst

    real, allocatable, dimension(:)   :: upfrst ,dnfrst

    real  :: mstdir

    real, allocatable, dimension(:)   :: ruwfst1,uwflst1
    real, allocatable, dimension(:)   :: rvwfst1,vwflst1

    real, allocatable, dimension(:,:) :: ruwfst ,uwflst
    real, allocatable, dimension(:,:) :: rvwfst ,vwflst

    real, allocatable, dimension(:)   :: mruwst ,muwfst
    real, allocatable, dimension(:)   :: mrvwst ,mvwfst

    real  :: mstruw, mstuwf
    real  :: mstrvw, mstvwf

    contains

    subroutine initvarspec(ndat  ,nrst  )

    implicit none

    integer, intent(in) :: ndat  ,nrst

    allocate(normt  (1:ndat))        ; normt  (1:ndat) = 0.0
    allocate(zst    (1:nrst))        ; zst    (1:nrst) = 0.0

    allocate(avgstet(1:ndat))        ; avgstet(1:ndat) = 0.0

    allocate(avgstn (1:ndat))        ; avgstn (1:ndat) = 0.0

    allocate(ustprt1(1:nrst))        ; ustprt1(1:nrst) = 0.0
    allocate(vstprt1(1:nrst))        ; vstprt1(1:nrst) = 0.0
    allocate(tnstprt(1:nrst))        ; tnstprt(1:nrst) = 0.0
    allocate(tnsthil(1:nrst))        ; tnsthil(1:nrst) = 0.0
    allocate(ustbas1(1:nrst))        ; ustbas1(1:nrst) = 0.0
    allocate(vstbas1(1:nrst))        ; vstbas1(1:nrst) = 0.0

    allocate(ubst   (1:nrst))        ; ubst(1:nrst) = 0.0
    allocate(ubdgst (1:nrst))        ; ubdgst(1:nrst) = 0.0

    allocate(ustprt (1:ndat,1:nrst)) ; ustprt (1:ndat,1:nrst) = 0.0
    allocate(vstprt (1:ndat,1:nrst)) ; vstprt (1:ndat,1:nrst) = 0.0
    allocate(tstprt (1:ndat,1:nrst)) ; tstprt (1:ndat,1:nrst) = 0.0
    allocate(ustbas (1:ndat,1:nrst)) ; ustbas (1:ndat,1:nrst) = 0.0
    allocate(vstbas (1:ndat,1:nrst)) ; vstbas (1:ndat,1:nrst) = 0.0
    allocate(tstbas (1:ndat,1:nrst)) ; tstbas (1:ndat,1:nrst) = 0.0
    allocate(nstbas (1:ndat,1:nrst)) ; nstbas (1:ndat,1:nrst) = 0.0

    allocate(wavmst (1:nrst))        ; wavmst (1:nrst) = 0.0
    allocate(rstpsd (1:nrst,1:ndat)) ; rstpsd (1:nrst,1:ndat) = 0.0

    allocate(rstpsd1(1:nrst))        ; rstpsd1(1:nrst) = 0.0

    allocate(phist  (1:ndat))        ; phist  (1:ndat) = 0.0
    allocate(degst  (1:ndat))        ; degst  (1:ndat) = 0.0
    allocate(degnst (1:ndat))        ; degnst (1:ndat) = 0.0
    allocate(w_fst  (1:ndat))        ; w_fst  (1:ndat) = 0.0
    allocate(dfst   (1:ndat))        ; dfst   (1:ndat) = 0.0
    allocate(mbst   (1:ndat))        ; mbst   (1:ndat) = 0.0
    allocate(kbst   (1:ndat))        ; kbst   (1:ndat) = 0.0
    allocate(mbast  (1:ndat))        ; mbast  (1:ndat) = 0.0
    allocate(kbast  (1:ndat))        ; kbast  (1:ndat) = 0.0
    allocate(czst   (1:ndat))        ; czst   (1:ndat) = 0.0
    allocate(cist   (1:ndat))        ; cist   (1:ndat) = 0.0
    allocate(cixst  (1:ndat))        ; cixst  (1:ndat) = 0.0
    allocate(ciyst  (1:ndat))        ; ciyst  (1:ndat) = 0.0
    allocate(cgxst  (1:ndat))        ; cgxst  (1:ndat) = 0.0
    allocate(cgyst  (1:ndat))        ; cgyst  (1:ndat) = 0.0
    allocate(wgrst  (1:ndat))        ; wgrst  (1:ndat) = 0.0
    allocate(cgrst  (1:ndat))        ; cgrst  (1:ndat) = 0.0
    allocate(cxst   (1:ndat))        ; cxst   (1:ndat) = 0.0
    allocate(cyst   (1:ndat))        ; cyst   (1:ndat) = 0.0
    allocate(cpxst  (1:ndat))        ; cpxst  (1:ndat) = 0.0
    allocate(cpyst  (1:ndat))        ; cpyst  (1:ndat) = 0.0
  
    allocate(angs   (1:12))          ; angs   (1:12)   = 0.0
    allocate(angst  (1:12))          ; angst  (1:12)   = 0.0

    allocate(cclst  (1:ndat))        ; cclst  (1:ndat) = 0.0
    allocate( clst  (1:ndat))        ;  clst  (1:ndat) = 0.0

    allocate(upfrst (1:ndat))        ; upfrst (1:ndat) = 0.0
    allocate(dnfrst (1:ndat))        ; dnfrst (1:ndat) = 0.0

    allocate(ruwfst1(1:nrst))        ; ruwfst1(1:nrst) = 0.0
    allocate(uwflst1(1:nrst))        ; uwflst1(1:nrst) = 0.0
    allocate(rvwfst1(1:nrst))        ; rvwfst1(1:nrst) = 0.0
    allocate(vwflst1(1:nrst))        ; vwflst1(1:nrst) = 0.0

    allocate(ruwfst (1:ndat,1:nrst)) ; ruwfst (1:ndat,1:nrst) = 0.0
    allocate(uwflst (1:ndat,1:nrst)) ; uwflst (1:ndat,1:nrst) = 0.0
    allocate(rvwfst (1:ndat,1:nrst)) ; rvwfst (1:ndat,1:nrst) = 0.0
    allocate(vwflst (1:ndat,1:nrst)) ; vwflst (1:ndat,1:nrst) = 0.0

    allocate(mruwst (1:ndat))        ; mruwst (1:ndat) = 0.0
    allocate(muwfst (1:ndat))        ; muwfst (1:ndat) = 0.0
    allocate(mrvwst (1:ndat))        ; mrvwst (1:ndat) = 0.0
    allocate(mvwfst (1:ndat))        ; mvwfst (1:ndat) = 0.0

    allocate(rstbas (1:ndat,1:nrst)) ; rstbas (1:ndat,1:nrst) = 0.0

    allocate(rstbas1(1:nrst))        ; rstbas1(1:nrst) = 0.0

    return
   
    end subroutine initvarspec

    end module varspec
