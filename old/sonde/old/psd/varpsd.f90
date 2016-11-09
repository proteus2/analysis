!2345
    module varpsd

    real, allocatable, dimension(:)   :: normt  
    real, allocatable, dimension(:)   :: zst    ,ztr  
    real, allocatable, dimension(:,:) :: ustprt ,vstprt ,tstprt ,tstbas ,nstbas
    real, allocatable, dimension(:,:) :: utrprt ,vtrprt ,ttrprt ,ttrbas ,ntrbas
    real, allocatable, dimension(:)   :: tnstprt,tntrprt
    real, allocatable, dimension(:)   :: avgstn ,avgtrn

    real, allocatable, dimension(:)   :: wavnst ,wavlst
    real, allocatable, dimension(:)   :: wavntr ,wavltr

    real, allocatable, dimension(:)   :: tstpsd1 ,tstpsdc1 ,ustpsd1 ,vstpsd1
    real, allocatable, dimension(:)   :: tstpsdm1,tstpsdcm1,ustpsdm1,vstpsdm1
    real, allocatable, dimension(:)   :: ttrpsd1 ,          utrpsd1 ,vtrpsd1
    real, allocatable, dimension(:)   :: ttrpsdm1,          utrpsdm1,vtrpsdm1
    real, allocatable, dimension(:)   :: tstsat1 ,ustsat1
    real, allocatable, dimension(:)   :: ttrsat1 ,utrsat1
    real, allocatable, dimension(:)   :: sterat1 ,trerat1
    real, allocatable, dimension(:)   :: tsterr1 ,ttrerr1

    real, allocatable, dimension(:,:) :: tstpsd ,tstpsdc ,ustpsd ,vstpsd
    real, allocatable, dimension(:,:) :: tstpsdm,tstpsdcm,ustpsdm,vstpsdm
    real, allocatable, dimension(:,:) :: ttrpsd ,         utrpsd ,vtrpsd
    real, allocatable, dimension(:,:) :: ttrpsdm,         utrpsdm,vtrpsdm
    real, allocatable, dimension(:,:) :: tstsat ,ustsat 
    real, allocatable, dimension(:,:) :: ttrsat ,utrsat 
    real, allocatable, dimension(:,:) :: sterat ,trerat

    real, allocatable, dimension(:)   :: tstpsd_mean ,tstpsdc_mean,tstpsdm_mean
    real, allocatable, dimension(:)   :: ustpsd_mean ,ustpsdm_mean
    real, allocatable, dimension(:)   :: vstpsd_mean ,vstpsdm_mean
    real, allocatable, dimension(:)   :: tstsat_mean ,ustsat_mean ,sterat_mean
    real, allocatable, dimension(:)   :: ttrpsd_mean ,             ttrpsdm_mean
    real, allocatable, dimension(:)   :: utrpsd_mean ,utrpsdm_mean
    real, allocatable, dimension(:)   :: vtrpsd_mean ,vtrpsdm_mean
    real, allocatable, dimension(:)   :: ttrsat_mean ,utrsat_mean ,trerat_mean


    contains

    subroutine initvarpsd(nt    ,nzst   ,nztr   )

    implicit none

    integer, intent(in) :: nt    ,nzst   ,nztr

    allocate(normt(1:nt))         ; normt(1:nt) = 0.0

    allocate(zst(1:nzst))         ; zst(1:nzst) = 0.0
    allocate(ztr(1:nztr))         ; ztr(1:nztr) = 0.0

    allocate(avgstn(1:nt))        ; avgstn(1:nt) = 0.0
    allocate(avgtrn(1:nt))        ; avgtrn(1:nt) = 0.0

    allocate(ustprt(1:nt,1:nzst)) ; ustprt(1:nt,1:nzst) = 0.0
    allocate(vstprt(1:nt,1:nzst)) ; vstprt(1:nt,1:nzst) = 0.0
    allocate(tstprt(1:nt,1:nzst)) ; tstprt(1:nt,1:nzst) = 0.0
    allocate(tstbas(1:nt,1:nzst)) ; tstbas(1:nt,1:nzst) = 0.0
    allocate(nstbas(1:nt,1:nzst)) ; nstbas(1:nt,1:nzst) = 0.0
    allocate(utrprt(1:nt,1:nztr)) ; utrprt(1:nt,1:nztr) = 0.0
    allocate(vtrprt(1:nt,1:nztr)) ; vtrprt(1:nt,1:nztr) = 0.0
    allocate(ttrprt(1:nt,1:nztr)) ; ttrprt(1:nt,1:nztr) = 0.0
    allocate(ttrbas(1:nt,1:nztr)) ; ttrbas(1:nt,1:nztr) = 0.0
    allocate(ntrbas(1:nt,1:nztr)) ; ntrbas(1:nt,1:nztr) = 0.0

    allocate(tnstprt(1:nzst))     ; tnstprt(1:nzst) = 0.0
    allocate(tntrprt(1:nztr))     ; tntrprt(1:nztr) = 0.0

    allocate(wavnst(1:nzst/2+1))      ; wavnst(1:nzst/2+1) = 0.0
    allocate(wavntr(1:nztr/2+1))      ; wavntr(1:nztr/2+1) = 0.0
    allocate(wavlst(1:nzst/2  ))      ; wavlst(1:nzst/2)   = 0.0
    allocate(wavltr(1:nztr/2  ))      ; wavltr(1:nztr/2)   = 0.0

    allocate(tstpsd1  (1:nzst/2+1))    ; tstpsd1  (1:nzst/2+1) = 0.0
    allocate(tstpsdc1 (1:nzst/2+1))    ; tstpsdc1 (1:nzst/2+1) = 0.0
    allocate(ustpsd1  (1:nzst/2+1))    ; ustpsd1  (1:nzst/2+1) = 0.0
    allocate(vstpsd1  (1:nzst/2+1))    ; vstpsd1  (1:nzst/2+1) = 0.0
    allocate(tstpsdm1 (1:nzst/2+1))    ; tstpsdm1 (1:nzst/2+1) = 0.0
    allocate(tstpsdcm1(1:nzst/2+1))    ; tstpsdcm1(1:nzst/2+1) = 0.0
    allocate(ustpsdm1 (1:nzst/2+1))    ; ustpsdm1 (1:nzst/2+1) = 0.0
    allocate(vstpsdm1 (1:nzst/2+1))    ; vstpsdm1 (1:nzst/2+1) = 0.0
    allocate(ttrpsd1  (1:nztr/2+1))    ; ttrpsd1  (1:nztr/2+1) = 0.0
    allocate(utrpsd1  (1:nztr/2+1))    ; utrpsd1  (1:nztr/2+1) = 0.0
    allocate(vtrpsd1  (1:nztr/2+1))    ; vtrpsd1  (1:nztr/2+1) = 0.0
    allocate(ttrpsdm1 (1:nztr/2+1))    ; ttrpsdm1 (1:nztr/2+1) = 0.0
    allocate(utrpsdm1 (1:nztr/2+1))    ; utrpsdm1 (1:nztr/2+1) = 0.0
    allocate(vtrpsdm1 (1:nztr/2+1))    ; vtrpsdm1 (1:nztr/2+1) = 0.0
    allocate(tstsat1  (1:nzst/2+1))    ; tstsat1  (1:nzst/2+1) = 0.0
    allocate(ustsat1  (1:nzst/2+1))    ; ustsat1  (1:nzst/2+1) = 0.0
    allocate(ttrsat1  (1:nztr/2+1))    ; ttrsat1  (1:nztr/2+1) = 0.0
    allocate(utrsat1  (1:nztr/2+1))    ; utrsat1  (1:nztr/2+1) = 0.0
    allocate(sterat1  (1:nzst/2+1))    ; sterat1  (1:nzst/2+1) = 0.0
    allocate(trerat1  (1:nzst/2+1))    ; trerat1  (1:nzst/2+1) = 0.0
    allocate(tsterr1  (1:nzst/2+1))    ; tsterr1  (1:nzst/2+1) = 0.0
    allocate(ttrerr1  (1:nztr/2+1))    ; ttrerr1  (1:nztr/2+1) = 0.0

    allocate(tstpsd  (1:nzst/2+1,1:nt))  ; tstpsd  (1:nzst/2+1,1:nt) = 0.0
    allocate(tstpsdc (1:nzst/2+1,1:nt))  ; tstpsdc (1:nzst/2+1,1:nt) = 0.0
    allocate(ustpsd  (1:nzst/2+1,1:nt))  ; ustpsd  (1:nzst/2+1,1:nt) = 0.0
    allocate(vstpsd  (1:nzst/2+1,1:nt))  ; vstpsd  (1:nzst/2+1,1:nt) = 0.0
    allocate(tstpsdm (1:nzst/2+1,1:nt))  ; tstpsdm (1:nzst/2+1,1:nt) = 0.0
    allocate(tstpsdcm(1:nzst/2+1,1:nt))  ; tstpsdcm(1:nzst/2+1,1:nt) = 0.0
    allocate(ustpsdm (1:nzst/2+1,1:nt))  ; ustpsdm (1:nzst/2+1,1:nt) = 0.0
    allocate(vstpsdm (1:nzst/2+1,1:nt))  ; vstpsdm (1:nzst/2+1,1:nt) = 0.0
    allocate(ttrpsd  (1:nztr/2+1,1:nt))  ; ttrpsd  (1:nztr/2+1,1:nt) = 0.0
    allocate(utrpsd  (1:nztr/2+1,1:nt))  ; utrpsd  (1:nztr/2+1,1:nt) = 0.0
    allocate(vtrpsd  (1:nztr/2+1,1:nt))  ; vtrpsd  (1:nztr/2+1,1:nt) = 0.0
    allocate(ttrpsdm (1:nztr/2+1,1:nt))  ; ttrpsdm (1:nztr/2+1,1:nt) = 0.0
    allocate(utrpsdm (1:nztr/2+1,1:nt))  ; utrpsdm (1:nztr/2+1,1:nt) = 0.0
    allocate(vtrpsdm (1:nztr/2+1,1:nt))  ; vtrpsdm (1:nztr/2+1,1:nt) = 0.0
    allocate(tstsat  (1:nzst/2+1,1:nt))  ; tstsat  (1:nzst/2+1,1:nt) = 0.0
    allocate(ustsat  (1:nzst/2+1,1:nt))  ; ustsat  (1:nzst/2+1,1:nt) = 0.0
    allocate(ttrsat  (1:nztr/2+1,1:nt))  ; ttrsat  (1:nztr/2+1,1:nt) = 0.0
    allocate(utrsat  (1:nztr/2+1,1:nt))  ; utrsat  (1:nztr/2+1,1:nt) = 0.0
    allocate(sterat  (1:nzst/2+1,1:nt))  ; sterat  (1:nzst/2+1,1:nt) = 0.0
    allocate(trerat  (1:nztr/2+1,1:nt))  ; trerat  (1:nztr/2+1,1:nt) = 0.0

    allocate( tstpsd_mean(1:nzst/2+1))   ;  tstpsd_mean(1:nzst/2+1) = 0.0
    allocate(tstpsdc_mean(1:nzst/2+1))   ; tstpsdc_mean(1:nzst/2+1) = 0.0
    allocate(tstpsdm_mean(1:nzst/2+1))   ; tstpsdm_mean(1:nzst/2+1) = 0.0
    allocate( ustpsd_mean(1:nzst/2+1))   ;  ustpsd_mean(1:nzst/2+1) = 0.0
    allocate(ustpsdm_mean(1:nzst/2+1))   ; ustpsdm_mean(1:nzst/2+1) = 0.0
    allocate( vstpsd_mean(1:nzst/2+1))   ;  vstpsd_mean(1:nzst/2+1) = 0.0
    allocate(vstpsdm_mean(1:nzst/2+1))   ; vstpsdm_mean(1:nzst/2+1) = 0.0
    allocate( tstsat_mean(1:nzst/2+1))   ;  tstsat_mean(1:nzst/2+1) = 0.0
    allocate( ustsat_mean(1:nzst/2+1))   ;  ustsat_mean(1:nzst/2+1) = 0.0
    allocate( sterat_mean(1:nzst/2+1))   ;  sterat_mean(1:nzst/2+1) = 0.0
    allocate( ttrpsd_mean(1:nztr/2+1))   ;  ttrpsd_mean(1:nztr/2+1) = 0.0
    allocate(ttrpsdm_mean(1:nztr/2+1))   ; ttrpsdm_mean(1:nztr/2+1) = 0.0
    allocate( utrpsd_mean(1:nztr/2+1))   ;  utrpsd_mean(1:nztr/2+1) = 0.0
    allocate(utrpsdm_mean(1:nztr/2+1))   ; utrpsdm_mean(1:nztr/2+1) = 0.0
    allocate( vtrpsd_mean(1:nztr/2+1))   ;  vtrpsd_mean(1:nztr/2+1) = 0.0
    allocate(vtrpsdm_mean(1:nztr/2+1))   ; vtrpsdm_mean(1:nztr/2+1) = 0.0
    allocate( ttrsat_mean(1:nztr/2+1))   ;  ttrsat_mean(1:nztr/2+1) = 0.0
    allocate( utrsat_mean(1:nztr/2+1))   ;  utrsat_mean(1:nztr/2+1) = 0.0
    allocate( trerat_mean(1:nztr/2+1))   ;  trerat_mean(1:nztr/2+1) = 0.0

    return
    end subroutine initvarpsd

    end module varpsd

