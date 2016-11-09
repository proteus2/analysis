#undef  CONTROL1
#undef  CONTROL2
#undef  SGWDCLM0
#undef  SGWDCLM1
#define SGWDCLM2

    module var

    use grid
    use interpolation

    implicit none

    integer, dimension(nm) :: daymon = &
             (/31,28,31,30,31,30,31,31,30,31,30,31/)
 
    real, allocatable, dimension(:,:,:,:) :: udy   ,vdy   ,tdy   ,omgdy

    real, dimension(ny )        :: lat
    real, dimension(nz )        :: lev   ,logz
    real, dimension(nym)        :: latm
    real, dimension(nzm)        :: levm  ,logzm
    real, dimension(nx ,ny ,nz) :: u     ,v     ,t     ,omga
    real, dimension(nx ,ny ,nz) :: utm   ,vtm   ,ttm   ,omgatm
    real, dimension(nym,nzm)    :: vs    ,ws
    real, dimension(nym,nzm)    :: cor   ,cur   ,advyu ,advzu
    real, dimension(nym,nzm)    :: epy   ,epyst ,epytr
    real, dimension(nym,nzm)    :: epz   ,epzst ,epztr
    real, dimension(nym,nzm)    :: epd   ,epdst ,epdtr
    real, dimension(nym,nzm)    :: advyt ,advzt
    real, dimension(nym,nzm)    :: ft    ,ftst  ,fttr
    real, dimension(ny ,nz )    :: vsn   ,wsn
    real, dimension(ny ,nz )    :: corn  ,curn  ,advyun,advzun
    real, dimension(ny ,nz )    :: epyn  ,epystn,epytrn
    real, dimension(ny ,nz )    :: epzn  ,epzstn,epztrn
    real, dimension(ny ,nz )    :: epdn  ,epdstn,epdtrn
    real, dimension(ny ,nz )    :: advytn,advztn
    real, dimension(ny ,nz )    :: ftn   ,ftstn ,fttrn
    real, dimension(ny ,nz )    :: hfy   ,hfyst ,hfytr 
    real, dimension(ny ,nz )    :: hfz   ,hfzst ,hfztr 
    real, dimension(ny ,nz )    :: mfy   ,mfyst ,mfytr 
    real, dimension(ny ,nz )    :: mfz   ,mfzst ,mfztr 
 
    real, dimension(nym,nzm)    :: vs1   ,ws1
    real, dimension(nym,nzm)    :: cor1  ,cur1  ,advyu1,advzu1
    real, dimension(nym,nzm)    :: epy1  ,epyst1,epytr1
    real, dimension(nym,nzm)    :: epz1  ,epzst1,epztr1
    real, dimension(nym,nzm)    :: epd1  ,epdst1,epdtr1
    real, dimension(nym,nzm)    :: advyt1,advzt1
    real, dimension(nym,nzm)    :: ft1   ,ftst1 ,fttr1
    real, dimension(ny ,nz )    :: hfy1  ,hfyst1,hfytr1
    real, dimension(ny ,nz )    :: hfz1  ,hfzst1,hfztr1
    real, dimension(ny ,nz )    :: mfy1  ,mfyst1,mfytr1
    real, dimension(ny ,nz )    :: mfz1  ,mfzst1,mfztr1

    contains

!------------------------------------------------------------------------------------

    subroutine initvar
 
    vs    (1:nym,1:nzm) = 0.0; vsn   (1:ny ,1:nz ) = 0.0
    ws    (1:nym,1:nzm) = 0.0; wsn   (1:ny ,1:nz ) = 0.0
    cor   (1:nym,1:nzm) = 0.0; corn  (1:ny ,1:nz ) = 0.0
    cur   (1:nym,1:nzm) = 0.0; curn  (1:ny ,1:nz ) = 0.0
    advyu (1:nym,1:nzm) = 0.0; advyun(1:ny ,1:nz ) = 0.0
    advzu (1:nym,1:nzm) = 0.0; advzun(1:ny ,1:nz ) = 0.0
    epy   (1:nym,1:nzm) = 0.0; epyn  (1:ny ,1:nz ) = 0.0
    epyst (1:nym,1:nzm) = 0.0; epystn(1:ny ,1:nz ) = 0.0
    epytr (1:nym,1:nzm) = 0.0; epytrn(1:ny ,1:nz ) = 0.0
    epz   (1:nym,1:nzm) = 0.0; epzn  (1:ny ,1:nz ) = 0.0
    epzst (1:nym,1:nzm) = 0.0; epzstn(1:ny ,1:nz ) = 0.0
    epztr (1:nym,1:nzm) = 0.0; epztrn(1:ny ,1:nz ) = 0.0
    epd   (1:nym,1:nzm) = 0.0; epdn  (1:ny ,1:nz ) = 0.0
    epdst (1:nym,1:nzm) = 0.0; epdstn(1:ny ,1:nz ) = 0.0
    epdtr (1:nym,1:nzm) = 0.0; epdtrn(1:ny ,1:nz ) = 0.0
    advyt (1:nym,1:nzm) = 0.0; advytn(1:ny ,1:nz ) = 0.0
    advzt (1:nym,1:nzm) = 0.0; advztn(1:ny ,1:nz ) = 0.0
    ft    (1:nym,1:nzm) = 0.0; ftn   (1:ny ,1:nz ) = 0.0
    ftst  (1:nym,1:nzm) = 0.0; ftstn (1:ny ,1:nz ) = 0.0
    fttr  (1:nym,1:nzm) = 0.0; fttrn (1:ny ,1:nz ) = 0.0
    hfy   (1:ny ,1:nz ) = 0.0
    hfyst (1:ny ,1:nz ) = 0.0
    hfytr (1:ny ,1:nz ) = 0.0
    hfz   (1:ny ,1:nz ) = 0.0
    hfzst (1:ny ,1:nz ) = 0.0
    hfztr (1:ny ,1:nz ) = 0.0
    mfy   (1:ny ,1:nz ) = 0.0
    mfyst (1:ny ,1:nz ) = 0.0
    mfytr (1:ny ,1:nz ) = 0.0
    mfz   (1:ny ,1:nz ) = 0.0
    mfzst (1:ny ,1:nz ) = 0.0
    mfztr (1:ny ,1:nz ) = 0.0

    return
    end subroutine initvar

    subroutine initvartmp

    vs1   (1:nym,1:nzm) = 0.0
    ws1   (1:nym,1:nzm) = 0.0
    cor1  (1:nym,1:nzm) = 0.0
    cur1  (1:nym,1:nzm) = 0.0
    advyu1(1:nym,1:nzm) = 0.0
    advzu1(1:nym,1:nzm) = 0.0
    epy1  (1:nym,1:nzm) = 0.0
    epyst1(1:nym,1:nzm) = 0.0
    epytr1(1:nym,1:nzm) = 0.0
    epz1  (1:nym,1:nzm) = 0.0
    epzst1(1:nym,1:nzm) = 0.0
    epztr1(1:nym,1:nzm) = 0.0
    epd1  (1:nym,1:nzm) = 0.0
    epdst1(1:nym,1:nzm) = 0.0
    epdtr1(1:nym,1:nzm) = 0.0
    advyt1(1:nym,1:nzm) = 0.0
    advzt1(1:nym,1:nzm) = 0.0
    ft1   (1:nym,1:nzm) = 0.0
    ftst1 (1:nym,1:nzm) = 0.0
    fttr1 (1:nym,1:nzm) = 0.0
    hfy1  (1:ny ,1:nz ) = 0.0
    hfyst1(1:ny ,1:nz ) = 0.0
    hfytr1(1:ny ,1:nz ) = 0.0
    hfz1  (1:ny ,1:nz ) = 0.0
    hfzst1(1:ny ,1:nz ) = 0.0
    hfztr1(1:ny ,1:nz ) = 0.0
    mfy1  (1:ny ,1:nz ) = 0.0
    mfyst1(1:ny ,1:nz ) = 0.0
    mfytr1(1:ny ,1:nz ) = 0.0
    mfz1  (1:ny ,1:nz ) = 0.0
    mfzst1(1:ny ,1:nz ) = 0.0
    mfztr1(1:ny ,1:nz ) = 0.0

    return
    end subroutine initvartmp

    subroutine addresult

    vs     = vs     + vs1
    ws     = ws     + ws1
    cor    = cor    + cor1
    cur    = cur    + cur1
    advyu  = advyu  + advyu1
    advzu  = advzu  + advzu1
    epy    = epy    + epy1
    epyst  = epyst  + epyst1
    epytr  = epytr  + epytr1
    epz    = epz    + epz1
    epzst  = epzst  + epzst1
    epztr  = epztr  + epztr1
    epd    = epd    + epd1
    epdst  = epdst  + epdst1
    epdtr  = epdtr  + epdtr1
    advyt  = advyt  + advyt1
    advzt  = advzt  + advzt1
    ft     = ft     + ft1
    ftst   = ftst   + ftst1
    fttr   = fttr   + fttr1
    hfy    = hfy    + hfy1
    hfyst  = hfyst  + hfyst1
    hfytr  = hfytr  + hfytr1
    hfz    = hfz    + hfz1
    hfzst  = hfzst  + hfzst1
    hfztr  = hfztr  + hfztr1
    mfy    = mfy    + mfy1
    mfyst  = mfyst  + mfyst1
    mfytr  = mfytr  + mfytr1
    mfz    = mfz    + mfz1
    mfzst  = mfzst  + mfzst1
    mfztr  = mfztr  + mfztr1

    return
    end subroutine addresult

    subroutine avgresult(im)

    integer, intent(in) :: im

    vs     = vs   /float(daymon(im))
    ws     = ws   /float(daymon(im))
    cor    = cor  /float(daymon(im))
    cur    = cur  /float(daymon(im))
    advyu  = advyu/float(daymon(im))
    advzu  = advzu/float(daymon(im))
    epy    = epy  /float(daymon(im))
    epyst  = epyst/float(daymon(im))
    epytr  = epytr/float(daymon(im))
    epz    = epz  /float(daymon(im))
    epzst  = epzst/float(daymon(im))
    epztr  = epztr/float(daymon(im))
    epd    = epd  /float(daymon(im))
    epdst  = epdst/float(daymon(im))
    epdtr  = epdtr/float(daymon(im))
    advyt  = advyt/float(daymon(im))
    advzt  = advzt/float(daymon(im))
    ft     = ft   /float(daymon(im))
    ftst   = ftst /float(daymon(im))
    fttr   = fttr /float(daymon(im))
    hfy    = hfy  /float(daymon(im))
    hfyst  = hfyst/float(daymon(im))
    hfytr  = hfytr/float(daymon(im))
    hfz    = hfz  /float(daymon(im))
    hfzst  = hfzst/float(daymon(im))
    hfztr  = hfztr/float(daymon(im))
    mfy    = mfy  /float(daymon(im))
    mfyst  = mfyst/float(daymon(im))
    mfytr  = mfytr/float(daymon(im))
    mfz    = mfz  /float(daymon(im))
    mfzst  = mfzst/float(daymon(im))
    mfztr  = mfztr/float(daymon(im))

    return
    end subroutine avgresult

    subroutine makeintpol

    call int2d(nym   ,nzm   ,latm  ,logzm ,vs    ,ny    ,nz    ,lat   ,logz  ,vsn   ,1,1,1.e+20)
    call int2d(nym   ,nzm   ,latm  ,logzm ,ws    ,ny    ,nz    ,lat   ,logz  ,wsn   ,1,1,1.e+20)
    call int2d(nym   ,nzm   ,latm  ,logzm ,cor   ,ny    ,nz    ,lat   ,logz  ,corn  ,1,1,1.e+20)
    call int2d(nym   ,nzm   ,latm  ,logzm ,cur   ,ny    ,nz    ,lat   ,logz  ,curn  ,1,1,1.e+20)
    call int2d(nym   ,nzm   ,latm  ,logzm ,advyu ,ny    ,nz    ,lat   ,logz  ,advyun,1,1,1.e+20)
    call int2d(nym   ,nzm   ,latm  ,logzm ,advzu ,ny    ,nz    ,lat   ,logz  ,advzun,1,1,1.e+20)
    call int2d(nym   ,nzm   ,latm  ,logzm ,epy   ,ny    ,nz    ,lat   ,logz  ,epyn  ,1,1,1.e+20)
    call int2d(nym   ,nzm   ,latm  ,logzm ,epyst ,ny    ,nz    ,lat   ,logz  ,epystn,1,1,1.e+20)
    call int2d(nym   ,nzm   ,latm  ,logzm ,epytr ,ny    ,nz    ,lat   ,logz  ,epytrn,1,1,1.e+20)
    call int2d(nym   ,nzm   ,latm  ,logzm ,epz   ,ny    ,nz    ,lat   ,logz  ,epzn  ,1,1,1.e+20)
    call int2d(nym   ,nzm   ,latm  ,logzm ,epzst ,ny    ,nz    ,lat   ,logz  ,epzstn,1,1,1.e+20)
    call int2d(nym   ,nzm   ,latm  ,logzm ,epztr ,ny    ,nz    ,lat   ,logz  ,epztrn,1,1,1.e+20)
    call int2d(nym   ,nzm   ,latm  ,logzm ,epd   ,ny    ,nz    ,lat   ,logz  ,epdn  ,1,1,1.e+20)
    call int2d(nym   ,nzm   ,latm  ,logzm ,epdst ,ny    ,nz    ,lat   ,logz  ,epdstn,1,1,1.e+20)
    call int2d(nym   ,nzm   ,latm  ,logzm ,epdtr ,ny    ,nz    ,lat   ,logz  ,epdtrn,1,1,1.e+20)
    call int2d(nym   ,nzm   ,latm  ,logzm ,advyt ,ny    ,nz    ,lat   ,logz  ,advytn,1,1,1.e+20)
    call int2d(nym   ,nzm   ,latm  ,logzm ,advzt ,ny    ,nz    ,lat   ,logz  ,advztn,1,1,1.e+20)
    call int2d(nym   ,nzm   ,latm  ,logzm ,ft    ,ny    ,nz    ,lat   ,logz  ,ftn   ,1,1,1.e+20)
    call int2d(nym   ,nzm   ,latm  ,logzm ,ftst  ,ny    ,nz    ,lat   ,logz  ,ftstn ,1,1,1.e+20)
    call int2d(nym   ,nzm   ,latm  ,logzm ,fttr  ,ny    ,nz    ,lat   ,logz  ,fttrn ,1,1,1.e+20)
    
    return
    end subroutine makeintpol

    subroutine read_timemean(iyear ,imonth)

    implicit none

    include 'netcdf.inc'

    integer, intent(in) :: iyear ,imonth

    character(len=100) :: namehead,fname
    integer :: idy
    integer :: istat ,ncid
    integer :: latid ,levid
    integer :: uid   ,vid   ,tid   ,omgid

#if ( defined CONTROL1 )
    namehead = '/export24/sis/waccm-control1/daily/control1.waccm.day.'
#elif ( defined CONTROL2 )
    namehead = '/export28/sis/waccm-control2/daily/control2.waccm.day.'
#elif ( defined SGWDCLM0 )
    namehead = '/export22/sis/waccm-sgwdclm0/daily/sgwdclm0.waccm.day.'
#elif ( defined SGWDCLM1 )
    namehead = '/export23/sis/waccm-sgwdclm1/daily/sgwdclm1.waccm.day.'
#elif ( defined SGWDCLM2 )
    namehead = '/export22/sis/waccm-sgwdclm2/daily/sgwdclm2.waccm.day.'
#endif

    write(fname,'(A,I4.4,A,I2.2,A)') trim(namehead),iyear,'-',imonth,'.nc'
    write(6,*) trim(fname)

    allocate(udy  (1:nx,1:ny,1:nz,1:daymon(imonth))) 
    allocate(vdy  (1:nx,1:ny,1:nz,1:daymon(imonth))) 
    allocate(tdy  (1:nx,1:ny,1:nz,1:daymon(imonth))) 
    allocate(omgdy(1:nx,1:ny,1:nz,1:daymon(imonth))) 

    istat = nf_open(trim(fname),nf_nowrite,ncid)
    istat = nf_inq_varid(ncid,'lat'  ,latid)
    istat = nf_inq_varid(ncid,'lev'  ,levid)
    istat = nf_inq_varid(ncid,'U'    ,uid) 
    istat = nf_inq_varid(ncid,'V'    ,vid) 
    istat = nf_inq_varid(ncid,'T'    ,tid) 
    istat = nf_inq_varid(ncid,'OMEGA',omgid)
    istat = nf_get_var_real(ncid,latid,lat   ); print *,istat
    istat = nf_get_var_real(ncid,levid,lev   ); print *,istat
    istat = nf_get_var_real(ncid,uid  ,udy   ); print *,istat
    istat = nf_get_var_real(ncid,vid  ,vdy   ); print *,istat
    istat = nf_get_var_real(ncid,tid  ,tdy   ); print *,istat
    istat = nf_get_var_real(ncid,omgid,omgdy ); print *,istat
    istat = nf_close(ncid)

    utm   (1:nx,1:ny,1:nz) = 0.0
    vtm   (1:nx,1:ny,1:nz) = 0.0
    ttm   (1:nx,1:ny,1:nz) = 0.0
    omgatm(1:nx,1:ny,1:nz) = 0.0

    do idy=1,daymon(imonth)
      utm    = utm    + udy  (1:nx,1:ny,1:nz,idy)
      vtm    = vtm    + vdy  (1:nx,1:ny,1:nz,idy)
      ttm    = ttm    + tdy  (1:nx,1:ny,1:nz,idy)
      omgatm = omgatm + omgdy(1:nx,1:ny,1:nz,idy)
    end do

    utm    = utm   /float(daymon(imonth))
    vtm    = vtm   /float(daymon(imonth))
    ttm    = ttm   /float(daymon(imonth))
    omgatm = omgatm/float(daymon(imonth))
   
    return
    end subroutine read_timemean

    subroutine get_day(iday  )

    integer, intent(in) :: iday

    u   (1:nx,1:ny,1:nz) = udy  (1:nx,1:ny,1:nz,iday)  
    v   (1:nx,1:ny,1:nz) = vdy  (1:nx,1:ny,1:nz,iday)  
    t   (1:nx,1:ny,1:nz) = tdy  (1:nx,1:ny,1:nz,iday)  
    omga(1:nx,1:ny,1:nz) = omgdy(1:nx,1:ny,1:nz,iday)  

    return
    end subroutine get_day
 
    subroutine dump(iyear ,imonth)

    include 'netcdf.inc'

    integer, intent(in) :: iyear ,imonth

    character(len=100) :: namehead,fname
    integer :: istat ,ncid
    integer :: latdid,levdid
    integer :: latid ,levid
    integer :: ivar
    integer,          dimension(32) :: varid 
    character(len=5), dimension(32) :: varnam

    varnam = (/'VS   ','WS   ','COR  ','CUR  ','ADVYU','ADVZU','EPY  ','EPYST','EPYTR','EPZ  ',  &
               'EPZST','EPZTR','EPD  ','EPDST','EPDTR','ADVYT','ADVZT','FT   ','FTST ','FTTR ',  &
               'HFY  ','HFYST','HFYTR','HFZ  ','HFZST','HFZTR','MFY  ','MFYST','MFYTR','MFZ  ',  & 
               'MFZST','MFZTR'/)

#if ( defined CONTROL1 )
    namehead = '/export18/sis/waccm1b_anal/epflux/control1/control1.'
#elif ( defined CONTROL2 )
    namehead = '/export18/sis/waccm1b_anal/epflux/control2/control2.'
#elif ( defined SGWDCLM0 )
    namehead = '/export18/sis/waccm1b_anal/epflux/sgwdclm0/sgwdclm0.'
#elif ( defined SGWDCLM1 )
    namehead = '/export18/sis/waccm1b_anal/epflux/sgwdclm1/sgwdclm1.'
#elif ( defined SGWDCLM2 )
    namehead = '/export18/sis/waccm1b_anal/epflux/sgwdclm2/sgwdclm2.'
#endif

    write(fname,'(A,I4.4,A,I2.2,A)') trim(namehead),iyear,'-',imonth,'.nc'
    write(6,*) trim(fname)

    istat = nf_create(trim(fname),nf_clobber,ncid)
    istat = nf_def_dim(ncid,'lat',ny,latdid)
    istat = nf_def_dim(ncid,'lev',nz,levdid)
    istat = nf_def_var(ncid,'lat',nf_real,1,latdid,latid)
    istat = nf_def_var(ncid,'lev',nf_real,1,levdid,levid)
    do ivar=1,32
      istat = nf_def_var(ncid,varnam(ivar),nf_real,2,(/latdid,levdid/),varid(ivar))
      istat = nf_put_att_real(ncid,varid(ivar),'_FillValue',nf_real,1,1.e+20)
    end do
    istat = nf_enddef(ncid)
    istat = nf_put_var_real(ncid,latid,lat) 
    istat = nf_put_var_real(ncid,levid,lev) 
    istat = nf_put_var_real(ncid,varid( 1),vsn   )
    istat = nf_put_var_real(ncid,varid( 2),wsn   )
    istat = nf_put_var_real(ncid,varid( 3),corn  )
    istat = nf_put_var_real(ncid,varid( 4),curn  )
    istat = nf_put_var_real(ncid,varid( 5),advyun)
    istat = nf_put_var_real(ncid,varid( 6),advzun)
    istat = nf_put_var_real(ncid,varid( 7),epyn  )
    istat = nf_put_var_real(ncid,varid( 8),epystn)
    istat = nf_put_var_real(ncid,varid( 9),epytrn)
    istat = nf_put_var_real(ncid,varid(10),epzn  )
    istat = nf_put_var_real(ncid,varid(11),epzstn)
    istat = nf_put_var_real(ncid,varid(12),epztrn)
    istat = nf_put_var_real(ncid,varid(13),epdn  )
    istat = nf_put_var_real(ncid,varid(14),epdstn)
    istat = nf_put_var_real(ncid,varid(15),epdtrn)
    istat = nf_put_var_real(ncid,varid(16),advytn)
    istat = nf_put_var_real(ncid,varid(17),advztn)
    istat = nf_put_var_real(ncid,varid(18),ftn   )
    istat = nf_put_var_real(ncid,varid(19),ftstn )
    istat = nf_put_var_real(ncid,varid(20),fttrn )
    istat = nf_put_var_real(ncid,varid(21),hfy   )
    istat = nf_put_var_real(ncid,varid(22),hfyst )
    istat = nf_put_var_real(ncid,varid(23),hfytr )
    istat = nf_put_var_real(ncid,varid(24),hfz   )
    istat = nf_put_var_real(ncid,varid(25),hfzst )
    istat = nf_put_var_real(ncid,varid(26),hfztr )
    istat = nf_put_var_real(ncid,varid(27),mfy   )
    istat = nf_put_var_real(ncid,varid(28),mfyst )
    istat = nf_put_var_real(ncid,varid(29),mfytr )
    istat = nf_put_var_real(ncid,varid(30),mfz   )
    istat = nf_put_var_real(ncid,varid(31),mfzst )
    istat = nf_put_var_real(ncid,varid(32),mfztr )
    istat = nf_close(ncid)

    deallocate(udy  ) 
    deallocate(vdy  ) 
    deallocate(tdy  ) 
    deallocate(omgdy) 

    return
    end subroutine dump

    end module var 
