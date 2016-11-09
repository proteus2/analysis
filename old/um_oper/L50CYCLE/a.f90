
    program residualsph

!==================================================================================
!
!   program residualsph
!
!   PURPOSE:
!
!   To calculate numerically residual mean meridional circulation.
!
!   In-Sun Song
!   Laboratory for Mesoscale Dynamics
!   Department of Atmospheric Sciences, Yonsei University, Seoul, Korea.
!
!   Jun 24, 2003
!   First written
!
!   Dec 21, 2004
!   Reviewed
!
!   Jan 15, 2005
!   Mud2crf90 is used as to utilize more general quasi-geostrophic
!   set of equations. This set of equations is equivalent to that used in
!   Holton (1975), Boyd (1976), and Garcia and Solomon (1983). 
!   Bottom boundary condition (formerly, w* = 0 at 1000 hPa) is set to w* 
!   at 100 hPa calculated using the steady-state formulation in Haynes
!   in Haynes et al. [1991, see their (2.7)].
!
!   APR  4, 2005
!   Bottom boundary can also be set to 1000 hPa. In this case, extrapolated
!   wave drag below maximum height of orography at a certain latitude is
!   set to zero.
!
!==================================================================================
!
    use regrid
    use mud2crf90
!
    implicit none
!
    include 'netcdf.inc'
!
    integer, parameter :: ny = 481, nz = 21, nt = 5
    integer, parameter :: iixp = 2, jjyq = 3          ! for MUD2CR
    integer, parameter :: iiex = 9, jjey = 6          ! for MUD2CR
    integer, parameter :: nny = iixp*(2**(iiex-1))+1  ! for MUD2CR
    integer, parameter :: nnz = jjyq*(2**(jjey-1))+1  ! for MUD2CR

!------------------------------------------------------------------------------------
!   Some physical constants
!------------------------------------------------------------------------------------

    real, parameter :: rd = 287.05
    real, parameter :: cp = 1005.0
    real, parameter :: grav = 9.80665
    real, parameter :: kappa = rd/cp
    real, parameter :: p0 = 100000.0
    real, parameter :: hscal = 7000.0
    real, parameter :: hscalp = 7500.0
    real, parameter :: arad = 6371229.0
    real, parameter :: pi = 3.141592!65358979323846
    real, parameter :: omega = 2.*pi/86400.

!------------------------------------------------------------------------------------
!   Grids and data needed for calculation
!------------------------------------------------------------------------------------

    real, dimension(ny)       :: lat   ,y     ,maxh
    real, dimension(nz)       :: z
    real, dimension(nt)       :: mon
    real, dimension(ny,nz,nt) :: zu    ,zt 
    real, dimension(ny,nt)    :: bnep  ,bng   ,bno   ,bnsg
    real, dimension(ny,nz,nt) :: epd   ,gwd   ,gwo   ,gwc
  
    real, dimension(nny)        :: latm  ,ym
    real, dimension(nnz)        :: zm    
    real, dimension(nnz,nt)     :: t0
    real, dimension(nnz,nt)     :: n0sq
    real, dimension(nny,nnz,nt) :: t1
    real, dimension(nny,nnz,nt) :: s
    real, dimension(nny,nnz,nt) :: zum
    real, dimension(nny,nnz,nt) :: ztm
    real, dimension(nny,nnz,nt) :: dudym
    real, dimension(nny,nnz,nt) :: dudzm
    real, dimension(nny,nnz,nt) :: dtdym
    real, dimension(nny,nt)     :: bnepm ,bngm  ,bnom  ,bnsgm
    real, dimension(nny,nnz,nt) :: epdm  ,gwdm  ,gwom  ,gwcm
    real, dimension(nny,nnz,nt) :: chiep ,chigw ,chigo ,chigwc
    real, dimension(nny,nnz,nt) :: vsepd ,vsgwd ,vsgwo ,vsgwc
    real, dimension(nny,nnz,nt) :: wsepd ,wsgwd ,wsgwo ,wsgwc
    real, dimension(nny,nnz,nt) :: dtvepd,dtvgwd,dtvgwo,dtvgwc
    real, dimension(nny,nnz,nt) :: dtwepd,dtwgwd,dtwgwo,dtwgwc

    integer :: n

    call readdata
    call setgrid
    call dnwardctrl
    call dataintpol
    call refatm
    do n=1,nt
    call elliptic(nny,nnz,n0sq(:,n),s(:,:,n),dtdym(:,:,n),zum(:,:,n), &
                  dudzm(:,:,n),dudym(:,:,n),                          &
                  bnepm(:,n),epdm(:,:,n),chiep(:,:,n),vsepd(:,:,n),   &
                  wsepd(:,:,n),dtvepd(:,:,n),dtwepd(:,:,n))
    call elliptic(nny,nnz,n0sq(:,n),s(:,:,n),dtdym(:,:,n),zum(:,:,n), &
                  dudzm(:,:,n),dudym(:,:,n),  &
                  bngm(:,n) ,gwdm(:,:,n),chigw(:,:,n),vsgwd(:,:,n),   &
                  wsgwd(:,:,n),dtvgwd(:,:,n),dtwgwd(:,:,n))
    call elliptic(nny,nnz,n0sq(:,n),s(:,:,n),dtdym(:,:,n),zum(:,:,n), &
                  dudzm(:,:,n),dudym(:,:,n),  &
                  bngm(:,n) ,gwom(:,:,n),chigo(:,:,n),vsgwo(:,:,n),   &
                  wsgwo(:,:,n),dtvgwo(:,:,n),dtwgwo(:,:,n))
    call elliptic(nny,nnz,n0sq(:,n),s(:,:,n),dtdym(:,:,n),zum(:,:,n), &
                  dudzm(:,:,n),dudym(:,:,n),  &
                  bnsgm(:,n),gwcm(:,:,n),chigwc(:,:,n),vsgwc(:,:,n),  &
                  wsgwc(:,:,n),dtvgwc(:,:,n),dtwgwc(:,:,n))
    end do
    call dump

    contains

!------------------------------------------------------------------------------

    subroutine readdata 
  
    character(len=100) :: rfn0  ,rfn1  ,rfn2  ,rfn3  ,rfn4  ,rfn5

    integer :: j     ,k     ,n     ,nn    ,temi
    integer :: istat ,ncid0 ,ncid1 ,ncid2 ,ncid3 ,ncid4 ,ncid5
    integer :: oroid ,latid ,tid   ,monid ,zid   ,zuid  ,ztid  
    integer :: epdid ,gwdid ,gwoid ,gwcid ,bldid

    real, dimension(ny,nz,11,62) :: zu00  ,zt00
    real, dimension(ny,nz,nt,62) :: epd00 ,gwd00 ,gwo00 ,gwc00, bld00

    write(rfn1,'(A)') '/data8/kyh/UM_OPER/anal/L50CYCLE/GWDC/pp/gwdc.uvtz_ypt.2010.01.nc'
    write(rfn2,'(A)') '/data8/kyh/UM_OPER/anal/L50CYCLE/GWDC/pp/gwdc.tem_ypft.201001.nc'
    write(rfn3,'(A)') '/data8/kyh/UM_OPER/anal/L50CYCLE/GWDC/pp/gwdc.zgwd_ypft.2010.01.nc'

    istat = nf_open(trim(rfn1),nf_nowrite,ncid1)
    istat = nf_open(trim(rfn2),nf_nowrite,ncid2)
    istat = nf_open(trim(rfn3),nf_nowrite,ncid3)

    istat = nf_inq_varid(ncid2,'lat'   ,latid)
    istat = nf_inq_varid(ncid2,'zp'    ,zid  )
    istat = nf_inq_varid(ncid2,'fcst'  ,monid)
    istat = nf_inq_varid(ncid2,'t'     ,tid  )
    istat = nf_inq_varid(ncid1,'u'     ,zuid )
    istat = nf_inq_varid(ncid1,'temp'  ,ztid )
    istat = nf_inq_varid(ncid2,'epd'   ,epdid)
    istat = nf_inq_varid(ncid3,'fu_ussp',gwdid)
    istat = nf_inq_varid(ncid3,'fu_gwdo',gwoid)
    istat = nf_inq_varid(ncid3,'fu_gwdc',gwcid)
    istat = nf_inq_varid(ncid3,'fu_bldo',bldid)

    istat = nf_get_var_real(ncid2,latid,lat  )
    istat = nf_get_var_real(ncid2,zid  ,z    )
    istat = nf_get_var_real(ncid2,monid,mon  )
    istat = nf_get_var_real(ncid1,zuid ,zu00 )
    istat = nf_get_var_real(ncid1,ztid ,zt00 )
    istat = nf_get_var_real(ncid2,epdid,epd00)
    istat = nf_get_var_real(ncid3,gwdid,gwd00)
    istat = nf_get_var_real(ncid3,gwoid,gwo00)
    istat = nf_get_var_real(ncid3,gwcid,gwc00)
    istat = nf_get_var_real(ncid3,bldid,bld00)

    istat = nf_close(ncid1)
    istat = nf_close(ncid2)
    istat = nf_close(ncid3)

    zu = 0.  ;  zt = 0.
    do nn=1, 62
    do n=1, 5
      zu(:,:,n) = zu(:,:,n) + 0.5*zu00(:,:,n*2,nn) + &
                              0.25*(zu00(:,:,n*2-1,nn)+zu00(:,:,n*2+1,nn))
      zt(:,:,n) = zt(:,:,n) + 0.5*zt00(:,:,n*2,nn) + &
                              0.25*(zt00(:,:,n*2-1,nn)+zt00(:,:,n*2+1,nn))
    enddo
    enddo
    zu = zu / 62.
    zt = zt / 62.

    epd = 0.  ;  gwd = 0.  ;  gwo = 0.  ;  gwc = 0.
    temi = 0
    do nn=1, 62
      epd(2:ny-1,:,:) = epd(2:ny-1,:,:) + epd00(2:ny-1,:,:,nn)
      if ( gwc00(5,5,1,nn) < 0.9e+20 ) then
        temi = temi + 1
        gwd(:,:,:) = gwd(:,:,:) + gwd00(:,:,:,nn)
        gwo(:,:,:) = gwo(:,:,:) + gwo00(:,:,:,nn) + bld00(:,:,:,nn)
        gwc(:,:,:) = gwc(:,:,:) + gwc00(:,:,:,nn)
      end if
    enddo
    epd = epd / 62. / 86400.
    gwd = gwd / temi / 86400.
    gwo = gwo / temi / 86400.
    gwc = gwc / temi / 86400.

    return
    end subroutine readdata

!-----------------------------------------------------------------------

    subroutine setgrid

    integer :: j     ,k

    do j=1,ny
      y(j) = arad*lat(j)*pi/180.
    end do

    do j=1,nny
      latm(j) = -90.0 + float(j-1)*180.0/float(nny-1)
      ym(j) = arad*latm(j)*pi/180.
    end do
    do k=1,nnz
      zm(k) = 0.0 + float(k-1)*48354./float(nnz-1)
    end do

!gridcheck
    print *,'original grid'
    print *,z
    print *,'grid for mudpack'
    print *,zm(1), zm(nnz), zm(2)-zm(1)

    return
    end subroutine setgrid

!-----------------------------------------------------------------------

    subroutine dnwardctrl

    integer :: j     ,k     ,n
    integer ::  j15n, j15s

    real, dimension(ny) :: phi
    real, dimension(nz) :: rho
    real :: fac,dz
    real :: fac1,fac2
    real, parameter ::  latbdy = 15.

    do j=1,ny
      phi(j) = lat(j)*pi/180.
    end do
    do k=1,nz
      rho(k) = p0/(grav*hscal)*exp(-z(k)/hscal)
    end do
!
!   The formulation used here is based on Haynes et al's steady-state
!   quasi-geostrophic formulation. Therefore near the equator,
!   this formulation may not be valid.
!
    bnep = 0.
    bng  = 0.
    bno  = 0.
    bnsg = 0.

    do n=1,nt

!yh  : f - 1/cos(phi(j)) * d(u(j,k,n)*cos(phi(j)))/dy
    do k=2,nz    ! k = 1 is 1000 hPa.
      dz = z(k)-z(k-1)
      do j=1,ny
        fac = (rho(k-1)+rho(k))*0.5*cos(phi(j))/(2.*omega*sin(phi(j)))
        bnep(j,n) = bnep(j,n) + fac*(epd(j,k-1,n)+epd(j,k,n))*0.5*dz
        bng (j,n) = bng (j,n) + fac*(gwd(j,k-1,n)+gwd(j,k,n))*0.5*dz
        bno (j,n) = bno (j,n) + fac*(gwo(j,k-1,n)+gwo(j,k,n))*0.5*dz
        bnsg(j,n) = bnsg(j,n) + fac*(gwc(j,k-1,n)+gwc(j,k,n))*0.5*dz
      end do
    end do

    do j=1,ny
      bnep(j,n) = (-1./rho(1))*bnep(j,n)
      bng (j,n) = (-1./rho(1))*bng (j,n)
      bno (j,n) = (-1./rho(1))*bno (j,n)
      bnsg(j,n) = (-1./rho(1))*bnsg(j,n)
    end do

    do j=1, ny
      if (lat(j) >= latbdy) then
        j15n = j
        EXIT
      end if
    enddo
    do j=ny, 1, -1
      if (lat(j) <= -latbdy) then
        j15s = j
        EXIT
      end if
    enddo
    call regrid1d(6,(/y(j15s-2),y(j15s-1),y(j15s),y(j15n),y(j15n+1),y(j15n+2)/),&
            (/bnep(j15s-2,n),bnep(j15s-1,n),bnep(j15s  ,n),   &
              bnep(j15n  ,n),bnep(j15n+1,n),bnep(j15n+2,n)/), &
            j15n-j15s-1,y(j15s+1:j15n-1),3,0., bnep(j15s+1:j15n-1,n))
    call regrid1d(6,(/y(j15s-2),y(j15s-1),y(j15s),y(j15n),y(j15n+1),y(j15n+2)/),&
            (/bng (j15s-2,n),bng (j15s-1,n),bng (j15s  ,n),   &
              bng (j15n  ,n),bng (j15n+1,n),bng (j15n+2,n)/), &
            j15n-j15s-1,y(j15s+1:j15n-1),3,0., bng (j15s+1:j15n-1,n))
    call regrid1d(6,(/y(j15s-2),y(j15s-1),y(j15s),y(j15n),y(j15n+1),y(j15n+2)/),&
            (/bno (j15s-2,n),bno (j15s-1,n),bno (j15s  ,n),   &
              bno (j15n  ,n),bno (j15n+1,n),bno (j15n+2,n)/), &
            j15n-j15s-1,y(j15s+1:j15n-1),3,0., bno (j15s+1:j15n-1,n))
    call regrid1d(6,(/y(j15s-2),y(j15s-1),y(j15s),y(j15n),y(j15n+1),y(j15n+2)/),&
            (/bnsg(j15s-2,n),bnsg(j15s-1,n),bnsg(j15s  ,n),   &
              bnsg(j15n  ,n),bnsg(j15n+1,n),bnsg(j15n+2,n)/), &
            j15n-j15s-1,y(j15s+1:j15n-1),3,0., bnsg(j15s+1:j15n-1,n))

    end do

    return
    end subroutine dnwardctrl

!-----------------------------------------------------------------------

    subroutine dataintpol

    integer :: j     ,k     ,n
!intpolcheck
    integer :: istat,ncid
    integer :: latdid,levdid,latid,levid
    integer :: zuid,ztid,epdid,gwdid,gwoid,gwcid

!
!   When ixp = 2 and iex = 6 (or ixp = 3 and iex = 5), latitudinal extrapolation
!   is done only at poles because the largest value of the original latitude
!   is 87.86 and the second largest value of the new grid is 87.18 (86.25)
!
    do n=1,nt
      call regrid1d(ny,y,bnep(:,n),nny,ym,1,0., bnepm(:,n))
      call regrid1d(ny,y,bng (:,n),nny,ym,1,0., bngm (:,n))
      call regrid1d(ny,y,bno (:,n),nny,ym,1,0., bnom (:,n))
      call regrid1d(ny,y,bnsg(:,n),nny,ym,1,0., bnsgm(:,n))

      call regrid2d(ny,nz,y,z,zu (:,:,n),nny,nnz,ym,zm,(/1,3/),1.e20, zum (:,:,n))
      call regrid2d(ny,nz,y,z,zt (:,:,n),nny,nnz,ym,zm,(/1,3/),1.e20, ztm (:,:,n))
      call regrid2d(ny,nz,y,z,epd(:,:,n),nny,nnz,ym,zm,(/1,3/),0.   , epdm(:,:,n))
      call regrid2d(ny,nz,y,z,gwd(:,:,n),nny,nnz,ym,zm,(/1,3/),0.   , gwdm(:,:,n))
      call regrid2d(ny,nz,y,z,gwo(:,:,n),nny,nnz,ym,zm,(/1,3/),0.   , gwom(:,:,n))
      call regrid2d(ny,nz,y,z,gwc(:,:,n),nny,nnz,ym,zm,(/1,3/),0.   , gwcm(:,:,n))
    end do
!
!   In author's experience, cyy (actually, czz in problem of our interest)
!   becomes frequently negative due to term du/dy. In this case, the
!   differential equation is not elliptic any more.
!   To prevent this difficulty, zu at poles is simply set to that at the
!   closest meridional grid. In case of the meridional temperature gradient,
!   that is not included in coefficients of the equation, and thus special
!   treatments are not necessary, but zt at poles is determined in the same
!   way as in zu.
! 
    do n=1,nt
    do k=1,nnz
      zum(1,  k,n) = 0.0
      zum(nny,k,n) = 0.0
      ztm(1,  k,n) = ztm(2    ,k,n)
      ztm(nny,k,n) = ztm(nny-1,k,n)
    end do
    end do

!intpolcheck
    istat = nf_create('res/intpol_check.nc',nf_clobber,ncid)
    istat = nf_def_dim(ncid,'lat',nny,latdid)
    istat = nf_def_dim(ncid,'z'  ,nnz,levdid)
    istat = nf_def_var(ncid,'lat',nf_real,1,latdid,latid)
    istat = nf_def_var(ncid,'z'  ,nf_real,1,levdid,levid)
    istat = nf_def_var(ncid,'zu' ,nf_real,2,(/latdid,levdid/),zuid)
    istat = nf_def_var(ncid,'zt' ,nf_real,2,(/latdid,levdid/),ztid)
    istat = nf_def_var(ncid,'epd',nf_real,2,(/latdid,levdid/),epdid)
    istat = nf_def_var(ncid,'gwd',nf_real,2,(/latdid,levdid/),gwdid)
    istat = nf_def_var(ncid,'gwo',nf_real,2,(/latdid,levdid/),gwoid)
    istat = nf_enddef(ncid)
    istat = nf_put_var_real(ncid,latid,latm)
    istat = nf_put_var_real(ncid,levid,zm)
    istat = nf_put_var_real(ncid,zuid,zum(:,:,1))
    istat = nf_put_var_real(ncid,ztid,ztm(:,:,1))
    istat = nf_put_var_real(ncid,epdid,epdm(:,:,1))
    istat = nf_put_var_real(ncid,gwdid,gwdm(:,:,1))
    istat = nf_put_var_real(ncid,gwoid,gwom(:,:,1))
    istat = nf_close(ncid)

    return
    end subroutine dataintpol

!-----------------------------------------------------------------------
 
    subroutine refatm

    integer :: j     ,k     ,n

    real, dimension(nnz,nt)     :: tmpz
    real, dimension(nny,nnz,nt) :: tmpz1,dtdz
    real, dimension(nny,nnz,nt) :: tmpy
!bgvalcheck
    integer :: istat,ncid
    integer :: latdid,levdid
    integer :: latid,levid,nsqid,dtdyid,sid,dudyid,dudzid
!
!   Reference temperature profile
!
    do n=1,nt
    do k=1,nnz
      t0(k,n) = 0.0
      do j=1,nny
        t0(k,n) = t0(k,n) + ztm(j,k,n)
      end do
      t0(k,n) = t0(k,n)/float(nny)
    end do
    end do
!
!   Define deviation of temperature from its reference profile.
!
    do n=1,nt
    do k=1,nnz
      do j=1,nny
        t1(j,k,n) = ztm(j,k,n) - t0(k,n)
      end do
    end do
    end do
!
!   Brunt-Vaisala frequency of the reference temperature
! 
    do n=1,nt
    tmpz(1:nnz,n) = 0.0
    do k=1,nnz-1
      tmpz(k,n) = (rd/hscal)*( (t0(k+1,n)-t0(k,n))/(zm(k+1)-zm(k)) +  &
                  kappa*(t0(k+1,n)+t0(k,n))*0.5/hscal )
    end do
    do k=2,nnz-1
      n0sq(k,n) = (tmpz(k,n)+tmpz(k-1,n))*0.5
    end do
    n0sq(1,n)   = n0sq(2,n)
    n0sq(nnz,n) = n0sq(nnz-1,n)
    end do
!
!   Vertical gradient of the temperature deviation from its reference
!
    do n=1,nt
    tmpz1(1:nny,1:nnz,n) = 0.0
    do k=1,nnz-1
      do j=1,nny
        tmpz1(j,k,n) = (t1(j,k+1,n)-t1(j,k,n))/(zm(k+1)-zm(k))
      end do
    end do
    do k=2,nnz-1
      do j=1,nny
        dtdz(j,k,n) = (tmpz1(j,k,n)+tmpz1(j,k-1,n))*0.5
      end do
    end do
    do j=1,nny
      dtdz(j,1,n)   = dtdz(j,2,n)
      dtdz(j,nnz,n) = dtdz(j,nnz-1,n)
    end do
    end do
!
!   Stability term in zonal-mean temperature equation
!
    do n=1,nt
    do k=1,nnz
      do j=1,nny
        s(j,k,n) = (hscal/rd)*n0sq(k,n)+dtdz(j,k,n)
      end do
    end do
    end do
!
!   Meridional temperature gradient
!
    do n=1,nt
    tmpy(1:nny,1:nnz,n) = 0.0
    do k=1,nnz
      do j=1,nny-1
        tmpy(j,k,n) = (ztm(j+1,k,n)-ztm(j,k,n))/(ym(j+1)-ym(j))
      end do
    end do
    do k=1,nnz
      do j=2,nny-1
        dtdym(j,k,n) = (tmpy(j-1,k,n)+tmpy(j,k,n))*0.5
      end do
      dtdym(  1,k,n) = 0.0
      dtdym(nny,k,n) = 0.0
    end do
    end do
!
!   Meridional gradient of zonal-mean zonal wind
!
    do n=1,nt
    tmpy(1:nny,1:nnz,n) = 0.0
    do k=1,nnz
      do j=1,nny-1
        tmpy(j,k,n) = (zum(j+1,k,n)-zum(j,k,n))/(ym(j+1)-ym(j))
      end do
    end do
    do k=1,nnz
      do j=2,nny-1
        dudym(j,k,n) = (tmpy(j-1,k,n)+tmpy(j,k,n))*0.5
      end do
      dudym(  1,k,n) = 0.0
      dudym(nny,k,n) = 0.0
    end do
    end do
!
!   Vertical gradient of zonal-mean zonal wind
! 
    do n=1,nt
    tmpz1(1:nny,1:nnz,n) = 0.0
    do k=1,nnz-1
      do j=1,nny
        tmpz1(j,k,n) = (zum(j,k+1,n)-zum(j,k,n))/(zm(k+1)-zm(k))
      end do
    end do
    do k=2,nnz-1
      do j=1,nny
        dudzm(j,k,n) = (tmpz1(j,k,n)+tmpz1(j,k-1,n))*0.5
      end do
    end do
    do j=1,nny
      dudzm(j,1,n)   = 0.0
      dudzm(j,nnz,n) = 0.0
    end do
    end do

!bgvalcheck
    istat = nf_create('res/bgval_check.nc',nf_clobber,ncid)
    istat = nf_def_dim(ncid,'lat',nny,latdid)
    istat = nf_def_dim(ncid,'z'  ,nnz,levdid)
    istat = nf_def_var(ncid,'lat',nf_real,1,latdid,latid)
    istat = nf_def_var(ncid,'z'  ,nf_real,1,levdid,levid)
    istat = nf_def_var(ncid,'nsq',nf_real,1,levdid,nsqid)
    istat = nf_def_var(ncid,'dtdy',nf_real,2,(/latdid,levdid/),dtdyid)
    istat = nf_def_var(ncid,'stbl',nf_real,2,(/latdid,levdid/),sid)
    istat = nf_def_var(ncid,'dudy',nf_real,2,(/latdid,levdid/),dudyid)
    istat = nf_def_var(ncid,'dudz',nf_real,2,(/latdid,levdid/),dudzid)
    istat = nf_enddef(ncid)
    istat = nf_put_var_real(ncid,latid,latm)
    istat = nf_put_var_real(ncid,levid,zm)
    istat = nf_put_var_real(ncid,nsqid,n0sq(:,1))
    istat = nf_put_var_real(ncid,dtdyid,dtdym(:,:,1))
    istat = nf_put_var_real(ncid,sid,s(:,:,1))
    istat = nf_put_var_real(ncid,dudyid,dudym(:,:,1))
    istat = nf_put_var_real(ncid,dudzid,dudzm(:,:,1))
    istat = nf_close(ncid)
!
    return
    end subroutine refatm
!
!-----------------------------------------------------------------------
! 
    subroutine elliptic(nny1,nnz1,n0sqin,stblin,dtdyin,zuin,dudzin,dudyin,  &
                        bndyin,frcin,chiout,vsout,wsout,dtvout,dtwout)
!
    integer,                    intent(in) :: nny1  ,nnz1
    real, dimension(nnz1),      intent(in) :: n0sqin
    real, dimension(nny1,nnz1), intent(in) :: stblin,dtdyin
    real, dimension(nny1,nnz1), intent(in) :: zuin,dudzin,dudyin
    real, dimension(nny1),      intent(in) :: bndyin
    real, dimension(nny1,nnz1), intent(in) :: frcin
    real, dimension(nny1,nnz1), intent(out) :: chiout,vsout ,wsout
    real, dimension(nny1,nnz1), intent(out) :: dtvout,dtwout
!   
!-----------------------------------------------------------------------
!   Coefficients of elliptic partial differential equation (cxxusr, cxyusr,
!   cyyusr, cxusr, cyusr, ceusr) are declared in module mud2crf90.
!-----------------------------------------------------------------------
!
    real, dimension(nny1,nnz1) :: forcing
    real, dimension(nny1,nnz1) :: tmp   ,tmp1
!
    real(r8) :: fcor,dfdy,tanova,cosf,fsign,dx,dy,eps
    real(r8), dimension(nny1)  :: ymr8
    real(r8), dimension(nnz1)  :: zmr8
    real(r8), dimension(nny1,nnz1) :: frcr8,chir8
!
    integer :: j     ,k
    integer :: j1    ,j2    ,k1    ,k2
    integer :: jj    ,kk    ,l
    integer :: npos
!
    integer  :: ixp,iex,jyq,jey
    integer  :: nxa,nxb,nyc,nyd
    integer  :: iguess,maxcy,method
    integer, dimension(4) :: mgopt
!
!   xusr, yusr, cxxusr, cxyusr, cyyusr, cxusr, cyusr, ceusr are declared in
!   mud2crf90, and they should be allocated here.
!
    allocate(xusr  (1:nny1))
    allocate(yusr  (1:nnz1))
    allocate(cxxusr(1:nny1,1:nnz1))
    allocate(cxyusr(1:nny1,1:nnz1))
    allocate(cyyusr(1:nny1,1:nnz1))
    allocate(cxusr (1:nny1,1:nnz1))
    allocate(cyusr (1:nny1,1:nnz1))
    allocate(ceusr (1:nny1,1:nnz1))
!      
    do j=1,nny1
      xusr(j) = ym(j)
    end do
    do k=1,nnz1
      yusr(k) = zm(k)
    end do 
!
    do k=1,nnz1
      do j=1,nny1
        fcor   = 2._r8*omega*dsin(latm(j)*pi/180._r8)
        if ( fcor == 0.0 ) fcor = 0.2 * 2._r8*omega*dsin(latm(j+1)*pi/180._r8)
        dfdy   = (2._r8*omega/arad)*dcos(latm(j)*pi/180._r8)
        tanova = dtan(latm(j)*pi/180._r8)/arad
        cxxusr(j,k) = (rd/hscal)*stblin(j,k)
        cxyusr(j,k) = 2._r8*fcor*dudzin(j,k)
        cyyusr(j,k) = fcor*(fcor+zuin(j,k)*tanova-dudyin(j,k))
        cxusr (j,k) = (rd/hscal)*stblin(j,k)*tanova - (fcor/hscal)*dudzin(j,k)
        cyusr (j,k) = -(fcor/hscal)*(fcor+zuin(j,k)*tanova-dudyin(j,k)) + &
                       dudzin(j,k)*(2._r8*fcor*tanova+dfdy)
        ceusr (j,k) = 0._r8
      end do
    end do
!
!   Ellipticity check
!
    eps = 1.e-20
    dx = xusr(2)-xusr(1)
    dy = yusr(2)-yusr(1)
    do k=1,nnz1
      do j=1,nny1
!
!   In author's experience, in most cases cxxusr is positive because negative
!   temperature calculated from time-averaged temperature is not frequent in.
!   large-scale models. Therefore, cxx should be made positive in advance.
!
        fcor   = 2._r8*omega*dsin(latm(j)*pi/180._r8)
        if ( fcor == 0.0 ) fcor = 0.2 * 2._r8*omega*dsin(latm(j+1)*pi/180._r8)
        tanova = dtan(latm(j)*pi/180._r8)/arad
        if ( cxxusr(j,k) <= 0.0 .or. cyyusr(j,k) <= 0.0 .or. & 
                    4.*cxxusr(j,k)*cyyusr(j,k) <= cxyusr(j,k)**2      ) then
          write(6,'(2(A,1X,F7.2))') 'Ellipticity check failed at lat = ',latm(j),' and z (km) = ',zm(k)/1.e3
          if ( cxxusr(j,k) <= 0.0 ) write(6,*) 'Stability should be checked'
          if ( cyyusr(j,k) <= 0.0 ) write(6,*) 'du/dy should be checked'
          cxxusr(j,k) = n0sqin(k)
          cyyusr(j,k) = fcor*fcor
          cxyusr(j,k) = eps
          cxusr (j,k) = n0sqin(k)*tanova
          cyusr (j,k) = -fcor*fcor/hscal
          ceusr (j,k) = 0.0 
        end if
!       This warning is related to MUD2CR error code -4.
!       MUD2CR change cxx into max(cxx,0.5*|cx|*dx) to work around.'
!       However, in most problems of our interest, cxx is more important than
!       cx. Therefore in this program, cx will be modified.
!       cxusr(j,k) = sign(1.0,cxusr(j,k))*2.*abs(cxxusr(j,k))/dx*eps
        if ( abs(cxusr(j,k))*dx > 2.*abs(cxxusr(j,k)) ) then
          cxusr(j,k) = stblin(j,k)*tanova
        end if
        if ( abs(cyusr(j,k))*dy > 2.*abs(cyyusr(j,k)) ) then
          cyusr(j,k) = -fcor*fcor/hscal
        end if
      end do
    end do

!-----------------------------------------------------------------------
!   Determine forcing of the Elliptic equation
!   Initialize mass stream function
!-----------------------------------------------------------------------

    tmp(1:nny1,1:nnz1) = 0.0

    do k=1,nnz1-1
      do j=1,nny1
        tmp(j,k) = (2.0*omega*sin(latm(j)*pi/180.))*  &
                    (frcin(j,k+1)-frcin(j,k))/(zm(k+1)-zm(k))*  &
                     cos(latm(j)*pi/180.)
      end do
    end do

    do k=2,nnz1-1
      do j=1,nny1
        forcing(j,k) = (tmp(j,k)+tmp(j,k-1))*0.5
      end do
    end do

    do j=1,nny1
      forcing(j,1) = forcing(j,2)
      forcing(j,nnz1) = forcing(j,nnz1-1)
    end do

    do k=1,nnz1
      do j=1,nny1
        forcing(j,k) = forcing(j,k)
      end do
    end do

    do k=2,nnz1
      do j=1,nny1
        chiout(j,k) = 0.0
      end do 
    end do
    do j=1,nny1
      chiout(j,1) = bndyin(j)
    end do

    do k=1,nnz1
      zmr8(k) = zm(k)*1._r8
    end do
    do j=1,nny1
      ymr8(j) = ym(j)*1._r8
    end do
    do k=1,nnz1
      do j=1,nny1
        frcr8(j,k) = forcing(j,k)*1._r8
        chir8(j,k) = chiout(j,k)*1._r8
      end do
    end do

!-----------------------------------------------------------------------
!   Solve Elliptic equation using MULTIGRID PACKAGE
!-----------------------------------------------------------------------

    ixp = iixp; jyq = jjyq; iex = iiex; jey = jjey
    nxa = 1; nxb = 1; nyc = 1; nyd = 2
    iguess = 0; maxcy = 1; method = 0
    mgopt = (/2,2,1,3/)

    call sol2cr(ixp,jyq,iex,jey,nny1,nnz1,  &
                nxa,nxb,nyc,nyd,iguess,maxcy,method,  &
                mgopt,ymr8,zmr8,frcr8,chir8,50)

    do k=1,nnz1
      do j=1,nny1
        if ( abs(chir8(j,k)) < 1.e-20 ) then
          chir8(j,k) = 0.0
        end if 
        chiout(j,k) = chir8(j,k)
      end do
    end do

    deallocate(xusr)
    deallocate(yusr)
    deallocate(cxxusr)
    deallocate(cxyusr)
    deallocate(cyyusr)
    deallocate(cxusr)
    deallocate(cyusr)
    deallocate(ceusr)
      
!-----------------------------------------------------------------------
!   meridional flow
!-----------------------------------------------------------------------

    tmp (1:nny1,1:nnz1) = 0.0
    tmp1(1:nny1,1:nnz1) = 0.0
    vsout(1:nny1,1:nnz1) = 0.0
    dtvout(1:nny1,1:nnz1) = 0.0

    do k=1,nnz1-1 
    do j=1,nny1
      tmp(j,k) = (chiout(j,k+1)-chiout(j,k))/(zm(k+1)-zm(k)) 
    end do
    end do
    
    do k=2,nnz1-1
    do j=1,nny1
      tmp1(j,k) = (tmp(j,k)+tmp(j,k-1))*0.5 
    end do
    end do

    do k=2,nnz1-1
    do j=1,nny1
      vsout (j,k) = (-1./cos(latm(j)*pi/180.))*(tmp1(j,k)-chiout(j,k)/hscal)
      dtvout(j,k) = -vsout(j,k)*dtdyin(j,k)
    end do
    end do

    do j=1,nny1
      vsout(j,1) = 1.e+20
      vsout(j,nnz1) = 1.e+20
      dtvout(j,1) = 1.e+20
      dtvout(j,nnz1) = 1.e+20
    end do

    tmp (1:nny1,1:nnz1) = 0.0
    tmp1(1:nny1,1:nnz1) = 0.0
    wsout(1:nny1,1:nnz1) = 0.0
    dtwout(1:nny1,1:nnz1) = 0.0

    do k=1,nnz1
    do j=1,nny1-1
      tmp(j,k) = (chiout(j+1,k)-chiout(j,k))/(ym(j+1)-ym(j))
    end do
    end do

    do k=1,nnz1
    do j=2,nny1-1
      tmp1(j,k) = (tmp(j,k)+tmp(j-1,k))*0.5
    end do
    end do

    do k=1,nnz1
    do j=2,nny1-1
      wsout(j,k) = (1./cos(latm(j)*pi/180.))*tmp1(j,k)
      dtwout(j,k) = -wsout(j,k)*stblin(j,k)
    end do
    end do

    do k=1,nnz1
      wsout(1,k) = 1.e+20
      wsout(nny1,k) = 1.e+20
      dtwout(1,k) = 1.e+20
      dtwout(nny1,k) = 1.e+20
    end do

    return
    end subroutine elliptic

!-----------------------------------------------------------------------

    subroutine dump

    character(len=100) :: wfn1
  
    integer :: istat ,ncid
    integer :: latdid,zdid  ,mondid
    integer :: latid ,zid   ,monid
    integer :: frcid ,chiid ,dudtid,dtdtid   ! for PALMER case
    integer :: epdid1,gwdid1,gwoid1,gwcid1
    integer :: chiid1,chiid2,chiid3,chiid4
    integer :: vsid1 ,vsid2 ,vsid3 ,vsid4
    integer :: wsid1 ,wsid2 ,wsid3 ,wsid4
    integer :: dtvid1,dtvid2,dtvid3,dtvid4
    integer :: dtwid1,dtwid2,dtwid3,dtwid4

    write(wfn1,'(A)') 'res/sgwdclm2.nc'

    istat = nf_create(trim(wfn1),nf_clobber,ncid)
    istat = nf_def_dim(ncid,'lat',nny,latdid)
    istat = nf_def_dim(ncid,'z'  ,nnz,zdid  )
    istat = nf_def_dim(ncid,'mon',nt ,mondid)
    istat = nf_def_var(ncid,'lat'    ,nf_real,1,latdid,latid)
    istat = nf_def_var(ncid,'z'      ,nf_real,1,zdid  ,zid  )
    istat = nf_def_var(ncid,'mon'    ,nf_real,1,mondid,monid)
    istat = nf_def_var(ncid,'epd'    ,nf_real,3,(/latdid,zdid,mondid/),epdid1)
    istat = nf_def_var(ncid,'gwd'    ,nf_real,3,(/latdid,zdid,mondid/),gwdid1)
    istat = nf_def_var(ncid,'gwo'    ,nf_real,3,(/latdid,zdid,mondid/),gwoid1)
    istat = nf_def_var(ncid,'chiepd' ,nf_real,3,(/latdid,zdid,mondid/),chiid1)
    istat = nf_def_var(ncid,'chigwd' ,nf_real,3,(/latdid,zdid,mondid/),chiid2)
    istat = nf_def_var(ncid,'chigwo' ,nf_real,3,(/latdid,zdid,mondid/),chiid3)
    istat = nf_def_var(ncid,'vsepd'  ,nf_real,3,(/latdid,zdid,mondid/),vsid1)
    istat = nf_def_var(ncid,'vsgwd'  ,nf_real,3,(/latdid,zdid,mondid/),vsid2)
    istat = nf_def_var(ncid,'vsgwo'  ,nf_real,3,(/latdid,zdid,mondid/),vsid3)
    istat = nf_def_var(ncid,'wsepd'  ,nf_real,3,(/latdid,zdid,mondid/),wsid1)
    istat = nf_def_var(ncid,'wsgwd'  ,nf_real,3,(/latdid,zdid,mondid/),wsid2)
    istat = nf_def_var(ncid,'wsgwo'  ,nf_real,3,(/latdid,zdid,mondid/),wsid3)
    istat = nf_def_var(ncid,'dtvepd' ,nf_real,3,(/latdid,zdid,mondid/),dtvid1)
    istat = nf_def_var(ncid,'dtvgwd' ,nf_real,3,(/latdid,zdid,mondid/),dtvid2)
    istat = nf_def_var(ncid,'dtvgwo' ,nf_real,3,(/latdid,zdid,mondid/),dtvid3)
    istat = nf_def_var(ncid,'dtwepd' ,nf_real,3,(/latdid,zdid,mondid/),dtwid1)
    istat = nf_def_var(ncid,'dtwgwd' ,nf_real,3,(/latdid,zdid,mondid/),dtwid2)
    istat = nf_def_var(ncid,'dtwgwo' ,nf_real,3,(/latdid,zdid,mondid/),dtwid3)
    istat = nf_put_att_real(ncid,epdid1,'_FillValue',nf_real,1,1.e+20)
    istat = nf_put_att_real(ncid,gwdid1,'_FillValue',nf_real,1,1.e+20)
    istat = nf_put_att_real(ncid,gwoid1,'_FillValue',nf_real,1,1.e+20)
    istat = nf_put_att_real(ncid,vsid1 ,'_FillValue',nf_real,1,1.e+20)
    istat = nf_put_att_real(ncid,vsid2 ,'_FillValue',nf_real,1,1.e+20)
    istat = nf_put_att_real(ncid,vsid3 ,'_FillValue',nf_real,1,1.e+20)
    istat = nf_put_att_real(ncid,wsid1 ,'_FillValue',nf_real,1,1.e+20)
    istat = nf_put_att_real(ncid,wsid2 ,'_FillValue',nf_real,1,1.e+20)
    istat = nf_put_att_real(ncid,wsid3 ,'_FillValue',nf_real,1,1.e+20)
    istat = nf_put_att_real(ncid,dtvid1,'_FillValue',nf_real,1,1.e+20)
    istat = nf_put_att_real(ncid,dtvid2,'_FillValue',nf_real,1,1.e+20)
    istat = nf_put_att_real(ncid,dtvid3,'_FillValue',nf_real,1,1.e+20)
    istat = nf_put_att_real(ncid,dtwid1,'_FillValue',nf_real,1,1.e+20)
    istat = nf_put_att_real(ncid,dtwid2,'_FillValue',nf_real,1,1.e+20)
    istat = nf_put_att_real(ncid,dtwid3,'_FillValue',nf_real,1,1.e+20)
    istat = nf_def_var(ncid,'gwc'   ,nf_real,3,(/latdid,zdid,mondid/),gwcid1)
    istat = nf_def_var(ncid,'chigwc' ,nf_real,3,(/latdid,zdid,mondid/),chiid4)
    istat = nf_def_var(ncid,'vsgwc'  ,nf_real,3,(/latdid,zdid,mondid/),vsid4)
    istat = nf_def_var(ncid,'wsgwc'  ,nf_real,3,(/latdid,zdid,mondid/),wsid4)
    istat = nf_def_var(ncid,'dtvgwc' ,nf_real,3,(/latdid,zdid,mondid/),dtvid4)
    istat = nf_def_var(ncid,'dtwgwc' ,nf_real,3,(/latdid,zdid,mondid/),dtwid4)
    istat = nf_put_att_real(ncid,gwcid1,'_FillValue',nf_real,1,1.e+20)
    istat = nf_put_att_real(ncid,vsid4  ,'_FillValue',nf_real,1,1.e+20)
    istat = nf_put_att_real(ncid,wsid4  ,'_FillValue',nf_real,1,1.e+20)
    istat = nf_put_att_real(ncid,dtvid4 ,'_FillValue',nf_real,1,1.e+20)
    istat = nf_put_att_real(ncid,dtwid4 ,'_FillValue',nf_real,1,1.e+20)
    istat = nf_enddef(ncid)
    istat = nf_put_var_real(ncid,latid ,latm )
    istat = nf_put_var_real(ncid,zid   ,zm)
    istat = nf_put_var_real(ncid,monid ,mon)
    istat = nf_put_var_real(ncid,epdid1,epdm)
    istat = nf_put_var_real(ncid,gwdid1,gwdm)
    istat = nf_put_var_real(ncid,gwoid1,gwom)
    istat = nf_put_var_real(ncid,chiid1,chiep)
    istat = nf_put_var_real(ncid,chiid2,chigw)
    istat = nf_put_var_real(ncid,chiid3,chigo)
    istat = nf_put_var_real(ncid,vsid1,vsepd)
    istat = nf_put_var_real(ncid,vsid2,vsgwd)
    istat = nf_put_var_real(ncid,vsid3,vsgwo)
    istat = nf_put_var_real(ncid,wsid1,wsepd)
    istat = nf_put_var_real(ncid,wsid2,wsgwd)
    istat = nf_put_var_real(ncid,wsid3,wsgwo)
    istat = nf_put_var_real(ncid,dtvid1,dtvepd)
    istat = nf_put_var_real(ncid,dtvid2,dtvgwd)
    istat = nf_put_var_real(ncid,dtvid3,dtvgwo)
    istat = nf_put_var_real(ncid,dtwid1,dtwepd)
    istat = nf_put_var_real(ncid,dtwid2,dtwgwd)
    istat = nf_put_var_real(ncid,dtwid3,dtwgwo)
    istat = nf_put_var_real(ncid,gwcid1,gwcm)
    istat = nf_put_var_real(ncid,chiid4,chigwc)
    istat = nf_put_var_real(ncid,vsid4 ,vsgwc)
    istat = nf_put_var_real(ncid,wsid4 ,wsgwc)
    istat = nf_put_var_real(ncid,dtvid4,dtvgwc)
    istat = nf_put_var_real(ncid,dtwid4,dtwgwc)
    istat = nf_close(ncid)

    return
    end subroutine dump

    end program residualsph
