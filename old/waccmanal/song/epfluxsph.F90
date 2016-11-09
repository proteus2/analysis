!2345
#undef DEBUG
    module epfluxsph

    implicit none

    contains

    subroutine epfluxsttr(nx    ,ny    ,nz    ,nym   ,nzm   ,lat   ,lev   ,  &
                          u     ,v     ,t     ,omga  ,utm   ,vtm   ,ttm   ,  &
                          omgatm,latm  ,logzm ,logz  ,vsf   ,wsf   ,corif ,  &
                          curvf ,advyuf,advzuf,epyf  ,epystf,epytrf,epzf  ,  &
                          epzstf,epztrf,epdf  ,epdstf,epdtrf,advytf,advztf,  &
                          ftf   ,ftstf ,fttrf ,hfyf  ,hfystf,hfytrf,hfzf  ,  &
                          hfzstf,hfztrf,mfyf  ,mfystf,mfytrf,mfzf  ,mfzstf,  &
                          mfztrf)  

!-----------------------------------------------------------------------------
!
!   Calculate terms of transformed Eulerian mean equation for zonal-mean zonal
!   wind and potential temperature equations in spherical geometry.
!
!   Author
!   I.-S. Song
!
!   Suppose that certain variables A and B are functions of x, y, z, t,
!   they can be divided into time-averaged and time-varying components:
!
!                         A = [A](x,y,z) + {A}(x,y,z,t),
!                    and  B = [B](x,y,z) + {B}(x,y,z,t),
!
!   where [ ] and { } means the time-averaged and time-varying components
!   of a variable. Perturbations (A' and B') of the variables A and B are 
!   calculated by removing their zonal means, and the perturbations can be
!   represented as the sum of perturbations of their time-averaged and
!   time-varying components of the variables:
!
!                              A' = [A]' + {A}',              
!                         and  B' = [B]' + {B}'.
!
!   Flux term due to the two perturbation variables can be written as
!
!              A'B' = [A]'[B]' + [A]'{B}' + {A}'[B]' + {A}'{B}'.
!
!   The zonal mean of the flux term becomes
!
!         <A'B'> = <[A]'[B]'> + <[A]'{B}'> + <{A}'[B]'> + <{A}'{B}'>.
!
!   If the flux term is averaged over time, the time-averaged flux term
!   becomes
!
!                [<A'B'>] = [<[A]'[B]'>] + [<{A}'{B}'>],
!                         =  <[A]'[B]'>  + [<{A}'{B}'>].
!
!   Therefore, the time-averaged flux terms due to wave perturbations can
!   be represented as the sum of the flux terms due to the stationary
!   (time-averaged) wave perturbations and to the transient (time-varying)
!   wave perturbations.
!
!
!
!   IMPORTANT!!!!!!  (For consistency)
!
!   When we have monthly mean data, we can calculated EP flux due to stationary
!   waves and its divergence. In this case, factors needed to obtain the
!   EP flux (e.g., du/dz or dtheta/dz) are calculated using the monthly
!   mean data. On the other hand, when we calculate EP fluxes due to both the 
!   stationary and transient waves using daily or 6-hourly data, the factors
!   can be calculated using the daily or 6-hourly data. In the two cases,
!   however, EP fluxes by stationary waves are not consistent with each other.
!   If we calculate EP flux due to stationary wave using the monthly-averaged
!   factors, such an inconsistency disappears, and total EP flux becomes
!   rigorously the sum of EP fluxes of stationary and transient waves.
!   For comparison, total EP fluxes will also be calculated using time-varying
!   variables (e.g., epy, epz, epd, fyt, and fzt). The EP fluxes obtained
!   using the monthly-averaged factors are distinguished (e.g., epyst, epytr,
!   epzst, epztr, epdst, epdtr, fytst, fyttr, fztst, and fzttr).
!
!
!   NOV 22, 2004
!   Fixed minor bugs
!   Added documentations
!
!   DEC 15, 2004
!   Added calculation of terms in TEM potential temperature equation.
!   Refined documentaions and source codes.
!
!   USAGE
!
!   INPUT VARIABLES
!
!   nx     : The number of grids in the longitude
!   ny     : The number of grids in the latitude
!   nz     : The number of grids in the vertical
!   nym    : The number of grids in the latitude ( = ny-1 )
!   nzm    : The number of grids in the vertical ( = nz-1 )
!   lat    : Latitude (degree) in an increasing order  (-90 ... 0 ... 90)
!   lev    : Pressure (mb)     in an increasing order  (... 1 ... 100 ... 1000)
!   u      : Zonal wind (m/s)
!   v      : Meridional wind (m/s)
!   t      : Temperature (K)
!   omga   : Pressure velocity (pa/s)
!   utm    : Time-averaged (usually monthly mean) zonal wind (m/s)
!   vtm    : Time-averaged (usually monthly mean) meridional wind (m/s)
!   ttm    : Time-averaged (usually monthly mean) temperature (K)
!   omgatm : Time-averaged (usually monthly mean) pressure velocity (Pa/s) 
!
!   OUTPUT VARIABLES
!
!   latm   : Latitude (degree)
!   logzm  : Log-Pressure height (m)
!   logz   : Log-Pressure height (m)
!   vsf    : Meridional flow in residual mean meridional circulation
!   wsf    : vertical flow in residual mean meridional circulation
!   corif  : Coriolis effect term by residual mean meridional wind 
!   advyuf : Meridional advection of zonal-mean zonal flow due to vstar
!   advzuf : Vertical advection of zonal-mean zonal flow due to wstar
!   epyf   : Meridional component of EP flux
!   epystf : Meridional component of EP flux due to stationary wave
!   epytrf : Meridional component of EP flux due to transient wave
!   epzf   : Vertical component of EP flux
!   epzstf : Vertical component of EP flux due to stationary wave
!   epztrf : Vertical component of EP flux due to transient wave
!   epdf   : EP flux divergence
!   epdstf : EP flux divergence due to stationary wave
!   epdtrf : EP flux divergence due to transient wave
!   advytf : Meridional advection of zonal-mean potential temperature due to vstar
!   advztf : Vertical advection of zonal-mean potential temperature due to wstar
!   ftf    : Tendency of potential temperature 
!   ftstf  : Tendency of potential temperature due to stationary waves
!   fttrf  : Tendency of potential temperature due to transient waves
!   hfyf   : Meridional heatflux
!   hfystf : Meridional heatflux due to stationary wave
!   hfytrf : Meridional heatflux due to transient wave
!   hfzf   : Vertical heatflux
!   hfzstf : Vertical heatflux due to stationary wave
!   hfztrf : Vertical heatflux due to transient wave
!   mfyf   : Meridional momentum flux
!   mfystf : Meridional momentum flux due to stationary wave
!   mfytrf : Meridional momentum flux due to transient wave
!   mfzf   : Vertical momentum flux
!   mfzstf : Vertical momentum flux due to stationary wave
!   mfztrf : Vertical momentum flux due to transient wave
! 
!-----------------------------------------------------------------------------

    implicit none

    integer, parameter :: r8 = selected_real_kind(12)
    integer, parameter :: i8 = selected_int_kind(13)

    real(r8), parameter :: a = 6370000.
    real(r8), parameter :: h = 7000.
!
!   Here, h is set to the constant. However, this constant h is valid below
!   z = 100 km where the ratios of constituents of gases constituting the
!   atmosphere are constant (i.e., the homosphere). In fact, above the lower
!   thermosphere h increases rapidly in the altitude. Therefore the results
!   obtained using this subroutine is in principle valid only below z = 100 km.
!
    real(r8), parameter :: rhos = 1.125
    real(r8), parameter :: ps = 1000.0
    real(r8), parameter :: g = 9.81
    real(r8), parameter :: rd = 287.0
    real(r8), parameter :: cp = 1004.
    real(r8), parameter :: pi = 3.141592
    real(r8), parameter :: radang = 2.*pi/86400.

    integer :: i,j,k,it,jj,kk
    integer :: jsub1,jsub2,ksub1,ksub2,ndtdz
    real    :: dtdzm

    integer,                   intent(in)  :: nx    ,ny    ,nz
    integer,                   intent(in)  :: nym   ,nzm
    real, dimension(ny),       intent(in)  :: lat
    real, dimension(nz),       intent(in)  :: lev
    real, dimension(nx,ny,nz), intent(in)  :: u     ,v     ,t     ,omga
    real, dimension(nx,ny,nz), intent(in)  :: utm   ,vtm   ,ttm   ,omgatm

    real, dimension(nym),      intent(out) :: latm
    real, dimension(nzm),      intent(out) :: logzm
    real, dimension(nz)     ,  intent(out) :: logz
    real, dimension(nym,nzm),  intent(out) :: vsf   ,wsf   ,corif ,curvf
    real, dimension(nym,nzm),  intent(out) :: advyuf,advzuf
    real, dimension(nym,nzm),  intent(out) :: epyf  ,epzf  ,epdf
    real, dimension(nym,nzm),  intent(out) :: epystf,epzstf,epdstf
    real, dimension(nym,nzm),  intent(out) :: epytrf,epztrf,epdtrf
    real, dimension(nym,nzm),  intent(out) :: advytf,advztf
    real, dimension(nym,nzm),  intent(out) :: ftf   ,ftstf ,fttrf
    real, dimension(ny,nz)  ,  intent(out) :: hfyf  ,hfystf,hfytrf
    real, dimension(ny,nz)  ,  intent(out) :: hfzf  ,hfzstf,hfztrf
    real, dimension(ny,nz)  ,  intent(out) :: mfyf  ,mfystf,mfytrf
    real, dimension(ny,nz)  ,  intent(out) :: mfzf  ,mfzstf,mfztrf

    real(r8), dimension(ny)         :: sinph  ,cosph
    real(r8), dimension(nym)        :: sinphm ,cosphm
    real(r8), dimension(nz)         :: rho0
    real(r8), dimension(nzm)        :: rho0m

    real(r8), dimension(nx,ny,nz)   :: thtm  ,wtm
    real(r8), dimension(nx,ny,nz)   :: th    ,w
    real(r8), dimension(ny,nz)      :: zmtmu ,zmtmv ,zmtmth,zmtmw
    real(r8), dimension(ny,nz)      :: zmu   ,zmv   ,zmth  ,zmw

    real(r8), dimension(nx,ny,nz)   :: up    ,vp    ,wp    ,thp 
    real(r8), dimension(nx,ny,nz)   :: upst  ,vpst  ,wpst  ,thpst
    real(r8), dimension(nx,ny,nz)   :: uptr  ,vptr  ,wptr  ,thptr

    real(r8), dimension(ny,nz)      :: dtdz  ,dudz  ,dmtdz ,dmudz
    real(r8), dimension(ny, nz)     :: dtdy  ,dudy  ,dmtdy ,dmudy

    real(r8), dimension(ny,nz)      :: hfy   ,hfz   ,mfy   ,mfz
    real(r8), dimension(ny,nz)      :: hfyst ,hfzst ,mfyst ,mfzst
    real(r8), dimension(ny,nz)      :: hfytr ,hfztr ,mfytr ,mfztr

    real(r8), dimension(nym,nzm)    :: vs    ,ws
    real(r8), dimension(nym,nzm)    :: cori  ,advyu ,advzu ,curv
    real(r8), dimension(nym,nzm)    :: epy   ,epyst ,epytr   
    real(r8), dimension(nym,nzm)    :: epz   ,epzst ,epztr
    real(r8), dimension(nym,nzm)    :: epd   ,epdst ,epdtr
    real(r8), dimension(nym,nzm)    :: advyt ,advzt
    real(r8), dimension(nym,nzm)    :: ft   ,ftst   ,fttr

    real(r8), dimension(ny ,nz )    :: tmp1  ,tmp2  ,tmp3  ,tmp4  ,tmp5  ,tmp6
   
    real(r8) :: fac0,fac1,fac2,fac3,fac4,fac5
 
!------------------------------------------------------------------------
!   Define coordinate 
!------------------------------------------------------------------------

    do j=1,ny
      sinph(j) = dsin(lat(j)*pi/180.0)
      cosph(j) = dcos(lat(j)*pi/180.0)
    end do
    do j=1,nym
      latm(j) = (lat(j)+lat(j+1))*0.5
      sinphm(j) = dsin(latm(j)*pi/180.0)
      cosphm(j) = dcos(latm(j)*pi/180.0)
    end do

#ifdef DEBUG
    write(6,'(6a10)') 'LAT','LATM','SIN','SINM','COS','COSM'
    do j=1,ny
      write(6,'(f10.3,10x,f10.3,10x,f10.3)') lat(j),sinph(j),cosph(j)
      write(6,'(10x,f10.3,10x,f10.3,10x,f10.3)') latm(j),sinphm(j),cosphm(j)
    end do
    pause
#endif 

    do k=1,nz
      logz(k) = h*log(ps/lev(k))       ! logz(1) = the highest level
      rho0(k) = rhos*dexp(-logz(k)/h)
    end do 
    do k=1,nzm
      logzm(k) = (logz(k)+logz(k+1))*0.5
      rho0m(k) = rhos*dexp(-logzm(k)/h)
    end do 

#ifdef DEBUG
    write(6,'(4a10)') 'LOGZ','LOGZM','RHO0','RHO0M'
    do k=1,nz
      write(6,'(f10.3,10x,f10.3)') logz(k), rho0(k)
      write(6,'(10x,f10.3,10x,f10.3)') logzm(k), rho0m(k) 
    end do
    pause
#endif

!------------------------------------------------------------------------
!   Define time-averaged and time-varying theta and w 
!------------------------------------------------------------------------

    do k=1,nz
    do j=1,ny
    do i=1,nx
      thtm(i,j,k) = ttm(i,j,k)*(ps/lev(k))**(rd/cp)
      th  (i,j,k) =   t(i,j,k)*(ps/lev(k))**(rd/cp)
      wtm (i,j,k) = -(h*omgatm(i,j,k))/(lev(k)*100.0)
      w   (i,j,k) = -(h*omga  (i,j,k))/(lev(k)*100.0)
    end do 
    end do
    end do

!------------------------------------------------------------------------
!   Zonal- and time-averaged u, v, w, and theta 
!------------------------------------------------------------------------

    zmtmu (1:ny,1:nz) = 0.d0
    zmtmv (1:ny,1:nz) = 0.d0
    zmtmth(1:ny,1:nz) = 0.d0
    zmtmw (1:ny,1:nz) = 0.d0

    do k=1,nz
    do j=1,ny
      do i=1,nx
        zmtmu (j,k) = zmtmu (j,k) + dble(utm(i,j,k))
        zmtmv (j,k) = zmtmv (j,k) + dble(vtm(i,j,k))
        zmtmth(j,k) = zmtmth(j,k) + thtm(i,j,k)
        zmtmw (j,k) = zmtmw (j,k) + wtm (i,j,k)
      end do
      zmtmu (j,k) = zmtmu (j,k)/float(nx)
      zmtmv (j,k) = zmtmv (j,k)/float(nx)
      zmtmth(j,k) = zmtmth(j,k)/float(nx)
      zmtmw (j,k) = zmtmw (j,k)/float(nx)
    end do
    end do

!------------------------------------------------------------------------
!   Zonally averaged u, v, w, and theta 
!   u, v, w, and theta are function of y, z, and t
!------------------------------------------------------------------------

    zmu (1:ny,1:nz) = 0.d0
    zmv (1:ny,1:nz) = 0.d0
    zmth(1:ny,1:nz) = 0.d0
    zmw (1:ny,1:nz) = 0.d0

    do k=1,nz
    do j=1,ny
      do i=1,nx
        zmu (j,k) = zmu (j,k) + dble(u(i,j,k))
        zmv (j,k) = zmv (j,k) + dble(v(i,j,k))
        zmth(j,k) = zmth(j,k) + th(i,j,k)
        zmw (j,k) = zmw (j,k) + w (i,j,k)
      end do
      zmu (j,k) = zmu (j,k)/float(nx)
      zmv (j,k) = zmv (j,k)/float(nx)
      zmth(j,k) = zmth(j,k)/float(nx)
      zmw (j,k) = zmw (j,k)/float(nx)
    end do
    end do

#if ( defined DEBUG )
    do k = 1, nz
      print '(i4,3(1x,f10.3))', k, logz(k), zmtheta(1,k), zmu(1,k)
    end do
#endif

!------------------------------------------------------------------------
!   Stationary wave perturbation
!------------------------------------------------------------------------

    do k=1,nz
    do j=1,ny
    do i=1,nx
      upst (i,j,k) = dble(utm(i,j,k)) - zmtmu (j,k)
      vpst (i,j,k) = dble(vtm(i,j,k)) - zmtmv (j,k)
      wpst (i,j,k) =      wtm(i,j,k)  - zmtmw (j,k)
      thpst(i,j,k) =     thtm(i,j,k)  - zmtmth(j,k)
    end do
    end do
    end do

!------------------------------------------------------------------------
!   Total wave perturbation
!------------------------------------------------------------------------

    do k=1,nz
    do j=1,ny
    do i=1,nx
      up (i,j,k) = dble(u(i,j,k)) - zmu (j,k)
      vp (i,j,k) = dble(v(i,j,k)) - zmv (j,k)
      wp (i,j,k) =      w(i,j,k)  - zmw (j,k)
      thp(i,j,k) =     th(i,j,k)  - zmth(j,k)
    end do
    end do
    end do

!------------------------------------------------------------------------
!   Transient wave perturbation
!------------------------------------------------------------------------

    uptr  = up  - upst
    vptr  = vp  - vpst
    wptr  = wp  - wpst
    thptr = thp - thpst

!--------------------------------------------------------------------------------
!   Calculate vertical derivative of zonal mean potential temperature.
!   Finally, dtdz and dmtdz are defined at the vertical levels of lev.
!   dtdz  is time-varying   d<theta>/dz.
!   dmtdz is time-invariant d<[theta]>/dz.
!   dtdz and dmtdz are extrapolated at both the top and bottom boundaries.
!--------------------------------------------------------------------------------

    tmp1(1:ny,1:nz) = 0.d0
    tmp2(1:ny,1:nz) = 0.d0
    do k=1,nzm
      do j=1,ny
        fac0 = zmth  (j,k)-zmth  (j,k+1)
        fac1 = zmtmth(j,k)-zmtmth(j,k+1)
        fac2 = logz(k)-logz(k+1)
        tmp1(j,k) = fac0/fac2
        tmp2(j,k) = fac1/fac2
      end do
    end do
    do j=1,ny
      do k=2,nzm
        dtdz (j,k) = (tmp1(j,k)+tmp1(j,k-1))*0.5
        dmtdz(j,k) = (tmp2(j,k)+tmp2(j,k-1))*0.5
      end do
      dtdz (j,1)  = dtdz (j,2) 
      dtdz (j,nz) = dtdz (j,nz-1)
      dmtdz(j,1)  = dmtdz(j,2) 
      dmtdz(j,nz) = dmtdz(j,nz-1)
    end do

    do k=1,nz
    do j=1,ny
      if ( dtdz(j,k) <= 0.0 ) then
        print *,'Invalid vertical gradient of potential temperature'
!       jsub1 = j-2
!       jsub2 = j+2
!       ksub1 = k-2
!       ksub2 = k+2
!       if ( jsub1 < 1  ) jsub1 = 1
!       if ( jsub2 > ny ) jsub2 = ny
!       if ( ksub1 < 1  ) jsub1 = 1
!       if ( ksub2 > nz ) jsub2 = nz
!       dtdzm = 0.0
!       ndtdz = 0
!       do kk=ksub1,ksub2
!       do jj=jsub1,jsub2
!         if ( dtdz(jj,kk) > 0.0 ) then
!           dtdzm = dtdzm + dtdz(jj,kk)
!           ndtdz = ndtdz + 1
!         end if
!       end do
!       end do
!       dtdz(j,k) = dtdzm/float(ndtdz)
        dtdzm = 0.0
        ndtdz = 0
        do jj=1,ny
          if ( dtdz(jj,k) > 0.0 ) then
            dtdzm = dtdzm + dtdz(jj,k)
            ndtdz = ndtdz + 1
          end if
        end do
        dtdz(j,k) = dtdzm/float(ndtdz)
      end if
    end do
    end do 

    do k=1,nz
    do j=1,ny
      if ( dmtdz(j,k) <= 0.0 ) then
        print *,'Invalid vertical gradient of potential temperature'
!       jsub1 = j-2
!       jsub2 = j+2
!       ksub1 = k-2
!       ksub2 = k+2
!       if ( jsub1 < 1  ) jsub1 = 1
!       if ( jsub2 > ny ) jsub2 = ny
!       if ( ksub1 < 1  ) jsub1 = 1
!       if ( ksub2 > nz ) jsub2 = nz
!       dtdzm = 0.0
!       ndtdz = 0
!       do kk=ksub1,ksub2
!       do jj=jsub1,jsub2
!         if ( dmtdz(jj,kk) > 0.0 ) then
!           dtdzm = dtdzm + dmtdz(jj,kk)
!           ndtdz = ndtdz + 1
!         end if
!       end do
!       end do
!       dmtdz(j,k) = dtdzm/float(ndtdz)
        dtdzm = 0.0
        ndtdz = 0
        do jj=1,ny
          if ( dmtdz(jj,k) > 0.0 ) then
            dtdzm = dtdzm + dmtdz(jj,k)
            ndtdz = ndtdz + 1
          end if
        end do
        dmtdz(j,k) = dtdzm/float(ndtdz)
      end if
    end do
    end do

!--------------------------------------------------------------------------------
!   Calculate the vertical derivative of zonal mean zonal wind.
!   Finally, dudz and mdzudz are defined at the vertical levels of lev.
!   dudz  is time-varying   d<u>/dz.
!   dmudz is time-invariant d<[u]>/dz.
!   dudz and dmudz are extrapolated at both the top and bottom boundaries.
!--------------------------------------------------------------------------------

    tmp1(1:ny,1:nz) = 0.d0
    tmp2(1:ny,1:nz) = 0.d0
    do k=1,nzm
      do j=1,ny
        fac0 = zmu  (j,k)-zmu  (j,k+1)
        fac1 = zmtmu(j,k)-zmtmu(j,k+1)
        fac2 = logz(k)-logz(k+1)
        tmp1(j,k) = fac0/fac2
        tmp2(j,k) = fac1/fac2
      end do
    end do
    do j=1,ny
      do k=2,nzm
        dudz (j,k) = (tmp1(j,k)+tmp1(j,k-1))*0.5
        dmudz(j,k) = (tmp2(j,k)+tmp2(j,k-1))*0.5
      end do
      dudz (j, 1) = dudz (j,2) 
      dudz (j,nz) = dudz (j,nz-1) 
      dmudz(j, 1) = dmudz(j,2) 
      dmudz(j,nz) = dmudz(j,nz-1) 
    end do

#if (defined DEBUG)
    write(6,'(4a10)') 'd<TH>/dz', 'd[<TH>]/dz', 'd<U>/dz', 'd[<U>]/dz'
    do k=1,nz
      write(6,'(4f10.5)') dtdz(1,k), dmtdz(1,k), dudz(1,k), dmudz(1,k)
    end do
    pause
#endif

!--------------------------------------------------------------------------------
!   Calculate meridional derivative of zonal-mean potential temperature.
!   Finally, dtdy and dmtdy are defined at the grid of lat.
!   dtdy  is time-variant d<theta>/dphi.
!   dmtdy is time-invariant d<[theta]>/dphi.
!   dtdy and dmtdy are extrapolated at both poles.
!--------------------------------------------------------------------------------

    tmp1(1:ny,1:nz) = 0.d0
    tmp2(1:ny,1:nz) = 0.d0
    do k=1,nz
      do j=1,nym
        fac0 = zmth  (j+1,k)-zmth  (j,k)
        fac1 = zmtmth(j+1,k)-zmtmth(j,k)
        fac2 = (lat(j+1)-lat(j))*pi/180.
        tmp1(j,k) = fac0/fac2
        tmp2(j,k) = fac1/fac2
      end do
    end do
    do k=1,nz
      do j=2,nym
        dtdy (j,k) = (tmp1(j-1,k)+tmp1(j,k))*0.5
        dmtdy(j,k) = (tmp2(j-1,k)+tmp2(j,k))*0.5
      end do
      dtdy (1 ,k) = dtdy (2  ,k)
      dtdy (ny,k) = dtdy (nym,k)
      dmtdy(1 ,k) = dmtdy(2  ,k)
      dmtdy(ny,k) = dmtdy(nym,k)
    end do

!-------------------------------------------------------------------------------
!   Calculate meridional derivative of zonal-mean zonal wind.
!   Finally, dudy and dmudy are defined at the grid of lat.
!   dudy  is time-variant d<u>/dphi.
!   dmudy is time-invariant d<[u]>/dphi.
!   dudy and dmudy are extrapolated at both poles. 
!-------------------------------------------------------------------------------

    tmp1(1:ny,1:nz) = 0.d0
    tmp2(1:ny,1:nz) = 0.d0
    do k=1,nz
      do j=1,nym
        fac0 = zmu  (j+1,k)-zmu  (j,k)
        fac1 = zmtmu(j+1,k)-zmtmu(j,k)
        fac2 = (lat(j+1)-lat(j))*pi/180.
        tmp1(j,k) = fac0/fac2
        tmp2(j,k) = fac1/fac2
      end do
    end do
    do k=1,nz
      do j=2,nym
        dudy (j,k) = (tmp1(j-1,k)+tmp1(j,k))*0.5
        dmudy(j,k) = (tmp2(j-1,k)+tmp2(j,k))*0.5
      end do
      dudy (1 ,k) = dudy (2  ,k)
      dudy (ny,k) = dudy (nym,k)
      dmudy(1 ,k) = dmudy(2  ,k)
      dmudy(ny,k) = dmudy(nym,k)
    end do
 
!-------------------------------------------------------------------------------
!   Calculate zonally averaged meridional and vertical heat flux and
!   meridional and vertical momentum flux
!--------------------------------------------------------------------------------

    hfy  (1:ny,1:nz) = 0.d0; hfz  (1:ny,1:nz) = 0.d0
    mfy  (1:ny,1:nz) = 0.d0; mfz  (1:ny,1:nz) = 0.d0
    hfyst(1:ny,1:nz) = 0.d0; hfzst(1:ny,1:nz) = 0.d0
    mfyst(1:ny,1:nz) = 0.d0; mfzst(1:ny,1:nz) = 0.d0
    hfytr(1:ny,1:nz) = 0.d0; hfztr(1:ny,1:nz) = 0.d0
    mfytr(1:ny,1:nz) = 0.d0; mfztr(1:ny,1:nz) = 0.d0

    do k=1,nz
    do j=1,ny
      do i=1,nx
        hfy  (j,k) = hfy  (j,k) + vp  (i,j,k)*thp  (i,j,k)
        hfz  (j,k) = hfz  (j,k) + wp  (i,j,k)*thp  (i,j,k)
        mfy  (j,k) = mfy  (j,k) + vp  (i,j,k)*up   (i,j,k)
        mfz  (j,k) = mfz  (j,k) + wp  (i,j,k)*up   (i,j,k)
        hfyst(j,k) = hfyst(j,k) + vpst(i,j,k)*thpst(i,j,k)
        hfzst(j,k) = hfzst(j,k) + wpst(i,j,k)*thpst(i,j,k)
        mfyst(j,k) = mfyst(j,k) + vpst(i,j,k)*upst (i,j,k)
        mfzst(j,k) = mfzst(j,k) + wpst(i,j,k)*upst (i,j,k)
        hfytr(j,k) = hfytr(j,k) + vptr(i,j,k)*thptr(i,j,k)
        hfztr(j,k) = hfztr(j,k) + wptr(i,j,k)*thptr(i,j,k)
        mfytr(j,k) = mfytr(j,k) + vptr(i,j,k)*uptr (i,j,k)
        mfztr(j,k) = mfztr(j,k) + wptr(i,j,k)*uptr (i,j,k)
      end do
      hfy  (j,k) = hfy  (j,k)/float(nx)
      hfz  (j,k) = hfz  (j,k)/float(nx)
      mfy  (j,k) = mfy  (j,k)/float(nx)
      mfz  (j,k) = mfz  (j,k)/float(nx)
      hfyst(j,k) = hfyst(j,k)/float(nx)
      hfzst(j,k) = hfzst(j,k)/float(nx)
      mfyst(j,k) = mfyst(j,k)/float(nx)
      mfzst(j,k) = mfzst(j,k)/float(nx)
      hfytr(j,k) = hfytr(j,k)/float(nx)
      hfztr(j,k) = hfztr(j,k)/float(nx)
      mfytr(j,k) = mfytr(j,k)/float(nx)
      mfztr(j,k) = mfztr(j,k)/float(nx)
    end do
    end do

!---------------------------------------------------------------------------------
!
!   Calculate residual mean meridional circulation (v* and w*) vstar and wstar
!   are defined at mid points of the horizontal and vertical grid
!
!                     1     d           <u'theta'>
!      <v>* = <v> - ------ ---- ( rho0 ------------ )
!                    rho0   dz          d<theta>/dz
!
!                        1         d                 <u'theta'>
!      <w>* = <w> + ------------- ------ ( cos(phi) ------------ )
!                     a cos(phi)   dphi              d<theta>/dz
!
!      < > means the zonal averaging
! 
!--------------------------------------------------------------------------------

    tmp1(1:ny,1:nz) = 0.d0
    do k=1,nzm
    do j=1,ny
      fac0 = 1.0/rho0m(k)
      fac1 = rho0(k  )*hfy(j,k  )/dtdz(j,k  ) 
      fac2 = rho0(k+1)*hfy(j,k+1)/dtdz(j,k+1)
      fac3 = logz(k)-logz(k+1)
      tmp1(j,k) = fac0*((fac1-fac2)/fac3)
    end do
    end do
    do k=1,nzm
    do j=1,nym
      fac0 = (zmv(j,k)+zmv(j+1,k)+zmv(j,k+1)+zmv(j+1,k+1))*0.25
      fac1 = (tmp1(j,k)+tmp1(j+1,k))*0.5
      vs(j,k) = fac0-fac1
    end do
    end do

    tmp1(1:ny,1:nz) = 0.d0
    do k=1,nz
    do j=1,nym
      fac0 = 1.0/(a*cosphm(j))
      fac1 = cosph(j+1)*hfy(j+1,k)/dtdz(j+1,k)
      fac2 = cosph(j  )*hfy(j  ,k)/dtdz(j  ,k)
      fac3 = (lat(j+1)-lat(j))*pi/180.
      tmp1(j,k) = fac0*((fac1-fac2)/fac3)
    end do
    end do
    do k=1,nzm
    do j=1,nym
      fac0 = (zmw(j,k)+zmw(j+1,k)+zmw(j,k+1)+zmw(j+1,k+1))*0.25
      fac1 = (tmp1(j,k+1)+tmp1(j,k))*0.5
      ws(j,k) = fac0+fac1
    end do
    end do

!-------------------------------------------------------------------------------
!   Calculate meridional advection terms in TEM zonal wind equation
!
!                 1     d<u>
!     tmp1    =  --- ( ------ ) 
!                 a     dphi
!
!     advy    = - <v>* advy1
!
!                       d<u>
!     advz    = - <w>* ------
!                        dz
!--------------------------------------------------------------------------------

    tmp1(1:ny,1:nz) = 0.d0
    do k=1,nz
    do j=1,nym
      fac0 = 1.0/a
      fac1 = zmu(j+1,k)-zmu(j,k)
      fac2 = (lat(j+1)-lat(j))*pi/180.
      tmp1(j,k) = fac0*fac1/fac2
    end do
    end do
    do k=1,nzm
    do j=1,nym
      fac0 = (tmp1(j,k+1)+tmp1(j,k))*0.5
      fac1 = (dudz(j,k)+dudz(j+1,k)+dudz(j,k+1)+dudz(j+1,k+1))*0.25
      advyu(j,k) = -vs(j,k)*fac0
      advzu(j,k) = -ws(j,k)*fac1 
    end do
    end do

!--------------------------------------------------------------------------------
!   Calculate meridional temperature advection term
!
!               <v>* d<theta>
!   advyt   = - ---- --------
!                a     dphi
!
!                     d<theta>
!   advzt   = - <w>* ---------
!                        dz
!--------------------------------------------------------------------------------

    tmp1(1:ny,1:nz) = 0.d0
    do k=1,nz
    do j=1,nym
      fac0 = 1.0/a
      fac1 = zmth(j+1,k)-zmth(j,k)
      fac2 = (lat(j+1)-lat(j))*pi/180.
      tmp1(j,k) = fac0*fac1/fac2
    end do
    end do
    do k=1,nzm
    do j=1,nym
      fac0 = (tmp1(j,k+1)+tmp1(j,k))*0.5
      fac1 = (dtdz(j,k)+dtdz(j+1,k)+dtdz(j,k+1)+dtdz(j+1,k+1))*0.25
      advyt(j,k) = -vs(j,k)*fac0
      advzt(j,k) = -ws(j,k)*fac1
    end do
    end do

!-------------------------------------------------------------------------------
!   Curvature term
!-------------------------------------------------------------------------------

    do k=1,nzm
    do j=1,nym
      fac0 = 1.0/a
      fac1 = (zmu  (j,k)+zmu  (j+1,k)+zmu  (j,k+1)+zmu  (j+1,k+1))*0.25
      fac2 = (zmtmu(j,k)+zmtmu(j+1,k)+zmtmu(j,k+1)+zmtmu(j+1,k+1))*0.25
      fac3 = dtan(latm(j)*pi/180.)
      curv (j,k) = fac0*fac1*fac3*vs(j,k)
    end do
    end do

!-------------------------------------------------------------------------------
!   Corioris term
!-------------------------------------------------------------------------------

    do k=1,nzm
    do j=1,nym
      cori (j,k) = vs(j,k)*2.0*radang*sinphm(j)
    end do
    end do

!-------------------------------------------------------------------------------
!   Eliassen-Palm flux and its divergence
!
!                               d<u>   <v'theta'>
!     epy  = rho0*a*cos(phi)*( ------ ------------ - <v'u'> )
!                                dz    d<theta>/dz
!
!                                       1       d                     <v'theta'>
!     epz  = rho0*a*cos(phi)*{( f - ---------- ---- ( <u>cos(phi) ) ) -----------
!                                   a*cos(phi) dphi                   d<theta>/dz
!
!            -  <w'u'> }
! 
!
!     [epdiv] = m s^(-2) 
!--------------------------------------------------------------------------------

    tmp1(1:ny,1:nz) = 0.d0
    tmp2(1:ny,1:nz) = 0.d0
    tmp3(1:ny,1:nz) = 0.d0
    do k=1,nz
    do j=1,ny
      fac0 = rho0(k)*a*cosph(j)
      fac1 = dudz (j,k)*hfy  (j,k)/dtdz (j,k) - mfy  (j,k)
      fac2 = dmudz(j,k)*hfyst(j,k)/dmtdz(j,k) - mfyst(j,k)
      fac3 = dmudz(j,k)*hfytr(j,k)/dmtdz(j,k) - mfytr(j,k)
      tmp1(j,k) = fac0*fac1
      tmp2(j,k) = fac0*fac2
      tmp3(j,k) = fac0*fac3
    end do
    end do
    do k=1,nzm
    do j=1,nym
      epy  (j,k) = (tmp1(j,k)+tmp1(j+1,k)+tmp1(j,k+1)+tmp1(j+1,k+1))*0.25
      epyst(j,k) = (tmp2(j,k)+tmp2(j+1,k)+tmp2(j,k+1)+tmp2(j+1,k+1))*0.25
      epytr(j,k) = (tmp3(j,k)+tmp3(j+1,k)+tmp3(j,k+1)+tmp3(j+1,k+1))*0.25
    end do
    end do
    
    tmp4(1:ny,1:nz) = 0.d0
    tmp5(1:ny,1:nz) = 0.d0
    tmp6(1:ny,1:nz) = 0.d0
    do k=1,nz
    do j=1,nym
      fac0 = 1.0/(a*cosphm(j))
      fac1 = tmp1(j+1,k)*cosph(j+1) - tmp1(j,k)*cosph(j)
      fac2 = tmp2(j+1,k)*cosph(j+1) - tmp2(j,k)*cosph(j)
      fac3 = tmp3(j+1,k)*cosph(j+1) - tmp3(j,k)*cosph(j)
      fac4 = (lat(j+1)-lat(j))*pi/180.
      tmp4(j,k) = fac0*fac1/fac4
      tmp5(j,k) = fac0*fac2/fac4
      tmp6(j,k) = fac0*fac3/fac4
    end do
    end do

    epd  (1:nym,1:nzm) = 0.d0
    epdst(1:nym,1:nzm) = 0.d0
    epdtr(1:nym,1:nzm) = 0.d0
    do k=1,nzm
    do j=1,nym
      fac0 = 1.0/(rho0m(k)*a*cosphm(j))
      epd  (j,k) = fac0*(tmp4(j,k+1)+tmp4(j,k))*0.5
      epdst(j,k) = fac0*(tmp5(j,k+1)+tmp5(j,k))*0.5
      epdtr(j,k) = fac0*(tmp6(j,k+1)+tmp6(j,k))*0.5
    end do
    end do

    tmp1(1:ny,1:nz) = 0.d0
    tmp2(1:ny,1:nz) = 0.d0
    tmp3(1:ny,1:nz) = 0.d0
    do k=1,nz
    do j=1,ny
      fac0 = rho0(k)*a*cosph(j)
      fac1 = 2.0*radang*sinph(j)-(1.0/a)*dudz (j,k)+zmu(j,k)*dtan(lat(j)*pi/180.)/a
      fac2 = 2.0*radang*sinph(j)-(1.0/a)*dmudz(j,k)+zmu(j,k)*dtan(lat(j)*pi/180.)/a
      fac3 = fac1*hfy  (j,k)/dtdz (j,k) - mfz  (j,k)
      fac4 = fac2*hfyst(j,k)/dmtdz(j,k) - mfzst(j,k)
      fac5 = fac2*hfytr(j,k)/dmtdz(j,k) - mfztr(j,k)
      tmp1(j,k) = fac0*fac3
      tmp2(j,k) = fac0*fac4
      tmp3(j,k) = fac0*fac5
    end do
    end do
    do k=1,nzm
    do j=1,nym
      epz  (j,k) = (tmp1(j,k)+tmp1(j+1,k)+tmp1(j,k+1)+tmp1(j+1,k+1))*0.25
      epzst(j,k) = (tmp2(j,k)+tmp2(j+1,k)+tmp2(j,k+1)+tmp2(j+1,k+1))*0.25
      epztr(j,k) = (tmp3(j,k)+tmp3(j+1,k)+tmp3(j,k+1)+tmp3(j+1,k+1))*0.25
    end do
    end do

    tmp4(1:ny,1:nz) = 0.d0
    tmp5(1:ny,1:nz) = 0.d0
    tmp6(1:ny,1:nz) = 0.d0
    do k=1,nzm
    do j=1,ny
      fac0 = tmp1(j,k)-tmp1(j,k+1)
      fac1 = tmp2(j,k)-tmp2(j,k+1)
      fac2 = tmp3(j,k)-tmp3(j,k+1)
      fac3 = logz(k)-logz(k+1)
      tmp4(j,k) = fac0/fac3
      tmp5(j,k) = fac1/fac3
      tmp6(j,k) = fac2/fac3
    end do
    end do

    do k=1,nzm
    do j=1,nym
      fac0 = 1.0/(rho0m(k)*a*cosphm(j))
      epd  (j,k) = epd  (j,k) + fac0*(tmp4(j+1,k)+tmp4(j,k))*0.5
      epdst(j,k) = epdst(j,k) + fac0*(tmp5(j+1,k)+tmp5(j,k))*0.5
      epdtr(j,k) = epdtr(j,k) + fac0*(tmp6(j+1,k)+tmp6(j,k))*0.5
    end do
    end do

!-------------------------------------------------------------------------------
!   Temperature forcing
!
!              1   d    rho  <v'th'>  d<th>
!   term1 = - --- --- ( --- --------- ----- )
!             rho  dz    a   d<th>/dz   dz
!
!              1   d    
!   term2 = - --- --- ( rho <w'th'> )
!             rho  dz
!-------------------------------------------------------------------------------

    tmp1(1:ny,1:nz) = 0.d0
    tmp2(1:ny,1:nz) = 0.d0
    tmp3(1:ny,1:nz) = 0.d0
    do k=1,nz
    do j=1,ny
      fac0 = rho0(k)/a
      fac1 = hfy  (j,k)/dtdz (j,k)*dtdy (j,k)
      fac2 = hfyst(j,k)/dmtdz(j,k)*dmtdy(j,k)
      fac3 = hfytr(j,k)/dmtdz(j,k)*dmtdy(j,k)
      tmp1(j,k) = fac0*fac1
      tmp2(j,k) = fac0*fac2
      tmp3(j,k) = fac0*fac3
    end do
    end do

    tmp4(1:ny,1:nz) = 0.d0
    tmp5(1:ny,1:nz) = 0.d0
    tmp6(1:ny,1:nz) = 0.d0
    do k=1,nzm
    do j=1,ny
      fac0 = -1.0/rho0m(k)
      fac1 = tmp1(j,k)-tmp1(j,k+1)
      fac2 = tmp2(j,k)-tmp2(j,k+1)
      fac3 = tmp3(j,k)-tmp3(j,k+1)
      fac4 = logz(k)-logz(k+1)
      tmp4(j,k) = fac0*fac1/fac4
      tmp5(j,k) = fac0*fac2/fac4
      tmp6(j,k) = fac0*fac3/fac4
    end do
    end do
   
    ft  (1:nym,1:nzm) = 0.d0
    ftst(1:nym,1:nzm) = 0.d0
    fttr(1:nym,1:nzm) = 0.d0
    do k=1,nzm
    do j=1,nym
      ft  (j,k) = (tmp4(j,k)+tmp4(j+1,k))*0.5
      ftst(j,k) = (tmp5(j,k)+tmp5(j+1,k))*0.5
      fttr(j,k) = (tmp6(j,k)+tmp6(j+1,k))*0.5
    end do
    end do
 
    tmp1(1:ny,1:nz) = 0.d0
    tmp2(1:ny,1:nz) = 0.d0
    tmp3(1:ny,1:nz) = 0.d0
    do k=1,nz
    do j=1,ny
      tmp1(j,k) = rho0(k)*hfz  (j,k)
      tmp2(j,k) = rho0(k)*hfzst(j,k)
      tmp3(j,k) = rho0(k)*hfztr(j,k)
    end do
    end do

    tmp4(1:ny,1:nz) = 0.d0
    tmp5(1:ny,1:nz) = 0.d0
    tmp6(1:ny,1:nz) = 0.d0
    do k=1,nzm
    do j=1,ny
      fac0 = -1.0/rho0m(k)
      fac1 = tmp1(j,k)-tmp1(j,k+1)
      fac2 = tmp2(j,k)-tmp2(j,k+1)
      fac3 = tmp3(j,k)-tmp3(j,k+1)
      fac4 = logz(k)-logz(k+1)
      tmp4(j,k) = fac0*fac1/fac4
      tmp5(j,k) = fac0*fac2/fac4
      tmp6(j,k) = fac0*fac3/fac4
    end do
    end do

    do k=1,nzm
    do j=1,nym
      ft  (j,k) = ft  (j,k) + (tmp4(j,k)+tmp4(j+1,k))*0.5
      ftst(j,k) = ftst(j,k) + (tmp5(j,k)+tmp5(j+1,k))*0.5
      fttr(j,k) = fttr(j,k) + (tmp6(j,k)+tmp6(j+1,k))*0.5
    end do
    end do

    corif  = real(cori);  curvf  = real(curv)
    vsf    = real(vs);    wsf    = real(ws)
    advyuf = real(advyu); advzuf = real(advzu)
    epyf   = real(epy);   epystf = real(epyst); epytrf = real(epytr)
    epzf   = real(epz);   epzstf = real(epzst); epztrf = real(epztr)
    epdf   = real(epd);   epdstf = real(epdst); epdtrf = real(epdtr)
    advytf = real(advyt); advztf = real(advzt)
    ftf    = real(ft);    ftstf  = real(ftst);  fttrf  = real(fttr)
    hfyf   = real(hfy);   hfystf = real(hfyst); hfytrf = real(hfytr)
    hfzf   = real(hfz);   hfzstf = real(hfzst); hfztrf = real(hfztr)
    mfyf   = real(mfy);   mfystf = real(mfyst); mfytrf = real(mfytr)
    mfzf   = real(mfz);   mfzstf = real(mfzst); mfztrf = real(mfztr)

    return
    end subroutine epfluxsttr    

   subroutine epfluxst(nx    ,ny    ,nz    ,nym   ,nzm   ,lat   ,lev   , &
                       utm   ,vtm   ,ttm   ,omgatm,latm  ,logzm ,logz  , &
                       epystf,epzstf,epdstf,ftstf ,hfystf,hfzstf,mfystf,  &
                       mfzstf)

!-----------------------------------------------------------------------------
!
!   USAGE
!
!   INPUT VARIABLES
!
!   nx     : The number of grids in the longitude
!   ny     : The number of grids in the latitude
!   nz     : The number of grids in the vertical
!   nym    : The number of grids in the latitude ( = ny-1 )
!   nzm    : The number of grids in the vertical ( = nz-1 )
!   lat    : Latitude (degree) in an increasing order
!   lev    : Pressure (mb) in an increasing order
!   utm    : Time-averaged zonal wind (m/s)
!   vtm    : Time-averaged meridional wind (m/s)
!   ttm    : Time-averaged temperature (K)
!   omgatm : Time-averaged pressure velocity (Pa/s) 
!
!   OUTPUT VARIABLES
!
!   latm   : Latitude (degree)
!   logzm  : Log-Pressure height (m)
!   logz   : Log-Pressure height (m)
!   epystf : Meridional component of EP flux due to stationary wave
!   epzstf : Vertical component of EP flux due to stationary wave
!   epdstf : EP flux divergence due to stationary wave
!   ftstf  : Potential temperature forcing due to stationary wave
!   hfystf : Meridional heat flux due to stationary wave
!   hfzstf : Vertical heat flux due to stationary wave
!   mfystf : Meridional momentum flux due to stationary wave
!   mfzstf : Vertical momentum flux due to stationary wave
!
!-----------------------------------------------------------------------------

    implicit none

    integer, parameter :: r8 = selected_real_kind(12)
    integer, parameter :: i8 = selected_int_kind(13)

    real(r8), parameter :: a = 6370000.
    real(r8), parameter :: h = 7000.
    real(r8), parameter :: rhos = 1.125
    real(r8), parameter :: ps = 1000.0
    real(r8), parameter :: g = 9.81
    real(r8), parameter :: rd = 287.0
    real(r8), parameter :: cp = 1004.
    real(r8), parameter :: pi = 3.141592
    real(r8), parameter :: radang = 2.*pi/86400.

    integer :: i,j,k,it

    integer,                   intent(in)  :: nx    ,ny    ,nz
    integer,                   intent(in)  :: nym   ,nzm
    real, dimension(ny),       intent(in)  :: lat
    real, dimension(nz),       intent(in)  :: lev
    real, dimension(nx,ny,nz), intent(in)  :: utm   ,vtm   ,ttm   ,omgatm

    real, dimension(nym),      intent(out) :: latm
    real, dimension(nzm),      intent(out) :: logzm
    real, dimension(nz),       intent(out) :: logz
    real, dimension(nym,nzm),  intent(out) :: epystf,epzstf,epdstf,ftstf
    real, dimension(ny,nz),    intent(out) :: hfystf,hfzstf
    real, dimension(ny,nz),    intent(out) :: mfystf,mfzstf

    real(r8), dimension(ny)         :: sinph  ,cosph
    real(r8), dimension(nym)        :: sinphm ,cosphm
    real(r8), dimension(nz)         :: rho0
    real(r8), dimension(nzm)        :: rho0m

    real(r8), dimension(nx,ny,nz)   :: thtm  ,wtm
    real(r8), dimension(ny,nz)      :: zmtmu ,zmtmv ,zmtmth,zmtmw

    real(r8), dimension(nx,ny,nz)   :: upst  ,vpst  ,wpst  ,thpst
    real(r8), dimension(ny,nz)      :: dmtdz ,dmudz ,dmtdy ,dmudy
    real(r8), dimension(ny,nz)      :: hfyst ,hfzst ,mfyst ,mfzst

    real(r8), dimension(nym,nzm)    :: epyst ,epzst ,epdst
    real(r8), dimension(nym,nzm)    :: ftst 

    real(r8), dimension(ny,nz)      :: tmp1  ,tmp2  ,tmp3

    real(r8) :: fac0,fac1,fac2,fac3,fac4
    
!------------------------------------------------------------------------
!   Define coordinate 
!------------------------------------------------------------------------

    do j=1,ny
      sinph(j) = dsin(lat(j)*pi/180.0)
      cosph(j) = dcos(lat(j)*pi/180.0)
    end do
    do j=1,nym
      latm(j) = (lat(j)+lat(j+1))*0.5
      sinphm(j) = dsin(latm(j)*pi/180.0)
      cosphm(j) = dcos(latm(j)*pi/180.0)
    end do

    do k=1,nz
      logz(k) = h*log(ps/lev(k))       ! logz(1) = the highest level
      rho0(k) = rhos*dexp(-logz(k)/h)
    end do 
    do k=1,nzm
      logzm(k) = (logz(k)+logz(k+1))*0.5
      rho0m(k) = rhos*dexp(-logzm(k)/h)
    end do 

!------------------------------------------------------------------------
!   Define time-averaged theta and w 
!------------------------------------------------------------------------

    do k=1,nz
    do j=1,ny
    do i=1,nx
      thtm(i,j,k) = ttm(i,j,k)*(ps/lev(k))**(rd/cp)
      wtm (i,j,k) = -(h*omgatm(i,j,k))/(lev(k)*100.0)
    end do 
    end do
    end do

!------------------------------------------------------------------------
!   Zonal- and time-averaged u, v, w, and theta 
!------------------------------------------------------------------------

    zmtmu (1:ny,1:nz) = 0.d0
    zmtmv (1:ny,1:nz) = 0.d0
    zmtmth(1:ny,1:nz) = 0.d0
    zmtmw (1:ny,1:nz) = 0.d0
    do k=1,nz
    do j=1,ny
      do i=1,nx
        zmtmu (j,k) = zmtmu (j,k) + dble(utm(i,j,k))
        zmtmv (j,k) = zmtmv (j,k) + dble(vtm(i,j,k))
        zmtmth(j,k) = zmtmth(j,k) + thtm(i,j,k)
        zmtmw (j,k) = zmtmw (j,k) + wtm (i,j,k)
      end do
      zmtmu (j,k) = zmtmu (j,k)/float(nx)
      zmtmv (j,k) = zmtmv (j,k)/float(nx)
      zmtmth(j,k) = zmtmth(j,k)/float(nx)
      zmtmw (j,k) = zmtmw (j,k)/float(nx)
    end do
    end do

!------------------------------------------------------------------------
!   Stationary wave perturbation
!------------------------------------------------------------------------

    do k=1,nz
    do j=1,ny
    do i=1,nx
      upst (i,j,k) = dble(utm(i,j,k)) - zmtmu (j,k)
      vpst (i,j,k) = dble(vtm(i,j,k)) - zmtmv (j,k)
      wpst (i,j,k) =      wtm(i,j,k)  - zmtmw (j,k)
      thpst(i,j,k) =     thtm(i,j,k)  - zmtmth(j,k)
    end do
    end do
    end do

!--------------------------------------------------------------------------------
!   Calculate vertical derivative of zonal mean potential temperature.
!   Finally, dmtdz are defined at the vertical levels of lev.
!   dmtdz is time-invariant d<[theta]>/dz.
!   dmtdz is extrapolated at both the top and bottom boundaries.
!--------------------------------------------------------------------------------

    tmp1(1:ny,1:nz) = 0.d0
    do k=1,nzm
      do j=1,ny
        fac0 = zmtmth(j,k)-zmtmth(j,k+1)
        fac1 = logz(k)-logz(k+1)
        tmp1(j,k) = fac0/fac1
      end do
    end do
    do j=1,ny
      do k=2,nzm
        dmtdz(j,k) = (tmp1(j,k)+tmp1(j,k-1))*0.5
      end do
      dmtdz(j,1)  = dmtdz(j,2) 
      dmtdz(j,nz) = dmtdz(j,nz-1)
    end do

    do k=1,nz
    do j=1,ny
      if ( dmtdz(j,k) <= 0.0 ) then
        print *,'Invalid vertical gradient of potential temperature'
      end if
    end do
    end do 

!--------------------------------------------------------------------------------
!   Calculate the vertical derivative of zonal mean zonal wind.
!   Finally, dmudz are defined at the vertical levels of lev.
!   dmudz is time-invariant d<[u]>/dz.
!   dmudz is extrapolated at both the top and bottom boundaries.
!--------------------------------------------------------------------------------

    tmp1(1:ny,1:nz) = 0.d0
    do k=1,nzm
      do j=1,ny
        fac0 = zmtmu(j,k)-zmtmu(j,k+1)
        fac1 = logz(k)-logz(k+1)
        tmp1(j,k) = fac0/fac1
      end do
    end do
    do j=1,ny
      do k=2,nzm
        dmudz(j,k) = (tmp2(j,k)+tmp2(j,k-1))*0.5
      end do
      dmudz(j, 1) = dmudz(j,2) 
      dmudz(j,nz) = dmudz(j,nz-1) 
    end do

!--------------------------------------------------------------------------------
!   Calculate meridional derivative of zonal-mean potential temperature.
!   Finally, dmtdy are defined at the grid of lat.
!   dmtdy is time-invariant d<[theta]>/dphi.
!   dmtdy is extrapolated at both poles.
!--------------------------------------------------------------------------------

    tmp1(1:ny,1:nz) = 0.d0
    do k=1,nz
      do j=1,nym
        fac0 = zmtmth(j+1,k)-zmtmth(j,k)
        fac1 = (lat(j+1)-lat(j))*pi/180.
        tmp1(j,k) = fac0/fac1
      end do
    end do
    do k=1,nz
      do j=2,nym
        dmtdy(j,k) = (tmp1(j-1,k)+tmp1(j,k))*0.5
      end do
      dmtdy(1 ,k) = dmtdy(2  ,k)
      dmtdy(ny,k) = dmtdy(nym,k)
    end do

!-------------------------------------------------------------------------------
!   Calculate meridional derivative of zonal-mean zonal wind
!   Finally, dmudy are defined at the grid of lat
!   dmudy is time-invariant d<[u]>/dphi
!   dmudy is extrapolated at both poles.
!-------------------------------------------------------------------------------

    tmp1(1:ny,1:nz) = 0.d0
    do k=1,nz
      do j=1,nym
        fac0 = zmtmu(j+1,k)-zmtmu(j,k)
        fac1 = (lat(j+1)-lat(j))*pi/180.
        tmp2(j,k) = fac0/fac1
      end do
    end do
    do k=1,nz
      do j=2,nym
        dmudy(j,k) = (tmp1(j-1,k)+tmp1(j,k))*0.5
      end do
      dmudy(1 ,k) = dmudy(2  ,k)
      dmudy(ny,k) = dmudy(nym,k)
    end do

!-------------------------------------------------------------------------------
!   Calculate zonally averaged meridional and vertical heat flux and
!   meridional and vertical momentum flux
!--------------------------------------------------------------------------------

    hfyst(1:ny,1:nz) = 0.d0; hfzst(1:ny,1:nz) = 0.d0
    mfyst(1:ny,1:nz) = 0.d0; mfzst(1:ny,1:nz) = 0.d0
    do k=1,nz
    do j=1,ny
      do i=1,nx
        hfyst(j,k) = hfyst(j,k) + vpst(i,j,k)*thpst(i,j,k)
        hfzst(j,k) = hfzst(j,k) + wpst(i,j,k)*thpst(i,j,k)
        mfyst(j,k) = mfyst(j,k) + vpst(i,j,k)*upst (i,j,k)
        mfzst(j,k) = mfzst(j,k) + wpst(i,j,k)*upst (i,j,k)
      end do
      hfyst(j,k) = hfyst(j,k)/float(nx)
      hfzst(j,k) = hfzst(j,k)/float(nx)
      mfyst(j,k) = mfyst(j,k)/float(nx)
      mfzst(j,k) = mfzst(j,k)/float(nx)
    end do
    end do

!-------------------------------------------------------------------------------
!   Eliassen-Palm flux and its divergence
!
!                               d<u>   <v'theta'>
!     epy  = rho0*a*cos(phi)*( ------ ------------ - <v'u'> )
!                                dz    d<theta>/dz
!
!                                       1       d                     <v'theta'>
!     epz  = rho0*a*cos(phi)*{( f - ---------- ---- ( <u>cos(phi) ) ) -----------
!                                   a*cos(phi) dphi                   d<theta>/dz
!
!            -  <w'u'> }
! 
!
!     [epdiv] = m s^(-2) 
!--------------------------------------------------------------------------------

    tmp1(1:ny,1:nz) = 0.d0
    do k=1,nz
    do j=1,ny
      fac0 = rho0(k)*a*cosph(j)
      fac1 = dmudz(j,k)*hfyst(j,k)/dmtdz(j,k) - mfyst(j,k)
      tmp1(j,k) = fac0*fac1
    end do
    end do
    do k=1,nzm
    do j=1,nym
      epyst(j,k) = (tmp1(j,k)+tmp1(j+1,k)+tmp1(j,k+1)+tmp1(j+1,k+1))*0.25
    end do
    end do

    tmp2(1:ny,1:nz) = 0.d0
    do k=1,nz
    do j=1,nym
      fac0 = 1.0/(a*cosphm(j))
      fac1 = tmp1(j+1,k)*cosph(j+1) - tmp1(j,k)*cosph(j)
      fac2 = (lat(j+1)-lat(j))*pi/180.
      tmp2(j,k) = fac0*fac1/fac2
    end do
    end do

    epdst(1:nym,1:nzm) = 0.d0
    do k=1,nzm
    do j=1,nym
      fac0 = 1.0/(rho0m(k)*a*cosphm(j))
      epdst(j,k) = fac0*(tmp2(j,k+1)+tmp2(j,k))*0.5
    end do
    end do

    tmp1(1:ny,1:nz) = 0.d0
    do k=1,nz
    do j=1,ny
      fac0 = rho0(k)*a*cosph(j)
      fac1 = 2.0*radang*sinph(j)-(1.0/a)*dmudz(j,k)+zmtmu(j,k)*dtan(lat(j)*pi/180.)/a
      fac2 = fac1*hfyst(j,k)/dmtdz(j,k) - mfzst(j,k)
      tmp1(j,k) = fac0*fac2
    end do
    end do
    do k=1,nzm
    do j=1,nym
      epzst(j,k) = (tmp1(j,k)+tmp1(j+1,k)+tmp1(j,k+1)+tmp1(j+1,k+1))*0.25
    end do
    end do

    tmp2(1:ny,1:nz) = 0.d0
    do k=1,nzm
    do j=1,ny
      fac0 = tmp1(j,k)-tmp1(j,k+1)
      fac1 = logz(k)-logz(k+1)
      tmp2(j,k) = fac0/fac1
    end do
    end do

    do k=1,nzm
    do j=1,nym
      fac0 = 1.0/(rho0m(k)*a*cosphm(j))
      epdst(j,k) = epdst(j,k) + fac0*(tmp2(j+1,k)+tmp2(j,k))*0.5
    end do
    end do

!-------------------------------------------------------------------------------
!   Temperature forcing
!
!              1   d    rho  <v'th'>  d<th>
!   term1 = - --- --- ( --- --------- ----- )
!             rho  dz    a   d<th>/dz   dz
!
!              1   d    
!   term2 = - --- --- ( rho <w'th'> )
!             rho  dz
!-------------------------------------------------------------------------------

    tmp1(1:ny,1:nz) = 0.d0
    do k=1,nz
    do j=1,ny
      fac0 = rho0(k)/a
      fac1 = hfyst(j,k)/dmtdz(j,k)*dmtdy(j,k)
      tmp1(j,k) = fac0*fac1
    end do
    end do

    tmp2(1:ny,1:nz) = 0.d0
    do k=1,nzm
    do j=1,ny
      fac0 = -1.0/rho0m(k)
      fac1 = tmp1(j,k)-tmp1(j,k+1)
      fac2 = logz(k)-logz(k+1)
      tmp2(j,k) = fac0*fac1/fac2
    end do
    end do
   
    ftst(1:nym,1:nzm) = 0.d0

    do k=1,nzm
    do j=1,nym
      ftst(j,k) = (tmp2(j,k)+tmp2(j+1,k))*0.5
    end do
    end do
 
    tmp1(1:ny,1:nz) = 0.d0
    do k=1,nz
    do j=1,ny
      tmp1(j,k) = rho0(k)*hfzst(j,k)
    end do
    end do

    tmp2(1:ny,1:nz) = 0.d0
    do k=1,nzm
    do j=1,ny
      fac0 = -1.0/rho0m(k)
      fac1 = tmp1(j,k)-tmp1(j,k+1)
      fac2 = logz(k)-logz(k+1)
      tmp2(j,k) = fac0*fac1/fac2
    end do
    end do

    do k=1,nzm
    do j=1,nym
      ftst(j,k) = ftst(j,k) + (tmp2(j,k)+tmp2(j+1,k))*0.5
    end do
    end do

    epystf = real(epyst)
    epzstf = real(epzst)
    epdstf = real(epdst) 
    ftstf  = real(ftst)
    hfystf = real(hfyst)
    hfzstf = real(hfzst)
    mfystf = real(mfyst)
    mfzstf = real(mfzst)

    return
    end subroutine epfluxst  

    end module epfluxsph
