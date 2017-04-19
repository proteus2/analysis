MODULE mflx_ctop_sc05

  implicit none

  real, dimension(:), allocatable ::  u_sfc, v_sfc
  real, dimension(:), allocatable ::  u_ct, v_ct, u_cb, v_cb, t_ct,      &
                                      n_q, n_ct, rho_ct, cqx, cqy,       &
                                      heatmax, zcta, zcba
  integer, dimension(:), allocatable ::  kcta
 
  real, dimension(:,:,:), allocatable ::  mfs_ct

  real, dimension(:,:)  , allocatable ::  diag_znwcq

  real, dimension(:,:,:), allocatable ::  diag_spec_ctop_0

  real, dimension(:), allocatable ::  mflx_ct_0_east, mflx_ct_0_west,    &
                                      mflx_ct_0_north, mflx_ct_0_south


CONTAINS


SUBROUTINE args_sc05(ncol,nz,u,v,t,nbv,rho,z,heat,kcb,kct,z_ref)
! The all input here are at flux levels.

!
! PURPOSE:  To obtain the parameters used in the GWDC scheme of
!           Song and Chun (2005, JAS), with a modification of (Cqx, Cqy)
!           (Choi and Chun, 2011, JAS)
!
! METHOD:
!
!   1) Determine the heating center height and half depth, and maximum 
!      heating rate.
!   2) Estimate moving velocity of convective cells.
!
! HISTORY:
!

  USE switch_dump,  ONLY: l_diag_znwcq_o

  implicit none

! SUBROUTINE ARGUMENTS

  integer, intent(in) ::  ncol, nz

  ! data arrays
  real, dimension(ncol,nz), intent(in) ::  u, v, t, nbv, rho, z, heat
  integer, dimension(ncol), intent(in) ::  kcb, kct
  real, dimension(nz)     , intent(in), optional ::  z_ref

! LOCAL VARIABLES

  ! data arrays
  real   , dimension(ncol) ::  zm, zd
  integer, dimension(ncol) ::  kcba

  ! work arrays
  real ::  tmp, eta_3km, tmp1d(ncol)
  real, dimension(nz) ::  eta
  real, dimension(ncol,nz) ::  heatdz, dz
  integer ::  tmpi, cnt, k_3km
  integer ::  k,l   ! loop counters

  ! parameters and constants
  real, parameter ::  zt_eta = 100.e3  ! 100 km

  if ( allocated(kcta) )  deallocate( kcta )
  allocate( kcta(ncol) )
  if ( allocated(u_ct) )  deallocate( u_ct, v_ct, u_cb, v_cb, t_ct,      &
                                      zcta, zcba, n_q, n_ct, rho_ct,     &
                                      cqx, cqy, heatmax )
  allocate( u_ct(ncol), v_ct(ncol), u_cb(ncol), v_cb(ncol),              &
            t_ct(ncol), zcta(ncol), zcba(ncol), n_q(ncol), n_ct(ncol),   &
            rho_ct(ncol), cqx(ncol), cqy(ncol), heatmax(ncol) )
 
  if ( l_diag_znwcq_o ) then
    if ( allocated(diag_znwcq) )  deallocate( diag_znwcq )
    allocate( diag_znwcq(ncol,10) )
  end if

!
! Determine heating center height and half depth,
! following Song et al. (2007, JAS), and maximum heating.
!
  do k=2, nz-1
    dz(:,k) = 0.5*( z(:,k+1) - z(:,k-1) )
  enddo
  dz(:,1) = 0.5*( z(:,1) + z(:,2) )
  dz(:,nz) = 0.0  ! may not be used

  do l=1, ncol
    heatmax(l) = maxval( heat(l,kcb(l):kct(l)) )
  enddo

  do k=1, maxval(kct)
    heatdz(:,k) = heat(:,k)*dz(:,k)
  enddo

  tmp1d(:) = 0.0

  do l=1, ncol
    tmp = 0.0
    do k=kcb(l), kct(l)
      if (heat(l,k) > 0.0) then
        tmp = tmp + z(l,k)*heatdz(l,k)
        tmp1d(l) = tmp1d(l) + heatdz(l,k)
      end if
    enddo
    zm(l) = tmp/tmp1d(l)
  enddo

  do l=1, ncol
    tmp = 0.0
    cnt = 0
    do k=kcb(l), kct(l)
      if (heat(l,k) > 0.0) then
        tmp = tmp + ((z(l,k)-zm(l))**2)*heatdz(l,k)
        tmpi = k
        cnt  = cnt + 1
      end if
    enddo
    if (cnt > 1) then
      zd(l) = sqrt(tmp/tmp1d(l))
    else
      zd(l) = (z(l,tmpi+1)-z(l,tmpi-1))/2.
    end if
    if (zd(l) > zm(l))  zd(l) = zm(l) - 10.
    ! zcba must be larger than 0 (i.e., zcba >= 10 m)
  enddo

  ! cloud-top and bottom height
  zcta(:) = zm(:) + zd(:)
  zcba(:) = zm(:) - zd(:)

!
! Determine approximated vertical indices for cloud top and bottom.
!
  ! indices for cloud bottom
  do l=1, ncol
    kcba(l) = minloc(abs(z(l,:) - zcba(l)),1)
  enddo

  ! indices for cloud top
  do l=1, ncol
    kcta(l) = kcba(l) + minloc(abs(z(l,kcba(l)+1:) - zcta(l)),1)
  enddo

!
! Estimate moving velocity of convective cells (See Choi and Chun, 
! 2011, JAS)
!
  eta_3km = 3.2e3 / zt_eta

  if ( present( z_ref ) ) then

    eta(:) = z_ref(:)/zt_eta

    k_3km = nz-1
    do k=2, nz
      if ( eta(k) > eta_3km ) then
        k_3km = k-1  ;  EXIT
      end if
    enddo

    cqx(:) = 0.0  ;  cqy(:) = 0.0
    do k=1, k_3km
      cqx(:) = cqx(:) + u(:,k)*dz(:,k)
      cqy(:) = cqy(:) + v(:,k)*dz(:,k)
    enddo
    cqx(:) = cqx(:)/( 0.5*(z(:,k_3km) + z(:,k_3km+1)) )
    cqy(:) = cqy(:)/( 0.5*(z(:,k_3km) + z(:,k_3km+1)) )

  else

    do l=1, ncol

      k_3km = nz-1
      do k=2, nz
        if ( eta(k) > eta_3km ) then
          k_3km = k-1  ;  EXIT
        end if
      enddo

      cqx(l) = sum(u(l,:k_3km)*dz(l,:k_3km)) /                 &
               ( 0.5*(z(l,k_3km) + z(l,k_3km+1)) )
      cqy(l) = sum(v(l,:k_3km)*dz(l,:k_3km)) /                 &
               ( 0.5*(z(l,k_3km) + z(l,k_3km+1)) )

    enddo

  end if

!
! Obtain the basic-state profile assumed for the analytic solution
! by Song and Chun (2005, JAS).
!
  do l=1, ncol
    k = kcta(l)
    n_q(l) = sum(nbv(l,1:k)*dz(l,1:k))/z(l,k)
    n_ct  (l) = nbv(l,k)
    t_ct  (l) = t  (l,k)
    rho_ct(l) = rho(l,k)
    u_ct  (l) = u  (l,k)
    v_ct  (l) = v  (l,k)
    u_cb  (l) = u  (l,k)
    v_cb  (l) = v  (l,k)
  enddo

  if ( .not. allocated(u_sfc) ) then
    allocate( u_sfc(ncol), v_sfc(ncol) )
    u_sfc(:) = u(:,1)  ;  v_sfc(:) = v(:,1)
  end if
 
  if ( l_diag_znwcq_o ) then
    diag_znwcq(:,1 ) = zcba  (:)
    diag_znwcq(:,2 ) = zcta  (:)
    diag_znwcq(:,3 ) = rho_ct(:)
    diag_znwcq(:,4 ) = n_q   (:)
    diag_znwcq(:,5 ) = n_ct  (:)
    diag_znwcq(:,6 ) = t_ct  (:)
    diag_znwcq(:,7 ) = cqx   (:)
    diag_znwcq(:,8 ) = cqy   (:)
    diag_znwcq(:,9 ) = u_ct  (:)
    diag_znwcq(:,10) = v_ct  (:)
  end if

  RETURN

END SUBROUTINE args_sc05

SUBROUTINE calc_sc05(ncol,nz,shear_ct)

!
! PURPOSE:  To calculate the cloud-top momentum flux spectrum in the 
!           GWDC scheme based on Song and Chun (2005, JAS)
!
! METHOD:
!
!   1) Calculate theta function
!   2) Calculate |X|^2 and the wave-filtering and resonance factor.
!   3) Calculate the cloud-top momentum flux spectrum.
!
! HISTORY:
!

  USE param_gwp
  USE switch_dump,  ONLY: l_diag_znwcq_o

  implicit none

! SUBROUTINE ARGUMENTS

  integer                , intent(in) ::  ncol, nz
  real, dimension(ncol,2), intent(in), optional ::  shear_ct

! LOCAL VARIABLES

  integer, parameter :: dp  = kind(1.0d0)
  integer, parameter :: dpc = kind((1.0d0,1.0d0))

  ! data arrays
  real, dimension(-nc:nc,ncol,nphi) ::  th_ftn, wfrf
  real, dimension(-nc-1:nc+1) ::  x_sq, x2tmp
  real, dimension(ncol) ::  q0sqc
  real, dimension(ncol,nphi) ::  cqh
  real(dp), dimension(ncol) ::  zm, zd
  real(dp), dimension(ncol,nphi) ::  ub_sfc, ub_cb, ub_ct, zcsa
  real(dp), dimension(-nc:nc) ::  c_intr

  ! variables used to calculate |X|^2
  logical ::  singular1, singular2
  integer ::  ic_cl
  real(dp)     ::  shear, ri_cs, n2dn1, mu, ztm_zs, ztm_zb, c2drip2dzd2, &
                   ztqdzb, ztqdzt, ztqs, zalph, zb_mzalph, zs_mzalph,    &
                   ut_mc, lm1, ztsb, ztss, ztsdzb, ztsdzs, ztus, ztut,   &
                   ztudzs, ztudzt
  complex(dpc) ::  aimu, c05p_aimu, c05m_aimu, ailm1,                    &
                   bfac, coef1, coef2, coef3, coef1c, coef2c, coef3c,    &
                   x0, x1, x2, x3, x4, x5, x1c, x2c, x3c, x4c, xa, xac,  &
                   y1, y2, y3, y1c, y2c, y3c, xp, xn, numer, denom

  real, parameter ::  xsq_limit = 100.0

  ! work arrays
  real(dp) ::  tmp1d(ncol), tmp2d(ncol,nphi)
  real     ::  c0, factor, ri_ct, tmp
  integer  ::  ic_u0(-1:2)
  integer  ::  ic,l,iphi   ! loop counters

  ! parameters and constants
  real, parameter ::  g  = 9.80665
  real, parameter ::  cp = 1005.0

  complex(dpc), parameter ::  ai = (0.0d0,1.0d0)

  real(dp), parameter ::  pi = 3.14159265358979323846d0

  real(dp) ::  v_small

  v_small = aimag(exp(ai*pi))
  if (v_small .le. 0.d0) then
    write(6,*) 'v_small should be positive :', v_small  ;  STOP
  end if

  if ( allocated(mfs_ct) )  deallocate( mfs_ct )
  allocate( mfs_ct(-nc:nc,ncol,nphi) )

!
! Obtain the basic-state profile assumed for the analytic solution
! by Song and Chun (2005, JAS).
!
  do iphi=1, nphi
    ub_ct (:,iphi) = dble(u_ct (:)*cosphi(iphi) + v_ct (:)*sinphi(iphi))
    ub_cb (:,iphi) = dble(u_cb (:)*cosphi(iphi) + v_cb (:)*sinphi(iphi))
    ub_sfc(:,iphi) = dble(u_sfc(:)*cosphi(iphi) + v_sfc(:)*sinphi(iphi))
    cqh(:,iphi) = cqx(:)*cosphi(iphi) + cqy(:)*sinphi(iphi)
  enddo

  tmp2d(:,:) = ub_cb(:,:) - ub_sfc(:,:)
  do iphi=1, nphi
  do l=1, ncol
    zcsa(l,iphi) = dble(zcta(l))
    if ( abs(tmp2d(l,iphi)) > 0.1d0 ) then
      tmp1d(l) = (ub_ct(l,iphi)-ub_sfc(l,iphi))/tmp2d(l,iphi)
      if (tmp1d(l) > 0.0d0) then
        zcsa(l,iphi) = tmp1d(l)*dble(zcba(l))
        zcsa(l,iphi) = min( dble(zcta(l)), max( zcsa(l,iphi), dble(zcba(l)) ) )
      end if
    end if
  enddo
  enddo

! will not be used
!     ! revise u_cbottom after determining zcsa
!     do iphi=1, nphi
!       tmp2d(:,iphi) = zcba(:)
!     enddo
!     ub_cb(:,:) = ub_sfc(:,:) + (ub_ct(:,:)-ub_sfc(:,:))*tmp2d(:,:)/   &
!    &             zcsa(:,:)

!
! Calculate theta function that represents the horizontal and temporal 
! structure of diabatic forcing in phase speed and direction domain.
!
  q0sqc(:) = (cfactor*cp*heatmax(:))**2 *   &
                        ( (hscale*tscale) / (32.0*sngl(pi)**1.5) )**2

  c0 = hscale/tscale

  do iphi=1, nphi
  do l=1, ncol
    th_ftn(:,l,iphi) = q0sqc(l)/(1.0+((c_phase(:)-cqh(l,iphi))/c0)**2)
  enddo
  enddo

  zm(:) = 0.5d0*dble(zcta(:) + zcba(:))
  zd(:) = 0.5d0*dble(zcta(:) - zcba(:))

!
! Calculate |X|^2 and the wave-filtering and resonance factor.
!
  do iphi=1, nphi
  do l=1, ncol

    shear = (ub_ct(l,iphi)-ub_sfc(l,iphi))/zcsa(l,iphi)
    ri_cs = (n_q(l)/shear)**2
    n2dn1 = dble(n_ct(l)/n_q(l))

    c_intr(:) = dble(c_phase(:)) - ub_ct(l,iphi)
    ic_cl = minloc(abs(c_intr),1)-(nc+1)
    ic_u0(0) = minloc(abs(dble(c_phase(:))-ub_sfc(l,iphi)),1)-(nc+1)
    if (ic_u0(0) > ic_cl) then
      ic_u0 = ic_u0(0) + (/-1,0,1,2/)
      ic_u0(2) = min(nc+1, ic_u0(2))
    else if (ic_u0(0) < ic_cl) then
      ic_u0 = ic_u0(0) + (/1,0,-1,-2/)
      ic_u0(2) = max(-nc-1, ic_u0(2))
    end if

    if (ri_cs > 0.25d0) then

      x_sq(:) = 0.0   ! -nc-1:nc+1

      if (ri_cs <= 2.5d3) then

        ! |X|^2 for the sheared wind
        bfac        = (-1.0d0,0.0d0)-ai*sign(v_small, shear)
        mu          = sqrt(ri_cs-0.25d0)
        aimu        = ai*mu
        c05p_aimu   = 0.5d0 + aimu
        c05m_aimu   = 0.5d0 - aimu
        ztm_zs      = dble(zcta(l))-zcsa(l,iphi)
        c2drip2dzd2 = 2.0d0/(ri_cs+2.0d0)/zd(l)**2
        coef1       = 0.25d0/aimu - 0.5d0
        ztqdzb      = 2.0d0/zd(l)
        ztqdzt      = -ztqdzb
        ztqs        = 1.0d0 - ((zcsa(l,iphi)-zm(l))/zd(l))**2

        do ic=-nc, nc
          ! do not calculate X2 near two singular points
          ! (1) near critical level (at ic+-1 for dc = 2)
          !     X2 = 0
          ! (2) near ub_sfc (at 2 pts. for dc = 2)
          !     X2 will be interpolated later
          singular1 = ic >= ic_cl-1 .and.                        &
                      ic <= ic_cl+1
          singular2 = ic == ic_u0(0) .or. ic == ic_u0(1)
          if ( singular1 .or. singular2 ) then
            x_sq(ic) = 0.0
! pracitally may not happen
!         else if (dble(c_phase(ic)) == ub_cb(l,iphi)) then
!           x_sq(ic) = 0.0
          else
            zalph     = (dble(c_phase(ic))-ub_sfc(l,iphi))/shear
            zb_mzalph = dble(zcba(l)) - zalph
            zs_mzalph = zcsa(l,iphi) - zalph
            ut_mc     = (-1.0d0)*c_intr(ic)
            lm1       = dble(n_q(l))/abs(ut_mc)
            ailm1     = ai*lm1
            coef2     = ailm1*zs_mzalph
            coef3     = coef2 - c05p_aimu
            coef2     = coef2 - c05m_aimu
            coef1c    = conjg(coef1)
            coef2c    = conjg(coef2)
            coef3c    = conjg(coef3)
            ztsb      = c2drip2dzd2*zb_mzalph*zb_mzalph
            ztss      = c2drip2dzd2*zs_mzalph*zs_mzalph
            ztsdzb    = 2.0d0*c2drip2dzd2*zb_mzalph
            ztsdzs    = 2.0d0*c2drip2dzd2*zs_mzalph
            ztus      = 2.0d0/(lm1*zd(l))**2
            ztut      = ztus
            ztudzs    = 0.0d0
            ztudzt    = 0.0d0
            if (zb_mzalph > 0.0d0) then
              y1  = coef1 /c05m_aimu/zb_mzalph*                          &
                     zb_mzalph**c05m_aimu*                               &
                     ( c05m_aimu*ztsb - (ztqdzb+ztsdzb)*zb_mzalph )
              y1c = conjg(y1)
            else
              y1  = coef1 /c05m_aimu/zb_mzalph*                          &
                     ((-zb_mzalph)*bfac)**c05m_aimu*                     &
                     ( c05m_aimu*ztsb - (ztqdzb+ztsdzb)*zb_mzalph )
              y1c = coef1c/c05p_aimu/zb_mzalph*                          &
                     ((-zb_mzalph)*bfac)**c05p_aimu*                     &
                     ( c05p_aimu*ztsb - (ztqdzb+ztsdzb)*zb_mzalph )
            end if
            y2  = (ztqs+ztss)-(ztsdzs-ztudzs)*zs_mzalph-                 &
                   ailm1*zs_mzalph*(ztss-ztus)
            y3  = ztut + (1.0d0/ailm1)*(ztqdzt+ztudzt)
            y2c = conjg(y2)
            y3c = conjg(y3)
            if (zalph < 0.0d0) then
              x0 = (-zalph)**(2.0d0*aimu)
            else
              x0 = (zalph*bfac)**(2.0d0*aimu)
            end if
            if (zs_mzalph > 0.0d0) then
              x2  = coef2*zs_mzalph**c05p_aimu
              x3  = coef3*zs_mzalph**c05m_aimu
              x2c = conjg(x2)
              x3c = conjg(x3)
            else 
              x2  = coef2 *((-zs_mzalph)*bfac)**c05p_aimu
              x3  = coef3 *((-zs_mzalph)*bfac)**c05m_aimu
              x2c = coef2c*((-zs_mzalph)*bfac)**c05m_aimu
              x3c = coef3c*((-zs_mzalph)*bfac)**c05p_aimu
            end if
            xp  = 1.0d0+sign(1.0d0,ut_mc)*n2dn1
            xn  = 1.0d0-sign(1.0d0,ut_mc)*n2dn1
            x4  = 2.0d0*exp(ailm1*ztm_zs)
            x4c = conjg(x4)
            xa  = x2 /coef2 *(real(y2)-c05m_aimu*(ztss-ztus)) +          &
                   2.0d0*aimu*zs_mzalph*y1c
            xac = x2c/coef2c*(real(y2)-c05p_aimu*(ztss-ztus)) -          &
                   2.0d0*aimu*zs_mzalph*y1
            x5  = 4.0d0*(xa-xac*x0)
            x4  = x4 *(x2 -x3 *x0)
            x4c = x4c*(x3c-x2c*x0)
            numer = x4*y3c + x4c*y3 + x5
            denom = xn*x4 + xp*x4c
            x_sq(ic) = sngl( abs(numer/denom)**2 )
            x_sq(ic) = min(xsq_limit, x_sq(ic))
          end if  ! singular1 .or. singular2
        enddo  ! ic

      else  ! ri_cs > 2.5e3

        ! |X|^2 for the uniform wind
        ! (for numerically stable calculation)
        ztm_zb = zd(l)*2.0d0
        ztqdzb = 2.0d0/zd(l)

        do ic=-nc, nc
          ! do not calculate X2 near critical level (same as the 
          ! shear-wind case)
          singular1 = ic >= ic_cl-1 .and.                        &
                      ic <= ic_cl+1
          if ( singular1 ) then
            x_sq(ic) = 0.0
          else
            ut_mc = (-1.0d0)*c_intr(ic)
            lm1   = dble(n_q(l))/abs(ut_mc)
            ailm1 = ai*lm1
            y1    = 2.0d0/(lm1*zd(l))**2 + (1.0d0/ailm1)*ztqdzb
            y1c   = conjg(y1)
            x1    = 2.0d0*exp(ailm1*dble(zcba(l)))
            x1c   = conjg(x1)
            x2    = 2.0d0*exp(ailm1*ztm_zb)
            x2c   = conjg(x2)
            xp    = 1.0d0+sign(1.0d0,ut_mc)*n2dn1
            xn    = 1.0d0-sign(1.0d0,ut_mc)*n2dn1
            numer = x1*x2*(2.0d0*y1-x2c*y1c)+x1c*x2c*(2.0d0*y1c-x2*y1)
            denom = 2.0d0*(xn*x1*x2+xp*x1c*x2c)
            x_sq(ic) = sngl( abs(numer/denom)**2 )
            x_sq(ic) = min(xsq_limit, x_sq(ic))
          end if  ! singular1
        enddo  ! ic

      end if  ! ri_cs <= 2.5e3

      ! prevent from overestimation of the spectral "density" due to
      ! the coarse resolution (dc = 2), where peaks are sharp (i.e. 
      ! may not dense)
      ! x2 = a*sqrt(x2/a), if x2/a > 1, where a = x2(i-1)+x2(i+1)
      x2tmp(:) = x_sq(:)
      do ic=-nc, nc
        tmp = x2tmp(ic-1) + x2tmp(ic+1)
        if (x_sq(ic) > tmp) then
          x_sq(ic) = sqrt(tmp*x_sq(ic))
        end if
      enddo

      ! interpolate x2 near ub_sfc
      if (ic_u0(0) /= ic_cl) then
        if ( abs(ic_u0(0)-ic_cl) > 1 ) then
          x_sq(ic_u0(0)) = (x_sq(ic_u0(-1))*2.+x_sq(ic_u0(2)))/3.
        else
          x_sq(ic_u0(0)) = 0.0
        end if
        x_sq(ic_u0(1)) = 0.5*(x_sq(ic_u0(0))+x_sq(ic_u0(2)))
      end if

      ! Wave-filtering and resonance factor
      c_intr(ic_cl) = -999.d0   ! prevent from 0/0 in case
      wfrf(:,l,iphi) = x_sq(-nc:nc)/sngl(c_intr(:))*n_ct(l)

    else  ! ri_cs <= 0.25

      wfrf(:,l,iphi) = 0.0

    end if  ! ri_cs > 0.25

  enddo  ! l
  enddo  ! iphi

!
! Calculate cloud-top momentum flux spectrum.
!
  tmp = (2.0*(2.0*sngl(pi))**3)/(ah*lt)*(g/cp)**2
  do iphi=1, nphi
  do l=1, ncol
    factor = rho_ct(l)*tmp/(t_ct(l)*n_q(l)**2)**2
    mfs_ct(:,l,iphi) = factor*wfrf(:,l,iphi)*th_ftn(:,l,iphi)
  enddo
  enddo

  if ( present(shear_ct) ) then
    do iphi=1, nphi
    do l=1, ncol
      tmp = (shear_ct(l,1)*cosphi(iphi) + shear_ct(l,2)*sinphi(iphi))**2
      if (tmp < 1.e-16)  CYCLE
      ri_ct = n_ct(l)**2 / tmp
      if (ri_ct < 0.25)  mfs_ct(:,l,iphi) = 0.0
    enddo
    enddo
  end if

  RETURN

END SUBROUTINE calc_sc05

END module mflx_ctop_sc05

