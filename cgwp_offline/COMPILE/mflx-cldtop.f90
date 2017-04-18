MODULE mflx_cldtop

  implicit none

  integer, dimension(:)    , allocatable ::  kcta
  real   , dimension(:,:,:), allocatable ::  mfs_ct
  real   , dimension(:,:)  , allocatable ::  diag_znwcq

  real, dimension(:,:,:), allocatable ::  diag_spec_ctop_0

  real, dimension(:), allocatable ::  mflx_ct_0_east, mflx_ct_0_west,    &
                                      mflx_ct_0_north, mflx_ct_0_south


CONTAINS


SUBROUTINE gw_ctop(ncol,nz,  &
     &   eta_theta_levels,                                              &
     &   u_flev,v_flev,rho_flev,z_flev,nbv_flev,         &
     &   heat,heatmax,u_sfc,v_sfc,kcb,kct,                            &
     &   t_flev                                      &
     &   )

!
! PURPOSE:  To calculate the cloud-top momentum flux spectrum in the 
!           GWDC scheme based on Song and Chun (2005, JAS) with
!           the modification by Choi and Chun (2011, JAS)
!
! METHOD:
!
!   1) Determine the heating center height and half depth, and maximum 
!      heating rate.
!   2) Estimate moving velocity of convective cells.
!   3) Calculate theta function
!   4) Calculate |X|^2 and the wave-filtering and resonance factor.
!   5) Calculate the cloud-top momentum flux spectrum.
!
! HISTORY:
!

  USE param_gwp
  USE switch_dump,  ONLY: l_diag_znwcq_o

  implicit none

! SUBROUTINE ARGUMENTS

  integer, intent(in) ::  ncol, nz

  ! data arrays
  real, dimension(nz)     , intent(in) ::  eta_theta_levels
  real, dimension(ncol,nz), intent(in) ::  u_flev, v_flev, rho_flev,     &
                                           nbv_flev, t_flev, heat, z_flev
  real, dimension(ncol)   , intent(in) ::  heatmax, u_sfc, v_sfc
  integer, dimension(ncol), intent(in) ::  kcb, kct

! LOCAL VARIABLES

  ! data arrays
  real, dimension(ncol) ::  q0sqc, zm, zd, zcta, zcba, n_q, n_ct, t_ct,  &
                            rho_ct, cqx, cqy
  integer, dimension(ncol,nphi) ::  ic_cl
  integer, dimension(ncol) ::  kcba
  real, dimension(ncol,nphi) ::  cqh, zcsa, ub_sfc, ub_cb, ub_ct, ri_ct
  real, dimension(-nc:nc) ::  c_intr, mfs_ct0
  real, dimension(-nc-1:nc+1) ::  x2tmp, x_sq
  real, dimension(-nc:nc,ncol,nphi) ::  th_ftn, wfrf

  ! variables used to calculate |X|^2
  logical ::  singular1, singular2
  real    ::  shear, ri_cs, n2dn1, mu, ztm_zs, ztm_zb, c2drip2dzd2,      &
              ztqdzb, ztqdzt, ztqs, zalph, zb_mzalph, zs_mzalph, ut_mc,  &
              lm1, ztsb, ztss, ztsdzb, ztsdzs, ztus, ztut, ztudzs, ztudzt
  complex ::  aimu, c05p_aimu, c05m_aimu, ailm1,                         &
              bfac, coef1, coef2, coef3, coef1c, coef2c, coef3c,         &
              x0, x1, x2, x3, x4, x5, x1c, x2c, x3c, x4c, xa, xac,       &
              y1, y2, y3, y1c, y2c, y3c, xp, xn, numer, denom

  real, parameter ::  xsq_limit = 100.0

  ! work arrays
  real ::  tmp, eta_3km, factor, tmp1d(ncol), tmp2d(ncol,nphi)
  real, dimension(nz) ::  d_eta
  real, dimension(ncol,nz) ::  heatdh
  integer ::  tmpi, cnt, k_3km, ic_u0(-1:2)
  integer ::  ic,k,l,iphi   ! loop counters

  ! parameters and constants
  real, parameter ::  g = 9.80665
  real, parameter ::  cp = 1005.

  complex, parameter ::  ai = (0.0,1.0)

  include 'c_math.inc'

  real ::  c0, v_small

  v_small = aimag(exp(ai*pi))
print*,v_small, 'should be positive'

  if ( allocated(kcta) )  deallocate( kcta )
  allocate( kcta(ncol) )

  if ( allocated(mfs_ct) )  deallocate( mfs_ct )
  allocate( mfs_ct(-nc:nc,ncol,nphi) )

  if ( l_diag_znwcq_o ) then
    if ( allocated(diag_znwcq) )  deallocate( diag_znwcq )
    allocate( diag_znwcq(ncol,10) )
  end if

!
! Determine heating center height and half depth,
! following Song et al. (2007, JAS), and maximum heating.
!
      do k=2, nz-1
        d_eta(k) = 0.5*( eta_theta_levels(k+1) - eta_theta_levels(k-1) )
      enddo
      d_eta(1) = 0.5*( eta_theta_levels(1) + eta_theta_levels(2) )
      d_eta(nz) = 0.0  ! may not be used

      do k=minval(kcb), maxval(kct)
      do l=1, ncol
        heatdh(l,k) = heat(l,k)*d_eta(k)
      enddo
      enddo

      tmp1d(:) = 0.0

      do l=1, ncol
        tmp = 0.0
        do k=kcb(l), kct(l)
          if (heat(l,k) > 0.0) then
            tmp = tmp + z_flev(l,k)*heatdh(l,k)
            tmp1d(l) = tmp1d(l) + heatdh(l,k)
          end if
        enddo
        zm(l) = tmp/tmp1d(l)
      enddo

      do l=1, ncol
        tmp = 0.0
        cnt  = 0
        do k=kcb(l), kct(l)
          if (heat(l,k) > 0.0) then
            tmp = tmp + ((z_flev(l,k)-zm(l))**2)*heatdh(l,k)
            tmpi = k
            cnt  = cnt + 1
          end if
        enddo
        if (cnt > 1) then
          zd(l) = sqrt(tmp/tmp1d(l))
        else
          zd(l) = (z_flev(l,tmpi+1)-z_flev(l,tmpi-1))/2.
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
        if (zcba(l) <= z_flev(l,1)) then
          kcba(l) = 1
        else
          do k=1, nz-1
            if ( zcba(l) >  z_flev(l,k  ) .and.                   &
     &           zcba(l) <= z_flev(l,k+1) ) then
              if ( (zcba(l)-z_flev(l,k)) >                        &
     &             (z_flev(l,k+1)-zcba(l)) ) then
                kcba(l) = k+1
              else
                kcba(l) = k
              end if
              EXIT
            end if
          enddo
        end if
      enddo

      ! indices for cloud top
      kcta(:) = kcba(:)+1
      do l=1, ncol
      do k=kcba(l)+1, nz-1
        if ( zcta(l) >  z_flev(l,k  ) .and.                       &
     &       zcta(l) <= z_flev(l,k+1) ) then
          if ( (zcta(l)-z_flev(l,k)) >                            &
     &         (z_flev(l,k+1)-zcta(l)) ) then
            kcta(l) = k+1
          else
            kcta(l) = k
          end if
          EXIT
        end if
      enddo
      enddo

!
! Estimate moving velocity of convective cells (See Choi and Chun, 
! 2011, JAS)
!
      eta_3km = 3.2e3 / z_flev(1,nz)

      k_3km = nz-1
      do k=2, nz
        if ( eta_theta_levels(k) > eta_3km ) then
          k_3km = k-1
          EXIT
        end if
      enddo

      eta_3km = 0.5*(eta_theta_levels(k_3km)+eta_theta_levels(k_3km+1))

      cqx(:) = 0.0
      cqy(:) = 0.0

      do k=1, k_3km
      do l=1, ncol
        cqx(l) = cqx(l) + u_flev(l,k)*d_eta(k)
        cqy(l) = cqy(l) + v_flev(l,k)*d_eta(k)
      enddo
      enddo
      cqx(:) = cqx(:)/eta_3km
      cqy(:) = cqy(:)/eta_3km

      do iphi=1, nphi
      do l=1, ncol
        cqh(l,iphi) = cqx(l)*cosphi(iphi) + cqy(l)*sinphi(iphi)
      enddo
      enddo

!
! Obtain the basic-state profile assumed for the analytic solution
! by Song and Chun (2005, JAS).
!
      n_q(:) = 0.0
      do l=1, ncol
        do k=1, kcta(l)
          n_q(l) = n_q(l) + nbv_flev(l,k)*d_eta(k)
        enddo
        n_q(l) = n_q(l) / eta_theta_levels(kcta(l))
        k = kcta(l)
        n_ct  (l) = nbv_flev     (l,k)
        t_ct  (l) = t_flev       (l,k)
        rho_ct(l) = rho_flev(l,k)
      enddo

      do iphi=1, nphi
      do l=1, ncol
        ub_sfc(l,iphi) = u_sfc(l)*cosphi(iphi) + v_sfc(l)*sinphi(iphi)
        ub_ct(l,iphi) = u_flev(l,kcta(l))*cosphi(iphi) +                 &
                        v_flev(l,kcta(l))*sinphi(iphi)
        ub_cb(l,iphi) = u_flev(l,kcba(l))*cosphi(iphi) +                 &
                        v_flev(l,kcba(l))*sinphi(iphi)
      enddo
      enddo

      do iphi=1, nphi
      do l=1, ncol
        tmpi = kcta(l)
        tmp = ( (u_flev(l,tmpi+1) - u_flev(l,tmpi-1))*cosphi(iphi) +     &
                (v_flev(l,tmpi+1) - v_flev(l,tmpi-1))*sinphi(iphi) ) /   &
              ( z_flev(l,tmpi+1) - z_flev(l,tmpi-1) )
        tmp = tmp*tmp
        if (tmp > 1.e-16) then
          ri_ct(l,iphi) = n_ct(l)**2 / tmp
        else
          ri_ct(l,iphi) = 1.e10
        end if
      enddo
      enddo

      tmp2d(:,:) = ub_cb(:,:) - ub_sfc(:,:)
      do iphi=1, nphi
      do l=1, ncol
        zcsa(l,iphi) = zcta(l)
        if ( abs(tmp2d(l,iphi)) > 0.1 ) then
          tmp1d(l) = (ub_ct(l,iphi)-ub_sfc(l,iphi))/tmp2d(l,iphi)
          if (tmp1d(l) > 0.0) then
            zcsa(l,iphi) = tmp1d(l)*zcba(l)
            zcsa(l,iphi) = min( zcta(l), max( zcsa(l,iphi), zcba(l) ) )
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
      q0sqc(:) = (cfactor*cp*heatmax(:))**2

      tmp = ( (hscale*tscale) / (32.*pi**1.5) )**2
      q0sqc(:) = q0sqc(:) * tmp

      c0 = hscale/tscale

      do iphi=1, nphi
      do l=1, ncol
      do ic=-nc, nc
        th_ftn(ic,l,iphi) = q0sqc(l) /                                  &
     &                     (1.+((c_phase(ic)-cqh(l,iphi))/c0)**2)
      enddo
      enddo
      enddo

!
! Calculate |X|^2 and the wave-filtering and resonance factor.
!
      do iphi=1, nphi
      do l=1, ncol

        shear = (ub_ct(l,iphi)-ub_sfc(l,iphi))/zcsa(l,iphi)
        ri_cs = (n_q(l)/shear)**2
        n2dn1 = n_ct(l)/n_q(l)

        c_intr(:) = c_phase(:) - ub_ct(l,iphi)
        ic_cl(l,iphi) = minloc(abs(c_intr),1)-(nc+1)
        ic_u0(0) = minloc(abs(c_phase(:)-ub_sfc(l,iphi)),1)-(nc+1)
        if (ic_u0(0) > ic_cl(l,iphi)) then
          ic_u0 = ic_u0(0) + (/-1,0,1,2/)
          ic_u0(2) = min(nc+1, ic_u0(2))
        else if (ic_u0(0) < ic_cl(l,iphi)) then
          ic_u0 = ic_u0(0) + (/1,0,-1,-2/)
          ic_u0(2) = max(-nc-1, ic_u0(2))
        end if

        if (ri_cs > 0.25) then

          x_sq(:) = 0.0   ! -nc-1:nc+1

          if (ri_cs <= 2.5e3) then

            ! |X|^2 for the sheared wind
            bfac        = (-1.0,0.0)-ai*sign(v_small, shear)
            mu          = sqrt(ri_cs-0.25)
            aimu        = ai*mu
            c05p_aimu   = 0.5 + aimu
            c05m_aimu   = 0.5 - aimu
            ztm_zs      = zcta(l)-zcsa(l,iphi)
            c2drip2dzd2 = 2.0/(ri_cs+2.0)/zd(l)**2
            coef1       = 0.25/aimu - 0.5
            ztqdzb      = 2.0/zd(l)
            ztqdzt      = -ztqdzb
            ztqs        = 1.0 - ((zcsa(l,iphi)-zm(l))/zd(l))**2

            do ic=-nc, nc
              ! do not calculate X2 near two singular points
              ! (1) near critical level (at ic+-1 for dc = 2)
              !     X2 = 0
              ! (2) near u_sfc (at 2 pts. for dc = 2)
              !     X2 will be interpolated later
              singular1 = ic >= ic_cl(l,iphi)-1 .and.                   &
     &                    ic <= ic_cl(l,iphi)+1
              singular2 = ic == ic_u0(0) .or. ic == ic_u0(1)
              if ( singular1 .or. singular2 ) then
                x_sq(ic) = 0.0
! pracitally may not happen
!             else if (c_phase(ic) == ub_cb(l,iphi)) then
!               x_sq(ic) = 0.0
              else
                zalph     = (c_phase(ic)-ub_sfc(l,iphi))/shear
                zb_mzalph = zcba(l) - zalph
                zs_mzalph = zcsa(l,iphi) - zalph
                ut_mc     = (-1.0)*c_intr(ic)
                lm1       = n_q(l)/abs(ut_mc)
                ailm1     = ai*lm1
                coef2     = ailm1*zs_mzalph
                coef3     = coef2 - c05p_aimu
                coef2     = coef2 - c05m_aimu
                coef1c    = conjg(coef1)
                coef2c    = conjg(coef2)
                coef3c    = conjg(coef3)
                ztsb      = c2drip2dzd2*zb_mzalph*zb_mzalph
                ztss      = c2drip2dzd2*zs_mzalph*zs_mzalph
                ztsdzb    = 2.0*c2drip2dzd2*zb_mzalph
                ztsdzs    = 2.0*c2drip2dzd2*zs_mzalph
                ztus      = 2.0/(lm1*zd(l))**2
                ztut      = ztus
                ztudzs    = 0.0
                ztudzt    = 0.0
                if (zb_mzalph > 0.0) then
                  y1  = coef1 /c05m_aimu/zb_mzalph*                     &
     &                   zb_mzalph**c05m_aimu*                          &
     &                   ( c05m_aimu*ztsb - (ztqdzb+ztsdzb)*zb_mzalph )
                  y1c = conjg(y1)
                else
                  y1  = coef1 /c05m_aimu/zb_mzalph*                     &
     &                   ((-zb_mzalph)*bfac)**c05m_aimu*                &
     &                   ( c05m_aimu*ztsb - (ztqdzb+ztsdzb)*zb_mzalph )
                  y1c = coef1c/c05p_aimu/zb_mzalph*                     &
     &                   ((-zb_mzalph)*bfac)**c05p_aimu*                &
     &                   ( c05p_aimu*ztsb - (ztqdzb+ztsdzb)*zb_mzalph )
                end if
                y2  = (ztqs+ztss)-(ztsdzs-ztudzs)*zs_mzalph-            &
     &                 ailm1*zs_mzalph*(ztss-ztus)
                y3  = ztut + (1.0/ailm1)*(ztqdzt+ztudzt)
                y2c = conjg(y2)
                y3c = conjg(y3)
                if (zalph < 0.0) then
                  x0 = (-zalph)**(2.0*aimu)
                else
                  x0 = (zalph*bfac)**(2.0*aimu)
                end if
                if (zs_mzalph > 0.0) then
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
                xp  = 1.0+sign(1.0,ut_mc)*n2dn1
                xn  = 1.0-sign(1.0,ut_mc)*n2dn1
                x4  = 2.0*exp(ailm1*ztm_zs)
                x4c = conjg(x4)
                xa  = x2 /coef2 *(real(y2)-c05m_aimu*(ztss-ztus)) +     &
     &                 2.0*aimu*zs_mzalph*y1c
                xac = x2c/coef2c*(real(y2)-c05p_aimu*(ztss-ztus)) -     &
     &                 2.0*aimu*zs_mzalph*y1
                x5  = 4.0*(xa-xac*x0)
                x4  = x4 *(x2 -x3 *x0)
                x4c = x4c*(x3c-x2c*x0)
                numer = x4*y3c + x4c*y3 + x5
                denom = xn*x4 + xp*x4c
                x_sq(ic) = abs(numer/denom)**2
                x_sq(ic) = min(xsq_limit, x_sq(ic))
              end if  ! singular1 .or. singular2
            enddo  ! ic

          else  ! ri_cs > 2.5e3

            ! |X|^2 for the uniform wind
            ! (for numerically stable calculation)
            ztm_zb = zcta(l)-zcba(l)
            ztqdzb = 2.0/zd(l)

            do ic=-nc, nc
              ! do not calculate X2 near critical level (same as the 
              ! shear-wind case)
              singular1 = ic >= ic_cl(l,iphi)-1 .and.                   &
     &                    ic <= ic_cl(l,iphi)+1
              if ( singular1 ) then
                x_sq(ic) = 0.0
              else
                ut_mc = (-1.0)*c_intr(ic)
                lm1   = n_q (l)/abs(ut_mc)
                ailm1 = ai*lm1
                y1    = 2.0/(lm1*zd(l))**2 + (1.0/ailm1)*ztqdzb
                y1c   = conjg(y1)
                x1    = 2.0*exp(ailm1*zcba(l))
                x1c   = conjg(x1)
                x2    = 2.0*exp(ailm1*ztm_zb)
                x2c   = conjg(x2)
                xp    = 1.0+sign(1.0,ut_mc)*n2dn1
                xn    = 1.0-sign(1.0,ut_mc)*n2dn1
                numer = x1*x2*(2.0*y1-x2c*y1c)+x1c*x2c*(2.0*y1c-x2*y1)
                denom = 2.0*(xn*x1*x2+xp*x1c*x2c)
                x_sq(ic) = abs(numer/denom)**2
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

          ! interpolate x2 near u_sfc
          if (ic_u0(0) /= ic_cl(l,iphi)) then
            if ( abs(ic_u0(0)-ic_cl(l,iphi)) > 1 ) then
              x_sq(ic_u0(0)) = (x_sq(ic_u0(-1))*2.+x_sq(ic_u0(2)))/3.
            else
              x_sq(ic_u0(0)) = 0.0
            end if
            x_sq(ic_u0(1)) = 0.5*(x_sq(ic_u0(0))+x_sq(ic_u0(2)))
          end if

          ! Wave-filtering and resonance factor
          c_intr(ic_cl(l,iphi)) = -999.   ! prevent from 0/0 in case
          wfrf(:,l,iphi) = x_sq(-nc:nc)/c_intr(:)*n_ct(l)

        else  ! ri_cs <= 0.25

          wfrf(:,l,iphi) = 0.0

        end if  ! ri_cs > 0.25

      enddo  ! l
      enddo  ! iphi

!
! Calculate cloud-top momentum flux spectrum.
!
      tmp = (2.0*(2.0*pi)**3)/(ah*lt)*(g/cp)**2
      do iphi=1, nphi
      do l=1, ncol
        if (ri_ct(l,iphi) > 0.25) then
          factor = rho_ct(l)*tmp/(t_ct(l)*n_q(l)**2)**2
          mfs_ct(:,l,iphi) = factor*wfrf(:,l,iphi)*th_ftn(:,l,iphi)
        else
          mfs_ct(:,l,iphi) = 0.0
        end if
      enddo
      enddo

!
! Diag (z_cb,z_ct,rho_ct,n1,n_ct,t_ct,cqx,cqy, u_ct,v_ct)
!
      diag_znwcq(:,1 ) = zcba  (:)
      diag_znwcq(:,2 ) = zcta  (:)
      diag_znwcq(:,3 ) = rho_ct(:)
      diag_znwcq(:,4 ) = n_q   (:)
      diag_znwcq(:,5 ) = n_ct  (:)
      diag_znwcq(:,6 ) = t_ct  (:)
      diag_znwcq(:,7 ) = cqx   (:)
      diag_znwcq(:,8 ) = cqy   (:)
      do l=1, ncol
        diag_znwcq(l,9 ) = u_flev(l,kcta(l))
        diag_znwcq(l,10) = v_flev(l,kcta(l))
      enddo


      RETURN

      END SUBROUTINE gw_ctop

END module mflx_cldtop
