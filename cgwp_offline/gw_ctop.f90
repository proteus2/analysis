

!
!  subroutine gw_ctop to calculate cloud-top GW momentum flux spectrum
!
      SUBROUTINE gw_ctop(row_length,rows,nrows,levels,wet_levels,       &
     &   ncol,nphi,nc,hscale,tscale,cfactor,lt,ah,igwdc,jgwdc,          &
     &   eta_theta_levels,                                              &
     &   ut_gwdc,vt_gwdc,rho_thlev_gwdc,r_thlev_gwdc,nbv_thlev,         &
     &   heat_gwdc,heatmax,ub_thlev,c_phase,c_m05dc,kcb,kct,            &
     &   t_thlev,cosphi,sinphi,dc,                                      &
     &   kcta_thlev,mfs_ct,ic_cl,diag_znwcq                             &
     &   )

!
! PURPOSE:  To calculate the cloud-top momentum flux spectrum in the 
!           GWDC scheme based on Song and Chun (2005, JAS) with
!           the modification by Choi and Chun (2011, JAS)
!
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
!
! HISTORY:
!
!     Date    Comment
!   --------  -------
!   22/04/03  First written.                                  I.-S. Song
!   15/04/05  Bugs on AIX (IBM) are fixed.                    I.-S. Song
!   16/05/09  Reconstruction of the code for the UM.           Y.-H. Kim
!   27/05/10  Modification in the estimation of moving velocities of
!             convective cells based on Choi and Chun (2011, JAS).
!                                                              Y.-H. Kim
!
! CODE DESCRIPTION:
!
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.
!

      implicit none

!
! SUBROUTINE ARGUMENTS
!
!   model dimensions

      integer, intent(in) ::                                            &
     &  row_length                                                      &
                           ! number of points per row
     &, rows                                                            &
                           ! number of rows on theta and u grids
     &, nrows                                                           &
                           ! number of rows on v grid
     &, levels                                                          &
                           ! number of model levels
     &, wet_levels                                                      &
                           ! number of wet model levels
     &, ncol               ! number of columns GWDC is calculated at.

!   GWDC parameters

      integer, intent(in) ::                                            &
     &  nphi                                                            &
     &, nc

      real, intent(in) ::                                               &
     &  hscale                                                          &
     &, tscale                                                          &
     &, cfactor                                                         &
     &, lt                                                              &
     &, ah

!   data arrays

      real, intent(in) ::                                               &
     &  eta_theta_levels(0:levels)                                      &
     &, ut_gwdc(ncol, 0:levels)                                         &
     &, vt_gwdc(ncol, 0:levels)                                         &
     &, rho_thlev_gwdc(ncol, levels)                                    &
     &, nbv_thlev     (ncol, levels)                                    &
     &, t_thlev   (ncol, wet_levels)                                    &
     &, heat_gwdc (ncol, wet_levels)                                    &
     &, heatmax(ncol)                                                   &
     &, ub_thlev(ncol, 0:levels, nphi)                                  &
     &, c_phase(-nc:nc)                                                 &
     &, c_m05dc(-nc:nc+1)

      integer, intent(in) ::                                            &
     &  igwdc(ncol)                                                     &
     &, jgwdc(ncol)                                                     &
     &, kcb(ncol)                                                       &
     &, kct(ncol) 

      real, intent(in) ::                                               &
     &  dc                                                              &
     &, cosphi(nphi)                                                    &
     &, sinphi(nphi)

!   coordinate arrays

      real, intent(in) ::                                               &
     &  r_thlev_gwdc(ncol, levels)

!   output variables

      integer, intent(out) ::                                           &
     &  kcta_thlev(ncol)                                                &
     &, ic_cl(ncol, nphi)

      real, intent(out) ::                                              &
     &  mfs_ct(-nc:nc, ncol, nphi)                                      &
     &, diag_znwcq(ncol, levels)

!
! global VARIABLES
!
!*L------------------COMDECK C_G----------------------------------------
! G IS MEAN ACCEL DUE TO GRAVITY AT EARTH'S SURFACE

      Real, Parameter :: G = 9.80665

!*----------------------------------------------------------------------
!*L------------------COMDECK C_PI---------------------------------------
!LL
!LL 4.0 19/09/95  New value for PI. Old value incorrect
!LL               from 12th decimal place. D. Robinson
!LL 5.1 7/03/00   Fixed/Free format P.Selwood
!LL

      ! Pi
      Real, Parameter :: Pi                 = 3.14159265358979323846

      ! Conversion factor degrees to radians
      Real, Parameter :: Pi_Over_180        = Pi/180.0

      ! Conversion factor radians to degrees
      Real, Parameter :: Recip_Pi_Over_180  = 180.0/Pi

!*----------------------------------------------------------------------
!*L------------------COMDECK C_R_CP-------------------------------------
! History:
! Version  Date      Comment.
!  5.0  07/05/99  Add variable P_zero for consistency with
!                 conversion to C-P 'C' dynamics grid. R. Rawlins
!  5.1  07/03/00  Fixed/Free format conversion   P. Selwood

! R IS GAS CONSTANT FOR DRY AIR
! CP IS SPECIFIC HEAT OF DRY AIR AT CONSTANT PRESSURE
! PREF IS REFERENCE SURFACE PRESSURE

      Real, Parameter  :: R      = 287.05
      Real, Parameter  :: CP     = 1005.
      Real, Parameter  :: Kappa  = R/CP
      Real, Parameter  :: Pref   = 100000.

      ! Reference surface pressure = PREF
      Real, Parameter  :: P_zero = Pref
      Real, Parameter  :: sclht  = 6.8E+03  ! Scale Height H
!*----------------------------------------------------------------------

! 
! LOCAL VARIABLES
!
!   data arrays

      real ::                                                           &
     &  q0sqc  (ncol)                                                   &
     &, zm     (ncol)                                                   &
     &, zd     (ncol)                                                   &
     &, zcta   (ncol)                                                   &
     &, zcba   (ncol)                                                   &
     &, n_q    (ncol)                                                   &
     &, n_ct   (ncol)                                                   &
     &, t_ct   (ncol)                                                   &
     &, rho_ct (ncol)                                                   &
     &, cqh    (ncol, nphi)                                             &
     &, zcsa   (ncol, nphi)                                             &
     &, ub_sfc (ncol, nphi)                                             &
     &, ub_cb  (ncol, nphi)                                             &
     &, ub_ct  (ncol, nphi)                                             &
     &, ri_ct  (ncol, nphi)                                             &
     &, c_int  (-nc:nc)                                                 &
     &, x2tmp  (-nc-1:nc+1)                                             &
     &, th_ftn (-nc:nc, ncol, nphi)                                     &
     &, x_sq   (-nc-1:nc+1)                                             &
     &, wfrf   (-nc:nc, ncol, nphi)                                     &
     &, mfs_ct0(-nc:nc)


      integer ::                                                        &
     &  kcta_rholev(ncol)                                               &
     &, kcba_rholev(ncol)                                               &
     &, kcba_thlev (ncol)

      real ::                                                           &
     &  c0

      real ::                                                           &
     &  v_small

!   variables used to calculate |X|^2
      logical ::                                                        &
     &  singular1, singular2

      real ::                                                           &
     &  shear, ri_cs, n2dn1, mu, ztm_zs, ztm_zb, c2drip2dzd2            &
     &, ztqdzb, ztqdzt, ztqs                                            &
     &, zalph, zb_mzalph, zs_mzalph, ut_mc, lm1                         &
     &, ztsb, ztss, ztsdzb, ztsdzs, ztus, ztut, ztudzs, ztudzt

      complex ::                                                        &
     &  aimu, c05p_aimu, c05m_aimu, ailm1                               &
     &, bfac, coef1, coef2, coef3, coef1c, coef2c, coef3c               &
     &, x0, x1, x2, x3, x4, x5, x1c, x2c, x3c, x4c, xa, xac             &
     &, y1, y2, y3, y1c, y2c, y3c                                       &
     &, xp, xn, numer, denom

      real, parameter ::                                                &
     &  xsq_limit = 100.0

      complex, parameter ::                                             &
     &  ai = (0.0,1.0)

!   work arrays

      real ::                                                           &
     &  eta_3km                                                         &
     &, temp                                                            &
     &, factor                                                          &
     &, heatdh(ncol,wet_levels)                                         &
     &, cqx  (ncol)                                                     &
     &, cqy  (ncol)                                                     &
     &, tem1d(ncol)                                                     &
     &, tem2d(ncol,nphi)                                                &
     &, d_eta(levels)

      integer ::                                                        &
     &  temi, cnt                                                       &
     &, k_3km                                                           &
     &, ic_u0(-1:2)

      integer ::  ic,k,l,iphi   ! loop counters

      v_small = aimag(exp(ai*pi))


!
! Determine heating center height and half depth,
! following Song et al. (2007, JAS), and maximum heating.
!
      do k=2, levels-1
        d_eta(k) = 0.5*( eta_theta_levels(k+1) - eta_theta_levels(k-1) )
      enddo
      d_eta(1) = 0.5*( eta_theta_levels(1) + eta_theta_levels(2) )
      d_eta(levels) = 0.0  ! may not be used

      do k=minval(kcb), maxval(kct)
      do l=1, ncol
        heatdh(l,k) = heat_gwdc(l,k)*d_eta(k)
      enddo
      enddo

      tem1d(:) = 0.0

      do l=1, ncol
        temp = 0.0
        do k=kcb(l), kct(l)
          if (heat_gwdc(l,k) > 0.0) then
            temp = temp + r_thlev_gwdc(l,k)*heatdh(l,k)
            tem1d(l) = tem1d(l) + heatdh(l,k)
          end if
        enddo
        zm(l) = temp/tem1d(l)
      enddo

      do l=1, ncol
        temp = 0.0
        cnt  = 0
        do k=kcb(l), kct(l)
          if (heat_gwdc(l,k) > 0.0) then
            temp = temp + ((r_thlev_gwdc(l,k)-zm(l))**2)*heatdh(l,k)
            temi = k
            cnt  = cnt + 1
          end if
        enddo
        if (cnt > 1) then
          zd(l) = sqrt(temp/tem1d(l))
        else
          zd(l) = (r_thlev_gwdc(l,temi+1)-r_thlev_gwdc(l,temi-1))/2.
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
        if (zcba(l) <= r_thlev_gwdc(l,1)) then
          kcba_thlev(l) = 1
        else
          do k=1, levels-1
            if ( zcba(l) >  r_thlev_gwdc(l,k  ) .and.                   &
     &           zcba(l) <= r_thlev_gwdc(l,k+1) ) then
              if ( (zcba(l)-r_thlev_gwdc(l,k)) >                        &
     &             (r_thlev_gwdc(l,k+1)-zcba(l)) ) then
                kcba_thlev(l) = k+1
              else
                kcba_thlev(l) = k
              end if
              EXIT
            end if
          enddo
        end if
      enddo

      ! indices for cloud top
      kcta_thlev(:) = kcba_thlev(:)+1
      do l=1, ncol
      do k=kcba_thlev(l)+1, levels-1
        if ( zcta(l) >  r_thlev_gwdc(l,k  ) .and.                       &
     &       zcta(l) <= r_thlev_gwdc(l,k+1) ) then
          if ( (zcta(l)-r_thlev_gwdc(l,k)) >                            &
     &         (r_thlev_gwdc(l,k+1)-zcta(l)) ) then
            kcta_thlev(l) = k+1
          else
            kcta_thlev(l) = k
          end if
          EXIT
        end if
      enddo
      enddo

!
! Estimate moving velocity of convective cells (See Choi and Chun, 
! 2011, JAS)
!
      eta_3km = 3.2e3 / r_thlev_gwdc(1,levels)

      k_3km = levels-1
      do k=2, levels
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
        cqx(l) = cqx(l) + ut_gwdc(l,k)*d_eta(k)
        cqy(l) = cqy(l) + vt_gwdc(l,k)*d_eta(k)
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
        do k=1, kcta_thlev(l)
          n_q(l) = n_q(l) + nbv_thlev(l,k)*d_eta(k)
        enddo
        n_q(l) = n_q(l) / eta_theta_levels(kcta_thlev(l))
        k = kcta_thlev(l)
        n_ct  (l) = nbv_thlev     (l,k)
        t_ct  (l) = t_thlev       (l,k)
        rho_ct(l) = rho_thlev_gwdc(l,k)
      enddo

      do iphi=1, nphi
      do l=1, ncol
        ub_sfc(l,iphi) = ub_thlev(l,0            ,iphi)
        ub_cb (l,iphi) = ub_thlev(l,kcba_thlev(l),iphi)
        ub_ct (l,iphi) = ub_thlev(l,kcta_thlev(l),iphi)
      enddo
      enddo

      do iphi=1, nphi
      do l=1, ncol
        temi = kcta_thlev(l)
        temp = ( ub_thlev(l,temi+1,iphi)-ub_thlev(l,temi-1,iphi) ) /    &
     &         ( r_thlev_gwdc(l,temi+1)-r_thlev_gwdc(l,temi-1) )
        temp = temp*temp
        if (temp > 1.e-16) then
          ri_ct(l,iphi) = n_ct(l)**2 / temp
        else
          ri_ct(l,iphi) = 1.e10
        end if
      enddo
      enddo

      tem2d(:,:) = ub_cb(:,:) - ub_sfc(:,:)
      do iphi=1, nphi
      do l=1, ncol
        zcsa(l,iphi) = zcta(l)
        if ( abs(tem2d(l,iphi)) > 0.1 ) then
          tem1d(l) = (ub_ct(l,iphi)-ub_sfc(l,iphi))/tem2d(l,iphi)
          if (tem1d(l) > 0.0) then
            zcsa(l,iphi) = tem1d(l)*zcba(l)
            zcsa(l,iphi) = min( zcta(l), max( zcsa(l,iphi), zcba(l) ) )
          end if
        end if
      enddo
      enddo

! will not be used
!     ! revise u_cbottom after determining zcsa
!     do iphi=1, nphi
!       tem2d(:,iphi) = zcba(:)
!     enddo
!     ub_cb(:,:) = ub_sfc(:,:) + (ub_ct(:,:)-ub_sfc(:,:))*tem2d(:,:)/   &
!    &             zcsa(:,:)

!
! Calculate theta function that represents the horizontal and temporal 
! structure of diabatic forcing in phase speed and direction domain.
!
      q0sqc(:) = (cfactor*cp*heatmax(:))**2

      temp = ( (hscale*tscale) / (32.*pi**1.5) )**2
      q0sqc(:) = q0sqc(:) * temp

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

        c_int(:) = c_phase(:) - ub_ct(l,iphi)
        ic_cl(l,iphi) = minloc(abs(c_int),1)-(nc+1)
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
                ut_mc     = (-1.0)*c_int(ic)
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
                ut_mc = (-1.0)*c_int(ic)
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
            temp = x2tmp(ic-1) + x2tmp(ic+1)
            if (x_sq(ic) > temp) then
              x_sq(ic) = sqrt(temp*x_sq(ic))
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
          c_int(ic_cl(l,iphi)) = -999.   ! prevent from 0/0 in case
          wfrf(:,l,iphi) = x_sq(-nc:nc)/c_int(:)*n_ct(l)

        else  ! ri_cs <= 0.25

          wfrf(:,l,iphi) = 0.0

        end if  ! ri_cs > 0.25

      enddo  ! l
      enddo  ! iphi

!
! Calculate cloud-top momentum flux spectrum.
!
      temp = (2.0*(2.0*pi)**3)/(ah*lt)*(g/cp)**2
      do iphi=1, nphi
      do l=1, ncol
        if (ri_ct(l,iphi) > 0.25) then
          factor = rho_ct(l)*temp/(t_ct(l)*n_q(l)**2)**2
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
        diag_znwcq(l,9 ) = ut_gwdc(l,kcta_thlev(l))
        diag_znwcq(l,10) = vt_gwdc(l,kcta_thlev(l))
      enddo

      diag_znwcq(:,11:levels) = 1.e20


      RETURN

      END SUBROUTINE gw_ctop
