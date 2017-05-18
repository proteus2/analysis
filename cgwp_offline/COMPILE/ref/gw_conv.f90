SUBROUTINE gw_conv(  &
    nx, ny, nz, nz_wet, ncol,                  &
    nphi, nc, c_max, phi_dir,               &
    hscale, tscale, cfactor, lt, ah,         &
    mstar, p_wm, s_wm, t_wm, beta_wm,                                   &
        u,v,theta,exner_tlev,rho,rho_tlev,             &
        eta_tlev,r_tlev,sec_theta_latitude,             &
        scheat,schmax,kscbas,ksctop,igwdc,jgwdc,                        &
! diagnostics
        drag_u_grd, drag_v_grd, l_drag_u_on,l_drag_v_on,          &
        mflx_east,mflx_west,                                    &
        mflx_north,mflx_south,                                  &
        l_mflx_u_on,l_mflx_v_on,                   &
        mflx_e_ctop,mflx_w_ctop,                                &
        mflx_n_ctop,mflx_s_ctop,                                &
        l_mflx_u_ctop_on,l_mflx_v_ctop_on,                          &
        diag_znwcq,l_znwcq_on,                                       &
        diag_spec,l_spec_on                                          &
     &  )

!
! PURPOSE:  To calculate the convective gravity wave drag
!
! METHOD:
!
!   1) Interpolate data and gather at the GWDC columns into 2-D array.
!   2) Call <gw_ctop> to calculate the cloud-top momentum flux spectrum
!   3) Obtain momentum flux profiles (Warner and McIntyre, 1999, EPS; 
!      Song and Chun, 2006, JKMS)
!   4) Calculate convective gravity wave drag
!
! HISTORY:
!

      implicit none

!
! SUBROUTINE ARGUMENTS
!
!   model dimensions
 
      integer, intent(in) ::                                            &
     &  nx, ny, nz, nz_wet                                              &
     &, ncol               ! number of columns GWDC is calculated at.

!   model switches

      logical, intent(in) ::                                            &
     &  l_drag_u_on                                                  &
     &, l_drag_v_on                                                  &
     &, l_mflx_u_on                                                &
     &, l_mflx_v_on                                               &
     &, l_mflx_u_ctop_on                                              &
     &, l_mflx_v_ctop_on                                              &
     &, l_znwcq_on

!   GWDC parameters

      integer, intent(in) ::                                            &
     &  nphi                                                            &
     &, nc

      real, intent(in) ::                                               &
     &  c_max                                                           &
     &, hscale                                                          &
     &, tscale                                                          &
     &, cfactor                                                         &
     &, lt                                                              &
     &, ah                                                              &
     &, mstar                                                           &
     &, p_wm                                                            &
     &, s_wm                                                            &
     &, t_wm                                                            &
     &, beta_wm                                                         &
     &, phi_dir(nphi)

!   data arrays

      real, intent(in) ::                                               &
     &  u(nx, ny, nz)         &
     &, v(nx, ny, nz)        &
                           ! primary u, v field
     &, theta(nx, ny, nz)     &
                           ! primary theta field
     &, exner_tlev(nx, ny, nz)                  &
                           ! exner ftn field
     &, rho(nx, ny, nz)                               &
                           ! rho field
     &, rho_tlev(nx, ny, nz-1)                       &
                           ! rho on theta levels
     &, scheat(nx, ny, nz_wet)                            &
     &, schmax(nx, ny)

      integer, intent(in) ::                                            &
     &  igwdc(ncol)                                                     &
     &, jgwdc(ncol)                                                     &
     &, kscbas(nx, ny)                                        &
     &, ksctop(nx, ny)

!   coordinate arrays

      real, intent(in) ::                                               &
     &  eta_tlev(0:nz)                                      &
     &, r_tlev(nx,ny, 0:nz)                  &
                           ! r on theta levels
     &, sec_theta_latitude(nx,ny)

!   output variables

      real, intent(out) ::                                            &
     &  drag_u_grd(nx, ny, nz)                            &
     &, drag_v_grd(nx, ny, nz)

!   diagnostic variables

      real, intent(inout) ::                                            &
     &  mflx_east  (nx, ny, nz)                       &
     &, mflx_west  (nx, ny, nz)                       &
     &, mflx_north (nx, ny, nz)                       &
     &, mflx_south (nx, ny, nz)                       &
     &, mflx_e_ctop(nx, ny)                               &
     &, mflx_w_ctop(nx, ny)                               &
     &, mflx_n_ctop(nx, ny)                               &
     &, mflx_s_ctop(nx, ny)                               &
     &, diag_znwcq     (nx, ny, nz)                       &
     &, diag_spec      (nx, ny, nz, (nc*2+1)*3)

      Real, Parameter :: G = 9.80665
      Real, Parameter :: Pi = 3.14159265358979323846
      Real, Parameter :: Pi_Over_180 = Pi/180.0
      Real, Parameter :: two_omega = 2.*7.292116E-5

!
! LOCAL VARIABLES
!
!   data arrays

      real ::                                                           &
     &  ut_col        (ncol, 0:nz)                                  &
     &, vt_col        (ncol, 0:nz)                                  &
     &, rho_col       (ncol, nz)                                    &
     &, rho_tlev_col (ncol, nz)                                    &
     &, nbv_tlev     (ncol, nz)                                    &
     &, pt_col        (ncol, nz)                                    &
     &, t_tlev       (ncol, nz_wet)                                &
     &, heat_col      (ncol, nz_wet)                                &
     &, heatmax       (ncol)                                            &
     &, ub_tlev      (ncol, 0:nz, nphi)                            &
     &, c_phase       (-nc:nc)                                          &
     &, c_m05dc       (-nc:nc+1)                                        &
     &, c_int         (-nc:nc)

      real ::                                                           &
     &  mfs_ct(-nc:nc, ncol, nphi)                                      &
     &, mf_east (ncol, nz)                                          &
     &, mf_north(ncol, nz)                                          &
     &, drag_u  (ncol, nz)                                          &
     &, drag_v  (ncol, nz)                                          &
     &, fact_s  (ncol, nz)                                          &
     &, mf_pos(nz)                                                  &
     &, mf_neg(nz)                                                  &
     &, mfsp  (-nc:nc)                                                  &
     &, mfsp_s(-nc:nc)                                                  &
     &, mfct_pos                                                        &
     &, mfct_neg                                                        &
     &, diag_znwcq_col(ncol, nz)

      integer ::                                                        &
     &  kcb (ncol)                                                      &
     &, kct (ncol)                                                      &
     &, kcta(ncol)

      real ::                                                           &
     &  dc                                                              &
     &, cosphi(nphi)                                                    &
     &, sinphi(nphi)

      real ::                                                           &
     &  b0                                                              &
     &, w0                                                              &
     &, ome_min

      logical ::                                                        &
     &  l_mflx_on                                                     &
     &, l_mflx_on2                                                    &
     &, l_spec_on

!   coordinate arrays

      real ::                                                           &
     &  r_tlev_col (ncol, nz)

!   work arrays

      real ::  tmp, pm1, c2mp

      integer ::  temi, ipos, ineg

      integer ::  i,j,k,l,ic,iphi   ! loop counters

!   parameters and constants

      real, parameter ::                                                &
     &  n2bv_min = 1.e-6                                                &
     &, beta_eq  = 2.3e-11
!
! EXTERNAL ROUTINES
!
      external ::  gw_ctop

  if (ncol < 1)  RETURN

!------------------------------------------------------------------
!  1-1. INTERPOLATE WINDS TO THETA-GRID AND GATHER AT GWDC ARRAY.
!------------------------------------------------------------------

      ! Approximate theta-level winds.
      do k=1, nz-1
      do l=1, ncol
        i = igwdc(l)
        j = jgwdc(l)
        ut_col(l,k) = 0.5*(u(i,j,k) + u(i,j,k+1))
        vt_col(l,k) = 0.5*(v(i,j,k) + v(i,j,k+1))
      enddo
      enddo

      do l=1, ncol
        i = igwdc(l)
        j = jgwdc(l)
        ut_col(l,0 ) = u(i,j,1 )
        vt_col(l,0 ) = v(i,j,1 )
        ut_col(l,nz) = u(i,j,nz)
        vt_col(l,nz) = v(i,j,nz)
      enddo

!------------------------------------------------------------------
!  1-2. GATHER RHO, N, HEIGHT AND DCH AT GWDC ARRAY.
!------------------------------------------------------------------

      do k=1, nz
      do l=1, ncol
        i = igwdc(l)
        j = jgwdc(l)
        r_tlev_col(l,k) = r_tlev(i,j,k) - r_tlev(i,j,0)
        rho_col(l,k) = rho(i,j,k)
        pt_col (l,k) = theta  (i,j,k)
      enddo
      enddo

      do k=1, nz-1
      do l=1, ncol
        rho_tlev_col(l,k) = rho_tlev(igwdc(l),jgwdc(l),k)
      enddo
      enddo

      ! extrapolate rho at top (with constant scale height)
      rho_tlev_col(:,nz) = rho_tlev_col(:,nz-1)*            &
     &      ( rho_tlev_col(:,nz-1)/rho_tlev_col(:,nz-2) )** &
     &      ( (r_tlev_col(:,nz  )-r_tlev_col(:,nz-1))/      &
     &        (r_tlev_col(:,nz-1)-r_tlev_col(:,nz-2)) )

      do k=2, nz-1
      do l=1, ncol
        nbv_tlev(l,k) = g/pt_col(l,k)*(pt_col(l,k+1)-pt_col(l,k-1)) &
     &                    / (r_tlev_col(l,k+1)-r_tlev_col(l,k-1))
        nbv_tlev(l,k) = max(n2bv_min, nbv_tlev(l,k))
        nbv_tlev(l,k) = sqrt(nbv_tlev(l,k))
      enddo
      enddo
      nbv_tlev(:,1 ) = nbv_tlev(:,2   )
      nbv_tlev(:,nz) = nbv_tlev(:,nz-1)

      do k=1, nz_wet
      do l=1, ncol
        i = igwdc(l)
        j = jgwdc(l)
        t_tlev (l,k) = theta(i,j,k)*exner_tlev(i,j,k)
        heat_col(l,k) = scheat(i,j,k)
      enddo
      enddo

      do l=1, ncol
        i = igwdc(l)
        j = jgwdc(l)
        heatmax(l) = schmax(i,j)
        kcb    (l) = kscbas(i,j)
        kct    (l) = ksctop(i,j)
      enddo

!------------------------------------------------------------------
!  1-3. PREPARE SOME VARIABLES USED IN GWDC CALCULATIONS.
!------------------------------------------------------------------

      do iphi=1, nphi
      do k=0, nz
      do l=1, ncol
        ub_tlev(l,k,iphi) = ut_col(l,k)*cosphi(iphi) +                &
     &                       vt_col(l,k)*sinphi(iphi)
      enddo
      enddo
      enddo

!------------------------------------------------------------------
!  2. CALCULATE THE CLOUD-TOP MOMENTUM FLUX SPECTRUM
!------------------------------------------------------------------

      call gw_ctop(nx,ny,ny,nz,nz_wet,ncol,         &
     &   eta_tlev,                                              &
     &   ut_col ,vt_col ,rho_tlev_col ,r_tlev_col ,nbv_tlev,         &
     &   heat_col ,heatmax,ub_tlev,kcb,kct,             &
     &   t_tlev,                                                       &
     &   kcta,mfs_ct,diag_znwcq_col                              &
     &   )

  if ( l_znwcq_on ) then
    diag_znwcq(:,:,:) = 1.e20
    do k=1, 10
    do l=1, ncol
      diag_znwcq(igwdc(l),jgwdc(l),k) = diag_znwcq_col(l,k)
    enddo
    enddo
  end if

  RETURN

END SUBROUTINE gw_conv

