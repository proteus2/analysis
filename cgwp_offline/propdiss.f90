SUBROUTINE propdiss(  &
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

  USE subr_common

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
     &, l_mflx_v_ctop_on

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
     &, diag_spec      (nx, ny, nz, (nc*2+1)*3)

!      Real, Parameter :: g = 9.80665
      Real, Parameter :: two_omega = 2.*7.292116E-5

      Real, Parameter :: pi = 3.14159265358979323846
      Real, Parameter :: deg2rad = pi/180.0

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
     &, f2            (ncol)                                            &
     &, ub_tlev      (ncol, 0:nz, nphi)                            &
     &, c_phase       (-nc:nc)                                          &
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
     &, mfct_neg

      integer ::                                                        &
     &  kcta(ncol)                                                      &
     &, ic_cl(ncol, nphi)

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
     &  r_tlev_col(ncol, nz)

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

! lat_col
! ut_col, vt_col
! r_tlev_col
! rho_col
! rho_tlev_col
! nbv_tlev
! dc
! c_phase
! sec_theta_latitude_col
! kcta

!------------------------------------------------------------------
!  1-3. PREPARE SOME VARIABLES USED IN GWDC CALCULATIONS.
!------------------------------------------------------------------

  cosphi(:) = cos(phi_dir(:)*deg2rad)
  sinphi(:) = sin(phi_dir(:)*deg2rad)

  f2(:) = (two_omega*two_omega)*sin(lat_col(:)*deg2rad)**2

  call basic_u_phi(ut_col,vt_col,phi_dir, ub_tlev)

  do iphi=1, nphi
  do l=1, ncol
    ub_ct(l,iphi) = ub_tlev(l,kcta(l),iphi)
  enddo
  enddo
  call ind_c_critlayer(c_phase,-nc,ub_ct, ic_cl)

!------------------------------------------------------------------
!  2. CALCULATE THE CLOUD-TOP MOMENTUM FLUX SPECTRUM
!------------------------------------------------------------------

      call gw_ctop(nx,ny,ny,nz,nz_wet,             &
     &   ncol,nphi,nc,hscale,tscale,cfactor,lt,ah,igwdc,jgwdc,          &
     &   eta_tlev,                                              &
     &   ut_col ,vt_col ,rho_tlev_col ,r_tlev_col ,nbv_tlev,         &
     &   heat_col ,heatmax,ub_tlev,c_phase,c_m05dc,kcb,kct,            &
     &   t_tlev,cosphi,sinphi,dc,                                      &
     &   kcta,mfs_ct,ic_cl,diag_znwcq_col                              &
     &   )

!------------------------------------------------------------------
!  3. OPTAIN MOMENTUM FLUX PROFILES
!------------------------------------------------------------------
!  kcta,mfs_ct

  l_mflx_on = l_mflx_u_on .or. l_mflx_v_on .or.      &
              l_mflx_u_ctop_on .or. l_mflx_v_ctop_on

  l_mflx_on2 = l_mflx_u_on .or. l_mflx_v_on

  mf_east (:,:) = 0.0
  mf_north(:,:) = 0.0

  ! for calculating saturated spectrum
  tmp  = beta_wm/(sqrt(2.0)*pi)
  pm1  = p_wm - 1.0
  c2mp = 2.0 - p_wm
  do k=1, nz
  do l=1, ncol
    ome_min = sqrt(max(f2(l), nbv_tlev(l,k)*beta_eq/mstar))
    b0 = pm1*ome_min**pm1/(1.0-(ome_min/nbv_tlev(l,k))**pm1)
!    if (p_wm == 1.0)  b0 = 1.0/log(nbv_tlev(l,k)/ome_min)
    w0 = (nbv_tlev(l,k)**c2mp-ome_min**c2mp)/c2mp
    fact_s(l,k) = abs(tmp*b0*w0*rho_tlev_col(l,k))
  enddo
  enddo

      do iphi=1, nphi
      do l=1, ncol

        i = igwdc(l)
        j = jgwdc(l)

        do ic=-nc, nc
          mfsp(ic) = mfs_ct(ic,l,iphi)*dc
        enddo
        c_int(:) = c_phase(:) - ub_tlev(l,kcta(l),iphi)

        ineg = ic_cl(l,iphi) - 1
        ipos = ic_cl(l,iphi) + 1

        ! diagnostics - unfiltered cloud-top momentum flux
        ! saved at mflx_xxxx(:,:,1)
        if ( l_mflx_on2 ) then

          mfct_neg = 0.0
          mfct_pos = 0.0

          do ic=-nc, ineg    ! for components having negative MF
            mfct_neg = mfct_neg + mfsp(ic)
          enddo
          do ic=ipos, nc     ! for components having positive MF
            mfct_pos = mfct_pos + mfsp(ic)
          enddo

          if ( l_mflx_u_on ) then
            if (cosphi(iphi) > 0.0) then
              mflx_east(i,j,1) = mflx_east(i,j,1) +             &
     &                               mfct_pos*cosphi(iphi)
            else
              mflx_east(i,j,1) = mflx_east(i,j,1) +             &
     &                               mfct_neg*cosphi(iphi)
            end if
            if (cosphi(iphi) > 0.0) then
              mflx_west(i,j,1) = mflx_west(i,j,1) +             &
     &                               mfct_neg*cosphi(iphi)
            else
              mflx_west(i,j,1) = mflx_west(i,j,1) +             &
     &                               mfct_pos*cosphi(iphi)
            end if
          end if
          if ( l_mflx_v_on ) then
            mflx_north(i,j,1) =             &
     &                 mflx_north(i,j,1) + mfct_pos*sinphi(iphi)
            mflx_south(i,j,1) =             &
     &                 mflx_south(i,j,1) + mfct_neg*sinphi(iphi)
          end if

        end if  ! l_mflx_on2

!yhspec+ :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
        ! diagnostics - unfiltered cloud-top momentum flux spectrum
        ! saved at diag_spec(:,:,1)
        ! valid only for phi = [0,90] or [45,135]
        if ( l_spec_on .and. cosphi(iphi) /= 0.0 ) then
          if (cosphi(iphi) > 0.0) then
            do ic=ipos+(nc+1), nc+(nc+1)
              diag_spec(i,j,1,ic) =                                     &
     &           diag_spec(i,j,1,ic) + mfsp(ic-nc-1)
            enddo
            do ic=-nc+(nc*3+2), ineg+(nc*3+2)
              diag_spec(i,j,1,ic) =                                     &
     &           diag_spec(i,j,1,ic) + mfsp(ic-nc*3-2)
            enddo
          else
            do ic=-ineg+(nc+1), nc+(nc+1)
              diag_spec(i,j,1,ic) =                                     &
     &           diag_spec(i,j,1,ic) - mfsp(nc+1-ic)
            enddo
            do ic=-nc+(nc*3+2), -ipos+(nc*3+2)
              diag_spec(i,j,1,ic) =                                     &
     &           diag_spec(i,j,1,ic) - mfsp(nc*3+2-ic)
            enddo
          end if
        end if
!yhspec- :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! momentum flux profiles
        mf_neg(:) = 0.0
        mf_pos(:) = 0.0

        do k=kcta(l), nz

          ! calculate the saturation spectrum using Warner and
          ! McIntyre's method but as a function of phase speed,
          ! and apply the saturation condition

          do ic=-nc, ineg    ! for components having negative MF
            mfsp_s(ic) = fact_s(l,k)*c_int(ic)/nbv_tlev(l,k)*dc
            mfsp(ic) = max(mfsp_s(ic), mfsp(ic))
            mf_neg(k) = mf_neg(k) + mfsp(ic)
          enddo

          do ic=ipos, nc     ! for components having positive MF
            mfsp_s(ic) = fact_s(l,k)*c_int(ic)/nbv_tlev(l,k)*dc
            mfsp(ic) = min(mfsp_s(ic), mfsp(ic))
            mf_pos(k) = mf_pos(k) + mfsp(ic)
          enddo

          ! calculate eastward and northward MF
          tmp = mf_pos(k) + mf_neg(k)
          mf_east (l,k) = mf_east (l,k) + tmp*cosphi(iphi)
          mf_north(l,k) = mf_north(l,k) + tmp*sinphi(iphi)

!yhspec+ :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
        ! valid only for phi = [0,90] or [45,135]
        if ( l_spec_on .and. cosphi(iphi) /= 0.0 ) then

          if (k == kcta(l)) then
            if (cosphi(iphi) > 0.0) then
              do ic=ipos+(nc+1), nc+(nc+1)
                diag_spec(i,j,2,ic) =                                   &
     &             diag_spec(i,j,2,ic) + mfsp(ic-nc-1)
              enddo
              do ic=-nc+(nc*3+2), ineg+(nc*3+2)
                diag_spec(i,j,2,ic) =                                   &
     &             diag_spec(i,j,2,ic) + mfsp(ic-nc*3-2)
              enddo
            else
              do ic=-ineg+(nc+1), nc+(nc+1)
                diag_spec(i,j,2,ic) =                                   &
     &             diag_spec(i,j,2,ic) - mfsp(nc+1-ic)
              enddo
              do ic=-nc+(nc*3+2), -ipos+(nc*3+2)
                diag_spec(i,j,2,ic) =                                   &
     &             diag_spec(i,j,2,ic) - mfsp(nc*3+2-ic)
              enddo
            end if
          end if  ! k == kcta(l)

          if (cosphi(iphi) > 0.0) then
            do ic=ipos+(nc+1), nc+(nc+1)
              diag_spec(i,j,k,ic) =                                     &
     &           diag_spec(i,j,k,ic) + mfsp(ic-nc-1)
            enddo
            do ic=-nc+(nc*3+2), ineg+(nc*3+2)
              diag_spec(i,j,k,ic) =                                     &
     &           diag_spec(i,j,k,ic) + mfsp(ic-nc*3-2)
            enddo
          else
            do ic=-ineg+(nc+1), nc+(nc+1)
              diag_spec(i,j,k,ic) =                                     &
     &           diag_spec(i,j,k,ic) - mfsp(nc+1-ic)
            enddo
            do ic=-nc+(nc*3+2), -ipos+(nc*3+2)
              diag_spec(i,j,k,ic) =                                     &
     &           diag_spec(i,j,k,ic) - mfsp(nc*3+2-ic)
            enddo
          end if

        end if  ! l_spec_on .and. cosphi(iphi) /= 0.0
!yhspec- :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

          ! check the critical levels between k and k+1 levels,
          ! and prepare the next-level calculations
          if (k /= nz) then

            c_int(:) = c_phase(:) - ub_tlev(l,k+1,iphi)

            temi = ineg
            do ic=temi, -nc, -1
              if ( c_int(ic) >= 0.0 ) then
                mfsp(ic) = 0.0
                ineg = ic - 1
              else
                EXIT
              end if
            enddo

            temi = ipos
            do ic=temi, nc
              if ( c_int(ic) <= 0.0 ) then
                mfsp(ic) = 0.0
                ipos = ic + 1
              else
                EXIT
              end if
            enddo

          end if  ! k /= nz

        enddo  ! k

        ! diagnostics
        if ( l_mflx_on ) then

          if ( l_mflx_u_on ) then
            if (cosphi(iphi) > 0.0) then
              mflx_east(i,j,:) = mflx_east(i,j,:) +             &
     &                               mf_pos(:)*cosphi(iphi)
            else
              mflx_east(i,j,:) = mflx_east(i,j,:) +             &
     &                               mf_neg(:)*cosphi(iphi)
            end if
            if (cosphi(iphi) > 0.0) then
              mflx_west(i,j,:) = mflx_west(i,j,:) +             &
     &                               mf_neg(:)*cosphi(iphi)
            else
              mflx_west(i,j,:) = mflx_west(i,j,:) +             &
     &                               mf_pos(:)*cosphi(iphi)
            end if
          end if
          if ( l_mflx_v_on ) then
            mflx_north(i,j,:) =             &
     &                 mflx_north(i,j,:) + mf_pos(:)*sinphi(iphi)
            mflx_south(i,j,:) =             &
     &                 mflx_south(i,j,:) + mf_neg(:)*sinphi(iphi)
          end if

          k = kcta(l)
          if ( l_mflx_u_ctop_on ) then
            if (cosphi(iphi) > 0.0) then
              mflx_e_ctop(i,j) = mflx_e_ctop(i,j) +             &
     &                               mf_pos(k)*cosphi(iphi)
            else
              mflx_e_ctop(i,j) = mflx_e_ctop(i,j) +             &
     &                               mf_neg(k)*cosphi(iphi)
            end if
            if (cosphi(iphi) > 0.0) then
              mflx_w_ctop(i,j) = mflx_w_ctop(i,j) +             &
     &                               mf_neg(k)*cosphi(iphi)
            else
              mflx_w_ctop(i,j) = mflx_w_ctop(i,j) +             &
     &                               mf_pos(k)*cosphi(iphi)
            end if
          end if
          if ( l_mflx_v_ctop_on ) then
            mflx_n_ctop(i,j) =             &
     &                 mflx_n_ctop(i,j) + mf_pos(k)*sinphi(iphi)
            mflx_s_ctop(i,j) =             &
     &                 mflx_s_ctop(i,j) + mf_neg(k)*sinphi(iphi)
          end if

        end if  ! l_mflx_on

      enddo  ! l
      enddo  ! iphi

!yhspec+ :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  if ( l_spec_on ) then
    ! valid only for phi = [0,90] or [45,135]
    tmp = abs(cosphi(1))
    diag_spec(:,:,:,:) = diag_spec(:,:,:,:)*tmp
  end if
!yhspec- :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

!------------------------------------------------------------------
!  4. CALCULATE CONVECTIVE GRAVITY WAVE DRAG
!------------------------------------------------------------------

  if ( l_drag_u_on .or. l_drag_v_on ) then

    drag_u_grd(:,:,:) = 0.0
    drag_v_grd(:,:,:) = 0.0

    do l=1, ncol
      i = igwdc(l)
      j = jgwdc(l)
      do k=kcta(l)+1, nz
        tmp = 1.0/rho_col(l,k)/(r_tlev_col(l,k)-r_tlev_col(l,k-1))
        drag_u_grd(i,j,k) = -(mf_east (l,k)-mf_east (l,k-1))*tmp
        drag_v_grd(i,j,k) = -(mf_north(l,k)-mf_north(l,k-1))*tmp
      enddo
    enddo

    do l=1, ncol
      i = igwdc(l)
      j = jgwdc(l)
      k = kcta(l)
      tmp = 1.0/rho_col(l,k)/(r_tlev_col(l,k)-r_tlev_col(l,k-1))
      drag_u_grd(i,j,k) = -(mf_east (l,k)-mf_east (l,nz))*tmp
      drag_v_grd(i,j,k) = -(mf_north(l,k)-mf_north(l,nz))*tmp
    enddo
    ! as the outward flux of momentum at the model top is allowed,
    ! drags below the cloud top must be reduced for the momentum
    ! budget of the model (i.e., mf(nz)*tmp)

  end if

  RETURN

END SUBROUTINE propdiss

