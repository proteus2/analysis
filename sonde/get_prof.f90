PROGRAM get_prof

!  NAMELISTS
!-----------------------------------------------------------------------
! fi_dir   : base directory in which the input file exists
! fi_head  : head of the input file
! fi_tail  : tail of the input file
! nline_h  : number of header lines
! MIS_in   : missing value assigned in the input data
! fo_dir   : base directory in which the output file will be saved
!-----------------------------------------------------------------------
! limit_in : length used to filter the original data in data-reading process
!            when the height interval is too short (i.e., shorter than limit_in)
!            The filtered data is made so that it has intervals (at least)
!            longer than limit_in.
! itv      : itv*limin_in is used to cut-out the data in data-reading process.
!            It is needed when the data is bad (i.e., One of the intervals is 
!            longer than itv*limin_in) so that the interpolation may produce 
!            large error).
!            Too small itv may cause a short vertical profile.
!            Too large itv may cause large errors in the interpolation process.
!-----------------------------------------------------------------------
! dz       : interval for interpolating
!-----------------------------------------------------------------------
! zab      : bottom height for analysis. It should be a multiple of dz
! zat      : top height for analysis. It should be a multiple of dz
! order_p  : order of polynomials used to obtain basic states
! zbuf     : buffer height for polynomial fitting to obtain basic states.
!            It should be a multiple of dz
! p0_hydro : option for the basic-state pressure to be hydrostatically balanced.
!            If .false., it is defined by the polynomial regression of the
!            pressure profile, like the other variables.
!-----------------------------------------------------------------------

  use m_etc
  use easync
  
  implicit none

! namelist parameters
  integer            ::  nline_h, itv, order_p
  logical            ::  p0_hydro
  real               ::  MIS_in, limit_in, dz, zab, zat, zbuf
  character(len=256) ::  fi_dir, fo_dir
  character(len=64)  ::  fi_head, fi_tail

  namelist /IO_PARAM/  fi_dir, fi_head, fi_tail, nline_h, MIS_in, fo_dir
  namelist /FILTERING_PARAM/  limit_in, itv
  namelist /INTERPOLATION_PARAM/  dz
  namelist /BASE_PERT/  zab, zat, order_p, zbuf, p0_hydro

! main
  integer                         ::  nlev, nzt, nzw, nzp, nzf, nza,  &
                                      nbufb, nbuft
  integer                         ::  kzab, kzat, kzb2, kzt2
  real, dimension(:), allocatable ::  h, t, ws, deg, p
  real, dimension(:), allocatable ::  h_t, h_w, h_p, t1, u1, v1, p1,  &
                                      ws1, deg1
  real, dimension(:), allocatable ::  zf, tf, uf, vf, pf
  real, dimension(:), allocatable ::  za, ta, ua, va, pa,              &
                                      ta0, ua0, va0, pa0, ra0, nbva0,  &
                                      ta_prt, ua_prt, va_prt, pa0h
  real                            ::  ta0_bb, ta0_bt, pa0_bb, pa0_bt,  &
                                      zb2, zt2, lnpa0_avg, tmp
  character(len=256)              ::  ifn, ofn1, ofn2, ofn3
  character(len=256)              ::  tit1, tit2, tit3

  integer                         ::  i,k, tmpi

! constants
  real, parameter ::  rd = 287.0, cp = 1004.0, g = 9.806
  real, parameter ::  MIS = 1.e30
  real, parameter ::  deg2rad = 3.14159265358979323846/180.

!-----------------------------------------------------------------------
!  GET ARGUMENTS AND READ NAMELISTS
!-----------------------------------------------------------------------

  call input_arg
  call input_namelist

!-----------------------------------------------------------------------
!  INPUT DATA FILENAME
!-----------------------------------------------------------------------

  write(ifn,'(a,6a)') trim(fi_dir)//'/'//trim(stid)//'/',  &
       trim(fi_head),trim(year),month,date,hour,trim(fi_tail)
  write(6,'(/,a)') ' input : '//trim(ifn)

!-----------------------------------------------------------------------
!  OUTPUT FILENAME
!-----------------------------------------------------------------------

  write(tit1,'(a)') 'Temperature and wind profiles filtered from raw data'

  write(ofn1,'(a,6a)') trim(fo_dir)//'/'//trim(stid)//'/',  &
      'prof_filt_',trim(year),month,date,hour,'.nc'

  write(tit2,'(a)') 'Interpolated variables'

  write(ofn2,'(a,6a)') trim(fo_dir)//'/'//trim(stid)//'/',  &
      'prof_intp_',trim(year),month,date,hour,'.nc'

  write(tit3,'(a)') 'Basic-state and perturbation variables'

  write(ofn3,'(a,6a)') trim(fo_dir)//'/'//trim(stid)//'/',  &
      'prof_pert_',trim(year),month,date,hour,'.nc'

!-----------------------------------------------------------------------
!  READ AND FILTER INPUT DATA
!-----------------------------------------------------------------------

  call read_data
  ! out: h, p, t, and (ws, deg) with nlev.
  ! Data where h is missing are removed.
  ! missing value : MIS

  call filter_data
  ! Now, the data have intervals about limit_in.
  ! nlev is updated.

  ! remove missing data
  call remove_missv
  ! out: p1, t1, and (ws1, deg1) with h_p, h_t, and h_w, respectively

  deallocate(h, p, t, ws, deg)

!-----------------------------------------------------------------------
!  DEFINE VARIABLES
!-----------------------------------------------------------------------

  t1(:) = t1(:) + 273.15  ! [K]

  allocate(u1(nzw), v1(nzw))
! check the definition of the wind direction
  u1(:) = -ws1(:)*sin(deg1(:)*deg2rad)
  v1(:) = -ws1(:)*cos(deg1(:)*deg2rad)

  deallocate(ws1, deg1)

!-----------------------------------------------------------------------
!  DUMP THE FILTERED VARIABLES
!-----------------------------------------------------------------------

  call dump_1

!-----------------------------------------------------------------------
!  INTERPOLATION 
!-----------------------------------------------------------------------

  call interpol_vars
  ! out: zf, tf, uf, vf, pf (1:nzf)
  ! The extrapolated data below the bottom of data is usually not valid,
  ! but we include it for convenience (i.e., All profiles start from
  ! zf = 0).

  deallocate(h_t, h_w, h_p, t1, u1, v1, p1)

!-----------------------------------------------------------------------
!  DUMP THE INTERPOLATED VARIABLES
!-----------------------------------------------------------------------
  
  call dump_2

!=======================================================================
!  HERE, we have zf, tf, uf, vf, and pf (1:nzf)
!=======================================================================

  if (zf(nzf) < zat)  call stop_message(  &
        'The top height of data is lower than the critical value.\n'//  &
        'Analysis not performed')

!-----------------------------------------------------------------------
!  SET DATA RANGE FOR ANALYSIS
!-----------------------------------------------------------------------

  ! top and bottom of the range including the buffer
  zb2 = max(zf(1)  , zab - zbuf)
  zt2 = min(zf(nzf), zat + zbuf)

  ! indices
  ! We assume that zab, zat, and zbuf are multiples of dz
  kzab = nint(zab/dz)+1
  kzat = nint(zat/dz)+1
  kzb2 = nint(zb2/dz)+1
  kzt2 = nint(zt2/dz)+1
  nza   = kzat - kzab  + 1
  nbufb = kzab - kzb2
  nbuft = kzt2 - kzat

  write(6,*) 'nz_analysis : ', nza

  allocate(za(nza), ta(nza), ua(nza), va(nza), pa(nza))
  za(:) = zf(kzab:kzat)
  ta(:) = tf(kzab:kzat)
  ua(:) = uf(kzab:kzat)
  va(:) = vf(kzab:kzat)
  pa(:) = pf(kzab:kzat)

!-----------------------------------------------------------------------
!  OBTAIN THE BASIC STATE AND PERTURBATION
!-----------------------------------------------------------------------

  call regress_uvtp
  ! out: ua0, va0, ta0, pa0 (1:nza)

  if ( p0_hydro ) then
    allocate( pa0h(nza) )
    call cal_p_hydro
    ! out: pa0h
    pa0(:) = pa0h(:)
    deallocate( pa0h )
  end if

  ! rho0
  allocate(ra0(nza))
  ra0(:) = 100.*pa0(:)/(rd*ta0(:))

  ! Brunt-Vaisala freq.
  ! use ta0 and (ta0_bb, ta0_bt), and
  !     pa0 and (pa0_bb, pa0_bt) if p0_hydro is .false.
  call cal_nbv
  ! out: nbva0

  ! perturbation
  allocate(ta_prt(nza), ua_prt(nza), va_prt(nza))
  ta_prt(:) = ta(:) - ta0(:)
  ua_prt(:) = ua(:) - ua0(:)
  va_prt(:) = va(:) - va0(:)

!-----------------------------------------------------------------------
!  DUMP THE BASIC-STATE AND PERTURBATION VARIABLES
!-----------------------------------------------------------------------

  call dump_3

  deallocate(za, ta, ua, va, pa)
  deallocate(ta0, ua0, va0, pa0, ra0, nbva0)
  deallocate(ta_prt, ua_prt, va_prt)

  STOP


  CONTAINS


SUBROUTINE input_namelist

  open(10, file=trim(f_namelist), status='old')
  read(10, IO_PARAM)  ;  read(10, FILTERING_PARAM)
  read(10, INTERPOLATION_PARAM)  ;  read(10, BASE_PERT)
  close(10)

END subroutine input_namelist

SUBROUTINE read_data
 
  integer            ::  istat
  real               ::  dummy1, dummy2
  character(len=256) ::  tmpch

  open(11,file=ifn,action='read',form='formatted')
  do i=1, nline_h
    read(11,*) tmpch   ! header line
  enddo
  nlev = -1
  istat = 0
  do while (istat == 0)
    nlev = nlev + 1
    read(11,*,iostat=istat) tmpch
  enddo
  close(11)

  allocate( h(nlev), p(nlev), t(nlev), ws(nlev), deg(nlev) )
  h(:) = MIS_in  ;  p(:) = MIS_in  ;  t(:) = MIS_in
  ws(:) = MIS_in  ;  deg(:) = MIS_in

  open(11,file=ifn,action='read',form='formatted')
  do i=1, nline_h
    read(11,*) tmpch   ! header line
  enddo
  do i=1, nlev
    read(11,*, iostat=istat)  &
        dummy1, h(i), p(i), t(i), dummy2, ws(i), deg(i)
    ! delete below, if you have perfect datasets
!    if ( p(i) == h(i) .and. p(i) == t(i) ) then  ! i.e., absence of data in the line
!      h  (i) = MIS_in
!      p  (i) = MIS_in
!      t  (i) = MIS_in
!      ws (i) = MIS_in
!      deg(i) = MIS_in
!    end if
  enddo
  close(11)
  if ( all(deg(:) == deg(1)) )  call stop_message(  &
     'PROGRAM STOP - Check the input file format.')
 
  ! remove data where the height is missing
  tmpi = 0
  do i=1, nlev
    if (h(i) /= MIS_in) then
      tmpi = tmpi + 1
      h  (tmpi) = h  (i)
      p  (tmpi) = p  (i)
      t  (tmpi) = t  (i)
      ws (tmpi) = ws (i)
      deg(tmpi) = deg(i)
    end if
  enddo
  if (nlev > tmpi) then
    h  (tmpi+1:) = MIS
    p  (tmpi+1:) = MIS
    t  (tmpi+1:) = MIS
    ws (tmpi+1:) = MIS
    deg(tmpi+1:) = MIS
  end if

  nlev = tmpi

  ! change missing values
  where (p(1:nlev) == MIS_in)  p(1:nlev) = MIS
  where (t(1:nlev) == MIS_in)  t(1:nlev) = MIS
  where ( ws(1:nlev) == MIS_in .or. deg(1:nlev) == MIS_in )
    ws (1:nlev) = MIS
    deg(1:nlev) = MIS
  end where
 
END subroutine read_data

SUBROUTINE filter_data

  tmpi = 2
  do i=2, nlev

    ! use the data sparsely, if it's interval is too short.
    ! Reversed data (i.e., h(i) > h(i+1)) are also filtered.
    if ( (h(i)-h(tmpi-1)) >= limit_in ) then
      tmpi = tmpi + 1
      tmpi = min(i,tmpi)
    end if
    h  (tmpi) = h  (i)
    p  (tmpi) = p  (i)
    t  (tmpi) = t  (i)
    ws (tmpi) = ws (i)
    deg(tmpi) = deg(i)

    ! cut-out the data, if it's interval is too long.
    if ( (h(tmpi)-h(tmpi-1)) > (itv*limit_in) ) then
      write(6,*) 'DATA INTERVAL IS TOO LONG BETWEEN Z =',  &
                 h(tmpi-1), 'and', h(tmpi), 'm'
      tmpi = tmpi - 1
      EXIT
    end if

  enddo
  if ( (h(tmpi)-h(tmpi-1)) < limit_in ) tmpi = tmpi - 1

  if (nlev > tmpi) then
    h  (tmpi+1:) = MIS
    p  (tmpi+1:) = MIS
    t  (tmpi+1:) = MIS
    ws (tmpi+1:) = MIS
    deg(tmpi+1:) = MIS
  end if

  nlev = tmpi

  if (nlev < 20)  call stop_message(  &
     'The number of valid data is smaller than 20')

END subroutine filter_data

SUBROUTINE remove_missv

  call get_nz_nonmissing(t , nzt)
  call get_nz_nonmissing(ws, nzw)
  call get_nz_nonmissing(p , nzp)

  if (nzt == 0)  call stop_message('NO AVAILABLE DATA FOR TEMPERATURE')
  if (nzw == 0)  call stop_message('NO AVAILABLE DATA FOR WIND')
  if (nzp == 0)  call stop_message('NO AVAILABLE DATA FOR PRESSURE')

  allocate(h_t(nzt), t1(nzt))
  allocate(h_w(nzw), ws1(nzw), deg1(nzw))
  allocate(h_p(nzp), p1(nzp))

  call remove_missing(t  (1:nlev),h(1:nlev),nzt, t1  ,h_t)
  call remove_missing(ws (1:nlev),h(1:nlev),nzw, ws1 ,h_w)
  call remove_missing(deg(1:nlev),h(1:nlev),nzw, deg1,h_w)
  call remove_missing(p  (1:nlev),h(1:nlev),nzp, p1  ,h_p)
 
END subroutine remove_missv

SUBROUTINE get_nz_nonmissing(var,nz_nm)

  real,    dimension(:), intent(in ) ::  var
  integer,               intent(out) ::  nz_nm

  integer ::  ii

  nz_nm = 0
  do ii=1, size(var)
    if (var(ii) /= MIS)  nz_nm = nz_nm + 1
  end do

END subroutine get_nz_nonmissing

SUBROUTINE remove_missing(var,h_var,nz_nm,v1,h1)

  integer,                   intent(in ) ::  nz_nm
  real,    dimension(:),     intent(in ) ::  var, h_var
  real,    dimension(nz_nm), intent(out) ::  v1, h1

  integer ::  ii, jj

  jj = 0
  do ii=1, size(var)
    if (var(ii) /= MIS) then
      jj = jj + 1
      v1(jj) = var(ii)
      h1(jj) = h_var(ii)
    end if
  end do

END subroutine remove_missing

SUBROUTINE interpol_vars

  nzf = int( min(min(h_t(nzt),h_w(nzw)),h_p(nzp))/dz )+1

  allocate(zf(nzf), tf(nzf), uf(nzf), vf(nzf), pf(nzf))

  do k=1, nzf
    zf(k) = dz*float(k-1)
  end do

  ! interporation
  call cub_spl_intp0(nzt,h_t,t1,nzf,zf, tf)
  call cub_spl_intp0(nzw,h_w,u1,nzf,zf, uf)
  call cub_spl_intp0(nzw,h_w,v1,nzf,zf, vf)
  call cub_spl_intp0(nzp,h_p,log(p1),nzf,zf, pf)
  pf(:) = exp(pf(:))

END subroutine interpol_vars

SUBROUTINE regress_uvtp

  ! Chi^2 regression with basis functions poly_func

  use nr,  only: lfit

  integer                           ::  ma
  logical                           ::  maska(20)
  real                              ::  chisq, a(20),covar(20,20)
  real, dimension(:),   allocatable ::  sig
  real, dimension(:),   allocatable ::  zn
  real, dimension(:,:), allocatable ::  fzn

  INTERFACE
    SUBROUTINE poly_func(x,arr)
    use nrtype
    implicit none
    real(sp),               intent(in ) :: x
    real(sp), dimension(:), intent(out) :: arr
    END subroutine poly_func
  END interface

  allocate(ua0(nza), va0(nza), ta0(nza), pa0(nza))

  ma = order_p + 1

  allocate(zn(1-nbufb:nza+nbuft), sig(1-nbufb:nza+nbuft))
  allocate(fzn(ma,1-nbufb:nza+nbuft))

  maska(1:ma) = .true.
  sig(:) = 1.

  ! use normalized height (0~1) for the polynomial fitting
  zn(:) = (zf(kzb2:kzt2) - zb2)/(zt2 - zb2)
  do k=1-nbufb, nza+nbuft
    call poly_func(zn(k), fzn(:,k))
  enddo

  ! u0
  a(1:ma) = 0.
  call lfit(zn,uf(kzb2:kzt2),sig,a(1:ma),maska(1:ma),covar(1:ma,1:ma),  &
            chisq,poly_func)
  ua0(:) = sum(spread(a(1:ma),2,nza)*fzn(:,1:nza), dim=1)

  ! v0
  a(1:ma) = 0.
  call lfit(zn,vf(kzb2:kzt2),sig,a(1:ma),maska(1:ma),covar(1:ma,1:ma),  &
            chisq,poly_func)
  va0(:) = sum(spread(a(1:ma),2,nza)*fzn(:,1:nza), dim=1)

  ! T0
  a(1:ma) = 0.
  call lfit(zn,tf(kzb2:kzt2),sig,a(1:ma),maska(1:ma),covar(1:ma,1:ma),  &
            chisq,poly_func)
  ta0(:) = sum(spread(a(1:ma),2,nza)*fzn(:,1:nza), dim=1)
  if (nbufb /= 0)  ta0_bb = sum(a(1:ma)*fzn(:,0    ))
  if (nbuft /= 0)  ta0_bt = sum(a(1:ma)*fzn(:,nza+1))

  ! p0
  a(1:ma) = 0.
  call lfit(zn,pf(kzb2:kzt2),sig,a(1:ma),maska(1:ma),covar(1:ma,1:ma),  &
            chisq,poly_func)
  pa0(:) = sum(spread(a(1:ma),2,nza)*fzn(:,1:nza), dim=1)
  if (nbufb /= 0)  pa0_bb = sum(a(1:ma)*fzn(:,0    ))
  if (nbuft /= 0)  pa0_bt = sum(a(1:ma)*fzn(:,nza+1))

  deallocate(zn, fzn, sig)

END subroutine regress_uvtp

SUBROUTINE cal_p_hydro

  pa0h(1) = 0.
  do k=2, nza
    pa0h(k) = pa0h(k-1) - dz*(g/(rd*0.5*(ta0(k)+ta0(k-1))))
  enddo
  tmp = sum(pa0h(:))/float(nza)
  lnpa0_avg = sum(log(pa0(:)))/float(nza)
  pa0h(:) = exp(pa0h(:) - tmp + lnpa0_avg)
  if ( abs(pa0(1) - pa0h(1))/pa0h(1) > 0.1 .or.  &
       abs(pa0(nza) - pa0h(nza))/pa0h(nza) > 0.1 ) then
    write(6,'(/,a)')  &
        ' Difference btw. REGRESSED and HYDROSTATIC pressure > 10%.'
    write(6,*) '  bottom : ', pa0(1), pa0h(1)
    write(6,*) '  top    : ', pa0(nza), pa0h(nza)
  end if
 
END subroutine cal_p_hydro

SUBROUTINE cal_nbv

  allocate(nbva0(nza))

  if ( p0_hydro ) then
    do k=2, nza-1
      nbva0(k) = g/ta0(k)*( (ta0(k+1)-ta0(k-1))/(2.0*dz) + g/cp )
    end do
    nbva0(1) = nbva0(2)  ;  nbva0(nza) = nbva0(nza-1)
    if (nbufb /= 0)  nbva0(1) = g/ta0(1)*( (ta0(2)-ta0_bb)/(2.0*dz) + g/cp )
    if (nbuft /= 0)  nbva0(nza) = g/ta0(nza)*  &
        ( (ta0_bt-ta0(nza-1))/(2.0*dz) + g/cp )
  else
    do k=2, nza-1
      nbva0(k) = g*( (ta0(k+1)-ta0(k-1))/ta0(k) -  &
                     rd/cp*log(pa0(k+1)/pa0(k-1)) )/(2.0*dz)
    end do
    nbva0(1) = nbva0(2)  ;  nbva0(nza) = nbva0(nza-1)
    if (nbufb /= 0)  nbva0(1) =  &
        g*( (ta0(2)-ta0_bb)/ta0(1) - rd/cp*log(pa0(2)/pa0_bb) )/(2.0*dz)
    if (nbuft /= 0)  nbva0(nza) =  &
        g*( (ta0_bt-ta0(nza-1))/ta0(nza) -  &
            rd/cp*log(pa0_bt/pa0(nza-1)) )/(2.0*dz)
  end if

  do k=1, nza
    nbva0(k) = max(0., nbva0(k))
  enddo
  nbva0(:) = sqrt(nbva0(:))

END subroutine cal_nbv

SUBROUTINE dump_1

  call put_var('overwrite',ofn1,'Z_T',h_t)
  call put_var('append'   ,ofn1,'Z_W',h_w)
  call put_var('append'   ,ofn1,'Z_P',h_p)
  call put_var('append'   ,ofn1,'T'  ,t1 ,axis='Z_T')
  call put_var('append'   ,ofn1,'U'  ,u1 ,axis='Z_W')
  call put_var('append'   ,ofn1,'V'  ,v1 ,axis='Z_W')
  call put_var('append'   ,ofn1,'P'  ,p1 ,axis='Z_P')
  call put_att(ofn1,'global','title',tit1)

END subroutine dump_1

SUBROUTINE dump_2

  call put_var('overwrite',ofn2,'Z',zf)
  call put_var('append'   ,ofn2,'T',tf,axis='Z')
  call put_var('append'   ,ofn2,'U',uf,axis='Z')
  call put_var('append'   ,ofn2,'V',vf,axis='Z')
  call put_var('append'   ,ofn2,'P',pf,axis='Z')
  call put_att(ofn2,'global','title',tit2)

END subroutine dump_2

SUBROUTINE dump_3

  call put_var('overwrite',ofn3,'Z',za)
  call put_var('append'   ,ofn3,'T',ta,axis='Z')
  call put_var('append'   ,ofn3,'U',ua,axis='Z')
  call put_var('append'   ,ofn3,'V',va,axis='Z')
  call put_var('append'   ,ofn3,'P',pa,axis='Z')
  call put_var('append'   ,ofn3,'T_b' ,ta0  ,axis='Z')
  call put_var('append'   ,ofn3,'U_b' ,ua0  ,axis='Z')
  call put_var('append'   ,ofn3,'V_b' ,va0  ,axis='Z')
  call put_var('append'   ,ofn3,'P_b' ,pa0  ,axis='Z')
  call put_var('append'   ,ofn3,'RHO0',ra0  ,axis='Z')
  call put_var('append'   ,ofn3,'N0'  ,nbva0,axis='Z')
  call put_var('append'   ,ofn3,'T_prt',ta_prt,axis='Z')
  call put_var('append'   ,ofn3,'U_prt',ua_prt,axis='Z')
  call put_var('append'   ,ofn3,'V_prt',va_prt,axis='Z')
  call put_att(ofn3,'global','title',tit3)

END subroutine dump_3

END program get_prof

