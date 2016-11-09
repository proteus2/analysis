PROGRAM anal_gw

!  NAMELISTS
!-----------------------------------------------------------------------
! fi_dir  : base directory in which the input file exists
! fi_head : head of the input file
! fi_tail : tail of the input file
! fo_dir  : base directory in which the output file will be saved
! fo_head : head of the output file
! fo_tail : tail of the output file
!-----------------------------------------------------------------------
! zb      : bottom height for analysis.
! zt      : top height for analysis.
!-----------------------------------------------------------------------

  use m_etc
  use igw_char_prof
  use fft
  use easync
  
  implicit none

! namelist parameters
  real               ::  zb, zt
  character(len=256) ::  fi_dir, fo_dir
  character(len=64)  ::  fi_head, fi_tail, fo_head, fo_tail

  namelist /IO_PARAM/  fi_dir, fi_head, fi_tail, fo_dir, fo_head, fo_tail
  namelist /ANALYSIS_PARAM/  zb, zt

! main
  integer                            ::  nz, nwav
  real, dimension(:), allocatable    ::  z, tn_a, u_a, v_a,              &
                                         temp0, u0, v0, r0, n0,          &
                                         u0s_z, v0s_z, wgt, wf
  real, dimension(:), allocatable    ::  tn, u, v, tn_hil
  real                               ::  dz, ph_dir, r_ax, d_polar
  real                               ::  ome_i, m_wn, k_wn
  real                               ::  u0s, v0s, n0m
  real                               ::  sp_i, sp_d, sp_p, sp_q
  character(len=256)                 ::  ifn, ofn
  character(len=256)                 ::  tit
  character(len=16)                  ::  varn
  complex, dimension(:), allocatable ::  fc

  integer                         ::  k, iw

! constants
  real, parameter ::  rd = 287.0, cp = 1004.0, g = 9.806
  real, parameter ::  deg2rad = 3.14159265358979323846/180.

!-----------------------------------------------------------------------
!  GET ARGUMENTS AND READ NAMELISTS
!-----------------------------------------------------------------------

  call input_arg
  call input_namelist

!-----------------------------------------------------------------------
!  INPUT DATA FILENAME
!-----------------------------------------------------------------------

  write(ifn,'(a,7a)') trim(fi_dir)//'/'//trim(stid)//'/',  &
       trim(fi_head),trim(year),month,date,hour,trim(fi_tail),'.nc'
  write(6,'(/,a)') ' input : '//trim(ifn)

!-----------------------------------------------------------------------
!  OUTPUT FILENAME
!-----------------------------------------------------------------------

  write(tit,'(a)') 'IGW characteristics'

  write(ofn,'(a,7a)') trim(fo_dir)//'/'//trim(stid)//'/',  &
       trim(fo_head),trim(year),month,date,hour,trim(fo_tail),'.nc'
  write(6,'(/,a)') ' output : '//trim(ofn)

!-----------------------------------------------------------------------
!  READ INPUT DATA
!-----------------------------------------------------------------------

  call read_data
  ! out: z, tn_a, u_a, v_a, temp0, u0, v0, r0, n0 with dz and nz

!-----------------------------------------------------------------------
!  CALCULATE BACKGROUND VARIABLES
!-----------------------------------------------------------------------

  allocate( u0s_z(nz), v0s_z(nz) )
  u0s_z(2:nz-1) = (u0(3:nz) - u0(1:nz-2))/(z(3:nz) - z(1:nz-2))
  v0s_z(2:nz-1) = (v0(3:nz) - v0(1:nz-2))/(z(3:nz) - z(1:nz-2))

  ! vertical mean
  u0s = (u0(nz) - u0(1))/(z(nz) - z(1))
  v0s = (v0(nz) - v0(1))/(z(nz) - z(1))
  n0m = sum(n0)/float(nz)


!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  nwav = 1

  N_WAVE: DO iw=0, nwav

  if (iw == 0) then
    tn(:) = tn_a(:)
    u (:) = u_a (:)
    v (:) = v_a (:)
  end if

  if (iw > 0) then
    tn(:) = 0.  ;  u(:) = 0.  ;  v(:) = 0.
!-----------------------------------------------------------------------
!  DECOMPOSE MULTI-WAVES
!-----------------------------------------------------------------------
    ! for example
    if (iw == 1) then
      allocate( fc(nz) )
      call fft1d_f(tn_a, fc)
      fc(nz/10+2:nz-nz/10) = 0.
      call fft1d_b(fc, tn)
      call fft1d_f(u_a, fc)
      fc(nz/10+2:nz-nz/10) = 0.
      call fft1d_b(fc, u)
      call fft1d_f(v_a, fc)
      fc(nz/10+2:nz-nz/10) = 0.
      call fft1d_b(fc, v)
      deallocate( fc )
    end if
  end if

!-----------------------------------------------------------------------
!  CALCULATE MEAN VERTICAL WAVENUMBER
!-----------------------------------------------------------------------

  allocate( wf(nz) )
  call welch_window(nz, wf)
  call waveno_m_avg(u,v,dz,wf, m_wn)
  deallocate( wf )

!-----------------------------------------------------------------------
!  CALCULATE STOKES PARAMETERS AND DEGREE OF POLARIZATION
!-----------------------------------------------------------------------

  call stokes_param(u,v, sp_i,sp_d,sp_p,sp_q)

  call degree_polar(sp_i,sp_d,sp_p,sp_q, d_polar)

!-----------------------------------------------------------------------
!  OBTAIN THE WAVE PROPAGATION DIRECTION
!-----------------------------------------------------------------------

  call phase_dir(u,v,tn, ph_dir)

!-----------------------------------------------------------------------
!  CALCULATE VERTICAL-MEAN BACKGROUND VARIABLES
!-----------------------------------------------------------------------

  allocate( wgt(nz) )
  wgt(:) = 1.

  ! weighting by (horizontal) kinetic energy
!  wgt(:) = u(:)*u(:) + v(:)*v(:)

  ! weighting by absolute momentum flux
!  allocate( tn_hil(nz) )
!  call hilbert_transform(tn, tn_hil)
!  wgt(:) = sqrt((u(:)*u(:) + v(:)*v(:)))*abs(tn_hil(:))
!  deallocate( tn_hil )

  if ( any(wgt(:) /= 1.) ) then
    u0s = sum(u0s_z(2:nz-1)*wgt(2:nz-1))/sum(wgt(2:nz-1))
    v0s = sum(v0s_z(2:nz-1)*wgt(2:nz-1))/sum(wgt(2:nz-1))
    n0m = sum(n0(:)*wgt(:))/sum(wgt)
  end if

  deallocate( wgt )

!-----------------------------------------------------------------------
!  OBTAIN THE INTRINSIC FREQUENCY
!-----------------------------------------------------------------------

  call intr_freq(sp_d,sp_p,sp_q,ph_dir,1.e-4,u0s,v0s,n0m, ome_i)

!-----------------------------------------------------------------------
!  DUMP
!-----------------------------------------------------------------------

  call dump_1

  ENDDO  N_WAVE

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


  deallocate( z, temp0, u0, v0, r0, n0, tn_a, u_a, v_a )
  deallocate( u0s_z, v0s_z )
  deallocate( tn, u, v )
 
  STOP


  CONTAINS


SUBROUTINE input_namelist

  open(10, file=trim(f_namelist), status='old')
  read(10, IO_PARAM)  ;  read(10, ANALYSIS_PARAM)
  close(10)

END subroutine input_namelist

SUBROUTINE read_data
 
  integer ::  nza(1), kzb, kzt
  real, dimension(:), allocatable ::  za

  call inq_var(ifn,'Z',var_shape=nza)
  allocate( za(nza(1)) )
  call get_var(ifn,'Z',za)
  dz = za(2) - za(1)

  if ( za(1) > zb .or. za(nza(1)) < zt )  call stop_message(  &
     '[zb,zt] is not fully included in the input profile.')
  if ( any( (abs((za(3:nza(1))-za(2:nza(1)-1))-dz)/dz) > 0.01 ) )  &
     call stop_message('The data interval is not a constant.')
 
  do k=1, nza(1)-1
    if (za(k) <= zb)  kzb = k
    if (za(k) < zt )  kzt = k+1
  enddo
  nz = kzt - kzb + 1

  deallocate( za )

  allocate( z(nz) )
  allocate( temp0(nz), u0(nz), v0(nz), r0(nz), n0(nz) )
  allocate( tn_a(nz), u_a(nz), v_a(nz) )
  call get_var(ifn,'Z'    ,z    , start=(/kzb/), count=(/nz/))
  call get_var(ifn,'T_b'  ,temp0, start=(/kzb/), count=(/nz/))
  call get_var(ifn,'U_b'  ,u0   , start=(/kzb/), count=(/nz/))
  call get_var(ifn,'V_b'  ,v0   , start=(/kzb/), count=(/nz/))
  call get_var(ifn,'RHO0' ,r0   , start=(/kzb/), count=(/nz/))
  call get_var(ifn,'N0'   ,n0   , start=(/kzb/), count=(/nz/))
  call get_var(ifn,'T_prt',tn_a , start=(/kzb/), count=(/nz/))
  call get_var(ifn,'U_prt',u_a  , start=(/kzb/), count=(/nz/))
  call get_var(ifn,'V_prt',v_a  , start=(/kzb/), count=(/nz/))

  tn_a(:) = tn_a(:)/temp0(:)

  allocate( tn(nz), u(nz), v(nz) )

END subroutine read_data

SUBROUTINE dump_1

  if (iw == 0) then
    call put_var('overwrite',ofn,'WAVE',(/iw/), is_record=.True.)
    call put_var('append',ofn,'Z',z)
  else
    call put_var('append',ofn,'WAVE',(/iw/), is_record=.True.,  &
                 start=(/iw+1/))
  end if
  call put_var('append',ofn,'m_wn'   ,(/m_wn   /), axis='WAVE',  &
               is_record=.True., start=(/iw+1/))
  call put_var('append',ofn,'d_polar',(/d_polar/), axis='WAVE',  &
               is_record=.True., start=(/iw+1/))
  call put_var('append',ofn,'ph_dir' ,(/ph_dir /), axis='WAVE',  &
               is_record=.True., start=(/iw+1/))
  call put_var('append',ofn,'ome_i'  ,(/ome_i  /), axis='WAVE',  &
               is_record=.True., start=(/iw+1/))
  call put_var('append',ofn,'u0s'    ,(/u0s    /), axis='WAVE',  &
               is_record=.True., start=(/iw+1/))
  call put_var('append',ofn,'v0s'    ,(/v0s    /), axis='WAVE',  &
               is_record=.True., start=(/iw+1/))
  call put_var('append',ofn,'n0m'    ,(/n0m    /), axis='WAVE',  &
               is_record=.True., start=(/iw+1/))
  
  write(varn,'(a,i0)')  'u_', iw
  call put_var('append',ofn,varn,u, axis='Z')
  write(varn,'(a,i0)')  'v_', iw
  call put_var('append',ofn,varn,v, axis='Z')
  write(varn,'(a,i0)')  't_', iw
  call put_var('append',ofn,varn,tn(:)*temp0(:), axis='Z')

  call put_att(ofn,'global','title',tit)

END subroutine dump_1


END program anal_gw

