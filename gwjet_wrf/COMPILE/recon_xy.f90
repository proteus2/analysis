PROGRAM RECON_xy

  use kl_filter
  use fft
  use easync

  implicit none

  integer            ::  i_wave, nnn(3), n_kl(2), k_rng(2), l_rng(2)
  real               ::  xy0(2), lengths(2), lat_c
  character(len=16)  ::  var_name
  character(len=256) ::  f_namelist, file_i, file_o

  namelist /PARAM/ I_WAVE, VAR_NAME, NNN, XY0, LENGTHS, N_KL, K_RNG, L_RNG, LAT_C
  namelist /FILEIO/ FILE_I, FILE_O

  integer ::  n1, nx, ny, nt, nk, nl, kk0, kk9, ll0, ll9
  integer ::  i,j,n
  real    ::  dk, dl

  real,    dimension(:,:,:,:), allocatable ::  var_in
  real,    dimension(:,:,:),   allocatable ::  var3d
  real,    dimension(:),       allocatable ::  lon, lat, t
  complex, dimension(:,:,:),   allocatable ::  var_kl, tmp
  complex, dimension(:,:),     allocatable ::  var_k

! READ NAMELISTS

  call getarg(1,f_namelist)
  open(10, file=trim(f_namelist), status='old')
  read(10, PARAM)  ;  read(10, FILEIO)
  close(10)

! INITIALIZE ARRAYS and READ DATA

  call init_read

!!  tmp(:,:,:) = 0.
!!  tmp   (k_rng(1)+1:k_rng(2)+1,l_rng(1)+1:l_rng(2)+1,:) =  &
!!  var_kl(k_rng(1)+1:k_rng(2)+1,l_rng(1)+1:l_rng(2)+1,:)
!!  var_kl(:,:,:) = tmp(:,:,:)

  call kl_filter_case(nk,nl,dk,dl,i_wave,var_kl)

  ! calculate FFT
  do n=1, nt
    var_k(:,:) = (0.,0.)
    do i=1, nk
      call fft1d_b(ny,var_kl(i,:,n),var_k(i,:))
    enddo
    var_k(2:nk-1,:) = var_k(2:nk-1,:)*2.  ! assume nx is even.
    do j=1, ny
      call fft1d_b(nx,var_k(:,j),var3d(:,j,n))
    enddo
  enddo

! DUMP and FINALIZE

  call dump
  call finalize

  STOP


  CONTAINS


  SUBROUTINE init_read

  n1 = nnn(1)
  nk = n_kl(1)  ;  nl = n_kl(2)  ;  nt = (nnn(2) - n1)/nnn(3) + 1

  allocate( t(nt) )
  call get_var(file_i,'t',t,start=(/n1/), count=(/nt/), stride=(/nnn(3)/))

  nx = (nk-1)*2
  ny = nl
  allocate( lon(nx) )
  do i=1, nx
    lon(i) = xy0(1) + float(i-1)*(lengths(1)/nx)
  enddo
  allocate( lat(ny) )
  do j=1, ny
    lat(j) = xy0(2) + float(j-1)*(lengths(2)/ny)
  enddo

  dk = 1./(lengths(1)*111.*cos(lat_c*0.017453293))  ! [cyc/km]
  dl = 1./(lengths(2)*111.)

  allocate( var_kl(nk,nl,nt), tmp(nk,nl,nt) )
  allocate( var_in(nk,nl,nt,2) )

  call get_var(file_i,'fc_'//trim(var_name)//'_r',var_in(:,:,:,1),          &
               start=(/1,1,n1/), count=(/nk,nl,nt/), stride=(/1,1,nnn(3)/), &
               map=(/1,nk,nk*nl/))
  call get_var(file_i,'fc_'//trim(var_name)//'_i',var_in(:,:,:,2),          &
               start=(/1,1,n1/), count=(/nk,nl,nt/), stride=(/1,1,nnn(3)/), &
               map=(/1,nk,nk*nl/))
  var_kl = cmplx(var_in(:,:,:,1),var_in(:,:,:,2))

  deallocate( var_in )

  allocate( var_k(nx,nl) )
  allocate( var3d(nx,ny,nt) )

  END subroutine init_read

  SUBROUTINE dump

  call put_var('overwrite',file_o,'lon',lon,axis='lon')
  call put_var('append'   ,file_o,'lat',lat,axis='lat')
  call put_var('append'   ,file_o,'t'  ,t  ,axis='t'  )

  call put_var('append'   ,file_o,var_name,var3d, &
               axes=(/'lon','lat','t'/))

  write(6,*)  ;  write(6,*) trim(file_o)  ;  write(6,*)

  END subroutine dump

  SUBROUTINE finalize

  deallocate( var3d, var_kl, var_k, tmp )
  deallocate( lon, lat, t )

  END subroutine finalize

END program RECON_xy

