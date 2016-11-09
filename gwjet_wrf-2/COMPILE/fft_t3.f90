PROGRAM FFT_t3

  use kl_filter
  use fft
  use easync

  implicit none

  integer            ::  i_wave, nnn(3), k_rng(2), l_rng(2)
  real               ::  lengths(2), lat_c
  character(len=16)  ::  var_name(2)
  character(len=256) ::  f_namelist, file_i, file_o

  namelist /PARAM/ I_WAVE, VAR_NAME, NNN, K_RNG, L_RNG, LENGTHS, LAT_C
  namelist /FILEIO/ FILE_I, FILE_O

  integer ::  i1, j1, n1, nt, nk, nl, no
  integer ::  i,j,n
  real    ::  dk, dl

  real,    dimension(:,:,:,:), allocatable ::  var4d
  real,    dimension(:,:,:,:), allocatable ::  var_in
  real,    dimension(:),       allocatable ::  kwn, lwn, ome
  complex, dimension(:,:,:),   allocatable ::  var_in_c
  complex, dimension(:),       allocatable ::  fc_var

! READ NAMELISTS

  call getarg(1,f_namelist)
  open(10, file=trim(f_namelist), status='old')
  read(10, PARAM)  ;  read(10, FILEIO)
  close(10)

! INITIALIZE ARRAYS and READ DATA

  call init_read

  if ( i_wave /= 0 .and. i_wave /= -999 ) then
    dk = 1./(lengths(1)*111.*cos(lat_c*0.017453293))  ! [cyc/km]
    dl = 1./(lengths(2)*111.)
    call kl_filter_case(nk,nl,dk,dl,i_wave,var_in_c)
  end if

  ! calculate FFT
  allocate( fc_var(nt) )
  do j=1, nl
  do i=1, nk
    call fft1d_f(nt,var_in_c(i,j,:),fc_var)
    var4d(i,j,:,1) = real (fc_var(:))
    var4d(i,j,:,2) = aimag(fc_var(:))
  enddo
  enddo
  deallocate( fc_var )

! DUMP and FINALIZE

  call dump
  call finalize

  STOP


  CONTAINS


  SUBROUTINE init_read

  i1 = k_rng(1)+1  ;  j1 = l_rng(1)+1  ;  n1 = nnn(1)
  nk = k_rng(2)+1 - i1 + 1  ;  nl = l_rng(2)+1 - j1 + 1  ;  nt = (nnn(2) - n1)/nnn(3) + 1

  allocate( kwn(nk) )
  call get_var(file_i,'k_wn',kwn,start=(/i1/), count=(/nk/))
  if ( kwn(1) /= float(k_rng(1)) ) then
    print*, 'k_wn miss match'  ;  STOP
  end if

  allocate( lwn(nl) )
  call get_var(file_i,'l_wn',lwn,start=(/j1/), count=(/nl/))
  if ( lwn(1) /= float(l_rng(1)) ) then
    print*, 'l_wn miss match'  ;  STOP
  end if

  no = nt
  allocate( ome(no) )
  do n=1, no
    ome(n) = float(n-1)
  enddo

  allocate( var_in_c(nk,nl,nt) )
  allocate( var_in(nk,nl,nt,2) )

  call get_var(file_i,var_name(1),var_in(:,:,:,1),start=(/i1,j1,n1/), &
               count=(/nk,nl,nt/), stride=(/1,1,nnn(3)/), map=(/1,nk,nk*nl/))
  call get_var(file_i,var_name(2),var_in(:,:,:,2),start=(/i1,j1,n1/), &
               count=(/nk,nl,nt/), stride=(/1,1,nnn(3)/), map=(/1,nk,nk*nl/))
  var_in_c = cmplx(var_in(:,:,:,1),var_in(:,:,:,2))

  deallocate( var_in )

  allocate( var4d(nk,nl,no,2) )

  END subroutine init_read

  SUBROUTINE dump

  call put_var('overwrite',file_o,'k_wn'  ,kwn,axis='k_wn'  )
  call put_var('append'   ,file_o,'l_wn'  ,lwn,axis='l_wn'  )
  call put_var('append'   ,file_o,'ome_fr',ome,axis='ome_fr')

  call put_var('append'   ,file_o,var_name(1),var4d(:,:,:,1), &
               axes=(/'k_wn','l_wn','ome_fr'/))
  call put_var('append'   ,file_o,var_name(2),var4d(:,:,:,2), &
               axes=(/'k_wn','l_wn','ome_fr'/))

  write(6,*)  ;  write(6,*) trim(file_o)  ;  write(6,*)

  END subroutine dump

  SUBROUTINE finalize

  deallocate( var4d, var_in_c )
  deallocate( kwn, lwn, ome )

  END subroutine finalize

END program FFT_t3

