PROGRAM INTERPOL_DATA

  use fft
  use easync

  implicit none

  integer            ::  nd(4), it(2), nxp
  character(len=16)  ::  var_names(99)
  character(len=256) ::  f_namelist, file_i(2), file_o(99)

  namelist /PARAM/ VAR_NAMES, ND, IT, NXP
  namelist /FILEIO/ FILE_I, FILE_O

  integer ::  nx, ny, nz, nn, nn1, nn2, nxp2
  integer ::  j,k,n, iv

  real,    dimension(:,:,:,:), allocatable ::  var_in1, var_in2, var4d
  real,    dimension(:,:,:),   allocatable ::  var_in1p, var_in2p
  real,    dimension(:,:),     allocatable ::  coef
  real,    dimension(:),       allocatable ::  abs_fc1, abs_fc2, abs_fc, abs_fc0
  real,    dimension(:),       allocatable ::  var1d
  complex, dimension(:),       allocatable ::  fc_var1, fc_var2, fc_var

! READ NAMELISTS

  call getarg(1,f_namelist)
  open(10, file=trim(f_namelist), status='old')
  read(10, PARAM)  ;  read(10, FILEIO)
  close(10)

  nn1 = it(1) - 1  ;  nn2 = it(1) + it(2)
  allocate( coef(it(2),2) )
  do n=1, it(2)
    coef(n,2) = float(n)/float(it(2) + 1)
    coef(n,1) = 1. - coef(n,2)
  enddo
  if (nn1 == 0   )  nn1 = nd(4)
  if (nn2 > nd(4))  nn2 = nn2 - nd(4)

  VARS:  DO iv=1, 99

  if ( trim(var_names(iv)) == '-999' )  exit

! INITIALIZE ARRAYS and READ DATA

  call init_read

  ! make the data periodic in x direction
  if (mod(nx,nxp) == 0) then
    nxp2 = nx
    allocate( var_in1p(nx,ny,nz), var_in2p(nx,ny,nz) )
    var_in1p(:,:,:) = var_in1(:,:,:,1)
    var_in2p(:,:,:) = var_in2(:,:,:,1)
  else
    nxp2 = nx + (nxp - mod(nx,nxp))
    allocate( var_in1p(nxp2,ny,nz), var_in2p(nxp2,ny,nz) )
    var_in1p(1:nx,:,:) = var_in1(:,:,:,1)
    var_in2p(1:nx,:,:) = var_in2(:,:,:,1)
    var_in1p(nx+1:nxp2,:,:) = var_in1(nx+1-nxp:nxp2-nxp,:,:,1)
    var_in2p(nx+1:nxp2,:,:) = var_in2(nx+1-nxp:nxp2-nxp,:,:,1)
  end if
  deallocate( var_in1, var_in2 )

  ! calculate FFT
  allocate( fc_var1(nxp2), fc_var2(nxp2), fc_var(nxp2) )
  allocate( abs_fc1(nxp2), abs_fc2(nxp2), abs_fc(nxp2), abs_fc0(nxp2) )
  allocate( var1d(nxp2) )
  do k=1, nz
  do j=1, ny
    call fft1d_f(nxp2,var_in1p(:,j,k),fc_var1)
    call fft1d_f(nxp2,var_in2p(:,j,k),fc_var2)
    abs_fc1(:) = cabs(fc_var1(:))
    abs_fc2(:) = cabs(fc_var2(:))
    do n=1, it(2)
      fc_var(:) = coef(n,1)*fc_var1(:) + &
                  coef(n,2)*fc_var2(:)
      abs_fc(:) = coef(n,1)*abs_fc1(:) + &
                  coef(n,2)*abs_fc2(:)
      abs_fc0(:) = cabs(fc_var(:))
      where ( abs_fc0 /= 0. )
        fc_var(:) = (fc_var(:)/abs_fc0(:))*abs_fc(:)
      end where
      call fft1d_b(nxp2,fc_var,var1d)
      var4d(:,j,k,n) = var1d(1:nx)
    enddo
  enddo
  enddo
  deallocate( var1d )
  deallocate( abs_fc1, abs_fc2, abs_fc, abs_fc0 )
  deallocate( fc_var1, fc_var2, fc_var )
  deallocate( var_in1p, var_in2p )

! DUMP and FINALIZE

  call dump
  deallocate( var4d )

  ENDDO  VARS

  STOP


  CONTAINS


  SUBROUTINE init_read

  nx = nd(1)  ;  ny = nd(2)  ;  nz = nd(3)
  if (var_names(iv) == 'U') then
    nx = nx + 1
  end if
  if (var_names(iv) == 'V') then
    ny = ny + 1
  end if
  if (var_names(iv) == 'W') then
    nz = nz + 1
  end if
  if (var_names(iv) == 'PH') then
    nz = nz + 1
  end if

  allocate( var_in1(nx,ny,nz,1), var_in2(nx,ny,nz,1) )

  nn = it(1) - 1
  if (nn == 0)  nn = nd(4)
  call get_var(file_i(1),var_names(iv),var_in1,start=(/1,1,1,nn/), &
               count=(/nx,ny,nz,1/))
  nn = it(1) + it(2)
  if (nn > nd(4))  nn = nn - nd(4)
  call get_var(file_i(2),var_names(iv),var_in2,start=(/1,1,1,nn/), &
               count=(/nx,ny,nz,1/))

  allocate( var4d(nx,ny,nz,it(2)) )

  END subroutine init_read

  SUBROUTINE dump

  do n=1, it(2)
    nn = it(1) + (n-1)
    if (nn > nd(4))  nn = nn - nd(4)
    call put_var('append',file_o(n),trim(var_names(iv)), &
                 reshape(var4d(:,:,:,n),(/nx,ny,nz,1/)), &
                 start=(/1,1,1,nn/),count=(/nx,ny,nz,1/))
    write(6,*)  ;  write(6,*) trim(file_o(n))  ;  write(6,*)
  enddo

  END subroutine dump

END program INTERPOL_DATA

