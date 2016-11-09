PROGRAM VERT_INTERPOL_TO_Z

  use nr, only: spline, splint
  use easync

  implicit none

  integer            ::  n1, kk(2), z_out(3)
  real               ::  z_divided_by, h_scale  ! h_scale in [km]
  character(len=16)  ::  var_name, grid_i, z_name
  character(len=256) ::  f_namelist, file_i_h, file_i_t, file_z_h, file_z_t, &
                         file_d, file_o_h, file_o_t

  namelist /PARAM/ VAR_NAME, GRID_I, N1, KK, Z_NAME, Z_DIVIDED_BY, Z_OUT, &
                   H_SCALE
  namelist /FILEIO/ FILE_I_H, FILE_I_T, FILE_Z_H, FILE_Z_T, FILE_D, FILE_O_H, &
                    FILE_O_T

  integer ::  nx, ny, nz, nz0
  integer ::  i,j,k
  real    ::  t(1)
  character(len=256) ::  file_i, file_z, file_o

  real, dimension(:,:,:,:), allocatable ::  var_out
  real, dimension(:,:,:),   allocatable ::  var_in, z_in
  real, dimension(:),       allocatable ::  lon, lat, z
  real, dimension(:),       allocatable ::  work

!  real, parameter ::  r_earth = 6371.e3

! READ NAMELISTS

  call getarg(1,f_namelist)
  open(10, file=trim(f_namelist), status='old')
  read(10, PARAM)  ;  read(10, FILEIO)
  close(10)

! INITIALIZE ARRAYS and READ DATA

  call init_read

  if ( trim(var_name) == 'P' .or. trim(var_name) == 'T' ) then
    var_in(:,:,:) = log(var_in(:,:,:))
  end if

  if ( trim(grid_i) == 'W' ) then
    allocate( work(0:nz0) )
  else
    allocate( work(nz0) )
  end if

  Y_LOOP:  DO j=1, ny
  X_LOOP:  DO i=1, nx

  call spline(z_in(:,i,j),var_in(:,i,j),1.e32,1.e32,work)

  do k=1, nz
    var_out(i,j,k,1) = splint(z_in(:,i,j),var_in(:,i,j),work,z(k))
  enddo

  ENDDO  X_LOOP
  ENDDO  Y_LOOP

  deallocate( work )

  if ( trim(var_name) == 'P' .or. trim(var_name) == 'T' ) then
    var_out(:,:,:,1) = exp(var_out(:,:,:,1))
  end if
  ! for WRF
  if ( trim(var_name) == 'T' )  var_out(:,:,:,1) = var_out(:,:,:,1) - 300.

! DUMP and FINALIZE

  call dump
  call finalize

  STOP


  CONTAINS


  SUBROUTINE init_read

  integer ::  tmpi(3)
  real, dimension(:,:,:), allocatable ::  var_in0, z_in0
  real, dimension(:,:),   allocatable ::  tmp2d

  call inq_var(file_d,'XLONG', var_shape=tmpi)
  nx = tmpi(1)
  ny = tmpi(2)

  allocate( lon(nx), lat(ny) )
  allocate( tmp2d(nx,ny) )
  call get_var(file_d,'XLONG',tmp2d)
  lon(:) = tmp2d(:,1)
  call get_var(file_d,'XLAT' ,tmp2d)
  lat(:) = tmp2d(1,:)
  deallocate( tmp2d )
  t(1) = float(n1-1)

  nz0 = kk(2) - kk(1) + 1

  if ( trim(grid_i) == 'W' ) then
    allocate( var_in(0:nz0,nx,ny) )
  else
    allocate( var_in(nz0,nx,ny) )
  end if

  select case ( trim(grid_i) )
    case ( 'U' )
      allocate( var_in0(nx+1,ny,1) )
      do k=1, nz0
        write(file_i,'(a,i3.3,a)') trim(file_i_h), kk(1) + k - 1, trim(file_i_t)
        call get_var(file_i,var_name,var_in0, start=(/1,1,n1/), count=(/nx+1,ny,1/), &
                     map=(/1,nx+1,(nx+1)*ny/))
        var_in(k,:,:) = 0.5*(var_in0(1:nx,:,1) + var_in0(2:nx+1,:,1))
      enddo
      deallocate( var_in0 )
    case ( 'V' )
      allocate( var_in0(nx,ny+1,1) )
      do k=1, nz0
        write(file_i,'(a,i3.3,a)') trim(file_i_h), kk(1) + k - 1, trim(file_i_t)
        call get_var(file_i,var_name,var_in0, start=(/1,1,n1/), count=(/nx,ny+1,1/), &
                     map=(/1,nx,nx*(ny+1)/))
        var_in(k,:,:) = 0.5*(var_in0(:,1:ny,1) + var_in0(:,2:ny+1,1))
      enddo
      deallocate( var_in0 )
    case ( 'W' )
      allocate( var_in0(nx,ny,1) )
      do k=0, nz0
        write(file_i,'(a,i3.3,a)') trim(file_i_h), kk(1) + k - 1, trim(file_i_t)
        call get_var(file_i,var_name,var_in0, start=(/1,1,n1/), count=(/nx,ny,1/), &
                     map=(/1,nx,nx*ny/))
        var_in(k,:,:) = var_in0(:,:,1)
      enddo
      deallocate( var_in0 )
    case default
      allocate( var_in0(nx,ny,1) )
      do k=1, nz0
        write(file_i,'(a,i3.3,a)') trim(file_i_h), kk(1) + k - 1, trim(file_i_t)
        call get_var(file_i,var_name,var_in0, start=(/1,1,n1/), count=(/nx,ny,1/), &
                     map=(/1,nx,nx*ny/))
        var_in(k,:,:) = var_in0(:,:,1)
      enddo
      deallocate( var_in0 )
  end select
  if ( trim(var_name) == 'T' ) then
    var_in(:,:,:) = var_in(:,:,:) + 300.  ! for WRF
  end if

  allocate( z_in0(0:nz0,nx,ny) )
  allocate( var_in0(nx,ny,1) )
  do k=0, nz0
    write(file_i,'(a,i3.3,a)') trim(file_z_h), kk(1) + k - 1, trim(file_z_t)
    call get_var(file_i,z_name,var_in0, start=(/1,1,n1/), count=(/nx,ny,1/), &
                 map=(/1,nx,nx*ny/))
    z_in0(k,:,:) = var_in0(:,:,1)
  enddo
  deallocate( var_in0 )
  z_in0(:,:,:) = z_in0(:,:,:)/z_divided_by
  z_in0(:,:,:) = z_in0(:,:,:)*1.e-3

  if ( trim(grid_i) == 'W' ) then
    allocate( z_in(0:nz0,nx,ny) )
    z_in(:,:,:) = z_in0(:,:,:)
  else
    allocate( z_in(nz0,nx,ny) )
    z_in0(:,:,:) = exp(z_in0(:,:,:)/(-h_scale))
    z_in(:,:,:) = log(0.5*(z_in0(0:nz0-1,:,:) + z_in0(1:nz0,:,:)))*(-h_scale)
  end if
  deallocate( z_in0 )

  nz = (z_out(2) - z_out(1))/z_out(3) + 1
  allocate( z(nz) )
  do k=1, nz
    z(k) = float(z_out(1)) + float(z_out(3))*float(k-1)
  enddo
  z(:) = z(:)*1.e-3

  allocate( var_out(nx,ny,nz,1) )

  END subroutine init_read

  SUBROUTINE dump

  do k=1, nz

    file_o = ''
    if (int(z(k)*1.e3) >= 0) then
      write(file_o,'(a,i5.5,a)') trim(file_o_h), int(z(k)*1.e3), trim(file_o_t)
    else
      write(file_o,'(a4,i1.1)') '_-00', abs(int(z(k)*1.e3))
      if (z(k) <= -0.01)  write(file_o,'(a3,i2.2)') '_-0', abs(int(z(k)*1.e3))
      if (z(k) <= -0.1 )  write(file_o,'(a2,i3.3)') '_-' , abs(int(z(k)*1.e3))
      file_o = trim(file_o_h)//trim(file_o)//trim(file_o_t)
    end if
 
    call put_var('overwrite',file_o,'lon',lon,axis='lon')
    call put_var('append'   ,file_o,'lat',lat,axis='lat')
    call put_var('append'   ,file_o,'t'  ,t  ,axis='t'  ,is_record=.TRUE.)

    call put_var('append'   ,file_o,trim(var_name),var_out(:,:,k,:), &
                 is_record=.TRUE., axes=(/'lon','lat','t'/))

    if (k == 1) then
      write(6,*)  ;  write(6,*) trim(file_o)
    end if
    if (k == nz) then
      write(6,*) trim(file_o)  ;  write(6,*)
    end if

  enddo

  END subroutine dump

  SUBROUTINE finalize

  deallocate( var_in, var_out, z_in )
  deallocate( lon, lat, z )

  END subroutine finalize

END program VERT_INTERPOL_TO_Z

