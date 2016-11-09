PROGRAM PERTURBATION_EXTRACTION

  use easync

  implicit none

  integer            ::  nnn(3), jj(2), ii(2), n_dom
  real               ::  d_crit
  character(len=16)  ::  var_name, grid_i
  character(len=256) ::  f_namelist, file_i, file_d, file_o

  namelist /PARAM/ VAR_NAME, GRID_I, NNN, JJ, II, N_DOM, D_CRIT
  namelist /FILEIO/ FILE_I, FILE_D, FILE_O

  integer ::  i1, j1, n1, nx, ny, nt, jbuf, is, ie
  integer ::  i,j,n, id,j2, j3, nj2, j2_i
  real    ::  cosdang_crit, dlon, dlat

  real, dimension(:,:,:), allocatable ::  var_in, prt, mean
  real, dimension(:,:),   allocatable ::  cosdang
  real, dimension(:),     allocatable ::  lon, lat, t
  real, dimension(:),     allocatable ::  coslat, sinlat, cosdlon
  real, dimension(:),     allocatable ::  norm, wgt, wgtl
  integer, dimension(:,:), allocatable ::  i_dlon

  real, parameter ::  r_earth = 6371.e3
  real, parameter ::  deg2rad = 0.017453292519943295769

! READ NAMELISTS

  call getarg(1,f_namelist)
  open(10, file=trim(f_namelist), status='old')
  read(10, PARAM)  ;  read(10, FILEIO)
  close(10)

! INITIALIZE ARRAYS and READ DATA

  call init_read

  cosdang_crit = cos(d_crit*1.e3/r_earth)

  ! calculate mean
  allocate( coslat(ny), sinlat(ny), cosdlon(0:nx*n_dom/2) )
  allocate( i_dlon(1+jbuf:ny-jbuf,ny), cosdang(1+jbuf:ny-jbuf,ny) )

  coslat(:) = cos(lat(:)*deg2rad)
  sinlat(:) = sin(lat(:)*deg2rad)
  do i=0, nx*n_dom/2
    cosdlon(i) = cos(float(i)*dlon*deg2rad)
  enddo

  i_dlon = -999
  do id=0, nx*n_dom/2
    do j2=1, ny
      cosdang(:,j2) = coslat(1+jbuf:ny-jbuf)*coslat(j2)*cosdlon(id) + &
                      sinlat(1+jbuf:ny-jbuf)*sinlat(j2)
    enddo
    where ( cosdang >= cosdang_crit )
      i_dlon(:,:) = id
    end where
  enddo

  allocate( mean(nx,1+jbuf:ny-jbuf,nt), norm(1+jbuf:ny-jbuf) )
  allocate( wgt(nx*2), wgtl(1-nx*(n_dom/2+1):nx*(n_dom/2+2)) )

  j2_i = ny+1  ;  nj2 = 0
  do j2=1, ny
    if ( i_dlon(1+jbuf,j2) /= -999 ) then
      j2_i = j2  ;  EXIT
    end if
  enddo
  do j2=j2_i, ny
    if ( i_dlon(1+jbuf,j2) /= -999 )  nj2 = nj2 + 1
  enddo
  
  mean(:,:,:) = 0.
  norm(:) = 0.
  do j=1+jbuf, ny-jbuf
!  do j2=j2_i+(j-1-jbuf), j2_i+(j-1-jbuf)+nj2-1
  do j3=1, nj2
    j2 = j2_i + (j3-1) + (j-1-jbuf)
    is = 1 - i_dlon(j,j2)
    ie = 1 + i_dlon(j,j2)
    wgt = 0.  ;  wgtl = 0.
    wgtl(is:ie) = coslat(j2)
    do id=-n_dom/2+1, n_dom/2+1
      wgt(1:nx) = wgt(1:nx) + wgtl(1+nx*id:nx*(id+1))
    enddo
    wgt(nx+1:nx*2) = wgt(1:nx)
    do i=1, nx
      mean(i,j,:) = mean(i,j,:) +  &
          sum(var_in(1:nx,j2,:)*spread(wgt(nx+2-i:nx*2+1-i),2,nt), dim=1)
    enddo
    norm(j) = norm(j) + sum(wgt(1:nx))
  enddo
  enddo
  do j=1+jbuf, ny-jbuf
    mean(:,j,:) = mean(:,j,:)/norm(j)
  enddo
  prt(:,:,:) = var_in(:,1+jbuf:ny-jbuf,:) - mean(:,:,:)

  deallocate( i_dlon, cosdang )
  deallocate( coslat, sinlat, cosdlon )

  deallocate( norm )

! DUMP and FINALIZE

  call dump
  call finalize

  STOP


  CONTAINS


  SUBROUTINE init_read

  real, dimension(:,:,:), allocatable ::  var_in0
  real, dimension(:,:),   allocatable ::  tmp2d

  allocate( tmp2d(2,1) )
  call get_var(file_d,'XLONG',tmp2d, start=(/1,1/), count=(/2,1/))
  dlon = tmp2d(2,1) - tmp2d(1,1)
  call get_var(file_d,'XLAT', tmp2d, start=(/1,1/), count=(/1,2/), map=(/1,1/))
  dlat = tmp2d(2,1) - tmp2d(1,1)
  deallocate( tmp2d )
  jbuf = int(d_crit*1.e3/r_earth/deg2rad/dlat) + 1

  i1 = ii(1)  ;  n1 = nnn(1)
  nx = ii(2) - i1 + 1  ;  nt = (nnn(2) - n1)/nnn(3) + 1
  j1 = jj(1) - jbuf
  ny = jj(2) + jbuf - j1 + 1

  allocate( lon(nx), lat(ny), t(nt) )

  allocate( tmp2d(nx,ny) )
  call get_var(file_d,'XLONG',tmp2d, start=(/i1,j1/), count=(/nx,ny/))
  lon(:) = tmp2d(:,1)
  call get_var(file_d,'XLAT' ,tmp2d, start=(/i1,j1/), count=(/nx,ny/))
  lat(:) = tmp2d(1,:)
  deallocate( tmp2d )
  do n=1, nt
    t(n) = float(n1 + (n - 1)*nnn(3))
  enddo

  allocate( prt(nx,1+jbuf:ny-jbuf,nt), var_in(nx,ny,nt) )

  select case ( trim(grid_i) )
    case ( 'U' )
      allocate( var_in0(nx+1,ny,nt) )
      call get_var(file_i,var_name,var_in0, start=(/i1,j1,n1/), count=(/nx+1,ny,nt/), &
                   stride=(/1,1,nnn(3)/), map=(/1,nx+1,(nx+1)*ny/))
      var_in(:,:,:) = 0.5*(var_in0(1:nx,:,:) + var_in0(2:nx+1,:,:))
      deallocate( var_in0 )
    case ( 'V' )
      allocate( var_in0(nx,ny+1,nt) )
      call get_var(file_i,var_name,var_in0, start=(/i1,j1,n1/), count=(/nx,ny+1,nt/), &
                   stride=(/1,1,nnn(3)/), map=(/1,nx,nx*(ny+1)/))
      var_in(:,:,:) = 0.5*(var_in0(:,1:ny,:) + var_in0(:,2:ny+1,:))
      deallocate( var_in0 )
    case default
      call get_var(file_i,var_name,var_in, start=(/i1,j1,n1/), count=(/nx,ny,nt/), &
                   stride=(/1,1,nnn(3)/), map=(/1,nx,nx*ny/))
  end select

  END subroutine init_read

  SUBROUTINE dump

  call put_var('overwrite',file_o,'lon',lon,axis='lon')
  call put_var('append'   ,file_o,'lat',lat(1+jbuf:ny-jbuf),axis='lat')
  call put_var('append'   ,file_o,'t'  ,t  ,axis='t'  )

  call put_var('append'   ,file_o,'prt_'//trim(var_name),prt, &
               axes=(/'lon','lat','t'/))
  call put_var('append'   ,file_o,'mean_'//trim(var_name),mean, &
               axes=(/'lon','lat','t'/))

  write(6,*)  ;  write(6,*) trim(file_o)  ;  write(6,*)

  END subroutine dump

  SUBROUTINE finalize

  deallocate( prt, mean, var_in )
  deallocate( lon, lat, t )

  END subroutine finalize

END program PERTURBATION_EXTRACTION

