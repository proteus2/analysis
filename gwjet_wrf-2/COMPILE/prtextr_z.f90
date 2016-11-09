PROGRAM PERTURBATION_EXTRACTION

  use easync

  implicit none

  integer            ::  nnn(3), jj(2), ii(2)
  real               ::  d_crit
  character(len=16)  ::  var_name
  character(len=256) ::  f_namelist, file_i, file_o(2)

  namelist /PARAM/ VAR_NAME, NNN, JJ, II, D_CRIT
  namelist /FILEIO/ FILE_I, FILE_O

  integer ::  i1, j1, n1, nx, ny, nt, jbuf
  integer ::  i,j,n, is, i2, j2, j2_i, j2_e
  real    ::  dlon, dlat, norm

  real, dimension(:,:,:), allocatable ::  var_in, prt, mean
  real, dimension(:,:),   allocatable ::  wgt, wgt0
  real, dimension(:),     allocatable ::  lon, lat, t
  real, dimension(:),     allocatable ::  coslat, sinlat, sin2dlonp2

  real, parameter ::  d_kernel = 20.   ! [km]
  real, parameter ::  r_earth = 6371.  ! [km]
  real, parameter ::  deg2rad = 0.017453292519943295769
  real, parameter ::  pi05 = 1.5707963267948966192

! READ NAMELISTS

  call getarg(1,f_namelist)
  open(10, file=trim(f_namelist), status='old')
  read(10, PARAM)  ;  read(10, FILEIO)
  close(10)

! INITIALIZE ARRAYS and READ DATA

  call init_read

  ! calculate mean
  allocate( coslat(ny), sinlat(ny), sin2dlonp2(0:nx/2) )
  allocate( wgt0(0:nx/2,-jbuf:jbuf) )
  allocate( wgt(nx*2,-jbuf:jbuf) )

  coslat(:) = cos(lat(:)*deg2rad)
  sinlat(:) = sin(lat(:)*deg2rad)
  do i2=0, nx/2
    sin2dlonp2(i2) = sin(0.5*float(i2)*dlon*deg2rad)**2
  enddo

  mean(:,:,:) = 0.

  J_LAT:  DO j=1+jbuf, ny-jbuf

  do j2=-jbuf, jbuf
    wgt0(:,j2) = 2.*asin( sqrt( sin(0.5*abs(lat(j)-lat(j+j2))*deg2rad)**2 +  &
                                coslat(j)*coslat(j+j2)*sin2dlonp2(:) ) )
  enddo
  wgt0(:,:) = (d_crit - wgt0(:,:)*r_earth)/d_kernel*pi05
  
  wgt0(:,:) = max(-pi05, min(pi05 , wgt0(:,:)) )

  wgt0(:,:) = 0.5*(1.+sin(wgt0(:,:)))

  do j2=-jbuf, jbuf
    wgt0(:,j2) = wgt0(:,j2)*cos(lat(j+j2)*deg2rad)
  enddo
 
  wgt(1:nx/2,:) = wgt0(0:nx/2-1,:)
  wgt(nx/2+1:nx,:) = wgt0((nx+1)/2:1:-1,:)
  wgt(nx+1:nx*2,:) = wgt(1:nx,:)

  do j2=-jbuf, jbuf
  do n=1, nt
  do i=1, nx
    mean(i,j,n) = mean(i,j,n) +  &
        sum(var_in(1:nx,j+j2,n)*wgt(nx+2-i:nx*2+1-i,j2))
  enddo
  enddo
  enddo
  mean(:,j,:) = mean(:,j,:)/sum(wgt(1:nx,:))

  ENDDO  J_LAT

  ! 9-pt local smoothing (1-2-1 smoothing)
  N_SM:  DO is=1, 3

  prt = mean
  mean(2:nx-1,:,:) = 0.25*(2.*prt(2:nx-1,:,:) + prt(1:nx-2,:,:) + prt(3:nx,:,:))
  mean(1 ,:,:) = 0.25*(2.*prt(1 ,:,:) + prt(nx  ,:,:) + prt(2,:,:))
  mean(nx,:,:) = 0.25*(2.*prt(nx,:,:) + prt(nx-1,:,:) + prt(1,:,:))
  prt = mean
  do j=2+jbuf, ny-1-jbuf
    mean(:,j,:) = 0.25*(2.*prt(:,j,:) + prt(:,j-1,:) + prt(:,j+1,:))
  enddo

  ENDDO  N_SM

  ! calculate perturbation
  prt(:,:,:) = var_in(:,1+jbuf:ny-jbuf,:) - mean(:,:,:)

  deallocate( wgt, wgt0 )
  deallocate( coslat, sinlat, sin2dlonp2 )

! DUMP and FINALIZE

  call dump
  call finalize

  STOP


  CONTAINS


  SUBROUTINE init_read

  real ::  tmp1d(2)

  call get_var(file_i,'lon',tmp1d, count=(/2/))
  dlon = tmp1d(2) - tmp1d(1)
  call get_var(file_i,'lat',tmp1d, count=(/2/))
  dlat = tmp1d(2) - tmp1d(1)
  jbuf = int((d_crit+d_kernel)/r_earth/deg2rad/dlat) + 1

  i1 = ii(1)  ;  n1 = nnn(1)
  nx = ii(2) - i1 + 1  ;  nt = (nnn(2) - n1)/nnn(3) + 1
  j1 = jj(1) - jbuf
  ny = jj(2) + jbuf - j1 + 1

  allocate( lon(nx), lat(ny), t(nt) )

  call get_var(file_i,'lon',lon, start=(/i1/), count=(/nx/))
  call get_var(file_i,'lat',lat, start=(/j1/), count=(/ny/))
  call get_var(file_i,'t'  ,t  , start=(/n1/), count=(/nt/), stride=(/nnn(3)/))

  allocate( prt(nx,1+jbuf:ny-jbuf,nt), mean(nx,1+jbuf:ny-jbuf,nt) )
  allocate( var_in(nx,ny,nt) )

  call get_var(file_i,var_name,var_in, start=(/i1,j1,n1/), count=(/nx,ny,nt/), &
               stride=(/1,1,nnn(3)/), map=(/1,nx,nx*ny/))

  END subroutine init_read

  SUBROUTINE dump

  call put_var('overwrite',file_o(1),'lon',lon,axis='lon')
  call put_var('append'   ,file_o(1),'lat',lat(1+jbuf:ny-jbuf),axis='lat')
  call put_var('append'   ,file_o(1),'t'  ,t  ,axis='t'  )

  call put_var('append'   ,file_o(1),'prt_'//trim(var_name),prt, &
               axes=(/'lon','lat','t'/))

  call put_var('overwrite',file_o(2),'lon',lon,axis='lon')
  call put_var('append'   ,file_o(2),'lat',lat(1+jbuf:ny-jbuf),axis='lat')
  call put_var('append'   ,file_o(2),'t'  ,t  ,axis='t'  )

  call put_var('append'   ,file_o(2),'mean_'//trim(var_name),mean, &
               axes=(/'lon','lat','t'/))

  write(6,*)  ;  write(6,*) trim(file_o(1))
  write(6,*) trim(file_o(2))  ;  write(6,*)

  END subroutine dump

  SUBROUTINE finalize

  deallocate( prt, mean, var_in )
  deallocate( lon, lat, t )

  END subroutine finalize

END program PERTURBATION_EXTRACTION

