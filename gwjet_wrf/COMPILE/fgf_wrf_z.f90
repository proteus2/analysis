PROGRAM FGF_WRF_Z

  use fgf
  use easync

  implicit none

  integer            ::  nnn(3), jj(2), ii(2)
  character(len=16)  ::  var_name(3)
  character(len=256) ::  f_namelist, file_i(3), file_o_d, file_o_t

  namelist /PARAM/ VAR_NAME, NNN, JJ, II
  namelist /FILEIO/ FILE_I, FILE_O_D, FILE_O_T

  integer, parameter ::  nv = 6
  character(len=64),  dimension(nv) ::  var_o
  character(len=256), dimension(nv) ::  file_o
  integer ::  i1, j1, n1, nx, ny, nt
  integer ::  iv, n, ism

  real, dimension(:,:,:), allocatable ::  u, v, pt, varo
  real, dimension(:,:),   allocatable ::  pt_grad, ff, div, def,         &
                                          ff_div, ff_def
  real, dimension(:),     allocatable ::  lon, lat, t

  real, parameter ::  p0 = 1000.e2
  real, parameter ::  cp = 1004.
  real, parameter ::  rd = 287.
  real, parameter ::  rdrp0 = rd/p0
  real, parameter ::  cvrcp = 1. - rd/cp

! READ NAMELISTS

  call getarg(1,f_namelist)
  open(10, file=trim(f_namelist), status='old')
  read(10, PARAM)  ;  read(10, FILEIO)
  close(10)

! INITIALIZE ARRAYS

  call init

  N_TIME:  DO n=1, nt

! READ DATA

  call read_n

! CALCULATE FRONTOGENESIS FUNCTION

  call fgf_2d( nx,ny,lat,u,v,pt,0., &
               pt_grad,ff,div,def,ff_div,ff_def )

! DUMP

  call dump_n

  ENDDO  N_TIME

! FINALIZE

  call finalize

  STOP


  CONTAINS


  SUBROUTINE init

  i1 = ii(1)  ;  j1 = jj(1)  ;  n1 = nnn(1)
  nx = ii(2) - i1 + 1  ;  ny = jj(2) - j1 + 1  ;  nt = (nnn(2) - n1)/nnn(3) + 1

  allocate( lon(nx), lat(ny), t(nt) )

  call get_var(file_i(1),'lon',lon, start=(/i1/), count=(/nx/))
  call get_var(file_i(1),'lat',lat, start=(/j1/), count=(/ny/))
  call get_var(file_i(1),'t'  ,t  , start=(/n1/), count=(/nt/), stride=(/nnn(3)/))

  allocate( u(nx,ny,1), v(nx,ny,1), pt(nx,ny,1) )

  allocate( pt_grad(nx,ny), ff(nx,ny), div(nx,ny), def(nx,ny) )
  allocate( ff_div(nx,ny), ff_def(nx,ny) )
  allocate( varo(nx,ny,nv) )

  var_o(1:6) = varname_fgf_2d(:)

  ! dump coordinates
  do iv=1, nv
    file_o(iv) = trim(file_o_d)//'/'//trim(var_o(iv))//trim(file_o_t)
    call put_var('overwrite',file_o(iv),'lon',lon,axis='lon')
    call put_var('append'   ,file_o(iv),'lat',lat,axis='lat')
  enddo

  END subroutine init

  SUBROUTINE read_n

  integer ::  nn
  real, dimension(nx,ny,1) ::  tmp

  nn = n1 + (n-1)*nnn(3)

  call get_var(file_i(1),var_name(1),u, start=(/i1,j1,nn/),count=(/nx,ny,1/), &
               stride=(/1,1,nnn(3)/))

  call get_var(file_i(2),var_name(2),v, start=(/i1,j1,nn/),count=(/nx,ny,1/), &
               stride=(/1,1,nnn(3)/))

  call get_var(file_i(3),var_name(3),pt, start=(/i1,j1,nn/),count=(/nx,ny,1/), &
               stride=(/1,1,nnn(3)/))
  pt(:,:,:) = pt(:,:,:) + 300.

  END subroutine read_n

  SUBROUTINE dump_n

  varo(:,:,1) = pt_grad
  varo(:,:,2) = ff
  varo(:,:,3) = div
  varo(:,:,4) = def
  varo(:,:,5) = ff_div
  varo(:,:,6) = ff_def

  do iv=1, nv
    call put_var('append',file_o(iv),'t',(/t(n),t(n)/), is_record=.TRUE., &
                 axis='t', start=(/n/), count=(/1/))
    call put_var('append',file_o(iv),var_o(iv),                       &
                 reshape(varo(:,:,iv),(/nx,ny,1/)), is_record=.TRUE., &
                 axes=(/'lon','lat','t'/), start=(/1,1,n/))
  enddo
 
  END subroutine dump_n

  SUBROUTINE finalize

  write(6,*)
  do iv=1, nv
    write(6,*) trim(file_o(iv))
  enddo
  write(6,*)

  deallocate( u, v, pt )
  deallocate( lon, lat, t )

  END subroutine finalize

END program FGF_WRF_Z

