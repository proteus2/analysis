PROGRAM PZ07_WRF_Z

  use pz07
  use easync

  implicit none

  integer            ::  nnn(3), jj(2), ii(2), n_smooth9
  real               ::  dz_shear, dt_tend
  character(len=16)  ::  var_name(5)
  character(len=256) ::  f_namelist, file_i(5), file_v_i(2), file_o

  namelist /PARAM/ VAR_NAME, NNN, JJ, II, DZ_SHEAR, DT_TEND, N_SMOOTH9
  namelist /FILEIO/ FILE_I, FILE_V_I, FILE_O

  integer, parameter ::  nv = 3
  character(len=16), dimension(nv) ::  var_o =                           &
      (/'nbe','dpt_dt','dav_dt'/)
 
  integer ::  i1, j1, n1, nx, ny, nt
  integer ::  iv, n, ism

  real, dimension(:,:,:), allocatable ::  u, v, p, pt, inv_rho, w, pts,  &
                                          ut, vt, ptt, varo
  real, dimension(:,:),   allocatable ::  nbe, dptdt, davdt
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

! CALCULATE DIAGNOSTICS

  call nbe_z_pz07( nx,ny,lat,u,v,p,inv_rho,0.,n_smooth9, nbe )

  call dptdt_z_pz07( nx,ny,lat,pt,u,v,w,pts,ptt,0., dptdt )

  call davdt_z_pz07( nx,ny,lat,u,v,ut,vt,0., davdt )

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

  allocate( u(nx,ny,1), v(nx,ny,1), p(nx,ny,1), pt(nx,ny,1), w(nx,ny,1) )
  allocate( inv_rho(nx,ny,1) )
  allocate( pts(nx,ny,1), ut(nx,ny,1), vt(nx,ny,1), ptt(nx,ny,1) )

  allocate( nbe(nx,ny), dptdt(nx,ny), davdt(nx,ny) )
  allocate( varo(nx,ny,nv) )

  ! dump coordinates
  call put_var('overwrite',file_o,'lon',lon,axis='lon')
  call put_var('append'   ,file_o,'lat',lat,axis='lat')

  END subroutine init

  SUBROUTINE read_n

  integer ::  nn
  real, dimension(nx,ny,1) ::  tmp
  real, dimension(nx,ny,2) ::  tmp2

  nn = n1 + (n-1)*nnn(3)

  call get_var(file_i(1),var_name(1),u, start=(/i1,j1,nn/),count=(/nx,ny,1/), &
               stride=(/1,1,nnn(3)/))

  call get_var(file_i(2),var_name(2),v, start=(/i1,j1,nn/),count=(/nx,ny,1/), &
               stride=(/1,1,nnn(3)/))

  call get_var(file_i(3),var_name(3),p, start=(/i1,j1,nn/),count=(/nx,ny,1/), &
               stride=(/1,1,nnn(3)/))

  call get_var(file_i(4),var_name(4),pt, start=(/i1,j1,nn/),count=(/nx,ny,1/), &
               stride=(/1,1,nnn(3)/))

  call get_var(file_i(5),var_name(5),w, start=(/i1,j1,nn/),count=(/nx,ny,1/), &
               stride=(/1,1,nnn(3)/))

  pt(:,:,:) = pt(:,:,:) + 300.  ! pt := +300 for WRF

  inv_rho(:,:,:) = rdrp0*(pt(:,:,:))*(p0/p(:,:,:))**cvrcp

  call get_var(file_v_i(1),var_name(4),tmp, start=(/i1,j1,nn/),count=(/nx,ny,1/), &
               stride=(/1,1,nnn(3)/))
  call get_var(file_v_i(2),var_name(4),pts, start=(/i1,j1,nn/),count=(/nx,ny,1/), &
               stride=(/1,1,nnn(3)/))

  pts(:,:,:) = (pts(:,:,:) - tmp(:,:,:))/dz_shear

  call get_var(file_i(1),var_name(1),tmp2, start=(/i1,j1,nn-1/),count=(/nx,ny,2/), &
               stride=(/1,1,2/))

  ut(:,:,1) = (tmp2(:,:,2) - tmp2(:,:,1))/dt_tend

  call get_var(file_i(2),var_name(2),tmp2, start=(/i1,j1,nn-1/),count=(/nx,ny,2/), &
               stride=(/1,1,2/))

  vt(:,:,1) = (tmp2(:,:,2) - tmp2(:,:,1))/dt_tend

  call get_var(file_i(4),var_name(4),tmp2, start=(/i1,j1,nn-1/),count=(/nx,ny,2/), &
               stride=(/1,1,2/))

  ptt(:,:,1) = (tmp2(:,:,2) - tmp2(:,:,1))/dt_tend

  END subroutine read_n

  SUBROUTINE dump_n

  varo(:,:,1) = nbe
  varo(:,:,2) = dptdt
  varo(:,:,3) = davdt

  call put_var('append',file_o,'t',(/t(n),t(n)/), is_record=.TRUE.,      &
               axis='t', start=(/n/), count=(/1/))
  do iv=1, nv
    call put_var('append',file_o,var_o(iv),                              &
                 reshape(varo(:,:,iv),(/nx,ny,1/)), is_record=.TRUE.,    &
                 axes=(/'lon','lat','t'/), start=(/1,1,n/))
  enddo
 
  END subroutine dump_n

  SUBROUTINE finalize

  write(6,*)
  write(6,*) trim(file_o)
  write(6,*)

  deallocate( u, v, p, inv_rho, w )
  deallocate( pts, ut, vt, ptt )
  deallocate( lon, lat, t )

  END subroutine finalize

END program PZ07_WRF_Z

