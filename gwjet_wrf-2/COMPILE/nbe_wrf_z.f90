PROGRAM NBE_WRF_Z

  use nbe
  use easync

  implicit none

  integer            ::  nnn(3), jj(2), ii(2), n_smooth9
  real               ::  dz_shear
  character(len=16)  ::  var_name(5)
  character(len=256) ::  f_namelist, file_i(5), file_v_i(4), file_o_d,   &
                         file_o_t

  namelist /PARAM/ VAR_NAME, NNN, JJ, II, DZ_SHEAR, N_SMOOTH9
  namelist /FILEIO/ FILE_I, FILE_V_I, FILE_O_D, FILE_O_T

  integer, parameter ::  nv = 7
  character(len=16), dimension(nv) ::  var_o =                           &
      (/'d_nbe_maj','curv','nl_min','adv','ageo','fvor','div'/)
  character(len=16), dimension(nv) ::  file_o_h =                        &
      (/'d_nbe-maj','d_nbe-curv','d_nbe-nl_min','adv_div','d_nbe-ageo',  &
        'd_nbe-fvor','div'/)
 
  character(len=256), dimension(nv) ::  file_o
  integer ::  i1, j1, n1, nx, ny, nt
  integer ::  iv, n, ism

  real, dimension(:,:,:), allocatable ::  u, v, p, inv_rho, w, us, vs,   &
                                          varo
  real, dimension(:,:),   allocatable ::  div, d_nbe, adv, nbe_maj,      &
                                          curv, nl_min, ageo, two_jac,   &
                                          fvor
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

! CALCULATE D_NBE

  call nbe_z( nx,ny,lat,u,v,p,inv_rho,w,us,vs,0.,n_smooth9, &
              div,d_nbe,adv,nbe_maj,curv,nl_min,ageo,two_jac,fvor )

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

  allocate( u(nx,ny,1), v(nx,ny,1), p(nx,ny,1), inv_rho(nx,ny,1), w(nx,ny,1) )
  allocate( us(nx,ny,1), vs(nx,ny,1) )

  allocate( div(nx,ny), d_nbe(nx,ny), adv(nx,ny) )
  allocate( nbe_maj(nx,ny), curv(nx,ny), nl_min(nx,ny) )
  allocate( ageo(nx,ny), two_jac(nx,ny), fvor(nx,ny) )
  allocate( varo(nx,ny,nv) )

  ! dump coordinates
  do iv=1, nv
    file_o(iv) = trim(file_o_d)//'/'//trim(file_o_h(iv))//trim(file_o_t)
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

  call get_var(file_i(3),var_name(3),p, start=(/i1,j1,nn/),count=(/nx,ny,1/), &
               stride=(/1,1,nnn(3)/))

  call get_var(file_i(4),var_name(4),inv_rho, start=(/i1,j1,nn/),count=(/nx,ny,1/), &
               stride=(/1,1,nnn(3)/))

  call get_var(file_i(5),var_name(5),w, start=(/i1,j1,nn/),count=(/nx,ny,1/), &
               stride=(/1,1,nnn(3)/))

  inv_rho(:,:,:) = rdrp0*(inv_rho(:,:,:)+300.)*(p0/p(:,:,:))**cvrcp  ! pt := +300 for WRF

  call get_var(file_v_i(1),var_name(1),tmp, start=(/i1,j1,nn/),count=(/nx,ny,1/), &
               stride=(/1,1,nnn(3)/))
  call get_var(file_v_i(2),var_name(1),us,  start=(/i1,j1,nn/),count=(/nx,ny,1/), &
               stride=(/1,1,nnn(3)/))

  us(:,:,:) = (us(:,:,:) - tmp(:,:,:))/dz_shear

  call get_var(file_v_i(3),var_name(2),tmp, start=(/i1,j1,nn/),count=(/nx,ny,1/), &
               stride=(/1,1,nnn(3)/))
  call get_var(file_v_i(4),var_name(2),vs,  start=(/i1,j1,nn/),count=(/nx,ny,1/), &
               stride=(/1,1,nnn(3)/))

  vs(:,:,:) = (vs(:,:,:) - tmp(:,:,:))/dz_shear

  END subroutine read_n

  SUBROUTINE dump_n

  varo(:,:,1) = nbe_maj
  varo(:,:,2) = curv
  varo(:,:,3) = nl_min
  varo(:,:,4) = adv
  varo(:,:,5) = ageo
  varo(:,:,6) = fvor
  varo(:,:,7) = div

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
  do iv=1, 7
    write(6,*) trim(file_o(iv))
  enddo
  write(6,*)

  deallocate( u, v, p, inv_rho, w )
  deallocate( us, vs )
  deallocate( lon, lat, t )

  END subroutine finalize

END program NBE_WRF_Z

