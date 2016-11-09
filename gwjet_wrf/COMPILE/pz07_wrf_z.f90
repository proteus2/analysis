PROGRAM PZ07_WRF_Z

  use pz07
  use easync

  implicit none

  integer            ::  nnn(3), nnn0(3), jj(2), ii(2), n_smooth9
  real               ::  dz_shear, dt_tend
  character(len=16)  ::  var_name(3), var_name0(3)
  character(len=256) ::  f_namelist, file_i, file_v_i(2), file_i0(3),    &
                         file_o

  namelist /PARAM/ VAR_NAME, VAR_NAME0, NNN, NNN0, JJ, II, DZ_SHEAR,     &
                   DT_TEND, N_SMOOTH9
  namelist /FILEIO/ FILE_I, FILE_V_I, FILE_I0, FILE_O

  integer, parameter ::  nv = 6
  character(len=64), dimension(nv) ::  var_o
 
  integer ::  i1, j1, n1, n1_0, nx, ny, nt, nt0
  integer ::  iv, n, ism
  real    ::  pt0

  real, dimension(:,:,:), allocatable ::  nbes, nbest, dptdt, davdts,    &
                                          u, v, varo
  real, dimension(:,:),   allocatable ::  f_d, f_t, f_v, frc, f_d1, frc1
  real, dimension(:),     allocatable ::  lon, lat, t

  real, parameter ::  p0 = 1000.e2
  real, parameter ::  cp = 1004.
  real, parameter ::  rd = 287.
  real, parameter ::  rdrp0 = rd/p0
  real, parameter ::  cvrcp = 1. - rd/cp

  var_o(:) = varname_pz07_z(:)

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

  call total_z_pz07( nx,ny,lat,dptdt,nbes,davdts,nbest, &
                     pt0,u,v,0.,n_smooth9,              &
                     f_d,f_t,f_v,frc,f_d1,frc1 )

! DUMP

  call dump_n

  ENDDO  N_TIME

! FINALIZE

  call finalize

  STOP


  CONTAINS


  SUBROUTINE init

  i1 = ii(1)  ;  j1 = jj(1)  ;  n1 = nnn(1)  ;  n1_0 = nnn0(1)
  nx = ii(2) - i1 + 1  ;  ny = jj(2) - j1 + 1  ;  nt = (nnn(2) - n1)/nnn(3) + 1
  nt0 = (nnn0(2) - n1_0)/nnn0(3) + 1
  if (nt /= nt0) then
    print*, 'Check NNN0'  ;  STOP
  end if

  allocate( lon(nx), lat(ny), t(nt) )

  call get_var(file_i,'lon',lon, start=(/i1/), count=(/nx/))
  call get_var(file_i,'lat',lat, start=(/j1/), count=(/ny/))
  call get_var(file_i,'t'  ,t  , start=(/n1/), count=(/nt/), stride=(/nnn(3)/))

  allocate( nbes(nx,ny,1), nbest(nx,ny,1), dptdt(nx,ny,1), davdts(nx,ny,1) )
  allocate( u(nx,ny,1), v(nx,ny,1) )

  allocate( f_d(nx,ny), f_t(nx,ny), f_v(nx,ny), frc(nx,ny) )
  allocate( f_d1(nx,ny), frc1(nx,ny) )
  allocate( varo(nx,ny,nv) )

  ! dump coordinates
  call put_var('overwrite',file_o,'lon',lon,axis='lon')
  call put_var('append'   ,file_o,'lat',lat,axis='lat')

  END subroutine init

  SUBROUTINE read_n

  integer ::  nn, nn0
  real, dimension(nx,ny,1) ::  tmp
  real, dimension(nx,ny,2) ::  tmp2

  nn = n1 + (n-1)*nnn(3)
  nn0 = n1_0 + (n-1)*nnn0(3)

  call get_var(file_v_i(1),var_name(1),tmp, start=(/i1,j1,nn/),count=(/nx,ny,1/), &
               stride=(/1,1,nnn(3)/))

  call get_var(file_v_i(2),var_name(1),nbes, start=(/i1,j1,nn/),count=(/nx,ny,1/), &
               stride=(/1,1,nnn(3)/))

  nbes(:,:,:) = (nbes(:,:,:) - tmp(:,:,:))/dz_shear

  call get_var(file_v_i(1),var_name(3),tmp, start=(/i1,j1,nn/),count=(/nx,ny,1/), &
               stride=(/1,1,nnn(3)/))

  call get_var(file_v_i(2),var_name(3),davdts, start=(/i1,j1,nn/),count=(/nx,ny,1/), &
               stride=(/1,1,nnn(3)/))

  davdts(:,:,:) = (davdts(:,:,:) - tmp(:,:,:))/dz_shear

  call get_var(file_i,var_name(2),dptdt, start=(/i1,j1,nn/),count=(/nx,ny,1/), &
               stride=(/1,1,nnn(3)/))


  call get_var(file_i0(1),var_name0(1),u, start=(/i1,j1,nn0/),count=(/nx,ny,1/), &
               stride=(/1,1,nnn0(3)/))

  call get_var(file_i0(2),var_name0(2),v, start=(/i1,j1,nn0/),count=(/nx,ny,1/), &
               stride=(/1,1,nnn0(3)/))

  call get_var(file_i0(3),var_name0(3),tmp, start=(/i1,j1,nn0/),count=(/nx,ny,1/), &
               stride=(/1,1,nnn0(3)/))

  pt0 = sum(tmp(:,ny/2+1,1))/float(nx) + 300.  ! pt0 := +300 for WRF


  call get_var(file_v_i(1),var_name(1),tmp2, start=(/i1,j1,nn-1/),count=(/nx,ny,2/), &
               stride=(/1,1,2/))

  tmp(:,:,1) = (tmp2(:,:,2) - tmp2(:,:,1))/dt_tend

  call get_var(file_v_i(2),var_name(1),tmp2, start=(/i1,j1,nn-1/),count=(/nx,ny,2/), &
               stride=(/1,1,2/))

  nbest(:,:,1) = ( (tmp2(:,:,2) - tmp2(:,:,1))/dt_tend - tmp(:,:,1) )/dz_shear

  END subroutine read_n

  SUBROUTINE dump_n

  varo(:,:,1) = f_d
  varo(:,:,2) = f_t
  varo(:,:,3) = f_v
  varo(:,:,4) = frc
  varo(:,:,5) = f_d1
  varo(:,:,6) = frc1

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

  deallocate( nbes, nbest, dptdt, davdts )
  deallocate( u, v )
  deallocate( lon, lat, t )

  END subroutine finalize

END program PZ07_WRF_Z

