! for the operational output

PROGRAM analysis_UM

  use netio
  use um_anal
  use um_axis

  implicit none

  integer, parameter ::  dcase = 10010712
  integer, parameter ::  ftail1 = 048, ftail2 = 120, fitv = 24
  integer, parameter ::  nvar = 1
  integer, parameter ::  np = 21
  real,    parameter ::  missv = 1.e32
  real,    parameter ::  zh_int = 12.0
  real               ::  p(np)
  character(len=128) ::  ifdir, expname, vartype, ofdir
  character(len=64)  ::  ivfile(nvar), varname(nvar)

  ! pressure levels
  data  p / 1000., 925., 850., 700., 600., 500., 400., 300., 250., 200., &
            150., 100., 70., 50., 30., 20., 10., 5., 3., 2., 1. /
  ! files
  data                   ifdir   /'/data11/kyh/UM_case/'/
  data                   expname /'gwdc'/
  data                   vartype /'gwdc'/
  data                   ofdir   /'/data11/kyh/UM_case/'/ 
  ! var. names in input files
  data                   ivfile(1) /'fu'/
  ! input variables
  data                   varname(1) /'unspecified'/


  real, dimension(:,:,:), allocatable ::  var, exner
  real, dimension(:),     allocatable ::  lon, lat, zh, lonu, latv, zht
  real, dimension(:),     allocatable ::  time, t_fct, t
  real, dimension(np)                 ::  exner_o

  integer ::  nfct, fct
  integer ::  nx, ny, nz, nxu, nyv, nzt, nt, nt0, ifc
  integer ::  nxr, nyr, nzr, ntr, nd1a, nd2a, nd3a, nd4a
  integer ::  iz, it
  integer ::  i,j,k,n,iv, ncid, ncid2
  real    ::  temp
  character(len=32)  ::  c_axis(3,2)
  character(len=128) ::  fname, fnamep

  type ::  vset
    real,    dimension(:,:,:,:), allocatable ::  var_out
    real,    dimension(:),       allocatable ::  axis1, axis2, axis3, axis4
    integer, dimension(4)                    ::  nd
    character(len=32)                        ::  vname, axis(4)
  end type vset
  type(vset), dimension(1) ::  set


  exner_o(:) = (p(:)/1.e3)**kappa


  N_VAR:   DO iv=1, nvar
  N_TIM:   DO it=ftail1, ftail2, fitv


  write(fname,'(a,i8.8,a,i8.8,a,i3.3,a)')                  &
        trim(ifdir)//trim(expname)//'/20',dcase,'/'//      &
        trim(vartype)//'.'//trim(ivfile(iv))//'.20',dcase, &
        '+',it,'.nc'

  write(fnamep,'(a,i8.8,a,i8.8,a,i3.3,a)')                 &
        trim(ifdir)//trim(expname)//'/20',dcase,'/'//      &
        'std.p.20',dcase,'+',it,'.nc'

  call opennc(trim(fname),ncid)
  call opennc(trim(fnamep),ncid2)

  ! get dim. sizes and axis
  if (it == ftail1) then
    call diminfo(ncid,.TRUE.,zh_int, nx,ny,nz,nxu,nyv,nzt,c_axis)
    allocate( lonu(nxu), lat(ny), zh(nz) )
    allocate( lon(1), latv(1), zht(1) )
    call axisinfo(ncid,1,ny,nz,nxu,1,1,c_axis, lon,lat,zh,lonu,latv,zht)
    nx = nxu
  end if

  call timeinfo(ncid, nt)
  allocate( time(nt) )
  call get1d(ncid,'t',nt, time)

  ! get var
  nxr = NXU  ;   nd1a = NXU
  nyr = NY   ;   nd2a = NY
  nzr = NZ   ;   nd3a = NP
  ntr = NT   ;   nd4a = NT

  if (it == ftail1)  allocate( set(1)%var_out(nd1a,nd2a,nd3a,nd4a) )

  allocate( var(nxr,nyr,nzr), exner(nxr,nyr,nzr) )

  do n=1, nt

print*, 'read', time(n)
    call geta4d(ncid ,trim(varname(iv)),1,nxr,1,nyr,1,nzr,n,1, var  )
    call geta4d(ncid2,'p'              ,1,nx ,1,nyr,1,nzr,n,1, exner)
    exner(:,:,:) = (exner(:,:,:)/1.e5)**kappa

    do k=1, nz
    do j=1, ny
      temp = exner(1,j,k)
      do i=1, nx-1
        exner(i,j,k) = 0.5*(exner(i,j,k) + exner(i+1,j,k))
      enddo
      exner(nx,j,k) = 0.5*(exner(nx,j,k) + temp)
    enddo
    enddo

print*, 'interpolate'
    call vintp_p(nxu,ny,nz,exner,var,np,exner_o, set(1)%var_out(:,:,:,n))

  enddo

  deallocate( var, exner )

  if (it == ftail1) then
    set(1)%nd(1) = nd1a
    set(1)%nd(2) = nd2a
    set(1)%nd(3) = nd3a
    set(1)%axis = (/'longitude','latitude','p','t'/)
    set(1)%vname = varname(iv)
    allocate( set(1)%axis1(set(1)%nd(1)) )
    allocate( set(1)%axis2(set(1)%nd(2)) )
    allocate( set(1)%axis3(set(1)%nd(3)) )
    set(1)%axis1 = lonu
    set(1)%axis2 = lat
    set(1)%axis3 = p
  end if
  set(1)%nd(4) = nd4a
  allocate( set(1)%axis4(set(1)%nd(4)) )
  set(1)%axis4 = time

  call closenc(ncid)
  call closenc(ncid2)

  deallocate( time )

  ! dump
  write(fname,'(a,i8.8,a,i3.3,a)')                           &
        trim(ofdir)//trim(expname)//'/anal/'//               &
        trim(vartype)//'.'//trim(ivfile(iv))//'_p.20',dcase, &
        '+',it,'.nc'

  write(6,*)
  write(6,*) trim(fname)
  write(6,*)

  call outnc(trim(fname),nvar,set,'')


  ENDDO  N_TIM


  deallocate( lonu, lat, zh )
  deallocate( lon, latv, zht )

  deallocate( set(1)%var_out )


  ENDDO  N_VAR


  STOP

END program


SUBROUTINE vintp_p(nx,ny,nz,exner,var,np,exner_o,var_o)

  implicit none

  integer,                   intent(in) ::  nx, ny, nz, np
  real, dimension(np),       intent(in) ::  exner_o
  real, dimension(nx,ny,nz), intent(in) ::  var, exner

  real, dimension(nx,ny,np), intent(out) ::  var_o

  integer ::  i,j,k,kp
  integer ::  last, k0(nx,ny)


  last = 1


  N_LEV:   DO kp=1, np


  do j=1, ny
  do i=1, nx
    if (exner_o(kp) <= exner(i,j,last) ) then
      ! searching up
      k0(i,j) = -1  ! top level
      do k=last+1, nz
        if (exner_o(kp) > exner(i,j,k)) then
          k0(i,j) = k - 1
          EXIT
        end if
      enddo
    else
      ! searching down
      k0(i,j) = 0  ! bottom level
      do k=last-1, 1, -1
        if (exner_o(kp) <= exner(i,j,k)) then
          k0(i,j) = k
          EXIT
        end if
      enddo

    end if
    last = max(k0(i,j), 1)
  enddo
  enddo

  do j=1, ny
  do i=1, nx
    if (k0(i,j) == -1) then
      var_o(i,j,kp) = var(i,j,nz)
    else if (k0(i,j) == 0) then
      var_o(i,j,kp) = var(i,j,1)
    else
      var_o(i,j,kp) = ( (exner_o(kp)-exner(i,j,k0(i,j)))*           &
                        var(i,j,k0(i,j)+1) -                        &
                        (exner_o(kp)-exner(i,j,k0(i,j)+1))*         &
                        var(i,j,k0(i,j)) ) /                        &
                      ( exner(i,j,k0(i,j)+1) - exner(i,j,k0(i,j)) )
    end if
  enddo
  enddo


  ENDDO  N_LEV


  RETURN

END subroutine


