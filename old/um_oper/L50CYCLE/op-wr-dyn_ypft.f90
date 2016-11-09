PROGRAM W_RES_DYN

  use netio
  use um_anal
  use um_axis

  implicit none

  integer, parameter ::  nv = 5
  real,    parameter ::  h_scale = 7.e3
  real,    parameter ::  missv = 1.e32
  real,    parameter ::  zb_out = 0.e3
  integer, parameter ::  nvo = nv + 1
  character(len=64)  ::  ovarname(nvo)
  character(len=128) ::  f_namelist

  integer ::  yyyymm
  real    ::  lat_avg
  character(len=10)  ::  expname
  character(len=128) ::  file_i(2), file_o

  namelist /ANALCASE/ EXPNAME, YYYYMM
  namelist /PARAM/ LAT_AVG
  namelist /FILEIO/ FILE_I, FILE_O

  ! input variables
  data  f_namelist  /'./namelist/nl.input'/
  ! output variables
  data  ovarname(1  ) /'epd    '/
  data  ovarname(2  ) /'fu_ussp'/
  data  ovarname(3  ) /'fu_gwdo'/
  data  ovarname(4  ) /'fu_bldo'/
  data  ovarname(5  ) /'fu_gwdc'/
  data  ovarname(nvo) /'u_tend '/

  integer ::  nt, nfct, ny, nz, ilat_avg(2), nyo
  integer ::  nd1a, nd2a, nd3a, nd4a
  integer ::  i,j,k,n, iv, itim, ifct, it, ncid, ncid2, iz, ni
  logical ::  l_getdim, ex1, ex2
  character(len=32)  ::  c_axis(3,2)
  character(len=256) ::  fname, fname1, fname2

  real, dimension(:,:,:,:), allocatable ::  um, frc
  real, dimension(:,:),     allocatable ::  wmr
  real, dimension(:),       allocatable ::  lat, p, lato
  real, dimension(:),       allocatable ::  zp, t, t_fct

  type ::  vset
    real,    dimension(:,:,:,:), allocatable ::  var_out
    real,    dimension(:),       allocatable ::  axis1, axis2, axis3, axis4
    integer, dimension(4)                    ::  nd
    character(len=32)                        ::  vname, axis(4)
  end type vset
  type(vset), dimension(nvo) ::  set

! READ NAMELISTS

  open(10, file=trim(f_namelist), status='old')
  read(10, ANALCASE)  ;  read(10, PARAM)  ;  read(10, FILEIO)
  close(10)

  write(fname1,'(a,i6.6,a)') trim(file_i(1))//'/'// &
       trim(expname)//'.tem_ypft.',yyyymm,'.nc'
  write(fname2,'(a,i4.4,a,i2.2,a)') trim(file_i(1))//'/'// &
       trim(expname)//'.zgwd_ypft.',yyyymm/100,'.',mod(yyyymm,100),'.nc'
  call check_ex(fname1, ex1)  ;  call check_ex(fname2, ex2)
  if (.not. (ex1 .and. ex2))  STOP

  call opennc(fname1,ncid )
  call opennc(fname2,ncid2)

! DEFINE AXIS

  call getdim

  call get_um

! LOOP

  L_VAR:  DO iv=1, nvo


  ! allocate output
  nd1a = NY
  nd2a = NZ
  nd3a = NFCT
  nd4a = NT
  allocate( set(iv)%var_out(nd1a,nd2a,nd3a,nd4a) )

  ! get var
  allocate( frc(ny,nz,nfct,nt) )
  if (iv == 1) then
    call get4d(ncid,ovarname(iv),ny,nz,nfct,nt, frc)
  else if (iv == nvo) then
    call get4d(ncid,ovarname(iv),ny,nz,nfct,nt, frc)
    frc(:,:,:,:) = -frc(:,:,:,:)
  else
    call get4d(ncid2,ovarname(iv),ny,nz,nfct,nt, frc)
  end if
  frc(:,:,:,:) = frc(:,:,:,:)/86400.

  allocate( wmr(nyo,nz) )


  L_TIM:  DO itim=1, nt
  L_FCT:  DO ifct=1, nfct


  ! calculate forcing
  call wmr_dynbal_hydp_frc(ny,nz,1,lat,p,um(:,:,ifct,itim), &
           frc(:,:,ifct,itim),ilat_avg,nyo, lato,wmr)

  do k=1, nz
    set(iv)%var_out(1:ilat_avg(1) ,k,ifct,itim) = wmr(1:nyo/2    ,k)
    set(iv)%var_out(ilat_avg(2):ny,k,ifct,itim) = wmr(nyo/2+1:nyo,k)
    set(iv)%var_out(ilat_avg(1)+1:ilat_avg(2)-1,k,ifct,itim) = wmr(nyo/2,k)
  enddo


  ENDDO  L_FCT
  ENDDO  L_TIM


  deallocate( wmr )
  deallocate( frc )

  if (zb_out /= 0.) then
    iz = nd2a + 1
    do k=2, nd2a
      if (zp(k) > zb_out) then  ;  iz = k - 1  ;  EXIT  ;  end if
    enddo
    set(iv)%var_out(:,:iz-1,:,:) = missv
  end if


  set(iv)%vname = ovarname(iv)
  set(iv)%axis = (/'lat  ','zp ','fcst','t'/)
  set(iv)%nd(:) = (/nd1a,nd2a,nd3a,nd4a/)
  allocate( set(iv)%axis1(set(iv)%nd(1)) )
  allocate( set(iv)%axis2(set(iv)%nd(2)) )
  allocate( set(iv)%axis3(set(iv)%nd(3)) )
  allocate( set(iv)%axis4(set(iv)%nd(4)) )
  set(iv)%axis1 = lat
  set(iv)%axis2 = zp
  set(iv)%axis3 = t_fct
  set(iv)%axis4 = t


  ENDDO  L_VAR


  call closenc(ncid )
  call closenc(ncid2)

! DUMP

  write(6,*) '\n',trim(file_o),'\n'

  call outnc(trim(file_o),nvo,set,'W_RES_DYN')

! END

  deallocate( um )
  deallocate( t_fct, t )
  deallocate( lat, zp, p, lato )
  do iv=1, nvo
    deallocate( set(iv)%axis1, set(iv)%axis2, set(iv)%axis3, set(iv)%axis4 )
    deallocate( set(iv)%var_out )
  enddo

  STOP


  CONTAINS

  SUBROUTINE getdim

    implicit none

    real    ::  templ
    logical ::  l_err

    call dilen(ncid,'lat' , ny  ,l_err)
    call dilen(ncid,'zp'  , nz  ,l_err)
    call dilen(ncid,'fcst', nfct,l_err)
    call dilen(ncid,'t'   , nt  ,l_err)

    allocate( lat(ny), zp(nz), p(nz) )
    allocate( t_fct(nfct), t(nt) )
    call get1d(ncid,'lat' ,ny  , lat  )
    call get1d(ncid,'zp'  ,nz  , zp   )
    call get1d(ncid,'fcst',nfct, t_fct)
    call get1d(ncid,'t'   ,nt  , t    )
    p(:) = 1.e5*exp(-zp(:)/h_scale)

    templ = minval(abs(lat(:)+lat_avg))
    do j=1, ny
      if (abs(lat(j)+lat_avg) == templ)  ilat_avg(1) = j
    enddo
    ilat_avg(2) = ny + 1 - ilat_avg(1)
    nyo = ny - (ilat_avg(2)-ilat_avg(1)-1)
    allocate( lato(nyo) )

  END subroutine getdim

  SUBROUTINE get_um

    implicit none

    character(len=64)              ::  ivarname

    ! input variables
    data  ivarname /'u'/

    integer            ::  nx, tempy, tempz, nxu, nyv, nf, nf1, nf2, ncid3, ncid4
    integer            ::  nf_intg(nfct)
    real               ::  dayintg, t1(nfct), t2(nfct), fct_intg
    character(len=256) ::  fname3, fname4
    character(len=10)  ::  timec
    real, parameter    ::  v_small = 1.e-5

    real, dimension(:,:,:,:), allocatable ::  u
    real, dimension(:),       allocatable ::  fday


    allocate( um(ny,nz,nfct,nt) )

    fct_intg = t_fct(2) - t_fct(1)

    do itim=1, nt

      write(timec,'(i6.6,2i2.2)') yyyymm, int(t(itim)), &
            int(24.*(t(itim)-int(t(itim))))

      fname3 = trim(file_i(2))//'/'//timec//'/'// &
               trim(expname)//'.std_p.'//timec//'+024.nc'
      fname4 = trim(file_i(2))//'/'//timec//'/'// &
               trim(expname)//'.std_p.'//timec//'+120.nc'

      call check_ex(fname3, ex1)  ;  call check_ex(fname4, ex2)
      if (.not. (ex1 .and. ex2))  CYCLE

      if (itim == 1) then
        call opennc(fname3,ncid3)
        call opennc(fname4,ncid4)
        call diminfop(ncid3,.TRUE., nx,tempy,tempz,nxu,nyv,c_axis)
        call timeinfo(ncid3, nf1)
        call timeinfo(ncid4, nf2)
        nf = nf1 + nf2
        allocate( fday(nf) )
        call get1d(ncid3,'t',nf1, fday(1:nf1))
        call get1d(ncid4,'t',nf2, fday(nf1+1:nf))
        call closenc(ncid3)
        call closenc(ncid4)
        do ifct=1, nfct
          t1(ifct) = t_fct(ifct) - 0.5*fct_intg
          t2(ifct) = t_fct(ifct) + 0.5*fct_intg
          do n=1, nf
            if ( abs(fday(n)-t2(ifct)) < v_small )  nf_intg(ifct) = &
               int(fct_intg/(fday(n)-fday(n-1))+1)
          enddo
        enddo
      end if

      do ifct=1, nfct

        allocate( u(nx,ny,nz,nf_intg(ifct)) )

    NI:   DO ni=1, nf_intg(ifct)


    dayintg = t2(ifct) + fct_intg*real(ni-nf_intg(ifct))/(nf_intg(ifct)-1)

    do n=1, nf
      if ( abs(fday(n) - dayintg) < v_small )  it = n
    enddo

    fname = fname3

    if (dayintg > 1.0) then
      fname = fname4
      it = it - nf1
    end if

    call opennc(trim(fname),ncid3)
    call geta4d(ncid3,trim(ivarname),1,nxu,1,ny,1,nz,it,1, u(:,:,:,ni))
    call closenc(ncid3)


    ENDDO  NI


        do k=1, nz
        do j=1, ny
          um(j,k,ifct,itim) = sum( u(:,j,k,:) )/float(nxu*nf_intg(ifct))
        enddo
        enddo

        deallocate( u )

      enddo  ! ifct

    enddo  ! itim

    deallocate( fday )

  END subroutine get_um

END program W_RES_DYN


SUBROUTINE check_ex(fname,existence)
  
  implicit none
  
  character(len=*), intent(in) ::  fname

  logical ::  existence
  
  
  inquire(file=trim(fname), exist=existence)
  if (.not. existence)  print*, '    ',trim(fname),' not found. - passed'
  
  RETURN
  
END subroutine check_ex


SUBROUTINE wmr_dynbal_hydp_frc(ny,nz,nt,lat,p,um,frc,ilat_avg,nyo, &
                               lato,wmr_d)

! following Randel et al. (2002)

  use avg
  use deriv
  use integ
  use utiletc

  implicit none

  logical               ::  uniform_y = .TRUE., uniform_z = .FALSE.
  logical, dimension(2) ::  l_ybdy = (/.FALSE.,.FALSE./), &
                            l_zbdy = (/.FALSE.,.FALSE./)

  integer,                   intent(in) ::  ny, nz, nt
  integer,                   intent(in) ::  ilat_avg(2), nyo
  real, dimension(ny),       intent(in) ::  lat
  real, dimension(nz),       intent(in) ::  p
  real, dimension(ny,nz,nt), intent(in) ::  um
  real, dimension(ny,nz,nt), intent(in) ::  frc

  real, dimension(nyo),    intent(out) ::  lato
  real, dimension(nyo,nz), intent(out) ::  wmr_d

  integer                   ::  j,k,n,m, k1, kn, ki, j1, j2
  real                      ::  ztop, ypole1, ypole2
  real, dimension(ny)       ::  y, f, cosphi, rsinphi
  real, dimension(nz)       ::  zp, rho0s, rrdsinphi
  real, dimension(ny)       ::  trid_a, trid_b, trid_c, trid_d
  real, dimension(nz)       ::  tempz1, tempz2
  real, dimension(0:ny+1)   ::  fn_str2
  real, dimension(ny,nz)    ::  rho0, r_fac
  real, dimension(nyo,nz)   ::  temp2d
  real, dimension(ny,nz,nt) ::  f_hat, dumdzdfh, temp3d
  real, dimension(ny,nz,nt) ::  fn_str
  real, dimension(3,ny)     ::  dcoef
  real*8                    ::  dcoef8(3)
  logical                   ::  l_zrev

  include 'c_math.inc'
  include 'c_phys.inc'


  j1 = ilat_avg(1)  ;  j2 = ilat_avg(2)

  if (p(1) < p(2))  l_zrev = .TRUE.

! latr, zp, cosphir, rsinphir, rho0r, rho0o, f_hatr

  do k=1, nz
    zp(k) = -h_scale*log(p(k)/p0)
  enddo
  ztop = maxval(zp)

  y(:) = r_earth*lat(:)*deg2rad
  ypole1 = r_earth*90.*deg2rad
  ypole1 = sign(ypole1,y(1))
  ypole2 = -ypole1

  call coslat(ny,lat, cosphi)
  rsinphi(:) = r_earth*sin(lat(:)*deg2rad)

  rho0s(:) = rho_s*exp(-zp(:)/h_scale)

  ! for hydrostatic eqn.
  do k=1, nz
    rho0(:,k) = rho0s(k)
  enddo
  r_fac(:,:) = 1.

!  ! for rho0(y,z) in the nonhydrostatic eqn.
!  do k=1, nz
!    r_fac(:,k) = rho0(:,k) / rho0s(k)
!  enddo

  rrdsinphi(:) = 0.
  do j=j1, j2-1
    rrdsinphi(:) = rrdsinphi(:) + 0.5*(r_fac(j,:)+r_fac(j+1,:))* &
                   (rsinphi(j+1)-rsinphi(j))
  enddo


  call f_coriolis(ny,lat, f)
  call div_lat(1,ny,nz,nt,uniform_y,lat,um,l_ybdy,0., temp3d)
  do n=1, nt
  do k=1, nz
    f_hat(:,k,n) = f(:) - temp3d(:,k,n)
  enddo
  enddo

  call deriv1d((/ny,nz,nt,1/),2,uniform_z,zp,um,l_zbdy,0., dumdzdfh)
  dumdzdfh(:,:,:) = dumdzdfh(:,:,:) / f_hat(:,:,:)

  dcoef(:,:) = 0.
  do j=2, ny-1
    call fdcoef(2,3,dble(y(j)),dble(y(j-1:j+1)), dcoef8)
    dcoef(:,j) = dcoef8(:)
  enddo
  if (abs(lat(1)) /= 90.) then
    call fdcoef(2,3,dble(y(1)),dble((/ypole1,y(1),y(2)/)), dcoef8)
    dcoef(:,1) = dcoef8(:)
  end if
  if (abs(lat(ny)) /= 90.) then
    call fdcoef(2,3,dble(y(ny)),dble((/y(ny-1),y(ny),ypole2/)), dcoef8)
    dcoef(:,ny) = dcoef8(:)
  end if

  dcoef(1,j1) = -1./(y(j1)-y(j1-1))
  dcoef(2,j1) =  1./(y(j1)-y(j1-1))
  dcoef(3,j1) =  0.

  dcoef(1,j2) =  0.
  dcoef(2,j2) = -1./(y(j2+1)-y(j2))
  dcoef(3,j2) =  1./(y(j2+1)-y(j2))

dcoef = 0.

  if ( l_zrev ) then
    k1 = nz
    kn = 2
    ki = 1
  else
    k1 = 1
    kn = nz-1
    ki = -1
  end if

  do k=kn, k1, ki
    tempz1(k) = 2./(zp(k)-zp(k-ki)) - 1./h_scale
    tempz2(k) = 2./(zp(k)-zp(k-ki)) + 1./h_scale
  enddo


! for each forcing

  do n=1, nt
    temp3d(:,:,n) = frc(:,:,n)*r_fac(:,:)/f_hat(:,:,n)
  enddo

  fn_str(:,:,:) = 0.
  fn_str2(0) = 0.  ;  fn_str2(ny+1) = 0.

  do n=1, nt
  do k=kn, k1, ki

    trid_a(:) = dcoef(1,:)*dumdzdfh(:,k,n)
    trid_b(:) = dcoef(2,:)*dumdzdfh(:,k,n) + tempz1(k)
    trid_c(:) = dcoef(3,:)*dumdzdfh(:,k,n)

    fn_str2(1:ny) = fn_str(:,k-ki,n)
    do j=1, ny
      trid_d(j) = -sum( dcoef(:,j)*fn_str2(j-1:j+1) )* &
                  dumdzdfh(j,k-k1,n)
    enddo
    trid_d(:) = trid_d(:) + tempz2(k)*fn_str(:,k-ki,n) + &
                cosphi(:)*(temp3d(:,k-ki,n)+temp3d(:,k,n))

    ! boundary condition
    trid_a(1 ) = 0.
    trid_c(ny) = 0.

    ! solve
    call tridag(j1,trid_a(1:j1),trid_b(1:j1),trid_c(1:j1), &
                trid_d(1:j1))
    fn_str(1:j1,k,n) = trid_d(1:j1)

    call tridag(ny-j2+1,trid_a(j2:ny),trid_b(j2:ny),trid_c(j2:ny), &
                trid_d(j2:ny))
    fn_str(j2:ny,k,n) = trid_d(j2:ny)

  enddo
  enddo

  call avg_d3((/j1    ,nz,nt,1/),fn_str(1:j1 ,:,:),1., temp2d(1:j1    ,:))
  call avg_d3((/nyo-j1,nz,nt,1/),fn_str(j2:ny,:,:),1., temp2d(j1+1:nyo,:))

  call deriv1d((/j1,nz,1,1/),1,.False.,rsinphi(1:j1), &
               temp2d(1:j1,:),l_ybdy,0., wmr_d(1:j1,:))

  call deriv1d((/nyo-j1,nz,1,1/),1,.False.,rsinphi(j2:ny), &
               temp2d(j1+1:nyo,:),l_ybdy,0., wmr_d(j1+1:nyo,:))

  wmr_d(1:j1-1  ,:) = wmr_d(1:j1-1  ,:) / r_fac(1:j1-1 ,:)
  wmr_d(j1+2:nyo,:) = wmr_d(j1+2:nyo,:) / r_fac(j2+1:ny,:)

  wmr_d(j1  ,:) = (temp2d(j1+1,:)-temp2d(j1,:))/rrdsinphi(:)
  wmr_d(j1+1,:) = wmr_d(j1,:)

  lato(1   :j1 ) = lat(1 :j1)
  lato(j1+1:nyo) = lat(j2:ny)

  RETURN

END subroutine wmr_dynbal_hydp_frc

