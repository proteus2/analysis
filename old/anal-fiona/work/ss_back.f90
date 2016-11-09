program ko_spectrum_S

  use pwrspd
  use netcdfio

  implicit none

  integer, parameter :: nx =1202, nz = 102
  integer, parameter :: nzp0 = 350
  real,    parameter :: dx = 1.e3, dt = 60., dz = 300.
  real,    parameter :: h_scale = 6.3e3, zref = 60.e3 
  real,    parameter :: ptbar1 = 240.5123  !300.2588

  integer, parameter :: nst =  61,    nnd = 480
  real,    parameter :: xst =   0.e3, xnd =1200.e3
  real,    parameter :: zst =  60.e3, znd =  90.e3, zint =  5.e3
  real,    parameter :: x0 =    0.e3, z0 = 60.e3

  integer, parameter :: eval = 0
  integer, parameter :: nt = nnd-nst+1
  integer, parameter :: noo = nt/2*2+1
  integer :: i,k,n, ncid, tag1, tag2, kre, kfind
  integer :: nxp, ist, ind, nkk, nzre, tmpi
  real, dimension(nx,nz) :: fdata
  real                   :: x(nx), z(nz), t(nt), zfind
  real                   :: ubardat(nzp0), n2bardat(nzp0), n2int(nzp0)
  real                   :: oo(noo), oo2(nt)
  character*6   :: time
  character*128 :: fname

  integer, dimension(:),  allocatable :: kst, knd
  real, dimension(:),     allocatable :: zre, kint, kk, kk2
  real, dimension(:),     allocatable :: ubar, n2bar, ptbar, dudz, d2udz
  real, dimension(:,:),   allocatable :: fp, cc
  real, dimension(:,:,:), allocatable :: f_tot, d_ave, mch2
  double complex, dimension(:,:), allocatable :: coef
  complex, dimension(:,:,:), allocatable :: ss, s_tot, ss2


  ! domain (input)
  do i=1, nx
    x(i) = x0 + (i-1.5)*dx
  enddo
  do k=1, nz
    z(k) = z0 + (k-1.5)*dz
  enddo
  do n=1, nt
    t(n) = (nst-1+n)*dt
  enddo

  ! analysis domain (x) 
  tag1 = 0   ;   tag2 = 0
  do i=1, nx
    if (x(i) .ge. xst)  tag1 = tag1 + 1
    if (x(i) .gt. xnd)  tag2 = tag2 + 1
    if (tag1 .eq. 1)  ist = i
    if (tag2 .eq. 1)  ind = i - 1
  enddo
  nxp = ind - ist + 1
  nkk = nxp/2 + 1
  allocate(kk(nkk),cc(nkk,noo),kk2(nxp))

  ! analysis domain (z)
  if ((znd-zst)/zint-int((znd-zst)/zint) .ne. 0) then
    print *, '( ( zt - zb ) / interval ) must be an integer.'
    STOP
  end if

  nzre = int((znd-zst)/zint)
  allocate(zre(nzre), kst(nzre), knd(nzre), kint(nzre))
  zfind = zst   ;   kfind = 1
  do kre=1, nzre
    zre(kre) = zst + (kre-0.5) * zint
    tag1 = 0
    do k=kfind, nz
      if (z(k) .ge. zfind)  tag1 = tag1 + 1
      if (tag1 .eq. 1)  kst(kre) = k
    enddo
    zfind = zfind + zint
    kfind = kst(kre) + 1
  enddo
  do kre=1, nzre-1
    knd(kre) = kst(kre+1) - 1
  enddo
  tag2 = 0
  do k=kst(nzre), nz
    if (z(k) .ge. znd)  tag2 = tag2 + 1
    if (tag2 .eq. 1)  knd(nzre) = k - 1
  enddo
  kint(:) = real(knd(:) - kst(:) + 1)

  ! ubar, n2bar, ptbar
  allocate(ubar(kst(1):knd(nzre)),n2bar(kst(1):knd(nzre)),ptbar(kst(1):knd(nzre)))
  call opennc('/usr/users/kyh/result4/2ndgen/105jul/UNbar/unbar.nc',ncid)
  call geta2d(ncid, "ubar", 1,1, 2,nzp0+1, ubardat)
  call geta2d(ncid, "nbar", 1,1, 2,nzp0+1, n2bardat)
  call closenc(ncid)
!!!!!!!!!!!!!!!!!
ubardat = 0.
n2bardat = 0.02
!!!!!!!!!!!!!!!!!
  n2bardat(:) = n2bardat(:)**2
  n2int(1) = 0.
  do k=2, nzp0
    n2int(k) = n2int(k-1) + 0.5*(n2bardat(k-1)+n2bardat(k))*dz
  enddo
  do k=1, nzp0
    if ((k-0.5)*dz .ge. zst) then
      tmpi = k
      EXIT
    end if
  enddo
  ubar(:) = ubardat(tmpi:tmpi+knd(nzre)-kst(1))
  n2bar(:) = n2bardat(tmpi:tmpi+knd(nzre)-kst(1))
  ptbar(:) = ptbar1 * exp(n2int(tmpi:tmpi+knd(nzre)-kst(1))/9.8)
print*,ubar

  ! k, omega
  call wavnum(nxp,dx,kk2)
  call wavnum(nt,dt,oo2)

  do i=1, nxp/2+1
    kk(i) = kk2(i)
  end do
  do n=1, nt/2+1
    oo(n) = oo2(nt/2+2-n)
  end do
  do n=nt/2+2, nt/2*2+1
    oo(n) = oo2(nt+nt/2+2-n)
  end do
  oo(:) = -oo(:)

  cc(1,:) = 0.
  do n=1, noo
  do i=2, nkk
    cc(i,n) = oo(n) / kk(i)
  enddo
  enddo

  ! m2 * (U-c)^2
  allocate(dudz(kst(1):knd(nzre)),d2udz(kst(1):knd(nzre)))
  do k=kst(1)+1, knd(nzre)-1
    dudz(k) = (ubar(k+1)-ubar(k-1)) / (2.*dz)
  enddo
  dudz(kst(1)) = dudz(kst(1)+1)
  dudz(knd(nzre)) = dudz(knd(nzre)-1)
  do k=kst(1)+1, knd(nzre)-1
    d2udz(k) = (dudz(k+1)-dudz(k-1)) / (2.*dz)              ! smoothing
  enddo
  d2udz(kst(1)) = 0.
  d2udz(knd(nzre)) = 0.

  allocate(mch2(nkk,noo,kst(1):knd(nzre)))
  mch2(1,:,:) = 0.
  do k=kst(1), knd(nzre)
  do n=1, noo
  do i=2, nkk
    mch2(i,n,k) = n2bar(k) - (ubar(k)-cc(i,n))*(d2udz(k)+dudz(k)/h_scale) &
                 - (ubar(k)-cc(i,n))**2*(kk(i)**2+0.25/h_scale**2)
  enddo
  enddo
  enddo
!!
do k=kst(1), knd(nzre)
do n=1, noo
do i=2, nkk
if (mch2(i,n,k) .eq. 0.) then
print*, i,n,k
endif
enddo
enddo
enddo
!!

  ! calculate S
  allocate(f_tot(nxp,nt,kst(1):knd(nzre)))
  allocate(fp(nxp,nt))
  allocate(coef(nxp,nt))
  allocate(ss(nkk,noo,kst(1):knd(nzre)),s_tot(nkk,noo,kst(1):knd(nzre)))
  allocate(ss2(nkk,noo,kst(1):knd(nzre)))
  s_tot(:,:,:) = (0.,0.)

  ! calculate fw_hat
  do n=1, nt
    write(time,'(i6.6)') int((nst-1+n)*dt)
    open(2,file='/export32/kyh/ing/105jul/fw/fw.bin'//time, &
         form='unformatted',convert='big_endian')
    read(2) fdata
    close(2)

    do k=kst(1), knd(nzre)
    do i=1, nxp
      f_tot(i,n,k) = fdata(ist+i-1,k)
    enddo
    enddo
  enddo

  do k=kst(1), knd(nzre)
    fp(:,:) = f_tot(:,:,k)
    call fft2df(nxp,nt,fp,coef)
    do n=1, nt/2+1
    do i=1, nxp/2+1
      ss(i,n,k) = coef(i,nt/2+2-n)
    enddo
    enddo
    do n=nt/2+2, nt/2*2+1
    do i=1, nxp/2+1
      ss(i,n,k) = coef(i,nt+nt/2+2-n)
    enddo
    enddo
  enddo

  ss(1,:,:) = 0.
  do k=kst(1), knd(nzre)
  do n=1, noo
  do i=2, nkk
    ss(i,n,k) = ss(i,n,k) * (0.,1.)*kk(i)*(ubar(k)-cc(i,n))
  enddo
  enddo
  enddo

  s_tot = s_tot + ss

  ! calculate fu_hat
  do n=1, nt
    write(time,'(i6.6)') int((nst-1+n)*dt)
    open(2,file='/export32/kyh/ing/105jul/fu/fu.bin'//time, &
         form='unformatted',convert='big_endian')
    read(2) fdata
    close(2)

    do k=kst(1), knd(nzre)
    do i=1, nxp
      f_tot(i,n,k) = fdata(ist+i-1,k)
    enddo
    enddo
  enddo

  do k=kst(1), knd(nzre)
    fp(:,:) = f_tot(:,:,k)
    call fft2df(nxp,nt,fp,coef)
    do n=1, nt/2+1
    do i=1, nxp/2+1
      ss2(i,n,k) = coef(i,nt/2+2-n)
    enddo
    enddo
    do n=nt/2+2, nt/2*2+1
    do i=1, nxp/2+1
      ss2(i,n,k) = coef(i,nt+nt/2+2-n)
    enddo
    enddo
  enddo

  do k=kst(1)+1, knd(nzre)-1
    ss(:,:,k) = (ss2(:,:,k+1)-ss2(:,:,k-1)) / (2.*dz)
  enddo
  ss(:,:,kst(1)) = (ss2(:,:,kst(1)+1)-ss2(:,:,kst(1))) / dz
  ss(:,:,knd(nzre)) = (ss2(:,:,knd(nzre))-ss2(:,:,knd(nzre)-1)) / dz

  ss(1,:,:) = 0.
  do k=kst(1), knd(nzre)
  do n=1, noo
  do i=2, nkk
    ss(i,n,k) = -ss(i,n,k) * (ubar(k)-cc(i,n))
  enddo
  enddo
  enddo

  s_tot = s_tot + ss

  ! calculate ft_hat
  do n=1, nt
    write(time,'(i6.6)') int((nst-1+n)*dt)
    open(2,file='/export32/kyh/ing/105jul/ft/ft.bin'//time, &
         form='unformatted',convert='big_endian')
    read(2) fdata
    close(2)

    do k=kst(1), knd(nzre)
    do i=1, nxp
      f_tot(i,n,k) = fdata(ist+i-1,k)
    enddo
    enddo
  enddo

  do k=kst(1), knd(nzre)
    fp(:,:) = f_tot(:,:,k)
    call fft2df(nxp,nt,fp,coef)
    do n=1, nt/2+1
    do i=1, nxp/2+1
      ss(i,n,k) = coef(i,nt/2+2-n)
    enddo
    enddo
    do n=nt/2+2, nt/2*2+1
    do i=1, nxp/2+1
      ss(i,n,k) = coef(i,nt+nt/2+2-n)
    enddo
    enddo
  enddo

  ss(1,:,:) = 0.
  do k=kst(1), knd(nzre)
  do n=1, noo
  do i=2, nkk
    ss(i,n,k) = ss(i,n,k) * 9.8/ptbar(k)
  enddo
  enddo
  enddo

  s_tot = s_tot + ss

  do k=kst(1), knd(nzre)
    s_tot(:,:,k) = s_tot(:,:,k) * exp(-(z(k)-zref)/(2.*h_scale))
  enddo

  ! S / (U-c)^2 / m2
  s_tot(1,:,:) = 0.
  s_tot(2:,:,:) = s_tot(2:,:,:) / mch2(2:,:,:)

  ! PSD
  allocate(d_ave(nkk,noo,nzre))
  d_ave = 0.
  do kre=1, nzre
  do k=kst(kre), knd(kre)
    do n=1, noo
    do i=1, nkk
      d_ave(i,n,kre) = d_ave(i,n,kre) + 2. * cabs(s_tot(i,n,k))**2 / kint(kre)
    enddo
    enddo
  enddo
  enddo


  ! output
  do k=1, nzre
  do n=1, noo
  do i=1, nkk
    if (d_ave(i,n,k) .le. 1.e-32)  d_ave(i,n,k) = 0.
  enddo
  enddo
  enddo

  write(fname,'(a)') '../skoz.nc'
  call out3d(trim(fname),1,(/'PSD'/),d_ave,'k',nkk,kk,'o',noo,oo, &
                         'z',nzre,zre/1.e3,'PSD2D k,o')


  STOP

end
