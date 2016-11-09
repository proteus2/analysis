program ko_spectrum

  use pwrspd
  use netcdfio

  implicit none

  integer, parameter :: nx =1402, nz = 472
  real,    parameter :: dx = 1.e3, dt = 60., dz = 300., zref = 15.e3

  integer, parameter :: nst =   1,    nnd = 300
  real,    parameter :: xst =   0.e3, xnd =1200.e3
  real,    parameter :: zst =   0.e3, znd = 105.e3, zint = 0.3e3
  real,    parameter :: x0 = -100.e3, z0 =  0.e3

  integer, parameter :: eval = 0
  integer, parameter :: nt = nnd-nst+1
  integer, parameter :: noo = nt/2*2+1
  integer :: i,k,n, tag1, tag2, kre, kfind
  integer :: nxp, ist, ind, nkk, nzre
  real, dimension(nx,nz) :: wprt
  real                   :: x(nx), z(nz), t(nt)
  real                   :: rhobar(472), rhort(nz), rhostd, zfind
  real                   :: oo(noo)
  character*6   :: time
  character*128 :: fname

  integer, dimension(:),  allocatable :: kst, knd
  real, dimension(:),     allocatable :: zre, kint, kk
  real, dimension(:,:),   allocatable :: w, d
  real, dimension(:,:,:), allocatable :: w_tot, d_ave

  
  open(1,file='/export23/kyh/105jul/mfrc/frcj6/rhobar')
  read(1,*) rhobar
  rhostd = rhobar(2+int(zref/dz))
  do k=1, nz
    rhort(k) = sqrt(rhobar(int(z0/dz)+k)/rhostd)
  enddo


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
  allocate(kk(nkk))
  allocate(w(nxp,nt))
  allocate(d(nkk,noo))

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

  ! calculate PSD
  allocate(w_tot(nxp,nt,kst(1):knd(nzre)))
  do n=1, nt
    write(time,'(i6.6)') int((nst-1+n)*dt)
    open(2,file='/export23/kyh/105jul/mfrc/frcj6/bin/w.bin'//time, &
         form='unformatted',convert='big_endian')
    read(2) wprt
    close(2)

    do k=kst(1), knd(nzre)
    do i=1, nxp
      w_tot(i,n,k) = wprt(ist+i-1,k)
    enddo
    enddo
  enddo

  allocate(d_ave(nkk,noo,nzre))
  d_ave = 0.
  do kre=1, nzre
  do k=kst(kre), knd(kre)
    w(:,:) = w_tot(:,:,k) * rhort(k)
    kk = 0.  ;  oo = 0.
    d  = 0.
    call psd2d(nxp,nt,dx,dt,w,kk,oo,d)

    do n=1, noo
    do i=1, nkk
      d_ave(i,n,kre) = d_ave(i,n,kre) + d(i,noo+1-n) / kint(kre)
    enddo
    enddo
  enddo
  enddo

  print*, oo


  ! output
  do k=1, nzre
  do n=1, noo
  do i=1, nkk
    if (d_ave(i,n,k) .le. 1.e-32)  d_ave(i,n,k) = 0.
  enddo
  enddo
  enddo

  write(fname,'(a)') '../koz.nc'
  call out3d(trim(fname),1,(/'PSD'/),d_ave*1.e12,'k',nkk,kk,'o',noo,oo, &
                         'z',nzre,zre/1.e3,'PSD2D k,o')


  STOP

end
