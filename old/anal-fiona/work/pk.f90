program k_spectrum

  use pwrspd
  use netcdfio

  implicit none

  integer, parameter :: tkz = 0, kzt = 1
  integer, parameter :: nx =1002, nz = 472
  real,    parameter :: dx = 1.e3, dt = 60., dz = 300., zref = 95.e3

  integer, parameter :: nst =  50,    nnd =  50,    ntint =  1
  real,    parameter :: xst = 250.e3, xnd = 600.e3
  real,    parameter :: zst =  90.e3, znd = 100.e3, zint = 10.e3
  real,    parameter :: x0 = -100.e3, z0 =  0.e3

  integer, parameter :: eval = 0
  integer :: i,k,n, tag1, tag2, nre, kre, kfind
  integer :: nxp, ist, ind, nkk, ntre, nzre

  real                   :: wprt(nx,nz), x(nx), z(nz)
  real                   :: rhobar(472), rhort(nz), rhostd, zfind
  character*6   :: time
  character*128 :: fname

  integer, dimension(:),  allocatable :: kst, knd
  real, dimension(:),     allocatable :: w, d, kk, tre, zre, kint
  real, dimension(:,:,:), allocatable :: d_ave, d_out

  
  open(1,file='/usr/users/kyh/ideal/work/id0/rhobar')
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
  allocate(w(nxp), d(nkk), kk(nkk))

  ! analysis domain (z)
  if ((znd-zst)/zint-int((znd-zst)/zint) .ne. 0) then
    print *, '( ( zt - zb ) / interval ) must be an integer.'
    STOP
  end if

  ntre = int((nnd-nst+1)/ntint)
  allocate(tre(ntre))
  do nre=1, ntre
    tre(nre) = (nst-1 + nre*ntint - (ntint/2.-0.5)) * dt
  enddo
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
  allocate(d_ave(nkk,nzre,ntre))                  ;  d_ave(:,:,:)=0.
  do nre=1, ntre
  do n=1, ntint
    write(time,'(i6.6)') int((nst-1+(nre-1)*ntint+n)*dt)
    open(2,file='/export24/kyh/ideal/id0/bin/w.bin'//time, &
         form='unformatted',convert='big_endian')
    read(2) wprt
    close(2)

    do kre=1, nzre
    do k=kst(kre), knd(kre)
      do i=1, nxp
        w(i) = wprt(ist+i-1,k) * rhort(k)
      enddo
      kk(:)=0.
      d(:)=0.
      call psd1d(nxp,3,x(:nxp),dx,w,kk,d,1,eval)

      d_ave(:,kre,nre) = d_ave(:,kre,nre) + d(:) / (kint(kre)*ntint)
    enddo
    enddo
  enddo
  enddo

  print*, kk

  do n=1, ntre
  do k=1, nzre
  do i=1, nkk
    if (d_ave(i,k,n) .lt. 1.e-32)   d_ave(i,k,n) = 0.
  enddo
  enddo
  enddo


  ! output
  if (kzt .eq. 1) then
    allocate(d_out(nkk,nzre,ntre))                   ; d_out = 0.
    d_out = d_ave

    write(fname,'(a)') '../kzt.nc'
    call out3d(trim(fname),1,(/'PSD'/),d_out,'k',nkk,kk,'z',nzre,zre/1.e3, &
                           't',ntre,tre/60.,'PSD1D k')
    deallocate(d_out)
  end if

  if (tkz .eq. 1) then
    allocate(d_out(ntre,nkk,nzre))                   ; d_out = 0.
    do k=1, nzre
    do i=1, nkk
    do n=1, ntre
      d_out(n,i,k) = d_ave(i,k,n)
    enddo
    enddo
    enddo

    write(fname,'(a)') '../tkz.nc'
    call out3d(trim(fname),1,(/'PSD'/),d_out,'t',ntre,tre/60.,'k',nkk,kk, &
                           'z',nzre,zre/1.e3,'PSD1D k')
    deallocate(d_out)
  end if


  STOP

end
