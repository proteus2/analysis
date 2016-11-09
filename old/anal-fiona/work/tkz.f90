program k_tz

  use pwrspd
  use netcdfio

!--------------------------------------------------------
!  ntp  : number of data to "analyze"
!  dt   : time interval of datum to analyze
!  tst  : index to skip analysis
!         ( analyze [tst+1,tst+ntp] )
!         if tst=-1, PSD=0 at t=0
!--------------------------------------------------------  

  implicit none

  integer, parameter :: nx =1002, nz = 472, tst = -1, ntp = 300
  real,    parameter :: dx = 1.e3, dt = 60., dz = 300., zref = 15.e3

  real,    parameter :: xst = 250.e3, xnd = 600.e3
  real,    parameter :: zst = 20.e3,  znd = 100.e3, zint = 10.e3
  real,    parameter :: x0 = -100.e3, z0 = 0.e3
  integer, parameter :: eval = 0

  integer :: i,k,n, tag1, tag2, kre, kfind, nst
  integer :: nxp, ist, ind, nkk, nzre

  real    :: wprt(nx,nz), x(nx), z(nz), t_out(ntp)
  real    :: rhort(nz), rhostd, zfind
  character*6   :: time
  character*128 :: fname

  real, dimension(:),     allocatable :: w, d, kk, zre, kst, knd, kint
  real, dimension(:,:,:), allocatable :: d_ave, d_out

  
  open(1,file='/data/kyh/model/ideal/work/idn/rhobar')
  read(1,*) rhort
  rhostd = rhort(2+int(zref/dz))
  rhort(:) = sqrt(rhort(:)/rhostd)


  ! domain (input)
  do i=1, nx
    x(i) = x0 + (i-1.5)*dx
  enddo
  do k=1, nz
    z(k) = z0 + (k-1.5)*dz
  enddo
  do n=1, ntp
    t_out(n) = dt * (n+tst) / 60.
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
  kint(:) = knd(:) - kst(:) + 1

  ! calculate PSD
  allocate(d_ave(nkk,nzre,ntp))                  ;  d_ave(:,:,:)=0.
  nst = 1
  if (tst .eq. -1)  nst = 2
  do n=nst, ntp
    write(time,'(i6.6)') int(t_out(n)*60.)
    open(2,file='/data/kyh/result4/idn/bin/w.bin'//time, &
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
        call psd1d(nxp,3,x,dx,w,kk,d,1,eval)

        d_ave(:,kre,n) = d_ave(:,kre,n) + d(:)
      enddo
      d_ave(:,kre,n) = d_ave(:,kre,n) / kint(kre)
    enddo
  enddo

  print*, kk

  ! output
  allocate(d_out(ntp,nkk,nzre))
  do n=1, ntp
    d_out(n,:,:) = d_ave(:,:,n)
  enddo

  write(fname,'(a)') '../tkz.nc'
  call out3d(trim(fname),1,(/'PSD'/),d_out,'t',ntp,t_out,'k',nkk,kk, &
                         'z',nzre,zre,'PSD1D k')


  STOP

end
