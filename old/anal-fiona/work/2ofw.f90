 program oo_z

  use netcdfio
  use pwrspd

  implicit none

  integer, parameter :: nx =1202, nz = 302, nxp =1200, nzp = 300, nt = 360
  integer, parameter :: nkk = nxp/2+1, noo = nt/2*2+1        !!
  integer, parameter :: tst = 7200, xst =   0, zst = 0
  integer, parameter :: nzall = 472
  real,    parameter :: zbottom = 30.e3, zref = 75.e3
  real,    parameter :: dx = 1000., dt = 60., dz = 300.

  integer :: i,k,n, zini
  real    :: fw(nxp,nzp,nt), f1prt(nx,nz), f2prt(nx,nz)
  real    :: x(nxp), t(nt), fwxt(nxp,nt)
  real    :: kk(nkk), oo(noo), d(nkk,noo), d_ave(noo,nzp), tem(nkk,noo)
  real    :: d_out(noo,nzp+2), z_out(nzp+2)
  real    :: rhort(nzall), ri(1402,472), ric
  character*6   :: time
  character*128 :: fname


  open(1,file='/usr/users/kyh/anal/work/rhoalljul')
  read(1,*) rhort 

  do i=1, nxp
    x(i) = (i-0.5)*dx
  enddo
  do n=1, nt
    t(n) = tst + n*dt
  enddo
  do k=1, nzp+2
    z_out(k) = (dz*(k-1.5+zst)+zbottom) / 1000.
  enddo
  zini = int(zbottom/dz)
  ric = 1./3.


  do n=1, nt                              ! time range
    write(time,'(i6.6)') int(t(n))
    open(2,file='/export10/kyh/ing/105jul/fw/fw1.bin'//time, &
         form='unformatted',convert='big_endian')
    read(2) f1prt
    close(2)
    open(3,file='/export10/kyh/ing/105jul/fw/fw2.bin'//time, &
         form='unformatted',convert='big_endian')
    read(3) f2prt
    close(3)

    open(4,file='/export10/kyh/ing/105jul/bin/ri.bin'//time, &
         form='unformatted',convert='big_endian')
    read(4) ri
    close(4)

    do k=1, nzp
      do i=1, nxp
        fw(i,k,n) = f1prt(i+1+xst,k+1+zst) + f2prt(i+1+xst,k+1+zst)
!!! only in the breaking region !!!!!!
        if (ri(i+1+100+xst,k+1+zini+zst) .ge. ric)  fw(i,k,n) = 0.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      enddo
      fw(:,k,n) = fw(:,k,n) * rhort(k+1+zst+int(zbottom/dz)) / rhort(2+int(zref/dz))
    enddo
  enddo

  d_ave(:,:)=0.
  do k=1, nzp
    do n=1, nt
      do i=1, nxp
        fwxt(i,n) = fw(i,k,n)
      enddo
    enddo
    kk(:)=0.
    oo(:)=0.
    d(:,:)=0.
    call psd2d(nxp,nt,dx,dt,fwxt,kk,oo,d)

    tem(:,:) = d(:,:)
    do n=1, noo
      d(:,n) = tem(:,noo+1-n)
    enddo 

    do n=1, noo
      do i=1, nkk
        d_ave(n,k) = d_ave(n,k) + d(i,n)
      enddo
      d_ave(n,k) = d_ave(n,k) / (nxp*dx)
    enddo

  enddo

  print*, oo

  do k=2, nzp+1
    d_out(:,k) = d_ave(:,k-1)
  enddo
  d_out(:,1) = d_out(:,2)
  d_out(:,nzp+2) = d_out(:,nzp+1)


  write(fname,'(a)') '/usr/users/kyh/anal/work/ooz.nc'
  call out2d(trim(fname),1,(/'PSD'/),d_out,'o',noo,oo,'z',nzp+2,z_out,       &
                         'PSD2D omega')

  stop
 end
