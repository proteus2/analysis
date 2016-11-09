 program ko_z

  use pwrspd
  use netcdfio

!--------------------------------------------------------
!  ntp  : number of data time to "analysis"
!  dt   : time interval of datum to analysis
!  tst  : starting "time" to analysis
!  xst, zst : starting "index"
!  zbottom  : bottom height of the domain in binary data
!  zref : normalizing height
!--------------------------------------------------------

  implicit none

  integer, parameter :: nx =1202, nz = 302 , nxp =1200, nzp = 300, ntp = 360
  integer, parameter :: nkk = nxp/2+1, noo = ntp/2*2+1        !!
  integer, parameter :: xst =   0, zst = 0, tst = 7200
  integer, parameter :: nzall = 472, nzre = 300
  real,    parameter :: zbottom = 30.e3, zref = 75.e3
  real,    parameter :: dx = 1000., dt = 60., dz = 300.

  integer :: i,k,n, kre, nzave, zini
  real    :: fw(nxp,nzp,ntp), f1prt(nx,nz), f2prt(nx,nz)
  real    :: x(nxp), t(ntp), fwxt(nxp,ntp), z_out(nzre)
  real    :: kk(nkk), oo(noo), d(nkk,noo), d_ave(nkk,noo,nzre), tem(nkk,noo)
  real    :: rhort(nzall), ri(1402,472), ric
  character*6   :: time
  character*128 :: fname


  open(1,file='/usr/users/kyh/anal/work/rhoalljul')
  read(1,*) rhort 

  nzave = nzp/nzre

  do i=1, nxp
    x(i) = (i-0.5)*dx
  enddo
  do n=1, ntp
    t(n) = tst + n*dt
  enddo
  do k=1, nzre
    z_out(k) = (dz*(zst-0.5+nzave*k)+zbottom) / 1000.
  enddo
  zini = int(zbottom/dz)
  ric = 1./3.


  do n=1, ntp                              ! time range
    write(time,'(i6.6)') int(t(n))
    open(2,file='/export10/kyh/ing/105jul/fw/fw1.bin'//time, &
         form='unformatted',convert='big_endian')
    read(2) f1prt
    close(2)
    open(3,file='/export10/kyh/ing/105jul/fw/fw2.bin'//time, &
         form='unformatted',convert='big_endian')
    read(3) f2prt
    close(3)
    open(4,file='/export10/kyh/ing/105jul/fw/wtmix.bin'//time, &
         form='unformatted',convert='big_endian')
    read(4) f2prt
    close(4)
    open(5,file='/export10/kyh/ing/105jul/bin/ri.bin'//time, &
         form='unformatted',convert='big_endian')
    read(5) ri
    close(5)

    do k=1, nzp
      do i=1, nxp
        fw(i,k,n) = f1prt(i+1+xst,k+1+zst)+f2prt(i+1+xst,k+1+zst)
!!! only in the breaking region !!!!!!
        if (ri(i+1+100+xst,k+1+zini+zst) .ge. ric)  fw(i,k,n) = 0.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      enddo
      fw(:,k,n) = fw(:,k,n) * rhort(k+1+zst+int(zbottom/dz)) / rhort(2+int(zref/dz))
    enddo
  enddo

  d_ave(:,:,:)=0.
  do kre=1, nzre

    do k=1, nzave
      do i=1, nxp
        do n=1, ntp
          fwxt(i,n) = fw(i,int(k+(kre-1)*nzave),n)
        enddo
      enddo
      kk(:)=0.
      oo(:)=0.
      d(:,:)=0.
      call psd2d(nxp,ntp,dx,dt,fwxt,kk,oo,d)

      tem(:,:) = d(:,:)
      do n=1, noo
        d(:,n) = tem(:,noo+1-n)
      enddo

      d_ave(:,:,kre) = d_ave(:,:,kre) + d(:,:)
    enddo

  enddo
  d_ave(:,:,:) = d_ave(:,:,:) / nzave          ! / total...

! output
  write(fname,'(a,i5.5,a)') '/usr/users/kyh/anal/work/koz.nc'
  call out3d(trim(fname),1,(/'PSD'/),d_ave,'k',nkk,kk,'o',noo,oo, &
             'z',nzre,z_out,'PSD2D k,omega')


  stop

 end
