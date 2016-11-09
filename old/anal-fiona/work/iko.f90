 program ko_z

  use pwrspd
  use netcdfio

!--------------------------------------------------------
!  ntp  : number of data time to "analysis"
!  dt   : time interval of datum to analysis
!  tst  : starting "time" to analysis
!  xst, zst : starting "index"
!--------------------------------------------------------

  implicit none

  integer, parameter :: nx =1002, nz = 352 , nxp = 750, nzp = 350, ntp = 180
  integer, parameter :: xst = 150, zst = 0, tst = 7200, nzre = 350
  integer, parameter :: nkk = nxp/2+1, noo = ntp/2*2+1        !!
  real,    parameter :: dx = 1000., dt = 60., dz = 300.
  integer :: i,k,n, kre, nzave
  real    :: wprt(nxp,nzp,ntp), temp(nx,nz)
  real    :: x(nxp), t(ntp), w(nxp,ntp), z_out(nzre)
  real    :: kk(nkk), oo(noo), d(nkk,noo), d_ave(nkk,noo,nzre), tem(nkk,noo)
  real    :: rhort(nz)
  character*6   :: time
  character*128 :: fname


  open(1,file='/usr/users/kyh/anal/work/rhoallid0')
  read(1,*) rhort 

  nzave = nzp/nzre

  do i=1, nxp
    x(i) = (i-0.5)*dx
  enddo
  do n=1, ntp
    t(n) = tst + n*dt
  enddo
  do k=1, nzre
    z_out(k) = dz * (zst-0.5+nzave*k) / 1000.
  enddo


  do n=1, ntp                              ! time range
    write(time,'(i6.6)') int(t(n))
    open(2,file='/export22/kyh/ing/idu/bin/w.bin'//time, &
         form='unformatted',convert='big_endian')
    read(2) temp
    close(2)
    do k=1, nzp
      do i=1, nxp
        wprt(i,k,n) = temp(i+1+xst,k+1+zst)*rhort(k+1+zst)
      enddo
    enddo
  enddo

  d_ave(:,:,:)=0.
  do kre=1, nzre

    do k=1, nzave
      do i=1, nxp
        do n=1, ntp
          w(i,n) = wprt(i,int(k+(kre-1)*nzave),n)
        enddo
      enddo
      kk(:)=0.
      oo(:)=0.
      d(:,:)=0.
      call psd2d(nxp,ntp,dx,dt,w,kk,oo,d)

      tem(:,:) = d(:,:)
      do n=1, noo
        d(:,n) = tem(:,noo+1-n)
      enddo

      d_ave(:,:,kre) = d_ave(:,:,kre) + d(:,:)
    enddo

  enddo
  d_ave(:,:,:) = d_ave(:,:,:) / nzave          ! / total...

! output
!  open(10,file='/usr/users/kyh/anal/work/ko',form='unformatted')
!  write(10) kk
!  write(10) oo
!  write(10) d_ave

  write(fname,'(a,i5.5,a)') '/usr/users/kyh/anal/work/koz.nc'
  call out3d(trim(fname),1,(/'PSD'/),d_ave,'k',nkk,kk,'o',noo,oo, &
             'z',nzre,z_out,'PSD2D k,omega')


  stop

 end
