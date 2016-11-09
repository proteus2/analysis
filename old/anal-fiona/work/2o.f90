 program oo_z

  use netcdfio
  use pwrspd

  implicit none

  integer, parameter :: nx =1402, nz = 472, nxp =1200, nzp = 350, nt = 360
  integer, parameter :: nkk = nxp/2+1, noo = nt/2*2+1        !!
  integer, parameter :: tst = 7200, xst = 100, zst = 0
  real,    parameter :: dx = 1000., dt = 60., dz = 300.
  integer :: i,k,n
  real    :: wprt(nxp,nzp,nt), temp(nx,nz)
  real    :: x(nxp), t(nt), w(nxp,nt)
  real    :: kk(nkk), oo(noo), d(nkk,noo), d_ave(noo,nzp), tem(nkk,noo)
  real    :: d_out(noo,nzp+2), z_out(nzp+2)
  real    :: rhort(nz)
  character*6   :: time
  character*128 :: fname


  open(1,file='/usr/users/kyh/anal/work/rhoalljult1')
  read(1,*) rhort 

  do i=1, nxp
    x(i) = (i-0.5)*dx
  enddo
  do n=1, nt
    t(n) = tst + n*dt
  enddo
  do k=1, nzp+2
    z_out(k) = dz * (k-1.5+zst) / 1000.
  enddo


  do n=1, nt                              ! time range
    write(time,'(i6.6)') int(t(n))
    open(2,file='/export22/kyh/ing/105jult1/bin/w.bin'//time, &
         form='unformatted',convert='big_endian')
    read(2) temp
    close(2)
    do k=1, nzp
      do i=1, nxp
        wprt(i,k,n) = temp(i+1+xst,k+1+zst) * rhort(k+1+zst)
      enddo
    enddo
  enddo

  d_ave(:,:)=0.
  do k=1, nzp
    do n=1, nt
      do i=1, nxp
        w(i,n) = wprt(i,k,n)
      enddo
    enddo
    kk(:)=0.
    oo(:)=0.
    d(:,:)=0.
    call psd2d(nxp,nt,dx,dt,w,kk,oo,d)

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
