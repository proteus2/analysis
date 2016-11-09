 program ubar_z

  use pwrspd
  implicit none

  integer, parameter :: nx = 902, nz = 472 , nxp = 700, nzp = 300, nt = 180
  integer, parameter :: nkk = nxp/2+1, nmm = nzp/2*2+1        !!
  real,    parameter :: dx = 1000., dz = 300., dt = 120
  integer :: i,k,n
  real    :: wprt(nx,nz,nt), temp(nx,nz)
  real    :: x(nxp), z(nzp), w(nxp,nzp), t(nt)
  real    :: kk(nkk), mm(nmm), d(nkk,nmm), d_ave(nkk,nmm)
  real    :: rhort(nzp)
  character*6 :: time


  open(10,file='km')
  open(1,file='rho')
  read(1,*) rhort 

  do i=1, nxp
    x(i) = (i-100-0.5)*dx
  enddo
  do k=1, nzp 
    z(k) = (k+ 50-0.5)*dz
  enddo
  do n=1, nt
    t(n) = (n+60)*dt
  enddo

  do n=1, nt                              ! time range
    write(time,'(i6.6)') int(t(n))
    open(2,file='/usr/users/kyh/result4/105t0/bin/w.bin'//time, &
         form='unformatted',convert='big_endian')
    read(2) temp
    close(2)
    do k=1, nzp
      do i=1, nxp
        wprt(i,k,n) = temp(i+101,k+51)*rhort(k)
      enddo
    enddo
  enddo

  d_ave(:,:)=0.
  do n=1, nt 
    do i=1, nxp
      do k=1, nzp
        w(i,k) = wprt(i,k,n)
      enddo
    enddo
    kk(:)=0.
    mm(:)=0.
    d(:,:)=0.
    call psd2d(nxp,nzp,dx,dz,w,kk,mm,d)

    if (n==1) then                                    !!
      write(10,*) kk
      write(10,*) mm
    end if
    write(10,*) d

    d_ave(:,:) = d_ave(:,:) + d(:,:)
  enddo
  d_ave(:,:) = d_ave(:,:) / nt           ! / total...


!  write(10,*) kk
!  write(10,*) ome
!  write(10,*) d_ave


  stop
 end
