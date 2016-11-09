 program o_z

  use netcdfio
  use pwrspd

!--------------------------------------------------------
!  nt   : number of data time to "analysis"
!  dt   : time interval of datum to analysis
!  nxre : number of averaged x-group (result data)
!  tst  : starting time to analysis
!--------------------------------------------------------  

  implicit none

  integer, parameter :: nx =1002, nz = 422 , nxp = 800, nzp = 300, nt = 240
  integer, parameter :: noo = nt/2+1, tst = 3600, xst = 100, zst = 0
  integer, parameter :: eval = 0, nxre = 800
  real,    parameter :: dt = 60., dx = 1000., dz = 300.
  integer :: i,k,n, ire, nxave
  real    :: wprt(nxp,nzp,nt), temp(nx,nz)
  real,dimension(nt)   :: w, t
  real,dimension(noo)  :: oo, d
  real, dimension(noo,nzp,nxre) :: d_ave
  real    :: d_out(noo,nzp+2,nxre), z_out(nzp+2), x_out(nxre)
  real    :: rhort(nz)
  character*6 :: time
  character*128 :: fname


  open(1,file='/usr/users/kyh/anal/work/rhoallid1')
  read(1,*) rhort 

  nxave = nxp/nxre

  do i=1, nxre
    x_out(i) = (dx*nxave*i + xst-100) / 1000.
  enddo
  do k=1, nzp+2 
    z_out(k) = dz * (k-1.5+zst) / 1000. + 15.
  enddo 
  do n=1, nt
    t(n) = tst + n*dt
  enddo

  do n=1, nt                         
    write(time,'(i6.6)') int(t(n))
    open(2,file='/export22/kyh/ing/id1/bin/w.bin'//time, &
         form='unformatted',convert='big_endian')
    read(2) temp
    close(2)
    do k=1, nzp
      do i=1, nxp
        wprt(i,k,n) = temp(i+1+xst,k+1+zst) * rhort(k+1+zst)
      enddo
    enddo
  enddo

  d_ave(:,:,:)=0.
  do ire=1, nxre

    do k=1, nzp
     do i=1, nxave
      do n=1, nt
        w(n) = wprt(int(i+(ire-1)*nxave),k,n)
      enddo
      oo(:)=0.
      d(:)=0.
      call psd1d(nt,4,t,dt,w,oo,d,1,eval)

      d_ave(:,k,ire) = d_ave(:,k,ire) + d(:)
     enddo
    enddo

  enddo
  d_ave(:,:,:) = d_ave(:,:,:) / nxave

  print*, oo

! output
!  open(10,file='o', form='unformatted')
!  write(10) oo
!  write(10) d_ave

  do k=2, nzp+1
    d_out(:,k,:) = d_ave(:,k-1,:)
  enddo
  d_out(:,1,:) = d_out(:,2,:)
  d_out(:,nzp+2,:) = d_out(:,nzp+1,:)


  write(fname,'(a)') '/usr/users/kyh/anal/work/ozx.nc'
  call out3d(trim(fname),1,(/'PSD'/),d_out,'o',noo,oo,'z',nzp+2,z_out,       &
                         'x',nxre,x_out,'PSD1D omega')


  stop
 end
