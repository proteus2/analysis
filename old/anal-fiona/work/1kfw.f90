 program k_zt

  use pwrspd
  use netcdfio

!--------------------------------------------------------
!  ntp  : number of data time to "analysis"
!  ntre : number of averaged time (result data)
!  dt   : time interval of datum to analysis
!  tst  : starting "time" to analysis
!  xst, zst : starting "index"
!  zbottom  : bottom height of the domain in binary data
!  zref : normalizing height
!--------------------------------------------------------  

  implicit none

  integer, parameter :: nx =1202, nz = 302, nxp =1200, nzp = 300, ntp = 480
  integer, parameter :: nkk = nxp/2+1, tst = 0, xst =   0, zst = 0
  integer, parameter :: nzall = 472, eval = 0, ntre = 8
  real,    parameter :: dx = 1000., dt = 60., dz = 300.
  real,    parameter :: zbottom = 30.e3, zref = 75.e3

  integer :: i,k,n,nre, ntave, zini
  real    :: f1prt(nx,nz), f2prt(nx,nz)
  real, dimension(nxp)            :: fw, x
  real, dimension(nkk)            :: kk, d
  real, dimension(nkk,nzp,ntre)   :: d_ave
  real    :: d_out(nkk,nzp+2,ntre), z_out(nzp+2), t_out(ntre)
  real    :: rhort(nzall), ri(1402,472), ric
  character*6   :: time
  character*128 :: fname


  open(1,file='/usr/users/kyh/anal/work/rhoalljul')
  read(1,*) rhort

  ntave = ntp/ntre

  do i=1, nxp
    x(i) = (i-0.5)*dx
  enddo
  do k=1, nzp+2
    z_out(k) = (dz*(k-1.5+zst)+zbottom) / 1000.
  enddo
  do n=1, ntre
    t_out(n) = (dt*ntave*n + tst) / 60.
  enddo
  zini = int(zbottom/dz)
  ric = 1./3.


  d_ave(:,:,:)=0.
  do nre=1, ntre                   

    do n = 1, ntave
      write(time,'(i6.6)') int(tst + (ntave*(nre-1)+n)*dt)
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
          fw(i) = f1prt(i+1+xst,k+1+zst)+f2prt(i+1+xst,k+1+zst)
!!! only in the breaking region !!!!!!
          if (ri(i+1+100+xst,k+1+zini+zst) .ge. ric)  fw(i) = 0.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        enddo
        fw(:) = fw(:) * rhort(k+1+zst+int(zbottom/dz)) / rhort(2+int(zref/dz))
        kk(:)=0.
        d(:)=0.
        call psd1d(nxp,3,x,dx,fw,kk,d,1,eval)

        d_ave(:,k,nre) = d_ave(:,k,nre) + d(:)
      enddo
    enddo
   
  enddo
  d_ave(:,:,:) = d_ave(:,:,:) / ntave     

  print*, kk

! output
  do k=2, nzp+1
    d_out(:,k,:) = d_ave(:,k-1,:)
  enddo
  d_out(:,1,:) = d_out(:,2,:)
  d_out(:,nzp+2,:) = d_out(:,nzp+1,:)

  write(fname,'(a)') '/usr/users/kyh/anal/work/kzt.nc'
  call out3d(trim(fname),1,(/'PSD'/),d_out,'k',nkk,kk,'z',nzp+2,z_out, &
                         't',ntre,t_out,'PSD1D k')


  stop

 end
