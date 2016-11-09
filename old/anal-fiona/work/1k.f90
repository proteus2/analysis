 program k_zt

  use pwrspd
  use netcdfio

!--------------------------------------------------------
!  ntp  : number of data time to "analysis"
!  ntre : number of averaged time (result data)
!  dt   : time interval of datum to analysis
!  tst  : starting time to analysis
!  xst, zst : starting index
!--------------------------------------------------------  

  implicit none

  integer, parameter :: nx =1402, nz = 472, nxp =1200, nzp = 350, ntp = 480
  integer, parameter :: nkk = nxp/2+1, tst = 0, xst = 100, zst = 0
  integer, parameter :: eval = 0, ntre = 8
  real,    parameter :: dx = 1000., dt = 60., dz = 300.
  integer :: i,k,n,nre, ntave
  real    :: wprt(nx,nz)
  real, dimension(nxp)            :: w, x
  real, dimension(nkk)            :: kk, d
  real, dimension(nkk,nzp,ntre)   :: d_ave
  real    :: d_out(nkk,nzp+2,ntre), z_out(nzp+2), t_out(ntre)
  real    :: rhort(nz)
  character*6   :: time
  character*128 :: fname


  open(1,file='/usr/users/kyh/anal/work/rhoalljul')
  read(1,*) rhort

  ntave = ntp/ntre

  do i=1, nxp
    x(i) = (i-0.5)*dx
  enddo
  do k=1, nzp+2
    z_out(k) = dz * (k-1.5+zst) / 1000.
  enddo
  do n=1, ntre
    t_out(n) = (dt*ntave*n + tst) / 60.
  enddo


  d_ave(:,:,:)=0.
  do nre=1, ntre                   

    do n = 1, ntave
      write(time,'(i6.6)') int(tst + (ntave*(nre-1)+n)*dt)
      open(2,file='/export10/kyh/ing/105jul/bin/w.bin'//time, &
           form='unformatted',convert='big_endian')
      read(2) wprt 
      close(2)
      do k=1, nzp
        do i=1, nxp
          w(i) = wprt(i+1+xst,k+1+zst) * rhort(k+1+zst)
        enddo
        kk(:)=0.
        d(:)=0.
        call psd1d(nxp,3,x,dx,w,kk,d,1,eval)

        d_ave(:,k,nre) = d_ave(:,k,nre) + d(:)
      enddo
    enddo
   
  enddo
  d_ave(:,:,:) = d_ave(:,:,:) / ntave     

  print*, kk

! output
!  open(10,file='/usr/users/kyh/anal/work/k', form='unformatted')
!  write(10) kk
!  write(10) d_ave

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
