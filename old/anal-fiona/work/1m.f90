 program ubar_z

  use pwrspd
  implicit none

  integer, parameter :: nx = 902, nz = 472 , nxp = 700, nzp = 300, nt = 180
  integer, parameter :: nmm = nzp/2+1
  integer, parameter :: eval = 0
  real,    parameter :: dz = 300.
  integer :: i,k,n
  real    :: wprt(nx,nz)
  real,dimension(nzp)    :: w,z
  real,dimension(nmm)     :: mm,d,d_ave
  real    :: rhort(nz)
  character*6 :: time


  open(10,file='m')
  open(20,file='mave')
  open(1,file='rhoall')
  read(1,*) rhort 

  do k=1,nzp
    z(k) = (k+50-0.5)*dz
  enddo

  d_ave(:)=0.
  do n=1,nt                               ! time range
    write(time,'(i6.6)') 7200 + n*120
    open(2,file='../../result4/105t0/bin/w.bin'//time, &
         form='unformatted',convert='big_endian')
    read(2) wprt 
    close(2)
 
    do i=1, nxp
      do k=1, nzp
        w(k) = wprt(i+101,k+51) * rhort(k+51)
      enddo
      mm(:)=0.
      d(:)=0.
      call psd1d(nzp,4,z,dz,w,mm,d,eval)

      if (n==1.and.i==1) then
!       write(10,*) mm
      end if
!     write(10,*) d

      d_ave(:) = d_ave(:) + d(:)
    enddo
  enddo
  d_ave(:) = d_ave(:) / (nxp*nt)          ! / total...

  do k=2,nmm
    write(20,*) mm(k), d_ave(k)
  enddo


  stop
 end
