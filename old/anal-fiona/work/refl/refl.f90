program reflection


 implicit none
!==========================================================

! There must has the data between boundaries of a layer.

!==========================================================
 integer, parameter :: ikk= 28, ioo= 213
 integer, parameter :: n_layer = 5       ! change the bnd too
 integer, parameter :: nzr=202, nxp=700, nt=360
 real,    parameter :: dx=1000., dt=60., dz=300.
 integer, parameter :: nx=902, nz=352
!==========================================================
 integer, parameter :: nkk=nxp/2+1, noo = nt/2*2+1
 real    :: dkk=1./(nxp*dx), doo=1./(nt*dt)
 real    :: kk, oo
 real, dimension(n_layer+1) :: bnd
! real, dimension(nz) :: ooint
 integer :: i,j,n, step
 integer :: jl(n_layer)

 real*8 :: ubardat(nx,nz) 
 real :: m2(nkk,noo,nzr)
 real :: pi
 complex, dimension(n_layer) :: m2c, m
 complex :: ii
 
 pi = acos(-1.)
 ii = (0., 1.)

! specify the altitudes of the boundaries. ==================
 bnd = (/15., 24., 33., 42., 51., 60./)! + 0.15
!============================================================
 bnd = bnd * 1000.

 do j=1, n_layer
   jl(j) = int(((bnd(j)+bnd(j+1))*0.5 + 1.5*dz) / dz)
 enddo

! read m2
 open(1,file='/usr/users/kyh/anal/work/sep60/data/m2.bin',           &
        form='unformatted')
 read(1) m2
 close(1)
 do j=1, n_layer                           ! m sign.... think again
   m2c(j) = m2(ikk,ioo,jl(j))
   m(j) = sqrt(m2c(j)) * 2*pi
 enddo

! read U and determin intrinsic freq.
 open(2,file='/usr/users/kyh/result4/60src/bin/ubar.bin',            &
        form='unformatted',convert='big_endian')
 read(2) ubardat
 close(2)
 if ( (mod(int(bnd(1)),int(dz))) .eq. 0 ) then
  print*,1

 endif
! ubnd


! step = 1
! if (startn .gt. endn)  step = -1

! m = -1.
 
 


 end
