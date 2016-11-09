!2345
    module fft

!----------------------------------------------------------------------------------
!
!   Purpose:
!   To calculate Fourier coefficients of one-dimensional or two-dimensional 
!   data
! 
!   Author:
!   In-Sun Song
!   Laboratory for Mesoscale Dynamics
!   Department of Atmospheric Sciences, Yonsei University, Seoul, Korea
!
!   USAGE:
!   This module should be compiled with cfftpack.f90
!   First of all. cfftpack.f90 should be compiled prior to the compilation 
!   of this module
!   In this module, all varables are defined as double precisions.
!
!   f90 -c cfftpack.f90
!   f90 -c fft.f90 cfftpack.o
!
!   You will see two object file (cfftpack.o and fft.o) and module file
!   (cfftpack.mod)
!
!----------------------------------------------------------------------------------

    contains

    subroutine fft1df(nx,dx,data,freqx,coef)

    integer                            ,intent(in)  :: nx
    real                               ,intent(in)  :: dx
    real,            dimension(nx)     ,intent(in)  :: data
    real,            dimension(nx)     ,intent(out) :: freqx
    double complex  ,dimension(nx)     ,intent(out) :: coef

    integer                             :: i
    double precision, dimension(4*nx+15) :: work

    external cffti, cfftf

    work(1:4*nx+15) = 0.d0
 
    do i=1,nx/2+1
      freqx(i) = float(i-1)/(nx*dx)
    end do
    do i=nx/2+2,nx 
      freqx(i) = -freqx(nx-i+2)
    end do
 
    do i=1,nx
      coef(i) = dble(data(i))
    end do

    call CFFTI(nx,work)
    call CFFTF(nx,coef,work)

    return 
    end subroutine fft1df

    subroutine fft1db(nx,coef,data)

    integer                            ,intent(in)  :: nx
    real,            dimension(nx)     ,intent(out) :: data
    double complex  ,dimension(nx)     ,intent(in)  :: coef

    integer                             :: i
    double precision, dimension(4*nx+15) :: work

    external cffti, cfftb

    work(1:4*nx+15) = 0.d0

    call CFFTI(nx,work)
    call CFFTB(nx,coef,work)

    do i=1,nx
      data(i) = real( (conjg(coef(i))+coef(i))*0.5 )
    end do

    return
    end subroutine fft1db

    subroutine fft2df(nx,ny,dx,dy,data,freqx,freqy,coef)

    integer                            ,intent(in)  :: nx     ,ny
    real                               ,intent(in)  :: dx     ,dy
    real,            dimension(nx,ny)  ,intent(in)  :: data
    real,            dimension(nx)     ,intent(out) :: freqx
    real,            dimension(ny)     ,intent(out) :: freqy
    double complex  ,dimension(nx,ny)  ,intent(out) :: coef

    integer                              :: i      ,j
    double precision, dimension(4*nx+15) :: workx
    double precision, dimension(4*ny+15) :: worky
    double complex  , dimension(nx)      :: coefx
    double complex  , dimension(ny)      :: coefy

    external cffti, cfftf
 
    do i=1,nx/2+1
      freqx(i) = float(i-1)/(nx*dx)
    end do
    do i=nx/2+2,nx
      freqx(i) = -freqx(nx-i+2)
    end do

    do j=1,ny/2+1
      freqy(j) = float(j-1)/(ny*dy)
    end do
    do j=ny/2+2,ny
      freqy(j) = -freqy(ny-j+2)
    end do

    do j=1,ny
      do i=1,nx
        coef(i,j) = dble(data(i,j))
      end do
    end do

    do j=1,ny
      workx(1:4*nx+15) = 0.d0
      coefx(1:nx) = coef(1:nx,j)
      call cffti(nx,workx)
      call cfftf(nx,coefx,workx)
      coef(1:nx,j) = coefx(1:nx)
    end do

    do i=1,nx
      worky(1:4*ny+15) = 0.d0
      coefy(1:ny) = coef(i,1:ny)
      call cffti(ny,worky)
      call cfftf(ny,coefy,worky)
      coef(i,1:ny) = coefy(1:ny)
    end do 
 
    return 
    end subroutine fft2df

    subroutine fft2db(nx,ny,coef,data)

    integer                            ,intent(in)  :: nx     ,ny
    double complex  ,dimension(nx,ny)  ,intent(in)  :: coef
    real            ,dimension(nx,ny)  ,intent(out) :: data

    integer                              :: i      ,j
    double precision, dimension(4*nx+15) :: workx
    double precision, dimension(4*ny+15) :: worky
    double complex  , dimension(nx)      :: coefx
    double complex  , dimension(ny)      :: coefy
    double complex  , dimension(nx,ny)   :: coef2

    external cffti, cfftb

    coef2 = coef

    do j=1,ny
      workx(1:4*nx+15) = 0.d0
      coefx(1:nx) = coef2(1:nx,j)
      call cffti(nx,workx)
      call cfftb(nx,coefx,workx)
      coef2(1:nx,j) = (1./float(nx))*coefx(1:nx)
    end do

    do i=1,nx
      worky(1:4*ny+15) = 0.d0
      coefy(1:ny) = coef2(i,1:ny)
      call cffti(ny,worky)
      call cfftf(ny,coefy,worky)
      coef2(i,1:ny) = (1./float(ny))*coefy(1:ny)
    end do

    do i=1,nx
      do j=1,ny
        data(i,j) = real(coef2(i,j))
      end do
    end do 

    return
    end subroutine fft2db


    end module fft
