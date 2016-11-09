    module pwrspd

!----------------------------------------------------------------------------------
!
!   PURPOSE:
!
!   To calculate power spectral density, cospectrum, and phase speed spectrum.
!
!
!   AUTHOR:
!
!   In-Sun Song.
!   Laboratory for Mesoscale Dynamics.
!   Department of Atmospheric Sciences, Yonsei University, Seoul, Korea.
!
!
!   TESTED SYSTEMS:
!
!   Compaq workstation: Compaq Unix 4.0 or 5.1.
!
!
!   VERSION HISTORY:
!
!   19 NOV, 2001: First wrote this module.
!   16 Feb, 2004: Reviewed and added descriptions.
!
!
!   REQUIRED MODULE:
!
!   fft, regress
!
!----------------------------------------------------------------------------------

    use fft,     only: fft1df,fft2df,wavnum
    use regress, only: lfit  ,funcs

    public :: psd1d 

    contains 

    subroutine psd1d(ndata ,ma    ,x     ,dx    ,data  ,wavk  ,pwsd  ,eval  )

!-------------------------------------------------------------------------------
!
!   SUBROUTINE psd1d
!
!   To calculate one-dimensional power spectrum density.
!
!   INPUT DATA
!
!   ndata : The number of elements of time series vector to be analysized.
!   ma    : The number of coefficients of polynomial to be used to remove
!           a trend in time series
!           ma = 2 : Linear trend
!           ma = 3 : Second order polynomial trend
!           ma cannot exceed 50
!   x     : Time ( or distance )
!   dx    : dt ( or dx )
!   data  : The time series vector
!
!   OUTPUT DATA
!
!   wavk  : Wavenumber vector      : wavk(ndata/2+1)
!   pwsd  : Power spectral density : pwsd(ndata/2+1)
!
!
!   DESCRIPTION:
!
!   1. Mean of given time series is removed
!   2. Trend is removed in the time series
!   3. Welch window function is applied to the modified time series to prevent
!      energy leakage in each wavenumber bin
!   4. Calculate power spectral density
!
!-------------------------------------------------------------------------------

    implicit none

    integer,                       intent(in)  :: ndata ,ma    ,eval
    real,                          intent(in)  :: dx
    real,    dimension(ndata),     intent(in)  :: x     ,data
    real,    dimension(ndata/2+1), intent(out) :: wavk  ,pwsd 

    integer                  :: i     ,j
    integer, dimension(ma)   :: ia
    real                     :: chisq ,summ  ,an    ,bn    ,wss   ,sqr   ,sqf
    real, dimension(ndata)   :: wavk2
    real, dimension(ma)      :: a
    real, dimension(ndata)   :: data1 ,sig   ,func  ,trend
    real, dimension(ndata)   :: window
    real, dimension(ma+10,ma+10) :: covar
    real, dimension(50)      :: afunc

    double precision                    :: pi
    double complex                      :: ci
    double complex, dimension(ndata)    :: coef
   
    pi = 4.d0*datan(1.d0)
    ci = (0.d0,1.d0) 


    data1(:) = data(:)
!++ 1. Removing mean            (when not de-trending)
!  
!   summ = 0.0
!
!   if (ma .eq. 0) then
!     do i=1,ndata
!       summ  = summ  + data(i)
!     end do 
!     summ = summ/float(ndata)
!   end if
!
!   do i=1,ndata
!     data1(i) = data(i) - summ
!   end do

!++ 2. Removing trend in data 

    do i=1,ma
     ia(i) = 1           ! calculate the coeff.(=a) if ia.ne.0
      a(i) = 0.0         ! coeff. for ia=0
    end do
    do i=1,ndata
      sig(i) = 1.0
    end do
   
    call lfit(x,data1,sig,ndata,a,ia,ma,covar,ma+10,chisq)
!    print*, (a(i),i=1,ma)

    do i=1,ndata
      call funcs(x(i),afunc,ma)
      trend(i)=0.
      do j=1,ma
        trend(i) = trend(i) + a(j)*afunc(j)
      end do 
      data1(i) = data1(i) - trend(i)

      if (eval .eq. 1) then
        open(777,file='../work/trend')
        write(777,*) x(i),data(i),data1(i),trend(i)
      end if

    end do 

!++ 3. Define Welch window

    an = ( ndata - 1.0 )/2.
    bn = ( ndata + 1.0 )/2.

    wss = 0.0
    do i=1,ndata
      window(i) = 1. - (( i - an )/bn )**2
      wss = wss + window(i)**2
      sig(i) = 1.0
      trend(i) = 0.0
    end do 

!++ 4. Applying Welch window to data

    do i=1,ndata
      data1(i) = data1(i)*window(i)
    end do

!++ 5. Forward Fourier transform

    call fft1df(ndata,data1,coef)
    call wavnum(ndata,dx,wavk2)

    do i=1,ndata/2+1
      wavk(i) = wavk2(i)
    end do

!++ 6. Calculate power spectral density

    pwsd(1) = (dx/wss)*real(cdabs(coef(1))**2)
    do i=2,ndata/2+1
      pwsd(i) = (dx/wss)*real(cdabs(coef(i))**2 +    &
                              cdabs(coef(ndata-i+2))**2 )
    end do 

!++ 7. Evaluating PSD calculation

    sqr = 0.
    do i=1,ndata
      sqr = sqr + (abs(data1(i))**2)*dx
    end do
    sqf = 0.
    do i=1,ndata
      sqf = sqf + (dx/ndata)*real(cdabs(coef(i))**2)
    end do

    if (eval .eq. 1) then
      print *,'SQR should be almost the same as SQF'
      print *,'SQR = ',sqr
      print *,'SQF = ',sqf
    end if

    return
    end subroutine psd1d 

!-----------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine psd2d(nx, ny, dx, dy, data, wavk, wavl, pwsd)

    implicit none

    integer,     intent(in)  :: nx, ny
    real,        intent(in)  :: dx, dy
    real,        intent(in)  :: data(nx,ny)
    real,        intent(out) :: wavk(nx/2+1), wavl(ny/2*2+1)
    real,        intent(out) :: pwsd(nx/2+1,ny/2*2+1)

    integer                  :: i, j
    real                     :: sqr, sqf
    real                     :: wavk2(nx), wavl2(ny), data1(nx,ny)

    double precision                    :: pi
    double complex                      :: ci
    double complex, dimension(nx,ny)    :: coef

    pi = 4.d0*datan(1.d0)
    ci = (0.d0,1.d0)


    data1(:,:) = data(:,:)

!++ FFT2D

    call fft2df(nx,ny,data1,coef)
    call wavnum(nx,dx,wavk2)
    call wavnum(ny,dy,wavl2)

    do i=1, nx/2+1
      wavk(i) = wavk2(i)
    end do
    do j=ny/2+1, ny/2*2+1
      wavl(j) = wavl2(j-ny/2)
    end do
    do j=1, ny/2
      wavl(j) = (-1.) * wavl2(ny/2+2-j)
    end do

!++ PSD ( & sorting )

    do j=1, ny/2
      pwsd(1,j) = (dx*dy/nx/ny)*real( cdabs(coef(1,(ny+1)/2+j))**2 )
      do i=2,nx/2+1
        pwsd(i,j) = (dx*dy/nx/ny)*real( cdabs(coef(i,(ny+1)/2+j))**2 &
                            + cdabs(coef(nx-i+2,ny/2+2-j))**2 )
      end do
    end do

    pwsd(1,ny/2+1) = (dx*dy/nx/ny)*real( cdabs(coef(1,1))**2 )
    do i=2,nx/2+1
      pwsd(i,ny/2+1) = (dx*dy/nx/ny)*real( cdabs(coef(i,1))**2 &
                            + cdabs(coef(nx-i+2,1))**2 )
    end do 

    do j=ny/2+2, ny/2*2+1
      pwsd(1,j) = (dx*dy/nx/ny)*real( cdabs(coef(1,j-ny/2))**2 )
      do i=2,nx/2+1
        pwsd(i,j) = (dx*dy/nx/ny)*real( cdabs(coef(i,j-ny/2))**2 &
                            + cdabs(coef(nx-i+2,ny+ny/2+2-j))**2 ) 
      end do 
    end do

!++ Evaluating

    sqr = 0.
    sqf = 0.
    do i=1, nx
     do j=1, ny
      sqr = sqr + (abs(data1(i,j))**2)*dx*dy
      sqf = sqf + (dx*dy/nx/ny)*real(cdabs(coef(i,j))**2)
     end do
    end do

!    print *,'SQR should be almost the same as SQF'
!    print *,'SQR = ',sqr
!    print *,'SQF = ',sqf

    return
    end subroutine psd2d


    end module pwrspd
