!2345
    module window 
!------------------------------------------------------------------------------
!
!   Define window function to prohibit energy leakage in FFT.    
!
!------------------------------------------------------------------------------

    contains

!------------------------------------------------------------------------------
!   Welch window
!------------------------------------------------------------------------------

    subroutine welch(n,windowfcn,wss)

    implicit none

    integer,            intent(in)  :: n
    real, dimension(n), intent(out) :: windowfcn
    real,               intent(out) :: wss

    integer             :: k
    real                :: an    ,bn

    an = ( float(n) - 1.0 )/2.
    bn = ( float(n) + 1.0 )/2.
    
    do k=1,n  
      windowfcn(k) = 1. - (( k - an )/bn )**2
    end do

    wss = 0.0
    do k=1,n
      wss = wss + windowfcn(k)**2
    end do

    return
    
    end subroutine welch

!------------------------------------------------------------------------------
!   Bartlett window
!------------------------------------------------------------------------------

    subroutine bartlett(n,windowfcn,wss)

    implicit none

    integer,            intent(in)  :: n
    real, dimension(n), intent(out) :: windowfcn
    real,               intent(out) :: wss

    integer             :: k

    do k=1,n
      windowfcn(k)=1-abs(float(k-int(n*.5)))/(int(n*.5))
    end do

    wss = 0.0
    do k=1,n
      wss = wss + windowfcn(k)**2
    end do

    return

    end subroutine bartlett

!------------------------------------------------------------------------------
!   Taper band
!------------------------------------------------------------------------------

    subroutine taper(n,percent,windowfcn,wss)

    implicit none

    integer,            intent(in)  :: n    ,percent
    real, dimension(n), intent(out) :: windowfcn
    real,               intent(out) :: wss

    integer             :: k    ,dk
    real                :: pi

    pi = 4.0*atan(1.0)
    dk = int(n*(percent/100.))

    do k=1,dk
      windowfcn(k) = sin(5*pi*(k-1)/n)**2
    end do
    do k=n,n-dk+1,-1
      windowfcn(k) = sin(5*pi*(k-n)/n)**2
    end do
    do k=dk+1,n-dk
      windowfcn(k) = 1.0
    end do

    wss = 0.0
    do k=1,n
      wss = wss + windowfcn(k)**2
    end do

    return
    
    end subroutine taper

    end module window
