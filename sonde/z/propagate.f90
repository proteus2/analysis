!2345
    module propagate 
!----------------------------------------------------------------------------------
!
!   MODULE PROPAGATE
!
!   Calculate wave propagation characteristics
! 
!----------------------------------------------------------------------------------

    use window
    
    contains

!----------------------------------------------------------------------------------
!
!   SUBROUTINE ROTARY
!
!   Calculate fractions of upward and downward propagating waves 
! 
!----------------------------------------------------------------------------------

    subroutine rotary(nz,f,dz,uprt,vprt,wavm,rotpsd,cclock,clock,upfr,dnfr)  
   
    integer,             intent(in)  :: nz
    real,                intent(in)  :: dz     ,f
    real, dimension(nz), intent(in)  :: uprt   ,vprt
    real, dimension(nz), intent(out) :: wavm   ,rotpsd
    real,                intent(out) :: cclock ,clock
    real,                intent(out) :: upfr   ,dnfr

    integer             :: k
    real                :: sum   ,wss
    real, dimension(nz) :: winfcn,uprt1  ,vprt1
    double precision, dimension(4*nz+15) :: work
    double complex                       :: ci
    double complex,   dimension(nz)      :: seq   ,coef

    ci = (0.d0,1.d0)

    call welch(nz,winfcn,wss) 

    sum = 0.0
    do k=1,nz
      sum = sum + uprt(k)
    end do

    do k=1,nz
      uprt1(k) = uprt(k) - sum/float(nz)
    end do

    sum = 0.0
    do k=1,nz
      sum = sum + vprt(k)
    end do

    do k=1,nz
      vprt1(k) = vprt(k) - sum/float(nz)
    end do

    do k=1,nz
      seq(k) = dble(uprt1(k)*winfcn(k)) + ci*dble(vprt1(k)*winfcn(k))
      coef(k) = seq(k)
    end do

!----------------------------------------------------------------------------------
!   Construct rotary spectra
!----------------------------------------------------------------------------------

    work(1:4*nz+15) = 0.d0
    call cffti(nz,work)
    call cfftf(nz,coef,work)

!----------------------------------------------------------------------------------
!   Negative vertical wavenumber region
!   This region represents counter-clockwisely rotating hodograph.
!   In northern hemisphere, this part represents upward wave energy transport. 
!----------------------------------------------------------------------------------

    do k=1,nz/2
      rotpsd(k) = (dz/wss)*real(cdabs(coef(k+nz/2+1))**2)
    end do

!----------------------------------------------------------------------------------
!   Zero vertical wavenumber
!----------------------------------------------------------------------------------

    rotpsd(nz/2+1) = (dz/wss)*real(cdabs(coef(1))**2)

!----------------------------------------------------------------------------------
!   Positive vertical wavenumber region
!   This region represents clockwisely rotating hodograph.
!   In northern hemisphere, this part represents downward wave energy transport. 
!----------------------------------------------------------------------------------

    do k=nz/2+2,nz
      rotpsd(k) = (dz/wss)*real(cdabs(coef(k-nz/2))**2)
    end do

!----------------------------------------------------------------------------------
!   Wavenumber
!----------------------------------------------------------------------------------

    do k=1,nz/2
      wavm(k) = (-1)*float(nz/2-k+1)/(dz*nz)
    end do

    wavm(nz/2+1) = 0.0

    do k=nz/2+2,nz
      wavm(k) = float(k-nz/2-1)/(dz*nz)
    end do

!----------------------------------------------------------------------------------
!   Calculate the fraction of counter-clockwise and clockwise to total
!----------------------------------------------------------------------------------

    cclock = 0.0
    do k=1,nz/2-1
      cclock = cclock + 0.5*(rotpsd(k)+rotpsd(k+1))*(1./(dz*nz))
    end do

    clock = 0.0
    do k=nz/2+1,nz-1
      clock = clock + 0.5*(rotpsd(k)+rotpsd(k+1))*(1./(dz*nz))
    end do

    if ( f > 0.0 ) then
      upfr = cclock/(cclock+clock)
      dnfr =  clock/(cclock+clock) 
    else
      upfr =  clock/(cclock+clock)
      dnfr = cclock/(cclock+clock) 
    end if

    return

    end subroutine rotary

    end module propagate
