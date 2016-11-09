!2345
    module hilbert

!----------------------------------------------------------------------------------
!
!   This module requires double precision version of CFFTPACK
!
!----------------------------------------------------------------------------------

    contains

!----------------------------------------------------------------------------------
!
!   SUBROUTINE INVHIL 
!
!   This subroutine makes data whose phase is shifted by +90 degree from original
!   data.
!
!----------------------------------------------------------------------------------

    subroutine invhil(n     ,dz     ,var    ,var90 )

    implicit none

    integer,            intent(in)  :: n
    real,               intent(in)  :: dz
    real, dimension(n), intent(in)  :: var
    real, dimension(n), intent(out) :: var90

    integer                             :: k
    double complex                      :: ci 
    double complex, dimension(n)        :: seq   ,coef
    double precision                    :: sqr   ,sqf
    double precision, dimension(4*n+15) :: work


    ci = (0.d0,1.d0)

    do k=1,n
      seq(k)  = dble(var(k)) + ci*0.d0
      coef(k) = seq(k) 
    end do

!----------------------------------------------------------------------------------
!   Fourier transformation
!----------------------------------------------------------------------------------

    work(1:4*n+15) = 0.d0
    call cffti(n    ,work   )
    call cfftf(n    ,coef   ,work  )
 
!----------------------------------------------------------------------------------
!   Variance check
!----------------------------------------------------------------------------------

    sqr = 0.d0
    sqf = 0.d0 

    do k=1,n
      sqr = sqr + (cdabs(seq(k))**2)*dz
    end do

    do k=1,n
      sqf = sqf + dz*(cdabs(coef(k))**2)/dble(n)
    end do

    write(6,'(a,e17.10)') '(INVHIL): TOTAL VARIANCE OF DATA               : ',sqr
    write(6,'(a,e17.10)') '(INVHIL): TOTAL VARIANCE OF FOURI COEFFICIENTS : ',sqf

!----------------------------------------------------------------------------------
!   Hilbert transform : h(f) ==> i*sgn(f)*h(f)
!   Inverse Hilbert transform in positive wavenumber (frequency)
!   h(f) ==> (-i)*sgn(f)*h(f) = (-i)*h(f)
!----------------------------------------------------------------------------------

    do k=2,n/2+1
      coef(k)=(-1)*ci*coef(k)
    end do

!----------------------------------------------------------------------------------
!   Hilbert transform : h(f) ==> i*sgn(f)*h(f)
!   Inverse Hilbert transform in negative wavenumber (frequency)
!   h(f) ==> ( i)*sgn(f)*h(f) = (-i)*h(f)
!----------------------------------------------------------------------------------

    do k=n/2+2,n
      coef(k)=ci*coef(k)
    end do

!----------------------------------------------------------------------------------
!   Inverse Fourier transformation
!----------------------------------------------------------------------------------

    work(1:4*n+15) = 0.d0

    call cffti(n    ,work  )
    call cfftb(n    ,coef  ,work  )

!----------------------------------------------------------------------------------
!   Extract real part from the inverse Fourier transformed coefficients 
!----------------------------------------------------------------------------------

    do k=1,n
      var90(k)=real(coef(k))/float(n)
    end do

    return
    
    end subroutine invhil

    end module hilbert
