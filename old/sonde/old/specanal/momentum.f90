!2345
    module momentum

!-------------------------------------------------------------------------------
!
!   Vincent and Alexander (2000)
!
!   Measurements of w' are not possible in rawinsonde observation, but
!   indirect estimates of u'w' are possible using Hilbert transform.
!   From the polarization relations between u', T', and w' and
!   assuming a monochromatic wave of intrinsic frequency \omega, it is
!   found that
!                  ____     \omega g  ____
!                  u'w'  =  --------- u'T'_{+90},
!                              N^2
!
!   where T'_{+90} = i T' and is obtained by first Hilbert transforming
!   T' in the appropriate sense and then reverse transformed to produce
!   90 degree phase shift in T' at all spatial frequencies without
!   changing the amplitude (Vincent et al. 1997). The vertical flux of
!   meridional momentum v'w' is found in a similar manner. In practice,
!   the quantity u'w' \delta_(\omega) is required, where \delta_(\omega)
!   = ( 1 - f^2/(\omega)^2 ) (Fritts and Vincent 1987), which shows that
!   low-frequency waves are less effective at transporting momentum.
!
!   References:
!
!   Fritts, D. C., and R. A. Vincent, 1987: Mesospheric momentum flux
!      studies at Adelaide. Australia: Observations and a gravity wave-
!      tidal interaction model. J. Atmos. Sci., 44, 605-619.
!   Vincent, R. A., S. J. Allen, and S. D. Eckermann, 1997: Gravity-
!      wave parameters in the lower stratosphere, Gravity Wave Processes:
!      Their Parameterization in Global Climate Models. K. Hamilton ed.
!      7-25.
!   Vincent, R. A., and M. J. Alexander, 2000: Gravity waves in the tro-
!      pical lower stratosphere: An observational study of seasonal and
!      interannual variability. J. Geophys. Res., 105, 17971-17982.
!
!-------------------------------------------------------------------------------

    contains

    subroutine uwvw(nz,avgn,axr,f,uprt,vprt,tnhil,rhobas,rhouw,rhovw,uw,vw,  &
                    mruw,mrvw,muw,mvw)
  
    integer,             intent(in)  :: nz
    real,                intent(in)  :: avgn ,axr  ,f
    real, dimension(nz), intent(in)  :: uprt ,vprt ,tnhil,rhobas
    real, dimension(nz), intent(out) :: rhouw,rhovw,uw   ,vw
    real,                intent(out) :: mruw ,mrvw ,muw  ,mvw

    integer             :: k
    real                :: omegabar,grav
    real, dimension(nz) :: ut90    ,vt90

    grav     = 9.806
    omegabar = axr*f

    do k=1,nz
      x        = omegabar*grav/(avgn**2)
      deltam   = 1.0 - (f/omegabar)**2
      ut90(k)  = uprt(k)*tnhil(k)
      vt90(k)  = vprt(k)*tnhil(k)
      uw(k)    = x*ut90(k)*deltam
      vw(k)    = x*vt90(k)*deltam
      rhouw(k) = rhobas(k)*uw(k)
      rhovw(k) = rhobas(k)*vw(k)
    end do

    mruw = 0.0
    mrvw = 0.0
    muw  = 0.0
    mvw  = 0.0
    do k=1,nz
      mruw = mruw + rhouw(k)
      mrvw = mrvw + rhovw(k)
      muw  = muw  + uw(k)
      mvw  = mvw  + vw(k)
    end do 
    mruw = mruw/float(nz)
    mrvw = mrvw/float(nz)
    muw  = muw /float(nz)
    mvw  = mvw /float(nz)  

    return
    
    end subroutine uwvw

    end module momentum
