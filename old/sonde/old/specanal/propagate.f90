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
!   SUBROUTINE PROPDIR
!
!   Calculate wave propagation direction 
! 
!----------------------------------------------------------------------------------

    subroutine propdir(nz   ,uprt  ,vprt  ,tnhil  ,phi   ,deg   ,degn  )

    implicit none

    integer, intent(in)              :: nz
    real, dimension(nz), intent(in)  :: uprt  ,vprt  ,tnhil
    real,                intent(out) :: phi   ,deg   ,degn 

    integer                          :: k
    real                             :: x     ,y     ,pi
    real, dimension(nz)              :: ut90  ,vt90
   
    pi = 4.0*atan(1.0)
 
    do k=1,nz
      ut90(k) = uprt(k)*tnhil(k)
      vt90(k) = vprt(k)*tnhil(k)
    end do
    
    x = 0.0
    y = 0.0
    do k=1,nz
      x=x+ut90(k)
      y=y+vt90(k)
    end do
    x=x/float(nz)
    y=y/float(nz)

    if ((abs(x) < 1.e-10).and.(x >= 0.0)) then
      x=1.0e-10
    else if ((abs(x) < 1.e-10).and.(x < 0.0)) then
      x=-1.0e-10
    end if

    if ((abs(y) < 1.e-10).and.(y >= 0.0)) then
      y=1.0e-10
    else if ((abs(y) < 1.e-10).and.(y < 0.0)) then
      y=-1.0e-10
    end if

    phi=atan(y/x)

!----------------------------------------------------------------------------------
!   deg measures angle counterclockwisely from positive x axis.
!   Therefore, deg is different from the angle measuring wind direction
!----------------------------------------------------------------------------------

    if ((x >= 0.0).and.(y >= 0.0)) then
      deg=phi*180.0/pi
    else if ((x < 0.0).and.(y >= 0.0)) then
      deg=180.0+phi*180.0/pi
    else if ((x < 0.0).and.(y < 0.0)) then
      deg=180.0+phi*180.0/pi
    else
      deg=360.0+phi*180.0/pi
    end if  

!----------------------------------------------------------------------------------
!   degn measures the angle in the same way as measuring wind direction
!----------------------------------------------------------------------------------

    degn = 270.0 - deg
    if (degn < 0.0) then
      degn = 360.0 + degn
    end if

    return
  
    end subroutine propdir

!----------------------------------------------------------------------------------
!
!   SUBROUTINE STOKE
!
!   Calculate axial ratio of wind perturbation hodograph using Stoke's parameter
!   methods.
! 
!----------------------------------------------------------------------------------

    subroutine stoke(nz   ,dz    ,uprt  ,vprt  ,ubas  ,vbas  ,deg   ,avgn   ,axr   ,df    )

    integer,             intent(in)  :: nz
    real,                intent(in)  :: dz    ,deg   ,avgn
    real, dimension(nz), intent(in)  :: uprt  ,vprt  ,ubas  ,vbas
    real,                intent(out) :: axr   ,df

    integer                 :: k

    real                    :: its   ,dts   ,pts   ,qts
    real                    :: xi    ,axrcor,pi
    real                    :: ubartb,ubartt,mshear
    real, dimension(nz/2+1) :: it    ,dt    ,pt    ,qt
    real, dimension(nz)     :: ubart ,ubart_dz
 
    double precision                     :: sqr   ,sqf
    double precision, dimension(4*nz+15) :: work
    double complex                       :: ci
    double complex, dimension(nz)        :: ut    ,vt

    pi = 4.0*atan(1.0)
    ci = (0.d0,1.d0)

    do k=1,nz
      ut(k) = dble(uprt(k)) + ci*0.0
    end do

    do k=1,nz
      vt(k) = dble(vprt(k)) + ci*0.0
    end do

    work(1:4*nz+15) = 0.d0
    call cffti(nz,work)
    call cfftf(nz,ut,work)

    work(1:4*nz+15) = 0.d0
    call cffti(nz,work)
    call cfftf(nz,vt,work)

    qts=0.0
    dts=0.0
    its=0.0
    pts=0.0

    do k=2,nz/2+1
      it(k)=real(ut(k)*conjg(ut(k))+vt(k)*conjg(vt(k)))
      dt(k)=real(ut(k)*conjg(ut(k))-vt(k)*conjg(vt(k)))
      pt(k)=2.0*real(conjg(ut(k))*vt(k))
      qt(k)=real(2.0*aimag(conjg(ut(k))*vt(k)))
      qts=qts+qt(k)
      dts=dts+dt(k)
      its=its+it(k)
      pts=pts+pt(k)
    enddo

    df  = sqrt( qts**2 + pts**2 + dts**2 )/its
    xi  = asin( abs( qts/(df*its) ) )/2.0
    axr = 1.0/tan(xi)

    do k=1,nz
      ubart(k)=ubas(k)*sin(deg*pi/180.0) - vbas(k)*cos(deg*pi/180.0)
    enddo
    ubartb = 2.0*ubart(1) - ubart(2)
    ubartt = 2.0*ubart(nz) - ubart(nz-1)

    ubart_dz(1)=(ubart(2)-ubartb)/(2.*dz)
    ubart_dz(nz)=(ubartt-ubart(nz-1))/(2.*dz)

    do k=2,nz-1
      ubart_dz(k)=(ubart(k+1)-ubart(k-1))/(2.*dz)
    end do
  
    mshear = 0.0
    do k=1,nz
      mshear = mshear + ubart_dz(k)
    end do
    mshear = mshear/float(nz)

    axrcor = mshear/avgn
  
    write(8,*) 'shear', axrcor  

!    axr = axr - axrcor
    axr = 1./(1./axr-axrcor)

    return

    end subroutine stoke
 
!----------------------------------------------------------------------------------
!
!   SUBROUTINE PHASE
!
!   Calculate phase and group velocities 
! 
!----------------------------------------------------------------------------------

    subroutine phase(nz  ,f   ,axr ,kba ,mba ,deg ,ubas,vbas,  &
                     cz  ,ci  ,cix ,ciy ,cgx ,cgy ,wgr ,cgr ,  &
                     cx  ,cy  ,cpx ,cpy ,ubmm,dgmm)  

!----------------------------------------------------------------------------------
!
!   Input
!
!   f      : Coriolis parameter
!   axr    : Axial ratio (omega/f)
!   deg    : Wave propagation direction
!   kba    : Mean horizontal wavenumber (angular wavenumber)
!   mba    : Mean vertical wavenumber (angular wavenumber)
!   ubas   : Basic-state zonal wind
!   vbas   : Basic-state meridional wind
!
!   Output
!
!   cz     : Vertical phase speed
!   ci     : Intrinsic horizontal phase speed in the direction of propagation
!   cix    : Intrinsic zonal phase speed   
!   ciy    : Intrinsic meridional phase speed
!   cgx    : Zonal group velocity
!   cgy    : Meridional group velocity
!   wgr    : Ground-relative wave frequency
!   cgr    : Ground-relative wave phase speed in the direction of wawe propagation
!   cx     : Ground-relative zonal phase speed (not vector component)
!   cy     : Ground-relative meridional phase speed (not vector component)
!   cpx    : Ground-relative zonal phase speed (vector component)
!   cpy    : Ground-relative meridional phase speed (vector component)
!   ubmm   : Vertically averaged basic-state wind speed
!   dgmm   : Vertically averaged basic-state wind direction
!
!----------------------------------------------------------------------------------

    integer,             intent(in)  :: nz
    real,                intent(in)  :: f   ,axr ,deg
    real,                intent(in)  :: kba ,mba   
    real, dimension(nz), intent(in)  :: ubas,vbas
    real,                intent(out) :: cz  ,ci  ,cix ,ciy
    real,                intent(out) :: cgx ,cgy ,wgr ,cgr
    real,                intent(out) :: cx  ,cy  ,cpx ,cpy

    real                             :: x   ,y
    real                             :: umm ,vmm
    real                             :: ubmm,thmm,dgmm

    pi = 4.0*atan(1.0)

    cz = f*axr/mba
    ci = f*axr/kba
    cix = ci*cos(pi*deg/180.0)
    ciy = ci*sin(pi*deg/180.0)

    x = 0 
    y = 0
    do k=1,nz
      x = x + ubas(k)
      y = y + vbas(k)
    end do
    umm = x/float(nz)
    vmm = y/float(nz)

    ubmm = sqrt(umm**2 + vmm**2)
    thmm = atan(vmm/umm)

    if ((umm >= 0.0).and.(vmm >= 0.0)) then
      dgmm = thmm*180.0/pi
    else if ((umm < 0.0).and.(vmm >= 0.0)) then
      dgmm = 180.0+thmm*180.0/pi
    else if ((umm < 0.0).and.(vmm < 0.0)) then
      dgmm = 180.0+thmm*180.0/pi
    else
      dgmm = 360.0+thmm*180.0/pi
    end if

!----------------------------------------------------------------------------------
!   Caution! Horizontal group velocity calculation in this program is just
!   approximately valid (not exact calculation). Actually this approximation
!   is exact result when the assumption of hydrostatic internal gravity 
!   wave is valid
!----------------------------------------------------------------------------------
!
!   cgx = cix + umm
!   cgy = ciy + vmm
!
!----------------------------------------------------------------------------------
!   Exact group velocity
!----------------------------------------------------------------------------------

    cgx = cix*(1.-(1./axr)**2) + umm
    cgy = ciy*(1.-(1./axr)**2) + vmm

    wgr = f*axr + kba*ubmm*cos((dgmm-deg)*pi/180.0)
    cgr = ci + ubmm*cos((dgmm-deg)*pi/180.0) 
 
!----------------------------------------------------------------------------------
!   Phase speed relative to the ground is not vector. So, this program divide
!   absolute value of c by sine or cosine of phase propagation direction to
!   calculate phase speed component. However, you should use cp_x and cp_y
!   in order to represent the phase speed relative to the ground as vector
!----------------------------------------------------------------------------------
 
    cx  = cgr/cos(deg*pi/180.0)
    cy  = cgr/sin(deg*pi/180.0)
    cpx = cgr*cos(deg*pi/180.0)
    cpy = cgr*sin(deg*pi/180.0)

    return
 
    end subroutine phase 

!----------------------------------------------------------------------------------
!
!   SUBROUTINE PROPHIST
!
!   Calculate energy weighted wave propagation histogram 
! 
!----------------------------------------------------------------------------------

    subroutine prophist(ndat  ,avget ,deg   ,angs  ,anghist,meandir,msgv) 

    implicit none

    integer,               intent(in)  :: ndat
    real, dimension(ndat), intent(in)  :: avget
    real, dimension(ndat), intent(in)  :: deg
    real,                  intent(in)  :: msgv
    real,                  intent(out) :: meandir
    real, dimension(12),   intent(out) :: angs  ,anghist

    integer :: k    ,it
    real    :: x    ,y    ,pi    ,sumet
    

    pi = 4.0*atan(1.0)

    do k=1,12
      angs(k) = float(k)*30. - 1.0
    end do

    sumet = 0.0
    do it=1,ndat
      if ( avget(it) /= msgv ) sumet = sumet + avget(it)
    end do

    do it=1,ndat
      if ( deg(it) /= msgv ) then
        if      ( deg(it) >=   0.0 .and. deg(it) <  30.0 ) then
          anghist( 1) = anghist( 1) + avget(it)/sumet
        else if ( deg(it) >=  30.0 .and. deg(it) <  60.0 ) then
          anghist( 2) = anghist( 2) + avget(it)/sumet
        else if ( deg(it) >=  60.0 .and. deg(it) <  90.0 ) then
          anghist( 3) = anghist( 3) + avget(it)/sumet
        else if ( deg(it) >=  90.0 .and. deg(it) < 120.0 ) then
          anghist( 4) = anghist( 4) + avget(it)/sumet
        else if ( deg(it) >= 120.0 .and. deg(it) < 150.0 ) then
          anghist( 5) = anghist( 5) + avget(it)/sumet
        else if ( deg(it) >= 150.0 .and. deg(it) < 180.0 ) then
          anghist( 6) = anghist( 6) + avget(it)/sumet
        else if ( deg(it) >= 180.0 .and. deg(it) < 210.0 ) then
          anghist( 7) = anghist( 7) + avget(it)/sumet
        else if ( deg(it) >= 210.0 .and. deg(it) < 240.0 ) then
          anghist( 8) = anghist( 8) + avget(it)/sumet
        else if ( deg(it) >= 240.0 .and. deg(it) < 270.0 ) then
          anghist( 9) = anghist( 9) + avget(it)/sumet
        else if ( deg(it) >= 270.0 .and. deg(it) < 300.0 ) then
          anghist(10) = anghist(10) + avget(it)/sumet
        else if ( deg(it) >= 300.0 .and. deg(it) < 330.0 ) then
          anghist(11) = anghist(11) + avget(it)/sumet
        else if ( deg(it) >= 330.0 .and. deg(it) < 360.0 ) then
          anghist(12) = anghist(12) + avget(it)/sumet
        end if
      end if
    enddo

    x = 0.0
    y = 0.0
    do it=1,ndat
      if ( deg(it) /= msgv ) then
        x = x + avget(it)*cos(deg(it)*pi/180.0)
        y = y + avget(it)*sin(deg(it)*pi/180.0)
      end if
    end do
    meandir = atan(y/x)

    if      ( x >= 0.0 .and. y >= 0.0 ) then
      meandir = meandir*180.0/pi
    else if ( x <  0.0 .and. y >= 0.0 ) then
      meandir = 180.0 + meandir*180.0/pi
    else if ( x <  0.0 .and. y <  0.0 ) then
      meandir = 180.0 + meandir*180.0/pi
    else
      meandir = 360.0 + meandir*180.0/pi
    end if

    return
 
    end subroutine prophist

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
