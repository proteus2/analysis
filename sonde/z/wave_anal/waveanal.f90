    module wave_analysis 
!----------------------------------------------------------------------------------
!
!   MODULE wave_analysis 
!
!   Calculate characteristics parameters of inertia-gravity waves
!   
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
!----------------------------------------------------------------------------------
!
!   Input
!
!   nz      : the number of vertical grid point 
!   uprt    : zonal wind perturbation
!   vprt    : meridional wind perturbation
!   tnhil   : hilbert tranformed T/Tbas(temperature perturbation/basic temperature)
!   
!   Ouput
! 
!   deg and degn  : propagation direction (deg)
!
!----------------------------------------------------------------------------------
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
!----------------------------------------------------------------------------------
!
!   Input
!
!   nz      : the number of vertical grid point
!   dz      : vertical grid spacing 
!   uprt    : zonal wind perturbation
!   vprt    : meridional wind perturbation
!   ubas    : basic-state zonal wind 
!   vbas    : basic-state meridional wind
!   deg     : propagation direction
!   avgn    : verticallly averaged stability
!
!   Ouput
! 
!   axr     : axial ration (intrinsic frequency/f) 
!   df      : degree of polarization
!----------------------------------------------------------------------------------
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

    do k=1,nz/2+1
      it(k)=real(ut(k)*conjg(ut(k))+vt(k)*conjg(vt(k)))
      dt(k)=real(ut(k)*conjg(ut(k))-vt(k)*conjg(vt(k)))
      pt(k)=2.0*real(conjg(ut(k))*vt(k))
      qt(k)=2.0*real(aimag(conjg(ut(K))*vt(k)))
      qts=qts+qt(k)
      dts=dts+dt(k)
      its=its+it(k)
      pts=pts+pt(k)
    enddo

    df  = sqrt( qts**2 + pts**2 + dts**2 )/its
    xi  = asin( abs( qts/(df*its) ) )/2.0
    axr = 1.0/tan(xi)

    do k=1,nz
      ubart(k)=- ubas(k)*sin(deg*pi/180.0) + vbas(k)*cos(deg*pi/180.0)
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

    axr = 1/(1/axr-axrcor)
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
                     cz  ,ci  ,cix ,ciy ,cgx ,cgy ,cgix,cgiy,wgr ,cgr ,  &
                     cx  ,cy  ,cpx ,cpy ,ubmm,dgmm,cgz, avgn)  

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
    real,                intent(in)  :: f   ,axr ,deg, avgn
    real,                intent(in)  :: kba ,mba   
    real, dimension(nz), intent(in)  :: ubas,vbas
    real,                intent(out) :: cz  ,ci  ,cix ,ciy
    real,                intent(out) :: cgx ,cgy ,cgix,cgiy,wgr ,cgr, cgz
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

    cgix = cix*(1.-(1./axr)**2)
    cgiy = ciy*(1.-(1./axr)**2)

    wgr = f*axr + kba*ubmm*cos((dgmm-deg)*pi/180.0)
    cgr = ci + ubmm*cos((dgmm-deg)*pi/180.0) 

    cgz = (avgn**2)*(kba**2)/(f*axr*mba**3) 
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
!   SUBROUTINE meanmk 
!
!   Calculate mean wave numbers (i.e., m and k) using the dispersion relation
!   of inertia gravity waves
! 
!----------------------------------------------------------------------------------
    subroutine meanmk(nz,dz,f,wf,nbas,ustprt,vstprt,mbar  ,mbaran,kbar   ,kbaran)
!----------------------------------------------------------------------------------
!
!   Input
!
!   nz      : the number of vertical grid point 
!   dz      : vertical grid spacing 
!   f       : Coriolis parameter
!   wf      : Axial ratio (omega/f)
!   nbas    : mean stability
!   ustprt  : zonal wind perturbation
!   vstprt  : meridional wind perturbation
!
!
!   Output
!  
!   mbar    : mean vertical wavenumber ( = 1/vertical wavelength) 
!   mbaran  : mean vertical wavenumber ( = 2*pi/vertical wavelength)
!   kbar    : mean horizontal wavenumber ( = 1/horizontal wavelength)
!   kbaran  : mean horizontal wavenumber ( = 2*pi/horizontal wavelength)
!
!----------------------------------------------------------------------------------

    implicit none

    integer,             intent(in)  :: nz
    real,                intent(in)  :: dz    ,f     ,wf     ,nbas
    real, dimension(nz), intent(in)  :: ustprt,vstprt
    real,                intent(out) :: mbar  ,mbaran,kbar   ,kbaran

    integer                          :: k
    real                             :: pi    ,s1    ,s2
    real                             :: umean ,vmean ,wss
    real, dimension(nz)              :: winfcn
    real, dimension(nz)              :: us    ,vs    ,usvs
    real, dimension(nz/2+1)          :: m
    double precision, dimension(4*nz+15) :: work
    double complex                   :: ci
    double complex, dimension(nz)    :: seq, coef

    pi = 4.0*atan(1.0)

    do k=1,nz/2+1
      m(k) = float(k-1)/(nz*dz)
    end do

    call welch(nz,winfcn,wss)

    umean = 0.0
    vmean = 0.0
    do k=1,nz
      umean = umean + ustprt(k)
      vmean = vmean + vstprt(k)
    end do
    umean = umean/float(nz)
    vmean = vmean/float(nz)

 
    do k=1,nz
      seq(k) = dble((ustprt(k)-umean)*winfcn(k)) + ci*0.d0
      coef(k) = seq(k)
    end do

    work(1:4*nz+15) = 0.d0
    call cffti(nz,work)
    call cfftf(nz,coef,work)

    do k=1,nz
      us(k) = real(cdabs(coef(k)*conjg(coef(k))))
    end do

    do k=1,nz
      seq(k) = dble((vstprt(k)-vmean)*winfcn(k)) + ci*0.d0
      coef(k) = seq(k)
    end do

    work(1:4*nz+15) = 0.d0
    call cffti(nz,work)
    call cfftf(nz,coef,work)

    do k=1,nz
      vs(k) = real(cdabs(coef(k)*conjg(coef(k))))
    end do

    do k=1,nz
      usvs(k) = (us(k)+vs(k))*dz/wss
    end do

    s1 = 0.0
    s2 = 0.0
    do k=2,nz/2
      s1 = s1 + usvs(k)*m(k)
      s2 = s2 + usvs(k)
    end do

    mbar   = s1/s2
    mbaran = 2.0*pi*mbar
    kbaran = f*mbaran*sqrt(wf**2-1)/nbas
    kbar   = kbaran/(2.0*pi)

    return

    end subroutine meanmk    


    end module wave_analysis 
