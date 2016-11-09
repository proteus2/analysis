    
    module cal_xsq

    use params
    use quadpack

    contains

!===============================================================================
!
!   Subroutine inti 
!
!   initialization 
!
!===============================================================================

    subroutine inti
   
    implicit none

    integer     :: l

    do l= 1,nc
      cpgrd(l) = -60.+ dble(l-1)*dc
    end do 

    end subroutine inti

!===============================================================================
!
!   Subroutine xsqun
!
!   Calculates the term |X|^2 for the basic-state wind without vertical shear.
!
!===============================================================================

    subroutine xsqun(ut    ,x2out )

    implicit none

    double precision, intent(in)                 :: ut
    double precision, intent(out), dimension(nc) :: x2out
    double precision, external                   :: xsqu

    integer                 :: l
    logical                 :: check1,reglar
    double precision        :: cpa   ,cpb   
    double precision        :: xsq   ,abserr 
    integer                 :: ier 
!-------------------------------------------------------------------------------

    do l= 1,nc
      cpa = cpgrd(l) - dc*0.5
      cpb = cpgrd(l) + dc*0.5
      reglar = .true.
      check1 =  cpa <= ut .and. ut <= cpb 
      if ( check1 ) reglar = .false.
      if ( reglar ) then
        call ninteg1d ( xsqu , cpa, cpb, epsabs, epsrel, key, xsq  , abserr, ier )
        if (ier /= 0) print*, abserr,ier,'error'
        x2out(l) = xsq/dc
      else
        x2out(l) = 0.0
      end if
    end do

    return
    end subroutine xsqun
  
!===============================================================================
!
!   Subroutine xsqsh
!
!   Calculates the term |X|^2 for the basic-state wind with vertical shear.
!
!===============================================================================

    subroutine xsqsh(ut    ,ub    ,u0    ,x2out )
                     
                     
    implicit none
    
    double precision, intent(in)                 :: ut    ,ub    ,u0
    double precision, intent(out), dimension(nc) :: x2out
    double precision, external                   :: xsqs

    integer                 :: l
    logical                 :: check1,check2,check3,reglar
    double precision        :: cpa   ,cpb   
    double precision        :: xsq   ,abserr 
    integer                 :: ier 
!-------------------------------------------------------------------------------

    do l=1,nc
      cpa = cpgrd(l) - dc*0.5
      cpb = cpgrd(l) + dc*0.5
      reglar = .true.
      check1 =  cpa <= ut .and. ut <= cpb 
      check2 =  cpa <= ub .and. ub <= cpb
      if ( check1 .or. check2 ) reglar = .false.
      if ( .not. reglar ) then
        x2out(l) = 0.0
      else
        call ninteg1d ( xsqs , cpa, cpb, epsabs, epsrel, key, xsq  , abserr, ier )
        x2out(l) = xsq/dc
      end if
    end do
    
    return
    end subroutine xsqsh

    end module cal_xsq

!-------------------------------------------------------------------------------------------
!
! external function
!
!-------------------------------------------------------------------------------------------

    function xsqu(cpx)
   
    use params

    implicit none

    double precision, intent(in) :: cpx
    double precision             :: xsqu

    double precision             :: gwcpi
    double precision             :: lm1   ,lm2   ,zd    ,zm
    double complex               :: ai    ,xu1   ,xu1c  ,xu2   ,xu2c  ,xp    ,xn    ,yu1   ,yu1c
    double complex               :: numer ,denom ,xxc

    gwcpi  = 4.d0*atan(1.d0)

    ai = (0.d0,1.d0)
    zm = (zt+zb)*0.5
    zd = (zt-zb)*0.5

    lm2   = nst/abs(ut-cpx)
    lm1   = ntr/abs(ut-cpx)
    xu1   = 2.d0*exp( ai*lm1*zb)
    xu1c  = 2.d0*exp(-ai*lm1*zb)
    xu2   = 2.d0*exp( ai*lm1*(zt-zb))
    xu2c  = 2.d0*exp(-ai*lm1*(zt-zb))
    xp    = 1.0+sign(1.0,ut-cpx)*nst/ntr
    xn    = 1.0-sign(1.0,ut-cpx)*nst/ntr
    yu1   = 2.0/(lm1*zd)**2 + (1.d0/(ai*lm1))*(2.0/zd)
    yu1c  = conjg(yu1)
    numer = xu1*xu2*(2.0*yu1-xu2c*yu1c)+xu1c*xu2c*(2.0*yu1c-xu2*yu1)
    denom = 2.0*(xn*xu1*xu2+xp*xu1c*xu2c)
    xxc   = numer/denom
    xsqu  = (abs(xxc)**2)

    return
    end

    function xsqs(cpx)
  
    use params

    implicit none

    double precision, intent(in) :: cpx
    double precision             :: xsqs

    double precision     :: lm1   ,lm2   ,ri    ,mu    ,zm    ,zd
    double precision     :: zalph ,alph
    double precision     :: ztqdzt,ztqdzb,ztqs
    double precision     :: ztsb  ,ztsdzb,ztss  ,ztsdzs
    double precision     :: ztus  ,ztudzs,ztut  ,ztudzt
    double complex       :: ai    ,xs0   ,xs2   ,xs2c  ,xs3   ,xs3c  ,xs4   ,xs4c
    double complex       :: xs5   ,xa    ,xac   ,ys1   ,ys1c  ,ys2   ,ys2c  ,ys3   ,ys3c
    double complex       :: xp    ,xn
    double complex       :: bfac  ,coef1 ,coef2 ,coef3 ,coef1c,coef2c,coef3c
    double complex       :: numer ,denom ,xxc
    double precision     :: gwcpi

    gwcpi  = 4.d0*atan(1.d0)
    ai   = (0.d0,1.d0)
    alph = (ut-u0)/zs
    ri   = ntr*ntr/(alph*alph)
    mu   = sqrt(ri-0.25)
    zm   = (zt+zb)*0.5
    zd   = (zt-zb)*0.5

    !--- Set branch point factor
    bfac = exp(-sign(1.0,alph)*ai*gwcpi)

    zalph = (cpx-u0)/alph
    lm1 = ntr/abs(ut-cpx)
    lm2 = nst/abs(ut-cpx)
    coef1 = (0.5-ai*mu)/(2.*ai*mu)
    coef2  = ai*lm1*(zs-zalph)-(0.5-ai*mu)
    coef3  = ai*lm1*(zs-zalph)-(0.5+ai*mu)
    coef1c = conjg(coef1)
    coef2c = conjg(coef2)
    coef3c = conjg(coef3)
    ztqdzt = -2.0/zd
    ztqdzb = 2.0/zd
    ztqs   = 1.0 - (zs-zm)**2/zd**2
    ztsb   = 2.0/(ri+2.0)*(zb-zalph)**2/zd**2
    ztsdzb = 4.0/(ri+2.0)*(zb-zalph)/zd**2
    ztss   = 2.0/(ri+2.0)*(zs-zalph)**2/zd**2
    ztsdzs = 4.0/(ri+2.0)*(zs-zalph)/zd**2
    ztus   = 2.0/(lm1*zd)**2
    ztudzs = 0.0
    ztut   = 2.0/(lm1*zd)**2
    ztudzt = 0.0
    if ( zb-zalph > 0.0 ) then
      ys1  = coef1 /(0.5-ai*mu)/(zb-zalph)*(abs(zb-zalph))**(0.5-ai*mu)*  &
             ( (0.5-ai*mu)*ztsb - (ztqdzb + ztsdzb)*(zb-zalph) )
      ys1c = coef1c/(0.5+ai*mu)/(zb-zalph)*(abs(zb-zalph))**(0.5+ai*mu)*  &
             ( (0.5+ai*mu)*ztsb - (ztqdzb + ztsdzb)*(zb-zalph) )
    else 
      ys1  = coef1 /(0.5-ai*mu)/(zb-zalph)*  &
             (abs(zb-zalph)*bfac)**(0.5-ai*mu)*  &
             ( (0.5-ai*mu)*ztsb - (ztqdzb + ztsdzb)*(zb-zalph) )
      ys1c = coef1c/(0.5+ai*mu)/(zb-zalph)*  &
             (abs(zb-zalph)*bfac)**(0.5+ai*mu)*  &
             ( (0.5+ai*mu)*ztsb - (ztqdzb + ztsdzb)*(zb-zalph) )
    end if          
    ys2  = (ztqs+ztss)-(ztsdzs-ztudzs)*(zs-zalph)-  &
           ai*lm1*(zs-zalph)*(ztss-ztus)   
    ys3  = ztut + (1.0/(ai*lm1))*(ztqdzt+ztudzt)  
    ys2c = conjg(ys2)
    ys3c = conjg(ys3)
    if ( -zalph > 0.0 ) then
      xs0 = (abs(-zalph))**(2.0*ai*mu)
    else
      xs0 = (abs(-zalph)*bfac)**(2.0*ai*mu)
    end if
    if ( zs-zalph > 0.0 ) then
      xs2  = coef2 *(abs(zs-zalph))**(0.5+ai*mu)
      xs3  = coef3 *(abs(zs-zalph))**(0.5-ai*mu)
      xs2c = coef2c*(abs(zs-zalph))**(0.5-ai*mu)
      xs3c = coef3c*(abs(zs-zalph))**(0.5+ai*mu)
    else
      xs2  = coef2 *(abs(zs-zalph)*bfac)**(0.5+ai*mu)
      xs3  = coef3 *(abs(zs-zalph)*bfac)**(0.5-ai*mu)
      xs2c = coef2c*(abs(zs-zalph)*bfac)**(0.5-ai*mu)
      xs3c = coef3c*(abs(zs-zalph)*bfac)**(0.5+ai*mu)
    end if
    xp    = 1.0+sign(1.0,ut-cpx)*nst/ntr
    xn    = 1.0-sign(1.0,ut-cpx)*nst/ntr
    xs4   = 2.0*exp(ai*lm1*(zt-zs))
    xs4c  = conjg(xs4)
    xa    = xs2 /coef2 *(real(ys2)-(0.5-ai*mu)*(ztss-ztus))+  &
            2.0*ai*mu*(zs-zalph)*ys1c
    xac   = xs2c/coef2c*(real(ys2)-(0.5+ai*mu)*(ztss-ztus))-  &
            2.0*ai*mu*(zs-zalph)*ys1
    xs5   = 4.0*(xa-xac*xs0)
    numer = xs4*ys3c*(xs2-xs3*xs0)+xs4c*ys3*(xs3c-xs2c*xs0)+xs5
    denom = xn*xs4*(xs2-xs3*xs0)+xp*xs4c*(xs3c-xs2c*xs0)
    xxc   = numer/denom
    xsqs  = (abs(xxc)**2)

    return
    end
