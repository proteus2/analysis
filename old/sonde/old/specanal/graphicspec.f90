!2345
    module graphicspec

    use varspec

    contains

    subroutine pltspec (ndat,year,month,msgv)

    implicit none

    integer, intent(in) :: ndat  ,year ,month
    real,    intent(in) :: msgv

    integer :: i,j ,nff
    character(len=40), dimension(9) :: xlab
    character(len=40), dimension(9) :: ylab
    character(len=40), dimension(9) :: xlab1
    character(len=40), dimension(9) :: ylab1
    character(len=3) , dimension(12) :: mon
    character(len=100) :: titlest,titletr
    character(len=100) :: nmdst  ,nmdtr
    character(len=100) :: nupfrst,nupfrtr
    real :: x1a,x2a,y1a,y2a
    real :: x1b,x2b,y1b,y2b
    real :: x1c,x2c,y1c,y2c
    real :: x1d,x2d,y1d,y2d
    real :: x1e,x2e,y1e,y2e
    real :: x1f,x2f,y1f,y2f

    mon = (/'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'/)

    write(titlest,'(a,i4)') 'Wave propagation in the stratosphere in '//mon(month)//', ',year
    write(titletr,'(a,i4)') 'Wave propagation in the troposphere in '//mon(month)//', ',year

    xlab = (/'-40','-30','-20','-10','0','10','20','30','40'/)
    ylab = (/'-40','-30','-20','-10','0','10','20','30','40'/)
    xlab1 = (/'-10','-7.5','-5','-2.5','0','2.5','5','7.5','10'/)
    ylab1 = (/'-10','-7.5','-5','-2.5','0','2.5','5','7.5','10'/)

    x1a = 0.15; x2a = 0.40; y1a = 0.69; y2a = 0.94
    x1b = 0.15; x2b = 0.40; y1b = 0.37; y2b = 0.62
    x1c = 0.15; x2c = 0.40; y1c = 0.05; y2c = 0.30
    x1d = 0.50; x2d = 0.75; y1d = 0.69; y2d = 0.94
    x1e = 0.50; x2e = 0.75; y1e = 0.37; y2e = 0.62
    x1f = 0.50; x2f = 0.75; y1f = 0.05; y2f = 0.30

    call GACWK(1)

    call GSTXFP(-4,2)

    call GASETR('XMJ.',0.006)
    call GASETR('YMJ.',0.006)
    call GASETR('XMN.',0.003)
    call GASETR('YMN.',0.003)
    call GASETI('LTY.',1)
    call PCSETI('FN.',21)

    call SET(0.0,1.0,0.0,1.0,0.0,1.0,0.0,1.0,1)
    call PLCHHQ(0.5,0.985,':F22:'//trim(titlest),0.012,0.0,0.0)

    call SET(x1a,x2a,y1a,y2a,-40.0,40.0,-40.0,40.0,1)
    call GSPLCI(1)
    call GSLWSC(2.0)
    do i=1,ndat 
      if ( cixst(i) /= msgv ) then
        call POINT(cixst(i),ciyst(i))
      end if
    end do
    call GSLWSC(1.0)
    call GRIDAL(8,4,8,4,0,0,5,0.0,0.0)
    call SET(x1a-0.05,x2a+0.05,y1a-0.05,y2a+0.05,x1a-0.05,x2a+0.05,y1a-0.05,y2a+0.05,1)
    do i=1,9
      call PLCHHQ(x1a+(i-1)*(x2a-x1a)/8.,y1a-0.01,trim(xlab(i)),0.008,0.0,0.0)
    end do
    do j=1,9
      call PLCHHQ(x1a-0.005,y1a+(j-1)*(y2a-y1a)/8.,trim(ylab(j)),0.008,0.0,1.0)
    end do
    call PLCHHQ((x1a+x2a)*0.5,y2a+0.015,':F22:Intrinsic phase speed'  ,0.01, 0.0,0.0)
    call PLCHHQ((x1a+x2a)*0.5,y1a-0.03,':F22:c:B:ix:N: (m s:S:-1:N:)',0.01,0.0,0.0)
    call PLCHHQ(x1a-0.04,(y1a+y2a)*0.5,':F22:c:B:iy:N: (m s:S:-1:N:)',0.01,90.0,0.0)

    call SET(x1b,x2b,y1b,y2b,-40.0,40.0,-40.0,40.0,1)
    call GSPLCI(1)
    call GSLWSC(2.0)
    do i=1,ndat
      if ( cpxst(i) /= msgv ) then
        call POINT(cpxst(i),cpyst(i))
      end if
    end do
    call GSLWSC(1.0)
    call GRIDAL(8,4,8,4,0,0,5,0.0,0.0)
    call SET(x1b-0.05,x2b+0.05,y1b-0.05,y2b+0.05,x1b-0.05,x2b+0.05,y1b-0.05,y2b+0.05,1)
    do i=1,9
      call PLCHHQ(x1b+(i-1)*(x2b-x1b)/8.,y1b-0.01,trim(xlab(i)),0.008,0.0,0.0)
    end do
    do j=1,9
      call PLCHHQ(x1b-0.005,y1b+(j-1)*(y2b-y1b)/8.,trim(ylab(j)),0.008,0.0,1.0)
    end do
    call PLCHHQ((x1b+x2b)*0.5,y2b+0.015,':F22:Ground-relative phase velocity'  ,0.01, 0.0,0.0)
    call PLCHHQ((x1b+x2b)*0.5,y1b-0.03,':F22:c:B:px:N: (m s:S:-1:N:)',0.01,0.0,0.0)
    call PLCHHQ(x1b-0.04,(y1b+y2b)*0.5,':F22:c:B:py:N: (m s:S:-1:N:)',0.01,90.0,0.0)

    call SET(x1c,x2c,y1c,y2c,-40.0,40.0,-40.0,40.0,1)
    call GSPLCI(1)
    call GSLWSC(2.0)
    do i=1,ndat
      if ( cgxst(i) /= msgv ) then
        call POINT(cgxst(i),cgyst(i))
      end if
    end do
    call GSLWSC(1.0)
    call GRIDAL(8,4,8,4,0,0,5,0.0,0.0)
    call SET(x1c-0.05,x2c+0.05,y1c-0.05,y2c+0.05,x1c-0.05,x2c+0.05,y1c-0.05,y2c+0.05,1)
    do i=1,9
      call PLCHHQ(x1c+(i-1)*(x2c-x1c)/8.,y1c-0.01,trim(xlab(i)),0.008,0.0,0.0)
    end do
    do j=1,9
      call PLCHHQ(x1c-0.005,y1c+(j-1)*(y2c-y1c)/8.,trim(ylab(j)),0.008,0.0,1.0)
    end do
    call PLCHHQ((x1c+x2c)*0.5,y2c+0.015,':F22:Group velocity'  ,0.01, 0.0,0.0)
    call PLCHHQ((x1c+x2c)*0.5,y1c-0.03,':F22:c:B:gx:N: (m s:S:-1:N:)',0.01,0.0,0.0)
    call PLCHHQ(x1c-0.04,(y1c+y2c)*0.5,':F22:c:B:gy:N: (m s:S:-1:N:)',0.01,90.0,0.0)

    call SET(x1d,x2d,y1d,y2d,-0.5,0.5,-0.5,0.5,1)
    call GSFAIS(1)
    call GSFACI(2)
    call GSPLCI(1)
    call GASETC('XLF.','(f4.1)')
    call GASETC('YLF.','(f4.1)')
    call GASETR('XLO.',0.005)
    call GASETR('YLO.',0.005)
    call DRAWHODO(12,angs,angst)
    call GSPLCI(1)
    call GSFACI(1)
    call GRIDAL(10,5,10,5,1,1,10,0,0)
    call SET(x1d-0.05,x2d+0.05,y1d-0.05,y2d+0.05,x1d-0.05,x2d+0.05,y1d-0.05,y2d+0.05,1)
    call PLCHHQ((x1d+x2d)*0.5,y2d+0.015,':F22:Energy (Et) weighted propagation',0.01,0.0,0.0)

    call SET(x1e,x2e,y1e,y2e,-0.01,0.01,-0.01,0.01,1)
    call GSPLCI(1)
    call GSLWSC(2.0)
    do i=1,ndat
      if ( mruwst(i) /= msgv ) then
        call POINT(mruwst(i),mrvwst(i))
      end if
    end do
    call GSLWSC(1.0)
    call GRIDAL(8,4,8,4,0,0,5,0.0,0.0)
    call SET(x1e-0.05,x2e+0.05,y1e-0.05,y2e+0.05,x1e-0.05,x2e+0.05,y1e-0.05,y2e+0.05,1)
    do i=1,9
      call PLCHHQ(x1e+(i-1)*(x2e-x1e)/8.,y1e-0.01,trim(xlab1(i)),0.008,0.0,0.0)
    end do
    do j=1,9
      call PLCHHQ(x1e-0.005,y1e+(j-1)*(y2e-y1e)/8.,trim(ylab1(j)),0.008,0.0,1.0)
    end do
    call PLCHHQ((x1e+x2e)*0.5,y2e+0.015,':F22:Momentum flux'  ,0.01, 0.0,0.0)
    call PLCHHQ((x1e+x2e)*0.5,y1e-0.03,':F22:rho u'' w'' (x10:S:-3:N: N m:S:-2:N:)',0.01,0.0,0.0)
    call PLCHHQ(x1e-0.04,(y1e+y2e)*0.5,':F22:rho v'' w'' (x10:S:-3:N: N m:S:-2:N:)',0.01,90.0,0.0)

    call GDAWK(1)

    call GACWK(2)

    call SET(0.0,1.0,0.0,1.0,0.0,1.0,0.0,1.0,1)
    call PLCHHQ(0.5,0.985,':F22:'//trim(titletr),0.012,0.0,0.0)

    call SET(x1a,x2a,y1a,y2a,-40.0,40.0,-40.0,40.0,1)
    call GSPLCI(1)
    call GSLWSC(2.0)
    do i=1,ndat
      if ( cixtr(i) /= msgv ) then
        call POINT(cixtr(i),ciytr(i))
      end if
    end do
    call GSLWSC(1.0)
    call GRIDAL(8,4,8,4,0,0,5,0.0,0.0)
    call SET(x1a-0.05,x2a+0.05,y1a-0.05,y2a+0.05,x1a-0.05,x2a+0.05,y1a-0.05,y2a+0.05,1)
    do i=1,9
      call PLCHHQ(x1a+(i-1)*(x2a-x1a)/8.,y1a-0.01,trim(xlab(i)),0.008,0.0,0.0)
    end do
    do j=1,9
      call PLCHHQ(x1a-0.005,y1a+(j-1)*(y2a-y1a)/8.,trim(ylab(j)),0.008,0.0,1.0)
    end do
    call PLCHHQ((x1a+x2a)*0.5,y2a+0.015,':F22:Intrinsic phase speed'  ,0.01, 0.0,0.0)
    call PLCHHQ((x1a+x2a)*0.5,y1a-0.03,':F22:c:B:ix:N: (m s:S:-1:N:)',0.01,0.0,0.0)
    call PLCHHQ(x1a-0.04,(y1a+y2a)*0.5,':F22:c:B:iy:N: (m s:S:-1:N:)',0.01,90.0,0.0)

    call SET(x1b,x2b,y1b,y2b,-40.0,40.0,-40.0,40.0,1)
    call GSPLCI(1)
    call GSLWSC(2.0)
    do i=1,ndat
      if ( cpxtr(i) /= msgv ) then
        call POINT(cpxtr(i),cpytr(i))
      end if
    end do
    call GSLWSC(1.0)
    call GRIDAL(8,4,8,4,0,0,5,0.0,0.0)
    call SET(x1b-0.05,x2b+0.05,y1b-0.05,y2b+0.05,x1b-0.05,x2b+0.05,y1b-0.05,y2b+0.05,1)
    do i=1,9
      call PLCHHQ(x1b+(i-1)*(x2b-x1b)/8.,y1b-0.01,trim(xlab(i)),0.008,0.0,0.0)
    end do
    do j=1,9
      call PLCHHQ(x1b-0.005,y1b+(j-1)*(y2b-y1b)/8.,trim(ylab(j)),0.008,0.0,1.0)
    end do
    call PLCHHQ((x1b+x2b)*0.5,y2b+0.015,':F22:Ground-relative phase velocity'  ,0.01, 0.0,0.0)
    call PLCHHQ((x1b+x2b)*0.5,y1b-0.03,':F22:c:B:ix:N: (m s:S:-1:N:)',0.01,0.0,0.0)
    call PLCHHQ(x1b-0.04,(y1b+y2b)*0.5,':F22:c:B:iy:N: (m s:S:-1:N:)',0.01,90.0,0.0)

    call SET(x1c,x2c,y1c,y2c,-40.0,40.0,-40.0,40.0,1)
    call GSPLCI(1)
    call GSLWSC(2.0)
    do i=1,ndat
      if ( cgxtr(i) /= msgv ) then
        call POINT(cgxtr(i),cgytr(i))
      end if
    end do
    call GSLWSC(1.0)
    call GRIDAL(8,4,8,4,0,0,5,0.0,0.0)
    call SET(x1c-0.05,x2c+0.05,y1c-0.05,y2c+0.05,x1c-0.05,x2c+0.05,y1c-0.05,y2c+0.05,1)
    do i=1,9
      call PLCHHQ(x1c+(i-1)*(x2c-x1c)/8.,y1c-0.01,trim(xlab(i)),0.008,0.0,0.0)
    end do
    do j=1,9
      call PLCHHQ(x1c-0.005,y1c+(j-1)*(y2c-y1c)/8.,trim(ylab(j)),0.008,0.0,1.0)
    end do
    call PLCHHQ((x1c+x2c)*0.5,y2c+0.015,':F22:Group velocity'  ,0.01, 0.0,0.0)
    call PLCHHQ((x1c+x2c)*0.5,y1c-0.03,':F22:c:B:gx:N: (m s:S:-1:N:)',0.01,0.0,0.0)
    call PLCHHQ(x1c-0.04,(y1c+y2c)*0.5,':F22:c:B:gy:N: (m s:S:-1:N:)',0.01,90.0,0.0)

    call SET(x1d,x2d,y1d,y2d,-0.5,0.5,-0.5,0.5,1)
    call GSFAIS(1)
    call GSFACI(2)
    call GSPLCI(1)
    call GASETC('XLF.','(f4.1)')
    call GASETC('YLF.','(f4.1)')
    call GASETR('XLO.',0.005)
    call GASETR('YLO.',0.005)
    call DRAWHODO(12,angs,angtr)
    call GSPLCI(1)
    call GSFACI(1)
    call GRIDAL(10,5,10,5,1,1,10,0,0)
    call SET(x1d-0.05,x2d+0.05,y1d-0.05,y2d+0.05,x1d-0.05,x2d+0.05,y1d-0.05,y2d+0.05,1)
    call PLCHHQ((x1d+x2d)*0.5,y2d+0.015,':F22:Energy (Et) weighted propagation',0.01,0.0,0.0)

    call SET(x1e,x2e,y1e,y2e,-0.01,0.01,-0.01,0.01,1)
    call GSPLCI(1)
    call GSLWSC(2.0)
    do i=1,ndat
      if ( mruwtr(i) /= msgv ) then
        call POINT(mruwtr(i),mrvwtr(i))
      end if
    end do
    call GSLWSC(1.0)
    call GRIDAL(8,4,8,4,0,0,5,0.0,0.0)
    call SET(x1e-0.05,x2e+0.05,y1e-0.05,y2e+0.05,x1e-0.05,x2e+0.05,y1e-0.05,y2e+0.05,1)
    do i=1,9
      call PLCHHQ(x1e+(i-1)*(x2e-x1e)/8.,y1c-0.01,trim(xlab1(i)),0.008,0.0,0.0)
    end do
    do j=1,9
      call PLCHHQ(x1e-0.005,y1e+(j-1)*(y2e-y1e)/8.,trim(ylab1(j)),0.008,0.0,1.0)
    end do
    call PLCHHQ((x1e+x2e)*0.5,y2e+0.015,':F22:Momentum flux'  ,0.01, 0.0,0.0)
    call PLCHHQ((x1e+x2e)*0.5,y1e-0.03,':F22:rho u'' w'' (x10:S:-3:N: N m:S:-2:N:)',0.01,0.0,0.0)
    call PLCHHQ(x1e-0.04,(y1e+y2e)*0.5,':F22:rho v'' w'' (x10:S:-3:N: N m:S:-2:N:)',0.01,90.0,0.0)

    call GDAWK(2)

    call FRAME

    return

    end subroutine pltspec

    subroutine drawhodo(n,deg,amp)

    integer,            intent(in) :: n
    real, dimension(n), intent(in) :: deg  ,amp

    integer, parameter :: np=302
    integer :: i,j
    real :: pi,nip,nfp 
    real, dimension(np) :: x,y,phi     

    pi = 4.0*atan(1.0)

    do i=1,12
      print *,deg(i)
      nip=((deg(i)-29.0))*pi/180.0
      nfp=((deg(i)+ 1.0))*pi/180.0
      do j=1,np-2
        phi(j) = (j-1)*(nfp-nip)/float(np-3) + nip
        x(j) = amp(i)*cos(phi(j))
        y(j) = amp(i)*sin(phi(j))
      end do
      x(np-1) = 0.0
      y(np-1) = 0.0
      x(np) = x(1)
      y(np) = y(1)
      call gfa(np-1,x,y)
      call curve(x,y,np)
    end do

    return
    end subroutine drawhodo

    end module graphicspec
