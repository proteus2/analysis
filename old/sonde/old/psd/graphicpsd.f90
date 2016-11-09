!2345
    module graphicpsd

    use varpsd
    use varnfit

    contains

    subroutine pltpsd(it    ,date  ,nrst  ,nrtr  ,stanal,tranal)

    implicit none

    integer, intent(in)          :: it
    character(len=8), intent(in) :: date
    integer, intent(in)          :: nrst  ,nrtr
    logical, intent(in)          :: stanal,tranal

    character(len=40), dimension( 5) :: xlab
    character(len=40), dimension(12) :: ylab1
    character(len=40), dimension(11) :: ylab2

    integer :: i,j

    xlab  = (/'100','10','1','0.1','0.01'/)
    ylab1 = (/'10:S:-4:N:','10:S:-3:N:','10:S:-2:N:','10:S:-1:N:', &
              '10:S:0:N:' ,'10:S:1:N:' ,'10:S:2:N:' ,'10:S:2:N:' , &
              '10:S:3:N:' ,'10:S:4:N:' ,'10:S:5:N:' ,'10:S:6:N:'/)
    ylab2 = (/'10:S:-9:N:','10:S:-8:N:','10:S:-7:N:','10:S:-6:N:', &
              '10:S:-5:N:','10:S:-4:N:','10:S:-3:N:','10:S:-2:N:', &
              '10:S:-1:N:','10:S:0:N:' ,'10:S:1:N:' /)

    if ( stanal ) then
  
      call GACWK(1)

      call GSTXFP(-4,2)

      call GASETR('XMJ.',0.006)
      call GASETR('YMJ.',0.006)
      call GASETR('XMN.',0.003)
      call GASETR('YMN.',0.003)
      call GASETI('LTY.',1)
      call PCSETI('FN.',21)

      call SET(0.30,0.70,0.70,0.95,100.,0.01,1.e-04,1.e+06,4)
      call GSPLCI(1)
      call CURVED(wavlst(2:nrst/2+1),ustpsd1(2:nrst/2+1),nrst/2)
      call GSPLCI(2)
      call CURVED(wavlst(2:nrst/2+1),ustsat1(2:nrst/2+1),nrst/2)
      call GSPLCI(1)
      call GRIDAL(1,1,1,1,0,0,5,0.0,0.0)
      call SET(0.20,0.80,0.65,1.00,0.20,0.80,0.65,1.00,1)
      do i=1,5
        call PLCHHQ(0.30+(i-1)*0.4/4.,0.69,trim(xlab(i)),0.008,0.0,0.0)
      end do
      do j=1,12
        call PLCHHQ(0.29,0.70+(j-1)*0.25/11.,trim(ylab1(j)),0.008,0.0,1.0)
      end do
      call PLCHHQ(0.69,0.93,':F22:UPSD at '//date,0.01,0.0,1.0)
      call PLCHHQ(0.50,0.67,':F22:Vertical wavelength (km)',0.01,0.0,0.0)
      call PLCHHQ(0.25,0.825,':F22:PSD [(m:S:2:N: s:S:-2:N:)/(cycle m:S:-1:N:)]',0.01,90.0,0.0)

      call SET(0.30,0.70,0.40,0.65,100.,0.01,1.e-04,1.e+06,4)
      call CURVED(wavlst(2:nrst/2+1),vstpsd1(2:nrst/2+1),nrst/2)
      call GSPLCI(1)
      call GRIDAL(1,1,1,1,0,0,5,0.0,0.0)
      call SET(0.20,0.80,0.35,0.70,0.20,0.80,0.65,1.00,1)
      do i=1,5
        call PLCHHQ(0.30+(i-1)*0.4/4.,0.69,trim(xlab(i)),0.008,0.0,0.0)
      end do
      do j=1,12
        call PLCHHQ(0.29,0.70+(j-1)*0.25/11.,trim(ylab1(j)),0.008,0.0,1.0)
      end do
      call PLCHHQ(0.69,0.93,':F22:VPSD at '//date,0.01,0.0,1.0)
      call PLCHHQ(0.50,0.67,':F22:Vertical wavelength (km)',0.01,0.0,0.0)
      call PLCHHQ(0.25,0.825,':F22:PSD [(m:S:2:N: s:S:-2:N:)/(cycle m:S:-1:N:)]',0.01,90.0,0.0)

      call SET(0.30,0.70,0.10,0.35,100.,0.01,1.e-09,1.e+01,4)
      call GSPLCI(1)
      call CURVED(wavlst(2:nrst/2+1),tstpsd1 (2:nrst/2+1),nrst/2)
      call GSPLCI(3)
      call CURVED(wavlst(2:nrst/2+1),tstpsdc1(2:nrst/2+1),nrst/2)
      call GSPLCI(2)
      call CURVED(wavlst(2:nrst/2+1),tstsat1 (2:nrst/2+1),nrst/2)
      call GSPLCI(4)
      call CURVED(wavlst(2:nrst/2+1),fitst(2:nrst/2+1,it),nrst/2)
      call GSPLCI(1)
      call GRIDAL(1,1,1,1,0,0,5,0.0,0.0)
      call SET(0.20,0.80,0.05,0.40,0.20,0.80,0.65,1.00,1)
      do i=1,5
        call PLCHHQ(0.30+(i-1)*0.4/4.,0.69,trim(xlab(i)),0.008,0.0,0.0)
      end do
      do j=1,11
        call PLCHHQ(0.29,0.70+(j-1)*0.25/10.,trim(ylab2(j)),0.008,0.0,1.0)
      end do
      call PLCHHQ(0.69,0.93, ':F22:TPSD at '//date,0.01,0.0,1.0)
      call PLCHHQ(0.50,0.67, ':F22:Vertical wavelength (km)',0.01,0.0,0.0)
      call PLCHHQ(0.25,0.825,':F22:PSD [(cycle/m):S:-1:N:]',0.01,90.0,0.0)

      call GDAWK(1)

    end if

    if ( tranal ) then

      call GACWK(2)

      call GSTXFP(-4,2)

      call GASETR('XMJ.',0.006)
      call GASETR('YMJ.',0.006)
      call GASETR('XMN.',0.003)
      call GASETR('YMN.',0.003)
      call GASETI('LTY.',1)
      call PCSETI('FN.',21)

      call SET(0.30,0.70,0.70,0.95,100.,0.01,1.e-04,1.e+06,4)
      call GSPLCI(1)
      call CURVED(wavltr(2:nrtr/2+1),utrpsd1(2:nrtr/2+1),nrtr/2)
      call GSPLCI(2)
      call CURVED(wavltr(2:nrtr/2+1),utrsat1(2:nrtr/2+1),nrtr/2)
      call GSPLCI(1)
      call GRIDAL(1,1,1,1,0,0,5,0.0,0.0)
      call SET(0.20,0.80,0.65,1.00,0.20,0.80,0.65,1.00,1)
      do i=1,5
        call PLCHHQ(0.30+(i-1)*0.4/4.,0.69,trim(xlab(i)),0.008,0.0,0.0)
      end do
      do j=1,12
        call PLCHHQ(0.29,0.70+(j-1)*0.25/11.,trim(ylab1(j)),0.008,0.0,1.0)
      end do
      call PLCHHQ(0.69,0.93,':F22:UPSD at '//date,0.01,0.0,1.0)
      call PLCHHQ(0.50,0.67,':F22:Vertical wavelength (km)',0.01,0.0,0.0)
      call PLCHHQ(0.25,0.825,':F22:PSD [(m:S:2:N: s:S:-2:N:)/(cycle m:S:-1:N:)]',0.01,90.0,0.0)

      call SET(0.30,0.70,0.40,0.65,100.,0.01,1.e-04,1.e+06,4)
      call CURVED(wavltr(2:nrtr/2+1),vtrpsd1(2:nrtr/2+1),nrtr/2)
      call GSPLCI(1)
      call GRIDAL(1,1,1,1,0,0,5,0.0,0.0)
      call SET(0.20,0.80,0.35,0.70,0.20,0.80,0.65,1.00,1)
      do i=1,5
        call PLCHHQ(0.30+(i-1)*0.4/4.,0.69,trim(xlab(i)),0.008,0.0,0.0)
      end do
      do j=1,12
        call PLCHHQ(0.29,0.70+(j-1)*0.25/11.,trim(ylab1(j)),0.008,0.0,1.0)
      end do
      call PLCHHQ(0.69,0.93,':F22:VPSD at '//date,0.01,0.0,1.0)
      call PLCHHQ(0.50,0.67,':F22:Vertical wavelength (km)',0.01,0.0,0.0)
      call PLCHHQ(0.25,0.825,':F22:PSD [(m:S:2:N: s:S:-2:N:)/(cycle m:S:-1:N:)]',0.01,90.0,0.0)

      call SET(0.30,0.70,0.10,0.35,100.,0.01,1.e-09,1.e+01,4)
      call GSPLCI(1)
      call CURVED(wavltr(2:nrtr/2+1),ttrpsd1 (2:nrtr/2+1),nrtr/2)
      call GSPLCI(2)
      call CURVED(wavltr(2:nrtr/2+1),ttrsat1 (2:nrtr/2+1),nrtr/2)
      call GSPLCI(4)
      call CURVED(wavltr(2:nrtr/2+1),fittr(2:nrtr/2+1,it),nrtr/2)
      call GSPLCI(1)
      call GRIDAL(1,1,1,1,0,0,5,0.0,0.0)
      call SET(0.20,0.80,0.05,0.40,0.20,0.80,0.65,1.00,1)
      do i=1,5
        call PLCHHQ(0.30+(i-1)*0.4/4.,0.69,trim(xlab(i)),0.008,0.0,0.0)
      end do
      do j=1,11
        call PLCHHQ(0.29,0.70+(j-1)*0.25/10.,trim(ylab2(j)),0.008,0.0,1.0)
      end do
      call PLCHHQ(0.69,0.93, ':F22:TPSD at '//date,0.01,0.0,1.0)
      call PLCHHQ(0.50,0.67, ':F22:Vertical wavelength (km)',0.01,0.0,0.0)
      call PLCHHQ(0.25,0.825,':F22:PSD [(cycle/m):S:-1:N:]',0.01,90.0,0.0)

      call GDAWK(2)

    end if

    call FRAME

    return
    end subroutine pltpsd

    end module graphicpsd
