!2345
    program specanal

!----------------------------------------------------------------------------------
!
!   Purpose:
!
!   To estimate several spectral characteristics of gravity waves 
!
!   Author:
!
!   In-Sun Song
!   Laboratory for Mesoscale Dynamics
!   Department of Atmospheric Sciences, Yonsei University, Seoul, Korea 
!
!   Version history:
!
!   17 November, 2000
!   Revised from previous program
!
!   6 May, 2002
!   Added description and standardized source code  
!    
!   5 December, 2002
!   Combined tropospheric analysis and stratospheric analysis
!   Changed IO routine to read input data
!   Changed IO routine to dump results as a NetCDF data file
!
!   Required library:
!
!   NetCDF (higher version than 3.4)
!   NCAR GRAPHICS (higher version that 4.0)
!   HILBERT module for Hilbert transformation
!   CFFTPACK
!
!----------------------------------------------------------------------------------

    use varspec
    use meanwavnum
    use propagate
    use hilbert
    use momentum
    use netcdfio
    use graphicspec    

    implicit none 

    integer            :: ndat   ,nrst   ,nrtr 
    real               :: pi     ,f     
    real, parameter    :: grav=9.806 ,pp=5./3. ,omega=0.00007292 ,eps=1.e-5

!----------------------------------------------------------------------------------
!
!   Input variables
!
!----------------------------------------------------------------------------------
 
    character(len=100) :: rdir
    character(len=1)   :: leap
    integer            :: year  ,month  ,year2
    real               :: lat   ,tau

    character(len=8), dimension(1000) :: addmsg  

!----------------------------------------------------------------------------------
!
!   Miscellaneous variables
!
!----------------------------------------------------------------------------------

    integer                 :: id     ,ih      ,it     ,k     ,naddmsg

    logical                 :: stanal ,tranal
    logical                 :: iex1   ,iex2

    integer, dimension(12)  :: daymon

    character(len=100)      :: tmpnam1,tmpnam2
    character(len=100)      :: grimfst,grimftr

    integer                 :: yrmsg  ,monmsg ,daymsg ,hrmsg

    integer                 :: ncidst1,ncidst2,ncidst3,ncidst4
    integer                 :: ncidst5,ncidst6,ncidst7,ncidst8
    integer                 :: ncidst9
    integer                 :: ncidtr1,ncidtr2,ncidtr3,ncidtr4
    integer                 :: ncidtr5,ncidtr6,ncidtr7,ncidtr8
    integer                 :: ncidtr9

    character(len=100)      :: rfnst1 ,rfnst2 ,rfnst3 ,rfnst4
    character(len=100)      :: rfnst5 ,rfnst6 ,rfnst7 ,rfnst8
    character(len=100)      :: rfnst9
    character(len=100)      :: rfntr1 ,rfntr2 ,rfntr3 ,rfntr4
    character(len=100)      :: rfntr5 ,rfntr6 ,rfntr7 ,rfntr8
    character(len=100)      :: rfntr9

    character(len=100)      :: wfnst01,wfnst02,wfnst03,wfnst04
    character(len=100)      :: wfnst05,wfnst06,wfnst07,wfnst08
    character(len=100)      :: wfnst09,wfnst10,wfnst11,wfnst12
    character(len=100)      :: wfnst13
    character(len=100)      :: wfntr01,wfntr02,wfntr03,wfntr04
    character(len=100)      :: wfntr05,wfntr06,wfntr07,wfntr08
    character(len=100)      :: wfntr09,wfntr10,wfntr11,wfntr12
    character(len=100)      :: wfntr13

    character(len= 8)       :: date

    real                    :: msgv   ,dz

    namelist /specparm/  rdir   ,year   ,month   ,leap   ,tau   ,lat   ,addmsg  

!----------------------------------------------------------------------------------
!
!   Begin program
!
!----------------------------------------------------------------------------------

!----------------------------------------------------------------------------------
!   Define some constants
!----------------------------------------------------------------------------------

    do k=1,1000
      addmsg(k) = 'xxxxxxxx'
    end do

!----------------------------------------------------------------------------------
!   Read input parameters 
!----------------------------------------------------------------------------------

    read(5,specparm)

!----------------------------------------------------------------------------------
!   Define parameters
!----------------------------------------------------------------------------------

    pi = 4.0*atan(1.0)

    year2 = year
    if ( year >= 2000 ) then
      year = year - 2000
    else if ( year >= 1900 ) then
      year = year - 1900
    end if

    f = 2.0*omega*sin(lat*pi/180.0)  ! f is Coriolis parameter

    daymon = (/31,28,31,30,31,30,31,31,30,31,30,31/)
    if ( leap == 'T' ) then
      daymon(2) = 29
    end if

    do k=1,1000
      if ( addmsg(k) == 'xxxxxxxx' ) then
        naddmsg = k-1
        go to 10
      end if
    end do
10  continue

100 format(a,i2.2,i2.2,a)

!----------------------------------------------------------------------------------
!
!   Initialize NCAR graphics 
!
!----------------------------------------------------------------------------------

    call GOPKS(30,0)

    write(grimfst,100) trim(rdir)//'/st/stspec_',year,month,'.ps'
    write(grimftr,100) trim(rdir)//'/tr/trspec_',year,month,'.ps'

    call NGSETC('ME.',trim(grimfst))
    call GOPWK(1,0,20)

    call GSCR(1,0,1.000,1.000,1.000)
    call GSCR(1,1,0.000,0.000,0.000)
    call GSCR(1,2,0.700,0.000,0.000)
    call GSCR(1,3,0.000,0.000,0.700)
    call GSCR(1,4,0.000,0.700,0.000)

    call NGSETC('ME',trim(grimftr))
    call GOPWK(2,0,20)

    call GSCR(2,0,1.000,1.000,1.000)
    call GSCR(2,1,0.000,0.000,0.000)
    call GSCR(2,2,0.700,0.000,0.000)
    call GSCR(2,3,0.000,0.000,0.700)
    call GSCR(2,4,0.000,0.700,0.000)

!----------------------------------------------------------------------------------
!
!   Reading data
!
!----------------------------------------------------------------------------------

    write(rfnst1,100) trim(rdir)//'/st/avgetst',year,month, '.nc'
    write(rfnst2,100) trim(rdir)//'/st/bvstbas',year,month, '.nc'
    write(rfnst3,100) trim(rdir)//'/st/ustprt' ,year,month, '.nc'
    write(rfnst4,100) trim(rdir)//'/st/vstprt' ,year,month, '.nc'
    write(rfnst5,100) trim(rdir)//'/st/tstprt' ,year,month, '.nc'
    write(rfnst6,100) trim(rdir)//'/st/ustbas' ,year,month, '.nc'
    write(rfnst7,100) trim(rdir)//'/st/vstbas' ,year,month, '.nc'
    write(rfnst8,100) trim(rdir)//'/st/tstbas' ,year,month, '.nc'
    write(rfnst9,100) trim(rdir)//'/st/rstbas' ,year,month, '.nc'
    write(rfntr1,100) trim(rdir)//'/tr/avgettr',year,month, '.nc'
    write(rfntr2,100) trim(rdir)//'/tr/bvtrbas',year,month, '.nc'
    write(rfntr3,100) trim(rdir)//'/tr/utrprt' ,year,month, '.nc'
    write(rfntr4,100) trim(rdir)//'/tr/vtrprt' ,year,month, '.nc'
    write(rfntr5,100) trim(rdir)//'/tr/ttrprt' ,year,month, '.nc'
    write(rfntr6,100) trim(rdir)//'/tr/utrbas' ,year,month, '.nc'
    write(rfntr7,100) trim(rdir)//'/tr/vtrbas' ,year,month, '.nc'
    write(rfntr8,100) trim(rdir)//'/tr/ttrbas' ,year,month, '.nc'
    write(rfntr9,100) trim(rdir)//'/tr/rtrbas' ,year,month, '.nc'

    call opennc(trim(rfnst1),ncidst1)
    call opennc(trim(rfnst2),ncidst2)
    call opennc(trim(rfnst3),ncidst3)
    call opennc(trim(rfnst4),ncidst4)
    call opennc(trim(rfnst5),ncidst5)
    call opennc(trim(rfnst6),ncidst6)
    call opennc(trim(rfnst7),ncidst7)
    call opennc(trim(rfnst8),ncidst8)
    call opennc(trim(rfnst9),ncidst9)
    call opennc(trim(rfntr1),ncidtr1)
    call opennc(trim(rfntr2),ncidtr2)
    call opennc(trim(rfntr3),ncidtr3)
    call opennc(trim(rfntr4),ncidtr4)
    call opennc(trim(rfntr5),ncidtr5)
    call opennc(trim(rfntr6),ncidtr6)
    call opennc(trim(rfntr7),ncidtr7)
    call opennc(trim(rfntr8),ncidtr8)
    call opennc(trim(rfntr9),ncidtr9)

    call getdimlen(ncidst1,'TIME',ndat)
    call getdimlen(ncidst3,'Z'   ,nrst)
    call getdimlen(ncidtr3,'Z'   ,nrtr)

    call initvarspec(ndat ,nrst ,nrtr )

    call get1d(ncidst1,'TIME'  ,ndat    ,normt )
    call get1d(ncidst3,'Z'     ,nrst    ,zst   )
    call get1d(ncidtr3,'Z'     ,nrtr    ,ztr   )

    dz = zst(2)-zst(1)
    msgv = 1.e+32
    print *, 'dz', dz, 'nrst', nrst, 'nrtr', nrtr

    call get1d(ncidst1,'AVGETST',ndat     ,avgstet)
    call get2d(ncidst2,'BVBAS'  ,ndat,nrst,nstbas )
    call get2d(ncidst3,'UPRT'   ,ndat,nrst,ustprt )
    call get2d(ncidst4,'VPRT'   ,ndat,nrst,vstprt )
    call get2d(ncidst5,'TPRT'   ,ndat,nrst,tstprt )
    call get2d(ncidst6,'UBAS'   ,ndat,nrst,ustbas )
    call get2d(ncidst7,'VBAS'   ,ndat,nrst,vstbas )
    call get2d(ncidst8,'TBAS'   ,ndat,nrst,tstbas )
    call get2d(ncidst9,'RBAS'   ,ndat,nrst,rstbas )
    call get1d(ncidtr1,'AVGETTR',ndat     ,avgtret)
    call get2d(ncidtr2,'BVBAS  ',ndat,nrtr,ntrbas )
    call get2d(ncidtr3,'UPRT'   ,ndat,nrtr,utrprt )
    call get2d(ncidtr4,'VPRT'   ,ndat,nrtr,vtrprt )
    call get2d(ncidtr5,'TPRT'   ,ndat,nrtr,ttrprt )
    call get2d(ncidtr6,'UBAS'   ,ndat,nrtr,utrbas )
    call get2d(ncidtr7,'VBAS'   ,ndat,nrtr,vtrbas )
    call get2d(ncidtr8,'TBAS'   ,ndat,nrtr,ttrbas )
    call get2d(ncidtr9,'RBAS'   ,ndat,nrtr,rtrbas )


    do it=1,ndat
      do k=1,nrst
        avgstn(it) = avgstn(it) + nstbas(it,k)/float(nrst)
      end do
      do k=1,nrtr
        avgtrn(it) = avgtrn(it) + ntrbas(it,k)/float(nrtr)
      end do
      if (ntrbas(it,3).eq.msgv) then
        avgtrn(it) = msgv
      end if  
      write(8,*) it,'avgtrn', avgtrn(it)
    end do  



    call closenc(ncidst1)
    call closenc(ncidst2)
    call closenc(ncidst3)
    call closenc(ncidst4)
    call closenc(ncidst5)
    call closenc(ncidst6)
    call closenc(ncidst7)
    call closenc(ncidst8)
    call closenc(ncidst9)
    call closenc(ncidtr1)
    call closenc(ncidtr2)
    call closenc(ncidtr3)
    call closenc(ncidtr4)
    call closenc(ncidtr5)
    call closenc(ncidtr6)   
    call closenc(ncidtr7)   
    call closenc(ncidtr8)   
    call closenc(ncidtr9)   

!----------------------------------------------------------------------------------
!
!   Determine output file name
!
!----------------------------------------------------------------------------------

    write(wfnst01,100) trim(rdir)//'/result/st/stpdir_',year,month,'.nc'
    write(wfnst02,100) trim(rdir)//'/result/st/ststok_',year,month,'.nc'
    write(wfnst03,100) trim(rdir)//'/result/st/stmekm_',year,month,'.nc'
    write(wfnst04,100) trim(rdir)//'/result/st/stcpin_',year,month,'.nc'
    write(wfnst05,100) trim(rdir)//'/result/st/stgvel_',year,month,'.nc'
    write(wfnst06,100) trim(rdir)//'/result/st/stwcgr_',year,month,'.nc'
    write(wfnst07,100) trim(rdir)//'/result/st/stcpgr_',year,month,'.nc'
    write(wfnst08,100) trim(rdir)//'/result/st/strota_',year,month,'.nc'
    write(wfnst09,100) trim(rdir)//'/result/st/strpsd_',year,month,'.nc'
    write(wfnst10,100) trim(rdir)//'/result/st/stangh_',year,month,'.nc'
    write(wfnst11,100) trim(rdir)//'/result/st/stmdir_',year,month,'.nc'
    write(wfnst12,100) trim(rdir)//'/result/st/struvw_',year,month,'.nc'
    write(wfnst13,100) trim(rdir)//'/result/st/stuvwf_',year,month,'.nc'

    write(wfntr01,100) trim(rdir)//'/result/tr/trpdir_',year,month,'.nc'
    write(wfntr02,100) trim(rdir)//'/result/tr/trtrok_',year,month,'.nc'
    write(wfntr03,100) trim(rdir)//'/result/tr/trmekm_',year,month,'.nc'
    write(wfntr04,100) trim(rdir)//'/result/tr/trcpin_',year,month,'.nc'
    write(wfntr05,100) trim(rdir)//'/result/tr/trgvel_',year,month,'.nc'
    write(wfntr06,100) trim(rdir)//'/result/tr/trwcgr_',year,month,'.nc'
    write(wfntr07,100) trim(rdir)//'/result/tr/trcpgr_',year,month,'.nc'
    write(wfntr08,100) trim(rdir)//'/result/tr/trrota_',year,month,'.nc'
    write(wfntr09,100) trim(rdir)//'/result/tr/trrpsd_',year,month,'.nc'
    write(wfntr10,100) trim(rdir)//'/result/tr/trangh_',year,month,'.nc'
    write(wfntr11,100) trim(rdir)//'/result/tr/trmdir_',year,month,'.nc'
    write(wfntr12,100) trim(rdir)//'/result/tr/trruvw_',year,month,'.nc'
    write(wfntr13,100) trim(rdir)//'/result/tr/truvwf_',year,month,'.nc'

!----------------------------------------------------------------------------------
!
!   Loop
!
!----------------------------------------------------------------------------------

      do it=1,ndat 
        
        print *,'it',it
     
        stanal = .true.
        tranal = .true.

        if ( avgstn(it) .eq. msgv .or. nstbas(it,5) .eq. msgv ) stanal = .false.
        if ( avgtrn(it) .eq. msgv .or. ntrbas(it,5) .eq. msgv ) tranal = .false. 

        if (it.eq.17 .or. it.eq.29) then
          stanal = .false.                           ! OSAN-ERROR
          tranal = .false.
        end if

        if ( stanal ) then

!----------------------------------------------------------------------------------
!         Inverse Hilbert transform of normalized temperature perturbation
!         in the stratosphere
!----------------------------------------------------------------------------------

          do k=1,nrst
            ustprt1(k) = ustprt(it,k)
            vstprt1(k) = vstprt(it,k)
            ustbas1(k) = ustbas(it,k)
            vstbas1(k) = vstbas(it,k)
            rstbas1(k) = rstbas(it,k)
            tnstprt(k) = tstprt(it,k)/tstbas(it,k)
          end do

          call invhil(nrst  ,dz    ,tnstprt, tnsthil)     

!----------------------------------------------------------------------------------
!         Find wave propagation direction 
!----------------------------------------------------------------------------------

          call propdir(nrst       ,ustprt1    ,vstprt1    ,tnsthil    ,  &
                       phist(it)  ,degst(it)  ,degnst(it) )        

!----------------------------------------------------------------------------------
!         Calculate Stoke's parameter 
!----------------------------------------------------------------------------------

          call stoke  (nrst       ,dz         ,ustprt1    ,vstprt1    ,  &
                       ustbas1    ,vstbas1    ,degst(it)  ,avgstn(it) ,  &
                       w_fst(it)  ,dfst(it)   )

!----------------------------------------------------------------------------------
!         Calculate mean vertical and horizontal wavenumbers 
!----------------------------------------------------------------------------------

          call meanmk (nrst       ,dz         ,f          ,w_fst(it)  ,  &
                       avgstn(it) ,ustprt1    ,vstprt1    ,mbst(it)   ,  & 
                       mbast(it)  ,kbst(it)   ,kbast(it)  )

!----------------------------------------------------------------------------------
!         Calculate phase and group velocities 
!----------------------------------------------------------------------------------

          call phase  (nrst       ,f          ,w_fst(it)  ,kbast(it)  ,  &
                       mbast(it)  ,degst(it)  ,ustbas1    ,vstbas1    ,  & 
                       czst(it)   ,cist(it)   ,cixst(it)  ,ciyst(it)  ,  &
                       cgxst(it)  ,cgyst(it)  ,wgrst(it)  ,cgrst(it)  ,  &
                       cxst(it)   ,cyst(it)   ,cpxst(it)  ,cpyst(it)  ,  &
                       ubst(it)   ,ubdgst(it) )

!----------------------------------------------------------------------------------
!         Estimate a fraction of upward propagating waves by calculating
!         rotary spectra.
!----------------------------------------------------------------------------------

          call rotary (nrst       ,f          ,dz         ,ustprt1    ,  &
                       vstprt1    ,wavmst     ,rstpsd1    ,cclst(it)  ,  &
                       clst(it)   ,upfrst(it) ,dnfrst(it) )

          rstpsd(1:nrst,it) = rstpsd1(1:nrst)

!----------------------------------------------------------------------------------
!         Estimate wave momentum flux 
!----------------------------------------------------------------------------------

          call uwvw   (nrst       ,avgstn(it) ,w_fst(it)  ,f          ,  &
                       ustprt1    ,vstprt1    ,tnsthil    ,rstbas1    ,  &
                       ruwfst1    ,rvwfst1    ,uwflst1    ,vwflst1    ,  &
                       mruwst(it) ,mrvwst(it) ,muwfst(it) ,mvwfst(it) )

          ruwfst(it,1:nrst) = ruwfst1(1:nrst)
          rvwfst(it,1:nrst) = rvwfst1(1:nrst)
          uwflst(it,1:nrst) = uwflst1(1:nrst)
          vwflst(it,1:nrst) = vwflst1(1:nrst)

!----------------------------------------------------------------------------------
!   If mean axial ratio of hodograph is less than one, skip the data
!----------------------------------------------------------------------------------
!----------------------------------------------------------------------------------
!   Following Vincent and Alexander (2000), a small fraction of the frequencies  
!   inferred from soundings had values greater than 10f. These values are not
!   reliable since errors in the Digicora wind measuring system will limit the
!   accuracy to which the transverse velocity perturbation of the hodograph can
!   be determined. Digicora wind speed accuracy is believed to be about 0.4 - 0.5
!   m/s (Ivanov et al. 1991). In Vincent and Alexander (2000), rms wave perturbation
!   is about 4 m/s, which suggests that the axial ratios and hence intrinsic    
!   frequencies greater than 10f are not reliable and should be discarded when
!   computing other wave parameters, such as wavelength. 
!----------------------------------------------------------------------------------
 
          if (w_fst(it) < 1.0 .or. w_fst(it) > 10.0 ) then
            phist(it)     = msgv
            degst(it)     = msgv
            degnst(it)    = msgv
            w_fst(it)     = msgv
            dfst(it)      = msgv
            mbst(it)      = msgv
            kbst(it)      = msgv
            mbast(it)     = msgv
            kbast(it)     = msgv
            czst(it)      = msgv
            cist(it)      = msgv
            cixst(it)     = msgv
            ciyst(it)     = msgv
            cgxst(it)     = msgv
            cgyst(it)     = msgv
            wgrst(it)     = msgv
            cgrst(it)     = msgv
            cxst(it)      = msgv
            cyst(it)      = msgv
            cpxst(it)     = msgv
            cpyst(it)     = msgv
            cclst(it)     = msgv
            clst(it)      = msgv
            upfrst(it)    = msgv
            dnfrst(it)    = msgv
            rstpsd(1:nrst,it) = msgv
            ruwfst(it,1:nrst) = msgv
            rvwfst(it,1:nrst) = msgv
            uwflst(it,1:nrst) = msgv
            vwflst(it,1:nrst) = msgv
            mruwst(it) = msgv
            mrvwst(it) = msgv
            muwfst(it) = msgv
            mvwfst(it) = msgv
          end if

        else

          phist(it)     = msgv
          degst(it)     = msgv
          degnst(it)    = msgv
          w_fst(it)     = msgv
          dfst(it)      = msgv
          mbst(it)      = msgv
          kbst(it)      = msgv
          mbast(it)     = msgv
          kbast(it)     = msgv
          czst(it)      = msgv
          cist(it)      = msgv
          cixst(it)     = msgv
          ciyst(it)     = msgv
          cgxst(it)     = msgv
          cgyst(it)     = msgv
          wgrst(it)     = msgv
          cgrst(it)     = msgv
          cxst(it)      = msgv
          cyst(it)      = msgv
          cpxst(it)     = msgv
          cpyst(it)     = msgv
          cclst(it)     = msgv
          clst(it)      = msgv
          upfrst(it)    = msgv
          dnfrst(it)    = msgv
          rstpsd(1:nrst,it) = msgv
          ruwfst(it,1:nrst) = msgv
          rvwfst(it,1:nrst) = msgv
          uwflst(it,1:nrst) = msgv
          vwflst(it,1:nrst) = msgv
          mruwst(it) = msgv
          mrvwst(it) = msgv
          muwfst(it) = msgv
          mvwfst(it) = msgv
        end if

        if ( tranal ) then

!----------------------------------------------------------------------------------
!         Inverse Hilbert transform of normalized temperature perturbation
!         in the troposphere
!----------------------------------------------------------------------------------

          do k=1,nrtr
            utrprt1(k) = utrprt(it,k)
            vtrprt1(k) = vtrprt(it,k)
            utrbas1(k) = utrbas(it,k)
            vtrbas1(k) = vtrbas(it,k)
            rtrbas1(k) = rtrbas(it,k)
            tntrprt(k) = ttrprt(it,k)/ttrbas(it,k)
          end do

          call invhil(nrtr  ,dz    ,tntrprt, tntrhil)     

!----------------------------------------------------------------------------------
!         Find wave propagation direction 
!----------------------------------------------------------------------------------

          call propdir(nrtr       ,utrprt1    ,vtrprt1    ,tntrhil    ,  &
                       phitr(it)  ,degtr(it)  ,degntr(it) )        

!----------------------------------------------------------------------------------
!         Calculate Stoke's parameter 
!----------------------------------------------------------------------------------

          call stoke  (nrtr       ,dz         ,utrprt1    ,vtrprt1    ,  &
                       utrbas1    ,vtrbas1    ,degtr(it)  ,avgtrn(it) ,  &
                       w_ftr(it)  ,dftr(it)   )

!----------------------------------------------------------------------------------
!         Calculate mean vertical and horizontal wavenumbers 
!----------------------------------------------------------------------------------

          call meanmk (nrtr       ,dz         ,f          ,w_ftr(it)  ,  &
                       avgtrn(it) ,utrprt1    ,vtrprt1    ,mbtr(it)   ,  & 
                       mbatr(it)  ,kbtr(it)   ,kbatr(it)  )

!----------------------------------------------------------------------------------
!         Calculate phase and group velocities 
!----------------------------------------------------------------------------------

          call phase  (nrtr       ,f          ,w_ftr(it)  ,kbatr(it)  ,  &
                       mbatr(it)  ,degtr(it)  ,utrbas1    ,vtrbas1    ,  & 
                       cztr(it)   ,citr(it)   ,cixtr(it)  ,ciytr(it)  ,  &
                       cgxtr(it)  ,cgytr(it)  ,wgrtr(it)  ,cgrtr(it)  ,  &
                       cxtr(it)   ,cytr(it)   ,cpxtr(it)  ,cpytr(it)  ,  &
                       ubtr(it)   ,ubdgtr(it) )

!----------------------------------------------------------------------------------
!         Estimate a fraction of upward propagating waves by calculating
!         rotary spectra.
!----------------------------------------------------------------------------------

          call rotary (nrtr       ,f          ,dz         ,utrprt1    ,  &
                       vtrprt1    ,wavmtr     ,rtrpsd1    ,ccltr(it)  ,  &
                       cltr(it)   ,upfrtr(it) ,dnfrtr(it) )

          rtrpsd(1:nrtr,it) = rtrpsd1(1:nrtr)

!----------------------------------------------------------------------------------
!         Estimate wave momentum flux 
!----------------------------------------------------------------------------------

          call uwvw   (nrtr       ,avgtrn(it) ,w_ftr(it)  ,f          ,  &
                       utrprt1    ,vtrprt1    ,tntrhil    ,rtrbas1    ,  &
                       ruwftr1    ,rvwftr1    ,uwfltr1    ,vwfltr1    ,  &
                       mruwtr(it) ,mrvwtr(it) ,muwftr(it) ,mvwftr(it) )

          ruwftr(it,1:nrtr) = ruwftr1(1:nrtr)
          rvwftr(it,1:nrtr) = rvwftr1(1:nrtr)
          uwfltr(it,1:nrtr) = uwfltr1(1:nrtr)
          vwfltr(it,1:nrtr) = vwfltr1(1:nrtr)

!----------------------------------------------------------------------------------
!   If mean axial ratio of hodograph is less than one, skip the data
!----------------------------------------------------------------------------------
!----------------------------------------------------------------------------------
!   Following Vincent and Alexander (2000), a small fraction of the frequencies  
!   inferred from soundings had values greater than 10f. These values are not
!   reliable since errors in the Digicora wind measuring sytrem will limit the
!   accuracy to which the transverse velocity perturbation of the hodograph can
!   be determined. Digicora wind speed accuracy is believed to be about 0.4 - 0.5
!   m/s (Ivanov et al. 1991). In Vincent and Alexander (2000), rms wave perturbation
!   is about 4 m/s, which suggetrs that the axial ratios and hence intrinsic    
!   frequencies greater than 10f are not reliable and should be discarded when
!   computing other wave parameters, such as wavelength. 
!----------------------------------------------------------------------------------
 
          if (w_ftr(it) < 1.0 .or. w_ftr(it) > 10.0 ) then
            phitr(it)     = msgv
            degtr(it)     = msgv
            degntr(it)    = msgv
            w_ftr(it)     = msgv
            dftr(it)      = msgv
            mbtr(it)      = msgv
            kbtr(it)      = msgv
            mbatr(it)     = msgv
            kbatr(it)     = msgv
            cztr(it)      = msgv
            citr(it)      = msgv
            cixtr(it)     = msgv
            ciytr(it)     = msgv
            cgxtr(it)     = msgv
            cgytr(it)     = msgv
            wgrtr(it)     = msgv
            cgrtr(it)     = msgv
            cxtr(it)      = msgv
            cytr(it)      = msgv
            cpxtr(it)     = msgv
            cpytr(it)     = msgv
            ccltr(it)     = msgv
            cltr(it)      = msgv
            upfrtr(it)    = msgv
            dnfrtr(it)    = msgv
            rtrpsd(1:nrtr,it) = msgv
            ruwftr(it,1:nrtr) = msgv
            rvwftr(it,1:nrtr) = msgv
            uwfltr(it,1:nrtr) = msgv
            vwfltr(it,1:nrtr) = msgv
            mruwtr(it) = msgv
            mrvwtr(it) = msgv
            muwftr(it) = msgv
            mvwftr(it) = msgv
          end if

        else

          phitr(it)     = msgv
          degtr(it)     = msgv
          degntr(it)    = msgv
          w_ftr(it)     = msgv
          dftr(it)      = msgv
          mbtr(it)      = msgv
          kbtr(it)      = msgv
          mbatr(it)     = msgv
          kbatr(it)     = msgv
          cztr(it)      = msgv
          citr(it)      = msgv
          cixtr(it)     = msgv
          ciytr(it)     = msgv
          cgxtr(it)     = msgv
          cgytr(it)     = msgv
          wgrtr(it)     = msgv
          cgrtr(it)     = msgv
          cxtr(it)      = msgv
          cytr(it)      = msgv
          cpxtr(it)     = msgv
          cpytr(it)     = msgv
          ccltr(it)     = msgv
          cltr(it)      = msgv
          upfrtr(it)    = msgv
          dnfrtr(it)    = msgv
          rtrpsd(1:nrtr,it) = msgv
          ruwftr(it,1:nrtr) = msgv
          rvwftr(it,1:nrtr) = msgv
          uwfltr(it,1:nrtr) = msgv
          vwfltr(it,1:nrtr) = msgv
          mruwtr(it) = msgv
          mrvwtr(it) = msgv
          muwftr(it) = msgv
          mvwftr(it) = msgv
        end if

    end do

!----------------------------------------------------------------------------------
!   Calculate energy weighted wave propagation histogram 
!----------------------------------------------------------------------------------

    call prophist(ndat  ,avgstet,degst ,angs  ,angst  ,mstdir,msgv  )
    call prophist(ndat  ,avgtret,degtr ,angs  ,angtr  ,mtrdir,msgv  )

!----------------------------------------------------------------------------------
!
!   DUMP OUTPUT
!
!----------------------------------------------------------------------------------

    write(6,'(a)') 'DUMP RESULTS'

    call out1d3(trim(wfnst01),'PHIST','DEGST','DEGNST',phist ,degst ,degnst,  &
                'TIME',ndat,normt,'Wave propagation direction')
    call out1d2(trim(wfnst02),'AXRST','DPST',w_fst,dfst,'TIME',ndat,normt,  &
                'Axial ratio and degree of polarization')
    call out1d4(trim(wfnst03),'MBST','MBAST','KBST','KBAST',mbst,mbast,kbst,kbast,  &
                'TIME',ndat,normt,'Mean horizontal and vertical wavenumbers')
    call out1d4(trim(wfnst04),'CIST','CIZST','CIXST','CIYST',cist,czst,cixst,ciyst, &
                'TIME',ndat,normt,'Intrinsic phase velociites')
    call out1d2(trim(wfnst05),'CGXST','CGYST',cgxst,cgyst,'TIME',ndat,normt,  &
                'Group velocities')
    call out1d2(trim(wfnst06),'WGRST','CGRST',wgrst,cgrst,'TIME',ndat,normt,  &
                'Ground-relative frequency and phase speed')
    call out1d4(trim(wfnst07),'CXST','CYST','CPXST','CPYST',cxst,cyst,cpxst,cpyst,  &
                'TIME',ndat,normt,'Ground-relative phase speeds')
    call out1d4(trim(wfnst08),'CLST','CCLST','UPFRST','DNFRST',clst,cclst,upfrst,dnfrst,  &
                'TIME',ndat,normt,'Rotary spectra parameters')
    call out2d (trim(wfnst09),'RSTPSD',rstpsd,'WAVN',nrst,wavmst,'TIME',ndat,normt,  &
                'Rotary spectra')
    call out1d (trim(wfnst10),'ANGST' ,angst ,'DEG' ,12,angs    ,'Angular histogram')
    call out1d (trim(wfnst11),'MDIRST',mstdir,'TIME',1 ,normt(1),'Mean propagation direction')
    call out2d4(trim(wfnst12),'RUWST','RVWST','UWST','VWST',ruwfst,rvwfst,uwflst,vwflst, &
                'TIME',ndat,normt,'Z',nrst,zst,'Momentum flux profiles')
    call out1d4(trim(wfnst13),'MRUWST','MRVWST','MUWST','MVWST',mruwst,mrvwst,muwfst,mvwfst, &
                'TIME',ndat,normt,'Mean momentum flux')

    call out1d3(trim(wfntr01),'PHITR','DEGTR','DEGNTR',phitr ,degtr ,degntr,  &
                'TIME',ndat,normt,'Wave propagation direction')
    call out1d2(trim(wfntr02),'AXRTR','DPTR',w_ftr,dftr,'TIME',ndat,normt,  &
                'Axial ratio and degree of polarization')
    call out1d4(trim(wfntr03),'MBTR','MBATR','KBTR','KBATR',mbtr,mbatr,kbtr,kbatr,  &
                'TIME',ndat,normt,'Mean horizontal and vertical wavenumbers')
    call out1d4(trim(wfntr04),'CITR','CIZTR','CIXTR','CIYTR',citr,cztr,cixtr,ciytr, &
                'TIME',ndat,normt,'Intrinsic phase velociites')
    call out1d2(trim(wfntr05),'CGXTR','CGYTR',cgxtr,cgytr,'TIME',ndat,normt,  &
                'Group velocities')
    call out1d2(trim(wfntr06),'WGRTR','CGRTR',wgrtr,cgrtr,'TIME',ndat,normt,  &
                'Ground-relative frequency and phase speed')
    call out1d4(trim(wfntr07),'CXTR','CYTR','CPXTR','CPYTR',cxtr,cytr,cpxtr,cpytr,  &
                'TIME',ndat,normt,'Ground-relative phase speeds')
    call out1d4(trim(wfntr08),'CLTR','CCLTR','UPFRTR','DNFRTR',cltr,ccltr,upfrtr,dnfrtr,  &
                'TIME',ndat,normt,'Rotary spectra parameters')
    call out2d (trim(wfntr09),'RTRPSD',rtrpsd,'WAVN',nrtr,wavmtr,'TIME',ndat,normt,  &
                'Rotary spectra')
    call out1d (trim(wfntr10),'ANGTR',angtr,'DEG',12,angs,'Angular histogram')
    call out1d (trim(wfnst11),'MDIRTR',mtrdir,'TIME',1 ,normt(1),'Mean propagation direction')
    call out2d4(trim(wfntr12),'RUWTR','RVWTR','UWTR','VWTR',ruwftr,rvwftr,uwfltr,vwfltr, &
                'TIME',ndat,normt,'Z',nrtr,ztr,'Momentum flux profiles')
    call out1d4(trim(wfntr13),'MRUWTR','MRVWTR','MUWTR','MVWTR',mruwtr,mrvwtr,muwftr,mvwftr, &
                'TIME',ndat,normt,'Mean momentum flux')

!----------------------------------------------------------------------------------
!   PLOT
!----------------------------------------------------------------------------------

    write(6,'(a)') 'PLOT'   

    call pltspec(ndat,year2,month,msgv)

!----------------------------------------------------------------------------------
!   Close NCAR GKS
!----------------------------------------------------------------------------------

    call GCLWK(1)
    call GCLWK(2)
    call GCLKS

    stop
    end program specanal


