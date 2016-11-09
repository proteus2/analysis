!2345 
     
      PROGRAM OD

!--------------------------------------------------------------------------------------------
!       
!     Purpose
!      
!     1. To interpolate variables observed using rawinsonde data to user-specified regular 
!     vertical grids, define perturbation and basic-state variables and mean energy density
!     in the user-specified tropospheric and stratospheric data analysis region.
!     ( SONG's program(which program exist /data/sis/radiosonde/waveanal)
!      need to Vaisala data (good and regular). But this program operates well operational
!      and crude data.
!
!     2. 
!
!     Author:
!    
!     JUNG-SUK GOH and IN-SUN SONG
!     Laboratory for Mesoscale Dynamics
!     Department of Atmospheric Sciences, Yonsei University, Seoul, Korea 
!  
!     Strating day  : 2, JUN, 2003
!     Ending   day  :  
!      
!     Required Library
!     
!     NetCDF (higher version than 3.4)
!  
!--------------------------------------------------------------------------------------------

      use netcdfio
      use cubicspline
      use regress
                                
      implicit none
                                 
      Integer, Parameter                          :: nz=501, num=5000
      Real, Parameter                             :: dz=500.0
      Integer                                     :: i,j,k,zz,ww,nn,yy
      Integer                                     :: nlev,index,iargc
      Integer                                     :: pmid,wmid,pid,wid
      Integer                                     :: id,ct,istat
      Real                                        :: wavl,minw,maxw,rd,cp,zsfc,grav
      Real                                        :: MIS,zbuf
      Real, Allocatable, Dimension(:)             :: pmt,pmw,tm,tmw
      Real, Allocatable, Dimension(:)             :: hep,hew,um,vm
      Real, Allocatable, Dimension(:)             :: tempt,tempu,tempv
      Real, Allocatable, Dimension(:)             :: fp,ft,fu,fv
                                     
      Real, Allocatable, Dimension(:)             :: zst,tst,pst,ust,vst
      Real, Allocatable, Dimension(:)             :: rst,thest,bvst
      Real, Allocatable, Dimension(:)             :: tstbas,pstbas,ustbas,vstbas
      Real, Allocatable, Dimension(:)             :: rstbas,thestbas,bvstbas
      Real, Allocatable, Dimension(:)             :: tstprt,pstprt,ustprt,vstprt
                                         
      Real, Allocatable, Dimension(:)             :: zst2,tst2,pst2,ust2,vst2
      Real, Allocatable, Dimension(:)             :: rst2,thest2,bvst2
      Real, Allocatable, Dimension(:)             :: tstbas2,pstbas2,ustbas2,vstbas2
      Real, Allocatable, Dimension(:)             :: rstbas2,thestbas2,bvstbas2
      Real, Allocatable, Dimension(:)             :: tstprt2,pstprt2,ustprt2,vstprt2
                                        
      Real, Allocatable, Dimension(:)             :: ztr,ttr,ptr,utr,vtr
      Real, Allocatable, Dimension(:)             :: rtr,thetr,bvtr
      Real, Allocatable, Dimension(:)             :: ttrbas,ptrbas,utrbas,vtrbas
      Real, Allocatable, Dimension(:)             :: rtrbas,thetrbas,bvtrbas
      Real, Allocatable, Dimension(:)             :: ttrprt,ptrprt,utrprt,vtrprt
                                          
      Real, Allocatable, Dimension(:)             :: ztr2,ttr2,ptr2,utr2,vtr2
      Real, Allocatable, Dimension(:)             :: rtr2,thetr2,bvtr2
      Real, Allocatable, Dimension(:)             :: ttrbas2,ptrbas2,utrbas2,vtrbas2
      Real, Allocatable, Dimension(:)             :: rtrbas2,thetrbas2,bvtrbas2
      Real, Allocatable, Dimension(:)             :: ttrprt2,ptrprt2,utrprt2,vtrprt2
                                    
      Real, Allocatable, Dimension(:)             :: ek,ep,epbas,et,etbas
      Real, Dimension(nz)                         :: fz             
      Character (len=256)                         :: title,wfname
      Character (len=256)                         :: ifname
      Character (len=4)                           :: arg1,arg2,arg3,arg4
      External                                       getarg, iargc  
                                       
!--------------------------------------------------------------------------------------------
! SECOND POLYNOMIAL
!--------------------------------------------------------------------------------------------

      Integer                                     :: nrst,kzbst,kztst,dist
      Integer                                     :: nrst2,kzbst2,kztst2,kstbdz,ksttdz
      Integer                                     :: nrtr,kzbtr,kzttr,ditr
      Integer                                     :: nrtr2,kzbtr2,kzttr2,ktrbdz,ktrtdz
      Real                                        :: zbst,ztst,zbst2,ztst2,stbdz,sttdz
      Real                                        :: zbtr,zttr,zbtr2,zttr2,trbdz,trtdz

!--------------------------------------------------------------------------------------------
! DECLARATION INPUT VARIABLE
!--------------------------------------------------------------------------------------------

      Real                                         :: phi
      Real, Dimension(num)                         :: p,t,w,deg,u,v
 
!--------------------------------------------------------------------------------------------
! VARIABLES FOR GENERALIZED LINEAR REGRESSION
!--------------------------------------------------------------------------------------------

      Integer, Parameter                          :: ma=3, npc=3
      Integer, Dimension(ma)                      :: ia
      Real                                        :: chisq
      Real, Dimension(ma)                         :: a
      Real, Allocatable, Dimension(:)             :: sigst2,sigtr2
      Real, Dimension(npc,npc)                    :: covar
      Real, Dimension(50)                         :: afunc

!--------------------------------------------------------------------------------------------
! INITIATE VARIABLE   
!--------------------------------------------------------------------------------------------
                                    
      nlev=0
      pmid=0 ; wmid=0 ; pid=0 ; wid=0 
      p(:)=-999.0 ; t(:)=-999.0 ; deg(:)=-999.0 ; w(:)=-999.0
      phi=4*atan(1.0) ; MIS=1.e+32 ; minw=1500.0 ; maxw=6000.0
      zsfc=17.79 ; rd=287.0 ; cp=1004.0 ; grav=9.806
                       
!--------------------------------------------------------------------------------------------   
!  INPUT DATA NAME       
!--------------------------------------------------------------------------------------------  
                        
      index = iargc ( )
      if (index /= 3) then
        print *, 'PROGRAM  STOP'
        print *, 'analyod <month> <day> <hour>'
        print *, 'analyod  08 20 12'
!        stop
      else 
        print *, 'The number of arguments is ok'
      end if                                        
      call getarg(1,arg1)
      call getarg(2,arg2)        
      call getarg(3,arg3)        
                        
!--------------------------------------------------------------------------------------------
! OPEN AND READ INPUT DATA
!--------------------------------------------------------------------------------------------
                       
      write(ifname,'(a,a2,a2,a2,a)') &
       '/export5/data_export9/radiosonde/rusa/osan31_100/rawdata/asc/02',arg1,arg2,arg3,'.prn'
      print *, ifname
      open (11,file=ifname,form='formatted')
      do i=1,num
        read (11, '(3x,f5.1,12x,f4.1,13x,f3.1,5x,f3.1)', iostat=istat ) &
           p(i),t(i),deg(i),w(i)
        write(9,*) p(i),t(i),deg(i),w(i)
        if (istat .eq. 0) then
          nlev = nlev +1
        else 
          close(11) 
        end if
      end do   
      print *, nlev, dz
      if ( nlev .gt. 20 ) then      ! IF DATA NUMBER IS SMALLER THAN 20, STOP PROGRAM
        do i=1,nlev
          if (deg(i).eq.-99.9 .or. deg(i).eq.-999.0 .or. deg(i).eq.0.0 ) then
            deg(i)=MIS
          end if
          if (w(i).eq.-99.9 .or. w(i).eq.-999.0 .or. w(i).eq.0.0 ) then
            w(i)=MIS
          end if
          if (p(i).eq.-99.9 .or. p(i).eq.-999.0 .or. p(i).eq.0.0 ) then
            p(i)=MIS
          end if
          if (t(i).eq.-99.9 .or. t(i).eq.-999.0 .or. t(i).eq.0.0 ) then
            t(i)=MIS
          end if
        end do
      print*, deg                    
!-------------------------------------------------------------------------------------------- 
! REMOVE MISSING VALUE
!--------------------------------------------------------------------------------------------
                               
        do i=1,nlev
          if (t(i).ne.MIS ) then
            pmid=pmid+1
          end if
        end do
        zz = pmid                                 ! NEED VARIABLE ALLOCATE
        if ( zz.eq.0 ) then
          print *, 'TEMPERATURE DATA IS NOT EXIST' 
          STOP       
        end if
        allocate(pmt (1:zz)) ; pmt (1:zz) = 0.0    
        allocate(tm (1:zz)) ; tm (1:zz) = 0.0    
        allocate(hep (1:zz)) ; hep (1:zz) = 0.0    
        allocate(tempt (1:zz)) ; tempt (1:zz) = 0.0    
        do i=1,nlev
          if (t(i).ne.MIS) then
            pid = pid+1
            tm(pid)=t(i)
            pmt(pid)=p(i)
          end if
        end do
                       
!-------------------------------------------------------------------------------------------- 
! Because problem of pressure do not exist, don't have to modify 
!--------------------------------------------------------------------------------------------
                                  
!--------------------------------------------------------------------------------------------
! PRESSURE CORRESPONDING WITH WIND
!--------------------------------------------------------------------------------------------
                          
        do i=1,nlev
          if ( w(i).ne.MIS ) then
            wmid=wmid+1
          end if
        end do
        ww = wmid
        if ( ww.eq.0 ) then
          print *, 'WIND DATA IS NOT EXIST' 
          STOP       
        end if
        allocate(pmw (1:ww)) ; pmw (1:ww) = 0.0
        allocate(um (1:ww)) ; um (1:ww) = 0.0
        allocate(vm (1:ww)) ; vm (1:ww) = 0.0
        allocate(hew (1:ww)) ; hew (1:ww) = 0.0
        allocate(tmw (1:ww)) ; tmw (1:ww) = 0.0
        allocate(tempu (1:ww)) ; tempu (1:ww) = 0.0 
        allocate(tempv (1:ww)) ; tempv (1:ww) = 0.0
        do i=1,nlev
          if ( w(i).ne.MIS ) then
            wid = wid+1
            w(wid)=w(i)
            deg(wid)=deg(i) 
            um(wid)=-w(i)*10.*0.5*(sin(deg(i)*phi*10./180.))
            vm(wid)=-w(i)*10.*0.5*(cos(deg(i)*phi*10./180.))
            pmw(wid)=p(i)
          end if
        end do
               
!--------------------------------------------------------------------------------------------
! CHECK PRESSURE IS RIGHT
!-------------------------------------------------------------------------------------------- 
             
        do j=1,zz-1
          if (pmt(j) .le. pmt(j+1)) then
            print *, 'PRESSURE pmt IS NOT RIGHT'
            STOP   
          end if
        end do
        do j=1,ww-1
          if (pmw(j) .le. pmw(j+1)) then
            print *, 'PRESSURE pmw IS NOT RIGHT'
            STOP   
          end if
        end do  
                   
!--------------------------------------------------------------------------------------------
! CALCULATE HEIGHT CORRESPONDING TEMPERATURE
!--------------------------------------------------------------------------------------------
                    
        do i=1,zz
          hep(i) = (287.0*(tm(i)+273.15)/9.806*log(pmt(1)/pmt(i))) + zsfc
        end do
                 
!--------------------------------------------------------------------------------------------
! INTERPOLATION TEMPERATURE
!--------------------------------------------------------------------------------------------
                  
        do k=1,nz
          fz(k) = (k-1)*dz
        end do
                     
        nn = int( hep(zz) / dz )
                 
        allocate(fp (1:nn)) ; fp (1:nn)=0.0    
        allocate(ft (1:nn)) ; ft (1:nn)=0.0    
        allocate(ek (1:nn)) ; ek(1:nn)=0.0
        allocate(ep (1:nn)) ; ep(1:nn)=0.0
        allocate(epbas (1:nn)) ; epbas(1:nn)=0.0
        allocate(et (1:nn)) ; et(1:nn)=0.0
        allocate(etbas (1:nn)) ; etbas(1:nn)=0.0
                                      
        do k=1,zz
          tempt(k) = 0.0
        end do
                                  
        do k=1,zz
          tm(k)=tm(k)+273.15
        end do
                            
        call spline(hep,tm,zz,1.e+32,1.e+32,tempt)
        do k=1,nn 
          call splint(hep,tm,tempt,zz,fz(k),ft(k))  
        end do
                                                    
        do k=1,nn
          if (fz(k) .gt. 10000.0) then
            if (ft(k) .gt. 270.0 ) then
              print *, 'C-S IS ERROR'
              print *, 'T IS OVERESTIMATED' 
              STOP
            end if
          end if
        end do
                      
!--------------------------------------------------------------------------------------------
! ESTIMATE PRESSURE CORRESPONDING WITH HEIGHT 
!-------------------------------------------------------------------------------------------- 
                         
        do k=1,nn
          fp(k) = pmt(1)*exp( -(fz(k)-zsfc) / (287.0*ft(k)/9.806) )  
        end do 
                   
!--------------------------------------------------------------------------------------------
! EVALUATE HEIGHT CORRESPODING WIND PRESSURE
!--------------------------------------------------------------------------------------------
                         
        do j=1,ww
          if ( pmw(j).gt.fp(1) ) then
            tmw(j) = ft(1)
          end if
          if ( pmw(j).le.fp(nn) )  then
            tmw(j) = ft(nn)
          end if
        end do       
        do j=1,ww
          do k=1,nn-1
            if ( pmw(j).le.fp(k) .and. pmw(j).gt.fp(k+1) ) then
              tmw(j) = ft(k)
            end if
          end do
        end do
                                
        do j=1,ww
          hew(j) = (287.0*tmw(j)/9.806*log(pmw(1)/pmw(j))) + zsfc
        end do
                           
!--------------------------------------------------------------------------------------------
! INTERPOLATION UWIND,VWIND
!--------------------------------------------------------------------------------------------
                  
        yy = int( hew(ww) / dz )
                       
        allocate(fu (1:yy)) ; fu (1:yy)=0.0    
        allocate(fv (1:yy)) ; fv (1:yy)=0.0    
                        
        do j=1,ww 
          tempu(j) = 0.0
          tempv(j) = 0.0
        end do 
                  
        call spline(hew,um,ww,1.e+32,1.e+32,tempu)
        call spline(hew,vm,ww,1.e+32,1.e+32,tempv)
                  
        do k=1,yy
          call splint(hew,um,tempu,ww,fz(k),fu(k))
          call splint(hew,vm,tempv,ww,fz(k),fv(k))
        end do
                     
        print *, 'TEMPERATURE AND WIND DATA NUMBER', nn, yy
                   
!--------------------------------------------------------------------------------------------
! OUTPUT ORIGINAL AND INTERPOLATED RAWDATA
!-------------------------------------------------------------------------------------------- 
                        
        write(title,'(a)') 'Original Temperature '
        write(wfname, '(a,a2,a2,a2,a)') &
             '/export9/radiosonde/rusa/osan24_500/rawdata/Net/T/ort',arg1,arg2,arg3,'.nc'
        call out1d(trim(wfname),'T',tm,'Z',zz,hep,trim(title))
                                
        write(title,'(a)') 'Original U-WIND'
        write(wfname, '(a,a2,a2,a2,a)') & 
             '/export9/radiosonde/rusa/osan24_500/rawdata/Net/U/oru',arg1,arg2,arg3,'.nc'
        call out1d(trim(wfname),'U',um,'Z',ww,hew,trim(title))
                           
        write(title,'(a)') 'Oroginal V-WIND'
        write(wfname, '(a,a2,a2,a2,a)')  &
             '/export9/radiosonde/rusa/osan24_500/rawdata/Net/V/orv',arg1,arg2,arg3,'.nc'
        call out1d(trim(wfname),'V',vm,'Z',ww,hew,trim(title))
                          
        write(title,'(a)') ' Interpolated Temperature '
        write(wfname, '(a,a2,a2,a2,a)') &
             '/export9/radiosonde/rusa/osan24_500/rawdata/Net/T/int',arg1,arg2,arg3,'.nc'
        call out1d(trim(wfname),'T',ft,'Z',nn,fz,trim(title))
                                
        write(title,'(a)') 'Interpolated U-WIND'
        write(wfname, '(a,a2,a2,a2,a)') & 
             '/export9/radiosonde/rusa/osan24_500/rawdata/Net/U/inu',arg1,arg2,arg3,'.nc'
        call out1d(trim(wfname),'U',fu,'Z',yy,fz,trim(title))
                           
        write(title,'(a)') 'Intepolated V-WIND'
        write(wfname, '(a,a2,a2,a2,a)')  &
             '/export9/radiosonde/rusa/osan24_500/rawdata/Net/V/inv',arg1,arg2,arg3,'.nc'
        call out1d(trim(wfname),'V',fv,'Z',yy,fz,trim(title))
                       
!--------------------------------------------------------------------------------------------
! JURDGE PROPER TEMPERATURE AND WIND HEIGHT
!--------------------------------------------------------------------------------------------
               
!        if (hep(zz) .lt. 31000.) then
        if (hep(zz) .lt. 24500.) then
          print *, 'TEM. HEIGHT IS SMALLER THAN CRITICAL VALUE'
          stop
        else 
          print *, 'TEM. PERTURBATION HEIGHT SATISFY'
        end if
                   
!        if (hew(ww) .lt. 31000.) then
        if (hew(ww) .lt. 24500.) then
          print *, 'WIND. HEIGHT IS SMALLER THAN CRITICAL VALUE'
          stop
        else 
          print *, 'WIND. PERTURBATION HEIGHT SATISFY'
        end if
                    
!--------------------------------------------------------------------------------------------
! START ANALYSIS : 1. DEFINE PERTURBATION AND OBTAIN MEAN VALUES 
!-------------------------------------------------------------------------------------------- 
                           
!--------------------------------------------------------------------------------------------
! DEFINE HEIGHT RANGE AND OBTAIN TOTAL T, P, U, V, RHO 
!--------------------------------------------------------------------------------------------
                       
!--------------------------------------------------------------------------------------------
! FIRST, GET VARIABLES OF STRATOSPHERE
!--------------------------------------------------------------------------------------------
                                          
        zbuf   = 500.
               
!        zbst   = 17500.
!        ztst   = 30500.
               
        zbst   = 17000.
        ztst   = 24000.
 
        zbst2  = zbst - zbuf
        ztst2  = ztst + zbuf
                   
        stbdz  = zbst - zbst2   
        sttdz  = ztst2 - ztst
                
        kstbdz = int(stbdz/dz + 0.001)
        ksttdz = int(sttdz/dz + 0.001)
                                                
        kzbst  = int(zbst/dz + 0.001) 
        kztst  = int(ztst/dz + 0.001)
                 
        kzbst2 = int(zbst2/dz + 0.001)
        kztst2 = int(ztst2/dz + 0.001)
                 
        nrst   = kztst - kzbst   + 1
        nrst2  = kztst2 - kzbst2 + 1
         
        print *, 'nrst',nrst,'nrst2',nrst2  
            
        allocate(zst(1:nrst))       ; zst(:)=0.0
        allocate(tst(1:nrst))       ; tst(:)=0.0
        allocate(pst(1:nrst))       ; pst(:)=0.0
        allocate(ust(1:nrst))       ; ust(:) = 0.0
        allocate(vst(1:nrst))       ; vst(:) = 0.0
        allocate(rst(1:nrst))       ; rst(:) = 0.0
        allocate(thest(1:nrst))     ; thest(:) = 0.0
        allocate(bvst(1:nrst))      ; bvst(:) = 0.0
              
        allocate(tstbas(1:nrst))   ; tstbas(:) = 0.0
        allocate(pstbas(1:nrst))   ; pstbas(:) = 0.0
        allocate(ustbas(1:nrst))   ; ustbas(:) = 0.0
        allocate(vstbas(1:nrst))   ; vstbas(:) = 0.0
        allocate(rstbas(1:nrst))   ; rstbas(:) = 0.0
        allocate(thestbas(1:nrst)) ; thestbas(:) = 0.0
        allocate(bvstbas(1:nrst))  ; bvstbas(:) = 0.0
                 
        allocate(tstprt(1:nrst))   ; tstprt(:) = 0.0
        allocate(pstprt(1:nrst))   ; pstprt(:) = 0.0
        allocate(ustprt(1:nrst))   ; ustprt(:) = 0.0
        allocate(vstprt(1:nrst))   ; vstprt(:) = 0.0

        allocate(zst2(1:nrst2))     ; zst2(:) = 0.0
        allocate(tst2(1:nrst2))     ; tst2(:) = 0.0
        allocate(pst2(1:nrst2))     ; pst2(:) = 0.0
        allocate(ust2(1:nrst2))     ; ust2(:) = 0.0
        allocate(vst2(1:nrst2))     ; vst2(:) = 0.0
        allocate(rst2(1:nrst2))     ; rst2(:) = 0.0
        allocate(thest2(1:nrst2))   ; thest2(:) = 0.0
        allocate(bvst2(1:nrst2))    ; bvst2(:) = 0.0
           
        allocate(tstbas2(1:nrst2))   ; tstbas2(:) = 0.0
        allocate(pstbas2(1:nrst2))   ; pstbas2(:) = 0.0
        allocate(ustbas2(1:nrst2))   ; ustbas2(:) = 0.0
        allocate(vstbas2(1:nrst2))   ; vstbas2(:) = 0.0
        allocate(rstbas2(1:nrst2))   ; rstbas2(:) = 0.0
        allocate(thestbas2(1:nrst2)) ; thestbas2(:) = 0.0
        allocate(bvstbas2(1:nrst2))  ; bvstbas2(:) = 0.0
             
        allocate(tstprt2(1:nrst2))   ; tstprt2(:) = 0.0
        allocate(pstprt2(1:nrst2))   ; pstprt2(:) = 0.0
        allocate(ustprt2(1:nrst2))   ; ustprt2(:) = 0.0
        allocate(vstprt2(1:nrst2))   ; vstprt2(:) = 0.0
                           
        allocate(sigst2(1:nrst2))    ; sigst2(:)=0.0
                        
        do k=1,nrst
          zst(k) = (k-1)*dz + zbst
        end do  
                   
        do k=1,nrst2
          zst2(k) = (k-1)*dz + zbst2 
          tst2(k) = ft(k-1+kzbst2)
          pst2(k) = fp(k-1+kzbst2)
          ust2(k) = fu(k-1+kzbst2)
          vst2(k) = fv(k-1+kzbst2)
!          print *, k, ust2(k)
        end do       
        do k=1,nrst2
          rst2(k) = 100.0*pst2(k)/(rd*tst2(k))
        end do     
      
        PRINT *, 'ZST2 MAXVAL:: ', maxval(zst2) 
        PRINT *, 'ZST2 MINVAL:: ', minval(zst2) 
                     
        PRINT *, 'TST2 MAXVAL:: ', maxval(tst2) 
        PRINT *, 'TST2 MINVAL:: ', minval(tst2) 
                              
        PRINT *, 'PST2 MAXVAL:: ', maxval(pst2) 
        PRINT *, 'PST2 MINVAL:: ', minval(pst2) 
                                                   
        PRINT *, 'UST2 MAXVAL:: ', maxval(ust2) 
        PRINT *, 'UST2 MINVAL:: ', minval(ust2) 
                           
        PRINT *, 'VST2 MAXVAL:: ', maxval(vst2) 
        PRINT *, 'VST2 MINVAL:: ', minval(vst2) 
                       
!--------------------------------------------------------------------------------------------
! SECOND, GET VARIABLES OF TROPOSPHERE
!--------------------------------------------------------------------------------------------
                                  
        zbtr   = 2500.
        zttr   = 9500.
              
        zbtr2  = zbtr - zbuf
        zttr2  = zttr + zbuf
                             
        trbdz  = zbtr  - zbtr2
        trtdz  = zttr2 - zttr
                             
        ktrbdz = int(trbdz/dz + 0.001)
        ktrtdz = int(trtdz/dz + 0.001)
                              
        kzbtr  = int(zbtr/dz + 0.001)
        kzttr  = int(zttr/dz + 0.001) 
                             
        kzbtr2 = int(zbtr2/dz + 0.001)
        kzttr2 = int(zttr2/dz + 0.001) 
                          
        nrtr   = kzttr - kzbtr   +1
        nrtr2  = kzttr2 - kzbtr2 +1
                     
        allocate(ztr(1:nrtr))       ; ztr(:)=0.0
        allocate(ttr(1:nrtr))       ; ttr(:)=0.0
        allocate(ptr(1:nrtr))       ; ptr(:)=0.0
        allocate(utr(1:nrtr))       ; utr(:) = 0.0
        allocate(vtr(1:nrtr))       ; vtr(:) = 0.0
        allocate(rtr(1:nrtr))       ; rtr(:) = 0.0
        allocate(thetr(1:nrtr))     ; thetr(:) = 0.0
        allocate(bvtr(1:nrtr))      ; bvtr(:) = 0.0
  
        allocate(ttrbas(1:nrtr))   ; ttrbas(:) = 0.0
        allocate(ptrbas(1:nrtr))   ; ptrbas(:) = 0.0
        allocate(utrbas(1:nrtr))   ; utrbas(:) = 0.0
        allocate(vtrbas(1:nrtr))   ; vtrbas(:) = 0.0
        allocate(rtrbas(1:nrtr))   ; rtrbas(:) = 0.0
        allocate(thetrbas(1:nrtr)) ; thetrbas(:) = 0.0
        allocate(bvtrbas(1:nrtr))  ; bvtrbas(:) = 0.0
                 
        allocate(ttrprt(1:nrtr))   ; ttrprt(:) = 0.0
        allocate(ptrprt(1:nrtr))   ; ptrprt(:) = 0.0
        allocate(utrprt(1:nrtr))   ; utrprt(:) = 0.0
        allocate(vtrprt(1:nrtr))   ; vtrprt(:) = 0.0
                            
        allocate(ztr2(1:nrtr2))      ; ztr2(:) = 0.0
        allocate(ttr2(1:nrtr2))      ; ttr2(:) = 0.0
        allocate(ptr2(1:nrtr2))      ; ptr2(:) = 0.0
        allocate(utr2(1:nrtr2))      ; utr2(:) = 0.0
        allocate(vtr2(1:nrtr2))      ; vtr2(:) = 0.0
        allocate(rtr2(1:nrtr2))      ; rtr2(:) = 0.0
        allocate(thetr2(1:nrtr2))    ; thetr2(:) = 0.0
        allocate(bvtr2(1:nrtr2))     ; bvtr2(:) = 0.0
                 
        allocate(ttrbas2(1:nrtr2))   ; ttrbas2(:) = 0.0
        allocate(ptrbas2(1:nrtr2))   ; ptrbas2(:) = 0.0
        allocate(utrbas2(1:nrtr2))   ; utrbas2(:) = 0.0
        allocate(vtrbas2(1:nrtr2))   ; vtrbas2(:) = 0.0
        allocate(rtrbas2(1:nrtr2))   ; rtrbas2(:) = 0.0
        allocate(thetrbas2(1:nrtr2)) ; thetrbas2(:) = 0.0
        allocate(bvtrbas2(1:nrtr2))  ; bvtrbas2(:) = 0.0
                 
        allocate(ttrprt2(1:nrtr2))   ; ttrprt2(:) = 0.0
        allocate(ptrprt2(1:nrtr2))   ; ptrprt2(:) = 0.0
        allocate(utrprt2(1:nrtr2))   ; utrprt2(:) = 0.0
        allocate(vtrprt2(1:nrtr2))   ; vtrprt2(:) = 0.0
                            
        allocate(sigtr2(1:nrtr2))    ; sigtr2(:)=0.0
         
        do k=1,nrtr
          ztr(k) = (k-1)*dz + zbtr
        end do  
        do k=1,nrtr2
          ztr2(k) = (k-1)*dz + zbtr2 
          ttr2(k) = ft(k-1+kzbtr2)
          ptr2(k) = fp(k-1+kzbtr2)
          utr2(k) = fu(k-1+kzbtr2)
          vtr2(k) = fv(k-1+kzbtr2)
        end do       
        do k=1,nrtr2
          rtr2(k) = 100.0*ptr2(k)/(rd*ttr2(k))
        end do     
                    
        PRINT *, 'ZTR2 MAXVAL:: ', maxval(ztr2) 
        PRINT *, 'ZTR2 MINVAL:: ', minval(ztr2) 
                     
        PRINT *, 'TTR2 MAXVAL:: ', maxval(ttr2) 
        PRINT *, 'TTR2 MINVAL:: ', minval(ttr2) 
                              
        PRINT *, 'PTR2 MAXVAL:: ', maxval(ptr2) 
        PRINT *, 'PTR2 MINVAL:: ', minval(ptr2) 
                                                   
        PRINT *, 'UTR2 MAXVAL:: ', maxval(utr2) 
        PRINT *, 'UTR2 MINVAL:: ', minval(utr2) 
                           
        PRINT *, 'VTR2 MAXVAL:: ', maxval(vtr2) 
        PRINT *, 'VTR2 MINVAL:: ', minval(vtr2) 
                                                                                
!--------------------------------------------------------------------------------------------
! OBTAIN 2nd POLINOMIAL MEAN AND PERTURBATION IN THE STRATOSPHERE
!--------------------------------------------------------------------------------------------
                           
        ia(1:ma)        = 1
        a(1:ma)         = 1.0
        sigst2(1:nrst2) = 1.0
        call lfit(zst2,tst2,sigst2,nrst2,a,ia,ma,covar,npc,chisq)
        do k=1,nrst2
          tstbas2(k) = 0.0
          call funcs(zst2(k),afunc,ma)
          do i=1,ma
            tstbas2(k)=tstbas2(k)+a(i)*afunc(i)
          end do
          tstprt2(k)=tst2(k)-tstbas2(k)
        end do 
                
        ia(1:ma)        = 1
        a(1:ma)         = 1.0
        sigst2(1:nrst2) = 1.0
        call lfit(zst2,pst2,sigst2,nrst2,a,ia,ma,covar,npc,chisq)
        do k=1,nrst2
          pstbas2(k) = 0.0
          call funcs(zst2(k),afunc,ma)
          do i=1,ma
            pstbas2(k)=pstbas2(k)+a(i)*afunc(i)
          end do
          pstprt2(k)=pst2(k)-pstbas2(k)
        end do 
                        
        ia(1:ma)        = 1
        a(1:ma)         = 1.0
        sigst2(1:nrst2) = 1.0
        call lfit(zst2,ust2,sigst2,nrst2,a,ia,ma,covar,npc,chisq)
        do k=1,nrst2
          ustbas2(k) = 0.0
          call funcs(zst2(k),afunc,ma)
          do i=1,ma
            ustbas2(k)=ustbas2(k)+a(i)*afunc(i)
          end do
          ustprt2(k)=ust2(k)-ustbas2(k)
        end do 
                             
        ia(1:ma)        = 1
        a(1:ma)         = 1.0
        sigst2(1:nrst2) = 1.0
        call lfit(zst2,vst2,sigst2,nrst2,a,ia,ma,covar,npc,chisq)
        do k=1,nrst2
          vstbas2(k) = 0.0
          call funcs(zst2(k),afunc,ma)
          do i=1,ma
            vstbas2(k)=vstbas2(k)+a(i)*afunc(i)
          end do
          vstprt2(k)=vst2(k)-vstbas2(k)
        end do 
                          
        do k=1,nrst2
          thest2(k) = tst2(k)*(1000./pst2(k))**(rd/cp)
          thestbas2(k) = tstbas2(k)*(1000./pstbas2(k))**(rd/cp)
          rstbas2(k) = 100.0*pstbas2(k)/(rd*tstbas2(k))
        end do            
        do k=2,nrst2-1
          bvst2(k)    = grav*((tst2(k+1)-tst2(k-1))/(2.0*dz)+grav/cp) &
                        / ((tst2(k+1)+tst2(k-1))*0.5)            
          bvstbas2(k) = grav*((tstbas2(k+1)-tstbas2(k-1))/(2.0*dz)+grav/cp) &
                        / ((tstbas2(k+1)+tstbas2(k-1))*0.5)
        end do
        do k=2,nrst2-1
          if (bvst2(k) .lt. 0.0) then
            bvst2(k) = 0.0
          else
            bvst2(k) = sqrt(bvst2(k))
          end if
          if (bvstbas2(k) .lt. 0.0) then
            bvstbas2(k) = 0.0
          else 
            bvstbas2(k) = sqrt(bvstbas2(k))
          end if
        end do
        bvst2(1)        = bvst2(2)         
        bvst2(nrst2)    = bvst2(nrst2-1)
        bvstbas2(1)     = bvstbas2(2)       
        bvstbas2(nrst2) = bvstbas2(nrst2-1)
  
        print *, 'TST2 MAXVAL:: ' , maxval(tstbas2(:))  
        print *, 'TST2 MINVAL:: ' , minval(tstbas2(:))  
                                      
        print *, 'PST2 MAXVAL:: ' , maxval(pstbas2(:))  
        print *, 'PST2 MINVAL:: ' , minval(pstbas2(:))  
                                      
        print *, 'UST2 MAXVAL:: ' , maxval(ustbas2(:))  
        print *, 'UST2 MINVAL:: ' , minval(ustbas2(:))  
                                   
        print *, 'VST2 MAXVAL:: ' , maxval(vstbas2(:))  
        print *, 'VST2 MINVAL:: ' , minval(vstbas2(:))  
                                    
        print *, 'THEST2 MAXVAL:: ' , maxval(thestbas2(:))  
        print *, 'THEST2 MINVAL:: ' , minval(thestbas2(:))  
                              
        print *, 'BVST2 MAXVAL:: ' , maxval(bvstbas2(:))  
        print *, 'BVST2 MINVAL:: ' , minval(bvstbas2(:))  
                                     
!--------------------------------------------------------------------------------------------
! OBTAIN 2nd POLINOMIAL MEAN AND PERTURBATION IN THE TROPOSPHERE
!--------------------------------------------------------------------------------------------
                              
        ia(1:ma)        = 1
        a(1:ma)         = 1.0
        sigtr2(1:nrtr2) = 1.0
        call lfit(ztr2,ttr2,sigtr2,nrtr2,a,ia,ma,covar,npc,chisq)
        do k=1,nrtr2
          ttrbas2(k) = 0.0
          call funcs(ztr2(k),afunc,ma)
          do i=1,ma
            ttrbas2(k)=ttrbas2(k)+a(i)*afunc(i)
          end do
          ttrprt2(k)=ttr2(k)-ttrbas2(k)
        end do 
                
        ia(1:ma)        = 1
        a(1:ma)         = 1.0
        sigtr2(1:nrtr2) = 1.0
        call lfit(ztr2,ptr2,sigtr2,nrtr2,a,ia,ma,covar,npc,chisq)
        do k=1,nrtr2
          ptrbas2(k) = 0.0
          call funcs(ztr2(k),afunc,ma)
          do i=1,ma
            ptrbas2(k)=ptrbas2(k)+a(i)*afunc(i)
          end do
          ptrprt2(k)=ptr2(k)-ptrbas2(k)
        end do 
                           
        ia(1:ma)        = 1
        a(1:ma)         = 1.0
        sigtr2(1:nrtr2) = 1.0
        call lfit(ztr2,utr2,sigtr2,nrtr2,a,ia,ma,covar,npc,chisq)
        do k=1,nrtr2
          utrbas2(k) = 0.0
          call funcs(ztr2(k),afunc,ma)
          do i=1,ma
            utrbas2(k)=utrbas2(k)+a(i)*afunc(i)
          end do
          utrprt2(k)=utr2(k)-utrbas2(k)
        end do 
                          
        ia(1:ma)        = 1
        a(1:ma)         = 1.0
        sigtr2(1:nrtr2) = 1.0
        call lfit(ztr2,vtr2,sigtr2,nrtr2,a,ia,ma,covar,npc,chisq)
        do k=1,nrtr2
          vtrbas2(k) = 0.0
          call funcs(ztr2(k),afunc,ma)
          do i=1,ma
            vtrbas2(k)=vtrbas2(k)+a(i)*afunc(i)
          end do
          vtrprt2(k)=vtr2(k)-vtrbas2(k)
        end do 
                         
        do k=1,nrtr2
          thetr2(k)    = ttr2(k)*(1000./ptr2(k))**(rd/cp)
          thetrbas2(k) = ttrbas2(k)*(1000./ptrbas2(k))**(rd/cp)
          rtrbas2(k) = 100.0*ptrbas2(k)/(rd*ttrbas2(k))
        end do            
        do k=2,nrtr2-1
          bvtr2(k)    = grav*((ttr2(k+1)-ttr2(k-1))/(2.0*dz)+grav/cp) &
                        / ((ttr2(k+1)+ttr2(k-1))*0.5)            
          bvtrbas2(k) = grav*((ttrbas2(k+1)-ttrbas2(k-1))/(2.0*dz)+grav/cp) &
                        / ((ttrbas2(k+1)+ttrbas2(k-1))*0.5)
        end do
        do k=2,nrtr2-1
          if (bvtr2(k) .lt. 0.0) then
            bvtr2(k) = 0.0
          else
            bvtr2(k) = sqrt(bvtr2(k))
          end if
          if (bvtrbas2(k) .lt. 0.0) then
            bvtrbas2(k) = 0.0
          else 
            bvtrbas2(k) = sqrt(bvtrbas2(k))
          end if
        end do
        bvtr2(1)        = bvtr2(2)         
        bvtr2(nrtr2)    = bvtr2(nrtr2-1)
        bvtrbas2(1)     = bvtrbas2(2)       
        bvtrbas2(nrtr2) = bvtrbas2(nrtr2-1)
                                             
        print *, 'TTR2 MAXVAL:: ' , maxval(ttrbas2(:))  
        print *, 'TTR2 MINVAL:: ' , minval(ttrbas2(:))  
                                      
        print *, 'PTR2 MAXVAL:: ' , maxval(ptrbas2(:))  
        print *, 'PTR2 MINVAL:: ' , minval(ptrbas2(:))  
                                      
        print *, 'UTR2 MAXVAL:: ' , maxval(utrbas2(:))  
        print *, 'UTR2 MINVAL:: ' , minval(utrbas2(:))  
                                   
        print *, 'VTR2 MAXVAL:: ' , maxval(vtrbas2(:))  
        print *, 'VTR2 MINVAL:: ' , minval(vtrbas2(:))  
                                    
        print *, 'THETR2 MAXVAL:: ' , maxval(thetrbas2(:))  
        print *, 'THETR2 MINVAL:: ' , minval(thetrbas2(:))  
                              
        print *, 'BVTR2 MAXVAL:: ' , maxval(bvtrbas2(:))  
        print *, 'BVTR2 MINVAL:: ' , minval(bvtrbas2(:))  
                                     
!--------------------------------------------------------------------------------------------
! CHANGE WITH EXACT HEIGHT RANGE OF STRATOSPHERE AND TROPOSPHERE
!--------------------------------------------------------------------------------------------
                  
        do k=1,nrst2
          if (zst2(k) .eq. zst(1)) then
            dist = k
          end if
        end do
        print *, 'nrst=', nrst,',nrst2=',nrst2  
 
        do k=1,nrst
          tst(k)      = tst2(k+dist-1)      
          pst(k)      = pst2(k+dist-1)
          ust(k)      = ust2(k+dist-1)
          vst(k)      = vst2(k+dist-1)
          rst(k)      = rst2(k+dist-1)
          thest(k)    = thest2(k+dist-1)
          bvst(k)     = bvst2(k+dist-1)
          tstbas(k)   = tstbas2(k+dist-1)
          pstbas(k)   = pstbas2(k+dist-1)
          ustbas(k)   = ustbas2(k+dist-1)
          vstbas(k)   = vstbas2(k+dist-1)
          rstbas(k)   = rstbas2(k+dist-1)
          thestbas(k) = thestbas2(k+dist-1)
          bvstbas(k)  = bvstbas2(k+dist-1)
          tstprt(k)   = tstprt2(k+dist-1)
          pstprt(k)   = pstprt2(k+dist-1)
          ustprt(k)   = ustprt2(k+dist-1)
          vstprt(k)   = vstprt2(k+dist-1)
        end do
   
        do k=1,nrtr2
          if (ztr2(k) .eq. ztr(1)) then
            ditr = k
          end if
        end do
        do k=1,nrtr
          ttr(k)   = ttr2(k+ditr-1)
          ptr(k)   = ptr2(k+ditr-1)
          utr(k)   = utr2(k+ditr-1)
          vtr(k)   = vtr2(k+ditr-1)
          rtr(k)   = rtr2(k+ditr-1)
          thetr(k) = thetr2(k+ditr-1)
          bvtr(k)  = bvtr2(k+ditr-1)
          ttrbas(k)   = ttrbas2(k+ditr-1)
          ptrbas(k)   = ptrbas2(k+ditr-1)
          utrbas(k)   = utrbas2(k+ditr-1)
          vtrbas(k)   = vtrbas2(k+ditr-1)
          rtrbas(k)   = rtrbas2(k+ditr-1)
          thetrbas(k) = thetrbas2(k+ditr-1)
          bvtrbas(k)  = bvtrbas2(k+ditr-1)
          ttrprt(k)   = ttrprt2(k+ditr-1)
          ptrprt(k)   = ptrprt2(k+ditr-1)
          utrprt(k)   = utrprt2(k+ditr-1)
          vtrprt(k)   = vtrprt2(k+ditr-1)
        end do   
                           
!--------------------------------------------------------------------------------------------
! OUPUT DATA OF STRATOSPHERE
!--------------------------------------------------------------------------------------------
                                 
        write(title, '(a)') 'TOTAL TEMPERATURE' 
        write(wfname, '(a,a2,a2,a2,a)') &
             '/export9/radiosonde/rusa/osan24_500/T/tsttot',arg1,arg2,arg3,'.nc'
        call out1d(trim(wfname),'T',tst,'Z',nrst,zst,trim(title))

        write(title, '(a)') 'TOTAL PRESSURE' 
        write(wfname, '(a,a2,a2,a2,a)') &
            '/export9/radiosonde/rusa/osan24_500/P/psttot',arg1,arg2,arg3,'.nc'
        call out1d(trim(wfname),'P',pst,'Z',nrst,zst,trim(title))
    
        write(title, '(a)') 'TOTAL U-WIND' 
        write(wfname, '(a,a2,a2,a2,a)') &
             '/export9/radiosonde/rusa/osan24_500/U/usttot',arg1,arg2,arg3,'.nc'
        call out1d(trim(wfname),'U',ust,'Z',nrst,zst,trim(title))
                       
        write(title, '(a)') 'TOTAL V-WIND' 
        write(wfname, '(a,a2,a2,a2,a)') &
             '/export9/radiosonde/rusa/osan24_500/V/vsttot',arg1,arg2,arg3,'.nc'
        call out1d(trim(wfname),'V',vst,'Z',nrst,zst,trim(title))
                   
        write(title, '(a)') 'TOTAL RHO' 
        write(wfname, '(a,a2,a2,a2,a)') &
             '/export9/radiosonde/rusa/osan24_500/RHO/rsttot',arg1,arg2,arg3,'.nc'
        call out1d(trim(wfname),'RHO',rst,'Z',nrst,zst,trim(title))
                          
        write(title, '(a)') 'TOTAL THETA' 
        write(wfname, '(a,a2,a2,a2,a)') &
             '/export9/radiosonde/rusa/osan24_500/THE/thesttot',arg1,arg2,arg3,'.nc'
        call out1d(trim(wfname),'THE',thest,'Z',nrst,zst,trim(title))
                                 
        write(title, '(a)') 'TOTAL BVF' 
        write(wfname, '(a,a2,a2,a2,a)') &
             '/export9/radiosonde/rusa/osan24_500/BV/bvsttot',arg1,arg2,arg3,'.nc'
        call out1d(trim(wfname),'BV',bvst,'Z',nrst,zst,trim(title))
        
        write(title, '(a)') 'BASIC TEMPERATURE' 
        write(wfname, '(a,a2,a2,a2,a)') &
             '/export9/radiosonde/rusa/osan24_500/T/tstbas',arg1,arg2,arg3,'.nc'
        call out1d(trim(wfname),'TBAS',tstbas,'Z',nrst,zst,trim(title))
                                                              
        write(title, '(a)') 'BASIC PRESSURE' 
        write(wfname, '(a,a2,a2,a2,a)') &
             '/export9/radiosonde/rusa/osan24_500/P/pstbas',arg1,arg2,arg3,'.nc'
        call out1d(trim(wfname),'PBAS',pstbas,'Z',nrst,zst,trim(title))
                                      
        write(title, '(a)') 'BASIC U-WIND' 
        write(wfname, '(a,a2,a2,a2,a)') &
             '/export9/radiosonde/rusa/osan24_500/U/ustbas',arg1,arg2,arg3,'.nc'
        call out1d(trim(wfname),'UBAS',ustbas,'Z',nrst,zst,trim(title))
                       
        write(title, '(a)') 'BASIC V-WIND' 
        write(wfname, '(a,a2,a2,a2,a)') &
             '/export9/radiosonde/rusa/osan24_500/V/vstbas',arg1,arg2,arg3,'.nc'
        call out1d(trim(wfname),'VBAS',vstbas,'Z',nrst,zst,trim(title))
                    
        write(title, '(a)') 'BASIC RHO' 
        write(wfname, '(a,a2,a2,a2,a)') &
             '/export9/radiosonde/rusa/osan24_500/RHO/rstbas',arg1,arg2,arg3,'.nc'
        call out1d(trim(wfname),'RHOBAS',rstbas,'Z',nrst,zst,trim(title))
                          
        write(title, '(a)') 'BASIC THETA' 
        write(wfname, '(a,a2,a2,a2,a)') &
             '/export9/radiosonde/rusa/osan24_500/THE/thestbas',arg1,arg2,arg3,'.nc'
        call out1d(trim(wfname),'THEBAS',thestbas,'Z',nrst,zst,trim(title))
                                 
        write(title, '(a)') 'BASIC BVF' 
        write(wfname, '(a,a2,a2,a2,a)') &
             '/export9/radiosonde/rusa/osan24_500/BV/bvstbas',arg1,arg2,arg3,'.nc'
        call out1d(trim(wfname),'BVBAS',bvstbas,'Z',nrst,zst,trim(title))
                    
        write(title, '(a)') 'PERTURBATION TEMPERATURE' 
        write(wfname, '(a,a2,a2,a2,a)') &
             '/export9/radiosonde/rusa/osan24_500/T/tstprt',arg1,arg2,arg3,'.nc'
        call out1d(trim(wfname),'TPRT',tstprt,'Z',nrst,zst,trim(title))
                                                              
        write(title, '(a)') 'PERTURBATION PRESSURE' 
        write(wfname, '(a,a2,a2,a2,a)') &
             '/export9/radiosonde/rusa/osan24_500/P/pstprt',arg1,arg2,arg3,'.nc'
        call out1d(trim(wfname),'PPRT',pstprt,'Z',nrst,zst,trim(title))
                                      
        write(title, '(a)') 'PERTURBATION U-WIND' 
        write(wfname, '(a,a2,a2,a2,a)') &
             '/export9/radiosonde/rusa/osan24_500/U/ustprt',arg1,arg2,arg3,'.nc'
        call out1d(trim(wfname),'UPRT',ustprt,'Z',nrst,zst,trim(title))
                       
        write(title, '(a)') 'PERTURBATION V-WIND' 
        write(wfname, '(a,a2,a2,a2,a)') &
             '/export9/radiosonde/rusa/osan24_500/V/vstprt',arg1,arg2,arg3,'.nc'
        call out1d(trim(wfname),'VPRT',vstprt,'Z',nrst,zst,trim(title))
                    
!--------------------------------------------------------------------------------------------
! OUTPUT DATA OF TROPOSPHERE
!--------------------------------------------------------------------------------------------
                                 
        write(title, '(a)') 'TOTAL TEMPERATURE' 
        write(wfname, '(a,a2,a2,a2,a)') &
             '/export9/radiosonde/rusa/osan24_500/T/ttrtot',arg1,arg2,arg3,'.nc'
        call out1d(trim(wfname),'T',ttr,'Z',nrtr,ztr,trim(title))
                                                              
        write(title, '(a)') 'TOTAL PRESSURE' 
        write(wfname, '(a,a2,a2,a2,a)') &
             '/export9/radiosonde/rusa/osan24_500/P/ptrtot',arg1,arg2,arg3,'.nc'
        call out1d(trim(wfname),'P',ptr,'Z',nrtr,ztr,trim(title))
                                      
        write(title, '(a)') 'TOTAL U-WIND' 
        write(wfname, '(a,a2,a2,a2,a)') &
             '/export9/radiosonde/rusa/osan24_500/U/utrtot',arg1,arg2,arg3,'.nc'
        call out1d(trim(wfname),'U',utr,'Z',nrtr,ztr,trim(title))
                       
        write(title, '(a)') 'TOTAL V-WIND' 
        write(wfname, '(a,a2,a2,a2,a)') &
             '/export9/radiosonde/rusa/osan24_500/V/vtrtot',arg1,arg2,arg3,'.nc'
        call out1d(trim(wfname),'V',vtr,'Z',nrtr,ztr,trim(title))
                    
        write(title, '(a)') 'TOTAL RHO' 
        write(wfname, '(a,a2,a2,a2,a)') &
             '/export9/radiosonde/rusa/osan24_500/RHO/rtrtot',arg1,arg2,arg3,'.nc'
        call out1d(trim(wfname),'RHO',rtr,'Z',nrtr,ztr,trim(title))
                          
        write(title, '(a)') 'TOTAL THETA' 
        write(wfname, '(a,a2,a2,a2,a)') &
             '/export9/radiosonde/rusa/osan24_500/THE/thetrtot',arg1,arg2,arg3,'.nc'
        call out1d(trim(wfname),'THE',thetr,'Z',nrtr,ztr,trim(title))
                                 
        write(title, '(a)') 'TOTAL BVF' 
        write(wfname, '(a,a2,a2,a2,a)') &
             '/export9/radiosonde/rusa/osan24_500/BV/bvtrtot',arg1,arg2,arg3,'.nc'
        call out1d(trim(wfname),'BV',bvtr,'Z',nrtr,ztr,trim(title))
                              
                                
        write(title, '(a)') 'BASIC TEMPERATURE' 
        write(wfname, '(a,a2,a2,a2,a)') &
             '/export9/radiosonde/rusa/osan24_500/T/ttrbas',arg1,arg2,arg3,'.nc'
        call out1d(trim(wfname),'TBAS',ttrbas,'Z',nrtr,ztr,trim(title))
                                                              
        write(title, '(a)') 'BASIC PRESSURE' 
        write(wfname, '(a,a2,a2,a2,a)') &
             '/export9/radiosonde/rusa/osan24_500/P/ptrbas',arg1,arg2,arg3,'.nc'
        call out1d(trim(wfname),'PBAS',ptrbas,'Z',nrtr,ztr,trim(title))
                                      
        write(title, '(a)') 'BASIC U-WIND' 
        write(wfname, '(a,a2,a2,a2,a)') &
             '/export9/radiosonde/rusa/osan24_500/U/utrbas',arg1,arg2,arg3,'.nc'
        call out1d(trim(wfname),'UBAS',utrbas,'Z',nrtr,ztr,trim(title))
                       
        write(title, '(a)') 'BASIC V-WIND' 
        write(wfname, '(a,a2,a2,a2,a)') &
             '/export9/radiosonde/rusa/osan24_500/V/vtrbas',arg1,arg2,arg3,'.nc'
        call out1d(trim(wfname),'VBAS',vtrbas,'Z',nrtr,ztr,trim(title))
                    
        write(title, '(a)') 'BASIC RHO' 
        write(wfname, '(a,a2,a2,a2,a)') &
             '/export9/radiosonde/rusa/osan24_500/RHO/rtrbas',arg1,arg2,arg3,'.nc'
        call out1d(trim(wfname),'RHOBAS',rtrbas,'Z',nrtr,ztr,trim(title))
                          
        write(title, '(a)') 'BASIC THETA' 
        write(wfname, '(a,a2,a2,a2,a)') &
             '/export9/radiosonde/rusa/osan24_500/THE/thetrbas',arg1,arg2,arg3,'.nc'
        call out1d(trim(wfname),'THEBAS',thetrbas,'Z',nrtr,ztr,trim(title))
                                 
        write(title, '(a)') 'BASIC BVF' 
        write(wfname, '(a,a2,a2,a2,a)') &
             '/export9/radiosonde/rusa/osan24_500/BV/bvtrbas',arg1,arg2,arg3,'.nc'
        call out1d(trim(wfname),'BVBAS',bvtrbas,'Z',nrtr,ztr,trim(title))
                    
                             
        write(title, '(a)') 'PERTURBATION TEMPERATURE' 
        write(wfname, '(a,a2,a2,a2,a)') &
             '/export9/radiosonde/rusa/osan24_500/T/ttrprt',arg1,arg2,arg3,'.nc'
        call out1d(trim(wfname),'TPRT',ttrprt,'Z',nrtr,ztr,trim(title))
                                                              
        write(title, '(a)') 'PERTURBATION PRESSURE' 
        write(wfname, '(a,a2,a2,a2,a)') &
             '/export9/radiosonde/rusa/osan24_500/P/ptrprt',arg1,arg2,arg3,'.nc'
        call out1d(trim(wfname),'PPRT',ptrprt,'Z',nrtr,ztr,trim(title))
                                      
        write(title, '(a)') 'PERTURBATION U-WIND' 
        write(wfname, '(a,a2,a2,a2,a)') &
             '/export9/radiosonde/rusa/osan24_500/U/utrprt',arg1,arg2,arg3,'.nc'
        call out1d(trim(wfname),'UPRT',utrprt,'Z',nrtr,ztr,trim(title))
                       
        write(title, '(a)') 'PERTURBATION V-WIND' 
        write(wfname, '(a,a2,a2,a2,a)') &
             '/export9/radiosonde/rusa/osan24_500/V/vtrprt',arg1,arg2,arg3,'.nc'
        call out1d(trim(wfname),'VPRT',vtrprt,'Z',nrtr,ztr,trim(title))
                    
!--------------------------------------------------------------------------------------------
! CALCULATE EK, EP, ET
!-------------------------------------------------------------------------------------------- 
!       do k=1,nn
!         ek(k) = 0.5*(uprt(k)**2 + vprt(k)**2)
!         if (fbv(k) .le. 0.0 ) then
!           ep(k) = 0.0
!         else  
!           ep(k) = 0.5*(9.806*(tprt(k)/tbas(k))/fbv(k))**2
!         end if
!         if (bvbas(k) .le. 0.0 ) then
!           epbas(k) = 0.0
!         else
!           epbas(k) = 0.5*(9.806*(tprt(k)/tbas(k))/bvbas(k))**2
!         end if
!         et(k) = ek(k) + ep(k)
!         etbas(k) = ek(k) + epbas(k)
!       end do
!
!       print *, nn,zz,ww
!
!--------------------------------------------------------------------------------------------  
      else 
        stop
      end if

      STOP
      END

