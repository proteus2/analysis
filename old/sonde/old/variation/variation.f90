!23456
      
      PROGRAM variation

      Use netcdfio

      Implicit None

      Integer                              ::   j
      Integer                              ::   nt,nlev
      Integer                              ::   ncid
      Character(len=255)                   ::   fnamet,fnameu,fnamev 
      Real                                 ::   stn,avg,vac,mis
      Integer, Dimension(68)               ::   date
      Real, Allocatable, Dimension(:)      ::   tprt,vprt,uprt,mid
 
      date = (/082000, 082006, 082012, 082018, 082100, 082106, 082112, 082118, &
               082200, 082206, 082212, 082218, 082300, 082306, 082312, 082318, &
               082400, 082406, 082412, 082418, 082500, 082506, 082512, 082518, &
               082600, 082606, 082612, 082618, 082700, 082706, 082712, 082718, &
               082800, 082806, 082812, 082818, 082900, 082906, 082912, 082918, &
               083000, 083006, 083012, 083018, 083100, 083106, 083112, 083118, &
               090100, 090106, 090112, 090118, 090200, 090206, 090212, 090218, &
               090300, 090306, 090312, 090318, 090400, 090406, 090412, 090418, &
               090500, 090506, 090512, 090518                                  /)

      nlev=0 ; stn=0.0 ; avg=0.0 ; vac=0.0
      
      mis = 1.e+32

      OPEN(10, FILE='osan24_500.txt', FORM='FORMATTED')

      write(fnamet,'(a)') '/export9/radiosonde/rusa/osan24_500/add_data/st/tstprtmaxabs.nc'    
      call opennc(fnamet,ncid)
      call getdimlen(ncid,'TIME',nt) 
 
      allocate(tprt(1:nt)) ; tprt(:)=0.0 

      call get1d(ncid,'TPRTMAX',nt,tprt) 
      call closenc(ncid)


      write(fnameu,'(a)') '/export9/radiosonde/rusa/osan24_500/add_data/st/ustprtmaxabs.nc'    
      call opennc(fnameu,ncid)
   
      allocate(uprt(1:nt)) ; uprt(:)=0.0 

      call get1d(ncid,'UPRTMAX',nt,uprt) 
      call closenc(ncid)


      write(fnamev,'(a)') '/export9/radiosonde/rusa/osan24_500/add_data/st/vstprtmaxabs.nc'    
      call opennc(fnamev,ncid)
   
      allocate(vprt(1:nt)) ; vprt(:)=0.0 

      call get1d(ncid,'VPRTMAX',nt,vprt) 
      call closenc(ncid)

      do j=1,nt
        if (tprt(j) .ne. mis) then
          nlev = nlev+1
        end if
      end do
 
      allocate(mid(1:nlev)) ; mid(:)=0.0
 
      nlev=0 

      do j=1,nt
        if (tprt(j) .ne. mis) then
          nlev = nlev+1
          mid(nlev) = tprt(j)
          write(10,*) 'EXIST T-DATE'
          write(10,*) j, date(j) 
        end if
      end do
   
      do j=1,nlev
        avg = avg + mid(j)/float(nlev)
        vac = vac + mid(j)**2/float(nlev) 
      end do
      stn = (vac - avg**2)
      stn = sqrt(stn)
 
      do j=1,nlev
        if (mid(j) .gt. avg+3*stn) then
          write(10,*) 'T-ERROR'
          write(10,*) j,avg,stn,mid(j)
        end if
      end do 

      nlev=0 ; stn=0.0 ; avg=0.0 ; vac=0.0
      mid(:)=0.0 

      do j=1,nt
        if (uprt(j) .ne. mis) then
          nlev = nlev+1
          mid(nlev) = uprt(j)
        end if
      end do
   
      do j=1,nlev
        avg = avg + mid(j)/float(nlev)
        vac = vac + mid(j)**2/float(nlev) 
      end do
      stn = (vac - avg**2)
      stn = sqrt(stn)
 
      do j=1,nlev
        if (mid(j) .gt. avg+3*stn) then
          write(10,*) 'U-ERROR'
          write(10,*) j,avg,stn,mid(j)
        end if
      end do 

      nlev=0 ; stn=0.0 ; avg=0.0 ; vac=0.0
      mid(:)=0.0 

      do j=1,nt
        if (vprt(j) .ne. mis) then
          nlev = nlev+1
          mid(nlev) = vprt(j)
        end if
      end do
   
      do j=1,nlev
        avg = avg + mid(j)/float(nlev)
        vac = vac + mid(j)**2/float(nlev) 
      end do
      stn = (vac - avg**2)
      stn = sqrt(stn)
 
      do j=1,nlev
        if (mid(j) .gt. avg+3*stn) then
          write(10,*) 'V-ERROR'
          write(10,*) j,avg,stn,mid(j)
        end if
      end do 

      stop
      end

