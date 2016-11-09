!23456
   
      PROGRAM add_data

      Use netcdfio

      Implicit None

      integer,Parameter                     :: nday=68, nlev=151 
      Integer                               :: i,j 
      Integer                               :: ncid,it
      Integer,Dimension(nday)               :: date
      Real                                  :: MIS        
      Real, Dimension(nday)                 :: day 
      Real, Dimension(nday)                 :: tprtmax,uprtmax,vprtmax,bvprtmax
      Real, Dimension(nday)                 :: tprtmin,uprtmin,vprtmin,bvprtmin
      Real, Dimension(nlev)                 :: z,tprt,uprt,vprt,bvprt
      Real, Dimension(nday,nlev)            :: allt,allu,allv,allbv
      Character(len=255)                    :: fnamet,fnameu,fnamev,fnamebv
      Character(len=255)                    :: wfname,title
      Logical                               :: iex

      date = (/082000, 082006, 082012, 082018, 082100, 082106, 082112, 082118, &
               082200, 082206, 082212, 082218, 082300, 082306, 082312, 082318, &
               082400, 082406, 082412, 082418, 082500, 082506, 082512, 082518, &
               082600, 082606, 082612, 082618, 082700, 082706, 082712, 082718, &
               082800, 082806, 082812, 082818, 082900, 082906, 082912, 082918, &
               083000, 083006, 083012, 083018, 083100, 083106, 083112, 083118, &
               090100, 090106, 090112, 090118, 090200, 090206, 090212, 090218, &
               090300, 090306, 090312, 090318, 090400, 090406, 090412, 090418, &
               090500, 090506, 090512, 090518                                 /)
                           
      MIS = 1.e+32
              
      do i=1,nday

        write(fnamet, '(a,i6.6,a)')  &
                      '/export9/radiosonde/rusa/osan/T/tstprt', date(i), '.nc'
        write(fnameu, '(a,i6.6,a)')  &
                      '/export9/radiosonde/rusa/osan/U/ustprt', date(i), '.nc'
        write(fnamev, '(a,i6.6,a)')  &
                      '/export9/radiosonde/rusa/osan/V/vstprt', date(i), '.nc'
!        write(fnamebv, '(a,i6.6,a)')  &
                        '/export9/radiosonde/rusa/osan/BV/bvstprt', date(i), '.nc'
   
        inquire(file=fnamet,exist=iex) 
   
        if (iex) then  
          call opennc(fnamet,ncid) 
          call get1d(ncid, 'Z', nlev, z)
          call get1d(ncid, 'TPRT', nlev, tprt) 
          call closenc(ncid)     

          call opennc(fnameu,ncid) 
          call get1d(ncid, 'Z', nlev, z)
          call get1d(ncid, 'UPRT', nlev, uprt) 
          call closenc(ncid)     

          call opennc(fnamev,ncid) 
          call get1d(ncid, 'Z', nlev, z)
          call get1d(ncid, 'VPRT', nlev, vprt) 
          call closenc(ncid)     
     
!          call opennc(fnamebv,ncid) 
!          call get1d(ncid, 'Z', nlev, z)
!          call get1d(ncid, 'BVPRT', nlev, bvprt) 
!          call closenc(ncid)     
      
            tprtmax(i)  = maxval(abs(tprt(:))) 
            uprtmax(i)  = maxval(abs(uprt(:)))
            vprtmax(i)  = maxval(abs(vprt(:)))
!            bvprtmax(i) = maxval(abs(bvprt(:)))
            tprtmin(i)  = minval(abs(tprt(:)))
            uprtmin(i)  = minval(abs(uprt(:)))
            vprtmin(i)  = minval(abs(vprt(:)))
!            bvprtmin(i) = minval(abs(bvprt(:)))
    
          do j=1,nlev
            allt(i,j)   = tprt(j)
            allu(i,j)   = uprt(j)
            allv(i,j)   = vprt(j)
!            allbv(i,j)  = bvprt(j)
          end do
        else 
            tprtmax(i)  = MIS
            uprtmax(i)  = MIS
            vprtmax(i)  = MIS
!            bvprtmax(i) = MIS
            tprtmin(i)  = MIS
            uprtmin(i)  = MIS
            vprtmin(i)  = MIS
!            bvprtmin(i) = MIS
          do j=1,nlev
            allt(i,j)   = MIS
            allu(i,j)   = MIS
            allv(i,j)   = MIS 
!            allbv(i,j)  = MIS
          end do
        end if  
      end do 

      do i=1,nday
        day(i) = i
      end do

      write(title,  '(a)')  'TPRTMAX'
      write(wfname, '(a)')  '/export9/radiosonde/rusa/osan/add_data/tstprtmaxabs.nc'
      call out1d(trim(wfname), 'TPRTMAX', tprtmax, 'TIME', nday, day, trim(title))
    
      write(title,  '(a)')  'TPRTMIN'
      write(wfname, '(a)')  '/export9/radiosonde/rusa/osan/add_data/tstprtminabs.nc'
      call out1d(trim(wfname), 'TPRTMIN', tprtmin, 'TIME', nday, day, trim(title))

      write(title,  '(a)')  'UPRTMAX'
      write(wfname, '(a)')  '/export9/radiosonde/rusa/osan/add_data/ustprtmaxabs.nc'
      call out1d(trim(wfname), 'UPRTMAX', uprtmax, 'TIME', nday, day, trim(title))
    
      write(title,  '(a)')  'UPRTMIN'
      write(wfname, '(a)')  '/export9/radiosonde/rusa/osan/add_data/ustprtminabs.nc'
      call out1d(trim(wfname), 'UPRTMIN', uprtmin, 'TIME', nday, day, trim(title))
    
      write(title,  '(a)')  'VPRTMAX'
      write(wfname, '(a)')  '/export9/radiosonde/rusa/osan/add_data/vstprtmaxabs.nc'
      call out1d(trim(wfname), 'VPRTMAX', vprtmax, 'TIME', nday, day, trim(title))
    
      write(title,  '(a)')  'VPRTMIN'
      write(wfname, '(a)')  '/export9/radiosonde/rusa/osan/add_data/vstprtminabs.nc'
      call out1d(trim(wfname), 'VPRTMIN', vprtmin, 'TIME', nday, day, trim(title))

!      write(title,  '(a)')  'BVPRTMAX'
!      write(wfname, '(a)')  '/export9/radiosonde/rusa/osan/add_data/bvstprtmaxabs.nc'
!      call out1d(trim(wfname), 'BVPRTMAX', bvprtmax, 'TIME', nday, day, trim(title))
    
!      write(title,  '(a)')  'BVPRTMIN'
!      write(wfname, '(a)')  '/export9/radiosonde/rusa/osan/add_data/bvstprtminabs.nc'
!      call out1d(trim(wfname), 'BVPRTMIN', bvprtmin, 'TIME', nday, day, trim(title))


      write(title,  '(a)')  'all TPRT'
      write(wfname, '(a)')  '/export9/radiosonde/rusa/osan/add_data/tstprt.nc'
      call out2d(trim(wfname), 'TPRT', allt, 'TIME', nday, day,   &
                                'Z', nlev, Z, trim(title)              )
 
      write(title,  '(a)')  'all UPRT'
      write(wfname, '(a)')  '/export9/radiosonde/rusa/osan/add_data/ustprt.nc'
      call out2d(trim(wfname), 'UPRT', allu, 'TIME', nday, day,   &
                                'Z', nlev, Z, trim(title)              )

      write(title,  '(a)')  'all VPRT'
      write(wfname, '(a)')  '/export9/radiosonde/rusa/osan/add_data/vstprt.nc'
      call out2d(trim(wfname), 'VPRT', allv, 'TIME', nday, day,   &
                                'Z', nlev, Z, trim(title)              )

!      write(title,  '(a)')  'all BVPRT'
!      write(wfname, '(a)')  '/export9/radiosonde/rusa/osan/add_data/bvstprt.nc'
!      call out2d(trim(wfname), 'BVPRT', allbv, 'TIME', nday, day,   &
!                                'Z', nlev, Z, trim(title)              )

      end
     

