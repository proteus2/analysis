!23456
   
      PROGRAM add_data

      Use netcdfio

      Implicit None

      integer,Parameter                 :: nday=68, nlev=15
      Integer                           :: i,j 
      Integer                           :: ncid,it
      Integer,Dimension(nday)           :: date
      Real                              :: MIS        
      Real, Dimension(nday)             :: day 
      Real, Dimension(nday)             :: tbasmax,ubasmax,vbasmax,bvbasmax
      Real, Dimension(nday)             :: tbasmin,ubasmin,vbasmin,bvbasmin
      Real, Dimension(nlev)             :: z,tbas,ubas,vbas,bvbas,rbas
      Real, Dimension(nday,nlev)        :: allt,allu,allv,allbv,allr
      Character(len=255)                :: fnamet,fnameu,fnamev,fnamebv,fnamer
      Character(len=255)                :: wfname,title
      Logical                           :: iex

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
                      '/export9/radiosonde/rusa/osan24_500/T/tstbas', date(i), '.nc'
        write(fnameu, '(a,i6.6,a)')  &
                      '/export9/radiosonde/rusa/osan24_500/U/ustbas', date(i), '.nc'
        write(fnamev, '(a,i6.6,a)')  &
                      '/export9/radiosonde/rusa/osan24_500/V/vstbas', date(i), '.nc'
        write(fnamebv, '(a,i6.6,a)')  &
                        '/export9/radiosonde/rusa/osan24_500/BV/bvstbas', date(i), '.nc'
        write(fnamer, '(a,i6.6,a)')  &
                        '/export9/radiosonde/rusa/osan24_500/RHO/rstbas', date(i), '.nc'
   
        inquire(file=fnamet,exist=iex) 
   
        if (iex) then  
          call opennc(fnamet,ncid) 
          call get1d(ncid, 'Z', nlev, z)
          call get1d(ncid, 'TBAS', nlev, tbas) 
          call closenc(ncid)     

          call opennc(fnameu,ncid) 
          call get1d(ncid, 'Z', nlev, z)
          call get1d(ncid, 'UBAS', nlev, ubas) 
          call closenc(ncid)     

          call opennc(fnamev,ncid) 
          call get1d(ncid, 'Z', nlev, z)
          call get1d(ncid, 'VBAS', nlev, vbas) 
          call closenc(ncid)     
     
          call opennc(fnamebv,ncid) 
          call get1d(ncid, 'Z', nlev, z)
          call get1d(ncid, 'BVBAS', nlev, bvbas) 
          call closenc(ncid)     
       
          call opennc(fnamer,ncid) 
          call get1d(ncid, 'Z', nlev, z)
          call get1d(ncid, 'RHOBAS', nlev, rbas) 
          call closenc(ncid)     

            tbasmax(i)  = maxval(abs(tbas(:))) 
            ubasmax(i)  = maxval(abs(ubas(:)))
            vbasmax(i)  = maxval(abs(vbas(:)))
            bvbasmax(i) = maxval(abs(bvbas(:)))
            tbasmin(i)  = minval(abs(tbas(:)))
            ubasmin(i)  = minval(abs(ubas(:)))
            vbasmin(i)  = minval(abs(vbas(:)))
            bvbasmin(i) = minval(abs(bvbas(:)))
    
          do j=1,nlev
            allt(i,j)   = tbas(j)
            allu(i,j)   = ubas(j)
            allv(i,j)   = vbas(j)
            allbv(i,j)  = bvbas(j)
            allr(i,j)   = rbas(j)
          end do
        else 
            tbasmax(i)  = MIS
            ubasmax(i)  = MIS
            vbasmax(i)  = MIS
            bvbasmax(i) = MIS
            tbasmin(i)  = MIS
            ubasmin(i)  = MIS
            vbasmin(i)  = MIS
            bvbasmin(i) = MIS
          do j=1,nlev
            allt(i,j)   = MIS
            allu(i,j)   = MIS
            allv(i,j)   = MIS 
            allbv(i,j)  = MIS
            allr(i,j)   = MIS
          end do
        end if  
      end do 

      do i=1,nday
        day(i) = i
      end do

      write(title,  '(a)')  'TBASMAX'
      write(wfname, '(a)')  '/export9/radiosonde/rusa/osan24_500/add_data/st/tstbasmaxabs.nc'
      call out1d(trim(wfname), 'TBASMAX', tbasmax, 'TIME', nday, day, trim(title))
    
      write(title,  '(a)')  'TBASMIN'
      write(wfname, '(a)')  '/export9/radiosonde/rusa/osan24_500/add_data/st/tstbasminabs.nc'
      call out1d(trim(wfname), 'TBASMIN', tbasmin, 'TIME', nday, day, trim(title))

      write(title,  '(a)')  'UBASMAX'
      write(wfname, '(a)')  '/export9/radiosonde/rusa/osan24_500/add_data/st/ustbasmaxabs.nc'
      call out1d(trim(wfname), 'UBASMAX', ubasmax, 'TIME', nday, day, trim(title))
    
      write(title,  '(a)')  'UBASMIN'
      write(wfname, '(a)')  '/export9/radiosonde/rusa/osan24_500/add_data/st/ustbasminabs.nc'
      call out1d(trim(wfname), 'UBASMIN', ubasmin, 'TIME', nday, day, trim(title))
    
      write(title,  '(a)')  'VBASMAX'
      write(wfname, '(a)')  '/export9/radiosonde/rusa/osan24_500/add_data/st/vstbasmaxabs.nc'
      call out1d(trim(wfname), 'VBASMAX', vbasmax, 'TIME', nday, day, trim(title))
    
      write(title,  '(a)')  'VBASMIN'
      write(wfname, '(a)')  '/export9/radiosonde/rusa/osan24_500/add_data/st/vstbasminabs.nc'
      call out1d(trim(wfname), 'VBASMIN', vbasmin, 'TIME', nday, day, trim(title))

      write(title,  '(a)')  'BVBASMAX'
      write(wfname, '(a)')  '/export9/radiosonde/rusa/osan24_500/add_data/st/bvstbasmaxabs.nc'
      call out1d(trim(wfname), 'BVBASMAX', bvbasmax, 'TIME', nday, day, trim(title))
    
      write(title,  '(a)')  'BVBASMIN'
      write(wfname, '(a)')  '/export9/radiosonde/rusa/osan24_500/add_data/st/bvstbasminabs.nc'
      call out1d(trim(wfname), 'BVBASMIN', bvbasmin, 'TIME', nday, day, trim(title))


      write(title,  '(a)')  'all TBAS'
      write(wfname, '(a)')  '/export9/radiosonde/rusa/osan24_500/add_data/st/tstbas0208.nc'
      call out2d(trim(wfname), 'TBAS', allt, 'TIME', nday, day,   &
                                'Z', nlev, Z, trim(title)              )
 
      write(title,  '(a)')  'all UBAS'
      write(wfname, '(a)')  '/export9/radiosonde/rusa/osan24_500/add_data/st/ustbas0208.nc'
      call out2d(trim(wfname), 'UBAS', allu, 'TIME', nday, day,   &
                                'Z', nlev, Z, trim(title)              )

      write(title,  '(a)')  'all VBAS'
      write(wfname, '(a)')  '/export9/radiosonde/rusa/osan24_500/add_data/st/vstbas0208.nc'
      call out2d(trim(wfname), 'VBAS', allv, 'TIME', nday, day,   &
                                'Z', nlev, Z, trim(title)              )

      write(title,  '(a)')  'all BVBAS'
      write(wfname, '(a)')  '/export9/radiosonde/rusa/osan24_500/add_data/st/bvstbas0208.nc'
      call out2d(trim(wfname), 'BVBAS', allbv, 'TIME', nday, day,   &
                                'Z', nlev, Z, trim(title)              )

      write(title,  '(a)')  'all RHOBAS'
      write(wfname, '(a)')  '/export9/radiosonde/rusa/osan24_500/add_data/st/rstbas0208.nc'
      call out2d(trim(wfname), 'RBAS', allr, 'TIME', nday, day,   &
                               'Z', nlev, Z, trim(title)              )
      end
     

