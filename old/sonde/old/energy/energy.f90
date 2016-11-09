!23456

      PROGRAM wave_energy

!----------------------------------------------------------------------------------------
! THIS PROGRAM OBTAIN WAVE ENERGY
! WAVE ENERGY CONSIST OF PIOTENTIAL ENERGY, KINETIC ENERGY, TOTAL ENERGY
!----------------------------------------------------------------------------------------
      Use Netcdfio

      IMPLICIT NONE 

      Integer,Parameter                       :: nrst=15, nrtr=15, nday=68
      Integer                                 :: i,j
      Integer                                 :: ncid,errday1,errday2
      Integer,Dimension(nday)                 :: date
      Real                                    :: MIS,grav 
      Real,Dimension(nday)                    :: day
      Real,Dimension(nrst)                    :: zst
      Real,Dimension(nrtr)                    :: ztr
      Real,Dimension(nday,nrst)               :: tstbas,tstprt,ustprt,vstprt,nstbas
      Real,Dimension(nday,nrtr)               :: ttrbas,ttrprt,utrprt,vtrprt,ntrbas
      Real,Dimension(nday,nrst)               :: etst,epst,ekst
      Real,Dimension(nday,nrtr)               :: ettr,eptr,ektr
      Real,Dimension(nday)                    :: avgnstbas,avgntrbas
      Real,Dimension(nday)                    :: avgetst,avgepst,avgekst
      Real,Dimension(nday)                    :: avgettr,avgeptr,avgektr
      Character(len=255)                      :: wfname,fname,chfile,title
      Logical                                 :: iex
 
      grav = 9.806  ;  MIS = 1.e+32 
      tstbas(:,:)=0.0 ; tstprt(:,:)=0.0 ; ustprt(:,:)=0.0 
      vstprt(:,:)=0.0 ; nstbas(:,:)=0.0  
      ttrbas(:,:)=0.0 ; ttrprt(:,:)=0.0 ; utrprt(:,:)=0.0 
      vtrprt(:,:)=0.0 ; ntrbas(:,:)=0.0   
      etst(:,:)=0.0   ;  epst(:,:)=0.0  ;  ekst(:,:)=0.0    
      ettr(:,:)=0.0   ;  eptr(:,:)=0.0  ;  ektr(:,:)=0.0   
      avgnstbas(:)=0.0; avgntrbas(:)=0.0
      avgetst(:)=0.0  ; avgepst(:)=0.0  ; avgekst(:)=0.0 
      avgettr(:)=0.0  ; avgeptr(:)=0.0  ; avgektr(:)=0.0 

      date = (/082000, 082006, 082012, 082018, 082100, 082106, 082112, 082118, &
               082200, 082206, 082212, 082218, 082300, 082306, 082312, 082318, &
               082400, 082406, 082412, 082418, 082500, 082506, 082512, 082518, &
               082600, 082606, 082612, 082618, 082700, 082706, 082712, 082718, &
               082800, 082806, 082812, 082818, 082900, 082906, 082912, 082918, &
               083000, 083006, 083012, 083018, 083100, 083106, 083112, 083118, &
               090100, 090106, 090112, 090118, 090200, 090206, 090212, 090218, &
               090300, 090306, 090312, 090318, 090400, 090406, 090412, 090418, &
               090500, 090506, 090512, 090518                                 /)

       errday1 = 082400  ! OSAN-ERROR
       errday2 = 082700  ! OSAN-ERROR

!----------------------------------------------------------------------------------------
! OPEN INPUT FILE IN STRATOSPHERE
!----------------------------------------------------------------------------------------

          write(fname, '(a)' )  &
                             '/export9/radiosonde/rusa/osan24_500/add_data/st/tstprt0208.nc'
          call opennc(fname,ncid)  
          call get1d(ncid,'Z',nrst,zst) 
          call get1d(ncid,'TIME',nday,day) 
          call get2d(ncid,'TPRT',nday,nrst,tstprt) 
          call closenc(ncid)     
                      
          write(fname, '(a)' )  &
                             '/export9/radiosonde/rusa/osan24_500/add_data/st/tstbas0208.nc'
          call opennc(fname,ncid)  
          call get2d(ncid,'TBAS',nday,nrst,tstbas) 
          call closenc(ncid)     
        
          write(fname, '(a)' )  &
                             '/export9/radiosonde/rusa/osan24_500/add_data/st/ustprt0208.nc'
          call opennc(fname,ncid)  
          call get2d(ncid,'UPRT',nday,nrst,ustprt) 
          call closenc(ncid)     
                  
          write(fname, '(a)' )  &
                             '/export9/radiosonde/rusa/osan24_500/add_data/st/vstprt0208.nc'
          call opennc(fname,ncid)  
          call get2d(ncid,'VPRT',nday,nrst,vstprt) 
          call closenc(ncid)     

          write(fname, '(a)' )  &
                             '/export9/radiosonde/rusa/osan24_500/add_data/st/bvstbas0208.nc'
          call opennc(fname,ncid)  
          call get2d(ncid,'BVBAS',nday,nrst,nstbas)
          call closenc(ncid)     

!----------------------------------------------------------------------------------------
! CALCULATE POTENTIAL ENERGY, KINETIC ENERGY, TOTAL ENERGY IN ST
!----------------------------------------------------------------------------------------
                 
      do j=1,nday
        write(chfile,'(a,i6.6,a)' ) & 
                           '/export9/radiosonde/rusa/osan24_500/T/tstprt',date(j),'.nc'
        inquire(file=chfile,exist=iex)
        if (iex) then
          do i=1,nrst
            avgnstbas(j) = avgnstbas(j) + nstbas(j,i)
          end do
          avgnstbas(j)=avgnstbas(j)/float(nrst)      
       
          do i=1,nrst
            ekst(j,i) = 0.5*(ustprt(j,i))**2 + 0.5*(vstprt(j,i))**2
            epst(j,i) = 0.5*(grav*(tstprt(j,i)/tstbas(j,i))/avgnstbas(j))**2
            etst(j,i) = ekst(j,i) + epst(j,i)
          end do
       
          do i=1,nrst
            avgekst(j) = avgekst(j) + ekst(j,i)      
            avgepst(j) = avgepst(j) + epst(j,i)      
            avgetst(j) = avgetst(j) + etst(j,i)
          end do
          avgekst(j) = avgekst(j)/float(nrst) 
          avgepst(j) = avgepst(j)/float(nrst)
          avgetst(j) = avgetst(j)/float(nrst) 
        else
          do i=1,nrst
            ekst(j,i) = MIS
            epst(j,i) = MIS
            etst(j,i) = MIS
          end do         
          avgekst(j) = MIS
          avgepst(j) = MIS
          avgetst(j) = MIS
        end if  
        if (date(j).eq.errday1 .or. date(j).eq.errday2) then
          do i=1,nrst
            ekst(j,i) = MIS
            epst(j,i) = MIS
            etst(j,i) = MIS
          end do         
          avgekst(j) = MIS
          avgepst(j) = MIS
          avgetst(j) = MIS
        end if 
      end do
 
          write(title, '(a)' ) 'AVGETST'
          write(wfname,'(a)' ) &
                       '/export9/radiosonde/rusa/osan24_500/add_data/st/avgetst0208.nc'
          call out1d(trim(wfname),'AVGETST',avgetst,'TIME',nday,day,trim(title))
       
          write(title, '(a)' ) 'AVGEPST'
          write(wfname,'(a)' ) &
                       '/export9/radiosonde/rusa/osan24_500/add_data/st/avgepst0208.nc'
          call out1d(trim(wfname),'AVGEPST',avgepst,'TIME',nday,day,trim(title))
    
          write(title, '(a)' ) 'AVGEKST'
          write(wfname,'(a)' ) &
                       '/export9/radiosonde/rusa/osan24_500/add_data/st/avgekst0208.nc'
          call out1d(trim(wfname),'AVGEKST',avgekst,'TIME',nday,day,trim(title))

          write(title, '(a)' ) 'ETST'
          write(wfname,'(a)' ) &
                       '/export9/radiosonde/rusa/osan24_500/add_data/st/etst0208.nc'
          call out2d(trim(wfname),'ETST',etst,'TIME',nday,day,'Z',nrst,zst,trim(title))
                                    
          write(title, '(a)' ) 'EKST'
          write(wfname,'(a)' ) &
                       '/export9/radiosonde/rusa/osan24_500/add_data/st/ekst0208.nc'
          call out2d(trim(wfname),'EKST',ekst,'TIME',nday,day,'Z',nrst,zst,trim(title))
                                   
          write(title, '(a)' ) 'EPST'
          write(wfname,'(a)' ) &
                       '/export9/radiosonde/rusa/osan24_500/add_data/st/epst0208.nc'
          call out2d(trim(wfname),'EPST',epst,'TIME',nday,day,'Z',nrst,zst,trim(title))
                              
   
!----------------------------------------------------------------------------------------
! CALCULATE POTENTIAL ENERGY, KINETIC ENERGY, TOTAL ENERGY IN TR
!----------------------------------------------------------------------------------------

          write(fname, '(a)' )  &
                             '/export9/radiosonde/rusa/osan24_500/add_data/tr/ttrprt0208.nc'
          call opennc(fname,ncid)  
          call get1d(ncid,'Z',nrtr,ztr) 
          call get2d(ncid,'TPRT',nday,nrtr,ttrprt) 
          call closenc(ncid)     
                  
          write(fname, '(a)' )  &
                             '/export9/radiosonde/rusa/osan24_500/add_data/tr/ttrbas0208.nc'
          call opennc(fname,ncid)  
          call get2d(ncid,'TBAS',nday,nrtr,ttrbas) 
          call closenc(ncid)     
                
          write(fname, '(a)' )  &
                             '/export9/radiosonde/rusa/osan24_500/add_data/tr/utrprt0208.nc'
          call opennc(fname,ncid)  
          call get2d(ncid,'UPRT',nday,nrtr,utrprt) 
          call closenc(ncid)     
                    
          write(fname, '(a)' )  &
                             '/export9/radiosonde/rusa/osan24_500/add_data/tr/vtrprt0208.nc'
          call opennc(fname,ncid)  
          call get2d(ncid,'VPRT',nday,nrtr,vtrprt) 
          call closenc(ncid)     
                          
          write(fname, '(a)' )  &
                             '/export9/radiosonde/rusa/osan24_500/add_data/tr/bvtrbas0208.nc'
          call opennc(fname,ncid)  
          call get2d(ncid,'BVBAS',nday,nrtr,ntrbas)
          call closenc(ncid)     
     
      do j=1,nday
        write(chfile,'(a,i6.6,a)' ) & 
                         '/export9/radiosonde/rusa/osan24_500/T/tstprt',date(j),'.nc'
        inquire(file=chfile,exist=iex)
        if (iex) then
          do i=1,nrtr
            avgntrbas(j) = avgntrbas(j) + ntrbas(j,i)
          end do
  
          avgntrbas(j)=avgntrbas(j)/float(nrtr)      
   
          do i=1,nrtr
            ektr(j,i) = 0.5*(utrprt(j,i))**2 + 0.5*(vtrprt(j,i))**2
            eptr(j,i) = 0.5*(grav*(ttrprt(j,i)/ttrbas(j,i))/avgntrbas(j))**2
            ettr(j,i) = ektr(j,i) + eptr(j,i)
          end do
          do i=1,nrtr
            avgektr(j) = avgektr(j) + ektr(j,i)      
            avgeptr(j) = avgeptr(j) + eptr(j,i)      
            avgettr(j) = avgettr(j) + ettr(j,i)
          end do
          avgektr(j) = avgektr(j)/float(nrtr) 
          avgeptr(j) = avgeptr(j)/float(nrtr)
          avgettr(j) = avgettr(j)/float(nrtr) 
        else 
          do i=1,nrtr
            ektr(j,i) = MIS
            eptr(j,i) = MIS
            ettr(j,i) = MIS
          end do         
          avgektr(j) = MIS
          avgeptr(j) = MIS
          avgettr(j) = MIS
        end if
        if (date(j).eq.errday1 .or. date(j).eq.errday2) then
          do i=1,nrtr
            ektr(j,i) = MIS
            eptr(j,i) = MIS
            ettr(j,i) = MIS
          end do         
          avgektr(j) = MIS
          avgeptr(j) = MIS
          avgettr(j) = MIS
        end if  
      end do     
      
!        do i=1,nrtr
!          ektr(37,i) = MIS
!          eptr(37,i) = MIS
!          ettr(37,i) = MIS
!        end do
!          avgektr(37) = MIS
!          avgeptr(37) = MIS
!          avgettr(37) = MIS
!

          write(title, '(a)' ) 'AVGETTR'
          write(wfname,'(a)' ) &
                   '/export9/radiosonde/rusa/osan24_500/add_data/tr/avgettr0208.nc'
          call out1d(trim(wfname),'AVGETTR',avgettr,'TIME',nday,day,trim(title))
       
          write(title, '(a)' ) 'AVGEPTR'
          write(wfname,'(a)' ) &
                       '/export9/radiosonde/rusa/osan24_500/add_data/tr/avgeptr0208.nc'
          call out1d(trim(wfname),'AVGEPTR',avgeptr,'TIME',nday,day,trim(title))
    
          write(title, '(a)' ) 'AVGEKTR'
          write(wfname,'(a)' ) &
                       '/export9/radiosonde/rusa/osan24_500/add_data/tr/avgektr0208.nc'
          call out1d(trim(wfname),'AVGEKTR',avgektr,'TIME',nday,day,trim(title))
    
          write(title, '(a)' ) 'ETTR'
          write(wfname,'(a)' ) &
                       '/export9/radiosonde/rusa/osan24_500/add_data/tr/ettr0208.nc'
          call out2d(trim(wfname),'ETTR',ettr,'TIME',nday,day,'Z',nrtr,ztr,trim(title))
                                
          write(title, '(a)' ) 'EKTR'
          write(wfname,'(a)' ) &
                       '/export9/radiosonde/rusa/osan24_500/add_data/tr/ektr0208.nc'
          call out2d(trim(wfname),'EKTR',ektr,'TIME',nday,day,'Z',nrtr,ztr,trim(title))
                               
          write(title, '(a)' ) 'EPTR'
          write(wfname,'(a)' ) &
                       '/export9/radiosonde/rusa/osan24_500/add_data/tr/eptr0208.nc'
          call out2d(trim(wfname),'EPTR',eptr,'TIME',nday,day,'Z',nrtr,ztr,trim(title))

     end


