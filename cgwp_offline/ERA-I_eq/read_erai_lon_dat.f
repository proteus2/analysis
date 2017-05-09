      subroutine read_erai_lon_data (nf,fname,
     &                               lon,zkm,data,nlon,nlev) 
      character*80 fname,label1,label2,label3
      real lon(nlon),zkm(nlev),data(nlon,nlev)  
      if (nlon.ne.360) then 
         write(6,*) 'nlon must be 360' 
         stop
      endif 
      if (nlev.ne.60) then 
         write(6,*) 'nlev must be 60' 
         stop
      endif 
      open (nf,file=fname,status='old')
      read (nf,5000) label1
      read (nf,5010) lon
      read (nf,5000) label2
      read (nf,5020) zkm
      read (nf,5000) label3
      do l=1,nlev 
         read (nf,5020) (data(i,l),i=1,nlon)
      enddo
      close (nf)         
c      write(6,*) label1
c      write(6,*) lon
c      write(6,*) label2
c      write(6,*) zkm
c      write(6,*) label3
c      do l=1,nlev 
c         write(6,*) (data(i,l),i=1,nlon)
c      enddo
 5000 format(a80)
 5010 format(360f8.2)
 5020 format(360f13.8)
      return
      end
