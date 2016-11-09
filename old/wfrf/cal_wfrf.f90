 program cal_wfrf

 use params
 use cal_xsq

 implicit none

 include 'netcdf.inc'

 double precision                          :: pi
 double precision                          :: shear  ,ricld
 double precision   , dimension(nd)        :: wpd 
 double precision   , dimension(nc,nd)     :: xsq    ,wfrf   ,theta  ,frt
 double precision ::  thcoef, c0

 character(len=120)               :: fname
 character(len=120)               :: casename 
 integer                          :: i      ,j   
 integer                          :: indx   ,iargc
 integer                          :: istat  ,ncid 
 integer                          :: phadid ,phaid ,wpddid ,wpdid ,xsqid  ,wfrfid ,thetid ,frtid
 integer, dimension(2)            :: var2id  
 character(len=100)               :: arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10,arg11,arg12,arg13
 external getarg

 pi = 4.d0*atan(1.d0)

 indx = iargc()
 if (indx /= 13) then 
   print*, 'insert input parameters'
   stop 
 else
   call getarg(1 , arg1 ); read(arg1 ,*) zb
   call getarg(2 , arg2 ); read(arg2 ,*) zt
   call getarg(3 , arg3 ); read(arg3 ,*) ntr
   call getarg(4 , arg4 ); read(arg4 ,*) nst
   call getarg(5 , arg5 ); read(arg5 ,*) ubt
   call getarg(6 , arg6 ); read(arg6 ,*) vbt
   call getarg(7 , arg7 ); read(arg7 ,*) ubb
   call getarg(8 , arg8 ); read(arg8 ,*) vbb
   call getarg(9 , arg9 ); read(arg9 ,*) ub0
   call getarg(10, arg10); read(arg10,*) vb0
   call getarg(11, arg11); read(arg11,*) cqx
   call getarg(12, arg12); read(arg12,*) cqy
   call getarg(13, arg13); read(arg13,*) casename 
 end if

 call inti

 thcoef = (delh*delt/32./pi**1.5)**2.
 c0 = delh/delt

 do i = 1, nd 
   wpd(i) = 0. + dd*(i-1) 
 end do

 do i = 1, nd
   ut = ubt*cos(wpd(i)*pi/180.) + vbt*sin(wpd(i)*pi/180.)
   ub = ubb*cos(wpd(i)*pi/180.) + vbb*sin(wpd(i)*pi/180.)
   u0 = ub0*cos(wpd(i)*pi/180.) + vb0*sin(wpd(i)*pi/180.)
   cq = cqx*cos(wpd(i)*pi/180.) + cqy*sin(wpd(i)*pi/180.)

   shear = 0.0
   ricld = 1.e20

   if (abs(ub-u0) > 0.1) then  ! shear case
     zs = (ut-u0)/(ub-u0)*zb
     zs = min( zt, max( zs, zb ) )
     shear = (ut-u0) / zs
     ricld = (nst/shear)**2.
   else
     zs = zt
   end if

   if (ricld > 0.25) then
     if (ricld > 1.e+4) then
       print*,'uniform'
       call xsqun(ut ,xsq(:,i))
     else
       print*,'shear'
       call xsqsh(ut ,ub ,u0 ,xsq(:,i))
     end if
   end if

   do j = 1, nc
     if (xsq(j,i) .ne. 0.) then
       wfrf(j,i) = xsq(j,i) * nst/abs(cpgrd(j)-ut)
     else
       wfrf(j,i) = 0.
     end if
     
     ! critical level filtering
     if ( ut-1.5*dc <= cpgrd(j) .and. cpgrd(j) <= ut+1.5*dc ) wfrf(j,i) = 0.
   end do

   do j=1, nc
     theta(j,i) = thcoef / (1.+((cpgrd(j)-cq)/c0)**2.)
   enddo

 end do

 frt = wfrf * theta

 write(fname,'(A,A,A)') '/data3/kyh/analy/wfrf/',trim(casename),'.nc'
 istat = nf_create(trim(fname),nf_clobber,ncid)
 print*, fname

 istat  = nf_def_dim(ncid,'PHASE',nc ,phadid)
 istat  = nf_def_dim(ncid,'WPD'  ,nd ,wpddid)
 var2id = (/phadid,wpddid/)
 istat  = nf_def_var(ncid,'PHASE',nf_real,1,phadid,phaid )
 istat  = nf_def_var(ncid,'WPD'  ,nf_real,1,wpddid,wpdid )
 istat  = nf_def_var(ncid,'XSQ'  ,nf_real,2,var2id,xsqid )
 istat  = nf_def_var(ncid,'WFRF' ,nf_real,2,var2id,wfrfid)
 istat  = nf_def_var(ncid,'THETA',nf_real,2,var2id,thetid)
 istat  = nf_def_var(ncid,'FRT'  ,nf_real,2,var2id,frtid )

 istat  = nf_enddef(ncid)

 istat = nf_put_var_real(ncid,phaid ,real(cpgrd))
 istat = nf_put_var_real(ncid,wpdid ,real(wpd)  )
 istat = nf_put_var_real(ncid,xsqid ,real(xsq)  )
 istat = nf_put_var_real(ncid,wfrfid,real(wfrf) )
 istat = nf_put_var_real(ncid,thetid,real(theta))
 istat = nf_put_var_real(ncid,frtid ,real(frt)  )
 istat = nf_close(ncid)

 end
