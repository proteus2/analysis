program uw_test              ! check units

  use netcdfio
  use fft,     only: fft2df, fft2db, wavnum
  use regress, only: lfit, funcs


 implicit none

 integer, parameter :: detbyuw = 1, pwanal = 0, k0int = 0 ! 2
 real,    parameter :: cmin = 25., cmax = 50.
 integer, parameter :: nx =1402, nz = 202, ntall = 480
 integer, parameter :: nxp = 600, nzp = 200, nt = 480
 integer, parameter :: xst = 700, zst = 0
 integer, parameter :: ma = 3
 real,    parameter :: dx = 1000., dt = 60., t0 = 0, dz = 300.
!-----------------------------------------------------
 integer, parameter :: nta = nt+nt/2
 integer :: i,j,n, ncid, ii
 integer :: ia(ma)
 real    :: x(nxp), t(nt), z(nzp), kk(nxp), oo(nta), cc, ccint, co
 real    :: rhobar(nz)
 real    :: a(ma), sig(nxp), afunc(50), trend(nxp), covar(ma+10,ma+10), chisq
 real, dimension(nzp)            :: ubar, ubar1, ubar2
 real, dimension(nx,nz)          :: udata, wdata, pdata
 real, dimension(0:ntall+2,nzp+2)  :: ubardat
 real, dimension(nxp,nzp,nt)  :: ua, wa, pa
 real, dimension(nxp,nta)     :: u, w, p, w_up, w_dn, u_up, u_dn, p_up, p_dn
 real, dimension(nt,nzp)      :: uw_up, uw_dn, pw_up, pw_dn
 double complex, dimension(nxp,nta) :: cu, cw, cp, cw_up, cw_dn, cu_up, cu_dn, cp_up, cp_dn
 double precision ::  rhobar8(nz)
 character*6   :: time
 character*128 :: fname, ftitle



 open(5,file='/usr/users/kyh/arps/work/2ndgen/x12jul/rhobar')
 read(5,*) rhobar8
 close(5)
 rhobar = real(rhobar8)

! x,t define
 do i=1, nxp
   x(i) = (i+xst-0.5) * dx - 100.e3   ! western sponge layer
 enddo
 do n=1, nt
   t(n) = n * dt + t0
 enddo
 do j=1, nzp
   z(j) = (j+zst-0.5) * dz
 enddo

! k, omega define
 call wavnum(nxp,dx,kk)
 call wavnum(nta,dt,oo)
 oo = (-1.) * oo

! read ubar or umean
 call opennc('/export23/kyh/srcjul/unbar.nc',ncid)
 call get2d(ncid, 'ubar', ntall+3, nzp+2, ubardat)
 call closenc(ncid)
 do j=1, nzp
   ubar1(j) = ubardat(0,zst+1+j)
   ubar2(j) = ubardat(ntall+2,zst+1+j)     ! avg. over selected time
 enddo

 ubar = ubar1
 print*,ubar,'    : u0'

! read u, w
print*, 'data reading'
 do n=1, nt
   write(time,'(i6.6)') int(t(n))
   open(1,file='/export23/kyh/srcjul/bin/uprt.bin'//time,    &
          form='unformatted',convert='big_endian')
   read(1) udata
   close(1)
   open(2,file='/export23/kyh/srcjul/bin/w.bin'//time,       &
          form='unformatted',convert='big_endian')
   read(2) wdata
   close(2)

   do j=1, nzp
   do i=1, nxp
     ua(i,j,n) = udata(xst+1+i,zst+1+j)
     wa(i,j,n) = wdata(xst+1+i,zst+1+j)
   enddo
   enddo
 enddo
 if (pwanal .ne. 0) then
   do n=1, nt
     write(time,'(i6.6)') int(t(n))
     open(3,file='/export23/kyh/srcjul/bin/pprt.bin'//time,    &
            form='unformatted',convert='big_endian')
     read(3) pdata
     close(3)

     do j=1, nzp
     do i=1, nxp
       pa(i,j,n) = pdata(xst+1+i,zst+1+j)
     enddo
     enddo
   enddo
 end if
print*, 'OK...'


 uw_up = 0.  ;  uw_dn = 0.
 pw_up = 0.  ;  pw_dn = 0.

 u = 0.   ;   w = 0.   ;   p = 0.
 
 do j=1, nzp

! 2-D FFT & determine up/downward
   do n=1, nt

     u(:,n) = ua(:,j,n)
     w(:,n) = wa(:,j,n)

     do i=1, ma
      ia(i) = 1           ! calculate the coeff.(=a) if ia.ne.0
       a(i) = 0.0         ! coeff. for ia=0
     end do
     do i=1, nxp
       sig(i) = 1.0
     end do
     call lfit(x,u(:,n),sig,nxp,a,ia,ma,covar,ma+10,chisq)
     do i=1, nxp
       call funcs(x(i),afunc,ma)
       trend(i)=0.
       do ii=1,ma
         trend(i) = trend(i) + a(ii)*afunc(ii)
       end do
       u(i,n) = u(i,n) - trend(i)
     end do

     do i=1, ma
      ia(i) = 1           ! calculate the coeff.(=a) if ia.ne.0
       a(i) = 0.0         ! coeff. for ia=0
     end do
     do i=1, nxp
       sig(i) = 1.0
     end do
     call lfit(x,w(:,n),sig,nxp,a,ia,ma,covar,ma+10,chisq)
     do i=1, nxp
       call funcs(x(i),afunc,ma)
       trend(i)=0.
       do ii=1,ma
         trend(i) = trend(i) + a(ii)*afunc(ii)
       end do
       w(i,n) = w(i,n) - trend(i)
     end do

   enddo

   if (pwanal .ne. 0) then
     do n=1, nt

       p(:,n) = pa(:,j,n)

       do i=1, ma
        ia(i) = 1           ! calculate the coeff.(=a) if ia.ne.0
         a(i) = 0.0         ! coeff. for ia=0
       end do
       do i=1, nxp
         sig(i) = 1.0
       end do
       call lfit(x,p(:,n),sig,nxp,a,ia,ma,covar,ma+10,chisq)
       do i=1, nxp
         call funcs(x(i),afunc,ma)
         trend(i)=0.
         do ii=1,ma
           trend(i) = trend(i) + a(ii)*afunc(ii)
         end do
         p(i,n) = p(i,n) - trend(i)
       end do

     enddo
   end if

   cw_up = (0.d0,0.d0)   ;   cw_dn = (0.d0,0.d0)
   cu_up = (0.d0,0.d0)   ;   cu_dn = (0.d0,0.d0)
   cp_up = (0.d0,0.d0)   ;   cp_dn = (0.d0,0.d0)

   call fft2df(nxp,nta,u,cu)
   call fft2df(nxp,nta,w,cw)
   call fft2df(nxp,nta,p,cp)

   if (detbyuw.eq.1 .and. pwanal.eq.0) then
     do n=1, nta
     do i=2+k0int, nxp
       co = real(cu(i,n))*real(cw(i,n))+aimag(cu(i,n))*aimag(cw(i,n))
       cc = oo(n)/kk(i)
       ccint = cc - ubar(j)
       if (cc.ge.cmin .and. cc.le.cmax) then
         if ((co*ccint) .gt. 0.0) then
           cw_up(i,n) = cw(i,n)
           cu_up(i,n) = cu(i,n)
         else if ((co*ccint) .lt. 0.0) then
           cw_dn(i,n) = cw(i,n)
           cu_dn(i,n) = cu(i,n)
         end if
       end if
     enddo
     enddo
   else if (detbyuw.eq.1 .and. pwanal.ne.0) then
     do n=1, nta
     do i=2+k0int, nxp
       co = real(cu(i,n))*real(cw(i,n))+aimag(cu(i,n))*aimag(cw(i,n))
       cc = oo(n)/kk(i)
       ccint = cc - ubar(j)
       if (cc.ge.cmin .and. cc.le.cmax) then
         if ((co*ccint) .gt. 0.0) then
           cw_up(i,n) = cw(i,n)
           cu_up(i,n) = cu(i,n)
           cp_up(i,n) = cp(i,n)
         else if ((co*ccint) .lt. 0.0) then
           cw_dn(i,n) = cw(i,n)
           cu_dn(i,n) = cu(i,n)
           cp_dn(i,n) = cp(i,n)
         end if
       end if
     enddo
     enddo
   else if (detbyuw.ne.1 .and. pwanal.ne.0) then
     do n=1, nta
     do i=2+k0int, nxp
       co = real(cp(i,n))*real(cw(i,n))+aimag(cp(i,n))*aimag(cw(i,n))
       cc = oo(n)/kk(i)
       if (cc.ge.cmin .and. cc.le.cmax) then
         if (co .gt. 0.0) then
           cw_up(i,n) = cw(i,n)
           cu_up(i,n) = cu(i,n)
           cp_up(i,n) = cp(i,n)
         else if (co .lt. 0.0) then
           cw_dn(i,n) = cw(i,n)
           cu_dn(i,n) = cu(i,n)
           cp_dn(i,n) = cp(i,n)
         end if
       end if
     enddo
     enddo
   else
     print*, ' ERROR. Check the option detbyuw and pwanal.'
   end if


! backward FFT
   call fft2db(nxp,nta,cw_up,w_up)
   call fft2db(nxp,nta,cw_dn,w_dn)
   call fft2db(nxp,nta,cu_up,u_up)
   call fft2db(nxp,nta,cu_dn,u_dn)
   if (pwanal .ne. 0) then
     call fft2db(nxp,nta,cp_up,p_up)
     call fft2db(nxp,nta,cp_dn,p_dn)
   end if

   do i=1, nxp
   do n=1, nt
     uw_up(n,j) = uw_up(n,j) + u_up(i,n)*w_up(i,n)*rhobar(zst+1+j)
     uw_dn(n,j) = uw_dn(n,j) + u_dn(i,n)*w_dn(i,n)*rhobar(zst+1+j)
   enddo
   enddo
   uw_up(:,j) = uw_up(:,j) / nxp
   uw_dn(:,j) = uw_dn(:,j) / nxp

   if (pwanal .ne. 0) then
     do i=1, nxp
     do n=1, nt
       pw_up(n,j) = pw_up(n,j) + p_up(i,n)*w_up(i,n)
       pw_dn(n,j) = pw_dn(n,j) + p_dn(i,n)*w_dn(i,n)
     enddo
     enddo
     pw_up(:,j) = pw_up(:,j) / nxp
     pw_dn(:,j) = pw_dn(:,j) / nxp
   end if

 enddo


! output
 if (detbyuw .eq. 1) then
   write(ftitle,'(a)') 'up/downward flux determined by uw, cp_hat'
 else
   write(ftitle,'(a)') 'up/downward flux determined by pw'
 end if
 write(fname,'(a)') '/usr/users/kyh/anal/work/mom_updn_tz.nc'

 if (pwanal .ne. 0) then
   call out2d(trim(fname),4,(/'mom_up','mom_dn','pw_up','pw_dn'/),  &
                            (/uw_up,uw_dn,pw_up,pw_dn/),            &
               't',nt,t/60.,'z',nzp,z/1000.,trim(ftitle))
 else
   call out2d(trim(fname),2,(/'mom_up','mom_dn'/),                  &
                            (/uw_up,uw_dn/),                        &
               't',nt,t/60.,'z',nzp,z/1000.,trim(ftitle))
 end if


 STOP

end

