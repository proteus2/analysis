program uw_test              ! check units

  use netcdfio
  use fft,     only: fft2df, fft2db, wavnum, fft2dcb

!-----------------------------------------------------------------------
!
! OUTPUT : upward w'                    (wup)
!          downward w'                  (wdn)
!          k_ome cospectrum             (kozco)
!          upward   k_ome cospec.       (kozcoup)
!          downward k_ome cospec.       (kozcodn)
!          k_ome PSD                    (koz)
!          upward   k_ome PSD           (kozup)
!          downward k_ome PSD           (kozdn)
!
! INPUT  : ubar, uprt, w binary files
!
! norm   : perform the normalization
!
! fieldopt : separation of w' field
! coopt    : separation of cospectum
! psopt    : separation of w' PSD 
!            (do not use psopt and fieldopt
!             when partial reflection portion is significant.)
!
!-----------------------------------------------------

 implicit none

 integer, parameter :: fieldopt = 0, coopt = 1, psopt = 0, norm = 1
 integer, parameter :: nx =1402, nz = 352, ntall = 480
 integer, parameter :: nxp =1200, nzp = 200, nt = 360
 integer, parameter :: xst = 100, zst = 0
 real,    parameter :: dx = 1000., dt = 60., t0 = 7200, dz = 300.
!-----------------------------------------------------
 integer, parameter :: nkknew = nxp/2+1, noonew = nt/2*2+1
 integer :: i,j,n, ncid, ii, nn
 real    :: x(nxp), t(nt), z(nzp), kk(nxp), oo(nt)
 real    :: ccint, mabs, co, pi
 real    :: kknew(nkknew), oonew(noonew), rhort(nz), m2(nkknew,noonew,nzp)
 real*8  :: rhobar(nz)
 real, dimension(nzp)            :: ubar, ubar1, ubar2, hinv
 real, dimension(nx,nz)          :: udata, wdata
 real, dimension(0:ntall+2,nzp+2)  :: ubardat
 real, dimension(nxp,nzp,nt)  :: ua, wa
 real, dimension(nxp,nt)      :: u, w, w_up, w_dn, co_up, co_resi
 real, dimension(:,:,:), allocatable :: wa_up, wa_dn
 real, dimension(:,:,:), allocatable :: pwsd1, pwsd2, cosd1, cosd2

 complex :: iii
 double complex :: cw_pos, cw_neg, cu_pos, cu_neg, dd
 double complex, dimension(nxp,nt) :: cu, cw, cw_up, cw_dn

 character*6   :: time
 character*128 :: fname

 data pi/3.141593/, iii/(0.,1.)/


 open(4,file='/usr/users/kyh/anal/work/rhoalljul')
 read(4,*) rhort
 close(4)
 open(5,file='/usr/users/kyh/arps/work/2ndgen/x12jul/rhobar')
 read(5,*) rhobar
 close(5)
 if (norm .eq. 0) then
   rhort = 1.
 end if
 do j=1, nzp
   hinv(j) = real(-(dlog(rhobar(j+2+zst))-dlog(rhobar(j+zst))) / (2.*dz))
 enddo
print*,1./hinv
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
 call wavnum(nt,dt,oo)
 kk = 2.*pi * kk
 oo = 2.*pi * (-oo)

! read ubar or umean
 call opennc('/usr/users/kyh/result4/srcjul/UNbar/unbar.nc',ncid)
 call get2d(ncid, 'ubar', ntall+3, nzp+2, ubardat)
 call closenc(ncid)
 do j=1, nzp
   ubar1(j) = ubardat(0,zst+1+j)
   ubar2(j) = ubardat(ntall+2,zst+1+j)     ! avg. over selected time
 enddo

 ubar = ubar1
 print*,ubar,'    : u0'

! read m^2
 call opennc('/export22/kyh/ing/x12jul/spec/msq/m2.nc',ncid)
 call geta3d(ncid,'m2',1,nkknew,1,noonew,zst+2,nzp,m2)
 call closenc(ncid)

! read u, w
print*, 'data reading'
 do n=1, nt
   write(time,'(i6.6)') int(t(n))
   open(1,file='/export24/kyh/ing/srcjul/bin/uprt.bin'//time,    &
          form='unformatted',convert='big_endian')
   read(1) udata
   close(1)
   open(2,file='/export24/kyh/ing/srcjul/bin/w.bin'//time,       &
          form='unformatted',convert='big_endian')
   read(2) wdata
   close(2)

   do j=1, nzp
   do i=1, nxp
     ua(i,j,n) = udata(xst+1+i,zst+1+j) * rhort(zst+1+j)
     wa(i,j,n) = wdata(xst+1+i,zst+1+j) * rhort(zst+1+j)
   enddo
   enddo
 enddo
print*, 'OK...'

 if (fieldopt .eq. 1) then
   allocate(wa_up(nxp,nzp,nt))    ; wa_up = 0.
   allocate(wa_dn(nxp,nzp,nt))    ; wa_dn = 0.
 end if
 if (coopt .eq. 1) then
   allocate(cosd1(nkknew,noonew,nzp))    ; cosd1 = 0.
   allocate(cosd2(nkknew,noonew,nzp))    ; cosd2 = 0.
 end if
 if (psopt .eq. 1) then
   allocate(pwsd1(nkknew,noonew,nzp))    ; pwsd1 = 0.
   allocate(pwsd2(nkknew,noonew,nzp))    ; pwsd2 = 0.
 end if

 do j=1, nzp

! 2-D FFT & determine up/downward
   do n=1, nt
   do i=1, nxp
     u(i,n) = ua(i,j,n)
     w(i,n) = wa(i,j,n)
   enddo
   enddo

   cw_up = (0.d0,0.d0)   ;   cw_dn = (0.d0,0.d0)
   co_up = 0.            ;   co_resi = 0.

   call fft2df(nxp,nt,u,cu)
   call fft2df(nxp,nt,w,cw)

   do n=1, nt
   do i=2, nxp

     ! find m2
     if (i .le. (nxp+1)/2) then
       ii = i
       nn = nt/2+2 - n
       if (n .gt. (nt+1)/2)  nn = nt + (nt/2) - n + 2
     else
       ii = nxp - i + 2
       nn = nt/2 + n
       if (n .gt. (nt+1)/2)  nn = n - (nt+1)/2
     end if

     ! up/dn
     if (m2(ii,nn,j) .gt. 0. .and. m2(ii,nn,j).le.(0.5e-3)**2) then
       mabs = 2.*pi*sqrt(m2(ii,nn,j))

       co = real(real(cu(i,n))*real(cw(i,n))+aimag(cu(i,n))*aimag(cw(i,n)))
       dd = (cu(i,n)*kk(i) + iii*hinv(j)*cw(i,n)) / sign(mabs,kk(i))

       cw_pos = (cw(i,n) + dd) / 2.
       cw_neg = (cw(i,n) - dd) / 2.

       cu_pos = cw_pos * (sign(mabs,kk(i)) - iii*hinv(j)) / kk(i)
       cu_neg = cw_neg * (-sign(mabs,kk(i)) - iii*hinv(j)) / kk(i)

       if (i.ne.1) then
         ccint = oo(n)/kk(i)-ubar(j)
       else if (oo(n).gt.0.0) then
         ccint = 10000.
       else if (oo(n).lt.0.0) then
         ccint = -10000.
       else
         ccint = 0.            ! k,o = 0 ,  and not calculate
       end if
       if (ccint .gt. 0.0) then
         cw_up(i,n) = cw_pos
         cw_dn(i,n) = cw_neg
         co_up(i,n) = &
            real(real(cu_pos)*real(cw_pos)+aimag(cu_pos)*aimag(cw_pos))
       else if (ccint .lt. 0.0) then
         cw_up(i,n) = cw_neg
         cw_dn(i,n) = cw_pos
         co_up(i,n) = &
            real(real(cu_neg)*real(cw_neg)+aimag(cu_neg)*aimag(cw_neg))
       end if
       co_resi(i,n) = co - co_up(i,n)
     end if

   enddo
   enddo

! backward FFT
   if (fieldopt .eq. 1) then
     call fft2db(nxp,nt,cw_up,w_up)
     call fft2db(nxp,nt,cw_dn,w_dn)

     do n=1, nt
     do i=1, nxp
       wa_up(i,j,n) = w_up(i,n) / rhort(zst+1+j)
       wa_dn(i,j,n) = w_dn(i,n) / rhort(zst+1+j)
     enddo
     enddo
   end if

   if (coopt .eq. 1) then
     call cospec(nxp,dx,nt,dt,co_up,  cosd1(:,:,j))
     call cospec(nxp,dx,nt,dt,co_resi,cosd2(:,:,j))
   end if

   if (psopt .eq. 1) then
     call powerspec(nxp,dx,nt,dt,cw_up,pwsd1(:,:,j))
     call powerspec(nxp,dx,nt,dt,cw_dn,pwsd2(:,:,j))
   end if

 enddo

! output
 if (fieldopt .eq. 1) then
   write(fname,'(a)') '/usr/users/kyh/anal/work/wup.nc'
   call out3d(trim(fname),1,(/'w_up'/),wa_up,                          &
               'x',nxp,x/1000.,'z',nzp,z/1000.,'t',nt,t/60.,           &
               'w upward')
   write(fname,'(a)') '/usr/users/kyh/anal/work/wdn.nc'
   call out3d(trim(fname),1,(/'w_dn'/),wa_dn,                          &
               'x',nxp,x/1000.,'z',nzp,z/1000.,'t',nt,t/60.,           &
               'w downward')
 end if

 do i=1, nkknew
   kknew(i) = kk(i)
 end do
 do n=(noonew-1)/2+1, noonew
   oonew(n) = (-1.) * oo(n-(noonew-1)/2)
 end do
 do n=1, (noonew-1)/2
   oonew(n) = oo((noonew-1)/2+2-n)
 end do

 if (coopt .eq. 1) then
   write(fname,'(a)') '/usr/users/kyh/anal/work/kozcoup.nc'
   call out3d(trim(fname),1,(/'Co_up'/),cosd1,                         &
               'k',nkknew,kknew,'o',noonew,oonew,'z',nzp,z/1000.,      &
               'cospectrum upward')
   write(fname,'(a)') '/usr/users/kyh/anal/work/kozcodn.nc'
   call out3d(trim(fname),1,(/'Co_dn'/),cosd2,                         &
               'k',nkknew,kknew,'o',noonew,oonew,'z',nzp,z/1000.,      &
               'cospectrum downward')
   cosd1 = cosd1 + cosd2
   write(fname,'(a)') '/usr/users/kyh/anal/work/kozco.nc'
!   call out3d(trim(fname),1,(/'Co'/),cosd1,                            &
!               'k',nkknew,kknew,'o',noonew,oonew,'z',nzp,z/1000.,      &
!               'cospectrum')
 end if

 if (psopt .eq. 1) then
   write(fname,'(a)') '/usr/users/kyh/anal/work/kozup.nc'
   call out3d(trim(fname),1,(/'PSD_up'/),pwsd1,                        &
               'k',nkknew,kknew,'o',noonew,oonew,'z',nzp,z/1000.,      &
               'PSD upward') 
   write(fname,'(a)') '/usr/users/kyh/anal/work/kozdn.nc' 
   call out3d(trim(fname),1,(/'PSD_dn'/),pwsd2,                        &
               'k',nkknew,kknew,'o',noonew,oonew,'z',nzp,z/1000.,      &
               'PSD downward')
   pwsd1 = pwsd1 + pwsd2
   write(fname,'(a)') '/usr/users/kyh/anal/work/koz.nc'
   call out3d(trim(fname),1,(/'PSD'/),pwsd1,                           &
               'k',nkknew,kknew,'o',noonew,oonew,'z',nzp,z/1000.,      &
               'PSD')
 end if


end


subroutine powerspec(nxp,dx,nt,dt,coeff,ddd)

 implicit none

 integer :: i,n
 integer :: nxp, nt
 real    :: dx, dt
 double complex, dimension(nxp,nt) :: coeff
 real, dimension(nxp/2+1,nt/2*2+1) :: ddd, tem


 do n=1, nt/2
   ddd(1,n) = (dx*dt/nxp/nt)*real( cdabs(coeff(1,(nt+1)/2+n))**2 ) 
   do i=2,nxp/2+1
     ddd(i,n) = (dx*dt/nxp/nt)*real( cdabs(coeff(i,(nt+1)/2+n))**2 & 
                             + cdabs(coeff(nxp-i+2,nt/2+2-n))**2 )
   end do
 end do

 ddd(1,nt/2+1) = (dx*dt/nxp/nt)*real( cdabs(coeff(1,1))**2 ) 
 do i=2,nxp/2+1
   ddd(i,nt/2+1) = (dx*dt/nxp/nt)*real( cdabs(coeff(i,1))**2 & 
                             + cdabs(coeff(nxp-i+2,1))**2 )
 end do

 do n=nt/2+2, nt/2*2+1
   ddd(1,n) = (dx*dt/nxp/nt)*real( cdabs(coeff(1,n-nt/2))**2 ) 
   do i=2,nxp/2+1
     ddd(i,n) = (dx*dt/nxp/nt)*real( cdabs(coeff(i,n-nt/2))**2 & 
                             + cdabs(coeff(nxp-i+2,nt+nt/2+2-n))**2 )
   end do
 end do

 if (nxp .eq. nxp/2*2)  ddd(nxp/2+1,:) = ddd(nxp/2+1,:) / 2.
 if (nt .eq. nt/2*2) then
   ddd(:,1) = ddd(:,1) / 2.
   ddd(:,nt+1) = ddd(:,nt+1) / 2.
 end if
 if (nxp.eq.nxp/2*2 .and. nt.eq.nt/2*2) then
   ddd(nxp/2+1,1) = ddd(nxp/2+1,1) * 2.
   ddd(nxp/2+1,nt+1) = ddd(nxp/2+1,nt+1) * 2.
 end if

 ! sorting
 tem(:,:) = ddd(:,:)
 do n=1, nt/2*2+1
   ddd(:,n) = tem(:,nt/2*2+1+1-n)
 enddo


return
end


subroutine cospec(nxp,dx,nt,dt,co_in,co_out)

 implicit none

 integer :: i, n
 integer :: nxp, nt
 real    :: dx, dt
 real, dimension(nxp,nt) :: co_in
 real, dimension(nxp/2+1,nt/2*2+1) :: co_out, tem


 do n=1, nt/2
   co_out(1,n) = (dx*dt/nxp/nt) * co_in(1,(nt+1)/2+n)
   do i=2,nxp/2+1
     co_out(i,n) = (dx*dt/nxp/nt) * ( co_in(i,(nt+1)/2+n)       &
                             + co_in(nxp-i+2,nt/2+2-n) )
   end do
 end do

 co_out(1,nt/2+1) = (dx*dt/nxp/nt) * co_in(1,1)
 do i=2,nxp/2+1
   co_out(i,nt/2+1) = (dx*dt/nxp/nt) * ( co_in(i,1)             &
                             + co_in(nxp-i+2,1) )
 end do

 do n=nt/2+2, nt/2*2+1
   co_out(1,n) = (dx*dt/nxp/nt) * co_in(1,n-nt/2)
   do i=2,nxp/2+1
     co_out(i,n) = (dx*dt/nxp/nt) * ( co_in(i,n-nt/2)           &
                             + co_in(nxp-i+2,nt+nt/2+2-n) )
   end do
 end do

 if (nxp .eq. nxp/2*2)  co_out(nxp/2+1,:) = co_out(nxp/2+1,:) / 2.
 if (nt .eq. nt/2*2) then
   co_out(:,1) = co_out(:,1) / 2.
   co_out(:,nt+1) = co_out(:,nt+1) / 2.
 end if
 if (nxp.eq.nxp/2*2 .and. nt.eq.nt/2*2) then
   co_out(nxp/2+1,1) = co_out(nxp/2+1,1) * 2.
   co_out(nxp/2+1,nt+1) = co_out(nxp/2+1,nt+1) * 2.
 end if

 ! sorting
 tem(:,:) = co_out(:,:)
 do n=1, nt/2*2+1
   co_out(:,n) = tem(:,nt/2*2+1+1-n)
 enddo


return
end


