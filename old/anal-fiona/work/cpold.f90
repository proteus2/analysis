program phase_speed_spectrum_wrt_z

 use netcdfio

!-------------------------------------------------------------------------
!
! OUTPUT : Cx, spectral density         (cpx)
!          summation number w.r.t. Cx for determining dcc  (checkno)
!
! INPUT  : k, omega, 
!          k_omega 2D PSD w.r.t. z
!
! power  : consider only absolute value of power spectra
!-------------------------------------------------------------------------

 implicit none

 integer, parameter :: nxp =1200, nzp = 350, nt = 360, power = 0
 real,    parameter :: dx = 1000., dt = 60.
 real,    parameter :: dcc = 1.0, cmin = -200., cmax = 200.
! real*8,  parameter :: coef_norm = 0.201648354312960   ! Jan
! real*8,  parameter :: coef_norm = 0.201555570315203   ! Apr
 real*8,  parameter :: coef_norm = 0.201349469516055   ! Jul line 52
! real*8,  parameter :: coef_norm = 0.203587889232740   ! Oct
!-------------------------------------------------------------------------
 integer, parameter :: ncc = (cmax-cmin)/dcc+1 
 integer, parameter :: nkk = nxp/2+1, noo = nt/2*2+1
 integer :: i,j,n,p, ncid
 integer :: cindex(nkk,noo,2)
 real    :: kk(nkk), oo(noo), cc(ncc), z_km(nzp), zout(nzp+2)
 real    :: dkk, doo, ccmid, summ
 real, dimension(nkk,noo,nzp)   :: pwsd
 real, dimension(ncc,nzp)       :: pwsdc
 real, dimension(ncc,nzp+2)     :: spout
 character*128 :: ifname, fname, ivarname, varname 


 ifname   = '/export22/kyh/ing/105jul/spec/koz/koz.nc'
! ifname   = '/export22/kyh/ing/kozcodn.nc'
! ifname = '/export24/kyh/kozcoresi.nc'
 ivarname = 'PSD'
 varname  = 'PSD'


 call opennc(trim(ifname),ncid)
 call get3d(ncid, trim(ivarname), nkk, noo, nzp, pwsd)
 call get1d(ncid, 'k', nkk, kk)
 call get1d(ncid, 'o', noo, oo)
 call get1d(ncid, 'z', nzp, z_km)
 call closenc(ncid)

 if (power.eq.1) pwsd = abs(pwsd)

 do p=1, ncc
   cc(p) = cmin + dcc*(p-1)
 enddo
 dkk = kk(2)-kk(1)
 doo = oo(2)-oo(1)
 
 cindex(:,:,1) = 999
 cindex(:,:,2) = 0
 do n=1, noo 
 do i=2, nkk
   ccmid = oo(n)/kk(i)
   do p=1, ncc
     if ( ccmid .lt. (cc(p)-dcc/2.) ) EXIT

     if ((ccmid.gt.(cc(p)-dcc/2.)).and.(ccmid.lt.(cc(p)+dcc/2.))) then
       cindex(i,n,1) = p
       cindex(i,n,2) = p
     else if ( ccmid .eq. (cc(p)-dcc/2.) ) then
       cindex(i,n,1) = p-1
       cindex(i,n,2) = p
     else if ( ccmid .eq. (cc(p)+dcc/2.) ) then
       cindex(i,n,1) = p
       cindex(i,n,2) = p+1
     end if
   enddo
 enddo
 enddo


 pwsdc = 0.
 do j=1, nzp
 do n=1, noo
 do i=2, nkk
   do p=cindex(i,n,1), cindex(i,n,2)
     pwsdc(p,j) = pwsdc(p,j) + pwsd(i,n,j)/(cindex(i,n,2)-cindex(i,n,1)+1)
   enddo
 enddo
 enddo
 enddo
 pwsdc(:,:) = pwsdc(:,:)*dkk*doo/dcc
 pwsdc = real(pwsdc * coef_norm)

 summ = 0.
 do p=1, ncc
   summ = summ + pwsdc(p,55)
 enddo
 print*, summ

! output
 do j=2, nzp+1
   spout(:,j) = pwsdc(:,j-1)
   zout(j) = z_km(j-1)
 enddo
 spout(:,1) = spout(:,2)
 spout(:,nzp+2) = spout(:,nzp+1)
 zout(1) = zout(2)*2. - zout(3)
 zout(nzp+2) = zout(nzp+1)*2. - zout(nzp)


 write(fname,'(a)') '/usr/users/kyh/anal/work/cpx2.nc'
 call out2d(trim(fname),1,(/trim(varname)/),spout,'c',ncc,cc,'z',nzp+2,zout, &
             'phase speed spectrum')


end
