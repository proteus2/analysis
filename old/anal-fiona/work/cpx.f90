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
!
! k0int  : integrate from Lx / 2^(k0int)
!-------------------------------------------------------------------------

 implicit none

 integer, parameter :: nxp =1200, nzp = 350, nt = 300, power = 0, k0int = 2
 integer, parameter :: intopt = 1, zst = 0, ntall = 480 ! for reading basic wind
 real,    parameter :: dx = 1000., dt = 60.
 real,    parameter :: dcc =2.0, cmin = -200., cmax = 200.
! real*8,  parameter :: coef_norm = 0.201648354312960   ! Jan
! real*8,  parameter :: coef_norm = 0.201555570315203   ! Apr
 real*8,  parameter :: coef_norm = 0.201349469516055   ! Jul line 52
! real*8,  parameter :: coef_norm = 0.203587889232740   ! Oct
!-------------------------------------------------------------------------
 integer, parameter :: ncc = (cmax-cmin)/dcc+1 
 integer, parameter :: nkk = nxp/2+1, noo = nt/2*2+1
 integer :: i,j,n,p, ncid
 integer :: swt, cindex(nkk,noo,2), pp1, pp2
 real    :: kk(nkk), oo(noo), cc(ncc), z_km(nzp), zout(nzp+2)
 real    :: dkk, doo, dkkh, dooh, dcch, gcmid
 real, dimension(nkk,noo)       :: gcmin, gcmax
 real, dimension(nkk,noo,nzp)   :: pwsd
 real, dimension(ncc,nzp)       :: psdc
 real, dimension(ncc,nzp+2)     :: spout
 real, dimension(0:ntall+2,nzp+2)  :: ubardat
 real, dimension(nzp)              :: ubar, ubar1, ubar2

 character*128 :: idname, ifname, ivarname, dname, fname, varname 


 idname   = '../'
 ifname   = 'kozcodn' !'kozcosrcup'
 ivarname = 'Co_dn' !'Co_up'
 dname   = '../' !trim(idname)
 fname   = 'cpxmomdn' !'cpxmomsrcup' 
 varname = 'mom_dn' !'mom_up' !trim(ivarname)


 call opennc(trim(idname)//trim(ifname)//'.nc',ncid)
 call get3d(ncid, trim(ivarname), nkk, noo, nzp, pwsd)
 call get1d(ncid, 'k', nkk, kk)
 call get1d(ncid, 'o', noo, oo)
 call get1d(ncid, 'z', nzp, z_km)
 call closenc(ncid)

! read ubar or umean
 call opennc('/export23/kyh/105jul/unbar.nc',ncid)
 call get2d(ncid, 'ubar', ntall+3, nzp+2, ubardat)
 call closenc(ncid)
 do j=1, nzp
   ubar1(j) = ubardat(0,zst+1+j)
   ubar2(j) = ubardat(ntall+2,zst+1+j)     ! avg. over selected time
 enddo
 ubar = ubar1


 if (power.eq.1) pwsd = abs(pwsd)

 do p=1, ncc 
   cc(p) = cmin + dcc*(p-1)
 enddo
 dcch = dcc/2.

 dkk  = kk(2)-kk(1)
 doo  = oo(2)-oo(1)
 dkkh = dkk/2.
 dooh = doo/2.
 
 do n=1, noo/2
 do i=2, nkk
   gcmin(i,n) = (oo(n)-dooh) / (kk(i)-dkkh)
   gcmax(i,n) = (oo(n)+dooh) / (kk(i)+dkkh)
 enddo
 enddo
 do i=2, nkk
   gcmin(i,noo/2+1) = -dooh / (kk(i)-dkkh)
   gcmax(i,noo/2+1) = dooh / (kk(i)-dkkh)
 enddo
 do n=noo/2+2, noo
 do i=2, nkk
   gcmin(i,n) = (oo(n)-dooh) / (kk(i)+dkkh)
   gcmax(i,n) = (oo(n)+dooh) / (kk(i)-dkkh)
 enddo
 enddo

 
 cindex(:,:,1) = ncc
 cindex(:,:,2) = 1
 do n=1, noo 
 do i=2, nkk
   gcmid = oo(n) / kk(i)
   swt = 0
   do p=1, ncc
     if ( (cc(p)+dcch).gt.gcmax(i,n) .and. gcmid.lt.(cc(p)-dcch) )  EXIT

     if ( (cc(p).ge.gcmin(i,n) .and. cc(p).le.gcmax(i,n)) .or. &
          (gcmid.ge.(cc(p)-dcch) .and. gcmid.lt.(cc(p)+dcch)) ) then
       swt = swt + 1
       if (swt .eq. 1)  cindex(i,n,1) = p
     end if
   enddo
   cindex(i,n,2) = cindex(i,n,1) - 1 + swt
   if (cindex(i,n,2) .gt. ncc)  cindex(i,n,2) = ncc
 enddo
 enddo


 psdc = 0.
 do j=1, nzp
 do n=1, noo
 do i=k0int+2, nkk

   gcmid = oo(n) / kk(i)
   pp1 = cindex(i,n,1)
   pp2 = cindex(i,n,2)

   ! near frequency = 0
   if (intopt.eq.1 .and. pp2.ge.pp1) then
     if (gcmid .gt. ubar(j)) then
       pp1 = max(cindex(i,n,1),int((ubar(j)-0.0001-cmin)/dcc+1)+1)
       pp1 = min(pp1,cindex(i,n,2))
       pwsd(i,n,j) = pwsd(i,n,j) * (cindex(i,n,2)-pp1+1) / (cindex(i,n,2)-cindex(i,n,1)+1)
     else if (gcmid .lt. ubar(j)) then
       pp2 = min(cindex(i,n,2),int((ubar(j)-cmin)/dcc+1))
       pp2 = max(pp2,cindex(i,n,1))
       pwsd(i,n,j) = pwsd(i,n,j) * (pp2-cindex(i,n,1)+1) / (cindex(i,n,2)-cindex(i,n,1)+1)
     else
       pwsd(i,n,j) = 0.
     end if
   end if

   do p=pp1, pp2
     psdc(p,j) = psdc(p,j) + pwsd(i,n,j) / (pp2-pp1+1)
   enddo
 enddo
 enddo
 enddo
 psdc(:,:) = psdc(:,:)*dkk*doo/dcc
 psdc = real(psdc * coef_norm)


! output
 do j=2, nzp+1
   spout(:,j) = psdc(:,j-1)
   zout(j) = z_km(j-1)
 enddo
 spout(:,1) = spout(:,2)
 spout(:,nzp+2) = spout(:,nzp+1)
 zout(1) = zout(2)*2. - zout(3)
 zout(nzp+2) = zout(nzp+1)*2. - zout(nzp)


 call out2d(trim(dname)//trim(fname)//'.nc',1,(/trim(varname)/),spout, &
            'c',ncc,cc,'z',nzp+2,zout,'phase speed spectrum')


end
