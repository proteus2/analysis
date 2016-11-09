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
 integer :: swt, cindex(nkk,noo,2)
 real    :: kk(nkk), oo(noo), cc(ncc), z_km(nzp), zout(nzp+2)
 real    :: dkk, doo, dkkh, dooh, dcch, summ
 real, dimension(nkk,noo)       :: gcmin, gcmax
 real, dimension(nkk,noo,nzp)   :: pwsd
 real, dimension(ncc,nzp)       :: psdc
 real, dimension(ncc,nzp+2)     :: spout
 character*128 :: idname, ifname, ivarname, dname, fname, varname 


 idname   = '/export22/kyh/ing/105jul/spec/koz/'
 ifname   = 'koz'
 ivarname = 'PSD'
 dname   = '/export22/kyh/ing/105jul/spec/cpx/' !trim(idname)
 fname   = trim(ifname)//'_cpx' 
 varname = trim(ivarname)


 call opennc(trim(idname)//trim(ifname)//'.nc',ncid)
 call get3d(ncid, trim(ivarname), nkk, noo, nzp, pwsd)
 call get1d(ncid, 'k', nkk, kk)
 call get1d(ncid, 'o', noo, oo)
 call get1d(ncid, 'z', nzp, z_km)
 call closenc(ncid)

 if (power.eq.1) pwsd = abs(pwsd)

 do p=1, ncc 
   cc(p) = cmin + dcc*(p-1)
 enddo
 dcch = 0. !dcc/2.

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
   gcmin(i,noo/2+1) = (oo(noo/2+1)-dooh) / (kk(i)-dkkh)
   gcmax(i,noo/2+1) = (oo(noo/2+1)+dooh) / (kk(i)-dkkh)
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
   swt = 0
   do p=1, ncc
     if ( (cc(p)+dcch) .gt. gcmax(i,n) )  EXIT

     if ( (cc(p)-dcch).ge.gcmin(i,n) .and. (cc(p)+dcch).le.gcmax(i,n) ) then
       swt = swt + 1
       if (swt .eq. 1)  cindex(i,n,1) = p
     end if
   enddo
   cindex(i,n,2) = cindex(i,n,1) - 1 + swt
 enddo
 enddo


 psdc = 0.
 do j=1, nzp
 do n=1, noo
 do i=2, nkk
   do p=cindex(i,n,1), cindex(i,n,2)
     psdc(p,j) = psdc(p,j) + pwsd(i,n,j)/(cindex(i,n,2)-cindex(i,n,1)+1)
   enddo
 enddo
 enddo
 enddo
 psdc(:,:) = psdc(:,:)*dkk*doo/dcc
 psdc = real(psdc * coef_norm)

 summ = 0.
 do p=1, ncc
   summ = summ + psdc(p,55)
 enddo
 print *, summ

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
