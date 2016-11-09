undef("a1dto2d")
function a1dto2d(u:numeric)
local uu, u, ndim
begin

 ndim = dimsizes(u)
 if ( dimsizes(ndim) .eq. 2 ) then
   uu = u
   uu@nz = ndim(0)
 else
   uu = onedtond(u,(/1,ndim/))
   uu@nz = 1
 end if

return(uu)
end

;-------------------------------------------------------------------------------

undef("a2dto1d")
function a2dto1d(uu:numeric)
local u, uu, ndim
begin

 ndim = dimsizes(uu)
 if ( ndim(0) .eq. 1 ) then
   u = uu(0,:)
 else
   u = uu
 end if

return(u)
end

;-------------------------------------------------------------------------------

undef("seasonal")
function seasonal(u:numeric, n_harmonic:integer)
local us, u, n_harmonic, uu, uus0, uus, fc, nz
begin

 uu = a1dto2d(u)
 nz = uu@nz

 uus0 = new((/nz,12/),float)
 do n=0, 11
   uus0(:,n) = dim_avg(uu(:,n::12))
 end do

 if (n_harmonic .eq. 0) then
   uus = uus0
 else
   fc = ezfftf(uus0)
   fc(:,:,n_harmonic:) = 0.
   uus = ezfftb(fc,fc@xbar)
 end if

 us = a2dto1d(uus)

return(us)
end

;-------------------------------------------------------------------------------

undef("deseasonalize")
function deseasonalize(u:numeric, n_harmonic:integer, period:numeric)
local ud, u, n_harmonic, period, uu, t, uud, uus, p1, p2, p2d,  \
      np, ndim, nz, nt, k, n, ip, i0_p
begin

 np = dimsizes(period)-1

 uu = a1dto2d(u)
 ndim = dimsizes(uu)
 nz = ndim(0)
 nt = ndim(1)
 t = uu&$uu!1$
 t!0 = "t"
 t&t = t

 uud = uu
 if (np .eq. 0) then
   if (period .eq. 0) then
     uus = seasonal(uu,n_harmonic)
     do k=0, nz-1
     do n=0, 11
       uud(k,n::12) = uu(k,n::12) - uus(k,n)
     end do
     end do
   else
     print("Check 'period' in the function 'deseasonalize'.")
     exit
   end if
 else
   do ip=0, np-1
     p1 = period(ip)
     p2 = period(ip+1) - 1.e-6
     p2d = p2
     if ( ip .eq. np-1 .and. p2 .lt. t(nt-1) ) then
       p2d = t(nt-1)
     end if
     uus = seasonal(uu(:,{p1:p2}),n_harmonic)
     i0_p = ind( t .ge. p1 .and. t .le. p2d )
     do k=0, nz-1
     do n=0, 11
       uud(k,i0_p(n::12)) = uu(k,i0_p(n::12)) - uus(k,n)
     end do
     end do
     if ( ip .eq. 0 .and. p1 .gt. t(0) ) then
       do k=0, nz-1
       do n=0, 11
         uud(k,i0_p(0)-1-n:0:12) = uu(k,i0_p(0)-1-n:0:12) - uus(k,11-n)
       end do
       end do
     end if
     delete(i0_p)
   end do
 end if

 ud = a2dto1d(uud)

return(ud)
end

