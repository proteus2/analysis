load "$UTIL/ncl/header.ncl"

begin

 f0 = addfile("../dat/xyz_ocn.nc","r")
 f1 = addfile("../dat/mask_ocn.nc","r")

 x  = f0->lon
 y  = f0->lat
 z  = f0->lev
 yb = f0->lat_bnds
 to = f1->thetao(0,:,:,:)

 nx = dimsizes(x)
 ny = dimsizes(y)
 nz = dimsizes(z)

 wa4grds = new((/nz,ny/),float)
 wa4grds(:,:) = 0.

 wa4grds(0,:) = (yb(:,1)-yb(:,0))*(x(1)-x(0))*cos(y(:)/180.*3.141592)
 wa4grds(0,0   ) = 0.
 wa4grds(0,ny-1) = 0.
 do k=1, nz-1
   wa4grds(k,:) = wa4grds(0,:)
 end do

 nmask = new((/nz,ny/),integer)
 nmask(:,:) = nx
 do k=0, nz-1
 do j=0, ny-1
 do i=0, nx-1
   if ( ismissing( to(k,j,i) ) ) then
     nmask(k,j) = nmask(k,j) - 1
   end if
 end do
 end do
 end do

 do k=0, nz-1
   tot = dim_sum(wa4grds(k,:)*int2flt(nmask(k,:)))
   wa4grds(k,:) = wa4grds(k,:)/tot
 end do

 wa4xavg = wa4grds(:,:)*int2flt(nmask(:,:))

 system("mv area_ocn.nc area_ocn.nc.old")
 f10 = addfile("area_ocn.nc","c")

 y!0 = "latitude"
 z!0 = "depth"

 wa4grds!1 = "latitude"
 wa4grds!0 = "depth"
 wa4xavg!1 = "latitude"
 wa4xavg!0 = "depth"

 wa4grds&latitude = y
 wa4grds&depth    = z
 wa4xavg&latitude = y
 wa4xavg&depth    = z

 f10->latitude = y
 f10->depth    = z
 f10->wa4grds  = wa4grds
 f10->wa4xavg  = wa4xavg

end

