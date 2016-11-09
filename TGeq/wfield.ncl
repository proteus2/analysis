load "/data3/kyh/util/ncl/header.ncl"

xspace = 170  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "/data3/kyh/util/ncl/figsize.ncl"
load "/data3/kyh/util/ncl/attribute.ncl"

axisfontratio = 1.2

begin

 wks = cr_wks("wfield",600,0)
 gsn_define_colormap(wks,"hotcold_18lev")
 icol_center = 11
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- parameters
 x0 = 0.
 x1 = 1000.
 z0 = 20.
 nx = 301
 lx = 500.

 ns = 2
 nw = 3
 
;--------------------------------------- read data
 f1 = addfile("ctl_0.nc" ,"r")
 f2 = addfile("npol_0.nc","r")

 z   = f1->z + z0
 nz  = dimsizes(z)

 wc = new((/ns,nw,2,nz/),float)
 wc(0,0,:,:) = (/ f1->wr,    f1->wi    /)
 wc(0,1,:,:) = (/ f1->wr_nm, f1->wi_nm /)
 wc(0,2,:,:) = (/ f1->wr_pm, f1->wi_pm /)
 wc(1,0,:,:) = (/ f2->wr,    f2->wi    /)
 wc(1,1,:,:) = (/ f2->wr_nm, f2->wi_nm /)
 wc(1,2,:,:) = (/ f2->wr_pm, f2->wi_pm /)

 x = fspan(x0,x1,nx)
 kx = x*(2.*acos(-1.)/lx)

 w = new((/ns,nw,nz,nx/),float)
 do is=0, ns-1
 do iw=0, nw-1
 do k=0, dimsizes(z)-1
   w(is,iw,k,:) = wc(is,iw,0,k)*cos(kx(:))-wc(is,iw,1,k)*sin(kx(:))
 end do
 end do
 end do

;--------------------------------------- line attribute
cnlevw = (/-1.5,-1.2,-0.9,-0.6,-0.3,0.3,0.6,0.9,1.2,1.5/)
cnfilw = icol_center + (/-10,-8,-6,-5,-4,0,2,3,4,6,9/)

;--------------------------------------- plot

 wi = 43 * (x1-x0)/2200.
 he = 43
 xf = 10 + (/0,1,2/)*28 ;+ (/0,6,6/)
 yf = 10 + (/0,1/)*68

 mstr = (/"w","uw","u","w"/)
 mtail = (/""," (smooth N:S:2:N:)"/)
 ystr = (/"z [km]",""/)

do is=0, ns-1
do iw=0, 0

 contour = cnshadeplot(wks,x/1.e2,z,w(is,iw,:,:), \
           cnlevw,cnfilw,"")

 set_size(contour,xf(is),yf(0),wi,he)
 axis_range(contour,-999,999,-999,999,False,False)
 axis_str(contour,mstr(iw)+mtail(is),"x [10:S:2:N: km]",ystr(is))

 cnfinalize(contour,wks)

 draw(contour)

end do
end do

;dr_txt(wks,"(a)",xf(0)-2,yf(0)-4,25,8,0.)
;dr_txt(wks,"(b)",xf(1)-2,yf(0)-4,25,8,0.)
;dr_txt(wks,"(c)",xf(0)-2,yf(1)-4,25,8,0.)
;dr_txt(wks,"(d)",xf(1)-2,yf(1)-4,25,8,0.)


 frame(wks)

end

