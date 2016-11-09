load "$UTIL/ncl/header.ncl"

xspace = 83   ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 1.2

begin

 wks = cr_wks("fig6",600,0)
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
 f1 = addfile("../../TGeq/dat/ctl_0.nc" ,"r")
 f2 = addfile("../../TGeq/dat/npol_0.nc","r")
 f0 = addfile("../res/uvn_background2_p1_region.nc","r")

 zn = f0->height_m/1.e3
 nb = f0->nbavg
 nbs= (/nb,nb-dtrend(nb,0)/)
 
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

 f3 = addfile("../../TGeq/dat/ctl_2.nc" ,"r")

 phi   = f3->R_phi
 r_abs = f3->R_abs

 nphi = dimsizes(phi)
 nr   = dimsizes(r_abs)

 rb = new((/nr,nphi/),float)
 do i=0, nr-1
   rb(i,:) = f3->r(0,:,i)
 end do
print(rb(0,(/3,10/)))

;--------------------------------------- line attribute
cnlevw = (/-1.5,-1.2,-0.9,-0.6,-0.3,0.3,0.6,0.9,1.2,1.5/)
cnfilw = icol_center + (/-10,-8,-6,-5,-4,0,2,3,4,6,9/)

cnlevr = (/0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16/)*0.05

;--------------------------------------- plot

 wi = 43 * (x1-x0)/2200.
 he = 43
 xf = 22 + (/0,1,0/)*40 ;+ (/0,6,6/)
 yf = 20 + (/0,0,1/)*62
 wi2 = 28+wi
 he2 = 30
 wi0 = 10

 mstr = (/"N","smooth N","w","w"/)
 ystr = (/"z [km]","","",""/)
 ylab = (/True,False,False,False/)


do is=0, ns-1
do iw=0, 0

 plot = xylineplot(wks,nbs(is,:)*1.e2,zn,"D",0,"T",1.2,"C",1)

 set_size(plot,xf(is)-13,yf(is),wi0,he)
 axis_range(plot,1.8,2.7,20,45,False,False)
 axis_str(plot,mstr(is),"[10:S:-2:N: s:S:-1:N:]",ystr(is))
 setvalues plot
  "tmYLLabelsOn" : ylab(is)
  "tmXBMaxTicks" : 2
 end setvalues

 draw(plot)


 contour = cnshadeplot(wks,x/1.e2,z,w(is,iw,:,:), \
           cnlevw,cnfilw,"")

 set_size(contour,xf(is),yf(is),wi,he)
 axis_range(contour,-999,999,-999,999,False,False)
 axis_str(contour,mstr(is+2),"x [10:S:2:N: km]",ystr(is+2))
 setvalues contour
  "tmYLLabelsOn" : ylab(is+2)
 end setvalues

 cnfinalize(contour,wks)

 draw(contour)

end do
end do


 contour3 = cnlineplot(wks,phi,r_abs,rb, \
           cnlevr,1)

 set_size(contour3,xf(2),yf(2),wi2,he2)
 axis_range(contour3,0,360,0,0.6,False,False)
 axis_str(contour3,"|R| (20 km)",":F33:Dq:F22: (deg)","|R| (45 km)")
 setvalues contour3
  "cnLineLabelFontHeightF" : p2n(3.6)
  "tmXBMode"         : "Manual"
  "tmXBTickSpacingF" : 90
 end setvalues
 cnfinalize(contour3,wks)

; draw(contour3)


dr_txt(wks,"(a)",xf(0)-12,yf(0)-9,25,8,0.)
dr_txt(wks,"(b)",xf(1)-12,yf(1)-9,25,8,0.)
;dr_txt(wks,"(c)",xf(2)-6,yf(2)-4,25,8,0.)
;dr_txt(wks,"(d)",xf(1)-2,yf(1)-4,25,8,0.)


 frame(wks)

end

