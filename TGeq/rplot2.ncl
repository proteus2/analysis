load "$UTIL/ncl/header.ncl"

xspace = 170  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 1.2

begin

 wks = cr_wks("rplot2",600,0)
 gsn_define_colormap(wks,"hotcold_18lev")
 icol_center = 11
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- parameters
 ns = 1
 
;--------------------------------------- read data
 f1 = addfile("ctl_2.nc" ,"r")

 phi   = f1->R_phi
 r_abs = f1->R_abs

 nphi = dimsizes(phi)
 nr   = dimsizes(r_abs)

 rb = new((/ns,nr,nphi/),float)
 do i=0, nr-1
   rb(0,i,:) = f1->r(0,:,i)
 end do

;--------------------------------------- line attribute
cnlevw = (/0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16/)*0.05
cnfilw = icol_center + (/-10,-8,-6,-5,-4,0,2,3,4,6,9/)
; labels = "|R:B:top:N:| = "+(/0.2,0.25,0.3/)
; ldsh = (/0,1,2/)

;--------------------------------------- plot

 wi = 43
 he = 43
 xf = 10 + (/0,1,2/)*28 ;+ (/0,6,6/)
 yf = 10 + (/0,1/)*68

 mstr = (/"w","uw","u","w"/)
 mtail = (/""," (smooth N:S:2:N:)"/)
 ystr = (/"z [km]",""/)

do is=0, ns-1

 contour = cnlineplot(wks,phi,r_abs,rb(is,:,:), \
           cnlevw,1)

 set_size(contour,xf(is),yf(is),wi,he)
 axis_range(contour,0,360,0,0.6,False,False)
 axis_str(contour,"|R_20 km|","phase diff. (deg)","|R_45 km|")
 setvalues contour
  "cnLineLabelFontHeightF" : p2n(3.6)
  "tmXBMode"         : "Manual"
  "tmXBTickSpacingF" : 90
 end setvalues
 cnfinalize(contour,wks)

 draw(contour)

end do

;dr_txt(wks,"(a)",xf(0)-2,yf(0)-4,25,8,0.)
;dr_txt(wks,"(b)",xf(1)-2,yf(0)-4,25,8,0.)
;dr_txt(wks,"(c)",xf(0)-2,yf(1)-4,25,8,0.)
;dr_txt(wks,"(d)",xf(1)-2,yf(1)-4,25,8,0.)


 frame(wks)

end

