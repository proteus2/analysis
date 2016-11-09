load "/data3/kyh/util/ncl/header.ncl"

xspace = 170  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "/data3/kyh/util/ncl/figsize.ncl"
load "/data3/kyh/util/ncl/attribute.ncl"

axisfontratio = 1.2

begin

 wks = cr_wks("rplot",600,0)
; gsn_define_colormap(wks,"hotcold_18lev")
; icol_center = 11
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- parameters
 ns = 1
 i_r_abs = (/4,5,6/)
 
;--------------------------------------- read data
 f1 = addfile("ctl_2.nc" ,"r")

 phi   = f1->R_phi
 r_abs = f1->R_abs(i_r_abs)

 nphi = dimsizes(phi)
 nr   = dimsizes(i_r_abs)

 rbi = new((/ns,nphi,nr/),float)
 rbi(0,:,:) = f1->r(0,:,i_r_abs)

 rb = new((/ns,nr,nphi/),float)
 do i=0, nr-1
   rb(:,i,:) = rbi(:,:,i)
 end do

;--------------------------------------- line attribute
 labels = "|R:B:top:N:| = "+(/0.2,0.25,0.3/)
 ldsh = (/0,1,2/)

;--------------------------------------- plot

 wi = 43
 he = 43
 xf = 10 + (/0,1,2/)*28 ;+ (/0,6,6/)
 yf = 10 + (/0,1/)*68

 mstr = (/"w","uw","u","w"/)
 mtail = (/""," (smooth N:S:2:N:)"/)
 ystr = (/"z [km]",""/)

do is=0, ns-1

 plot = xylineplot(wks,phi,rb(is,:,:),"D",ldsh,"T",1.0,"C",1)

 set_size(plot,xf(is),yf(is),wi,he)
 axis_range(plot,0,360,-999,999,False,False)
 axis_str(plot,"","phase diff. (deg)","|R:B:bottom:N:|")
; tick_exp(plot,"XB",xval,xlab,time)

 draw(plot)


 lg1 = lgplot(wks,labels,"D",ldsh,"T",1.0,"C",1,1)

 set_size_lg(lg1,xf(is)+12,yf(is),10,12)

 draw(lg1)

end do

;dr_txt(wks,"(a)",xf(0)-2,yf(0)-4,25,8,0.)
;dr_txt(wks,"(b)",xf(1)-2,yf(0)-4,25,8,0.)
;dr_txt(wks,"(c)",xf(0)-2,yf(1)-4,25,8,0.)
;dr_txt(wks,"(d)",xf(1)-2,yf(1)-4,25,8,0.)


 frame(wks)

end

