load "$UTIL/ncl/header.ncl"

xspace = 170  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 1.0

begin

 wks = cr_wks("fig4",600,0)
 colind = NhlNewColor(wks,0.6,0.6,0.6)
 colind2 = NhlNewColor(wks,0.8,0.8,0.8)

;--------------------------------------- parameters
 nf = 1
 nl = 4
;--------------------------------------- read data

 ; read variables to plot
 f1x = addfile("../res/wind/mfa/mfa-klf_p1_f1.nc","r")

 f1k = addfile("../res/temp2/mfa/mfa-klf_p1_f1.nc","r")


 z  = f1k->height({20000:}) / 1.e3
 nz = dimsizes(z)
;;;
 var_dat = new((/nf,nl,nz/),float,1.e32)
 var_dat(0,0,:) = f1x->mfx_e({20000:}) * 1.e3

 var_dat(0,1,:) = f1k->mfx_e({20000:}) * 1.e3

 var_dat(0,2,:) = f1k->mfx_up_e({20000:}) * 1.e3

 var_dat(0,3,:) = f1k->mfx_dn_e({20000:}) * 1.e3
;;;

print(-var_dat(0,1,:)/var_dat(0,0,:)*100)
print(-sum(var_dat(0,1,:))/sum(var_dat(0,0,:))*100)

;--------------------------------------- process 1
 var = var_dat

;--------------------------------------- line attribute
 labels = (/"F:B:W:N:","F:B:T total:N:","F:B:T up:N:","F:B:T down:N:"/)
 lcol = (/1,2,2,2/)
 ltck = (/0.6,1.3,1.0,1.0/)
 ldsh = (/0,0,2,16/)

;--------------------------------------- plot

 wi = 33
 he = 33
 xf = 10 + (/0,1/)*(wi+6)
 yf = 10 + (/0,0/)*(he+10)

 ystr = (/"z [km]","","","z [km]","",""/)

;do ifig=0, nf-1

 plot = xylineplot(wks,var(0,:,:),z,"D",ldsh,"T",ltck,"C",lcol)

 xyaddline(plot,0,z,"D",0,"T",0.6,"C",colind2)

 set_size(plot,xf(0),yf(0),wi,he)
 axis_range(plot,-999,2,-999,45,False,False)
 axis_str(plot,"","F [mPa]",ystr(0))
; tick_exp(plot,"XB",xval,xlab,time)

 draw(plot)


 lg1 = lgplot(wks,labels,"D",ldsh, \
                  "T",ltck,"C",lcol,1)

 set_size_lg(lg1,xf(0)+wi-10,yf,10,14)

 draw(lg1)


 plot2 = xylineplot(wks,var(0,3,:),z,"D",ldsh(3),"T",ltck(3),"C",lcol(3))

 xyaddline(plot2,0,z,"D",0,"T",0.6,"C",colind2)

 set_size(plot2,xf(1),yf(1),wi,he)
 axis_range(plot2,-0.12,0.01,-999,45,False,False)
 axis_str(plot2,"","F [mPa]",ystr(1))
; tick_exp(plot2,"XB",xval,xlab,time)

 draw(plot2)


; dr_txt(wks,"date",xf+wi/2,yf+he+10,22,6,0.)

;end do

dr_txt(wks,"(a)",xf(0)-2,yf(0)-4,25,8,0.)
dr_txt(wks,"(b)",xf(1)-2,yf(1)-4,25,8,0.)


 frame(wks)

end

