load "/data3/kyh/util/ncl/header.ncl"

xspace = 170  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "/data3/kyh/util/ncl/figsize.ncl"
load "/data3/kyh/util/ncl/attribute.ncl"

axisfontratio = 1.4

begin

 wks = cr_wks("mf-wt_z",300,0)
 colind = NhlNewColor(wks,0.6,0.6,0.6)
 colind2 = NhlNewColor(wks,0.8,0.8,0.8)

;--------------------------------------- parameters
 nf = 3
 nl = 2
;--------------------------------------- read data

 ; read variables to plot
 f1x = addfile("../res/wind/mfa/mfa-klf_p1_f1.nc","r")
 f1k = addfile("../res/temp/mfa/mfa-klf_p1_f1.nc","r")
 f2x = addfile("../res/wind/mfa/mfa-klf_p2_f1.nc","r")
 f2k = addfile("../res/temp/mfa/mfa-klf_p2_f1.nc","r")
 f3x = addfile("../res/wind/mfa/mfa-klf_p3_f1.nc","r")
 f3k = addfile("../res/temp/mfa/mfa-klf_p3_f1.nc","r")


 z  = f1k->height({20000:}) / 1.e3
 nz = dimsizes(z)
;;;
 var_dat = new((/nf,nl,nz/),float,1.e32)
 var_dat(0,0,:) = f1x->mfx({20000:}) * 1.e3
 var_dat(0,1,:) = f1k->mfx({20000:}) * 1.e3
 var_dat(1,0,:) = f2x->mfx({20000:}) * 1.e3
 var_dat(1,1,:) = f2k->mfx({20000:}) * 1.e3
 var_dat(2,0,:) = f3x->mfx({20000:}) * 1.e3
 var_dat(2,1,:) = f3k->mfx({20000:}) * 1.e3
;;;

;--------------------------------------- process 1
 var = var_dat

;--------------------------------------- line attribute
 labels = (/"F:B:W:N:","F:B:T:N:"/)
 lcol = (/1,2,2,2/)
 ltck = (/1.3,1.3,1.0,1.0/)
 ldsh = (/0,0,2,16/)

;--------------------------------------- plot

 wi = 43
 he = 43
 xf = 12 + (/0,1,2/)*(wi+10)
 yf = 10 + (/0,0,0/)*(he+10)

 ystr = (/"z [km]","","","z [km]","",""/)

do ifig=0, nf-1

 plot = xylineplot(wks,var(ifig,:,:),z,"D",ldsh,"T",ltck,"C",lcol)

 xyaddline(plot,0,z,"D",0,"T",0.6,"C",colind2)

 set_size(plot,xf(ifig),yf(ifig),wi,he)
 axis_range(plot,-999,999,-999,45,False,False)
 axis_str(plot,"","F [mPa]",ystr(ifig))
; tick_exp(plot,"XB",xval,xlab,time)

 draw(plot)

 lg1 = lgplot(wks,labels,"D",ldsh, \
                  "T",ltck,"C",lcol,1)

 set_size_lg(lg1,xf(ifig)+wi-10,yf,10,8)

 if (ifig .eq. 0) then
   draw(lg1)
 end if

; dr_txt(wks,"date",xf+wi/2,yf+he+10,22,6,0.)

end do

;dr_txt(wks,"(a)",xf(0)-2,yf(0)-4,25,8,0.)
;dr_txt(wks,"(b)",xf(1)-2,yf(1)-4,25,8,0.)


 frame(wks)

end

