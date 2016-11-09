load "$UTIL/ncl/header.ncl"

xspace = 83   ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 1.2

begin

 wks = cr_wks("fig1",600,0)
 colind = NhlNewColor(wks,0.6,0.6,0.6)
 colind2 = NhlNewColor(wks,0.8,0.8,0.8)

;--------------------------------------- parameters
 nf = 2
 nl = 2
;--------------------------------------- read data

 ; read variables to plot
 f1x = addfile("../res/wind/mfa/mfa-klf_p1_f1.nc","r")

 f1k = addfile("../res/temp/mfa/mfa-klf_p1_f1.nc","r")


 z  = f1k->height({20000:}) / 1.e3
 nz = dimsizes(z)
;;;
 var_dat = new((/nf,nl,nz/),float,1.e32)

 var_dat(0,0,:) = f1k->mfx({20000:}) * 1.e3
 var_dat(0,1,:) = f1x->mfx({20000:}) * 1.e3

 var_dat(1,0,:) = f1k->mfx_e({20000:}) * 1.e3
 var_dat(1,1,:) = f1k->mfx_w({20000:}) * 1.e3
;;;

;--------------------------------------- process 1
 var = var_dat

err = var(0,0,:)-var(0,1,:)
print(sprintf("%7.3f",err(:))+"   "+sprintf("%6.2f",err(:)/var(0,1,:)*100.))
print(avg(err)+"   "+avg(err)/avg(var(0,1,:))*100.)
;--------------------------------------- line attribute
 ldsh = (/0,0,0,1/)
 ltck = (/2.0,0.6,1.2,1.2/)

;--------------------------------------- plot

 wi = 43
 he = 43
 xf = 12 + (/0,0/)*(wi+7)
 yf = 12 + (/0,1/)*(he+12)

 labels = (/"F:B:T:N:","F:B:W:N:","F:B:T:N: (east)","F:B:T:N: (west)"/)
 ystr = (/"z [km]","z [km]"/)

do ifig=0, nf-1

 isub1 = ifig*nl
 isub2 = (ifig+1)*nl-1

 plot = xylineplot(wks,var(ifig,:,:),z,"D",ldsh(isub1:isub2), \
                   "T",ltck(isub1:isub2),"C",1)

 if (ifig .eq. 1) then
 xyaddline(plot,0,z,"D",0,"T",0.8,"C",colind2)
 end if

 set_size(plot,xf(ifig),yf(ifig),wi,he)
 axis_range(plot,-999,999,-999,45,False,False)
 axis_str(plot,"","F [mPa]",ystr(ifig))
; tick_exp(plot,"XB",xval,xlab,time)

 draw(plot)


 lg1 = lgplot(wks,labels(isub1:isub2),"D",ldsh(isub1:isub2), \
                  "T",ltck(isub1:isub2),"C",1,1)

 set_size_lg(lg1,xf(ifig)+wi-10,yf(ifig),10,8)

 draw(lg1)


end do

dr_txt(wks,"(a)",xf(0)-8,yf(0)-2,25,8,0.)
dr_txt(wks,"(b)",xf(1)-8,yf(1)-2,25,8,0.)


 frame(wks)

end

