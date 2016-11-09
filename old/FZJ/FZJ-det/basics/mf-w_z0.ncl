load "/data3/kyh/util/ncl/header.ncl"

xspace = 170  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "/data3/kyh/util/ncl/figsize.ncl"
load "/data3/kyh/util/ncl/attribute.ncl"

axisfontratio = 1.2

begin

 wks = cr_wks("mf-w_z0",300,0)
; gsn_define_colormap(wks,"posneg_2")
 colind = NhlNewColor(wks,0.5,0.5,0.5)

;--------------------------------------- parameters
 nf = 6
 nl = 6
;--------------------------------------- read data

 ; read variables to plot
 f1x = addfile("../../res/wind/mean-mf/mean-mfx_p1.nc","r")
 f2x = addfile("../../res/wind/mean-mf/mean-mfx_p2.nc","r")
 f3x = addfile("../../res/wind/mean-mf/mean-mfx_p3.nc","r")
 f1y = addfile("../../res/wind/mean-mf/mean-mfy_p1.nc","r")
 f2y = addfile("../../res/wind/mean-mf/mean-mfy_p2.nc","r")
 f3y = addfile("../../res/wind/mean-mf/mean-mfy_p3.nc","r")

 f1a = addfile("../res/wind/mean-mf/mean-mfx_p1.nc","r")
 f2a = addfile("../res/wind/mean-mf/mean-mfx_p2.nc","r")
 f3a = addfile("../res/wind/mean-mf/mean-mfx_p3.nc","r")
 f1b = addfile("../res/wind/mean-mf/mean-mfy_p1.nc","r")
 f2b = addfile("../res/wind/mean-mf/mean-mfy_p2.nc","r")
 f3b = addfile("../res/wind/mean-mf/mean-mfy_p3.nc","r")

 z  = f1x->height_m({20000:}) / 1.e3
 nz = dimsizes(z)
;;;
 var_dat = new((/nf,nl,nz/),float,1.e32)
 var_dat(0,0,:) = f1x->mfx({20000:}) * 1.e3
 var_dat(1,0,:) = f2x->mfx({20000:}) * 1.e3
 var_dat(2,0,:) = f3x->mfx({20000:}) * 1.e3
 var_dat(3,0,:) = f1y->mfy({20000:}) * 1.e3
 var_dat(4,0,:) = f2y->mfy({20000:}) * 1.e3
 var_dat(5,0,:) = f3y->mfy({20000:}) * 1.e3

 var_dat(0,1,:) = f1a->mfx({20000:}) * 1.e3
 var_dat(1,1,:) = f2a->mfx({20000:}) * 1.e3
 var_dat(2,1,:) = f3a->mfx({20000:}) * 1.e3
 var_dat(3,1,:) = f1b->mfy({20000:}) * 1.e3
 var_dat(4,1,:) = f2b->mfy({20000:}) * 1.e3
 var_dat(5,1,:) = f3b->mfy({20000:}) * 1.e3

 var_dat(0,2,:) = f1x->mfx_p({20000:}) * 1.e3
 var_dat(1,2,:) = f2x->mfx_p({20000:}) * 1.e3
 var_dat(2,2,:) = f3x->mfx_p({20000:}) * 1.e3
 var_dat(3,2,:) = f1y->mfy_p({20000:}) * 1.e3
 var_dat(4,2,:) = f2y->mfy_p({20000:}) * 1.e3
 var_dat(5,2,:) = f3y->mfy_p({20000:}) * 1.e3

 var_dat(0,3,:) = f1a->mfx_p({20000:}) * 1.e3
 var_dat(1,3,:) = f2a->mfx_p({20000:}) * 1.e3
 var_dat(2,3,:) = f3a->mfx_p({20000:}) * 1.e3
 var_dat(3,3,:) = f1b->mfy_p({20000:}) * 1.e3
 var_dat(4,3,:) = f2b->mfy_p({20000:}) * 1.e3
 var_dat(5,3,:) = f3b->mfy_p({20000:}) * 1.e3

 var_dat(0,4,:) = f1x->mfx_n({20000:}) * 1.e3
 var_dat(1,4,:) = f2x->mfx_n({20000:}) * 1.e3
 var_dat(2,4,:) = f3x->mfx_n({20000:}) * 1.e3
 var_dat(3,4,:) = f1y->mfy_n({20000:}) * 1.e3
 var_dat(4,4,:) = f2y->mfy_n({20000:}) * 1.e3
 var_dat(5,4,:) = f3y->mfy_n({20000:}) * 1.e3

 var_dat(0,5,:) = f1a->mfx_n({20000:}) * 1.e3
 var_dat(1,5,:) = f2a->mfx_n({20000:}) * 1.e3
 var_dat(2,5,:) = f3a->mfx_n({20000:}) * 1.e3
 var_dat(3,5,:) = f1b->mfy_n({20000:}) * 1.e3
 var_dat(4,5,:) = f2b->mfy_n({20000:}) * 1.e3
 var_dat(5,5,:) = f3b->mfy_n({20000:}) * 1.e3
;;;

;--------------------------------------- process 1
 var = var_dat

;--------------------------------------- line attribute
 lcol = (/4,2,4,2,4,2/)
 ltck = (/1.5,1.5,1.2,1.2,1.2,1.2/)
 ldsh = (/0,0,14,14,14,14/)

;--------------------------------------- plot

 wi = 40
 he = 35
 xf = 10 + (/0,1,2,0,1,2/)*(wi+10)
 yf = 10 + (/0,0,0,1,1,1/)*(he+10)

 xstr = (/"x","x","x","y","y","y"/)
 ystr = (/"z (km)","","","z (km)","",""/)

do ifig=0, nf-1

 plot = xylineplot(wks,var(ifig,:,:),z,"D",ldsh,"T",ltck,"C",lcol)

 xyaddline(plot,0,z,"D",0,"T",0.6,"C",colind)

 set_size(plot,xf(ifig),yf(ifig),wi,he)
 axis_range(plot,-999,999,-999,45,False,False)
 axis_str(plot,"","MF"+xstr(ifig)+" (mPa)",ystr(ifig))
; tick_exp(plot,"XB",xval,xlab,time)

 draw(plot)

; dr_txt(wks,"date",xf+wi/2,yf+he+10,22,6,0.)

end do


 frame(wks)

end

