load "/data3/kyh/util/ncl/header.ncl"

xspace = 170  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "/data3/kyh/util/ncl/figsize.ncl"
load "/data3/kyh/util/ncl/attribute.ncl"

axisfontratio = 1.2

begin

 wks = cr_wks("mf-w_z",300,0)
; gsn_define_colormap(wks,"posneg_2")
 colind = NhlNewColor(wks,0.5,0.5,0.5)

;--------------------------------------- parameters
 nf = 6
 nl = 9
;--------------------------------------- read data

 ; read variables to plot
 f1x = addfile("res/wind/mean-mf/mean-mfx_p1.nc","r")
 f2x = addfile("res/wind/mean-mf/mean-mfx_p2.nc","r")
 f3x = addfile("res/wind/mean-mf/mean-mfx_p3.nc","r")
 f1y = addfile("res/wind/mean-mf/mean-mfy_p1.nc","r")
 f2y = addfile("res/wind/mean-mf/mean-mfy_p2.nc","r")
 f3y = addfile("res/wind/mean-mf/mean-mfy_p3.nc","r")

 f1k = addfile("res/wind/mfa/mfa-klf_p1_f1.nc","r")
 f2k = addfile("res/wind/mfa/mfa-klf_p2_f1.nc","r")
 f3k = addfile("res/wind/mfa/mfa-klf_p3_f1.nc","r")

 f1c = addfile("res/wind/mfa/mfa-c_phi_p1_f1.nc","r")
 f2c = addfile("res/wind/mfa/mfa-c_phi_p2_f1.nc","r")
 f3c = addfile("res/wind/mfa/mfa-c_phi_p3_f1.nc","r")


 z  = f1k->height({20000:}) / 1.e3
 nz = dimsizes(z)
;;;
 var_dat = new((/nf,nl,nz/),float,1.e32)
 var_dat(0,0,:) = f1x->mfx({20000:}) * 1.e3
 var_dat(1,0,:) = f2x->mfx({20000:}) * 1.e3
 var_dat(2,0,:) = f3x->mfx({20000:}) * 1.e3
 var_dat(3,0,:) = f1y->mfy({20000:}) * 1.e3
 var_dat(4,0,:) = f2y->mfy({20000:}) * 1.e3
 var_dat(5,0,:) = f3y->mfy({20000:}) * 1.e3

 var_dat(0,1,:) = f1k->mfx({20000:}) * 1.e3
 var_dat(1,1,:) = f2k->mfx({20000:}) * 1.e3
 var_dat(2,1,:) = f3k->mfx({20000:}) * 1.e3
 var_dat(3,1,:) = f1k->mfy({20000:}) * 1.e3
 var_dat(4,1,:) = f2k->mfy({20000:}) * 1.e3
 var_dat(5,1,:) = f3k->mfy({20000:}) * 1.e3

 var_dat(0,2,:) = f1c->mfx({20000:}) * 1.e3
 var_dat(1,2,:) = f2c->mfx({20000:}) * 1.e3
 var_dat(2,2,:) = f3c->mfx({20000:}) * 1.e3
 var_dat(3,2,:) = f1c->mfy({20000:}) * 1.e3
 var_dat(4,2,:) = f2c->mfy({20000:}) * 1.e3
 var_dat(5,2,:) = f3c->mfy({20000:}) * 1.e3

 var_dat(0,3,:) = f1x->mfx_p({20000:}) * 1.e3
 var_dat(1,3,:) = f2x->mfx_p({20000:}) * 1.e3
 var_dat(2,3,:) = f3x->mfx_p({20000:}) * 1.e3
 var_dat(3,3,:) = f1y->mfy_p({20000:}) * 1.e3
 var_dat(4,3,:) = f2y->mfy_p({20000:}) * 1.e3
 var_dat(5,3,:) = f3y->mfy_p({20000:}) * 1.e3

 var_dat(0,4,:) = f1k->mfx_p({20000:}) * 1.e3
 var_dat(1,4,:) = f2k->mfx_p({20000:}) * 1.e3
 var_dat(2,4,:) = f3k->mfx_p({20000:}) * 1.e3
 var_dat(3,4,:) = f1k->mfy_p({20000:}) * 1.e3
 var_dat(4,4,:) = f2k->mfy_p({20000:}) * 1.e3
 var_dat(5,4,:) = f3k->mfy_p({20000:}) * 1.e3

 var_dat(0,5,:) = f1c->mfx_p({20000:}) * 1.e3
 var_dat(1,5,:) = f2c->mfx_p({20000:}) * 1.e3
 var_dat(2,5,:) = f3c->mfx_p({20000:}) * 1.e3
 var_dat(3,5,:) = f1c->mfy_p({20000:}) * 1.e3
 var_dat(4,5,:) = f2c->mfy_p({20000:}) * 1.e3
 var_dat(5,5,:) = f3c->mfy_p({20000:}) * 1.e3

 var_dat(0,6,:) = f1x->mfx_n({20000:}) * 1.e3
 var_dat(1,6,:) = f2x->mfx_n({20000:}) * 1.e3
 var_dat(2,6,:) = f3x->mfx_n({20000:}) * 1.e3
 var_dat(3,6,:) = f1y->mfy_n({20000:}) * 1.e3
 var_dat(4,6,:) = f2y->mfy_n({20000:}) * 1.e3
 var_dat(5,6,:) = f3y->mfy_n({20000:}) * 1.e3

 var_dat(0,7,:) = f1k->mfx_m({20000:}) * 1.e3
 var_dat(1,7,:) = f2k->mfx_m({20000:}) * 1.e3
 var_dat(2,7,:) = f3k->mfx_m({20000:}) * 1.e3
 var_dat(3,7,:) = f1k->mfy_m({20000:}) * 1.e3
 var_dat(4,7,:) = f2k->mfy_m({20000:}) * 1.e3
 var_dat(5,7,:) = f3k->mfy_m({20000:}) * 1.e3

 var_dat(0,8,:) = f1c->mfx_n({20000:}) * 1.e3
 var_dat(1,8,:) = f2c->mfx_n({20000:}) * 1.e3
 var_dat(2,8,:) = f3c->mfx_n({20000:}) * 1.e3
 var_dat(3,8,:) = f1c->mfy_n({20000:}) * 1.e3
 var_dat(4,8,:) = f2c->mfy_n({20000:}) * 1.e3
 var_dat(5,8,:) = f3c->mfy_n({20000:}) * 1.e3
;;;

;--------------------------------------- process 1
 var = var_dat

;--------------------------------------- line attribute
 lcol = (/1,2,-1,1,2,-1,1,2,-1/)
 ltck = (/1.5,1.5,1.5,1.2,1.2,1.2,1.2,1.2,1.2/)
 ldsh = (/0,0,0,14,14,14,14,14,14/)

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

