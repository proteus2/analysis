load "/data3/kyh/util/ncl/header.ncl"

xspace = 170  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "/data3/kyh/util/ncl/figsize.ncl"
load "/data3/kyh/util/ncl/attribute.ncl"

axisfontratio = 1.2

begin

 wks = cr_wks("basic",300,0)
; gsn_define_colormap(wks,"posneg_2")
 colind = NhlNewColor(wks,0.5,0.5,0.5)

;--------------------------------------- parameters
 nf = 3
 nl = 3
;--------------------------------------- read data

 ; read variables to plot
 f1x = addfile("/data4/atmosdata/T_SPECTRUM/P1/n_profile_ave_p1.nc","r")
 f2x = addfile("/data4/atmosdata/T_SPECTRUM/P2/n_profile_ave_p2.nc","r")
 f3x = addfile("/data4/atmosdata/T_SPECTRUM/P3/n_profile_ave_p3.nc","r")

 f1k = addfile("/data4/atmosdata/T_SPECTRUM/P1/u_profile_ave_p1.nc","r")
 f2k = addfile("/data4/atmosdata/T_SPECTRUM/P2/u_profile_ave_p2.nc","r")
 f3k = addfile("/data4/atmosdata/T_SPECTRUM/P3/u_profile_ave_p3.nc","r")

 f1c = addfile("/data4/atmosdata/T_SPECTRUM/P1/v_profile_ave_p1.nc","r")
 f2c = addfile("/data4/atmosdata/T_SPECTRUM/P2/v_profile_ave_p2.nc","r")
 f3c = addfile("/data4/atmosdata/T_SPECTRUM/P3/v_profile_ave_p3.nc","r")


 z  = f1k->height({20000:}) / 1.e3
 nz = dimsizes(z)
;;;
 var_dat = new((/nf,nl,nz/),float,1.e32)
 var_dat(0,0,:) = f1x->n_ave({20000:}) * 1.e2
 var_dat(1,0,:) = f2x->n_ave({20000:}) * 1.e2
 var_dat(2,0,:) = f3x->n_ave({20000:}) * 1.e2

 var_dat(0,1,:) = f1k->u_ave({20000:})
 var_dat(1,1,:) = f2k->u_ave({20000:})
 var_dat(2,1,:) = f3k->u_ave({20000:})

 var_dat(0,2,:) = f1c->v_ave({20000:})
 var_dat(1,2,:) = f2c->v_ave({20000:})
 var_dat(2,2,:) = f3c->v_ave({20000:})
;;;

;--------------------------------------- process 1
 var = var_dat

;--------------------------------------- line attribute
 lcol = (/-1,1,1/)
 ltck = (/1.2,1.2,1.2/)
 ldsh = (/0,0,14/)

;--------------------------------------- plot

 wi = 40
 he = 35
 xf = 10 + (/0,1,2,0,1,2/)*(wi+10)
 yf = 10 + (/0,0,0,1,1,1/)*(he+10)

 ystr = (/"z (km)","","","z (km)","",""/)

do ifig=0, nf-1

 plot = xylineplot(wks,var(ifig,:,:),z,"D",ldsh,"T",ltck,"C",lcol)

; xyaddline(plot,0,z,"D",0,"T",0.6,"C",colind)

 set_size(plot,xf(ifig),yf(ifig),wi,he)
 axis_range(plot,-999,999,-999,45,False,False)
 axis_str(plot,"","U, V (m s:S:-1:N:)",ystr(ifig))
; tick_exp(plot,"XB",xval,xlab,time)
 setvalues plot
  "tmXTOn" : False
 end setvalues

 plot2 = xylineplot(wks,var(ifig,0,:),z,"D",0,"T",1.2,"C",2)

 set_size(plot2,xf(ifig),yf(ifig),wi,he)
 axis_range(plot2,0,3,-999,45,False,False)
 axis_str(plot2,"","N (10:S:-2:N: s:S:-1:N:)","")
 setvalues plot2
  "tmXBOn" : False
  "tmXTLabelsOn" : True
  "tiXAxisSide" : "Top"
 end setvalues


 draw(plot)
 draw(plot2)

; dr_txt(wks,"date",xf+wi/2,yf+he+10,22,6,0.)

end do


 frame(wks)

end

