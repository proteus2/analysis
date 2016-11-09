load "$UTIL/ncl/header.ncl"

xspace = 170  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 1.2

begin

;--------------------------------------- parameters
 v_s = 30.
 t_b = 1.5
 c_w = 30.
 nc = 65

;--------------------------------------- workstation
 wks = cr_wks("fig/gaussian",600,0)
; gsn_define_colormap(wks,"prcp_3")
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- read data

;--------------------------------------- process 0
 c = fspan(v_s-80.,v_s+80.,nc)

 tau = t_b*exp(-(((c-v_s)/c_w)^2))
 tau(nc/2) = 0.

;--------------------------------------- cn level

;--------------------------------------- contour
 wi = 60
 he = 30
 xf = 15
 yf = 10

 plot = xylineplot(wks,c,tau,"D",0,"T",1.0,"C",-1)
 xyaddmarker(plot,16,"S",1.2,"T",0.1,"C",1)

 set_size(plot,xf,yf,wi,he)
 axis_range(plot,min(c),max(c),0,999,False,False)
 axis_str(plot,"","c (m s:S:-1:N:)","M (mPa)")

 draw(plot)

 frame(wks)

end

