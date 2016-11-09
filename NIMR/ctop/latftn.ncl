load "$UTIL/ncl/header.ncl"

xspace = 170  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 1.2

begin

;--------------------------------------- parameters
 wks = cr_wks("fig/latftn",600,0)
; gsn_define_colormap(wks,"StepSeq25")
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- base data
 lat = fspan(-90,90,181)
 ftn = 1. + 0.20*(tanh((lat-50.)/10.) - tanh((lat+50.)/10.))

;--------------------------------------- plot

 xf = 20
 yf = 10
 wi = 50
 he = 25
 ystr = (/"E", "W"/)+"-MF (mPa)"
 llbl = (/"CGW","CGW+BGW","BGW (CTL)"/)

 plot = xylineplot(wks,lat,ftn,"D",0,"T",1.0,"C",1)

 set_size(plot,xf,yf,wi,he)
 axis_range(plot,-90,90,0,1.1,False,False)
 axis_str(plot,"Reduction factor","lat (deg)","")

 draw(plot)

 frame(wks)

end



