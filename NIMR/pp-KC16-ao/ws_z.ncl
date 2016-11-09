load "$UTIL/ncl/header.ncl"

xspace = 170  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 1.0

begin

;--------------------------------------- parameters
 cexp    = (/"CTL","CC"/)
 expname = (/"tk","tc"/)
 ppid    = (/"j" /)
 varname = "u"

 y0 = 1951
 y1 = 1960
 y2 = 1999
 ydif = 98
 lat0 = 10.

 wks = cr_wks("fig/ws_z",600,0)
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- read data
 ne = dimsizes(expname)*2
 n0 = (y1-y0)*12
 n9 = (y2-y0+1)*12
 n0f = n0 + ydif*12
 n9f = n9 + ydif*12

 nt = n9 - n0

   fnam1 = "$DATD/AOL60CTL/tem/uan"+expname(0)+".tem-wc_yzt.nc"
   fnam2 = "$DATD/AOL60CGW/tem/uan"+expname(1)+".tem-wc_yzt.nc"
   f1 = addfile(fnam1,"r")
   f2 = addfile(fnam2,"r")
   za = f1->z / 1.e3

   ws0 = transpose( dim_avg( f1->w_res(n0:n9-1,:,{-lat0:lat0}) ) )
   ws1 = transpose( dim_avg( f2->w_res(n0:n9-1,:,{-lat0:lat0}) ) )
   ws2 = transpose( dim_avg( f1->w_res(n0f:n9f-1,:,{-lat0:lat0}) ) )
   ws3 = transpose( dim_avg( f2->w_res(n0f:n9f-1,:,{-lat0:lat0}) ) )

 ws0a = dim_avg(ws0)
 ws1a = dim_avg(ws1)
 ws2a = dim_avg(ws2)
 ws3a = dim_avg(ws3)

;--------------------------------------- process 0
 ws0a = ws0a*1.e3
 ws1a = ws1a*1.e3
 ws2a = ws2a*1.e3
 ws3a = ws3a*1.e3
 wsa = (/ws0a,ws1a,ws2a,ws3a/)

;--------------------------------------- cn level
 lcol = (/1,2,1,2/)
 ldsh = (/14,14,0,0/)

;--------------------------------------- contour
 wi = 35
 he = 35
 xf = (/15,62/)
 yf = 10

 ystr = "height (km)"
 xstr = (/"W* (mm s:S:-1:N:)","Std. Dev. of :F18:s:F22:U/:F18:s:F22:z (m s:S:-1:N: / km)"/)
 txtstr = (/y1+"-"+y2+", CTL",y1+"-"+y2+", CC",(y1+ydif)+"-"+(y2+ydif)+", CTL",(y1+ydif)+"-"+(y2+ydif)+", CC"/)

do ie=0, ne-1

 plot = xylineplot(wks,wsa(ie,:),za,"D",ldsh(ie),"T",1.2,"C",lcol(ie))
; xyaddmarker(plot,16,"S",1.0,"T",1.0,"C",lcol(ie))
 xyaddline(plot,0.,za,"D",0,"T",0.6,"C",colind)

 set_size(plot,xf(0),yf,wi,he)
 axis_range(plot,0,0.7,17,33,False,False)
 axis_str(plot,"",xstr(0),ystr)

 draw(plot)

 ps = (/100,70,50,40,30,20,10,7,5,3,2,1/)
 zs = 18.6 + 6.6*log(70./ps)
 plotx = xylineplot(wks,(/-1,0,1/),(/-1,0,1/),"D",0,"T",0.001,"C",lcol(ie))
 set_size(plotx,xf(0)+wi+4,yf,1,he)
 axis_range(plotx,-999,999,17,33,False,False)
 axis_str(plotx,"","","p (hPa)")
 setvalues plotx
  "tmXBOn" : False
  "tmXTOn" : False
  "tmYLOn" : False
  "tmYROn" : True
  "tmXBBorderOn" : False
  "tmXTBorderOn" : False
  "tmYLBorderOn" : False
  "tiYAxisSide" : "Right"
  "tiYAxisAngleF" : 270
  "tmYRLabelsOn" : True
  "tmYLMode" : "Explicit"
  "tmYLValues" : zs
  "tmYLLabels" : ps
 end setvalues
 draw(plotx)

end do

 lg = lgplot(wks,txtstr,"D",ldsh,"T",1.2,"C",lcol,2)
 set_size_lg(lg,xf(0),yf,wi*0.25,14)
 draw(lg)

 frame(wks)

end

