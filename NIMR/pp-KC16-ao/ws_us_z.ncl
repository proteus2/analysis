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

 wks = cr_wks("fig/ws_us_z",600,0)
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- read data
 ne = dimsizes(expname)*2
 n0 = (y1-y0)*12
 n9 = (y2-y0+1)*12
 n0f = n0 + ydif*12
 n9f = n9 + ydif*12

 nt = n9 - n0

   fnam1 = "$DATD/AOL60CTL/xm.uan"+expname(0)+"a.pj.nc"
   fnam2 = "$DATD/AOL60CGW/xm.uan"+expname(1)+"a.pj.nc"
   f1 = addfile(fnam1,"r")
   f2 = addfile(fnam2,"r")
   zu = f1->hybrid_ht / 1.e3

   nz = dimsizes(zu)

   u0 = transpose( dim_avg( f1->u(n0:n9-1,:,{-lat0:lat0}) ) )
   u1 = transpose( dim_avg( f2->u(n0:n9-1,:,{-lat0:lat0}) ) )
   u2 = transpose( dim_avg( f1->u(n0f:n9f-1,:,{-lat0:lat0}) ) )
   u3 = transpose( dim_avg( f2->u(n0f:n9f-1,:,{-lat0:lat0}) ) )

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

 us0 = u0*0.
 us1 = u1*0.
 us2 = u2*0.
 us3 = u3*0.
 do k=1, nz-2
   us0(k,:) = (u0(k+1,:) - u0(k-1,:))/(zu(k+1) - zu(k-1))
   us1(k,:) = (u1(k+1,:) - u1(k-1,:))/(zu(k+1) - zu(k-1))
   us2(k,:) = (u2(k+1,:) - u2(k-1,:))/(zu(k+1) - zu(k-1))
   us3(k,:) = (u3(k+1,:) - u3(k-1,:))/(zu(k+1) - zu(k-1))
 end do

 uss0 = dim_stddev(us0)
 uss1 = dim_stddev(us1)
 uss2 = dim_stddev(us2)
 uss3 = dim_stddev(us3)
 uss = (/uss0,uss1,uss2,uss3/)

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

 plot = xylineplot(wks,uss(ie,:),zu,"D",ldsh(ie),"T",1.2,"C",lcol(ie))
; xyaddmarker(plot,16,"S",1.0,"T",1.0,"C",lcol(ie))
 xyaddline(plot,0.,zu,"D",0,"T",0.6,"C",colind)

 set_size(plot,xf(1),yf,wi,he)
 axis_range(plot,0.0,6.0,17,33,False,False)
 axis_str(plot,"",xstr(1),ystr)

 draw(plot)

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
 set_size(plotx,xf(1)+wi+4,yf,1,he)
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

 lg = lgplot(wks,txtstr,"D",ldsh,"T",1.2,"C",lcol,2)
 set_size_lg(lg,xf(1),yf,wi*0.25,14)
 draw(lg)

 dr_txt(wks,"(a)",xf(0)-7,yf(0)-2,25,8,0.)
 dr_txt(wks,"(b)",xf(1)-7,yf(0)-2,25,8,0.)

 frame(wks)

end

