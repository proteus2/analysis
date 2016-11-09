load "$UTIL/ncl/header.ncl"

xspace = 170  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 1.0

begin

;--------------------------------------- parameters
 cexp    = "CC"
 expname = "tc"
 ppid    = "e"
 varname = (/"dT_conv","dT_lsrain"/)

 y0 = 1951
 y1 = 1955
 y2 = 1994

 ydif = 100
 lat0 = 15.

 zb = 0
 zt = 21

 wks = cr_wks("fig/heat_z-a",600,0)
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- read data
 n0 = (y1-y0)*12
 n9 = (y2-y0+1)*12
 n0f = n0 + ydif*12
 n9f = n9 + ydif*12
 n0m = n0 + (ydif/2)*12
 n9m = n9 + (ydif/2)*12
 zbr = (zb-1.5)*1.e3
 ztr = (zt+1.5)*1.e3

 nt = n9 - n0

 fnam = "$DATD/AOL60CGW/xm.uan"+expname+"a.p"+ppid+".nc"

 f1 = addfile(fnam,"r")
 za = f1->hybrid_ht({zbr:ztr}) / 1.e3

 dt0 = dim_avg_n( dim_avg( f1->$varname(0)$(n0:n9-1,{zbr:ztr},{-lat0:lat0}) ), 0 )
 dt0f = dim_avg_n( dim_avg( f1->$varname(0)$(n0f:n9f-1,{zbr:ztr},{-lat0:lat0}) ), 0 )
 dt0m = dim_avg_n( dim_avg( f1->$varname(0)$(n0m:n9m-1,{zbr:ztr},{-lat0:lat0}) ), 0 )
 dt1 = dim_avg_n( dim_avg( f1->$varname(1)$(n0:n9-1,{zbr:ztr},{-lat0:lat0}) ), 0 )
 dt1f = dim_avg_n( dim_avg( f1->$varname(1)$(n0f:n9f-1,{zbr:ztr},{-lat0:lat0}) ), 0 )
 dt1m = dim_avg_n( dim_avg( f1->$varname(1)$(n0m:n9m-1,{zbr:ztr},{-lat0:lat0}) ), 0 )

;--------------------------------------- process 0
; dta = (/dt0,dt0f,dt1,dt1f/)*(3.*24.)
; dta = where(dta .eq. 0., 1.e-20,dta)

 dta = (/dt0+dt1,dt0m+dt1m,dt0f+dt1f/)*(3.*24.)

;--------------------------------------- cn level
 lcol = (/1,1,1,2,2,2,3,3,3/)
 ldsh = (/0,2,16,0,2,16,0,2,16/)

;--------------------------------------- contour
 wi = 35
 he = 35
 xf = 15
 yf = 10

 ystr = "height [km]"
 xstr = "moist heating [K day:S:-1:N:]"
; txtstr = (/y1+"-"+y2+", Conv.",(y1+ydif)+"-"+(y2+ydif)+", Conv.", \
;            y1+"-"+y2+", LS rain",(y1+ydif)+"-"+(y2+ydif)+", LS rain"/)
 txtstr = (/y1+"-"+y2,(y1+ydif/2)+"-"+(y2+ydif/2),(y1+ydif)+"-"+(y2+ydif)/)

;do il=0, dimsizes(varname)*2-1

 plot = xylineplot(wks,dta,za,"D",ldsh,"T",1.2,"C",lcol)
; xyaddmarker(plot,16,"S",1.0,"T",1.0,"C",lcol(il))
 xyaddline(plot,0.,za,"D",0,"T",0.6,"C",colind)

 set_size(plot,xf,yf,wi,he)
 axis_range(plot,-999,999,zb,zt,False,False)
 axis_str(plot,"",xstr,ystr)

 draw(plot)

 ps = (/500,300,200,150,100,70,50,40,30,20,10,7,5,3,2,1/)
 zs = 18.6 + 6.6*log(70./ps)
 plotx = xylineplot(wks,(/-1,0,1/),(/-1,0,1/),"D",0,"T",0.001,"C",-1)
 set_size(plotx,xf+wi+2,yf,1,he)
 axis_range(plotx,-999,999,zb,zt,False,False)
 axis_str(plotx,"","","p [hPa]")
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

;end do

 lg = lgplot(wks,txtstr,"D",ldsh,"T",1.2,"C",lcol,2)
 set_size_lg(lg,xf+wi*0.5,yf,wi*0.3,9)
 draw(lg)

 frame(wks)

end

