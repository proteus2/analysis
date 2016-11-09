load "$UTIL/ncl/header.ncl"

xspace = 170  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 1.0

begin

;--------------------------------------- parameters
 cexp    = (/"CC"/)
 expname = (/"uj","merra","era-int"/)
 ppid    = (/"j" /)
 varname = "u"

 y0  = 1953
 y0r = 1979
 y1 = 1979
 y2 = 2006
 lat0 = 5.
 lat0r = lat0 + 0.9

 wks = cr_wks("fig/ws_z-"+expname(0),600,0)
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- read data
 ne = dimsizes(expname)
 n0 = (y1-y0)*12
 n9 = (y2-y0+1)*12
 n0r = (y1-y0r)*12
 n9r = (y2-y0r+1)*12

   fnam1 = "$DATD/L60CGW-t/xm.uan"+expname(0)+"a.pj.nc"
   fnam1m = "$DATD/clim/MERRA/"+expname(1)+".u_ypt.1979-2006.nc"
   fnam1e = "$DATD/clim/ERA-inter/"+expname(2)+".u_ypt.1979-2006.nc"
   f1 = addfile(fnam1,"r")
   f1m = addfile(fnam1m,"r")
   f1e = addfile(fnam1e,"r")
   zu = f1->hybrid_ht / 1.e3
   zum = 18.6 + 6.6*log( 70./ f1m->levels )
   zue = 18.6 + 6.6*log( 70./ f1e->lv_ISBL1 )
print(f1->latitude({-lat0:lat0}))
print(f1m->latitude({-lat0r:lat0r}))
print(f1e->g0_lat_2({-lat0r:lat0r}))

   u0 = transpose( dim_avg( f1->u(n0:n9-1,:,{-lat0:lat0}) ) )
   u0m = transpose( dim_avg( f1m->u(n0r:n9r-1,:,{-lat0r:lat0r}) ) )
   u0e = transpose( dim_avg( f1e->U_GDS0_ISBL(n0r:n9r-1,:,{-lat0r:lat0r}) ) )

   fnam2 = "$DATD/L60CGW-t/tem/uan"+expname(0)+".tem-wc_yzt.nc"
   fnam2m = "$DATD/tem/MERRA/"+expname(1)+".tem_ypt.1979-2006.nc"
   fnam2e = "$DATD/tem/ERA-inter/"+expname(2)+".tem_ypt.1979-2006.nc"
   f2 = addfile(fnam2,"r")
   f2m = addfile(fnam2m,"r")
   f2e = addfile(fnam2e,"r")
   za = f2->z / 1.e3
   zam = 18.6 + 6.6*log( 70./ f2m->p )
   zae = 18.6 + 6.6*log( 70./ f2e->p )

   ws = transpose( dim_avg( f2->w_res(n0:n9-1,:,{-lat0:lat0}) ) )
   wsm = transpose( dim_avg( f2m->w_res(n0r:n9r-1,:,{-lat0r:lat0r}) ) )
   wse = transpose( dim_avg( f2e->w_res(n0r:n9r-1,:,{-lat0r:lat0r}) ) )

 wsa = dim_avg(ws)
 wsam = dim_avg(wsm)
 wsae = dim_avg(wse)

;--------------------------------------- process 0
 wsa = wsa*1.e3
 wsam = wsam*1.e3
 wsae = wsae*1.e3

;--------------------------------------- cn level
 ldsh = (/0,16,2,0/)
 ltck = (/1.8,0.8,0.8,0.8/)

;--------------------------------------- contour
 wi = 35
 he = 35
 xf = 15
 yf = 10

 ystr = "z (km)"
 xstr ="W* (mm s:S:-1:N:)"
 txtstr = (/"model","MERRA","ERA-I"/)

do ie=0, ne-1

 if (ie .eq. 0) then
   plot = xylineplot(wks,wsa,za,"D",ldsh(ie),"T",ltck(ie),"C",1)
 end if
 if (ie .eq. 1) then
   plot = xylineplot(wks,wsam,zam,"D",ldsh(ie),"T",ltck(ie),"C",1)
 end if
 if (ie .eq. 2) then
   plot = xylineplot(wks,wsae,zae,"D",ldsh(ie),"T",ltck(ie),"C",1)
 end if
 xyaddmarker(plot,16,"S",1.0,"T",1.0,"C",1)
 xyaddline(plot,0.,za,"D",0,"T",0.6,"C",colind)

 set_size(plot,xf,yf,wi,he)
 axis_range(plot,-0.2,1,17,36,False,False)
 axis_str(plot,"",xstr,ystr)

 draw(plot)

end do

 lg = lgplot(wks,txtstr,"D",ldsh,"T",ltck,"C",1,2)
 set_size_lg(lg,xf+wi*0.58,yf+he-10,wi*0.25,10.5)
 draw(lg)

 frame(wks)

end

