load "$UTIL/ncl/header.ncl"

xspace = 170  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 1.0

begin

;--------------------------------------- parameters
 cexp    = (/"CC"/)
 expname = (/"uj","merra","era-int","sonde_hr"/)
 ppid    = (/"j" /)
 varname = "u"

 y0  = 1953
 y0r = 1979
 y1 = 1987
 y2 = 2006
 lat0 = 5.
 lat0r = lat0 + 0.9

 wks = cr_wks("fig/us_z-"+expname(0),600,0)
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- read data
 ne = dimsizes(expname)
 n0 = (y1-y0)*12
 n9 = (y2-y0+1)*12
 n0r = (y1-y0r)*12
 n9r = (y2-y0r+1)*12

 nt = n9 - n0

   fnam1 = "$DATD/L60CGW-t/xm.uan"+expname(0)+"a.pj.nc"
   fnam1m = "$DATD/clim/MERRA/"+expname(1)+".u_ypt.1979-2006.nc"
   fnam1e = "$DATD/clim/ERA-inter/"+expname(2)+".u_ypt.1979-2006.nc"
   fnam1o = "/data18/kyh/qbo_obs/qbo_"+expname(3)+".nc"
   f1 = addfile(fnam1,"r")
   f1m = addfile(fnam1m,"r")
   f1e = addfile(fnam1e,"r")
   f1o = addfile(fnam1o,"r")
   zu = f1->hybrid_ht / 1.e3
   zum = tofloat(18.6 + 6.6*log( 70./ f1m->levels ))
   zue = tofloat(18.6 + 6.6*log( 70./ f1e->lv_ISBL1 ))
   zuo = 18.6 + 6.6*log( 70./ f1o->p )

   nz = dimsizes(zu)
   nzm = dimsizes(zum)
   nze = dimsizes(zue)
   nzo = dimsizes(zuo)

   u0 = transpose( dim_avg( f1->u(n0:n9-1,:,{-lat0:lat0}) ) )
   u0m = transpose( dim_avg( f1m->u(n0r:n9r-1,:,{-lat0r:lat0r}) ) )
   u0e = transpose( dim_avg( f1e->U_GDS0_ISBL(n0r:n9r-1,:,{-lat0r:lat0r}) ) )
   u0o = f1o->u(:,:nt-1)

 us = u0*0.
 usm = u0m*0.
 use = u0e*0.
 uso = u0o*0.
 do k=1, nz-2
   us(k,:) = (u0(k+1,:) - u0(k-1,:))/(zu(k+1) - zu(k-1))
 end do
 do k=1, nzm-2
   usm(k,:) = (u0m(k+1,:) - u0m(k-1,:))/(zum(k+1) - zum(k-1))
 end do
 do k=1, nze-2
   use(k,:) = (u0e(k+1,:) - u0e(k-1,:))/(zue(k+1) - zue(k-1))
 end do
 do k=1, nzo-2
   uso(k,:) = (u0o(k+1,:) - u0o(k-1,:))/(zuo(k+1) - zuo(k-1))
 end do

 uss = dim_stddev(us)
 ussm = dim_stddev(usm)
 usse = dim_stddev(use)
 usso = dim_stddev(uso)
 usso@_FillValue = 0.

;--------------------------------------- process 0

;--------------------------------------- cn level
 ldsh = (/0,16,2,0/)
 ltck = (/1.8,0.8,0.8,0.8/)

;--------------------------------------- contour
 wi = 35
 he = 35
 xf = 15
 yf = 10

 ystr = "z (km)"
 xstr ="Std. Dev. of :F18:s:F22:U/:F18:s:F22:z (m s:S:-1:N: / km)"
 txtstr = (/"model","MERRA","ERA-I","sonde"/)

do ie=0, ne-1

 if (ie .eq. 0) then
   plot = xylineplot(wks,uss,zu,"D",ldsh(ie),"T",ltck(ie),"C",1)
 end if
 if (ie .eq. 1) then
   plot = xylineplot(wks,ussm,zum,"D",ldsh(ie),"T",ltck(ie),"C",1)
 end if
 if (ie .eq. 2) then
   plot = xylineplot(wks,usse,zue,"D",ldsh(ie),"T",ltck(ie),"C",1)
 end if
 if (ie .eq. 3) then
   plot = xylineplot(wks,usso,zuo,"D",ldsh(ie),"T",ltck(ie),"C",1)
 end if
 xyaddmarker(plot,16,"S",1.0,"T",1.0,"C",1)
 xyaddline(plot,0.,zu,"D",0,"T",0.6,"C",colind)

 set_size(plot,xf,yf,wi,he)
 axis_range(plot,0,7,17,36,False,False)
 axis_str(plot,"",xstr,ystr)

 draw(plot)

end do

 lg = lgplot(wks,txtstr,"D",ldsh,"T",ltck,"C",1,2)
 set_size_lg(lg,xf,yf,wi*0.25,14)
 draw(lg)

 frame(wks)

end

