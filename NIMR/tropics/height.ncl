load "$UTIL/ncl/header.ncl"

xspace =  83  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 2.0

begin

;--------------------------------------- parameters
 fnam0 = "$DATD/L60CGW-t/xm.uanuja.pj.nc"
 fnam1 = "$DATD/clim/MERRA/merra.h_ypt.1979-2010.nc"

 lat0 = 5.

 wks = cr_wks("fig/height",600,0)
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- read data
   f0 = addfile(fnam0,"r")
   f1 = addfile(fnam1,"r")

   p0   = dim_avg( f0->p  (:,{14.e3:50.e3},{-lat0:lat0}) )/100.
   r0   = dim_avg( f0->rho(:,{14.e3:50.e3},{-lat0:lat0}) )
   lnp0 = log(p0)*(-1.)
   z0   = f0->$f0->p!1$({14.e3:50.e3}) / 1.e3

   p   = f1->levels({150.:1.})
   lnp = log(p)*(-1.)
   z   = dim_avg( f1->h(:,{150.:1.},{-lat0:lat0}) ) / 1.e3

   tmp = dimsizes(p0)
   nt0 = tmp(0)
   tmp = dimsizes(z)
   nt  = tmp(0)

   hs = 6.6
   zp = (/ \
18.5 + hs*(log(70.  ) + lnp), \
 1.0 + hs*(log(1000.) + lnp)  \
        /)

   rs = p*100./9.8/(hs*1.e3)

;--------------------------------------- process 0
 lnp0s = new((/12,dimsizes(z0)/),float)
 do n=0, 11
   lnp0s(n,:) = dim_avg_n(lnp0(ispan(n,nt0-1,12),:),0)
 end do
 p0s = exp(-lnp0s)

 r0s = new((/12,dimsizes(z0)/),float)
 do n=0, 11
   r0s(n,:) = dim_avg_n(r0(ispan(n,nt0-1,12),:),0)
 end do

 zs = new((/12,dimsizes(p)/),float)
 do n=0, 11
   zs(n,:) = dim_avg_n(z(ispan(n,nt-1,12),:),0)
 end do

;--------------------------------------- cn level
 lcol = (/1,3,4,6,7/)

;--------------------------------------- contour
 wi = 60
 he = 60
 xf = 15
 xf2= 85
 yf = 20

 do n=0, 11

 plot = xylineplot(wks,p,zs(n,:),"D",0,"T",1.0,"C",-1)
 xyaddmarker(plot,1,"S",2.0,"T",1.0,"C",1)

 set_size(plot,xf,yf,wi,he)
 axis_range(plot,0.7,200.,10.,50.,True,False)
 axis_str(plot,"","","")

 draw(plot)

 plot0 = xylineplot(wks,p0s(n,:),z0,"D",0,"T",1.0,"C",-1)
 xyaddmarker(plot0,1,"S",2.0,"T",1.0,"C",2)

 set_size(plot0,xf,yf,wi,he)
 axis_range(plot0,0.7,200.,10.,50.,True,False)
 axis_str(plot0,"","","")

 draw(plot0)

 plotb0 = xylineplot(wks,r0s(n,:),z0,"D",0,"T",1.0,"C",-1)
 xyaddmarker(plotb0,1,"S",2.0,"T",1.0,"C",2)

 set_size(plotb0,xf2,yf,wi,he)
 axis_range(plotb0,0.8e-3,0.4,10.,50.,True,False)
 axis_str(plotb0,"","","")

 draw(plotb0)

 end do

 plot1 = xylineplot(wks,p,zp,"D",0,"T",1.0,"C",lcol)

 set_size(plot1,xf,yf,wi,he)
 axis_range(plot1,0.7,200.,10.,50.,True,False)
 axis_str(plot1,"",":F21:p [hPa]",":F21:z [km]")

 draw(plot1)

 plotb1 = xylineplot(wks,rs,zp,"D",0,"T",1.0,"C",lcol)

 set_size(plotb1,xf2,yf,wi,he)
 axis_range(plotb1,0.8e-3,0.4,10.,50.,True,False)
 axis_str(plotb1,"",":F33:r:F21: [kg m:S:-3:N:]","")

 draw(plotb1)

 dr_txt2(wks,"dot ::  p, :F33:r:F25: vs z @ 5:F34:0:F25:N:F30:-:F25:5:F34:0:F25:S in HadGEM2 (red) / MERRA (black)",xf+2,yf-12,25,8,0.,"CenterLeft",1)
 dr_txt2(wks,"line ::  log-p height where H = "+hs+" km, z:B:70:N: = 18.5 km (or, z:B:1000:N: = 1.0 km)",xf+2,yf-6,25,8,0.,"CenterLeft",1)

 frame(wks)

end



