load "$UTIL/ncl/header.ncl"

xspace = 170  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 1.0

begin

;--------------------------------------- parameters
 expname = (/"era-int","merra","cfsr","jra55","era-int_f","sonde_hr"/)
 varname = "u"

 y0 = 1979
 y1 = 1979
 y2 = 2010
 lat0 = 10
 lata = lat0 + 2.5

 wks = cr_wks("fig/ws_us_z-ya"+lat0,600,0)
 colind = NhlNewColor(wks,0.6,0.6,0.6)

;--------------------------------------- read data
 ne = dimsizes(expname)
 n0 = (y1-y0)*12
 n9 = (y2-y0+1)*12

 nt = n9 - n0

   fnam1m = "$DATD/clim/MERRA/"+expname(1)+".u_ypt.nc"
   fnam1e = "$DATD/clim/ERA-I/"+expname(0)+".u_ypt.nc"
   fnam1c = "$DATD/clim/CFSR/"+expname(2)+".u_ypt.nc"
   fnam1j = "$DATD/clim/JRA55/"+expname(3)+".u_ypt.nc"
   fnam1f = "$DATD/clim/ERA-I-nr/"+expname(4)+".u_ypt.nc"
   fnam1o = "/data18/kyh/qbo_obs/qbo_"+expname(5)+".nc"
   f1m = addfile(fnam1m,"r")
   f1e = addfile(fnam1e,"r")
   f1c = addfile(fnam1c,"r")
   f1j = addfile(fnam1j,"r")
   f1f = addfile(fnam1f,"r")
   f1o = addfile(fnam1o,"r")
   fzf = addfile("/data7/ERA-I-nr/ERA-Interim_coordvars.nc","r")
   zf0 = (fzf->a_model_alt)/100. + (fzf->b_model_alt)*1.e3
   zum = tofloat(18.6 + 6.6*log( 70./ f1m->levels ))
   zue = tofloat(18.6 + 6.6*log( 70./ f1e->lv_ISBL1 ))
   zuc = tofloat(18.6 + 6.6*log( 70./ f1c->plevel ))
   zuj = tofloat(18.6 + 6.6*log( 70./ f1j->lv_ISBL1 ))
   zuf = tofloat(18.6 + 6.6*log( 70./ zf0( f1f->level - 1 ) ))
print(f1f->level)
   zuo = 18.6 + 6.6*log( 70./ f1o->p )

   nzm = dimsizes(zum)
   nze = dimsizes(zue)
   nzc = dimsizes(zuc)
   nzj = dimsizes(zuj)
   nzf = dimsizes(zuf)
   nzo = dimsizes(zuo)

   latm = tofloat(f1m->latitude({-lata:lata}))
   late = tofloat(f1e->g0_lat_2({-lata:lata}))
   latc = tofloat(f1c->latitude({-lata:lata}))
   latj = tofloat(f1j->g0_lat_2({-lata:lata}))
   latf = tofloat(f1f->latitude({-lata:lata}))

   u0m = new((/nzm,nt/),float)
   u0e = new((/nze,nt/),float)
   u0c = new((/nzc,nt/),float)
   u0j = new((/nzj,nt/),float)
   u0f = new((/nzf,nt/),float)
   u0o = new((/nzo,nt/),float)
   do n=0, nt-1
     do k=0, nzm-1
       u0m(k,n) = ftcurvi(-lat0,lat0,latm, tofloat(f1m->u(n0+n,k,{-lata:lata})))
     end do
     do k=0, nze-1
       u0e(k,n) = ftcurvi(-lat0,lat0,late, f1e->U_GDS0_ISBL(n0+n,k,{-lata:lata}))
     end do
     do k=0, nzc-1
       u0c(k,n) = ftcurvi(-lat0,lat0,latc, f1c->UGRD(n0+n,k,{-lata:lata}))
     end do
     do k=0, nzj-1
       u0j(k,n) = ftcurvi(-lat0,lat0,latj, f1j->U_GRD_GDS0_ISBL(n0+n,k,{-lata:lata}))
     end do
     do k=0, nzf-1
       u0f(k,n) = ftcurvi(-lat0,lat0,latf, tofloat(f1f->u(n0+n,k,{-lata:lata})))
     end do
;     u0o(:,n) = f1o->u(:,n)
   end do
   u0m = u0m / tofloat(lat0+lat0)
   u0e = u0e / tofloat(lat0+lat0)
   u0c = u0c / tofloat(lat0+lat0)
   u0j = u0j / tofloat(lat0+lat0)
   u0f = u0f / tofloat(lat0+lat0)

   fnam2m = "$DATD/tem/MERRA/"+expname(1)+".tem_ypt.nc"
   fnam2e = "$DATD/tem/ERA-I/"+expname(0)+".tem_ypt.nc"
   fnam2c = "$DATD/tem/CFSR/"+expname(2)+".tem_ypt.nc"
   fnam2j = "$DATD/tem/JRA55/"+expname(3)+".tem_ypt.nc"
   fnam2f = "$DATD/tem/ERA-I-nr/"+expname(4)+".tem_ypt.nc"
   f2m = addfile(fnam2m,"r")
   f2e = addfile(fnam2e,"r")
   f2c = addfile(fnam2c,"r")
   f2j = addfile(fnam2j,"r")
   f2f = addfile(fnam2f,"r")
   zam = 18.6 + 6.6*log( 70./ f2m->p )
   zae = 18.6 + 6.6*log( 70./ f2e->p )
   zac = 18.6 + 6.6*log( 70./ f2c->p )
   zaj = 18.6 + 6.6*log( 70./ f2j->p )
   zaf = 18.6 + 6.6*log( 70./ f2f->p )

   wsm = transpose( ftcurvi(-lat0,lat0,latm, f2m->w_res(n0:n9-1,:,{-lata:lata})) )
   wse = transpose( ftcurvi(-lat0,lat0,late, f2e->w_res(n0:n9-1,:,{-lata:lata})) )
   wsc = transpose( ftcurvi(-lat0,lat0,latc, f2c->w_res(n0:n9-1,:,{-lata:lata})) )
   wsj = transpose( ftcurvi(-lat0,lat0,latj, f2j->w_res(n0:n9-1,:,{-lata:lata})) )
   wsf = transpose( ftcurvi(-lat0,lat0,latf, f2f->w_res(n0:n9-1,:,{-lata:lata})) )
   wsm = wsm / tofloat(lat0+lat0)
   wse = wse / tofloat(lat0+lat0)
   wsc = wsc / tofloat(lat0+lat0)
   wsj = wsj / tofloat(lat0+lat0)
   wsf = wsf / tofloat(lat0+lat0)


 wsam = dim_avg(wsm)
 wsae = dim_avg(wse)
 wsac = dim_avg(wsc)
 wsaj = dim_avg(wsj)
 wsaf = dim_avg(wsf)


 usm = u0m*0.
 use = u0e*0.
 usc = u0c*0.
 usj = u0j*0.
 usf = u0f*0.
 uso = u0o*0.
 do k=1, nzm-2
   usm(k,:) = (u0m(k+1,:) - u0m(k-1,:))/(zum(k+1) - zum(k-1))
 end do
 do k=1, nze-2
   use(k,:) = (u0e(k+1,:) - u0e(k-1,:))/(zue(k+1) - zue(k-1))
 end do
 do k=1, nzc-2
   usc(k,:) = (u0c(k+1,:) - u0c(k-1,:))/(zuc(k+1) - zuc(k-1))
 end do
 do k=1, nzj-2
   usj(k,:) = (u0j(k+1,:) - u0j(k-1,:))/(zuj(k+1) - zuj(k-1))
 end do
 do k=1, nzf-2
   usf(k,:) = (u0f(k+1,:) - u0f(k-1,:))/(zuf(k+1) - zuf(k-1))
 end do
 do k=1, nzo-2
   uso(k,:) = (u0o(k+1,:) - u0o(k-1,:))/(zuo(k+1) - zuo(k-1))
 end do

 ussm = dim_stddev(usm)
 usse = dim_stddev(use)
 ussc = dim_stddev(usc)
 ussj = dim_stddev(usj)
 ussf = dim_stddev(usf)
 usso = dim_stddev(uso)
 usso@_FillValue = 0.

;--------------------------------------- process 0
 wsam = wsam*1.e3
 wsae = wsae*1.e3
 wsac = wsac*1.e3
 wsaj = wsaj*1.e3
 wsaf = wsaf*1.e3

;--------------------------------------- cn level
 lcol = (/4,2,3,13,1,colind/)  ; default
 ltck = (/0.8,0.8,0.8,0.8,0.8,2.0/)

;--------------------------------------- contour
 wi = 35
 he = 35
 xf = (/15,56/)
 yf = 10

 ystr = "height (km)"
 xstr = (/"W* (mm s:S:-1:N:)","Std. Dev. of :F18:s:F22:U/:F18:s:F22:z (m s:S:-1:N: / km)"/)
 txtstr = (/"ERA-I","MERRA","CFSR","JRA-55","ERA-I_ml","sonde"/)

do ie=0, ne-2 ;ne-1

 if (ie .eq. 1) then
   plot = xylineplot(wks,ussm,zum,"D",0,"T",ltck(ie),"C",lcol(ie))
 end if
 if (ie .eq. 0) then
   plot = xylineplot(wks,usse,zue,"D",0,"T",ltck(ie),"C",lcol(ie))
 end if
 if (ie .eq. 2) then
   plot = xylineplot(wks,ussc,zuc,"D",0,"T",ltck(ie),"C",lcol(ie))
 end if
 if (ie .eq. 3) then
   plot = xylineplot(wks,ussj,zuj,"D",0,"T",ltck(ie),"C",lcol(ie))
 end if
 if (ie .eq. 4) then
   plot = xylineplot(wks,ussf,zuf,"D",0,"T",ltck(ie),"C",lcol(ie))
 end if
 if (ie .eq. 5) then
   plot = xylineplot(wks,usso,zuo,"D",0,"T",ltck(ie),"C",lcol(ie))
 end if
 xyaddmarker(plot,16,"S",1.0,"T",1.0,"C",lcol(ie))
 xyaddline(plot,0.,zue,"D",0,"T",0.6,"C",colind)

 set_size(plot,xf(1),yf,wi,he)
 axis_range(plot,1.0,5.0,17,35,False,False)
 axis_str(plot,"",xstr(1),"")

 draw(plot)

 if ( ie .eq. ne-1 ) then
   continue
 end if

 if (ie .eq. 1) then
   plot = xylineplot(wks,wsam,zam,"D",0,"T",ltck(ie),"C",lcol(ie))
 end if
 if (ie .eq. 0) then
   plot = xylineplot(wks,wsae,zae,"D",0,"T",ltck(ie),"C",lcol(ie))
 end if
 if (ie .eq. 2) then
   plot = xylineplot(wks,wsac,zac,"D",0,"T",ltck(ie),"C",lcol(ie))
 end if
 if (ie .eq. 3) then
   plot = xylineplot(wks,wsaj,zaj,"D",0,"T",ltck(ie),"C",lcol(ie))
 end if
 if (ie .eq. 4) then
   plot = xylineplot(wks,wsaf,zaf,"D",0,"T",ltck(ie),"C",lcol(ie))
 end if
 xyaddmarker(plot,16,"S",1.0,"T",1.0,"C",lcol(ie))
 xyaddline(plot,0.,zae,"D",0,"T",0.6,"C",colind)

 set_size(plot,xf(0),yf,wi,he)
 axis_range(plot,-0.2,1.0,17,35,False,False)
 axis_str(plot,"",xstr(0),ystr)

 draw(plot)

 ps = (/100,70,50,40,30,20,10,7,5,3,2,1/)
 zs = 18.6 + 6.6*log(70./ps)
 plotx = xylineplot(wks,(/-1,0,1/),(/-1,0,1/),"D",0,"T",0.001,"C",-1)
 set_size(plotx,xf(1)+wi+4,yf,1,he)
 axis_range(plotx,-999,999,17,35,False,False)
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

; lg = lgplot(wks,txtstr,"D",0,"T",ltck,"C",lcol,2)
; set_size_lg(lg,xf(1),yf,wi*0.25,3.5*6)
; draw(lg)

 lg = lgplot(wks,txtstr(:ne-2),"D",0,"T",ltck,"C",lcol,2)
 set_size_lg(lg,xf(0)+wi*0.60,yf+he-3.0*5+0.5,wi*0.20,3.0*5)
 draw(lg)

; dr_txt(wks,"(a)",xf(0)-7,yf(0)-2,25,8,0.)
; dr_txt(wks,"(b)",xf(1)-7,yf(0)-2,25,8,0.)

 frame(wks)

end

