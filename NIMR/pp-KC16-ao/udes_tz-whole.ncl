load "$UTIL/ncl/header.ncl"
xspace =  83  ; 83, 170  -  QJRMS
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 0.85

begin
;--------------------------------------- parameters
 cexp    = (/"CGW","MERRA"/)
 expname = "tc"
 ppid    = "j"
 varname = "u"

 y1 = 1950 ;2058 ;1960
 y2 = 2098 ;2097 ;1999
 lat0 = 5.
 opt_ss = 3  ; 0, 1, 2, 3 : no, mean, linear, and quadratic seasonal cycle
 zdot = 21.0  ; p0 = 50.

 opt_ra = 0
 y1r = 1979
 y2r = y2

 wks = cr_wks("fig/udes_tz-whole-"+expname(0),600,0)
 gsn_define_colormap(wks,"StepSeq25")
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- read data

; HadGEM

 fdir  = "$DATD/AOL60"+cexp(0)
 fname = "xm.uan"+expname+"a.p"+ppid+".nc"

 f0 = addfile(fdir+"/"+fname,"r")
 tmp = f0->$varname$
 z = tmp&$tmp!1$ / 1.e3
 nz = dimsizes(z)
 t = tmp&$tmp!0$
 t_char = tochar(t@units)
 if ( t_char(0:10) .ne. "days since " ) then
   print(t@units)
   print("Check this. It should be 'days since ...'.")
   exit
 end if
 y0i = toint(tostring(t_char(11:14)))
 m0i = toint(tostring(t_char(16:17)))
 d0i = toint(tostring(t_char(19:20)))
 print("date0 : "+y0i+" "+m0i+" "+d0i)

 t = y0i + (m0i-1)/12. + (d0i-1)/360. + t/360.
 nt = dimsizes(t)
 print("data : "+t(0)+" - "+t(nt-1))

 u = transpose( dim_avg( tmp(:,:,{-lat0:lat0}) ) )

 delete(tmp)

if (opt_ra .ne. 0) then

; MERRA

 print("MERRA")
 f0 = addfile("$DATD/clim/MERRA/merra.u_ypt."+y1r+"-2012.nc","r")
 latnam = f0->$varname$!2
 lat = f0->$latnam$({-lat0-1.0:lat0+1.0})
 print(lat(0)+" - "+lat(dimsizes(lat)-1))
 wgt = tofloat(lat)*0. + 1.
 wgt(0) = 0.5
 wgt(dimsizes(wgt)-1) = 0.5
 pnam = f0->$varname$!1
 p = tofloat( f0->$pnam$({1:100}) )
 nt0 = dimsizes( f0->$varname$&$(f0->$varname$!0)$ )
 u0 = transpose( dim_avg_wgt( \
      tofloat( f0->$varname$(:,{1:100},{-lat0-1.0:lat0+1.0}) ), wgt, 0 ) )
 t0 = tofloat(y1r) + fspan(1./24.,nt0/12.-1./24.,nt0)
 zp = 18.6 + 6.6*log(70./p)

 fh = addfile("$DATD/clim/MERRA/merra.h_p."+y1r+"-2012.nc","r")
 zm = fh->h({1:100}) / 1.e3

; z0 = zp
 z0 = zm
 nz0 = dimsizes(z0)

end if

;--------------------------------------- process 0
 if (opt_ss .ne. 0) then
   udes = u
   if (opt_ss .eq. 1) then  ; mean seasonal cycle
     uss = new((/nz,12/),float)
     um = dim_avg(u)
     do n=0, 11
       itmps = ispan(n,nt-1,12)
       uss(:,n) = dim_avg(u(:,itmps))
       udes(:,itmps) = u(:,itmps) - conform_dims(dimsizes(u(:,itmps)),uss(:,n)-um,0)
       delete(itmps)
     end do
   else
     if (opt_ss .eq. 2) then  ; seasonal cycle having a linear trend
       do n=0, 11
         itmps = ispan(n,nt-1,12)
         udes(:,itmps) = dtrend(u(:,itmps),False)
         delete(itmps)
       end do
       uss = u - udes
       utmp = uss - dtrend(uss,False)
     else if (opt_ss .eq. 3) then  ; seasonal cycle having a quadratic trend
       do n=0, 11
         itmps = ispan(n,nt-1,12)
         udes(:,itmps) = dtrend_quadratic(u(:,itmps),0)
         delete(itmps)
       end do
       uss = u - udes
       utmp = uss - dtrend_quadratic(uss,0)
     end if
     end if
     udes = udes + utmp
   end if 

   u = udes

 end if


if (opt_ra .ne. 0) then

; MERRA

 u0des = u0
 u0ss = new((/nz0,12/),float)
; u0m = dim_avg(u0)
 u0m = 0.
 do n=0, 11
   itmps = ispan(n,nt0-1,12)
   u0ss(:,n) = dim_avg(u0(:,itmps))
   u0des(:,itmps) = u0(:,itmps) - conform_dims(dimsizes(u0(:,itmps)),u0ss(:,n)-u0m,0)
   delete(itmps)
 end do

 u0 = u0des

end if

 do k=0, nz-1
   if (abs(z(k)-zdot)/zdot .le. 0.02) then
     kk1 = k
   end if
 end do

 ts = new(300,float,1.e20)
 zz = new(300,float,1.e20)
 cnt = 0
 do n=0, nt-2
   if ( ismissing(u(kk1,n+1)) ) then
     break
   end if
   if ( u(kk1,n)*u(kk1,n+1) .lt. 0. ) then
     x = ( abs(u(kk1,n))*t(n+1) + abs(u(kk1,n+1))*t(n) )/ \
         (abs(u(kk1,n))+abs(u(kk1,n+1)))
     ts(cnt) = x
     zz(cnt) = z(kk1)
     cnt = cnt + 1
   end if
 end do

;--------------------------------------- cn level
 cnlev = fspan(-40,40,21)
; cnlev = fspan(-40,28,18)
; cnfil = (/2,3,4,5,6,7,8,9,10,11,0,0,26,25,24,23,22,21,20,19,18,17/)
 cnfil = (/17,18,19,20,21,22,23,24,25,26,0,0,11,10,9,8,7,6,5,4,3,2/)

;--------------------------------------- contour
 yl = 20
 yls = 1;2

 nfig = 1
 if (opt_ra .ne. 0) then
   nfig = 2
 end if

do ifig=0, nfig-1

 wi = 5*yl
 he = 14
 xf = 20
 yf = 10 + ispan(0,20,1)*(he+5)

 if (ifig .eq. 0) then
   varf = u
   varx = t
   vary = z
   ys = ispan(y1+yls,y2,yl)
 else
   varf = u0
   varx = t0
   vary = z0
   ys = ispan(y1r+yls,y2r,yl)
   yf = yf + nl*(he+5) + 7
 end if
; ys(0) = ys(0) - yls

 nl = dimsizes(ys)

 wi2 = tofloat(wi)/yl*yls

 xf0 = xf - wi2
 wi0 = wi + wi2

 xstr  = new(nl,string)
 lbopt = new(nl,string)

 mstr = ":F25:"+(/"(a)","(b)"/)
 mstr(:) = mstr(:)+"   :F22:"+cexp
 xstr(:) = ""
 xstr(nl-1) = "year"
 ystr = "height (km)"
 lbopt(:) = ""
 if (ifig .eq. nfig-1) then
   lbopt(nl-1) = "(m s:S:-1:N:)"
 end if


do il=0, nl-1

 contour1 = cnshadeplot(wks,varx,vary,varf,cnlev,cnfil,lbopt(il))
 cnaddline(contour1,wks,varx,vary,varf,(/0./),0,"D",0,"T",0.8,"C",1)
 setvalues contour1
;  "pmLabelBarOrthogonalPosF" : 0.18 ;0.12
;  "pmLabelBarParallelPosF"   : 0.50
  "pmLabelBarOrthogonalPosF" : -0.8
  "pmLabelBarParallelPosF"   : 1.90
  "pmLabelBarKeepAspect"     : False
  "pmLabelBarWidthF"         : 0.6
  "pmLabelBarHeightF"        : 0.2
  "lbLabelFontHeightF" : 0.015
  "lbTitlePosition"    : "Right"
  "lbTitleDirection"   : "Across"
  "lbTitleOffsetF"     : 0.015
  "lbTitleFontHeightF" : 0.018
 end setvalues

 plot1 = xylineplot(wks,ts,zz,"D",0,"T",1.0,"C",-1)
 xyaddmarker(plot1,16,"S",0.8,"T",1.0,"C",1)

 ye = ys(il)+yl
 xf1 = tofloat(xf)
 wi1 = tofloat(wi)
; if (il .eq. 0) then
;   ye = ye + yls
;   xf1 = xf0
;   wi1 = wi0
; end if
 if (ye .gt. y2+1) then
   ye = y2+1
   wi1 = wi*(tofloat(y2+1-ys(il))/yl)
 end if

 set_size(contour1,xf1,yf(il),wi1,he)
 axis_range(contour1,ys(il),ye,17,39,False,False)
 axis_str(contour1,"",xstr(il),ystr)

 set_size(plot1,xf1,yf(il),wi1,he)
 axis_range(plot1,ys(il),ye,17,39,False,False)
 axis_str(plot1,"","","")

 ps = (/100,70,50,30,20,10,7,5,3,2,1/)
 zs = 18.6 + 6.6*log(70./ps)
 plotx = xylineplot(wks,(/-1,0,1/),(/-1,0,1/),"D",0,"T",0.001,"C",-1)
 set_size(plotx,xf1+wi1+2,yf(il),1,he)
 axis_range(plotx,-999,999,17,39,False,False)
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

 xtck = ispan(ys(il),ye,1)
 setvalues contour1
  "tmXBLabelJust" : "TopLeft"
  "tmXBMode" : "Explicit"
  "tmXBValues" : xtck
  "tmXBLabels" : where(xtck .eq. ye, "", "  "+xtck)
 end setvalues
 delete(xtck)
 setvalues plot1
  "tmXBOn" : False
  "tmXTOn" : False
  "tmYLOn" : False
  "tmYROn" : False
 end setvalues


 if (il .gt. 0) then

  setvalues contour1
   "tmYLOn" : False
   "tiYAxisString" : ""
  end setvalues

  contour2 = cnshadeplot(wks,varx,vary,varf,cnlev,cnfil,"")
  cnaddline(contour2,wks,varx,vary,varf,(/0./),0,"D",0,"T",0.8,"C",1)

  plot2 = xylineplot(wks,ts,zz,"D",0,"T",1.0,"C",-1)
  xyaddmarker(plot2,16,"S",0.8,"T",1.0,"C",1)

  set_size(contour2,xf-wi2,yf(il),wi2,he)
  axis_range(contour2,ys(il)-yls,ys(il),17,39,False,False)
  axis_str(contour2,"","",ystr)

  set_size(plot2,xf-wi2,yf(il),wi2,he)
  axis_range(plot2,ys(il)-yls,ys(il),17,39,False,False)
  axis_str(plot2,"","","")

  setvalues contour2
   "tmXBLabelsOn" : False
   "tmXBMode" : "Manual"
   "tmXBTickSpacingF"  : 1
   "tmXBMinorPerMajor" : 0
   "tmYROn" : False
  end setvalues
  setvalues plot2
   "tmXBOn" : False
   "tmXTOn" : False
   "tmYLOn" : False
   "tmYROn" : False
  end setvalues
  getvalues plot1
   "xyCoordDataSpec" : spec
  end getvalues
  getvalues spec
   "xyMarkerSizeF" : msize
  end getvalues
  getvalues plot2
   "xyCoordDataSpec" : spec
  end getvalues
  setvalues spec
   "xyMarkerSizeF" : msize
  end setvalues

  cnfinalize(contour2,wks)

  draw(contour2)
;  draw(plot2)

 end if


 cnfinalize(contour1,wks)

 draw(contour1)
; draw(plot1)
 draw(plotx)

end do  ; il

 if (nfig .gt. 1) then
   dr_txt2(wks,mstr(ifig),xf-wi2-5,yf(0)-4,22,6,0.,"CenterLeft",1)
 end if

 delete(varx)
 delete(vary)
 delete(varf)
 delete(ys)
 delete(mstr)
 delete(xstr)
 delete(lbopt)

end do  ; ifig

 frame(wks)

end

