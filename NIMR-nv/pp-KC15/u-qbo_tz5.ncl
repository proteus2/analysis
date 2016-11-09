load "$UTIL/ncl/header.ncl"

xspace =  83  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 1.0

begin

;--------------------------------------- parameters
 cexp    = (/"CC"/)
 expname = (/"uj"/)
 varname = "u"

 y1 = 1953
 y2 = 2006
 lat0 = 5.
 z0 = 21.0  ; p0 = 50.

 wks = cr_wks("fig/u-qbo_tz5-"+expname(0),600,0)
 gsn_define_colormap(wks,"StepSeq25")
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- read data
 ne = dimsizes(expname)
 nt = (y2-y1+1)*12

 fdir  = "../dat/L60CGW-t"
 fname = fdir+"/uan"+expname+".u_qbo.nc"

 ; get z
 f0 = addfile(fname(0),"r")
 tmp = f0->$varname$
 zname = tmp!1
 z = tmp&$zname$ / 1.e3
 nz = dimsizes(z)
 delete(tmp)

 u = new((/ne,nz,nt/),float,1.e20)
 do ie=0, ne-1
   f1 = addfile(fname(ie),"r")
   tmp = f1->$varname$
   u(ie,:,:) = transpose( tmp )
 end do

 t = y1 + (fspan(1,nt,nt)-0.5)/12.

; u = u - conform_dims( dimsizes(u), dim_avg(u), (/0,1/) )

;--------------------------------------- process 0
 do k=0, nz-1
   if (abs(z(k)-z0)/z0 .le. 0.02) then
     kk1 = k
   end if
 end do

 ts = new((/ne,100/),float,1.e20)
 zz = new((/ne,100/),float,1.e20)
 do ie=0, ne-1
   cnt = 0
   do n=0, nt-2
     if ( ismissing(u(ie,kk1,n+1)) ) then
       break
     end if
     if ( u(ie,kk1,n)*u(ie,kk1,n+1) .lt. 0. ) then
       x = ( abs(u(ie,kk1,n))*t(n+1) + abs(u(ie,kk1,n+1))*t(n) )/ \
           (abs(u(ie,kk1,n))+abs(u(ie,kk1,n+1)))
       ts(ie,cnt) = x
       zz(ie,cnt) = z(kk1)
       cnt = cnt + 1
     end if
   end do
 end do

;--------------------------------------- cn level
; cnlev = fspan(-40,40,21)
 cnlev = fspan(-40,24,17)
; cnfil = (/2,3,4,5,6,7,8,9,10,11,0,0,26,25,24,23,22,21,20,19,18,17/)
 cnfil = (/17,18,19,20,21,22,23,24,25,26,0,0,11,10,9,8,7,6,5,4,3,2/)
print(max(u(0,1:,:))+"  "+min(u(0,1:,:)))

;--------------------------------------- contour
 nl = 5
 ys = (/1953,1964,1975,1986,1997/)
 yl = 11
 yls = 2

 wi = 60
 he = 15
 xf = 20
 yf = 10 + (/0,1,2,3,4/)*20

 wi2 = tofloat(wi)/yl*yls

 mstr  = new(nl,string)
 xstr  = new(nl,string)
 ystr  = new(nl,string)
 lbopt = new(nl,string)

 mstr(:) = ""
; mstr(0) = cexp(0)
 xstr(:) = ""
 xstr(nl-1) = "year"
 ystr(:) = "z (km)"
 lbopt(:) = ""
 lbopt(nl-1) = "(m s:S:-1:N:)"

do il=0, nl-1

 contour1 = cnshadeplot(wks,t,z,u(0,:,:),cnlev,cnfil,lbopt(il))
 cnaddline(contour1,wks,t,z,u(0,:,:),(/0./),0,"D",0,"T",0.8,"C",1)
 setvalues contour1
  "pmLabelBarOrthogonalPosF" : 0.12
  "pmLabelBarParallelPosF"   : 0.50+0.05
  "pmLabelBarKeepAspect"     : False
  "pmLabelBarWidthF"         : 0.6
  "pmLabelBarHeightF"        : 0.2
  "lbLabelFontHeightF" : 0.015
  "lbTitlePosition"    : "Right"
  "lbTitleDirection"   : "Across"
  "lbTitleOffsetF"     : 0.015
  "lbTitleFontHeightF" : 0.018
 end setvalues

 ye = ys(il)+yl
 wi1 = tofloat(wi)
 if (ye .gt. y2+1) then
   ye = y2+1
   wi1 = wi*(tofloat(y2+1-ys(il))/yl)
 end if
 set_size(contour1,xf,yf(il),wi1,he)
 axis_range(contour1,ys(il),ye,17,39,False,False)
 axis_str(contour1,mstr(il),xstr(il),ystr(il))
 setvalues contour1
  "tmXBMode" : "Manual"
  "tmXBTickSpacingF"  : 1
  "tmXBMinorPerMajor" : 0
 end setvalues

if (il .gt. 0) then
  setvalues contour1
   "tmYLOn" : False
   "tiYAxisString" : ""
  end setvalues

  contour2 = cnshadeplot(wks,t,z,u(0,:,:),cnlev,cnfil,"")
  cnaddline(contour2,wks,t,z,u(0,:,:),(/0./),0,"D",0,"T",0.8,"C",1)
  set_size(contour2,xf-wi2,yf(il),wi2,he)
  axis_range(contour2,ys(il)-yls,ys(il),17,39,False,False)
  axis_str(contour2,"","",ystr(il))
  setvalues contour2
   "tmXBLabelsOn" : False
   "tmXBMode" : "Manual"
   "tmXBTickSpacingF"  : 1
   "tmXBMinorPerMajor" : 0
   "tmYROn" : False
  end setvalues
  cnfinalize(contour2,wks)
  draw(contour2)
end if

 cnfinalize(contour1,wks)

 draw(contour1)


 plot = xylineplot(wks,ts(0,:),zz(0,:),"D",0,"T",1.0,"C",-1)
 xyaddmarker(plot,16,"S",0.8,"T",1.0,"C",1)

 set_size(plot,xf,yf(il),wi1,he)
 axis_range(plot,ys(il),ye,17,39,False,False)
 axis_str(plot,"","","")
 setvalues plot
  "tmXBOn" : False
  "tmXTOn" : False
  "tmYLOn" : False
  "tmYROn" : False
 end setvalues

if (il .gt. 0) then
  plot2 = xylineplot(wks,ts(0,:),zz(0,:),"D",0,"T",1.0,"C",-1)
  xyaddmarker(plot2,16,"S",0.8,"T",1.0,"C",1)
  set_size(plot2,xf-wi2,yf(il),wi2,he)
  axis_range(plot2,ys(il)-yls,ys(il),17,39,False,False)
  axis_str(plot2,"","","")
  setvalues plot2
   "tmXBOn" : False
   "tmXTOn" : False
   "tmYLOn" : False
   "tmYROn" : False
  end setvalues
  getvalues plot
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
;  draw(plot2)
end if

; draw(plot)

end do


 frame(wks)

end

