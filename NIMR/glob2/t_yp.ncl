load "$UTIL/ncl/header.ncl"

xspace =  83  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 1.2

begin

;--------------------------------------- parameters
 jobname = (/"tk","tc"/)
 expname = (/"AOL60CTL","AOL60CGW"/)
 ppid    = "j"
 varname = "T_p_uv"

 y0 = 1951
 y1 = 1951
 y2 = 2009
 mms = "djf"

 wks = cr_wks("fig/t_yp-"+mms+"_"+y1+"-"+y2+"-"+jobname(1),600,0)
 gsn_define_colormap(wks,"StepSeq25")
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- read data
 ne = dimsizes(jobname)
 nyr = y2 - y1 + 1

 mm3s = (/"djf","mam","jja","son"/)
 mm3i = (/ (/1,2,12/), (/3,4,5/), (/6,7,8/), (/9,10,11/) /)
 do is=0, dimsizes(mm3s)-1
   if ( mms .eq. mm3s(is) ) then
     mm = mm3i(is,:)
   end if
 end do
 nm = dimsizes(mm)

do ie=0, ne-1

 fname = "../dat/"+expname(ie)+"/xm.uan"+jobname(ie)+"a.p"+ppid+".nc"
 f0 = addfile(fname,"r")
 n0 = (y1 - y0)*12

 ; get y, z
 tmp = f0->$varname$(n0:n0+nyr*12-1,:,:)
 lat = tmp&$tmp!2$
 p = tmp&$tmp!1$
 z = 7.0*log(1.e3/p)
 ny = dimsizes(lat)
 nz = dimsizes(p)

 if (ie .eq. 0) then
   vara = new((/ne,nm,nyr,nz,ny/),float)
 end if
 do im=0, nm-1
   vara(ie,im,:,:,:) = tmp(ispan(mm(im)-1,nyr*12-1,12),:,:)
 end do

 delete(tmp)

end do  ; ie

;--------------------------------------- process 0
 var0 = dim_avg_n(dim_avg_n(vara,2),1)
 vard = var0(1,:,:) - var0(0,:,:)

 var = array_append_record( var0, \
          conform_dims((/1,nz,ny/),vard,(/1,2/)), 0 )

;--------------------------------------- cn level
 scl = 4.
 cnlev = ispan(-9,10,1)*scl+250.
; cnlev = fspan(-45,50,20)+250.
 cnlev2 = array_append_record( ispan(-20,-9,1), ispan(10,25,1), 0 )*scl+250.
 cnlevd = ispan(-10,10,1)*scl*0.05
 cnlevd2 = array_append_record( ispan(-25,-10,1), ispan(10,25,1), 0 )*scl*0.05

; cnfil = (/2,3,4,5,6,7,8,9,10,11,0,0,26,25,24,23,22,21,20,19,18,17/)
 cnfil = (/17,18,19,20,21,22,23,24,25,26,0,11,10,9,8,7,6,5,4,3,2/)
 cnfild = (/17,18,19,20,21,22,23,24,25,26,0,0,11,10,9,8,7,6,5,4,3,2/)

;--------------------------------------- contour
 wi = 32
 he = 40
 xf = 10 + (/0,1,2/)*40
 yf = 10

 ystr  = new(ne+1,string)

 mstr = (/"CTL","CC","CC - CTL"/)
 xstr = "latitude (deg)"
 ystr(:) = ""
 ystr(0) = "log-p height (km)"
 lbopt = (/"(K)","","(K)"/)

do ie=0, ne-1

 contour1 = cnshadeplot(wks,lat,z,var(ie,:,:),cnlev,cnfil,lbopt(ie))
 cnaddline(contour1,wks,lat,z,var(ie,:,:),cnlev2,0,"D",0,"T",0.5,"C",1)
 setvalues contour1
  "pmLabelBarOrthogonalPosF" : 0.08
  "pmLabelBarParallelPosF"   : 1.15
  "pmLabelBarKeepAspect"     : False
  "pmLabelBarWidthF"         : 1.2
  "pmLabelBarHeightF"        : 0.12
  "lbLabelFontHeightF" : 0.026
  "lbTitlePosition"    : "Right"
  "lbTitleDirection"   : "Across"
  "lbTitleOffsetF"     : 0.02
  "lbTitleFontHeightF" : 0.032
 end setvalues

 set_size(contour1,xf(ie),yf,wi,he)
 axis_range(contour1,-90,90,0,50,False,False)
 axis_str(contour1,mstr(ie),xstr,ystr(ie))

 cnfinalize(contour1,wks)

 draw(contour1)

end do  ; ie

 ie = ne

 contour2 = cnshadeplot(wks,lat,z,var(ie,:,:),cnlevd,cnfild,lbopt(ie))
 cnaddline(contour2,wks,lat,z,var(ie,:,:),cnlevd2,0,"D",0,"T",0.5,"C",1)
 setvalues contour2
  "pmLabelBarOrthogonalPosF" : 0.08
  "pmLabelBarKeepAspect"     : False
  "pmLabelBarWidthF"         : 0.8
  "pmLabelBarHeightF"        : 0.12
  "lbLabelFontHeightF" : 0.026
  "lbTitlePosition"    : "Right"
  "lbTitleDirection"   : "Across"
  "lbTitleOffsetF"     : 0.02
  "lbTitleFontHeightF" : 0.032
 end setvalues

 set_size(contour2,xf(ie),yf,wi,he)
 axis_range(contour2,-90,90,0,50,False,False)
 axis_str(contour2,mstr(ie),xstr,ystr(ie))

 cnfinalize(contour2,wks)

 draw(contour2)


 frame(wks)

end

