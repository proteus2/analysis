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
 varname = "u"

 y0 = 1951
 y1 = 1951
 y2 = 1952
 mms = "djf"

 wks = cr_wks("fig/rc-w_yz-lb",600,0)
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
 z = tmp&$tmp!1$ / 1.e3
 ny = dimsizes(lat)
 nz = dimsizes(z)

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

 var = var@_FillValue
;--------------------------------------- cn level
 scl = 0.4
 cnlev = ispan(-10,10,1)*scl
 cnlev2 = cnlev*10.
; cnlev2 = array_append_record( ispan(-25,-10,1), ispan(10,25,1), 0 )*scl
 cnlevd = cnlev*0.1
 cnlevd2 = cnlev2*0.1

; cnfil = (/2,3,4,5,6,7,8,9,10,11,0,0,26,25,24,23,22,21,20,19,18,17/)
 cnfil = (/17,18,19,20,21,22,23,24,25,26,0,0,11,10,9,8,7,6,5,4,3,2/)
 
;--------------------------------------- contour
 wi = 32
 he = 40
 xf = 10 + (/0,1,2/)*40
 yf = 10

 ystr  = new(ne+1,string)

 mstr = (/"CTL","CC","CC - CTL"/)
 xstr = "latitude (deg)"
 ystr(:) = ""
 ystr(0) = "hybrid height (km)"
 lbopt = (/"(mm s:S:-1:N:)","","(mm s:S:-1:N:)"/)

do ie=0, ne-1

 contour1 = cnshadeplot(wks,lat,z,var(ie,:,:),cnlev,cnfil,lbopt(ie))
; cnaddline(contour1,wks,lat,z,var(ie,:,:),cnlev2,0,"D",0,"T",0.5,"C",(/1,-1,1/))
 setvalues contour1
;  "pmLabelBarOrthogonalPosF" : 0.08
;  "pmLabelBarParallelPosF"   : 1.15
  "pmLabelBarKeepAspect"     : False
  "pmLabelBarWidthF"         : 0.10;1.2
  "pmLabelBarHeightF"        : 1.0;0.12
  "lbOrientation"      : "Vertical"
  "lbLabelFontHeightF" : 0.026
  "lbTitlePosition"    : "Right"
  "lbTitleDirection"   : "Across"
  "lbTitleAngleF"      : 270
  "lbTitleOffsetF"     : 0.25
  "lbTitleFontHeightF" : 0.032
 end setvalues

 set_size(contour1,xf(ie),yf,wi,he)
 axis_range(contour1,-90,90,0,50,False,False)
 axis_str(contour1,"","","")

 cnfinalize(contour1,wks)

 draw(contour1)

end do  ; ie

 ie = ne

 contour2 = cnshadeplot(wks,lat,z,var(ie,:,:),cnlevd,cnfil,lbopt(ie))
; cnaddline(contour2,wks,lat,z,var(ie,:,:),cnlevd2,0,"D",0,"T",0.5,"C",(/1,-1,1/))
 setvalues contour2
;  "pmLabelBarOrthogonalPosF" : 0.08
  "pmLabelBarKeepAspect"     : False
  "pmLabelBarWidthF"         : 0.1;0.8
  "pmLabelBarHeightF"        : 1.0;0.12
  "lbOrientation"      : "Vertical"
  "lbLabelFontHeightF" : 0.026
  "lbTitlePosition"    : "Right"
  "lbTitleDirection"   : "Across"
  "lbTitleAngleF"      : 270
  "lbTitleOffsetF"     : 0.25
  "lbTitleFontHeightF" : 0.032
 end setvalues

 set_size(contour2,xf(ie),yf,wi,he)
 axis_range(contour2,-90,90,0,50,False,False)
 axis_str(contour2,"","","")

 cnfinalize(contour2,wks)

 draw(contour2)


 frame(wks)

end

