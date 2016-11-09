load "$UTIL/ncl/header.ncl"

xspace =  83  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 0.9

begin

;--------------------------------------- parameters
 jobname = (/"tc","tk"/)
 expname = (/"AOL60CGW","AOL60CTL"/)
 varname = (/"v_res","v_res","w_res"/)

 y0 = 1951
 y1 = 1960
 y2 = 1999;2009
 mms = "djf"

 wks = cr_wks("fig/rc2v_yz-"+mms+"_chg-"+jobname(1),600,0)
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

 fname = "../dat/"+expname(ie)+"/tem/uan"+jobname(ie)+".tem-wc_yzt.nc"
 f0 = addfile(fname,"r")
 n0 = (y1 - y0)*12

 ; get y, z
 tmp = f0->$varname(0)$(n0:n0+nyr*12-1,:,:)
 lat = tmp&$tmp!2$
 z = tmp&$tmp!1$ / 1.e3
 ny = dimsizes(lat)
 nz = dimsizes(z)
tmp = f0->$varname(0)$(n0+12*98:n0+12*98+nyr*12-1,:,:) - tmp

 epf = (/ f0->$varname(1)$(n0+12*98:n0+12*98+nyr*12-1,:,:) - \
          f0->$varname(1)$(n0:n0+nyr*12-1,:,:), \
          f0->$varname(2)$(n0+12*98:n0+12*98+nyr*12-1,:,:) - \
          f0->$varname(2)$(n0:n0+nyr*12-1,:,:) /)

 if (ie .eq. 0) then
   vara = new((/ne,nm,nyr,nz,ny/),float)
   veca = new((/ne,2,nm,nyr,nz,ny/),float)
 end if
 do im=0, nm-1
   vara(ie,im,:,:,:) = tmp(ispan(mm(im)-1,nyr*12-1,12),:,:)
   veca(ie,:,im,:,:,:) = epf(:,ispan(mm(im)-1,nyr*12-1,12),:,:)
 end do

 delete(tmp)
 delete(epf)

end do  ; ie

;--------------------------------------- process 0
 var0 = dim_avg_n(dim_avg_n(vara,2),1)
 vard = var0(0,:,:) - var0(1,:,:)

 var = array_append_record( var0, \
          conform_dims((/1,nz,ny/),vard,(/1,2/)), 0 )
 var = var;*1.e3  ; [mm/s]

 vec0 = dim_avg_n(dim_avg_n(veca,3),2)
 vecd = vec0(0,:,:,:) - vec0(1,:,:,:)

 vec = array_append_record( vec0, \
          conform_dims((/1,2,nz,ny/),vecd,(/1,2,3/)), 0 )
 vec = vec*1.e3  ; [mm/s]
 do k=19, nz-1
   vec(:,:,k,:) = vec(:,:,k,:)*exp(-z(k)/(7.*2.))
 end do

 veci = vec(:,:,24+6:52:2,::6)
 lati = lat(::6)
 zi = z(24+6:52:2)

;--------------------------------------- cn level
 scl = 0.02
 cnlev = ispan(-10,10,1)*scl
 cnlev2 = cnlev*10.
; cnlev2 = array_append_record( ispan(-25,-10,1), ispan(10,25,1), 0 )*scl
 cnlevd = cnlev*0.5
 cnlevd2 = cnlev2*0.5

; cnfil = (/2,3,4,5,6,7,8,9,10,11,0,0,26,25,24,23,22,21,20,19,18,17/)
 cnfil = (/17,18,19,20,21,22,23,24,25,26,0,0,11,10,9,8,7,6,5,4,3,2/)

;--------------------------------------- contour
 wi = 32
 he = 40*(25./50.)
 xf = 10 + (/0,1,2/)*40
 yf = 10

 rxy2x = (22./(180.*111.))/(tofloat(he)/wi)
 veci(:,0,:,:) = veci(:,0,:,:)*rxy2x

 ystr  = new(ne+1,string)

 mstr = (/"CC","CTL","CC - CTL"/)
 xstr = "latitude (deg)"
 ystr(:) = ""
 ystr(0) = "hybrid height (km)"
 lbopt = (/"(m s:S:-1:N:)","","(m s:S:-1:N:)"/)

do ie=0, ne-1

 contour1 = cnshadeplot(wks,lat,z,var(ie,:,:),cnlev,cnfil,lbopt(ie))
 cnaddline(contour1,wks,lat,z,var(ie,:,:),cnlev2,0,"D",0,"T",0.5,"C",(/1,-1,1/))
 setvalues contour1
  "pmLabelBarOrthogonalPosF" : 0.08
  "pmLabelBarParallelPosF"   : 1.15
  "pmLabelBarKeepAspect"     : False
  "pmLabelBarWidthF"         : 1.4
  "pmLabelBarHeightF"        : 0.12
  "lbLabelFontHeightF" : 0.026
  "lbTitlePosition"    : "Right"
  "lbTitleDirection"   : "Across"
  "lbTitleOffsetF"     : 0.04
  "lbTitleFontHeightF" : 0.032
 end setvalues

 vf = create "epf_vec" vectorFieldClass defaultapp
  "vfDataArray" : veci(ie,:,:,:)
  "vfXArray"    : lati
  "vfYArray"    : zi
  "vfMissingUValueV" : 1.e32
  "vfMissingVValueV" : 1.e32
 end create

 vec1 = create "vector" vectorPlotClass wks
  "vcVectorFieldData" : vf
  "vcMapDirection"  : False
  "vcRefMagnitudeF" : 0.5*0.3
  "vcRefLengthF"    : 0.2
  "vcRefAnnoOn"     : False
;  "vcLineArrowThicknessF" :  thc(nn)
;  "vcLineArrowColor"      :  col(nn)
;  "vcRefAnnoString1" : ref*tmpc(nn)
  "vcRefAnnoString2On" : False
;  "vcRefAnnoParallelPosF" :  1. - nn*0.1
  "vcGlyphStyle" : "CurlyVector"
  "vcFillArrowEdgeThicknessF" : 0.1
  "vcFillArrowWidthF" : 0.1
  "vcPositionMode" : "ArrowTail"
  "vcLineArrowHeadMaxSizeF" : 0.03
  "vcLineArrowHeadMinSizeF" : 0.002
 end create

 set_size(contour1,xf(ie),yf,wi,he)
 axis_range(contour1,-90,90,18,40,False,False)
 axis_str(contour1,mstr(ie),xstr,ystr(ie))

 set_size(vec1,xf(ie),yf,wi,he)
 axis_range(vec1,-90,90,18,40,False,False)
 axis_str(vec1,"","","")

 cnfinalize(contour1,wks)

 overlay(contour1,vec1)
 draw(contour1)
 NhlRemoveOverlay(contour1,vec1,False)

end do  ; ie

 ie = ne

 contour2 = cnshadeplot(wks,lat,z,var(ie,:,:),cnlevd,cnfil,lbopt(ie))
 cnaddline(contour2,wks,lat,z,var(ie,:,:),cnlevd2,0,"D",0,"T",0.5,"C",(/1,-1,1/))
 setvalues contour2
  "pmLabelBarOrthogonalPosF" : 0.08
  "pmLabelBarKeepAspect"     : False
  "pmLabelBarWidthF"         : 0.6
  "pmLabelBarHeightF"        : 0.12
  "lbLabelFontHeightF" : 0.026
  "lbTitlePosition"    : "Right"
  "lbTitleDirection"   : "Across"
  "lbTitleOffsetF"     : 0.04
  "lbTitleFontHeightF" : 0.032
 end setvalues

 vf = create "epf_vec" vectorFieldClass defaultapp
  "vfDataArray" : veci(ie,:,:,:)
  "vfXArray"    : lati
  "vfYArray"    : zi
  "vfMissingUValueV" : 1.e32
  "vfMissingVValueV" : 1.e32
 end create

 vec2 = create "vector" vectorPlotClass wks
  "vcVectorFieldData" : vf
  "vcMapDirection"  : False
  "vcRefMagnitudeF" : 0.5 * 0.15
  "vcRefLengthF"    : 0.2
  "vcRefAnnoOn"     : False
;  "vcLineArrowThicknessF" :  thc(nn)
;  "vcLineArrowColor"      :  col(nn)
;  "vcRefAnnoString1" : ref*tmpc(nn)
  "vcRefAnnoString2On" : False
;  "vcRefAnnoParallelPosF" :  1. - nn*0.1
  "vcGlyphStyle" : "CurlyVector"
  "vcFillArrowEdgeThicknessF" : 0.1
  "vcFillArrowWidthF" : 0.1
  "vcPositionMode" : "ArrowTail"
  "vcLineArrowHeadMaxSizeF" : 0.03
  "vcLineArrowHeadMinSizeF" : 0.002
 end create

 set_size(contour2,xf(ie),yf,wi,he)
 axis_range(contour2,-90,90,18,40,False,False)
 axis_str(contour2,mstr(ie),xstr,ystr(ie))

 set_size(vec2,xf(ie),yf,wi,he)
 axis_range(vec2,-90,90,18,40,False,False)
 axis_str(vec2,"","","")

 cnfinalize(contour2,wks)

 overlay(contour2,vec2)
 draw(contour2)
 NhlRemoveOverlay(contour2,vec2,False)


 frame(wks)

end

