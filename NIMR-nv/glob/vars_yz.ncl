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
 ppid = (/"j","e"/)
 varname = (/"u","du_gw","epd","f_y","f_z"/)
 nv = 3
 iepf = nv + (/0,1/)

 y0 = 1951
 y1 = 1951
 y2 = 1999;2009
 mms = "jja"

 wks = cr_wks("fig/vars_yz-"+mms+"_"+y1+"-"+y2+"-"+jobname(1),600,0)
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

 fname0 = "../dat/"+expname(ie)+"/xm.uan"+jobname(ie)+"a.p"+ppid(0)+".nc"
 fname1 = "../dat/"+expname(ie)+"/xm.uan"+jobname(ie)+"a.p"+ppid(1)+".nc"
 fname9 = "../dat/"+expname(ie)+"/tem/uan"+jobname(ie)+".tem_yzt.nc"
 f0 = addfile(fname0,"r")
 f1 = addfile(fname1,"r")
 f9 = addfile(fname9,"r")
 n0 = (y1 - y0)*12

 ; get y, z
 tmp0 = f0->$varname(0)$(n0:n0+nyr*12-1,:,:)
 tmp1 = f1->$varname(1)$(n0:n0+nyr*12-1,:,:)
 tmp2 = f9->$varname(2)$(n0:n0+nyr*12-1,:,:)
 lat = tmp0&$tmp0!2$
 z = tmp0&$tmp0!1$ / 1.e3
 ny = dimsizes(lat)
 nz = dimsizes(z)

 tmp1 = tmp1*72.

 if ( isfilevar(f0,"dudt_gwconv") ) then
   tmp1 = tmp1 + f1->dudt_gwconv(n0:n0+nyr*12-1,:,:)*86400.
 end if

 epf = (/ f9->$varname(iepf(0))$(n0:n0+nyr*12-1,:,:), \
          f9->$varname(iepf(1))$(n0:n0+nyr*12-1,:,:) /)

 if (ie .eq. 0) then
   vara = new((/ne,nv,nm,nyr,nz,ny/),float)
   veca = new((/ne,2,nm,nyr,nz,ny/),float)
 end if
 do im=0, nm-1
   vara(ie,0,im,:,:,:) = (/ tmp0(ispan(mm(im)-1,nyr*12-1,12),:,:) /)
   vara(ie,1,im,:,:,:) = (/ tmp1(ispan(mm(im)-1,nyr*12-1,12),:,:) /)
   vara(ie,2,im,:,:,:) = (/ tmp2(ispan(mm(im)-1,nyr*12-1,12),:,:) /)
   veca(ie,:,im,:,:,:) = epf(:,ispan(mm(im)-1,nyr*12-1,12),:,:)
 end do

 delete(tmp0)
 delete(tmp1)
 delete(tmp2)
 delete(epf)

end do  ; ie

;--------------------------------------- process 0
 var0 = dim_avg_n(dim_avg_n(vara,3),2)
 vard = var0(1,:,:,:) - var0(0,:,:,:)

 var = array_append_record( var0, \
          conform_dims((/1,nv,nz,ny/),vard,(/1,2,3/)), 0 )

 vec0 = dim_avg_n(dim_avg_n(veca,3),2)
 vecd = vec0(1,:,:,:) - vec0(0,:,:,:)

 vec = array_append_record( vec0, \
          conform_dims((/1,2,nz,ny/),vecd,(/1,2,3/)), 0 )
 vec = vec/6371.*1.e3  ; [mPa]
 do k=19, nz-1
   vec(:,:,k,:) = vec(:,:,k,:)*exp(z(k)/(7.))
 end do
 vec(:,0,:,:) = vec(:,0,:,:)/111.

 veci = vec(:,:,24::2,::8)
 lati = lat(::8)
 zi = z(24::2)

;--------------------------------------- cn level
; cnfil = (/2,3,4,5,6,7,8,9,10,11,0,0,26,25,24,23,22,21,20,19,18,17/)
 cnfil = (/17,18,19,20,21,22,23,24,25,26,0,0,11,10,9,8,7,6,5,4,3,2/)

;--------------------------------------- contour
 wi = 32
 he = 32
 xf = 10 + (/0,1,2/)*(wi+8)
 yf = 10 + (/0,1,1/)*(he+16)

 ystr  = new(ne+1,string)

 mstr = (/"U","GWD","EPD"/)
 xstr = (/"latitude(deg)","","latitude (deg)"/)
 ystr(:) = ""
 ystr(0) = "hybrid height (km)"
 lbopt = (/"(m s:S:-1:N:)","(m s:S:-1:N: day:S:-1:N:)","(m s:S:-1:N: day:S:-1:N:)"/)

do iv=0, nv-1

 if (iv .eq. 0) then
 scl = 5.
 cnlev = ispan(-10,10,1)*scl
; cnlev2 = cnlev*10.
 cnlev2 = array_append_record( ispan(-25,-10,1), ispan(10,25,1), 0 )*scl
 cnlevd = cnlev*0.1
 cnlevd2 = cnlev2*0.1
 end if
 if (iv .ge. 1) then
 scl = 0.3
 cnlev = ispan(-10,10,1)*scl
 cnlev2 = cnlev*10.
; cnlev2 = array_append_record( ispan(-25,-10,1), ispan(10,25,1), 0 )*scl
 cnlevd = cnlev*0.2
 cnlevd2 = cnlev2*0.2
 end if

do ie=0, 0

 contour1 = cnshadeplot(wks,lat,z,var(ie,iv,:,:),cnlev,cnfil,lbopt(iv))
 cnaddline(contour1,wks,lat,z,var(ie,iv,:,:),cnlev2,0,"D",0,"T",0.5,"C",(/1,-1,1/))
 setvalues contour1
;  "pmLabelBarOrthogonalPosF" : -0.04 ;0.08
;  "pmLabelBarParallelPosF"   : 1.15
  "pmLabelBarKeepAspect"     : False
  "pmLabelBarWidthF"         : 0.55 ;1.4
  "pmLabelBarHeightF"        : 0.12
  "lbLabelFontHeightF" : 0.026
  "lbTitlePosition"    : "Right"
  "lbTitleDirection"   : "Across"
  "lbTitleOffsetF"     : 0.04
  "lbTitleFontHeightF" : 0.032
 end setvalues

 set_size(contour1,xf(iv),yf(ie),wi,he)
 axis_range(contour1,-90,90,0,50,False,False)
 axis_str(contour1,mstr(iv),xstr(ie),ystr(iv))

 cnfinalize(contour1,wks)

 if (iv .eq. nv-1) then

 vf = create "epf_vec" vectorFieldClass defaultapp
  "vfDataArray" : veci(ie,:,:,:)
  "vfXArray"    : lati
  "vfYArray"    : zi
  "vfMissingUValueV" : 1.e32
  "vfMissingVValueV" : 1.e32
 end create

 vec1 = create "vector" vectorPlotClass wks
  "vcVectorFieldData" : vf
  "vcRefMagnitudeF" : 1.e6
  "vcRefLengthF"    : 0.05
  "vcRefAnnoOn"     : False
;  "vcLineArrowThicknessF" :  thc(nn)
;  "vcLineArrowColor"      :  col(nn)
;  "vcRefAnnoString1" : ref*tmpc(nn)
  "vcRefAnnoString2On" : False
;  "vcRefAnnoParallelPosF" :  1. - nn*0.1
  "vcGlyphStyle" : "FillArrow"
  "vcFillArrowEdgeThicknessF" : 0.1
  "vcFillArrowWidthF" : 0.1
  "vcPositionMode" : "ArrowTail"
 end create

 set_size(vec1,xf(ie),yf,wi,he)
 axis_range(vec1,-90,90,0,50,False,False)
 axis_str(vec1,"","","")

 overlay(contour1,vec1)

 end if

 draw(contour1)

 if (iv .eq. nv-1) then
   NhlRemoveOverlay(contour1,vec1,False)
 end if

end do  ; ie

 ie = ne

 contour2 = cnshadeplot(wks,lat,z,var(ie,iv,:,:),cnlevd,cnfil,lbopt(iv))
 cnaddline(contour2,wks,lat,z,var(ie,iv,:,:),cnlevd2,0,"D",0,"T",0.5,"C",(/1,-1,1/))
 setvalues contour2
;  "pmLabelBarOrthogonalPosF" : 0.08
  "pmLabelBarKeepAspect"     : False
  "pmLabelBarWidthF"         : 0.55
  "pmLabelBarHeightF"        : 0.12
  "lbLabelFontHeightF" : 0.026
  "lbTitlePosition"    : "Right"
  "lbTitleDirection"   : "Across"
  "lbTitleOffsetF"     : 0.04
  "lbTitleFontHeightF" : 0.032
 end setvalues

 set_size(contour2,xf(iv),yf(ie),wi,he)
 axis_range(contour2,-90,90,0,50,False,False)
 axis_str(contour2,"",xstr(ie),ystr(iv))

 cnfinalize(contour2,wks)

 if (iv .eq. nv-1) then

 vf = create "epf_vec" vectorFieldClass defaultapp
  "vfDataArray" : veci(ie,:,:,:)
  "vfXArray"    : lati
  "vfYArray"    : zi
  "vfMissingUValueV" : 1.e32
  "vfMissingVValueV" : 1.e32
 end create

 vec2 = create "vector" vectorPlotClass wks
  "vcVectorFieldData" : vf
  "vcRefMagnitudeF" : 1.e6 * 0.1
  "vcRefLengthF"    : 0.05
  "vcRefAnnoOn"     : False
;  "vcLineArrowThicknessF" :  thc(nn)
;  "vcLineArrowColor"      :  col(nn)
;  "vcRefAnnoString1" : ref*tmpc(nn)
  "vcRefAnnoString2On" : False
;  "vcRefAnnoParallelPosF" :  1. - nn*0.1
  "vcGlyphStyle" : "FillArrow"
  "vcFillArrowEdgeThicknessF" : 0.1
  "vcFillArrowWidthF" : 0.1
  "vcPositionMode" : "ArrowTail"
 end create

 set_size(vec2,xf(iv),yf(ie),wi,he)
 axis_range(vec2,-90,90,0,50,False,False)
 axis_str(vec2,"","","")

 overlay(contour2,vec2)

 end if

 draw(contour2)

 if (iv .eq. nv-1) then
   NhlRemoveOverlay(contour2,vec2,False)
 end if

 delete(cnlev2)
 delete(cnlevd2)

end do  ; iv


 frame(wks)

end

