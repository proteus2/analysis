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
 varname = (/"u","f_y","f_z"/)

 y0 = 1951
 y1 = 1960+98
 y2 = 1999+98
 mms = "djf"

 wks = cr_wks("fig/irefr_yz-"+mms+"_"+y1+"-"+y2+"-"+jobname(1),600,0)
; gsn_define_colormap(wks,"StepSeq25")
 gsn_merge_colormaps(wks,"StepSeq25","prcp_3")
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

 fnamej = "../dat/"+expname(ie)+"/xm.uan"+jobname(ie)+"a.pj.nc"
 fname = "../dat/"+expname(ie)+"/tem/uan"+jobname(ie)+".tem-wc_yzt.nc"
 fj = addfile(fnamej,"r")
 f0 = addfile(fname,"r")
 n0 = (y1 - y0)*12

 ; get y, z
 tmp = fj->$varname(0)$(n0:n0+nyr*12-1,:,:)
 lat = tmp&$tmp!2$
 z = tmp&$tmp!1$ / 1.e3
 ny = dimsizes(lat)
 nz = dimsizes(z)

 rho = fj->rho(n0:n0+nyr*12-1,:,:)
 pt = fj->theta(n0:n0+nyr*12-1,:,:)

 epf = (/ f0->$varname(1)$(n0:n0+nyr*12-1,:,:), \
          f0->$varname(2)$(n0:n0+nyr*12-1,:,:) /)

 if (ie .eq. 0) then
   vara = new((/ne,nm,nyr,nz,ny/),float)
   varr = new((/ne,nm,nyr,nz,ny/),float)
   varn = new((/ne,nm,nyr,nz,ny/),float)
   veca = new((/ne,2,nm,nyr,nz,ny/),float)
 end if
 do im=0, nm-1
   vara(ie,im,:,:,:) = tmp(ispan(mm(im)-1,nyr*12-1,12),:,:)
   varr(ie,im,:,:,:) = rho(ispan(mm(im)-1,nyr*12-1,12),:,:)
   varn(ie,im,:,:,:) = pt(ispan(mm(im)-1,nyr*12-1,12),:,:)
   veca(ie,:,im,:,:,:) = epf(:,ispan(mm(im)-1,nyr*12-1,12),:,:)
 end do

 delete(tmp)
 delete(rho)
 delete(pt)
 delete(epf)

end do  ; ie

;--------------------------------------- process 0
 tmp = log(varn)
 do k=1, nz-2
   varn(:,:,:,k,:) = 9.81*(tmp(:,:,:,k+1,:) - tmp(:,:,:,k-1,:)) \
                     /(z(k+1) - z(k-1))/1.e3
 end do
 varn(:,:,:,0,:) = 999.
 varn(:,:,:,nz-1,:) = 999.
 delete(tmp)

 var0 = dim_avg_n(dim_avg_n(vara,2),1)
 vard = var0(0,:,:) - var0(1,:,:)

 varr0 = dim_avg_n(dim_avg_n(varr,2),1)
 varn0 = dim_avg_n(dim_avg_n(varn,2),1)

 var = array_append_record( var0, \
          conform_dims((/1,nz,ny/),vard,(/1,2/)), 0 )

 vec0 = dim_avg_n(dim_avg_n(veca,3),2)
 vecd = vec0(0,:,:,:) - vec0(1,:,:,:)

 vec = array_append_record( vec0, \
          conform_dims((/1,2,nz,ny/),vecd,(/1,2,3/)), 0 )
 vec = vec/6371.*1.e3  ; [mPa]
 do k=19, nz-1
   vec(:,:,k,:) = vec(:,:,k,:)*exp(z(k)/14.)
 end do

 veci = vec(:,:,21::2,::4)
 lati = lat(::4)
 zi = z(21::2)

 s_w = 2.
 lat_w = 60.
 c_w = 0.1

 a_e = 6370.e3
 k_w = (2./a_e/cos(lat_w*3.141592/180.))
 h_s = 7.0e3
 cosphi = cos(lat*3.141592/180.)
 f = 2.*7.2921e-5*sin(lat*3.141592/180.)
 d2phi = (lat(2) - lat(0))*3.141592/180.
 d2y = d2phi*a_e
 nbv = sqrt(varn0)
 q_phi = new(dimsizes(var0),float)
 tmp = new(dimsizes(var0),float)
 do j=1, ny-2
   tmp(:,:,j) = (var0(:,:,j+1)*cosphi(j+1) - var0(:,:,j-1)*cosphi(j-1))/d2y \
                /cosphi(j)
 end do
 do j=2, ny-3
   q_phi(:,:,j) = 2.*7.2921e-5*cosphi(j) - \
                  (tmp(:,:,j+1) - tmp(:,:,j-1))/d2phi
 end do
 do k=1, nz-2
   tmp(:,k,:) = (var0(:,k+1,:) - var0(:,k-1,:))/(z(k+1) - z(k-1))/1.e3 \
                *varr0(:,k,:)/varn0(:,k,:)
 end do
 do k=2, nz-3
   q_phi(:,k,:) = q_phi(:,k,:) - conform_dims((/ne,ny/),f^2,1)/varr0(:,k,:)* \
                  (tmp(:,k+1,:) - tmp(:,k-1,:))/(z(k+1) - z(k-1))/1.e3
 end do

 ir = a_e*q_phi/(var0 - c_w) - k_w^2 - \
      (0.5*conform_dims((/ne,nz,ny/),f,2)*a_e/nbv/h_s)^2

 varird = where( ir(0,:,:) .gt. 0. .and. ir(1,:,:) .gt. 0., \
          (ir(0,:,:) - ir(1,:,:))/exp(0.5*(log(ir(0,:,:))+log(ir(1,:,:)))), \
          -999.)*100.
 
 varir = array_append_record( ir, \
          conform_dims((/1,nz,ny/),varird,(/1,2/)), 0 )

 varir@_FillValue = 1.e20

 varir(0:1,:,:) = where(varir(0:1,:,:) .lt. 0., varir@_FillValue,varir(0:1,:,:))
 varir(2,:,:) = where(varir(2,:,:) .eq. -99900., varir@_FillValue,varir(2,:,:))

 varir(:,:,0:1) = varir@_FillValue
 varir(:,:,ny-2:ny-1) = varir@_FillValue
 varir(:,0:1,:) = varir@_FillValue
 varir(:,nz-2:nz-1,:) = varir@_FillValue
 
;--------------------------------------- cn level
 scl = 0.5
 cnlev = ispan(-10,10,1)*scl
 cnlev2 = cnlev*10.
; cnlev2 = array_append_record( ispan(-25,-10,1), ispan(10,25,1), 0 )*scl
 cnlevd = cnlev*0.2
 cnlevd2 = cnlev2*0.2

 delete(cnlev)
 cnlev = ispan(5,15,1)*0.2
 delete(cnlevd)
 cnlevd = ispan(-10,10,1)*2.

; cnfil = (/17,18,19,20,21,22,23,24,25,26,0,0,11,10,9,8,7,6,5,4,3,2/)
 cnfil = (/-25,5,6,7,10,12,14,15,16,17,19,21/)+25
 cnfild = (/17,18,19,20,21,22,23,24,25,26,0,0,11,10,9,8,7,6,5,4,3,2/)

;--------------------------------------- contour
 wi = 32
 he = 40*(32./50.)
 xf = 10 + (/0,1,2/)*40
 yf = 10

 rxy2x = (32./(180.*111.))/(tofloat(he)/wi)
 veci(:,0,:,:) = veci(:,0,:,:)*rxy2x

 ystr  = new(ne+1,string)

 mstr = (/"CC","CTL","CC - CTL"/)
 xstr = "latitude (deg)"
 ystr(:) = ""
 ystr(0) = "hybrid height (km)"
 lbopt = (/" ",""," "/)

do ie=0, ne-1

 contour1 = cnshadeplot(wks,lat,z,log10(varir(ie,:,:)),cnlev,cnfil,lbopt(ie))
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
  "vcRefMagnitudeF" : 1.e6 * 0.5
  "vcRefLengthF"    : 0.2
  "vcRefAnnoOn"     : False
;  "vcLineArrowThicknessF" :  thc(nn)
;  "vcLineArrowColor"      :  col(nn)
;  "vcRefAnnoString1" : ref*tmpc(nn)
  "vcRefAnnoString2On" : False
;  "vcRefAnnoParallelPosF" :  1. - nn*0.1
;  "vcGlyphStyle" : "FillArrow"
  "vcFillArrowEdgeThicknessF" : 0.1
  "vcFillArrowWidthF" : 0.1
  "vcPositionMode" : "ArrowTail"
 end create

 set_size(contour1,xf(ie),yf,wi,he)
 axis_range(contour1,0,90,8,40,False,False)
 axis_str(contour1,mstr(ie),xstr,ystr(ie))

 set_size(vec1,xf(ie),yf,wi,he)
 axis_range(vec1,0,90,8,40,False,False)
 axis_str(vec1,"","","")

 cnfinalize(contour1,wks)

 overlay(contour1,vec1)
 draw(contour1)
 NhlRemoveOverlay(contour1,vec1,False)

end do  ; ie

 ie = ne

 contour2 = cnshadeplot(wks,lat,z,varir(ie,:,:),cnlevd,cnfild,lbopt(ie))
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
  "vcRefMagnitudeF" : 1.e6 * 0.02
  "vcRefLengthF"    : 0.2
  "vcRefAnnoOn"     : False
;  "vcLineArrowThicknessF" :  thc(nn)
;  "vcLineArrowColor"      :  col(nn)
;  "vcRefAnnoString1" : ref*tmpc(nn)
  "vcRefAnnoString2On" : False
;  "vcRefAnnoParallelPosF" :  1. - nn*0.1
;  "vcGlyphStyle" : "FillArrow"
  "vcFillArrowEdgeThicknessF" : 0.1
  "vcFillArrowWidthF" : 0.1
  "vcPositionMode" : "ArrowTail"
 end create

 set_size(contour2,xf(ie),yf,wi,he)
 axis_range(contour2,0,90,8,40,False,False)
 axis_str(contour2,mstr(ie),xstr,ystr(ie))

 set_size(vec2,xf(ie),yf,wi,he)
 axis_range(vec2,0,90,8,40,False,False)
 axis_str(vec2,"","","")

 cnfinalize(contour2,wks)

 overlay(contour2,vec2)
 draw(contour2)
 NhlRemoveOverlay(contour2,vec2,False)


 frame(wks)

end

