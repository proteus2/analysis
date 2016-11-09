load "$UTIL/ncl/header.ncl"

xspace =  83  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 1.0

begin

;--------------------------------------- parameters
 jobname = (/"tk","tc"/)
; jobname = (/"ui","uj"/)
 expname = "L60"+(/"CTL","CGW"/)
; expname = "L60"+(/"CTL","CGW-t"/)
 year1 = 1960;+2
 year2 = 1999;+27

 wks = cr_wks("fig/eawm_u300cli2-"+jobname(1),600,0)
 gsn_define_colormap(wks,"StepSeq25")
 colind = NhlNewColor(wks,0.6,0.6,0.6)

;--------------------------------------- data reading
 nt = year2 - year1 + 1
 time = fspan(year1,year2,nt)

 eawmi = new((/2,nt/),float)
 eawmi = 0.

 mon3 = (/1,2,12/)

 n = 0
 do iy=year1, year2
   do im=0, 2
;     file_uv = "/hippo0/HG2AMIP/"+expname+"/pj/uan"+jobname+"a.pj_"+ \
     file_uv = "/hippo0/HG2CMIP/"+expname+"/pj/uan"+jobname+"a.pj_"+ \
               iy+sprinti("%0.2i",mon3(im))+".nc"
     file_uvf = "/hippo0/HG2CMIP/"+expname+"/pj/uan"+jobname+"a.pj_"+ \
               (iy+98)+sprinti("%0.2i",mon3(im))+".nc"
     fu0 = addfile(file_uv(0),"r")
     fu1 = addfile(file_uv(1),"r")
     fu0f = addfile(file_uvf(0),"r")
     fu1f = addfile(file_uvf(1),"r")
     if ( .not. isdefined("um") ) then
       um = (/fu0->u_p_uv(0,{300},:,:),fu1->u_p_uv(0,{300},:,:)/)
       vm = (/fu0->v_p_uv(0,{300},:,:),fu1->v_p_uv(0,{300},:,:)/)
       um = 0.
       vm = 0.
       lon = fu0->u&longitude_1
       lat = fu0->u&latitude_1
       ua = new((/2,dimsizes(lat),dimsizes(lon),nt/),float)
       ua = 0.
     end if
     ; calculate EAWMI
     u = (/fu0f->u_p_uv(0,{300},:,:),fu1f->u_p_uv(0,{300},:,:)/) - \
         (/fu0->u_p_uv(0,{300},:,:),fu1->u_p_uv(0,{300},:,:)/)
     v = (/fu0f->v_p_uv(0,{300},:,:),fu1f->v_p_uv(0,{300},:,:)/) - \
         (/fu0->v_p_uv(0,{300},:,:),fu1->v_p_uv(0,{300},:,:)/)
     u!1 = "lat"
     u!2 = "lon"
     u&lat = lat
     u&lon = lon
     do ie=0, 1
       eawmi(ie,n) = eawmi(ie,n) + ( avg(u(ie,{27.5:37.5},{110:170})) - \
                                     avg(u(ie,{50.0:60.0},{ 80:140})) )/3.
     end do
     um = um + u/(tofloat(nt*3))
     vm = vm + v/(tofloat(nt*3))
     ua(:,:,:,n) = ua(:,:,:,n) + u(0:1,:,:)/3.
   end do
   n = n + 1
 end do
print(dim_avg(eawmi))

; eawmi = eawmi - avg(eawmi)
 eawmi = eawmi - 41.50

;--------------------------------------- process 0
 ud = um(1,:,:) - um(0,:,:)
 vd = vm(1,:,:) - vm(0,:,:)

 um = where(um .eq. 0., 1.e15, um)
 vm = where(vm .eq. 0., 1.e15, vm)
 ud = where(ud .eq. 0., 1.e15, ud)
 vd = where(vd .eq. 0., 1.e15, vd)

 wnd = (/um(:,::3,::3),vm(:,::3,::3)/)
 wind = new(2,graphic)
 do ie=0, 1
 wind(ie) = create "f1" vectorFieldClass defaultapp
  "vfDataArray" : wnd(:,ie,:,:)
  "vfXArray" : lon(::3)
  "vfYArray" : lat(::3)
  "vfMissingUValueV" : 1.e15
  "vfMissingVValueV" : 1.e15
 end create
 end do

 wndd = (/ud(::3,::3),vd(::3,::3)/)
 windd = create "f1" vectorFieldClass defaultapp
  "vfDataArray" : wndd
  "vfXArray" : lon(::3)
  "vfYArray" : lat(::3)
  "vfMissingUValueV" : 1.e15
  "vfMissingVValueV" : 1.e15
 end create

 um@_FillValue = 1.e15
 vm@_FillValue = 1.e15
 ud@_FillValue = 1.e15
 vd@_FillValue = 1.e15

 uv = sqrt(um^2 + vm^2)
 uv@_FillValue = 1.e15

 ua = where(ua .eq. 0., 1.e15, ua)
 var0 = dim_variance_n(ua(0,:,:,:),2)
 var1 = dim_variance_n(ua(1,:,:,:),2)
 var0 = where(var0 .eq. 0., 1.e15, var0)
 var1 = where(var1 .eq. 0., 1.e15, var1)

 tt = (um(1,:,:) - um(0,:,:))/sqrt( (var0+var1)/tofloat(nt) )

 tt = student_t(tt,nt-1)

;--------------------------------------- cn level
; cnlev = fspan(4,72,18)
 cnlevd1 = fspan(-10,10,21);*0.3
 cnfil = (/0,11,10,9,6,5,4,3,2,26,25,24,23,22,21,20,19,18,17/)
 cnlevd = fspan(-10,10,21)*0.3
 cnfild = (/17,18,19,20,21,22,23,24,25,26,0,0,11,10,9,8,7,6,5,4,3,2/)

;--------------------------------------- drawing
 wi = (/nt*2,36/)
 he = (/20,30/)
 xf = (/10,10,50,90/)
 yf = (/10,42/)

; plot = xylineplot(wks,time,eawmi,"D",0,"T",1.2,"C",(/24,4/))
; xyaddmarker(plot,16,"S",0.5,"T",1.0,"C",1)
; xyaddline(plot,(/year1-6,year2+5/),0.,"D",0,"T",0.6,"C",colind)
; set_size(plot,xf(0),yf(0),wi(0),he(0))
; axis_range(plot,year1-1.0,year2+0.5,-999,999,False,False)
; setvalues plot
;  "tmXBMode" : "Manual"
;  "tmXBTickStartF"    : 1950
;  "tmXBTickSpacingF"  : 5
;  "tmXBMinorPerMajor" : 4
; end setvalues
; axis_str(plot,"East Asian Winter Monsoon Index","year","(m s:S:-1:N:)")
; draw(plot)

 map = create "mp" mapPlotClass wks
  "mpProjection"     : "Mollweide" ;"Robinson"
  "mpCenterLonF"     : 125.
  "mpOutlineOn"      : True
  "mpGridAndLimbOn"  : False
  "mpPerimOn"        : True
  "mpPerimDrawOrder" : "PostDraw"

  "mpShapeMode"     : "FixedAspectFitBB"
  "mpGeophysicalLineColor" : "Orange4"
  "mpDataBaseVersion" : "Ncarg4_0"

  "mpLimitMode"    : "LatLon"
  "mpMinLonF"      : 85.
  "mpMaxLonF"      : 165.
  "mpMinLatF"      : 10.
  "mpMaxLatF"      : 60.
 end create
 set_size_mp(map,xf(1),yf(1),wi(1),he(1))

 vec1 = create "vector" vectorPlotClass wks
  "vcVectorFieldData" : wind(0)
  "vcRefMagnitudeF" : 2;30.
  "vcRefLengthF"    : 0.02
  "vcRefAnnoOn"     : False
 end create
 set_size(vec1,xf(1),yf(1),wi(1),he(1))
; axis_str(vec1,"DJF climatology of 300-hPa wind","","")

 contour1 = cnshadeplot(wks,lon,lat,um(0,:,:),cnlevd1,cnfild,"u (m s:S:-1:N:)")
 setvalues contour1
  "cnMissingValFillColor" : colind
  "lbTitlePosition"    : "Right"
  "lbTitleDirection"   : "Across"
  "lbTitleOffsetF"     : 0.06
 end setvalues
 set_size(contour1,xf(1),yf(1),wi(1),he(1))
 axis_str(contour1,"CTL","","")

 resp = True
 resp@gsLineThicknessF = 1.2
; dum1 = gsn_add_polyline(wks,contour1,(/80,80,140,140,80/),(/50,60,60,50,50/),resp)
; dum2 = gsn_add_polyline(wks,contour1,(/110,110,170,170,110/),(/27.5,37.5,37.5,27.5,27.5/),resp)

 overlay(map,vec1)
 NhlAddOverlay(map,contour1,map)

 draw(map)

 NhlRemoveOverlay(map,contour1,False)
 NhlRemoveOverlay(map,vec1    ,False)

 set_size_mp(map,xf(2),yf(1),wi(1),he(1))

 vec2 = create "vector" vectorPlotClass wks
  "vcVectorFieldData" : wind(1)
  "vcRefMagnitudeF" : 2;30.
  "vcRefLengthF"    : 0.02
  "vcRefAnnoOn"     : False
 end create
 set_size(vec2,xf(2),yf(1),wi(1),he(1))

 axis_str(vec2,"CC","","")

 contour2 = cnshadeplot(wks,lon,lat,um(1,:,:),cnlevd1,cnfild,"u (m s:S:-1:N:)")
 setvalues contour2
  "cnMissingValFillColor" : colind
  "lbTitlePosition"    : "Right"
  "lbTitleDirection"   : "Across"
  "lbTitleOffsetF"     : 0.06
 end setvalues
 set_size(contour2,xf(2),yf(1),wi(1),he(1))

; dum1 = gsn_add_polyline(wks,contour2,(/80,80,140,140,80/),(/50,60,60,50,50/),resp)
; dum2 = gsn_add_polyline(wks,contour2,(/110,110,170,170,110/),(/27.5,37.5,37.5,27.5,27.5/),resp)

 overlay(map,vec2)
 NhlAddOverlay(map,contour2,map)

 draw(map)

 NhlRemoveOverlay(map,contour2,False)
 NhlRemoveOverlay(map,vec2    ,False)

 set_size_mp(map,xf(3),yf(1),wi(1),he(1))

 vec3 = create "vector" vectorPlotClass wks
  "vcVectorFieldData" : windd
  "vcRefMagnitudeF" : 2.
  "vcRefLengthF"    : 0.02
  "vcRefAnnoOn"     : False
 end create
 set_size(vec3,xf(3),yf(1),wi(1),he(1))

 axis_str(vec3,"CC - CTL","","")

 contour3 = cnshadeplot(wks,lon,lat,ud,cnlevd,cnfild,"u (m s:S:-1:N:)")
 cnaddline(contour3,wks,lon,lat,tt,(/0.05/),0,"D",0,"T",1,"C",0)
 setvalues contour3
  "cnMissingValFillColor" : colind
  "lbTitlePosition"    : "Right"
  "lbTitleDirection"   : "Across"
  "lbTitleOffsetF"     : 0.06
 end setvalues
 set_size(contour3,xf(3),yf(1),wi(1),he(1))

; dum1 = gsn_add_polyline(wks,contour3,(/80,80,140,140,80/),(/50,60,60,50,50/),resp)
; dum2 = gsn_add_polyline(wks,contour3,(/110,110,170,170,110/),(/27.5,37.5,37.5,27.5,27.5/),resp)

 overlay(map,vec3)
 NhlAddOverlay(map,contour3,map)

 draw(map)

 frame(wks)

 status_exit(10)

end

