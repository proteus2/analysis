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
 year2 = 1999;+26

 wks = cr_wks("fig/easm_raincli2-"+jobname(1),600,0)
 gsn_merge_colormaps(wks,"precip_11lev","precip_diff_12lev")
 colind = NhlNewColor(wks,0.6,0.6,0.6)

;--------------------------------------- data reading
 nt = year2 - year1 + 1
 time = fspan(year1,year2,nt)

 easmi = new((/2,nt/),float)
 easmi = 0.

 n = 0
 do iy=year1, year2
   do im=6, 8
;     file_uv = "/hippo0/HG2AMIP/"+expname+"/pj/uan"+jobname+"a.pj_"+ \
     file_uv = "/hippo0/HG2CMIP/"+expname+"/pj/uan"+jobname+"a.pj_"+ \
               iy+sprinti("%0.2i",im)+".nc"
     file_uvf = "/hippo0/HG2CMIP/"+expname+"/pj/uan"+jobname+"a.pj_"+ \
               (iy+98)+sprinti("%0.2i",im)+".nc"
     fu0 = addfile(file_uv(0),"r")
     fu1 = addfile(file_uv(1),"r")
     fu0f = addfile(file_uvf(0),"r")
     fu1f = addfile(file_uvf(1),"r")
     if ( .not. isdefined("um") ) then
       um = (/fu0->precip(0,0,:,:),fu1->precip(0,0,:,:)/)
       um = 0.
       lon = fu0->u&longitude
       lat = fu0->u&latitude
       ua = new((/2,dimsizes(lat),dimsizes(lon),nt/),float)
       ua = 0.
     end if
     ; calculate EASMI
     u = (/fu0f->precip(0,0,:,:),fu1f->precip(0,0,:,:)/)*86400. - \
         (/fu0->precip(0,0,:,:),fu1->precip(0,0,:,:)/)*86400.
     u!1 = "lat"
     u!2 = "lon"
     u&lat = lat
     u&lon = lon
     do ie=0, 1
       easmi(ie,n) = easmi(ie,n) + ( avg(u(ie,{22.5:32.5},{110:140})) - \
                                     avg(u(ie,{ 5.0:15.0},{ 90:130})) )/3.
     end do
     um = um + u/(tofloat(nt*3))
     ua(:,:,:,n) = ua(:,:,:,n) + u(0:1,:,:)/3.
   end do
   n = n + 1
 end do
print(dim_avg(easmi))

; easmi = easmi - avg(easmi)
 easmi = easmi + 3.96

;--------------------------------------- process 0
 ud = um(1,:,:) - um(0,:,:)

 um = where(um .eq. 0., 1.e15, um)
 ud = where(ud .eq. 0., 1.e15, ud)

 um@_FillValue = 1.e15
 ud@_FillValue = 1.e15

 ua = where(ua .eq. 0., 1.e15, ua)
 var0 = dim_variance_n(ua(0,:,:,:),2)
 var1 = dim_variance_n(ua(1,:,:,:),2)
 var0 = where(var0 .eq. 0., 1.e15, var0)
 var1 = where(var1 .eq. 0., 1.e15, var1)

 tt = (um(1,:,:) - um(0,:,:))/sqrt( (var0+var1)/tofloat(nt) )

 tt = student_t(tt,nt-1)

;--------------------------------------- cn level
 cnlev = fspan(2,18,9)
 cnfil = (/0,4,5,6,8,9,10,11,12,13/)

 cnlevd1 = fspan(-6,6,13)*0.5
 cnlevd = fspan(-6,6,13)*0.2
 cnfild = (/14,15,16,17,18,19,20,20,21,22,23,24,25,26/)

;--------------------------------------- drawing
 wi = (/nt*2,36/)
 he = (/20,30/)
 xf = (/10,10,50,90/)
 yf = (/10,42/)

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
  "mpMinLatF"      : 0.
  "mpMaxLatF"      : 50.
 end create
 set_size_mp(map,xf(1),yf(1),wi(1),he(1))

 contour1 = cnshadeplot(wks,lon,lat,um(0,:,:),cnlevd1,cnfild,"(mm day:S:-1:N:)")
 setvalues contour1
  "cnMissingValFillColor" : colind
  "lbTitlePosition"    : "Right"
  "lbTitleDirection"   : "Across"
  "lbTitleOffsetF"     : 0.12
 end setvalues
 set_size(contour1,xf(1),yf(1),wi(1),he(1))
 axis_str(contour1,"CTL","","")

 resp = True
 resp@gsLineThicknessF = 1.2
; dum1 = gsn_add_polyline(wks,contour1,(/90,90,130,130,90/),(/5,15,15,5,5/),resp)
; dum2 = gsn_add_polyline(wks,contour1,(/110,110,140,140,110/),(/22.5,32.5,32.5,22.5,22.5/),resp)

 overlay(map,contour1)

 draw(map)

 NhlRemoveOverlay(map,contour1,False)

 set_size_mp(map,xf(2),yf(1),wi(1),he(1))

 contour2 = cnshadeplot(wks,lon,lat,um(1,:,:),cnlevd1,cnfild,"(mm day:S:-1:N:)")
 setvalues contour2
  "cnMissingValFillColor" : colind
  "lbTitlePosition"    : "Right"
  "lbTitleDirection"   : "Across"
  "lbTitleOffsetF"     : 0.12
 end setvalues
 set_size(contour2,xf(2),yf(1),wi(1),he(1))
 axis_str(contour2,"CC","","")

; dum1 = gsn_add_polyline(wks,contour2,(/90,90,130,130,90/),(/5,15,15,5,5/),resp)
; dum2 = gsn_add_polyline(wks,contour2,(/110,110,140,140,110/),(/22.5,32.5,32.5,22.5,22.5/),resp)

 overlay(map,contour2)

 draw(map)

 NhlRemoveOverlay(map,contour2,False)

 set_size_mp(map,xf(3),yf(1),wi(1),he(1))

 contour3 = cnshadeplot(wks,lon,lat,ud,cnlevd,cnfild,"(mm day:S:-1:N:)")
 cnaddline(contour3,wks,lon,lat,tt,(/0.05/),0,"D",0,"T",1,"C",1)
 setvalues contour3
  "cnMissingValFillColor" : colind
  "lbTitlePosition"    : "Right"
  "lbTitleDirection"   : "Across"
  "lbTitleOffsetF"     : 0.12
 end setvalues
 set_size(contour3,xf(3),yf(1),wi(1),he(1))
 axis_str(contour3,"CC - CTL","","")

; dum1 = gsn_add_polyline(wks,contour3,(/90,90,130,130,90/),(/5,15,15,5,5/),resp)
; dum2 = gsn_add_polyline(wks,contour3,(/110,110,140,140,110/),(/22.5,32.5,32.5,22.5,22.5/),resp)

 overlay(map,contour3)

 draw(map)

 frame(wks)

 status_exit(10)

end

