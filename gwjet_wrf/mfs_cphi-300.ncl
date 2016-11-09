load "$UTIL/ncl/header.ncl"

xspace = 170  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 0.75

begin

;--------------------------------------- parameters
 cmax = (/25,15/)
 ctck = 5
 fdir = "$DATD/gwjet_wrf"
 var1 = (/"wu","wv"/)
 lev = 8000
 var1name = var1+"_phi_c"
;	Time = UNLIMITED ; // (1152 currently)
 
 wks = cr_wks("fig/mfs_cphi-300",600,0)
; gsn_define_colormap(wks,"prcp_3")
 gsn_define_colormap(wks,"StepSeq25")
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- read data
 fname0 = fdir+"/x6_5min/mfs_z/prt_d300km/mfs_prt_phi_2d__z"+sprinti("%5.5i",lev)+"-"+ispan(1,2,1)+"a.nc"
 f01 = addfile(fname0(0),"r")
 f02 = addfile(fname0(1),"r")

 cmax_i = max(cmax)+5
 vs0 = (/ (/ f01->$var1name(0)$({:cmax_i},:), f01->$var1name(1)$({:cmax_i},:) /), \
          (/ f02->$var1name(0)$({:cmax_i},:), f02->$var1name(1)$({:cmax_i},:) /) /)
 dir0= f01->dir
 c   = f01->c({:cmax_i})

;--------------------------------------- process 0
 rho = (/0.465,0.500/)  ; at 8 km in D1 and D2
 vs0(0,:,:,:) = vs0(0,:,:,:)*rho(0) * 1.e3  ; mPa
 vs0(1,:,:,:) = vs0(1,:,:,:)*rho(1) * 1.e3  ; mPa

 ndir0 = dimsizes(dir0)
 dir = new(ndir0+1,float)
 dir(:ndir0-1) = dir0(:)
 dir(ndir0) = dir0(0) + 360

 vs = new((/2,2,dimsizes(c),ndir0+1/),float)
 vs(:,:,:,:ndir0-1) = vs0(:,:,:,:)
 vs(:,:,:,ndir0) = vs0(:,:,:,0)

; vs = where( vs .eq. 0., -999., vs )
; vs@_FillValue = -999.
; vsmall = 1.e-12
; vs = where( vs .eq. 0., vsmall, vs )

; d2r = 0.017453292519943
; c2d   = conform(vs,c  ,0)
; dir2d = conform(vs,dir,1)
; x = c2d*cos(dir2d*d2r)
; y = c2d*sin(dir2d*d2r)

;--------------------------------------- cn level
 cnlev = (/-160,-80,-40,-20,-10,-5,-2,2,5,10,20,40,80,160/)*0.1
 cnfil = (/17,18,19,20,21,24,26,0,10,8,6,5,4,3,2/)

;--------------------------------------- contour
 wi = 35
 he = wi
 xf = 10 + (/0,1,2,3/)*(wi+6)
 yf = 10 + (/0,1/)*(he+6)

 lbopt = new((/2,4/),string)
 lbopt = ""
 lbopt(1,1) = "[10:S:-4:N: mPa / deg m s:S:-1:N:]"

do id=0, 1
do il=0, 1
if ( .not. ( id .eq. 1 .and. il .eq. 3 ) ) then
 map = create "mp" mapPlotClass wks
  "mpProjection"          : "AzimuthalEquidistant"
  "mpOutlineBoundarySets" : "NoBoundaries"
  "mpGridAndLimbOn"       : True
;  "mpGridLineColor"       : colind
  "mpGridLineThicknessF"  : 0.5
  "mpGridLineDashPattern" : 2
  "mpEllipticalBoundary"  : True

  "mpLimitMode"    : "LatLon"
  "mpCenterLonF"   : 270
  "mpCenterLatF"   : 90
  "mpMinLonF"      : 0
  "mpMaxLonF"      : 360
  "mpMinLatF"      : 90-cmax(id)
  "mpMaxLatF"      : 90
  "mpGridLatSpacingF" : ctck
  "mpGridLonSpacingF" : 45
  "mpGridPolarLonSpacingF" : 45
 end create
 set_size_mp(map,xf(il),yf(id),wi,he)

 contour1 = cnshadeplot(wks,dir,90-c,vs(id,il,:,:)*1.e4,cnlev,cnfil,lbopt(id,il))
 setvalues contour1
;  "cnFillMode" : "CellFill"
;  "trGridType" : "Spherical"
  "pmLabelBarOrthogonalPosF" : 0.1
  "pmLabelBarParallelPosF"   : -0.2
  "pmLabelBarKeepAspect"     : False
  "pmLabelBarWidthF"         : 1.0
  "pmLabelBarHeightF"        : 0.10
  "lbLabelAutoStride"  : False
  "lbLabelStride"      : 1
  "lbLabelFontHeightF" : 0.015
  "lbTitlePosition"    : "Right"
  "lbTitleDirection"   : "Across"
  "lbTitleOffsetF"     : 0.1
  "lbTitleFontHeightF" : 0.02
 end setvalues

 set_size0(contour1,xf(il),yf(id),wi,he)
 axis_range(contour1,-999,999,90-cmax(id),999,False,False)
 axis_str(contour1,"","","")

 cnfinalize(contour1,wks)

 overlay(map,contour1)

 draw(map)

 axc = xylineplot(wks,(/0,1/),(/0,1/),"D",0,"T",1,"C",-1)
 set_size(axc,xf(il)+wi*0.5,yf(id),wi*0.5,he)
 axis_range(axc,0,cmax(id),-999,999,False,False)
 axis_str(axc,"","","")
 setvalues axc
  "tmYLBorderOn" : False
  "tmYRBorderOn" : False
  "tmXTBorderOn" : False
  "tmYLOn"       : False
  "tmYROn"       : False
  "tmXTOn"       : False
  "tmXBMode" : "Manual"
  "tmXBTickSpacingF" : ctck
  "tmXBMinorPerMajor" : 0
 end setvalues
 draw(axc)
 dr_txt(wks,"c [m s:S:-1:N:]",xf(il)+wi-4,yf(id)+he-2,21,4.,0.)

 delete(contour1)
 delete(map)

end if
end do
end do

 dr_txt(wks,"(a)",xf(0)-3,yf(0)-0.5,25,6,0.)
 dr_txt(wks,"(b)",xf(0)-3,yf(1)-0.5,25,6,0.)
 dr_txt(wks,"D1",xf(0)+3,yf(0)-0.5,22,6,0.)
 dr_txt(wks,"D2",xf(0)+3,yf(1)-0.5,22,6,0.)

 dr_txt(wks,"Zonal flux"     ,xf(0)+wi/2,yf(0)-4.5,22,5.5,0.)
 dr_txt(wks,"Meridional flux",xf(1)+wi/2,yf(0)-4.5,22,5.5,0.)

 frame(wks)

end

