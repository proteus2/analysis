load "$UTIL/ncl/header.ncl"

xspace = 170  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

;axisfontratio = 1.0

begin

;--------------------------------------- parameters
 dd = 2
 cmax = 30
 ctck = 5
 fdir = "$DATD/gwjet_wrf"
 lev = 8000
;	Time = UNLIMITED ; // (1152 currently)
 
 wks = cr_wks("fig/mfys_z_cphi-300-"+(lev/1.e3)+"-"+dd,600,0)
; gsn_define_colormap(wks,"prcp_3")
 gsn_define_colormap(wks,"StepSeq25")
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- read data
 fname1 = fdir+"/x6_5min/mfs_z/prt_d300km/mfs_prt_phi_2d__z"+sprinti("%5.5i",lev)+"-"+dd+".nc"
 f1 = addfile(fname1,"r")
 vs0 = f1->wv_phi_c({:cmax},:)
 dir0= f1->dir
 c   = f1->c({:cmax})

;--------------------------------------- process 0
 rho = (/0.465,0.500/)  ; at 8 km in D1 and D2
 vs0 = vs0*rho(dd-1) * 1.e3  ; mPa

 ndir0 = dimsizes(dir0)
 dir = new(ndir0+1,float)
 dir(:ndir0-1) = dir0(:)
 dir(ndir0) = dir0(0) + 360

 vs = new((/dimsizes(c),ndir0+1/),float)
 vs(:,:ndir0-1) = vs0(:,:)
 vs(:,ndir0) = vs0(:,0)

; vs = where( vs .eq. 0., -999., vs )
; vs@_FillValue = -999.

; d2r = 0.017453292519943
; c2d   = conform(vs,c  ,0)
; dir2d = conform(vs,dir,1)
; x = c2d*cos(dir2d*d2r)
; y = c2d*sin(dir2d*d2r)

;--------------------------------------- cn level
; cnlev = array_append_record( 10^(fspan(-3.25,-5.75,11)*1.)*(-1.0), \
;                              10^(fspan(-5.75,-3.25,11)*1.), 0 )
; cnfil = (/21,19,17,16,15,14,12,10,7,6,5,0,5,6,7,10,12,14,15,16,17,19,21/)
; cnlev = ispan(-10,10,1)*5.
 cnlev = (/-160,-80,-40,-20,-10,-5,-4,-3,-2,2,3,4,5,10,20,40,80,160/)
 cnfil = (/17,18,19,20,21,22,24,25,26,0,11,10,9,7,6,5,4,3,2/)

;--------------------------------------- contour
 wi = 50
 he = wi
 xf = 12 + (/0,1,2,3,0,1,2,3/)*(wi+6)
 yf = 10 + (/0,0,0,0,1,1,1,1/)*(he+6)

 lbopt = "[10:S:-5:N: mPa / deg m s:S:-1:N:]"

do il=0, 0

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
  "mpMinLatF"      : 90-cmax
  "mpMaxLatF"      : 90
  "mpGridLatSpacingF" : ctck
  "mpGridLonSpacingF" : 45
  "mpGridPolarLonSpacingF" : 45
 end create
 set_size_mp(map,xf(0),yf(0),wi,he)

; contour1 = cnshadeplot(wks,x,y,log10(vs),cnlev,cnfil,lbopt(il))
 contour1 = cnshadeplot(wks,dir,90-c,vs*1.e5,cnlev,cnfil,lbopt(il))
; cnaddline(contour1,wks,dir,90-c,vs,10.^(-5.75)*(/-1,1/),0,"D",(/0,0,14/),"T",1,"C",1)
 setvalues contour1
;  "cnFillMode" : "CellFill"
;  "trGridType" : "Spherical"
;  "pmLabelBarOrthogonalPosF" : 0.08
;  "pmLabelBarParallelPosF"   : -1.20;+0.05 ;+ 0.55
  "pmLabelBarKeepAspect"     : False
  "pmLabelBarWidthF"         : 0.60
  "pmLabelBarHeightF"        : 0.08
  "lbLabelAutoStride"      : False
  "lbLabelStride"      : 1
  "lbLabelFontHeightF" : 0.015
;  "lbTitlePosition"    : "Right"
;  "lbTitleDirection"   : "Across"
;  "lbTitleOffsetF"     : 0.06
  "lbTitleFontHeightF" : 0.02
 end setvalues

 set_size0(contour1,xf(il),yf(il),wi,he)
 axis_range(contour1,-999,999,90-cmax,999,False,False)
 axis_str(contour1,"","","")

 cnfinalize(contour1,wks)

 overlay(map,contour1)
 
 draw(map)

 axc = xylineplot(wks,(/0,1/),(/0,1/),"D",0,"T",1,"C",-1)
 set_size(axc,xf(il)+wi*0.5,yf(il),wi*0.5,he)
 axis_range(axc,0,cmax,-999,999,False,False)
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
 dr_txt(wks,"c [m s:S:-1:N:]",xf(il)+wi-4,yf(il)+he-2,21,4.5,0.)

 delete(contour1)
 delete(map)

end do

 frame(wks)

end

