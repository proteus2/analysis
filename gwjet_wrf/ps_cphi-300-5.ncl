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
 var1 = "W"
 lev = 8000
 var1name = "prt_"+var1
;	Time = UNLIMITED ; // (1152 currently)
 
 wks = cr_wks("fig/ps_cphi-300"+var1name+"-5",600,0)
 gsn_define_colormap(wks,"prcp_3")
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- read data
 fname1 = fdir+"/x6_5min/W_z/prt_d300km/psd/psd_"+var1name+"_phi_2d__z"+sprinti("%5.5i",lev)+"-w"+ispan(1,5,1)+".nc"
 fname0 = fdir+"/x6_5min/W_z/prt_d300km/psd/psd_"+var1name+"_phi_2d__z"+sprinti("%5.5i",lev)+"-"+ispan(1,2,1)+"a.nc"
 f1 = addfile(fname1(0),"r")
 f2 = addfile(fname1(1),"r")
 f3 = addfile(fname1(2),"r")
 f4 = addfile(fname1(3),"r")
 f5 = addfile(fname1(4),"r")
 f01 = addfile(fname0(0),"r")
 f02 = addfile(fname0(1),"r")

 vs0 = (/ (/ f01->psd_phi_c({:max(cmax)+5},:), \
             f1->psd_phi_c({:max(cmax)+5},:), \
             f2->psd_phi_c({:max(cmax)+5},:), \
             f3->psd_phi_c({:max(cmax)+5},:) /), \
          (/ f02->psd_phi_c({:max(cmax)+5},:), \
             f4->psd_phi_c({:max(cmax)+5},:), \
             f5->psd_phi_c({:max(cmax)+5},:), \
             f5->psd_phi_c({:max(cmax)+5},:) /) /)
 dir0= f1->dir
 c   = f1->c({:max(cmax)+5})

;--------------------------------------- process 0
 ndir0 = dimsizes(dir0)
 dir = new(ndir0+1,float)
 dir(:ndir0-1) = dir0(:)
 dir(ndir0) = dir0(0) + 360

 vs = new((/2,4,dimsizes(c),ndir0+1/),float)
 vs(:,:,:,:ndir0-1) = vs0(:,:,:,:)
 vs(:,:,:,ndir0) = vs0(:,:,:,0)

; vs = where( vs .eq. 0., -999., vs )
; vs@_FillValue = -999.
 vsmall = 1.e-12
 vs = where( vs .eq. 0., vsmall, vs )

; d2r = 0.017453292519943
; c2d   = conform(vs,c  ,0)
; dir2d = conform(vs,dir,1)
; x = c2d*cos(dir2d*d2r)
; y = c2d*sin(dir2d*d2r)

;--------------------------------------- cn level
 if (var1 .eq. "W") then
;   cnlev = fspan(-12.0,-7.5,10)*1.
   cnlev = fspan(-9.5,-7.,11)*1.
 end if
 if (var1 .eq. "U") then
   cnlev = fspan(-4.5,0.0,10)*1.
 end if
 cnfil = (/0,5,6,7,10,12,14,15,16,17,19,21/)

;--------------------------------------- contour
 wi = 35
 he = wi
 xf = 10 + (/3,0,1,2/)*(wi+2) + (/5,0,0,0/)
 yf = 10 + (/0,1/)*(he+8)

 lbopt = new((/2,4/),string)
 lbopt = ""
; lbopt(0,3) = "log:B:10:N:(PSD [m:S:2:N: s:S:-2:N: / deg m s:S:-1:N:])"
 lbopt(1,2) = "log:B:10:N:(PSD [m:S:2:N: s:S:-2:N: / deg m s:S:-1:N:])"

 txtstr = (/ (/"D1","W1","W2","W3"/), (/"D2","W4","W5",""/) /)

do id=0, 1
do il=0, 3
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

 contour1 = cnshadeplot(wks,dir,90-c,log10(vs(id,il,:,:)),cnlev,cnfil,lbopt(id,il))
 setvalues contour1
;  "cnFillMode" : "CellFill"
;  "trGridType" : "Spherical"
  "pmLabelBarOrthogonalPosF" : 0.1 ;1.0
;  "pmLabelBarParallelPosF"   : -1.20;+0.05 ;+ 0.55
  "pmLabelBarKeepAspect"     : False
  "pmLabelBarWidthF"         : 0.75 ;0.5
  "pmLabelBarHeightF"        : 0.10
  "lbLabelStride"      : 2
;  "lbLabelFontHeightF" : 0.015
;  "lbTitlePosition"    : "Right"
;  "lbTitleDirection"   : "Across"
;  "lbTitleOffsetF"     : 0.06
;  "lbTitleFontHeightF" : 0.02
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

 dr_txt(wks,txtstr(id,il),xf(il)+wi*0.5,yf(id)-1.2,22,5.,0.)

 delete(contour1)
 delete(map)

end if
end do
end do

; dr_txt(wks,"(c)",xf(0)-3,yf(0)-1.0,25,6,0.)
; dr_txt(wks,"(d)",xf(0)-3,yf(1)-1.0,25,6,0.)
 dr_txt(wks,"(a)",xf(1)-3,yf(0)-1.0,25,6,0.)
 dr_txt(wks,"(b)",xf(1)-3,yf(1)-1.0,25,6,0.)

 frame(wks)

end

