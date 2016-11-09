load "$UTIL/ncl/header.ncl"

xspace = 170  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 0.75

begin

;--------------------------------------- parameters
 cmax = (/25,25,15/)
 ctck = 5
 fdir = "$DATD/gwjet_wrf"
 var1 = "W"
 lev = (/250,1500,3000,5000/)
 var1name = "prt_"+var1
;	Time = UNLIMITED ; // (1152 currently)
 
 wks = cr_wks("fig/ps_cphi-300-"+var1name+"-3-4z",600,0)
 gsn_define_colormap(wks,"prcp_3")
 colind = NhlNewColor(wks,0.4,0.4,0.4)

do k=0, dimsizes(lev)-1
;--------------------------------------- read data
; fname1 = fdir+"/x6_5min/W_z/prt_d300km/psd/psd_"+var1name+"_phi_2d__z"+sprinti("%5.5i",lev)+"-w"+ispan(1,5,1)+".nc"
 fname0 = fdir+"/x6_5min/W_z/prt_d300km/psd/psd_"+var1name+"_phi_2d__z"+sprinti("%5.5i",lev(k))+"-"+ispan(1,2,1)+"a.nc"
 fnamef = fdir+"/x6_5min/W_z/prt_d300km/psd/psd_"+var1name+"_phi_2d__z"+sprinti("%5.5i",lev(k))+"-1a-f150.nc"
; f1 = addfile(fname1(0),"r")
; f2 = addfile(fname1(1),"r")
; f3 = addfile(fname1(2),"r")
; f4 = addfile(fname1(3),"r")
; f5 = addfile(fname1(4),"r")
 f01 = addfile(fname0(0),"r")
 f02 = addfile(fname0(1),"r")
 ff1 = addfile(fnamef,"r")

; vs0 = (/ (/ f01->psd_phi_c({:max(cmax)+5},:), \
;             f1->psd_phi_c({:max(cmax)+5},:), \
;             f2->psd_phi_c({:max(cmax)+5},:), \
;             f3->psd_phi_c({:max(cmax)+5},:) /), \
;          (/ f02->psd_phi_c({:max(cmax)+5},:), \
;             f4->psd_phi_c({:max(cmax)+5},:), \
;             f5->psd_phi_c({:max(cmax)+5},:), \
;             f5->psd_phi_c({:max(cmax)+5},:) /) /)
 vs0 = (/ f01->psd_phi_c({:max(cmax)+5},:) - ff1->psd_phi_c({:max(cmax)+5},:), \
          ff1->psd_phi_c({:max(cmax)+5},:), \
          f02->psd_phi_c({:max(cmax)+5},:) /) * exp(-lev(k)/7.e3)
 dir0= f01->dir
 c   = f01->c({:max(cmax)+5})

;--------------------------------------- process 0
 ndir0 = dimsizes(dir0)
 dir = new(ndir0+1,float)
 dir(:ndir0-1) = dir0(:)
 dir(ndir0) = dir0(0) + 360

 vs = new((/3,dimsizes(c),ndir0+1/),float)
 vs(:,:,:ndir0-1) = vs0(:,:,:)
 vs(:,:,ndir0) = vs0(:,:,0)

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
 xf = 10 + (/0,1,2,3,4/)*(wi+2) + (/0,0,0,0,0/)
 yf = 10 + (/0,1,2,3/)*(he+6)

 lbopt = (/"","","log:B:10:N:(PSD [m:S:2:N: s:S:-2:N: / deg m s:S:-1:N:])"/)
 if (k .ne. 1) then
   lbopt = ""
 end if

 txtstr = "z = "+(lev(k)/1000.)+" km"
 if (lev(k) .lt. 1000) then
   txtstr = "z = "+lev(k)+" m"
 end if

do id=0, 2
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
 set_size_mp(map,xf(k),yf(id),wi,he)

 contour1 = cnshadeplot(wks,dir,90-c,log10(vs(id,:,:)),cnlev,cnfil,lbopt(id))
 setvalues contour1
;  "cnFillMode" : "CellFill"
;  "trGridType" : "Spherical"
  "pmLabelBarOrthogonalPosF" : 0.1 ;1.0
  "pmLabelBarParallelPosF"   : 1.00
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

 set_size0(contour1,xf(k),yf(id),wi,he)
 axis_range(contour1,-999,999,90-cmax(id),999,False,False)
 axis_str(contour1,"","","")

 cnfinalize(contour1,wks)

 overlay(map,contour1)
 
 draw(map)

 axc = xylineplot(wks,(/0,1/),(/0,1/),"D",0,"T",1,"C",-1)
 set_size(axc,xf(k)+wi*0.5,yf(id),wi*0.5,he)
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
 dr_txt(wks,"c [m s:S:-1:N:]",xf(k)+wi-4,yf(id)+he-2,21,4.,0.)

if (id .eq. 0) then
 dr_txt(wks,txtstr(id),xf(k)+wi*0.5,yf(id)-2.5,22,6.,0.)
end if

 delete(contour1)
 delete(map)

end do

end do

; dr_txt(wks,"(c)",xf(0)-3,yf(0)-1.0,25,6,0.)
; dr_txt(wks,"(d)",xf(0)-3,yf(1)-1.0,25,6,0.)
 dr_txt(wks,"(a)",xf(0)+3,yf(0)-1.0,25,6,0.)
 dr_txt(wks,"(b)",xf(0)+3,yf(1)-1.0,25,6,0.)
 dr_txt(wks,"(c)",xf(0)+3,yf(2)-1.0,25,6,0.)
; dr_txt(wks,"D1 (Lh > 150 km)",xf(0)+10,yf(0)+he-2.0,22,5,0.)
; dr_txt(wks,"D1 (Lh < 150 km)",xf(0)+10,yf(1)+he-2.0,22,5,0.)
; dr_txt(wks,"D2",xf(0)+4,yf(2)+he-2.0,22,5,0.)

 frame(wks)

end

