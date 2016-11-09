load "$UTIL/ncl/header.ncl"

xspace = 170  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 0.85

begin

;--------------------------------------- parameters
 wks = cr_wks("fig/lb-fu_tz0",600,0)
 gsn_define_colormap(wks,"StepSeq25")
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- read data

;--------------------------------------- process 0

;--------------------------------------- cn level
 cnlev = fspan(-10,10,21)*2.
 cnfil = (/17,18,19,20,21,22,23,24,25,26,0,0,11,10,9,8,7,6,5,4,3,2/)

;--------------------------------------- contour
 wi = 64 + (64./13.*2.)
 he = 15 - 14.99
 xf = 17
 yf = 10

 lbopt = "(m s:S:-1:N: month:S:-1:N:)"

 dmp = (/ (/1,2,3/), (/4,5,6/), (/7,8,9/) /)*0.
 dmp@_FillValue = 0.
 contour1 = cnshadeplot(wks,(/-1,0,1/),(/-1,0,1/),dmp,cnlev,cnfil,lbopt)
 setvalues contour1
  "pmLabelBarOrthogonalPosF" : 0.12
  "pmLabelBarParallelPosF"   : 0.50;+0.05 ;+ 0.55
  "pmLabelBarKeepAspect"     : False
  "pmLabelBarWidthF"         : 0.6
  "pmLabelBarHeightF"        : 0.20
  "lbLabelStride"      : 2
  "lbLabelFontHeightF" : 0.012
  "lbTitlePosition"    : "Right"
  "lbTitleDirection"   : "Across"
  "lbTitleOffsetF"     : 0.05
  "lbTitleFontHeightF" : 0.015
  "tmXBOn" : False
  "tmXTOn" : False
  "tmYLOn" : False
  "tmYROn" : False
  "cnNoDataLabelOn" : False
 end setvalues
 set_size(contour1,xf,yf,wi,he)
 axis_range(contour1,-999,999,-999,999,False,False)
 axis_str(contour1,"","","")

 cnfinalize(contour1,wks)

 draw(contour1)

 frame(wks)

end

