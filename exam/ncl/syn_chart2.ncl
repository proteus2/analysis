begin   ; synoptic chart

 setvalues NhlGetWorkspaceObjectId()
   "wsMaximumSize" : 55000000
 end setvalues

 load "/usr/local/ncarg/lib/ncarg/nclscripts/csm/gsn_code.ncl"

 wks = create "contour" psWorkstationClass defaultapp
     "wkPSFileName" : "syn2002.ps"
;     "wkColorMap" : "BlWhRe"
 end create

 klevt = 14      ; index for temperature

 inixf = 0.03
 iniyf = 1.00

 xf = (/0.,0.,0.,0.,0.5,0.5,0.5,0.5/) + inixf
 yf = iniyf - (/0.,0.25,0.5,0.75,0.,0.25,0.5,0.75/)
;-------------------------------------------------------------
 file1 = addfile("200807.nc","r")
 x = file1->lon
 y = file1->lat
 p = file1->lev

 psfc = file1->P_MSL/100.
 t = file1->T-273.


 print(" T lev. :  "+p(klevt))

;-------------------------------------------------------------
do i=0, 7
;-------------------------------------------------------------
 field1 = create "f1" scalarFieldClass defaultapp
  "sfDataArray" : psfc(i+39,:,:);psfc(i+8,:,:)
  "sfXArray" : x
  "sfYArray" : y
 end create

 field2 = create "f1" scalarFieldClass defaultapp
  "sfDataArray" : t(i+39,klevt,:,:);t(i+8,klevt,:,:) 
  "sfXArray" : x
  "sfYArray" : y
 end create 

;-------------------------------------------------------------
 map = create "mp" mapPlotClass wks
  "vpXF"      : xf(i)
  "vpYF"      : yf(i)
  "vpWidthF"  : 0.45
  "vpHeightF" : 0.24

  "mpProjection"    : "LambertConformal"
  "mpLambertMeridianF" : 115
  "mpFillOn"        : True;False
  "mpFillColors"    : (/1,0,11,11/)
  "mpLabelsOn"      : True;False
  "mpOutlineOn"     : True
  "mpGridAndLimbOn" : True
  "mpGridLineColor" : "grey"
  "mpGridLineDashPattern" : 1
  "mpGridLineDashSegLenF" : 0.2
  "mpPerimOn"       : True

  "mpShapeMode"     : "FixedAspectFitBB"
  "mpGeophysicalLineColor" : "Orange4"
  "mpDataBaseVersion" : "Ncarg4_0"

  "mpGridLonSpacingF" : 10
  "mpGridLatSpacingF" : 10
  "mpLimitMode"    : "LatLon"
  "mpMinLonF"      : 80.
  "mpMaxLonF"      : 150.
  "mpMinLatF"      : 20.
  "mpMaxLatF"      : 60.
 end create
 getvalues map
  "vpWidthF"  : wfg
  "vpHeightF" : hfg
  "vpXF"      : xfg
  "vpYF"      : yfg
 end getvalues

;------------------------------------------------------------
 cont1 = create "contour" contourPlotClass wks
  "cnScalarFieldData" : field1
  "cnLevelSelectionMode" : "ManualLevels"
  "cnLevelSpacingF" : 2
;  "cnFillOn"       :  True
;  "cnFillColors"   : cnfil
;  "cnLinesOn"      : False
  "cnInfoLabelOn" : False
  "cnLineLabelPlacementMode" : "Constant"
  "cnLineLabelFont" : 21
  "cnLineLabelFontHeightF" : 0.007
  "cnLineLabelsOn" : True

  "cnHighLabelsOn"         : True
  "cnHighLabelFontColor"   : "blue"
  "cnHighLabelString"      : "H"
  "cnHighLabelFontHeightF" : 0.01
  "cnHighLabelPerimSpaceF" : 0.05
  "cnLowLabelsOn"          : True
  "cnLowLabelFontColor"    : "red"
  "cnLowLabelString"       : "L"
  "cnLowLabelFontHeightF"  : 0.01
  "cnLowLabelPerimSpaceF"  : 0.05

  "cnSmoothingDistanceF"  : 0.005
  "cnSmoothingTensionF"   : -20.0
  "cnSmoothingOn" : True

  "tmXBOn" : False
  "tmYLOn" : False
  "tmXTOn" : False
  "tmYROn" : False

  "vpXF" : xfg
  "vpYF" : yfg
  "vpWidthF" : wfg
  "vpHeightF" : hfg

  "trXMinF" : 55
  "trXMaxF" : 175
  "trYMinF" : 10
  "trYMaxF" : 70

 end create

 cont2 = create "contour" contourPlotClass wks
  "cnScalarFieldData" : field2
  "cnLevelSelectionMode" : "ManualLevels"
  "cnLevelSpacingF" : 4
;  "cnFillOn"       :  True
;  "cnFillColors"   : cnfil
;  "cnLinesOn"      : False
  "cnInfoLabelOn" : False
  "cnLineLabelPlacementMode" : "Constant"
  "cnLineLabelFont" : 21
  "cnLineLabelFontHeightF" : 0.007
  "cnLineLabelsOn"  : True
  "cnLineLabelFontColor" : "red"
  "cnMonoLineColor" : True
  "cnLineColor"     : "red"
  "cnMonoLineDashPattern" : True
  "cnLineDashPattern" : 2

  "cnHighLabelsOn"         : True
  "cnHighLabelFontColor"   : "red"
  "cnHighLabelString"      : "W"
  "cnHighLabelFontHeightF" : 0.01
  "cnHighLabelPerimSpaceF" : 0.05
  "cnLowLabelsOn"          : True
  "cnLowLabelFontColor"    : "blue"
  "cnLowLabelString"       : "C"
  "cnLowLabelFontHeightF"  : 0.01
  "cnLowLabelPerimSpaceF"  : 0.05

  "cnSmoothingDistanceF"  : 0.001
;  "cnSmoothingTensionF"   : 1.0
  "cnSmoothingOn" : True

  "tmXBOn" : False
  "tmYLOn" : False
  "tmXTOn" : False
  "tmYROn" : False

  "vpXF" : xfg
  "vpYF" : yfg
  "vpWidthF" : wfg
  "vpHeightF" : hfg

  "trXMinF" : 55
  "trXMaxF" : 175
  "trYMinF" : 10
  "trYMaxF" : 70
 end create
 
 overlay(map,cont1)
 NhlAddOverlay(map,cont2,map)
 draw(map)


 delete(field1)
 delete(map)
 delete(cont1)
 delete(cont2)
;--------------------------------------------------------------
end do
;--------------------------------------------------------------

 frame(wks)

end
