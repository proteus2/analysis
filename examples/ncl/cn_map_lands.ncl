begin

 wks = create "contour" psWorkstationClass defaultapp
  "wkPSFileName"  : "horizontal.ps"
  "wkColorMap"    : "rainbow+white+gray"
  "wkOrientation"  : "landscape"
  "wkDeviceLowerX" : 0
  "wkDeviceUpperX" : 612
  "wkDeviceLowerY" : 90
  "wkDeviceUpperY" : 702
 end create

 fdir1 = "/home2/datatmp/GR1_CFLv/"
 fnam1 = "table2.yonust15_gr1cfl.jan1979-dec1988.nc"

 fdir2 = "/home3/DATA/NRA2/NetCDF/"
 fnam2 = "table2.nra2.jan1979-dec1988.nc"

; READ DATA ---------------------------------
 f1 = addfile(fdir1+fnam1,"r")
 var1 = f1->rlut*(-1.)
 lon1  = f1->longitude
 lat1  = f1->latitude
 tim1  = f1->time

 nlon1 = dimsizes(lon1)
 nlat1 = dimsizes(lat1)
 ntim1 = dimsizes(tim1)

 f2 = addfile(fdir2+fnam2,"r")
 var2 = f2->ulwrftoa
 lon2 = f2->longitude
 lat2 = f2->latitude
 tim2 = f2->time

 nlon2 = dimsizes(lon2)
 nlat2 = dimsizes(lat2)
 ntim2 = dimsizes(tim2)

; AVERAGE and Min., Max. --------------------------------
 olr_yonu = new((/nlat1,nlon1/),float)
 olr_yonu(:,:) = 0.

 do j=0, nlat1-1
 do i=0, nlon1-1
   do n=0, ntim1/12-1
     olr_yonu(j,i) = olr_yonu(j,i) + var1(n*12+1,j,i)   ; Jan.
   end do
 end do
 end do
 olr_yonu(:,:) = olr_yonu(:,:) / (ntim1/12)


 olr_nra = new((/nlat2,nlon2/),float)
 olr_nra(:,:) = 0.

 do j=0, nlat2-1
 do i=0, nlon2-1
   do n=0, ntim2/12-1
     olr_nra(j,i) = olr_nra(j,i) + var2(n*12+1,j,i)   ; Jan.
   end do
 end do
 end do
 olr_nra(:,:) = olr_nra(:,:) / (ntim2/12)


 max_yonu = max(olr_yonu)
 min_yonu = min(olr_yonu)
 max_nra  = max(olr_nra)
 min_nra  = min(olr_nra)

; DEFINE SCALAR FIELD --------------------------------------
 field1 = create "w" scalarFieldClass defaultapp
  "sfDataArray" : olr_yonu
  "sfYArray" : lat1
  "sfXArray" : lon1
 end create

 field2 = create "w" scalarFieldClass defaultapp
  "sfDataArray" : olr_nra
  "sfYArray" : lat2
  "sfXArray" : lon2
 end create

; DEFINE CONTOUR LEVELS and COLOR -------------------------

 cnlev = new(14,float)
 do i=0, dimsizes(cnlev)-1
   cnlev(i) = 105. + i*15.
 end do

 cnfil = new(dimsizes(cnlev)+1,integer)
 do i=0, dimsizes(cnfil)-1
   cnfil(i) = 230 - i*14
 end do

; MAP (1st fig.) -----------------------------------------
 map = create "mp" mapPlotClass wks
  "vpXF"      : 0.11
  "vpYF"      : 0.70
  "vpWidthF"  : 0.38
  "vpHeightF" : 0.20

  "mpProjection"    : "CylindricalEquidistant"
  "mpFillOn"        : False
  "mpLabelsOn"      : False
  "mpOutlineOn"     : True
  "mpGridAndLimbOn" : False
  "mpPerimOn"       : True
  "mpPerimDrawOrder" : "PostDraw"

  "mpShapeMode"     : "FixedAspectFitBB"
  "mpGeophysicalLineColor" : "black"
  "mpDataBaseVersion" : "Ncarg4_0"

  "mpLimitMode"    : "LatLon"
  "mpMinLonF"      :   0
  "mpMaxLonF"      : 360
  "mpMinLatF"      : -90
  "mpMaxLatF"      :  90
  "mpCenterLonF"   : 180
 end create
 getvalues map
  "vpWidthF"  : wfg
  "vpHeightF" : hfg
  "vpXF"      : xfg
  "vpYF"      : yfg
 end getvalues

; CONTOUR (1st fig.) ----------------------------------------
 cont1 = create "contour" contourPlotClass wks
  "cnScalarFieldData" : field1
  "cnLinesOn"         : False
  "cnFillOn"          : True
  "cnLevelSelectionMode" : "ExplicitLevels"
  "cnLevels"     : cnlev
  "cnFillColors" : cnfil
  "cnInfoLabelOn"  : False
  "cnLineLabelsOn" : False

  "tmXBMode" : "Explicit"
  "tmXBValues" : (/0,30,60,90,120,150,180,210,240,270,300,330,360/)
  "tmXBLabels" : (/"0","30E","60E","90E","120E","150E","180","150W","120W","90W","60W","30W","0"/)
  "tmXBMajorLengthF" : -0.004
  "tmXBMinorLengthF" : -0.003
  "tmXBLabelFont"        : 21
  "tmXBLabelFontHeightF" : 0.007

  "tmYLMode" : "Explicit"
  "tmYLValues" : (/-90,-70,-50,-30,-10,10,30,50,70,90/)
  "tmYLLabels" : (/"90S","70S","50S","30S","10S","10N","30N","50N","70N","90N"/)
  "tmYLMajorLengthF" : -0.004
  "tmYLMinorLengthF" : -0.003
  "tmYLLabelFont"        : 21
  "tmYLLabelFontHeightF" : 0.007

  "vpXF" : xfg
  "vpYF" : yfg
  "vpWidthF"  : wfg
  "vpHeightF" : hfg

  "tiMainString" : "YONU"
  "tiMainFont"   : 22
  "tiMainFontHeightF" : 0.017
  "tiMainOffsetYF"    : 0.02
  "tiXAxisString" : ""
  "tiYAxisString" : ""
 end create

 plotover = create "plot" logLinPlotClass wks
  "vpXF" : xfg
  "vpYF" : yfg
  "vpWidthF"  : wfg
  "vpHeightF" : hfg

  "trXMinF" :   0
  "trXMaxF" : 360
  "trYMinF" : -90
  "trYMaxF" :  90
 end create

 overlay(plotover,cont1)

; DRAW (1st fig.) ------------------------------------
 draw(plotover)

 draw(map)

; TYPE Max. and Min. --------------------------------
 tid = create "Text" textItemClass wks
  "txString" : "Max:: "+floattoint(max_yonu*10.)/10.
  "txPosXF"  : xfg + 0.05
  "txPosYF"  : yfg + 0.02
  "txFont"   : 21
  "txFontHeightF" : 0.013
 end create
 draw(tid)
 setvalues tid
  "txString" : "Min:: "+floattoint(min_yonu*10.)/10.
  "txPosXF"  : xfg+wfg - 0.05
 end setvalues
 draw(tid)



; REPEAT (2nd fig.) ----------------------------------
 NhlRemoveOverlay(plotover,cont1,False)

 setvalues map
  "vpXF" : 0.53
  "vpWidthF" : wfg
 end setvalues
 getvalues map
  "vpXF" : xfg
 end getvalues
 setvalues cont1
  "cnScalarFieldData" : field2
  "tmYLLabelsOn" : False
  "tmYRLabelsOn" : True
  "vpXF" : xfg
  "tiMainString" : "NRA"
 end setvalues
 setvalues plotover
  "vpXF" : xfg
 end setvalues

 overlay(plotover,cont1)
 draw(plotover)

 draw(map)

 setvalues tid
  "txString" : "Max:: "+floattoint(max_nra*10.)/10.
  "txPosXF"  : xfg + 0.05
 end setvalues
 draw(tid)
 setvalues tid
  "txString" : "Min:: "+floattoint(min_nra*10.)/10.
  "txPosXF"  : xfg+wfg - 0.05
 end setvalues
 draw(tid)


; LABEL BAR ---------------------------------------
 label = create "lab" labelBarClass wks
  "lbOrientation" : "horizontal"
  "lbBoxCount"    : dimsizes(cnfil)
  "lbFillColors"  : cnfil
  "lbMonoFillPattern" : True
  "lbLabelsOn"     : True
  "lbLabelStrings" : cnlev
  "lbAutoManage"  : False
  "lbLabelFont"   : 21
  "lbLabelFontHeightF" : 0.007
  "lbLabelPosition"  : "Bottom"
  "lbLabelAlignment" : "InteriorEdges"
  "lbPerimOn" : False
  "lbTitleString"    : " OLR [W m:S:-2:N:]"
  "lbTitlePosition"  : "Right"
  "lbTitleDirection" : "Across"
  "lbTitleFont" : 22
  "lbTitleFontHeightF" : 0.013
  "vpXF" : 0.28
  "vpYF" : 0.40
  "vpWidthF"  : 0.45
  "vpHeightF" : 0.03
 end create
 draw(label)

;-------------------------------------------------

 frame(wks)

end



