begin

 wks = create "contour" psWorkstationClass defaultapp
  "wkPSFileName"  : "vertical.ps"
  "wkColorMap"    : "rainbow+white+gray"
  "wkOrientation"  : "landscape"
  "wkDeviceLowerX" : 0
  "wkDeviceUpperX" : 612
  "wkDeviceLowerY" : 90
  "wkDeviceUpperY" : 702
 end create

 fdir1 = "";"/home2/datatmp/GR1_CFLv/"
 fnam1 = "table1.yonust15_gr1cfl.jan1979-dec1988.nc"

 fdir2 = "";"/home3/DATA/NRA2/NetCDF/"
 fnam2 = "table1a.nra2.jan1979-dec1988.nc"

; READ DATA -------------------------------------------
 f1 = addfile(fdir1+fnam1,"r")
 ta1  = f1->ta-273.15
 lon1 = f1->longitude
 lat1 = f1->latitude
 pre1 = f1->pressure
 tim1 = f1->time

 nlon1 = dimsizes(lon1)
 nlat1 = dimsizes(lat1)
 npre1 = dimsizes(pre1)
 ntim1 = dimsizes(tim1)

 f2 = addfile(fdir2+fnam2,"r")
 tmp2 = f2->tmp-273.15
 lon2 = f2->longitude
 lat2 = f2->latitude
 pre2 = f2->level
 tim2 = f2->time

 nlon2 = dimsizes(lon2)
 nlat2 = dimsizes(lat2)
 npre2 = dimsizes(pre2)
 ntim2 = dimsizes(tim2)

; AVERAGE and Min., Max. --------------------------------
 ta_yonu = new((/npre1,nlat1/),float)
 ta_yonu(:,:) = 0.

 do k=0, npre1-1
 do j=0, nlat1-1
   do n=0, ntim1/12-1
   do i=0, nlon1-1
     ta_yonu(k,j) = ta_yonu(k,j) + ta1(n*12+1,k,j,i)   ; Jan.
   end do
   end do
 end do
 end do
 ta_yonu(:,:) = ta_yonu(:,:) / (nlon1*(ntim1/12))


 tmp_nra = new((/npre2,nlat2/),float)
 tmp_nra(:,:) = 0.

 do k=0, npre2-1
 do j=0, nlat2-1
   do n=0, ntim2/12-1
   do i=0, nlon2-1
     tmp_nra(k,j) = tmp_nra(k,j) + tmp2(n*12+1,k,j,i)   ; Jan.
   end do
   end do
 end do
 end do
 tmp_nra(:,:) = tmp_nra(:,:) / (nlon2*(ntim2/12))


 max_yonu = max(ta_yonu)
 min_yonu = min(ta_yonu)
 max_nra  = max(tmp_nra)
 min_nra  = min(tmp_nra)

; DEFINE SCALAR FIELD --------------------------------------
 field1 = create "w" scalarFieldClass defaultapp
  "sfDataArray" : ta_yonu
  "sfYArray" : pre1
  "sfXArray" : lat1
 end create

 field2 = create "w" scalarFieldClass defaultapp
  "sfDataArray" : tmp_nra
  "sfYArray" : pre2
  "sfXArray" : lat2
 end create

; DEFINE CONTOUR LEVELS and COLOR ------------------------

 cnlev = (/-70,-60,-50,-40,-30,-20,-10,0,5,10,15,20,25,30,35/)

 cnfil = new(dimsizes(cnlev)+1,integer)
 do i=0, dimsizes(cnfil)-1
   cnfil(i) = 10 + i*14
 end do

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

  "tmXBMode" : "Explicit" ;"Manual"
;  "tmXBTickSpacingF"  : 20
;  "tmXBMinorPerMajor" : 0
  "tmXBValues" : (/-90,-70,-50,-30,-10,10,30,50,70,90/)
  "tmXBLabels" : (/"90S","70S","50S","30S","10S","10N","30N","50N","70N","90N"/)
  "tmXBMajorLengthF" : -0.006
  "tmXBMinorLengthF" : -0.005
  "tmXBLabelFont"        : 21
  "tmXBLabelFontHeightF" : 0.009

  "tmYLMode" : "Explicit"
  "tmYLValues" : (/10.,20.,50.,100.,200.,500.,1000./)
  "tmYLLabels" : (/10,20,50,100,200,500,1000/)
  "tmYLMajorLengthF" : -0.006
  "tmYLMinorLengthF" : -0.005
  "tmYLLabelFont" : 21
  "tmYLLabelFontHeightF" : 0.009

  "vpXF" : 0.11
  "vpYF" : 0.70
  "vpWidthF" : 0.38
  "vpHeightF" : 0.20

  "tiMainString" : "YONU"
  "tiMainFont"   : 22
  "tiMainFontHeightF" : 0.017
  "tiMainOffsetYF"    : 0.02

  "tiXAxisString" : "latitude (deg)"
  "tiXAxisFont"   : 22
  "tiXAxisFontHeightF" : 0.015

  "tiYAxisString" : "pressure (hPa)"
  "tiYAxisFont"   : 22
  "tiYAxisFontHeightF" : 0.015
 end create

 plotover = create "plot" logLinPlotClass wks
  "vpXF" : 0.11
  "vpYF" : 0.70
  "vpWidthF" : 0.38
  "vpHeightF" : 0.20

  "trXMinF" : -90
  "trXMaxF" :  90
  "trYMinF" :  10.
  "trYMaxF" :  1000.
  "trXReverse" : True
  "trYReverse" : True
  "trYLog"  : True
 end create

 overlay(plotover,cont1)

 draw(plotover)

; CONTOUR (2st fig.) ----------------------------------------
 NhlRemoveOverlay(plotover,cont1,False)

 setvalues cont1
  "cnScalarFieldData" : field2
  "tmYLLabelsOn" : False
  "tmYRLabelsOn" : True
  "vpXF" : 0.53
  "tiMainString" : "NRA"
  "tiYAxisString" : ""
 end setvalues

 setvalues plotover
  "vpXF" : 0.53
 end setvalues

 overlay(plotover,cont1)

 draw(plotover)


; TYPE Max. and Min. --------------------------------
 tid = create "Text" textItemClass wks
  "txString" : "Max:: "+floattoint(max_yonu*10.)/10.
  "txPosXF"  : 0.11 + 0.05 
  "txPosYF"  : 0.70 + 0.02
  "txFont"   : 21
  "txFontHeightF" : 0.013
 end create
 draw(tid)
 setvalues tid  
  "txString" : "Min:: "+floattoint(min_yonu*10.)/10.
  "txPosXF"  : 0.11+0.38 - 0.05
 end setvalues 
 draw(tid)

 setvalues tid
  "txString" : "Max:: "+floattoint(max_nra*10.)/10.
  "txPosXF"  : 0.53 + 0.05
 end setvalues
 draw(tid)
 setvalues tid
  "txString" : "Min:: "+floattoint(min_nra*10.)/10.
  "txPosXF"  : 0.53+0.38 - 0.05
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
  "lbLabelFontHeightF" : 0.009
  "lbLabelPosition"  : "Bottom"
  "lbLabelAlignment" : "InteriorEdges"
  "lbPerimOn" : False
  "lbTitleString"    : " T [C]"
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



