begin

 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

;--------------------------------------------

wks = create "1sttry" psWorkstationClass defaultapp
      "wkPSFileName" : "sa2b2.ps"
      "wkColorMap"   : "BlWhRe"
end create

;---------------------------------------

 file1 = "sa2b2pic.nc"
 f1    = addfile(file1,"r")
 lag   = f1->lag
 freq  = f1->freq
 lat   = f1->lat1
 cor   = f1->COR
 coh   = f1->COH
 phs   = f1->PHASE
 amp   = f1->AMP

 lat1 = new(dimsizes(lat)+2,float)
 lat1(0) = -60.
 lat1(dimsizes(lat)+1) = 60.
 lat1(1:dimsizes(lat)) = lat(:)
 cor1 = new((/dimsizes(lat)+2,dimsizes(lag)/),float)
 cor1(0,:) = cor(0,:)
 cor1(dimsizes(lat)+1,:) = cor(dimsizes(lat)-1,:)
 cor1(1:dimsizes(lat),:) = cor(:,:)
 coh1 = new((/dimsizes(lat)+2,dimsizes(freq)/),float)
 coh1(0,:) = coh(0,:)
 coh1(dimsizes(lat)+1,:) = coh(dimsizes(lat)-1,:)
 coh1(1:dimsizes(lat),:) = coh(:,:)
 phs1 = new((/dimsizes(lat)+2,dimsizes(freq)/),float)
 phs1(0,:) = phs(0,:)
 phs1(dimsizes(lat)+1,:) = phs(dimsizes(lat)-1,:)
 phs1(1:dimsizes(lat),:) = phs(:,:)
 amp1 = new((/dimsizes(lat)+2,dimsizes(freq)/),float)
 amp1(0,:) = amp(0,:)
 amp1(dimsizes(lat)+1,:) = amp(dimsizes(lat)-1,:)
 amp1(1:dimsizes(lat),:) = amp(:,:)


;---------------------------------- cn level

 cnlev1 = new(5,float)
 cnfil1 = new(6,integer)
 do i=0,4
  cnlev1(i) = (i-2) * 0.1
 end do
 cnfil1 = (/-12,-6,-2,0,6,12/) + 52

 cnlev2 = new(4,float)
 cnfil2 = new(5,integer)
 do i=0,3
  cnlev2(i) = (i+1) * 0.2
 end do
 cnfil2 = (/0,-2,-6,-10,-20/) + 52

 cnlev3 = new(11,float)
 cnfil3 = new(12,integer)
 do i=0,10
  cnlev3(i) = (i-5) * 30.
 end do
 cnfil3 = (/-40,-33,-20,-7,-4,-2,0,4,7,20,33,40/) + 52

 cnlev4 = new(5,float)
 cnfil4 = new(6,integer)
 do i=0,4
  cnlev4(i) = (i+1) * 50
 end do
 cnfil4 = (/0,-3,-6,-9,-14,-20/) + 52


;------------------------------------------

 field = new(4,graphic)
 field(0) = create "w" scalarFieldClass defaultapp
   "sfDataArray" : cor1
   "sfYArray" : lat1
   "sfXArray" : lag
 end create
 field(1) = create "w" scalarFieldClass defaultapp
   "sfDataArray" : coh1
   "sfYArray" : lat1
   "sfXArray" : freq
 end create
 field(2) = create "w" scalarFieldClass defaultapp
   "sfDataArray" : phs1
   "sfYArray" : lat1
   "sfXArray" : freq
 end create 
 field(3) = create "w" scalarFieldClass defaultapp
   "sfDataArray" : amp1
   "sfYArray" : lat1
   "sfXArray" : freq 
 end create


  cont1 = create "contour" contourPlotClass wks
   "cnScalarFieldData" : field(0)
   "cnFillOn"        :  True           ; color fill
   "cnLevelSelectionMode" : "ExplicitLevels"
   "cnLevels" : cnlev1
   "cnFillColors" : cnfil1
   "cnLinesOn" : False
   "cnInfoLabelOn" : False
   "cnLineLabelPlacementMode" : "Constant"
   "cnLineLabelFont" : 21
   "cnLineLabelFontHeightF" : 0.011
   "cnLineLabelsOn" : True

   "tmXBMode" : "Manual"
   "tmXBTickSpacingF" : 3
   "tmXBMinorPerMajor" : 2
   "tmXBMajorLengthF" : -0.007
   "tmXBMinorLengthF" : -0.004
   "tmXBLabelFont" : 21
   "tmXBLabelFontHeightF" : 0.011
   "tmXBFormat" : "@*+^sg"

   "tmYLMode" : "Manual"
   "tmYLTickSpacingF" : 20
   "tmYLMinorPerMajor" : 1
   "tmYLMinorLengthF" : -0.004
   "tmYLMajorLengthF" : -0.007
   "tmYLLabelFont" : 21
   "tmYLLabelFontHeightF" : 0.011
   "tmYLFormat" : "@*+^sg"

   "vpXF" : 0.1
   "vpYF" : 0.97
   "vpWidthF" : 0.65
   "vpHeightF" : 0.18

   "trXMinF" : -12
   "trXMaxF" : 12
   "trYMinF" : -60
   "trYMaxF" : 60

;   "tiMainString" : "Correlation"
   "tiMainFont"  : 22
   "tiMainFontHeightF" : 0.017

   "tiXAxisString" : "lag (mon)"
   "tiXAxisFont"  : 21
   "tiXAxisFontHeightF" : 0.015

   "tiYAxisString" : "lat"
   "tiYAxisFont"  : 21
   "tiYAxisFontHeightF" : 0.015
  end create
  draw(cont1)
  delete(cont1)

  cont1 = create "contour" contourPlotClass wks
   "cnScalarFieldData" : field(1)
   "cnFillOn"        :  True           ; color fill
   "cnLevelSelectionMode" : "ExplicitLevels"
   "cnLevels" : cnlev2
   "cnFillColors" : cnfil2
   "cnLinesOn" : False
   "cnInfoLabelOn" : False
   "cnLineLabelPlacementMode" : "Constant"
   "cnLineLabelFont" : 21
   "cnLineLabelFontHeightF" : 0.011
   "cnLineLabelsOn" : True

   "tmXBMode" : "Manual"
   "tmXBTickSpacingF" : 0.01
   "tmXBMinorPerMajor" : 4
   "tmXBMajorLengthF" : -0.007
   "tmXBMinorLengthF" : -0.004
   "tmXBLabelFont" : 21
   "tmXBLabelFontHeightF" : 0.011
   "tmXBFormat" : "@*+^sg"

   "tmYLMode" : "Manual"
   "tmYLTickSpacingF" : 20
   "tmYLMinorPerMajor" : 1
   "tmYLMinorLengthF" : -0.004
   "tmYLMajorLengthF" : -0.007
   "tmYLLabelFont" : 21
   "tmYLLabelFontHeightF" : 0.011
   "tmYLFormat" : "@*+^sg"

   "vpXF" : 0.1
   "vpYF" : 0.69
   "vpWidthF" : 0.65
   "vpHeightF" : 0.18

   "trXMinF" : 0
   "trXMaxF" : 0.05
   "trYMinF" : -60
   "trYMaxF" : 60

;   "tiMainString" : "Coherence"
   "tiMainFont"  : 22
   "tiMainFontHeightF" : 0.017

   "tiXAxisString" : ""
   "tiXAxisFont"  : 21
   "tiXAxisFontHeightF" : 0.015

   "tiYAxisString" : "lat"
   "tiYAxisFont"  : 21
   "tiYAxisFontHeightF" : 0.015
  end create
  draw(cont1)

 setvalues cont1
  "cnScalarFieldData" : field(2)
  "cnLevels" : cnlev3
  "cnFillColors" : cnfil3
  "vpYF" : 0.47
  "tiXAxisString" : ""
;  "tiMainString" : "Phase"
 end setvalues
 draw(cont1)

 setvalues cont1
  "cnScalarFieldData" : field(3)
  "cnLevels" : cnlev4
  "cnFillColors" : cnfil4
  "vpYF" : 0.25
  "tiXAxisString" : "frequency (cycle/mon)"
;  "tiMainString" : "Amplitude"
 end setvalues 
 draw(cont1) 


; 6,5,12,6

   label = create "lab" labelBarClass wks
     "lbOrientation" : "vertical"
     "lbBoxCount" : 6
     "lbFillColors" : cnfil1
     "lbMonoFillPattern" : True
     "lbLabelsOn" : True
     "lbLabelStrings" : cnlev1
     "lbAutoManage" : False
     "lbLabelFont" : 21
     "lbLabelFontHeightF" : 0.010
     "lbLabelPosition" : "Bottom";,Right"
     "lbLabelAlignment" : "InteriorEdges"
     "lbLabelAngleF" : 0.
     "lbPerimOn" : False
     "lbTitleString" : ""
     "lbTitlePosition" : "Right";,"Top"
     "lbTitleDirection" : "Across"
     "lbTitleFont" : 22
     "lbTitleFontHeightF" : 0.010
     "vpXF" : 0.80
     "vpYF" : 0.97
     "vpWidthF"  : 0.06
     "vpHeightF" : 0.18
   end create
   draw(label)
   setvalues label
     "lbBoxCount" : 5
     "lbFillColors" : cnfil2
     "lbLabelStrings" : cnlev2
     "vpYF" : 0.69
   end setvalues
   draw(label)
   setvalues label 
     "lbBoxCount" : 12
     "lbFillColors" : cnfil3 
     "lbLabelStrings" : cnlev3 
     "vpYF" : 0.47 
   end setvalues 
   draw(label)
   setvalues label 
     "lbBoxCount" : 6 
     "lbFillColors" : cnfil4
     "lbLabelStrings" : cnlev4 
     "vpYF" : 0.25 
   end setvalues 
   draw(label)


;-----------------------------------------------

frame(wks)

end



