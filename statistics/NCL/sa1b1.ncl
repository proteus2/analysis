begin

 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

;--------------------------------------------

wks = create "1sttry" psWorkstationClass defaultapp
      "wkPSFileName" : "sa1b1.ps"
      "wkColorMap"   : "BlWhRe"
end create

;---------------------------------------

 file1 = "res/sa1b1.nc"
 f1   = addfile(file1,"r")
 lon  = f1->lon
 lat  = f1->lat
 T1   = f1->TmeanCO2
 T2   = f1->Tmeanpic

;---------------------------------- cn level

 cnlev = new(8,float)
 cnfil = new(9,integer)

 do i=0,7
  cnlev(i) = 10 * i + 230
 end do
 cnfil = (/-18,-12,-7,-4,-2,0,4,7,12/)
 cnfil = cnfil + 52

 cnlev2 = new(15,float)
 do i=0, 14
   cnlev2(i) = 5 * i + 230
 end do
;------------------------------------------

 field = new(2,graphic)
 field(0) = create "w" scalarFieldClass defaultapp
   "sfDataArray" : T1
   "sfYArray" : lat
   "sfXArray" : lon
 end create
 field(1) = create "w" scalarFieldClass defaultapp
   "sfDataArray" : T2
   "sfYArray" : lat 
   "sfXArray" : lon 
 end create

 yf = (/0.9,0.45/)
 str = (/"CO:B:2:N:","pictl"/)

  do i=0, 1
  cont1 = create "contour" contourPlotClass wks
   "cnScalarFieldData" : field(i)
   "cnFillOn"        :  True           ; color fill
   "cnLevelSelectionMode" : "ExplicitLevels"
   "cnLevels" : cnlev
   "cnFillColors" : cnfil
   "cnLinesOn" : False
;   "cnLineColor" : 254
   "cnInfoLabelOn" : False
   "cnLineLabelPlacementMode" : "Constant"
   "cnLineLabelFont" : 21
   "cnLineLabelFontHeightF" : 0.011
   "cnLineLabelsOn" : False ;True

   "tmXBMode" : "Manual"
   "tmXBTickSpacingF" : 60
   "tmXBMinorPerMajor" : 5
   "tmXBMajorLengthF" : -0.007
   "tmXBMinorLengthF" : -0.004
   "tmXBLabelFont" : 21
   "tmXBLabelFontHeightF" : 0.011

   "tmYLMode" : "Manual"
   "tmYLTickSpacingF" : 20
   "tmYLMinorPerMajor" : 3
   "tmYLMinorLengthF" : -0.004
   "tmYLMajorLengthF" : -0.007
   "tmYLLabelFont" : 21
   "tmYLLabelFontHeightF" : 0.011
   "tmYLFormat" : "@*+^sg"

   "vpXF" : 0.15
   "vpYF" : yf(i)
   "vpWidthF" : 0.7
   "vpHeightF" : 0.30

   "trXMinF" : 0.
   "trXMaxF" : 358.125
   "trYMinF" : -80
   "trYMaxF" : 80

   "tiMainString" : "mean Ta_"+str(i)
   "tiMainFont"  : 22
   "tiMainFontHeightF" : 0.017

   "tiXAxisString" : "lon"
   "tiXAxisFont"  : 21
   "tiXAxisFontHeightF" : 0.015

   "tiYAxisString" : "lat"
   "tiYAxisFont"  : 21
   "tiYAxisFontHeightF" : 0.015
  end create
  draw(cont1)
  setvalues cont1
   "cnFillOn"  : False
   "cnLinesOn" : True ;False
   "cnLevels" : cnlev2
  end setvalues
  draw(cont1)

  delete(cont1)
  end do

   labels = cnlev
   label = create "lab" labelBarClass wks
     "lbOrientation" : "horizontal"
     "lbBoxCount" : 9
     "lbFillColors" : cnfil
     "lbMonoFillPattern" : True
     "lbLabelsOn" : True
     "lbLabelStrings" : labels
     "lbAutoManage" : False
     "lbLabelFont" : 21
     "lbLabelFontHeightF" : 0.012
     "lbLabelPosition" : "Bottom";,Right"
     "lbLabelAlignment" : "InteriorEdges"
     "lbLabelAngleF" : 0.
     "lbPerimOn" : False
     "lbTitleString" : "K"
     "lbTitlePosition" : "Right"
     "lbTitleDirection" : "Across"
     "lbTitleFont" : 21
     "lbTitleFontHeightF" : 0.015
     "vpXF" : 0.50
     "vpYF" : 0.08
     "vpWidthF" : 0.40
     "vpHeightF" : 0.06
   end create
   draw(label)

;-------------------------------------------------

frame(wks)

end



