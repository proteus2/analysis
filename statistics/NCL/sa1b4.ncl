begin

 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

 setvalues NhlGetWorkspaceObjectId()
   "wsMaximumSize" : 55000000
 end setvalues

;--------------------------------------------

wks = create "1sttry" psWorkstationClass defaultapp
      "wkPSFileName" : "sa1b4.ps"
      "wkColorMap"   : "hotres"
end create

;---------------------------------------

 file1 = "res/sa1b4.nc"
 f1    = addfile(file1,"r")
 lat   = f1->lat_CO2
 freq1 = f1->freq_CO2
 freq2 = f1->freq_pic
 p10   = f1->PS_CO2
 p20   = f1->PS_pic

 p1 = new((/dimsizes(lat),dimsizes(freq1)/),float)
 p2 = new((/dimsizes(lat),dimsizes(freq2)/),float)
 do i=0, dimsizes(lat)-1
   p1(i,:) = p10(:,i)
   p2(i,:) = p20(:,i)
 end do

;---------------------------------- cn level

 cnlev = new(7,float)
 cnfil = new(8,integer)

 do i=0,6
  cnlev(i) = 10.^(0.5 * i - 3)
 end do
 cnfil = (/1,4,9,16,36,64,81,100/)*2

;------------------------------------------

 field = new(2,graphic)
 field(0) = create "w" scalarFieldClass defaultapp
   "sfDataArray" : p1
   "sfYArray" : lat
   "sfXArray" : freq1
 end create
 field(1) = create "w" scalarFieldClass defaultapp
   "sfDataArray" : p2
   "sfYArray" : lat 
   "sfXArray" : freq2
 end create

 yf = (/0.95,0.50/)
 str = (/"CO:B:2:N:","pictl"/)

  do i=0, 1
  cont1 = create "contour" contourPlotClass wks
   "cnScalarFieldData" : field(i)
   "cnFillOn"        :  True           ; color fill
   "cnLevelSelectionMode" : "ExplicitLevels"
   "cnLevels" : cnlev
   "cnFillColors" : cnfil
   "cnLinesOn" : False
   "cnInfoLabelOn" : False
   "cnLineLabelPlacementMode" : "Constant"
   "cnLineLabelFont" : 21
   "cnLineLabelFontHeightF" : 0.011
   "cnLineLabelsOn" : True

   "tmXBMode" : "Manual"
   "tmXBTickSpacingF" : 0.1
   "tmXBMinorPerMajor" : 4
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

   "trXMinF" : 0
   "trXMaxF" : 0.5
   "trYMinF" : -80
   "trYMaxF" : 80

   "tiMainString" : "PSD of Ta_"+str(i)
   "tiMainFont"  : 22
   "tiMainFontHeightF" : 0.017

   "tiXAxisString" : "frequency (cycle/mon)"
   "tiXAxisFont"  : 21
   "tiXAxisFontHeightF" : 0.015

   "tiYAxisString" : "lat"
   "tiYAxisFont"  : 21
   "tiYAxisFontHeightF" : 0.015
  end create
  draw(cont1)
  delete(cont1)
  end do


   labels = log10(cnlev)
   label = create "lab" labelBarClass wks
     "lbOrientation" : "horizontal"
     "lbBoxCount" : 8
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
     "lbTitleString" : "Log:B:10:N:(PSD)"
     "lbTitlePosition" : "Right"
     "lbTitleDirection" : "Across"
     "lbTitleFont" : 21
     "lbTitleFontHeightF" : 0.015
     "vpXF" : 0.50
     "vpYF" : 0.10
     "vpWidthF" : 0.35
     "vpHeightF" : 0.08
   end create
   draw(label)

;-------------------------------------------------

frame(wks)

end



