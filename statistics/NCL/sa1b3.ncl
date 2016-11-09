begin

 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

 setvalues NhlGetWorkspaceObjectId()
   "wsMaximumSize" : 55000000
 end setvalues

;--------------------------------------------

wks = new(2,graphic)

;---------------------------------------

 file1 = "res/sa1b3co2.nc"
 f1    = addfile(file1,"r")
 lon   = f1->lon
 freq1 = f1->freq
 p11   = f1->EQ_CO2
 p12   = f1->40N_CO2
 p13   = f1->40S_CO2

 file2 = "res/sa1b3pic.nc"
 f2    = addfile(file2,"r")
 freq2 = f2->freq 
 p21   = f2->EQ_pic 
 p22   = f2->40N_pic
 p23   = f2->40S_pic

;---------------------------------- cn level

 cnlev = new((/2,4/),float)
 cnfil = new(5,integer)

 do i=0,3
  cnlev(0,i) = 10.^(0.5 * i - 0)
  cnlev(1,i) = 10.^(0.5 * i - 0)
 end do
 cnfil = (/1,5,16,36,70/)*2    ;5

;------------------------------------------

 field = new((/2,3/),graphic)
 field(0,0) = create "w" scalarFieldClass defaultapp
   "sfDataArray" : p11
   "sfXArray" : lon
   "sfYArray" : freq1
 end create
 field(0,1) = create "w" scalarFieldClass defaultapp
   "sfDataArray" : p12
   "sfXArray" : lon 
   "sfYArray" : freq1
 end create
 field(0,2) = create "w" scalarFieldClass defaultapp
   "sfDataArray" : p13
   "sfXArray" : lon
   "sfYArray" : freq1 
 end create 
 field(1,0) = create "w" scalarFieldClass defaultapp
   "sfDataArray" : p21
   "sfXArray" : lon
   "sfYArray" : freq2 
 end create 
 field(1,1) = create "w" scalarFieldClass defaultapp
   "sfDataArray" : p22
   "sfXArray" : lon
   "sfYArray" : freq2
 end create
 field(1,2) = create "w" scalarFieldClass defaultapp
   "sfDataArray" : p23 
   "sfXArray" : lon 
   "sfYArray" : freq2  
 end create

 xf = (/0.1,0.1,0.55/)
 yf = (/0.95,0.50,0.50/)
 str = (/"EQ","40N","40S"/)
 xst = (/"","lon","lon"/)
 yst = (/"frequency (cycle/mon)","frequency (cycle/mon)",""/)

  fnn = (/"sa1b3co2.ps","sa1b3pic.ps"/)
  do j=0, 1

  wks(j) = create "1sttry" psWorkstationClass defaultapp
      "wkPSFileName" : fnn(j)
      "wkColorMap"   : "hotres"
  end create

  do i=0, 2

  cont1 = create "contour" contourPlotClass wks(j)
   "cnScalarFieldData" : field(j,i)
   "cnFillOn"        :  True           ; color fill
   "cnLevelSelectionMode" : "ExplicitLevels"
   "cnLevels" : cnlev(j,:)
   "cnFillColors" : cnfil
   "cnLinesOn" : False
   "cnInfoLabelOn" : False
   "cnLineLabelPlacementMode" : "Constant"
   "cnLineLabelFont" : 21
   "cnLineLabelFontHeightF" : 0.011
   "cnLineLabelsOn" : True

   "tmXBMode" : "Manual"
   "tmXBTickSpacingF" : 90
   "tmXBMinorPerMajor" : 2
   "tmXBMajorLengthF" : -0.007
   "tmXBMinorLengthF" : -0.004
   "tmXBLabelFont" : 21
   "tmXBLabelFontHeightF" : 0.011

   "tmYLMode" : "Manual"
   "tmYLTickSpacingF" : 0.1
   "tmYLMinorPerMajor" : 4
   "tmYLMinorLengthF" : -0.004
   "tmYLMajorLengthF" : -0.007
   "tmYLLabelFont" : 21
   "tmYLLabelFontHeightF" : 0.011
   "tmYLFormat" : "@*+^sg"

   "vpXF" : xf(i)
   "vpYF" : yf(i)
   "vpWidthF" : 0.35
   "vpHeightF" : 0.35

   "trXMinF" : 0
   "trXMaxF" : 358
   "trYMinF" : 0
   "trYMaxF" : 0.3

   "tiMainString" : "PSD of ua_"+str(i)
   "tiMainFont"  : 22
   "tiMainFontHeightF" : 0.017

   "tiXAxisString" : xst(i)
   "tiXAxisFont"  : 21
   "tiXAxisFontHeightF" : 0.015

   "tiYAxisString" : yst(i)
   "tiYAxisFont"  : 21
   "tiYAxisFontHeightF" : 0.015
  end create
  draw(cont1)
  delete(cont1)

  end do

   labels = new(4,float)
   labels(:) = log10(cnlev(j,:))
   label = create "lab" labelBarClass wks(j)
     "lbOrientation" : "horizontal"
     "lbBoxCount" : 5
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
     "vpXF" : 0.25
     "vpYF" : 0.10
     "vpWidthF" : 0.35
     "vpHeightF" : 0.07
   end create
   draw(label)

   delete(label)
   delete(labels)

   frame(wks(j))

  end do


;-------------------------------------------------

end



