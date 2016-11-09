begin

 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

;--------------------------------------------

wks = create "1sttry" psWorkstationClass defaultapp
      "wkPSFileName" : "waa.ps"
      "wkColorMap"   : "hotres"
end create

;---------------------------------------

 file1 = "res/waa1_0.nc"
 f1    = addfile(file1,"r")
 ps1   = f1->PS
 xout  = f1->Xout
 peri  = f1->FScale

 file2 = "res/waa1_wn.nc"
 f2    = addfile(file2,"r")
 ps2   = f2->PS

 file3 = "res/waa2.nc"
 f3    = addfile(file3,"r") 
 ps3   = f3->PS

 psd = new((/3,dimsizes(peri),dimsizes(xout)/),float)
 x   = new(dimsizes(xout),float)
 do i=1, dimsizes(xout)-1
   x(i) = xout(i-1)
   psd(0,:,i) = ps1(i-1,:)
   psd(1,:,i) = ps2(i-1,:)
   psd(2,:,i) = ps3(i-1,:)
 end do
 x(0) = 0.
 psd(0,:,0) = ps1(dimsizes(xout)-1,:)
 psd(1,:,0) = ps2(dimsizes(xout)-1,:)
 psd(2,:,0) = ps3(dimsizes(xout)-1,:)

;---------------------------------- cn level

 cnlev1 = new(7,float)
 cnlev2 = new(7,float)
 cnfil = new(8,integer)
 do i=0,6
  cnlev1(i) = 10.^(0.5 * i + 0.0)
  cnlev2(i) = 10.^(0.5 * i - 2.0)
 end do
 cnfil = (/1,4,9,16,36,64,81,100/)*2


;------------------------------------------

 field = new(3,graphic)
 do i=0, 2
 field(i) = create "w" scalarFieldClass defaultapp
   "sfDataArray" : psd(i,:,:)
   "sfYArray" : peri
   "sfXArray" : x
 end create
 end do

 yf = (/0.95,0.63,0.31/)
 do i=0, 2
  cont1 = create "contour" contourPlotClass wks
   "cnScalarFieldData" : field(i)
   "cnFillOn"        :  True           ; color fill
   "cnLevelSelectionMode" : "ExplicitLevels"
   "cnLevels" : cnlev1
   "cnFillColors" : cnfil
   "cnLinesOn" : False
   "cnInfoLabelOn" : False
   "cnLineLabelPlacementMode" : "Constant"
   "cnLineLabelFont" : 21
   "cnLineLabelFontHeightF" : 0.011
   "cnLineLabelsOn" : True

   "tmXBMode" : "Manual"
   "tmXBTickSpacingF" : 200
   "tmXBMinorPerMajor" : 3
   "tmXBMajorLengthF" : -0.007
   "tmXBMinorLengthF" : -0.004
   "tmXBLabelFont" : 21
   "tmXBLabelFontHeightF" : 0.011
   "tmXBFormat" : "@*+^sg"

   "tmYLMode" : "Explicit"
   "tmYLValues" : 2.^(/1,2,3,4,5,6,7,8,9/)
   "tmYLLabels" : 2.^(/1,2,3,4,5,6,7,8,9/)
   "tmYLMinorLengthF" : -0.004
   "tmYLMajorLengthF" : -0.007
   "tmYLLabelFont" : 21
   "tmYLLabelFontHeightF" : 0.011
   "tmYLFormat" : "@*+^sg"

   "vpXF" : 0.1
   "vpYF" : yf(i)
   "vpWidthF" : 0.65
   "vpHeightF" : 0.25

   "trXMinF" : 0
   "trXMaxF" : 1000
   "trYMinF" : 2
   "trYMaxF" : 1000
   "trYReverse" : "True"

;   "tiMainString" : "Correlation"
   "tiMainFont"  : 22
   "tiMainFontHeightF" : 0.017

   "tiXAxisString" : "time"
   "tiXAxisFont"  : 21
   "tiXAxisFontHeightF" : 0.015

   "tiYAxisString" : "period"
   "tiYAxisFont"  : 21
   "tiYAxisFontHeightF" : 0.015
  end create
  if (i .eq. 2) then
    setvalues cont1
      "cnLevels" : cnlev2
    end setvalues
  end if
  draw(cont1)
  delete(cont1)
  end do



   label = create "lab" labelBarClass wks
     "lbOrientation" : "vertical"
     "lbBoxCount" : 8
     "lbFillColors" : cnfil
     "lbMonoFillPattern" : True
     "lbLabelsOn" : True
     "lbLabelStrings" : log10(cnlev1)
     "lbAutoManage" : False
     "lbLabelFont" : 21
     "lbLabelFontHeightF" : 0.010
     "lbLabelPosition" : "Bottom";,Right"
     "lbLabelAlignment" : "InteriorEdges"
     "lbLabelAngleF" : 0.
     "lbPerimOn" : False
     "lbTitleString" : "log:B:10:N:(PSD)"
     "lbTitlePosition" : "Bottom";,"Top"
     "lbTitleDirection" : "Across"
     "lbTitleFont" : 22
     "lbTitleFontHeightF" : 0.010
     "vpXF" : 0.80
     "vpYF" : 0.75
     "vpWidthF"  : 0.06
     "vpHeightF" : 0.25
   end create
   draw(label)

   setvalues label
     "lbLabelStrings" : log10(cnlev2)
     "vpYF" : 0.30
   end setvalues
   draw(label)

;-----------------------------------------------

frame(wks)

end



