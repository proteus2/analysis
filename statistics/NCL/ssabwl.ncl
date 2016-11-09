begin

 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

;--------------------------------------------

wks = create "1sttry" psWorkstationClass defaultapp
      "wkPSFileName" : "ssabwl.ps"
      "wkColorMap"   : "hotres"
end create

;---------------------------------------

 file1 = "res/ssabwlpic.nc"
 f1    = addfile(file1,"r")
 ps1   = f1->PS
 xout  = f1->Xout
 peri  = f1->FScale

 file3 = "res/wabpic.nc"
 f3    = addfile(file3,"r")
 ps3   = f3->PS

 file4 = "res/wabpicsig.nc"
 f4    = addfile(file4,"r")
 sig3  = f4->Sig
 ps_theor3 = f4->PS_theor
 sig_glob3 = f4->Sig_glob


 psd = new((/3,dimsizes(peri),dimsizes(xout)-1/),float)
 x   = new(dimsizes(xout)-1,float)
 do i=0, dimsizes(xout)-2
   x(i) = xout(i)
   psd(0,:,i) = ps1(i,:)
   psd(1,:,i) = ps3(i,:)
   psd(2,:,i) = psd(1,:,i)/sig3(:)
 end do


;---------------------------------- cn level

 cnlev = new(6,float)
 cnfil = new(7,integer)
 do i=0,5
  cnlev(i) = 2.^(i-5 -1)
 end do
 cnfil = (/1,4,9,16,36,64,90/)*2


;------------------------------------------

 field = new(3,graphic)
 do i=0, 2
 field(i) = create "w" scalarFieldClass defaultapp
   "sfDataArray" : psd(i,:,:)
   "sfYArray" : peri
   "sfXArray" : x
 end create
 end do


 yf = (/0.95,0.50,0.50/)
 str = (/"Wavelet PS of reconstructed Nino3.4 index (mode1-5)", \
         "Wavelet PS of Nino3.4 index","Wavelet PS of Nino3.4 index"/)

 do i=0, 2
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
   "cnLineLabelsOn" : False

   "tmXBMode" : "Manual"
   "tmXBTickSpacingF" : 20
   "tmXBMinorPerMajor" : 3
   "tmXBMajorLengthF" : -0.007
   "tmXBMinorLengthF" : -0.004
   "tmXBLabelFont" : 21
   "tmXBLabelFontHeightF" : 0.011
   "tmXBFormat" : "@*+^sg"

   "tmYLMode" : "Explicit"
   "tmYLValues" : 2.^(/0,1,2,3,4,5,6,7,8,9/)
   "tmYLLabels" : 2.^(/0,1,2,3,4,5,6,7,8,9/)
   "tmYLMinorLengthF" : -0.004
   "tmYLMajorLengthF" : -0.007
   "tmYLLabelFont" : 21
   "tmYLLabelFontHeightF" : 0.011
   "tmYLFormat" : "@*+^sg"

   "vpXF" : 0.1
   "vpYF" : yf(i)
   "vpWidthF" : 0.65
   "vpHeightF" : 0.30

   "trXMinF" : 0
   "trXMaxF" : max(x)
   "trYMinF" : 1
   "trYMaxF" : 32
   "trYReverse" : "True"

   "tiMainString" : str(i)
   "tiMainFont"  : 22
   "tiMainFontHeightF" : 0.017

   "tiXAxisString" : "time (yr)"
   "tiXAxisFont"  : 21
   "tiXAxisFontHeightF" : 0.015

   "tiYAxisString" : "period (yr)"
   "tiYAxisFont"  : 21
   "tiYAxisFontHeightF" : 0.015
  end create
  if (i .eq. 2) then
    setvalues cont1
      "cnFillOn"     : False
      "cnLevels"     : 1.
      "cnLinesOn"    : True
    end setvalues
  end if
  draw(cont1)
  delete(cont1)
  end do




   label = create "lab" labelBarClass wks
     "lbOrientation" : "vertical"
     "lbBoxCount" : 7
     "lbFillColors" : cnfil
     "lbMonoFillPattern" : True
     "lbLabelsOn" : True
     "lbLabelStrings" : log(cnlev)/log(2.)
     "lbAutoManage" : False
     "lbLabelFont" : 21
     "lbLabelFontHeightF" : 0.010
     "lbLabelPosition" : "Bottom";,Right"
     "lbLabelAlignment" : "InteriorEdges"
     "lbLabelAngleF" : 0.
     "lbPerimOn" : False
     "lbTitleString" : "log:B:2:N: (PSD)"
     "lbTitlePosition" : "Bottom";,"Top"
     "lbTitleDirection" : "Across"
     "lbTitleFont" : 22
     "lbTitleFontHeightF" : 0.010
     "vpXF" : 0.80
     "vpYF" : 0.95
     "vpWidthF"  : 0.06
     "vpHeightF" : 0.25
   end create
   draw(label)


;-----------------------------------------------

frame(wks)

end



