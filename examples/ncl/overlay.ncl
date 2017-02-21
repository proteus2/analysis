begin

 setvalues NhlGetWorkspaceObjectId()
   "wsMaximumSize" : 55000000
 end setvalues

;---------------------------------------------- wks

 load "/usr/local/ncarg/lib/ncarg/nclscripts/csm/gsn_code.ncl"
 load "/usr/local/ncarg/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

 wks = create "contour" psWorkstationClass defaultapp
     "wkPSFileName" : "prinsec.ps"
     "wkColorMap" : "hotres"
 end create


; tkz ==============================================================================

;-------------------------------- reading

 fnam10 = "nc/tkz_x250-600.nc"
 f10 = addfile(fnam10,"r")
 d1 = f10->PSD
 k1 = f10->k
 z1 = f10->z
 t1 = f10->t

 fnam11= "nc/tkz_x450-800.nc"
 f11 = addfile(fnam11,"r")
 d2 = f11->PSD
 k2 = f11->k
 z2 = f11->z
 t2 = f11->t

;---------------------------------- cn level

 cnlev = new(5,float)
 cnfil = new(6,integer)

 do i=0, 4
;  cnlev(i) = 2.^(i - 3.) ;7.*(i+1)
  cnlev(i) = 10.^(i*0.5 + 1.5)
 end do
 cnfil = (/2,5,20,36,64,100/)*2

;--------------------------------------------- field

  field1 = create "w" scalarFieldClass defaultapp
    "sfDataArray" : d1(7,:,:) / (2.*3.141593)
    "sfXArray" : t1 / 60.
    "sfYArray" : k1*10000. * (2.*3.141593)
  end create

  field2 = create "w" scalarFieldClass defaultapp
    "sfDataArray" : d1(3,:,:) / (2.*3.141593)
    "sfXArray" : t1 / 60.
    "sfYArray" : k1*10000. * (2.*3.141593)
  end create 
 
  field3 = create "w" scalarFieldClass defaultapp
    "sfDataArray" : d2(1,:,:) / (2.*3.141593)
    "sfXArray" : t2 / 60.
    "sfYArray" : k2*10000. * (2.*3.141593)
  end create

;--------------------------------------------- contour
; wprt1  
;----------------------------------------------------
  
  cont1 = create "contour" contourPlotClass wks

   "cnScalarFieldData" : field1
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
   "tmXBTickSpacingF" : 1
   "tmXBMinorPerMajor" : 1
   "tmXBMinorLengthF" : -0.004
   "tmXBMajorLengthF" : -0.007
   "tmXBLabelFont" : 21
   "tmXBLabelFontHeightF" : 0.012
   "tmXBFormat" : "@*+^sg"

   "tmYLMode" : "Manual"
   "tmYLTickSpacingF" : 1
   "tmYLMinorPerMajor" : 1
   "tmYLMinorLengthF" : -0.004
   "tmYLMajorLengthF" : -0.007
   "tmYLLabelFont" : 21
   "tmYLLabelFontHeightF" : 0.012

   "vpXF" : 0.1
   "vpYF" : 0.91
   "vpWidthF" : 0.25
   "vpHeightF" : 0.25

   "trXMinF" : 1;k1(0)*100000 
   "trXMaxF" : 4
   "trYMinF" : 0;15
   "trYMaxF" : 4;105

   "tiMainString" : "z = 90 km"
   "tiMainFont"  : 21
   "tiMainFontHeightF" : 0.018

   "tiXAxisString" : "t (hr)"
   "tiXAxisFont"  : 21
   "tiXAxisFontHeightF" : 0.016
   "tiYAxisString" : "k (10:S:-4:N: rad/m)"
   "tiYAxisFont"  : 21
   "tiYAxisFontHeightF" : 0.016
  end create
   

;----------------------------------------------------
;  draw(cont1)
;---------------------------------------------- text

setvalues cont1
   "cnScalarFieldData" : field2

   "vpXF" : 0.40
   "vpYF" : 0.91
   "tiMainString" : "z = 50 km"
   "tiYAxisString" : ""
end setvalues
;draw(cont1)
 
setvalues cont1
   "cnScalarFieldData" : field3

   "vpXF" : 0.35 ;0.70
   "vpYF" : 0.91
   "tiMainString" : "PSD";"z = 40 km"
   "tiYAxisString" : "k (10:S:-4:N: rad/m)";""

   "tmYUseLeft"  : False
   "tmYROn"      : True
   "tmYRLabelsOn": True
   "tmYRMode" : "Explicit"
   "tmYRValues" : 2.*3.141593 / ((/10,20,50,100/)*1000.) * 10000.
   "tmYRLabels" : (/10,20,50,100/)
   "tmYRMinorValues" : 2.*3.141593 / ((/30,40,60,70,80,90/)*1000.) * 10000.
   "tmYRMinorLengthF" : -0.005
   "tmYRMajorLengthF" : -0.008
   "tmYRLabelFont" : 21
   "tmYRLabelFontHeightF" : 0.012
   "tmYRFormat" : "@*+^sg"

end setvalues
draw(cont1)

;frame(wks)


   labels = log10(cnlev)
   label = create "lab" labelBarClass wks
     "lbOrientation" : "horizontal"
     "lbBoxCount" : 6 
     "lbFillColors" : cnfil
     "lbMonoFillPattern" : True
     "lbLabelsOn" : True
     "lbLabelStrings" : labels
     "lbAutoManage" : False
     "lbLabelFont" : 21
     "lbLabelFontHeightF" : 0.013
     "lbLabelPosition" : "Bottom";,Right"
     "lbLabelAlignment" : "InteriorEdges"
     "lbLabelAngleF" : 0.
     "lbPerimOn" : False
     "lbTitleString" : "log:B:10:N: (PSD)"
     "lbTitlePosition" : "Right" ;"Top"
     "lbTitleDirection" : "Across"
     "lbTitleFont" : 21
     "lbTitleFontHeightF" : 0.015
     "vpXF" : 0.35
     "vpYF" : 0.58
     "vpWidthF" : 0.30
     "vpHeightF" : 0.04
   end create
   draw(label)


; to ===============================================================================

; read ---------------------------------
 fnam2 = new(3,string)

 fnam2 = (/"nc/ot_x550z30.nc","nc/ot_x550z30.nc","nc/ot_x550z30.nc"/)

xf = (/0.25,0.35,0.7/)
mstr = (/"z = 80 km","Wavelet PSD","z = 40 km","z = 40 km"/)
ystr = (/"",":F33:w:F21: (rad/s)","",""/)
do nn=0, 1;2

 f2    = addfile(fnam2(nn),"r")
 ps2   = f2->PS
 xout  = f2->Xout
 peri  = f2->FScale

 freq  = (2.*3.141593) / (peri*60.)
 ps2   = ps2 * 60. / (2.*3.141593)

 psd = new((/dimsizes(peri),dimsizes(xout)-1/),float)
 x   = new(dimsizes(xout)-1,float)
 do i=0, dimsizes(xout)-2
   x(i) = xout(i)
   psd(:,i) = ps2(i,:)
 end do

;---------------------------------- cn level

 cnlev2 = new(5,float)
 cnfil2 = new(6,integer)
 do i=0,4
;  cnlev2(i) = 2.^(i-6 -2)
  cnlev2(i) = 10.^(i*0.5 - 0.0)
 end do
 cnfil2 = cnfil

;------------------------------------------

 field = create "w" scalarFieldClass defaultapp
   "sfDataArray" : psd
   "sfYArray" : freq
   "sfXArray" : x / 60.
 end create

  cont1 = create "contour" contourPlotClass wks
   "cnScalarFieldData" : field
   "cnFillOn"        :  True           ; color fill
   "cnLevelSelectionMode" : "ExplicitLevels"
   "cnLevels" : cnlev2
   "cnFillColors" : cnfil2
   "cnLinesOn" : False
   "cnInfoLabelOn" : False
   "cnLineLabelPlacementMode" : "Constant"
   "cnLineLabelFont" : 21
   "cnLineLabelFontHeightF" : 0.011
   "cnLineLabelsOn" : False

   "tmXBMode" : "Manual"
   "tmXBTickSpacingF" : 1
   "tmXBMinorPerMajor" : 1
   "tmXBMajorLengthF" : -0.007
   "tmXBMinorLengthF" : -0.004
   "tmXBLabelFont" : 21
   "tmXBLabelFontHeightF" : 0.012
   "tmXBFormat" : "@*+^sg"

   "tmYLMode" : "Manual"
   "tmYLTickSpacingF" : 0.01
   "tmYLMinorPerMajor" : 1
   "tmYLMajorLengthF" : -0.007
   "tmYLMinorLengthF" : -0.004
   "tmYLLabelFont" : 21
   "tmYLLabelFontHeightF" : 0.012
   "tmYLFormat" : "@*+^sg"

   "vpXF" : xf(nn);0.4
   "vpYF" : 0.40 ;0.37
   "vpWidthF" : 0.25;0.55
   "vpHeightF" : 0.25

   "trXMinF" : 1
   "trXMaxF" : 4
   "trYMinF" : 1
   "trYMaxF" : 0.05

   "tiMainString" : mstr(nn)
   "tiMainFont"  : 21
   "tiMainFontHeightF" : 0.018

   "tiXAxisString" : "t (hr)"
   "tiXAxisFont"  : 21
   "tiXAxisFontHeightF" : 0.016

   "tiYAxisString" : ystr(nn)
   "tiYAxisFont"  : 21
   "tiYAxisFontHeightF" : 0.016
  end create
  if (nn .eq. 2) then
    setvalues cont1
      "tmYUseLeft"  : False
      "tmYROn"      : True
      "tmYRLabelsOn": True
      "tmYRMode" : "Explicit"
      "tmYRValues" : 2.*3.141593/60. / (/1,5,10,60/)
      "tmYRLabels" : (/1,5,10,60/) ;2.^(/0,1,2,3,4,5,6,7,8,9/)
      "tmYRMinorValues" : 2.*3.141593/60. / (/2,3,4,6,7,8,9,20,30,40,50/)
      "tmYRMinorLengthF" : -0.005
      "tmYRMajorLengthF" : -0.008
      "tmYRLabelFont" : 21
      "tmYRLabelFontHeightF" : 0.012
      "tmYRFormat" : "@*+^sg"
    end setvalues
  end if
  plotover = create "plot" logLinPlotClass wks
    "vpXF" : xf(nn);0.4
    "vpYF" : 0.40 ;0.37
    "vpWidthF" : 0.25;0.55
    "vpHeightF" : 0.25

    "trXMinF" : 1
    "trXMaxF" : 4
    "trYMinF" : 0
    "trYMaxF" : 0.05
  end create

  overlay(plotover,cont1)
if (nn .eq. 1) then
  draw(plotover)
end if
  delete(plotover)
  delete(cont1)

end do

;-------------------------------------------------------------


  setvalues label
    "lbFillColors" : cnfil2
    "lbLabelStrings" : log10(cnlev2)
    "vpXF" : 0.35
    "vpYF" : 0.07 ;0.04
  end setvalues
  draw(label)


; ===============================================================================

  txt = create "text" textItemClass wks
    "txPosXF"               :  0.37
    "txPosYF"               :  0.98
    "txFontHeightF"         :  0.018
    "txFont"                :  21
    "txString"              :  "x = 250 ~ 600 km"
  end create
;  draw(txt)
  setvalues txt
    "txPosXF"               :  0.82
    "txString"              :  "x = 450 ~ 800 km"
  end setvalues
;  draw(txt)


  setvalues txt
    "txPosXF"               :  0.375 ;0.225
    "txPosYF"               :  0.44
    "txFontHeightF"         :  0.018
    "txFont"                :  21
    "txString"              :  "x = 440 km"
  end setvalues
;  draw(txt)
  setvalues txt
    "txPosXF"               :  0.825
    "txString"              :  "x = 550 km"
  end setvalues
;  draw(txt)


  frame(wks)

 end



