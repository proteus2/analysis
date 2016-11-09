begin

 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

;--------------------------------------------

wks = create "1sttry" psWorkstationClass defaultapp
      "wkPSFileName" : "sa2b1.ps"
;      "wkColorMap"   : "ViBlGrWhYeOrRe"
end create

;---------------------------------------

 file1 = "sa2b1pic.nc"
 f1    = addfile(file1,"r")
 lag   = f1->lag
 freq  = f1->freq
 cor1  = f1->COR
 coh1  = f1->COH
 phs1  = f1->PHASE
 amp1  = f1->AMP

;------------------------------------------

 data0 = new(4,graphic)
 data0(0)= create "00" coordArraysClass defaultapp
   "caXArray" : lag
   "caYArray" : cor1
 end create
 data0(1)= create "00" coordArraysClass defaultapp
   "caXArray" : freq
   "caYArray" : coh1 
 end create 
 data0(2)= create "00" coordArraysClass defaultapp
   "caXArray" : freq
   "caYArray" : phs1 
 end create 
 data0(3)= create "00" coordArraysClass defaultapp
   "caXArray" : freq
   "caYArray" : amp1
 end create 


 plot0 = create "plot0" xyPlotClass wks
  "xyCoordData" : data0(0)

  "trYMaxF" : 0.5
  "trYMinF" : -0.5
  "trXMaxF" : 300
  "trXMinF" : -300

  "tmXBMode" : "MANUAL"
  "tmXBTickSpacingF" : 100.
  "tmXBMinorPerMajor" : 4
  "tmXBLabelFont" : 21
  "tmXBLabelFontHeightF" : 0.011
  "tmXBMajorLengthF" :  -0.007
  "tmXBMinorLengthF" :  -0.004
  "tmXBFormat" : "@*+^sg"

  "tmYLMinorLengthF" : -0.004
  "tmYLMajorLengthF" : -0.008
  "tmYLMode" : "MANUAL"
  "tmYLTickSpacingF" : 0.5
  "tmYLMinorPerMajor" : 4
  "tmYLLabelFont" : 21
  "tmYLLabelFontHeightF" : 0.011
  "tmYLFormat" : "@*+^sg"

  "tiMainFont" : 22
  "tiMainFontHeightF" : 0.020
  "tiMainString" : ""
  "tiXAxisFont" : 21
  "tiXAxisFontHeightF" : 0.017
  "tiXAxisString" : "lag"
  "tiYAxisFont" : 21
  "tiYAxisFontHeightF" : 0.017
  "tiYAxisString" : "Correlation"

  "vpYF" : 0.97
  "vpXF" : 0.1
  "vpHeightF" : 0.18
  "vpWidthF" : 0.6
 end create
 draw(plot0)

 setvalues plot0
  "xyCoordData" : data0(1)
  "trYMaxF" : 1.2
  "trYMinF" : 0
  "trXMaxF" : 0.05
  "trXMinF" : 0
  "vpYF" : 0.69
  "tmXBTickSpacingF" : 0.01
  "tmXBMinorPerMajor" : 4
  "tmYLTickSpacingF" : 0.2
  "tmYLMinorPerMajor" : 1
  "tiXAxisString" : ""
  "tiYAxisString" : "Coherence"
 end setvalues
 draw(plot0)

 setvalues plot0
  "xyCoordData" : data0(2)
  "trYMaxF" : 180
  "trYMinF" : -180
  "trXMaxF" : 0.05
  "trXMinF" : 0
  "vpYF" : 0.47
  "tmYLTickSpacingF" : 90
  "tmYLMinorPerMajor" : 2
  "tiXAxisString" : ""
  "tiYAxisString" : "Phase"
 end setvalues
 draw(plot0)

 setvalues plot0
  "xyCoordData" : data0(3)
  "trYMaxF" : 30.
  "trYMinF" : 0 
  "trXMaxF" : 0.05
  "trXMinF" : 0 
  "vpYF" : 0.25
  "tmYLTickSpacingF" : 10.
  "tmYLMinorPerMajor" : 1
  "tiXAxisString" : "frequency"
  "tiYAxisString" : "Amplitude"
 end setvalues 
 draw(plot0) 


;-----------------------------------------------

frame(wks)

end



