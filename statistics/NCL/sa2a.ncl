begin

 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

;--------------------------------------------

wks = create "1sttry" psWorkstationClass defaultapp
      "wkPSFileName" : "sa2a.ps"
;      "wkColorMap"   : "ViBlGrWhYeOrRe"
end create

;---------------------------------------

 file1 = "sa2a_0.nc"
 f1    = addfile(file1,"r")
 lag   = f1->lag
 freq  = f1->freq
 cor1  = f1->COR
 coh1  = f1->COH
 phs1  = f1->PHASE
 amp1  = f1->AMP

 file2 = "sa2a_wn.nc"
 f2    = addfile(file2,"r")
 cor1w = f2->COR
 coh1w = f2->COH
 phs1w = f2->PHASE
 amp1w = f2->AMP 

 file3 = "sa2a_phs.nc"
 f3    = addfile(file3,"r")
 cor2  = f3->COR
 coh2  = f3->COH
 phs2  = f3->PHASE
 amp2  = f3->AMP

;------------------------------------------

 data0 = new((/2,4/),graphic)
 dataw = new(4,graphic)
 data0(0,0)= create "00" coordArraysClass defaultapp
   "caXArray" : lag
   "caYArray" : cor1
 end create
 data0(0,1)= create "00" coordArraysClass defaultapp
   "caXArray" : freq
   "caYArray" : coh1 
 end create 
 data0(0,2)= create "00" coordArraysClass defaultapp
   "caXArray" : freq
   "caYArray" : phs1 
 end create 
 data0(0,3)= create "00" coordArraysClass defaultapp
   "caXArray" : freq
   "caYArray" : amp1/1000. 
 end create 
 data0(1,0)= create "00" coordArraysClass defaultapp
   "caXArray" : lag 
   "caYArray" : cor2
 end create 
 data0(1,1)= create "00" coordArraysClass defaultapp
   "caXArray" : freq
   "caYArray" : coh2
 end create  
 data0(1,2)= create "00" coordArraysClass defaultapp
   "caXArray" : freq
   "caYArray" : phs2
 end create  
 data0(1,3)= create "00" coordArraysClass defaultapp
   "caXArray" : freq
   "caYArray" : amp2/1000.
 end create
 dataw(0)= create "00" coordArraysClass defaultapp
   "caXArray" : lag
   "caYArray" : cor1w
 end create
 dataw(1)= create "00" coordArraysClass defaultapp
   "caXArray" : freq
   "caYArray" : coh1w
 end create  
 dataw(2)= create "00" coordArraysClass defaultapp
   "caXArray" : freq
   "caYArray" : phs1w
 end create  
 dataw(3)= create "00" coordArraysClass defaultapp
   "caXArray" : freq
   "caYArray" : amp1w/1000.
 end create  


 xf = (/0.1,0.55/)
 mamp = (/10,100/)

 plot0 = create "plot0" xyPlotClass wks
  "xyCoordData" : dataw(0)
 
  "trYMaxF" : 1
  "trYMinF" : -1
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
  "tiMainString" : "";+(i+1)
  "tiXAxisFont" : 21
  "tiXAxisFontHeightF" : 0.017
  "tiXAxisString" : "lag"
  "tiYAxisFont" : 21
  "tiYAxisFontHeightF" : 0.017
  "tiYAxisString" : "Correlation"
 
  "vpYF" : 0.97
  "vpXF" : xf(0)
  "vpHeightF" : 0.18
  "vpWidthF" : 0.35
 end create
 getvalues plot0
  "xyCoordDataSpec" : spec
 end getvalues
 setvalues spec
  "xyLineColor" : "grey"
 end setvalues
 draw(plot0)

 setvalues plot0
  "xyCoordData" : dataw(1)
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
 getvalues plot0
  "xyCoordDataSpec" : spec
 end getvalues
 setvalues spec
  "xyLineColor" : "grey"
 end setvalues
 draw(plot0)

 setvalues plot0
  "xyCoordData" : dataw(2)
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
 getvalues plot0
  "xyCoordDataSpec" : spec
 end getvalues
 setvalues spec
  "xyLineColor" : "grey"
 end setvalues
 draw(plot0)

 setvalues plot0
  "xyCoordData" : dataw(3)
  "trYMaxF" : mamp(0)
  "trYMinF" : 0
  "trXMaxF" : 0.05
  "trXMinF" : 0
  "vpYF" : 0.25
  "tmYLTickSpacingF" : mamp(0)/4.
  "tmYLMinorPerMajor" : 4
  "tiXAxisString" : "frequency"
  "tiYAxisString" : "Amplitude (10:S:3:N:)"
 end setvalues
 getvalues plot0
  "xyCoordDataSpec" : spec
 end getvalues
 setvalues spec
  "xyLineColor" : "grey"
 end setvalues
 draw(plot0)

 delete(plot0)


 do i=0, 1
 plot0 = create "plot0" xyPlotClass wks
  "xyCoordData" : data0(i,0)

  "trYMaxF" : 1
  "trYMinF" : -1
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
  "vpXF" : xf(i)
  "vpHeightF" : 0.18
  "vpWidthF" : 0.35
 end create
 draw(plot0)

 setvalues plot0
  "xyCoordData" : data0(i,1)
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
  "xyCoordData" : data0(i,2)
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
  "xyCoordData" : data0(i,3)
  "trYMaxF" : mamp(i)
  "trYMinF" : 0 
  "trXMaxF" : 0.05
  "trXMinF" : 0 
  "vpYF" : 0.25
  "tmYLTickSpacingF" : mamp(i)/4.
  "tmYLMinorPerMajor" : 4
  "tiXAxisString" : "frequency"
  "tiYAxisString" : "Amplitude (10:S:3:N:)"
 end setvalues 
 draw(plot0) 


 delete(plot0)
 end do

;-----------------------------------------------

frame(wks)

end



