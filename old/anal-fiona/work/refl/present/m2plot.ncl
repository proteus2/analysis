begin

wks = create "1sttry" psWorkstationClass defaultapp
      "wkPSFileName" :   "m2p.ps"
      "wkColorMap"   : "gsdtol"
end create

;---------------------------------------

f1  = addfile("/usr/users/kyh/convnc/m2k/m20.nc","r")
m20 = f1->m2
z   = f1->z
o   = f1->o

f2  = addfile("/usr/users/kyh/convnc/m2k/m30.nc","r")
m30 = f2->m2

f3  = addfile("/usr/users/kyh/convnc/m2k/m100.nc","r")
m100= f3->m2

maxo1  = new(352,float)
maxo2  = new(352,float)
maxo3  = new(352,float)

do i=0, 351
  maxo1(i)  = o( maxind(m20(i,:)) )
  maxo2(i)  = o( maxind(m30(i,:)) )
  maxo3(i)  = o( maxind(m100(i,:)) )
end do


 toplength1 = (/-1,-5,-10,-20,-60,60,20,10,5,1/)
 minorl1    = (/-2,-3,-4,-6,-7,-8,-9,-30,-40,-50,-70,-80,-90,-100,-200, \
                100,200,90,80,70,50,40,30,9,8,7,6,4,3,2/)
 topno1     = 100./6./toplength1
 minorn1    = 100./6./minorl1
 toplength1 = abs(toplength1)


;----------------------------------------

  field1 = create "w" scalarFieldClass defaultapp
    "sfDataArray" : m20
    "sfYArray" : z
    "sfXArray" : o*1000.
  end create

  field2 = create "w" scalarFieldClass defaultapp
    "sfDataArray" : m30
    "sfYArray" : z
    "sfXArray" : o*1000.
  end create

  field3 = create "w" scalarFieldClass defaultapp
    "sfDataArray" : m100
    "sfYArray" : z
    "sfXArray" : o*1000.
  end create

  data1 = create "T" coordArraysClass defaultapp
      "caXArray"  : maxo1*1000.
      "caYArray"  : z
  end create

  data2 = create "T" coordArraysClass defaultapp
      "caXArray"  : maxo2*1000.
      "caYArray"  : z
  end create

  data3 = create "T" coordArraysClass defaultapp
      "caXArray"  : maxo3*1000.
      "caYArray"  : z
  end create



;----------------------------------------

  cont1 = create "contour" contourPlotClass wks

   "cnScalarFieldData" : field1
   "cnFillOn"        :  True          ; color fill
   "cnLevelSelectionMode" : "ExplicitLevels"
   "cnLevels" : 0.
   "cnFillColors" : (/25,0/)
   "cnLinesOn" : True
   "cnInfoLabelOn" : False
   "cnLineLabelPlacementMode" : "Constant"
   "cnLineLabelFont" : 21
   "cnLineLabelFontHeightF" : 0.011
   "cnLineLabelsOn" : False

   "tmXBMode" : "Manual"
   "tmXBTickSpacingF" : 1
   "tmXBMinorPerMajor" : 4
   "tmXBMajorLengthF" : -0.007
   "tmXBMinorLengthF" : -0.003
   "tmXBLabelFont" : 21
   "tmXBLabelFontHeightF" : 0.011
   "tmXBFormat" : "@*+^sg"

  "tmXUseBottom" : False
  "tmXTLabelsOn" : True
  "tmXTOn" : True
  "tmXTMode" : "EXPLICIT"
  "tmXTValues" : topno1
  "tmXTLabels" : toplength1
  "tmXTMinorValues" : minorn1
  "tmXTLabelFont" : 21
  "tmXTLabelFontHeightF" : 0.011
  "tmXTMajorLengthF" : -0.007
  "tmXTMinorLengthF" : -0.003


   "tmYLMode" : "Manual"
   "tmYLTickSpacingF" : 15
   "tmYLMinorPerMajor" : 2
   "tmYLMinorLengthF" : -0.003
   "tmYLMajorLengthF" : -0.007
   "tmYLLabelFont" : 21
   "tmYLLabelFontHeightF" : 0.011

   "vpXF" : 0.1
   "vpYF" : 0.93
   "vpWidthF" : 0.60
   "vpHeightF" : 0.20

   "trXMinF" : -6
   "trXMaxF" : 6
   "trYMinF" : 0
   "trYMaxF" : 105

   "tiMainString" : "period (min)"
   "tiMainFont"  : 21
   "tiMainFontHeightF" : 0.015

;   "tiXAxisString" : "frequency (10:S:-3:N: cycle s:S:-1:N:)"
   "tiXAxisFont"  : 21
   "tiXAxisFontHeightF" : 0.015
;   "tiXAxisOn" : True

   "tiYAxisString" : "z (km)"
   "tiYAxisFont"  : 21
   "tiYAxisFontHeightF" : 0.015
  end create

  draw(cont1)

  setvalues cont1
   "cnScalarFieldData" : field2
   "vpXF" : 0.1
   "vpYF" : 0.68
   "tiMainString" : ""
  end setvalues
  draw(cont1)

  setvalues cont1
   "cnScalarFieldData" : field3
   "vpXF" : 0.1
   "vpYF" : 0.43
   "tiMainString" : ""
   "tiXAxisString" : "frequency (10:S:-3:N: cycle s:S:-1:N:)"
   "tiXAxisFont"  : 21
   "tiXAxisFontHeightF" : 0.015
  end setvalues 
  draw(cont1)

;-------------------------

plot1 =  create "lplot" xyPlotClass wks

  "xyCoordData" : data1

   "tmXBMode" : "Manual"
   "tmXBTickSpacingF" : 1
   "tmXBMinorPerMajor" : 4
   "tmXBMajorLengthF" : -0.007
   "tmXBMinorLengthF" : -0.003
   "tmXBLabelFont" : 21
   "tmXBLabelFontHeightF" : 0.011
   "tmXBFormat" : "@*+^sg"

  "tmXUseBottom" : False
  "tmXTLabelsOn" : True
  "tmXTOn" : True
  "tmXTMode" : "EXPLICIT"
  "tmXTValues" : topno1
  "tmXTLabels" : toplength1
  "tmXTMinorValues" : minorn1
  "tmXTLabelFont" : 21
  "tmXTLabelFontHeightF" : 0.011
  "tmXTMajorLengthF" : -0.007
  "tmXTMinorLengthF" : -0.003


   "tmYLMode" : "Manual"
   "tmYLTickSpacingF" : 15
   "tmYLMinorPerMajor" : 2
   "tmYLMinorLengthF" : -0.003
   "tmYLMajorLengthF" : -0.007
   "tmYLLabelFont" : 21
   "tmYLLabelFontHeightF" : 0.011

   "vpXF" : 0.1
   "vpYF" : 0.93
   "vpWidthF" : 0.60
   "vpHeightF" : 0.20

   "trXMinF" : -6
   "trXMaxF" : 6
   "trYMinF" : 0
   "trYMaxF" : 105

;   "tiMainString" : "period (min)"
   "tiMainFont"  : 21
   "tiMainFontHeightF" : 0.015

;   "tiXAxisString" : "frequency (10:S:-3:N: cycle s:S:-1:N:)"
   "tiXAxisFont"  : 21
   "tiXAxisFontHeightF" : 0.015
;   "tiXAxisOn" : True

   "tiYAxisString" : "z (km)"
   "tiYAxisFont"  : 21
   "tiYAxisFontHeightF" : 0.015

end create
getvalues plot1
  "xyCoordDataSpec" : spec
end getvalues
setvalues spec
  "xyMonoDashPattern" : True
  "xyDashPattern" : 1
end setvalues
draw(plot1)

setvalues plot1
  "xyCoordData" : data2
  "vpXF" : 0.1
  "vpYF" : 0.68
end setvalues
getvalues plot1
  "xyCoordDataSpec" : spec
end getvalues
setvalues spec
  "xyMonoDashPattern" : True
  "xyDashPattern" : 1
end setvalues
draw(plot1)

setvalues plot1
  "xyCoordData" : data3
  "vpXF" : 0.1
  "vpYF" : 0.43
end setvalues
getvalues plot1
  "xyCoordDataSpec" : spec
end getvalues 
setvalues spec
  "xyMonoDashPattern" : True
  "xyDashPattern" : 1
end setvalues 
draw(plot1) 





frame(wks)


end
