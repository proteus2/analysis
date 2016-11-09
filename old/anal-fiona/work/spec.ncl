begin

 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

;---------------------------------------

 nm = 150
 nk = 350
 no = 90

 rd1 = asciiread ("mave",(/nm,2/),"float")
 rd2 = asciiread ("kave",(/nk,2/),"float")
 rd3 = asciiread ("oave",(/no,2/),"float")

 mm1 = new(nm,float)
 ave1= new(nm,float)
 kk2 = new(nk,float)
 ave2= new(nk,float)
 oo3 = new(no,float)
 ave3= new(no,float)
 
 mm1(:)  = rd1(:,0)
 kk2(:)  = rd2(:,0)
 oo3(:)  = rd3(:,0)
 ave1(:)= rd1(:,1)
 ave2(:)= rd2(:,1)
 ave3(:)= rd3(:,1)


 toplength0 = (/100,50,20,10,5/)
 minorl0    = (/90,80,70,60,40,30,9,8,7,6/)
 topno0     = 100./toplength0
 minorn0    = 100./minorl0
 toplength1 = (/1000,100,50,20,10/)
 minorl1    = (/900,800,700,600,500,400,300,200,90,80,70,60,40,30/)
 topno1     = 100./toplength1
 minorn1    = 100./minorl1
 toplength2 = (/100,50,20,10,5/)
 minorl2    = (/400,300,200,90,80,70,60,40,30,9,8,7,6/)
 topno2     = 100./6./toplength2
 minorn2    = 100./6./minorl2

;--------------------------------------------


wks = create "1sttry" psWorkstationClass defaultapp
	"wkPSFileName" :   "Spec.ps"
	"wkColorMap"   : "cosam" 
end create

;------------------------------------------

data1 = create "m2_8" coordArraysClass defaultapp
	"caXArray"	: mm1*100000
	"caYArray"	: ave1/1000000
end create

data2 = create "k2_8" coordArraysClass defaultapp
      "caXArray"  : kk2*100000
      "caYArray"  : ave2/1000000
end create 

data3 = create "o2_8" coordArraysClass defaultapp
      "caXArray"  : oo3*1000
      "caYArray"  : ave3/10000
end create 

;-------------------------------------------

plot1 =  create "plot0" xyPlotClass wks

  "xyCoordData" : data1

  "trYMaxF" : 1.5
  "trYMinF" : 0
  "trXMaxF" : 13
  "trXMinF" : 0.

  "tmXUseBottom" : False

  "tmXTLabelsOn" : True
  "tmXTOn" : True

  "tmXBMode" : "MANUAL"
  "tmXBTickSpacingF" : 5
  "tmXBMinorPerMajor" : 4
  "tmXBLabelFont" : 21
  "tmXBLabelFontHeightF" : 0.015
  "tmXBMajorLengthF" : -0.008
  "tmXBMinorLengthF" : -0.004
  "tmXBAutoPrecision" : False
;  "tmXBFormat" : "@*+^sg"

  "tmXTMode" : "EXPLICIT"
  "tmXTValues" : topno0
  "tmXTLabels" : toplength0
  "tmXTMinorValues" : minorn0
  "tmXTLabelFont" : 21
  "tmXTLabelFontHeightF" : 0.015
  "tmXTMajorLengthF" : -0.008
  "tmXTMinorLengthF" : -0.004

  "tmYLMinorLengthF" : -0.004
  "tmYLMajorLengthF" : -0.008
  "tmYLMode" : "MANUAL"
  "tmYLTickSpacingF" : 0.5
  "tmYLMinorPerMajor" : 4
  "tmYLLabelFont" : 21
  "tmYLLabelFontHeightF" : 0.015
;  "tmYLFormat" : "@*+^sg" 
 
  "tiMainFont" : 25
  "tiMainFontHeightF" : 0.021
  "tiXAxisFont" : 12
  "tiXAxisFontHeightF" : 0.016
  "tiXAxisString" : "vertical wave no. (10:S:-5 :N:cycle m:S:-1:N:)"
  "tiYAxisFont" : 12
  "tiYAxisFontHeightF" : 0.016
  "tiYAxisString" : "PSD (10:S:6:N:m:S:2:N:s:S:-2:N:/cycle m:S:-1:N:)"

  "vpYF" : 0.90
  "vpXF" : 0.1
  "vpHeightF" : 0.20
  "vpWidthF" : 0.65

end create
draw(plot1)

xtoptxt = create "xtop" textItemClass wks
      "txString" : "(km)"
      "txPosXF" :"0.72"
      "txPosYF" :"0.92"
      "txFont" : 12
      "txFontHeightF" : 0.016
end create
draw(xtoptxt)


;-------------------------------------------

plot2 =  create "plot1" xyPlotClass wks

  "xyCoordData" : data2

  "trYMaxF" : 3
  "trYMinF" : 0
  "trXMaxF" : 8
  "trXMinF" : 0.

  "tmXUseBottom" : False

  "tmXTLabelsOn" : True
  "tmXTOn" : True

  "tmXBMode" : "MANUAL"
  "tmXBTickSpacingF" : 2
  "tmXBMinorPerMajor" : 3
  "tmXBLabelFont" : 21
  "tmXBLabelFontHeightF" : 0.015
  "tmXBMajorLengthF" : -0.008
  "tmXBMinorLengthF" : -0.004
  "tmXBAutoPrecision" : False
;  "tmXBFormat" : "@*+^sg"

  "tmXTMode" : "EXPLICIT"
  "tmXTValues" : topno1
  "tmXTLabels" : toplength1
  "tmXTMinorValues" : minorn1
  "tmXTLabelFont" : 21
  "tmXTLabelFontHeightF" : 0.015
  "tmXTMajorLengthF" : -0.008
  "tmXTMinorLengthF" : -0.004

  "tmYLMinorLengthF" : -0.004
  "tmYLMajorLengthF" : -0.008
  "tmYLMode" : "MANUAL"
  "tmYLTickSpacingF" : 1
  "tmYLMinorPerMajor" : 9
  "tmYLLabelFont" : 21
  "tmYLLabelFontHeightF" : 0.015
  "tmYLFormat" : "@*+^sg"

  "tiMainFont" : 25
  "tiMainFontHeightF" : 0.021
  "tiXAxisFont" : 12
  "tiXAxisFontHeightF" : 0.016
  "tiXAxisString" : "horizontal wave no. (10:S:-5 :N:cycle m:S:-1:N:)"
  "tiYAxisFont" : 12
  "tiYAxisFontHeightF" : 0.016
;  "tiYAxisString" : "PSD (10:S:5:N:m:S:2:N:s:S:-2:N:/cycle m:S:-1:N:)"

  "vpYF" : 0.60
  "vpXF" : 0.1
  "vpHeightF" : 0.20
  "vpWidthF" : 0.65

end create
draw(plot2)

xtoptxt = create "xtop" textItemClass wks
      "txString" : "(km)"
      "txPosXF" :"0.72"
      "txPosYF" :"0.62"
      "txFont" : 12
      "txFontHeightF" : 0.016
end create
draw(xtoptxt)

;-------------------------------------------

plot3 =  create "plot2" xyPlotClass wks

  "xyCoordData" : data3

  "trYMaxF" : 5
  "trYMinF" : 0
  "trXMaxF" : 2.5
  "trXMinF" : 0.

  "tmXUseBottom" : False

  "tmXTLabelsOn" : True
  "tmXTOn" : True

  "tmXBMode" : "MANUAL"
  "tmXBTickSpacingF" : 0.5
  "tmXBMinorPerMajor" : 4
  "tmXBLabelFont" : 21
  "tmXBLabelFontHeightF" : 0.015
  "tmXBMajorLengthF" : -0.008
  "tmXBMinorLengthF" : -0.004
  "tmXBAutoPrecision" : False
;  "tmXBFormat" : "0@*+^sg"

  "tmXTMode" : "EXPLICIT"
  "tmXTValues" : topno2
  "tmXTLabels" : toplength2
  "tmXTMinorValues" : minorn2
  "tmXTLabelFont" : 21
  "tmXTLabelFontHeightF" : 0.015
  "tmXTMajorLengthF" : -0.008
  "tmXTMinorLengthF" : -0.004

  "tmYLMinorLengthF" : -0.004
  "tmYLMajorLengthF" : -0.008
  "tmYLMode" : "MANUAL"
  "tmYLTickSpacingF" : 1
  "tmYLMinorPerMajor" : 4
  "tmYLLabelFont" : 21
  "tmYLLabelFontHeightF" : 0.015
  "tmYLFormat" : "@*+^sg"

  "tiMainFont" : 25
  "tiMainFontHeightF" : 0.021
  "tiXAxisFont" : 12
  "tiXAxisFontHeightF" : 0.016
  "tiXAxisString" : "frequency (10:S:-3 :N:cycle s:S:-1:N:)"
  "tiYAxisFont" : 12
  "tiYAxisFontHeightF" : 0.016
  "tiYAxisString" : "PSD (10:S:4:N:m:S:2:N:s:S:-2:N:/cycle s:S:-1:N:)"

  "vpYF" : 0.30
  "vpXF" : 0.1
  "vpHeightF" : 0.20
  "vpWidthF" : 0.65

end create
draw(plot3)

xtoptxt = create "xtop" textItemClass wks
      "txString" : "(min)"
      "txPosXF" :"0.72"
      "txPosYF" :"0.32"
      "txFont" : 12
      "txFontHeightF" : 0.016
end create
draw(xtoptxt)

;------------------------------------------


frame(wks)

end




