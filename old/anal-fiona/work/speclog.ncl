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

 imm1 = new(nm,float)
 ave1= new(nm,float)
 ikk2 = new(nk,float)
 ave2= new(nk,float)
 ioo3 = new(no,float)
 ave3= new(no,float)
 
 imm1(:)  = 1./rd1(:,0)
 ikk2(:)  = 1./rd2(:,0)
 ioo3(:)  = 1./rd3(:,0)
 ave1(:)= rd1(:,1)
 ave2(:)= rd2(:,1)
 ave3(:)= rd3(:,1)



wks = create "1sttry" psWorkstationClass defaultapp
	"wkPSFileName" :   "Speclog.ps"
	"wkColorMap"   : "cosam" 
end create

;------------------------------------------

data1 = create "m2_8" coordArraysClass defaultapp
	"caXArray"	: imm1/1000
	"caYArray"	: ave1
end create

data2 = create "k2_8" coordArraysClass defaultapp
      "caXArray"  : ikk2/1000
      "caYArray"  : ave2
end create 

data3 = create "o2_8" coordArraysClass defaultapp
      "caXArray"  : ioo3/60
      "caYArray"  : ave3
end create 

;-------------------------------------------

plot1 =  create "plot0" xyPlotClass wks

  "xyCoordData" : data1

  "trYLog"  : True
  "trXLog"  : True
  "trXReverse" : True
  "trYMaxF" : 10.^4
  "trYMinF" : 10.^(-5)
  "trXMaxF" : 300
  "trXMinF" : 0.1

  "tmXBMode" : "MANUAL"
  "tmXBTickSpacingF" : 1
  "tmXBMinorPerMajor" : 8
  "tmXBLabelFont" : 21
  "tmXBLabelFontHeightF" : 0.012
  "tmXBMajorLengthF" : -0.007
  "tmXBMinorLengthF" : -0.004
  "tmXBAutoPrecision" : False
;  "tmXBFormat" : "@*+^sg"


  "tmYLMinorLengthF" : -0.004
  "tmYLMajorLengthF" : -0.007
  "tmYLMode" : "MANUAL"
  "tmYLTickSpacingF" : 2
  "tmYLMinorPerMajor" : 8
  "tmYLLabelFont" : 21
  "tmYLLabelFontHeightF" : 0.012
;  "tmYLFormat" : "@*+^sg" 
 
  "tiMainFont" : 25
  "tiMainFontHeightF" : 0.021
  "tiXAxisFont" : 12
  "tiXAxisFontHeightF" : 0.014
  "tiXAxisString" : "vertical wavelength (km)"
  "tiYAxisFont" : 12
  "tiYAxisFontHeightF" : 0.014
  "tiYAxisString" : ""

  "vpYF" : 0.90
  "vpXF" : 0.1
  "vpHeightF" : 0.20
  "vpWidthF" : 0.50

end create

draw(plot1)


;-------------------------------------------

plot2 =  create "plot1" xyPlotClass wks

  "xyCoordData" : data2

  "trYLog"  : True
  "trXLog"  : True
  "trXReverse" : True
  "trYMaxF" : 10.^4    ;4
  "trYMinF" : 10.^(-3)
  "trXMaxF" : 1000
  "trXMinF" : 1.

  "tmXBMode" : "MANUAL"
  "tmXBTickSpacingF" : 1
  "tmXBMinorPerMajor" : 8
  "tmXBLabelFont" : 21
  "tmXBLabelFontHeightF" : 0.012
  "tmXBMajorLengthF" : -0.007
  "tmXBMinorLengthF" : -0.004
  "tmXBAutoPrecision" : False
;  "tmXBFormat" : "@*+^sg"

  "tmYLMinorLengthF" : -0.004
  "tmYLMajorLengthF" : -0.007
  "tmYLMode" : "MANUAL"
  "tmYLTickSpacingF" : 2
  "tmYLMinorPerMajor" : 8
  "tmYLLabelFont" : 21
  "tmYLLabelFontHeightF" : 0.012
  "tmYLFormat" : "@*+^sg"

  "tiMainFont" : 25
  "tiMainFontHeightF" : 0.021
  "tiXAxisFont" : 12
  "tiXAxisFontHeightF" : 0.014
  "tiXAxisString" : "horizontal wavelength (km)"
  "tiYAxisFont" : 12
  "tiYAxisFontHeightF" : 0.014
  "tiYAxisString" : "PSD (m:S:2:N:s:S:-2:N:/cycle m:S:-1:N:)"

  "vpYF" : 0.60
  "vpXF" : 0.1
  "vpHeightF" : 0.20
  "vpWidthF" : 0.50

end create
draw(plot2)


;-------------------------------------------

plot3 =  create "plot2" xyPlotClass wks

  "xyCoordData" : data3

  "trYLog"  : True
  "trXLog"  : True
  "trXReverse" : True
  "trYMaxF" : 10.^(2 )
  "trYMinF" : 10.^(-1)
  "trXMaxF" : 1000
  "trXMinF" : 1.

  "tmXBMode" : "MANUAL"
  "tmXBTickSpacingF" : 1
  "tmXBMinorPerMajor" : 8
  "tmXBLabelFont" : 21
  "tmXBLabelFontHeightF" : 0.012
  "tmXBMajorLengthF" : -0.007
  "tmXBMinorLengthF" : -0.004
  "tmXBAutoPrecision" : False
;  "tmXBFormat" : "0@*+^sg"

  "tmYLMinorLengthF" : -0.004
  "tmYLMajorLengthF" : -0.007
  "tmYLMode" : "MANUAL"
  "tmYLTickSpacingF" : 1
  "tmYLMinorPerMajor" : 8
  "tmYLLabelFont" : 21
  "tmYLLabelFontHeightF" : 0.012
  "tmYLFormat" : "@*+^sg"

  "tiMainFont" : 25
  "tiMainFontHeightF" : 0.021
  "tiXAxisFont" : 12
  "tiXAxisFontHeightF" : 0.014
  "tiXAxisString" : "period (min)"
  "tiYAxisFont" : 12
  "tiYAxisFontHeightF" : 0.014
  "tiYAxisString" : ""

  "vpYF" : 0.30
  "vpXF" : 0.1
  "vpHeightF" : 0.20
  "vpWidthF" : 0.50

end create
draw(plot3)


;------------------------------------------


frame(wks)

end




