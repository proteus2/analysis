; copy .sound to snd.out , and change the values of the boundary.
; for more accuracy, copy data of .output directly.(especially, pressure data)

begin

 nz = 222

 data = asciiread ("snd.out",(/nz,8/),"float")

 z  = new(nz,float)
 p  = new(nz,float)
 pt = new(nz,float)
 u  = new(nz,float)
 T  = new(nz,float)
 N2 = new(nz,float)
 N  = new(nz,float)
 do k=0, nz-1
   z (k) = data(k,1)
   p (k) = data(k,2)
   pt(k) = data(k,3)
   u (k) = data(k,6)
 end do

 do k=0,nz-1
   T(k) = pt(k)*(p(k)/100000.)^0.287
 end do
 do k=1,nz-2
   N2(k) = 9.8/pt(k)*(pt(k-1)-pt(k+1))/(z(k-1)-z(k+1))
 end do
 N2(0) = 9.8/pt(0)*(pt(0)-pt(1))/(z(0)-z(1))
 N2(nz-1) = 9.8/pt(nz-1)*(pt(nz-2)-pt(nz-1))/(z(nz-2)-z(nz-1))
; N2(0) = N2(1)
; N2(nz-1) = N2(nz-2)

 N = sqrt(N2)
 do k=0,nz-1
   print (k+" "+N(k))
 end do

;--------------------------------------------

cmap = (/     (/1.000,1.000,1.000/), \
              (/0.000,0.000,0.000/), \
              (/1.000,0.000,0.000/), \
              (/0.000,0.353,0.765/), \
              (/0.000,0.706,0.294/), \
              (/0.294,1.000,0.000/), \
              (/0.588,1.000,0.000/), \
              (/0.765,1.000,0.000/), \
              (/0.882,1.000,0.000/), \
              (/1.000,0.941,0.000/)  /)

wks = create "1sttry" psWorkstationClass defaultapp
	"wkPSFileName" :   "snd.ps"
	"wkColorMap"   : cmap
end create

;------------------------------------------

dataleft = create "u" coordArraysClass defaultapp
	"caXArray"	: u
	"caYArray"	: z/1000.
end create

dataright = create "N2" coordArraysClass defaultapp
      "caXArray"  : N
;     "caXArray"  : N2 *1.e4
;     "caXArray"  : pt
      "caYArray"  : z/1000.
end create

datadownl = create "T" coordArraysClass defaultapp
      "caXArray"  : T
      "caYArray"  : z/1000.
end create

datadownr = create "pt" coordArraysClass defaultapp
      "caXArray"  : pt/10000.
      "caYArray"  : z/1000.
end create

 
;-------------------------------------------

plotleft =  create "lplot" xyPlotClass wks

  "xyCoordData" : dataleft

  "trYMaxF" : 110
  "trYMinF" : 0.
  "trXMaxF" : 70.
  "trXMinF" : -80.

  "tmXBMode" : "MANUAL"
  "tmXBTickSpacingF" : 20.
  "tmXBMinorPerMajor" : 3
  "tmXBLabelFont" : 21
  "tmXBLabelFontHeightF" : 0.015
  "tmXBMajorLengthF" : -0.01
  "tmXBMinorLengthF" : -0.005

  "tmYLMinorLengthF" : -0.005
  "tmYLMajorLengthF" : -0.01
  "tmYLMode" : "MANUAL"
  "tmYLTickSpacingF" : 50.
  "tmYLMinorPerMajor" : 4
  "tmYLLabelFont" : 21
  "tmYLLabelFontHeightF" : 0.015
 
 
  "tiMainFont" : 25
  "tiMainFontHeightF" : 0.021
  "tiMainString" : "mean wind"
  "tiXAxisFont" : 12
  "tiXAxisFontHeightF" : 0.016
  "tiXAxisString" : "U (m/s)"
  "tiYAxisFont" : 12
  "tiYAxisFontHeightF" : 0.016
;  "tiYAxisString" : "z (km)"
 
  "vpYF" : 0.95 
  "vpXF" : 0.1
  "vpHeightF" : 0.3
  "vpWidthF" : 0.3

end create

 
plotright =  create "rplot" xyPlotClass wks

  "xyCoordData" : dataright

  "trYMaxF" : 110
  "trYMinF" : 0.
  "trXMaxF" : 0.03
  "trXMinF" : 0.
; "trYMaxF" : 80.0
; "trYMinF" : 0.
; "trXMaxF" : 10000
; "trXMinF" : 300.

  "tmXBMode" : "MANUAL"
  "tmXBTickSpacingF" : 0.01
  "tmXBMinorPerMajor" : 4
  "tmXBLabelFont" : 21
  "tmXBLabelFontHeightF" : 0.015
  "tmXBMajorLengthF" : -0.01
  "tmXBMinorLengthF" : -0.005

  "tmYLMinorLengthF" : -0.005
  "tmYLMajorLengthF" : -0.01
  "tmYLMode" : "MANUAL"
  "tmYLTickSpacingF" : 50.
  "tmYLMinorPerMajor" : 4
  "tmYLLabelFont" : 21
  "tmYLLabelFontHeightF" : 0.015


  "tiMainFont" : 25
  "tiMainFontHeightF" : 0.021
  "tiMainString" : "Brunt-Vaisala freq."
  "tiXAxisFont" : 12
  "tiXAxisFontHeightF" : 0.016
  "tiXAxisString" : "N (/s)"
  "tiYAxisFont" : 12
  "tiYAxisFontHeightF" : 0.016
;  "tiYAxisString" : "z (km)"

  "vpYF" : 0.95
  "vpXF" : 0.55
  "vpHeightF" : 0.3 
  "vpWidthF" : 0.3 
 
end create 
 

plotdownl =  create "dlplot" xyPlotClass wks

  "xyCoordData" : datadownl

  "trYMaxF" : 110
  "trYMinF" : 0.
  "trXMaxF" : 350.
  "trXMinF" : 150.

  "tmXBMode" : "MANUAL"
  "tmXBTickSpacingF" : 50
  "tmXBMinorPerMajor" : 4
  "tmXBLabelFont" : 21
  "tmXBLabelFontHeightF" : 0.015
  "tmXBMajorLengthF" : -0.01
  "tmXBMinorLengthF" : -0.005

  "tmYLMinorLengthF" : -0.005
  "tmYLMajorLengthF" : -0.01
  "tmYLMode" : "MANUAL"
  "tmYLTickSpacingF" : 50.
  "tmYLMinorPerMajor" : 4
  "tmYLLabelFont" : 21
  "tmYLLabelFontHeightF" : 0.015


  "tiMainFont" : 25
  "tiMainFontHeightF" : 0.021
  "tiMainString" : "temperature"
  "tiXAxisFont" : 12
  "tiXAxisFontHeightF" : 0.016
  "tiXAxisString" : "T (K)"
  "tiYAxisFont" : 12
  "tiYAxisFontHeightF" : 0.016
;  "tiYAxisString" : "z (km)"

  "vpYF" : 0.45
  "vpXF" : 0.1
  "vpHeightF" : 0.3
  "vpWidthF" : 0.3

end create


plotdownr =  create "drplot" xyPlotClass wks

  "xyCoordData" : datadownr

  "trYMaxF" : 110
  "trYMinF" : 0.
  "trXMaxF" : 5
  "trXMinF" : 0.

  "tmXBMode" : "MANUAL"
  "tmXBTickSpacingF" : 1
  "tmXBMinorPerMajor" : 4
  "tmXBLabelFont" : 21
  "tmXBLabelFontHeightF" : 0.015
  "tmXBMajorLengthF" : -0.01
  "tmXBMinorLengthF" : -0.005

  "tmYLMinorLengthF" : -0.005
  "tmYLMajorLengthF" : -0.01
  "tmYLMode" : "MANUAL"
  "tmYLTickSpacingF" : 50.
  "tmYLMinorPerMajor" : 4
  "tmYLLabelFont" : 21
  "tmYLLabelFontHeightF" : 0.015


  "tiMainFont" : 25
  "tiMainFontHeightF" : 0.021
  "tiMainString" : "pot. temp."
  "tiXAxisFont" : 12
  "tiXAxisFontHeightF" : 0.016
  "tiXAxisString" : "theta (10000K)"
  "tiYAxisFont" : 12
  "tiYAxisFontHeightF" : 0.016
;  "tiYAxisString" : "z (km)"

  "vpYF" : 0.45
  "vpXF" : 0.55
  "vpHeightF" : 0.3
  "vpWidthF" : 0.3

end create



;------------------------------------------------

draw(plotleft)
draw(plotright)
draw(plotdownl)
draw(plotdownr)
frame(wks)

end




