begin

 wks = create "contour" psWorkstationClass defaultapp
   "wkPSFileName" : "qbo_eof1.eps"
   "wkPSFormat"   : "eps"
   "wkColorMap"   : "BlGrYeOrReVi200"
 end create
 ind_grey = NhlNewColor(wks,0.7,0.7,0.7)

;-----------------------------------------------------------

 fnam = "eof_qbo.nc"
 f1 = addfile(fnam,"r")
 eof = f1->EOF(0:1,0,:)
 pc  = f1->PC(0:1,:)
 p   = f1->x
 t   = f1->t + 1979.
 evc = f1->e_value_100(0:4)

 ; change the sign of EOF2 and PC2 (if, PC2 is leading)
 eof(1,:) = -eof(1,:)
 pc (1,:) = -pc (1,:)

 series1 = new(3,graphic)
 series2 = new(3,graphic)
 do i=1, 2
   series1(i) = create "T" coordArraysClass defaultapp
     "caXArray"  : eof(i-1,:)
     "caYArray"  : p
   end create
   series2(i) = create "T" coordArraysClass defaultapp
     "caXArray"  : t
     "caYArray"  : pc(i-1,:)
   end create
 end do
 series1(0) = create "T" coordArraysClass defaultapp
   "caXArray"  : p*0.
   "caYArray"  : p
 end create
 series2(0) = create "T" coordArraysClass defaultapp
   "caXArray"  : t
   "caYArray"  : t*0.
 end create

 series3 = create "T" coordArraysClass defaultapp
   "caXArray"  : (/1,2,3,4,5/)
   "caYArray"  : evc
 end create

;----------------------------------------------------
; plot
;----------------------------------------------------

 xstr = (/"EOF","year","mode"/)
 ystr = (/"p (hPa)","PC","(%)"/)
 mstr = (/"EOFs of U_EQ","Principal components","Normalized eigen-values"/)

 idash = (/0,0,11/)
 icolor = (/ind_grey,1,1/)

 plot1 = create "lplot" xyPlotClass wks
   "xyCoordData" : series1

   "tmXBMode" : "Automatic"
   "tmXBMajorLengthF" : -0.006
   "tmXBMinorLengthF" : -0.003
   "tmXBLabelFont"        : 21
   "tmXBLabelFontHeightF" : 0.013
   "tmXBFormat" : "@*+^sg"

   "tmYLMode"   : "Explicit"
   "tmYLValues" : p
   "tmYLLabels" : p
   "tmYLMajorLengthF" : -0.006
   "tmYLMinorLengthF" : -0.003
   "tmYLLabelFont"        : 21
   "tmYLLabelFontHeightF" : 0.013

   "tiMainString"      : mstr(0)
   "tiMainFont"        : 22
   "tiMainFontHeightF" : 0.017

   "tiXAxisString"      : xstr(0)
   "tiXAxisFont"        : 22
   "tiXAxisFontHeightF" : 0.015

   "tiYAxisString"      : ystr(0)
   "tiYAxisFont"        : 22
   "tiYAxisFontHeightF" : 0.015

   "vpXF" : 0.1
   "vpYF" : 0.9
   "vpWidthF"  : 0.25
   "vpHeightF" : 0.30

;   "trXMinF" : min(eof)
;   "trXMaxF" : max(eof)
   "trYMinF" : min(p)
   "trYMaxF" : max(p)
   "trYReverse"  : True
   "trYAxisType" : "LogAxis"
 end create
 getvalues plot1
   "xyCoordDataSpec" : spec1
 end getvalues
 do i=0, dimsizes(series1)-1
   setvalues spec1(i)
     "xyDashPatterns"   : idash(i)
     "xyLineColors"     : icolor(i)
     "xyLineThicknessF" : 1.8
   end setvalues
 end do

 draw(plot1)


 plot2 = create "lplot" xyPlotClass wks
   "xyCoordData" : series2

   "tmXBMode"          : "Manual"
   "tmXBTickStartF"    : 1975
   "tmXBTickSpacingF"  : 5
   "tmXBMinorPerMajor" : 4
   "tmXBMajorLengthF" : -0.006
   "tmXBMinorLengthF" : -0.003
   "tmXBLabelFont"        : 21
   "tmXBLabelFontHeightF" : 0.013

   "tmYLMode" : "Automatic"
   "tmYLMajorLengthF" : -0.006
   "tmYLMinorLengthF" : -0.003
   "tmYLLabelFont"        : 21
   "tmYLLabelFontHeightF" : 0.013

   "tiMainString"      : mstr(1)
   "tiMainFont"        : 22
   "tiMainFontHeightF" : 0.017

   "tiXAxisString"      : xstr(1)
   "tiXAxisFont"        : 22
   "tiXAxisFontHeightF" : 0.015

   "tiYAxisString"      : ystr(1)
   "tiYAxisFont"        : 22
   "tiYAxisFontHeightF" : 0.015

   "vpXF" : 0.1
   "vpYF" : 0.48
   "vpWidthF"  : 0.75
   "vpHeightF" : 0.15

;   "trXMinF" : min(t)
;   "trXMaxF" : max(t)
;   "trYMinF" : min(pc)
;   "trYMaxF" : max(pc)
 end create
 getvalues plot2
   "xyCoordDataSpec" : spec2
 end getvalues
 do i=0, dimsizes(series2)-1
   setvalues spec2(i)
     "xyDashPatterns"   : idash(i)
     "xyLineColors"     : icolor(i)
     "xyLineThicknessF" : 1.8
   end setvalues
 end do

 draw(plot2)


 plot3 = create "lplot" xyPlotClass wks
   "xyCoordData" : series3

   "tmXBMode"    : "Explicit"
   "tmXBValues" : (/1,2,3,4,5/)
   "tmXBLabels" : (/1,2,3,4,5/)
   "tmXBMajorLengthF" : -0.006
   "tmXBMinorLengthF" : -0.003
   "tmXBLabelFont"        : 21
   "tmXBLabelFontHeightF" : 0.013
   "tmXBFormat" : "@*+^sg"

   "tmYLMode" : "Automatic"
   "tmYLMajorLengthF" : -0.006
   "tmYLMinorLengthF" : -0.003
   "tmYLLabelFont"        : 21
   "tmYLLabelFontHeightF" : 0.013

   "tiMainString"      : mstr(2)
   "tiMainFont"        : 22
   "tiMainFontHeightF" : 0.017

   "tiXAxisString"      : xstr(2)
   "tiXAxisFont"        : 22
   "tiXAxisFontHeightF" : 0.015

   "tiYAxisString"      : ystr(2)
   "tiYAxisFont"        : 22
   "tiYAxisFontHeightF" : 0.015

   "vpXF" : 0.54
   "vpYF" : 0.9
   "vpWidthF"  : 0.30
   "vpHeightF" : 0.30

   "trXMinF" : 0.5
   "trXMaxF" : 5.5
 end create
 getvalues plot3
   "xyCoordDataSpec" : spec3
 end getvalues
 setvalues spec3
   "xyMarkLineMode" : "Markers"
   "xyMarker"       : 0
   "xyMarkerSizeF"  : 0.015
 end setvalues

 draw(plot3)


;------ text ------------------------------------------------

; tid = create "Text" textItemClass wks
;   "txString" : sprintf("%6.2f",evc(i))
;   "txPosXF"  : 0.81
;   "txPosYF"  : 0.875
;   "txFont"   : 22
;   "txFontHeightF" : 0.018
; end create
; draw(tid)

;-- legend ---------------------------------------

 lg1 = create "lg1" legendClass wks
   "lgAutoManage"   : False
   "lgItemCount"    : 2
   "lgLabelStrings" : (/"EOF1","EOF2"/)

   "lgMonoLineColor"  : True
   "lgDashIndexes"    : idash(1:2)
   "lgLineThicknessF" : 1.8
   "lgLabelFont" : 21
   "lgLabelFontHeightF" : 0.013
   "lgPerimOn"   : False

   "vpXF"        : 0.25
   "vpYF"        : 0.90
   "vpWidthF"    : 0.08
   "vpHeightF"   : 0.05
 end create
 draw(lg1)

;-------------------------------------------------

 frame(wks)

end

