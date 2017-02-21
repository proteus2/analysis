begin

 wks = create "contour" psWorkstationClass defaultapp
   "wkPSFileName" : "qbo.eps"
   "wkPSFormat"   : "eps"
   "wkColorMap"   : "posneg_2"
 end create

;-----------------------------------------------------------

 fnam = "u_eq.nc"
 f1 = addfile(fnam,"r")
 ueq0 = f1->U_eq
 p    = f1->lvl
 t    = f1->t + 1979.

 ueq = new((/dimsizes(p),dimsizes(t)/),float)
 do k=0, dimsizes(p)-1
   ueq(k,:) = ueq0(:,k)
 end do

 field1 = create "w" scalarFieldClass defaultapp
   "sfDataArray" : ueq
   "sfYArray"    : p
   "sfXArray"    : t
 end create

;----------------------------------------------------
; plot
;----------------------------------------------------

 cnlvlpos = nice_mnmxintvl(0,max(fabs(ueq)),7,False)

 ncol = 2*floattoint(cnlvlpos(1)/cnlvlpos(2)+1)

 cncol = new(ncol,integer)
 cncolbase = (/10,14/)
 do icol=0, ncol/2-1
   cncol(ncol/2-1-icol) = cncolbase(0) - icol*1
   cncol(ncol/2+icol)   = cncolbase(1) + icol*1
 end do

 contour = create "contour" contourPlotClass wks
   "cnScalarFieldData" : field1
   "cnLinesOn"    : True
   "cnFillOn"     : True
   "cnFillColors" : cncol
   "cnLevelSelectionMode" : "ManualLevels"
   "cnMinLevelValF"       : -cnlvlpos(1)
   "cnLevelSpacingF"      : cnlvlpos(2)
   "cnInfoLabelOn" : False
   "cnLineLabelsOn" : False

   "tmXBMode"          : "Manual"
   "tmXBTickStartF"    : 1975
   "tmXBTickSpacingF"  : 5
   "tmXBMinorPerMajor" : 4
   "tmXBMajorLengthF" : -0.006
   "tmXBMinorLengthF" : -0.003
   "tmXBLabelFont"        : 21
   "tmXBLabelFontHeightF" : 0.013

   "tmYLMode"   : "Explicit"
   "tmYLValues" : p
   "tmYLLabels" : p
   "tmYLMajorLengthF" : -0.006
   "tmYLMinorLengthF" : -0.003
   "tmYLLabelFont"        : 21
   "tmYLLabelFontHeightF" : 0.013

   "tiMainString"      : "U_EQ"
   "tiMainFont"        : 22
   "tiMainFontHeightF" : 0.017

   "tiXAxisString"      : "year"
   "tiXAxisFont"        : 22
   "tiXAxisFontHeightF" : 0.015

   "tiYAxisString"      : "p (hPa)"
   "tiYAxisFont"        : 22
   "tiYAxisFontHeightF" : 0.015

   "vpWidthF"  : 0.75
   "vpHeightF" : 0.15
 end create
 getvalues contour
   "cnLevels" : levels
 end getvalues
 ncn = dimsizes(levels)

; ndash = new(ncn,integer)
; ndash = 0
; do i=0, ncn-1
;   if (levels(i) .lt. 0.0) then
;     ndash(i) = 5
;   end if
; end do
; setvalues contour
;   "cnMonoLineDashPattern" : False
;   "cnLineDashPatterns"    : ndash
; end setvalues

 overl = create "plot" irregularPlotClass wks
   "vpXF" : 0.1
   "vpYF" : 0.9
   "vpWidthF"  : 0.75
   "vpHeightF" : 0.15

   "trXMinF" : min(t)
   "trXMaxF" : max(t)
   "trYMinF" : min(p)
   "trYMaxF" : max(p)
   "trYReverse"  : True
   "trYAxisType" : "LogAxis"
 end create

 overlay(overl,contour)

 draw(overl)

;------ text ------------------------------------------------

; tid =  create "Text" textItemClass wks
;   "txString" : " time = 24hours"
;   "txPosXF" :"0.5"
;   "txPosYF" :"0.96"
;   "txFont" : 22
;   "txFontHeightF" : 0.018
; end create
; draw(tid)

;-------- label bar ----------------------------

 label = create "lab" labelBarClass wks
   "lbAutoManage" : False
   "lbOrientation" : "horizontal"
   "lbBoxCount"   : ncol
   "lbFillColors" : cncol
   "lbMonoFillPattern" : True
;   "lbLabelsOn"     : True
   "lbLabelStrings"     : levels
   "lbLabelFont"        : 21
   "lbLabelFontHeightF" : 0.011
   "lbLabelPosition"    : "Bottom"
   "lbLabelAlignment"   : "InteriorEdges"
   "lbPerimOn" : False
   "lbTitleString"      : "(m s:S:-1:N:)"
   "lbTitlePosition"    : "Right"
   "lbTitleDirection"   : "Across"
   "lbTitleFont"        : 21
   "lbTitleFontHeightF" : 0.013
   "vpXF" : 0.25
   "vpYF" : 0.68
   "vpWidthF"  : 0.50
   "vpHeightF" : 0.03
 end create
 draw(label)

;-------------------------------------------------

 frame(wks)

end



