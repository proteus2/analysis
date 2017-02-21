begin

 wks = create "contour" psWorkstationClass defaultapp
   "wkPSFileName" : "qbo_eof2.eps"
   "wkPSFormat"   : "eps"
   "wkColorMap"   : "BlGrYeOrReVi200"
 end create
 ind_grey = NhlNewColor(wks,0.7,0.7,0.7)

;-----------------------------------------------------------

 fnam = "eof_qbo.nc"
 f1 = addfile(fnam,"r")
 pc  = f1->PC(0:1,:)
 p   = f1->x
 t   = f1->t + 1979.

 ; change the sign of PC2 (if, PC2 is leading)
 pc (1,:) = -pc (1,:)

 nt4col = 12
 ncol = (dimsizes(t)-2)/nt4col+1

 series_0 = new(2,graphic)
 series_1 = new(ncol,graphic)

 series_0(0) = create "T" coordArraysClass defaultapp
   "caXArray"  : pc(0,:)*2.
   "caYArray"  : pc(1,:)*0.
 end create
 series_0(1) = create "T" coordArraysClass defaultapp
   "caXArray"  : pc(0,:)*0.
   "caYArray"  : pc(1,:)*2.
 end create

 do icol=0, ncol-2
   series_1(icol) = create "T" coordArraysClass defaultapp
     "caXArray"  : pc(0,icol*nt4col:(icol+1)*nt4col)
     "caYArray"  : pc(1,icol*nt4col:(icol+1)*nt4col)
   end create
 end do
 series_1(ncol-1) = create "T" coordArraysClass defaultapp
   "caXArray"  : pc(0,(icol-1)*nt4col:)
   "caYArray"  : pc(1,(icol-1)*nt4col:)
 end create

;----------------------------------------------------
; plot
;----------------------------------------------------

 xstr = (/"PC1"/)
 ystr = (/"PC2"/)
 mstr = (/"Trajectory of PCs"/)

 ; define colors used for trajectory
 icol_itv = 170.0/ncol   ; 180 ~ 10
 tr_col = new(ncol,integer)
 do icol=0, ncol-1
   tr_col(icol) = 180 - floattoint(icol*icol_itv)
 end do

 plot_1 = create "lplot" xyPlotClass wks
   "xyCoordData" : series_1

   "tmXBMode" : "Automatic"
   "tmXBMajorLengthF" : -0.006
   "tmXBMinorLengthF" : -0.003
   "tmXBLabelFont"        : 21
   "tmXBLabelFontHeightF" : 0.013

   "tmYLMode" : "Automatic"
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

   "vpXF" : 0.47
   "vpYF" : 0.9
   "vpWidthF"  : 0.30
   "vpHeightF" : 0.30
 end create
 getvalues plot_1
   "xyCoordDataSpec" : spec3_1
   "trXMinF" : pc1xmin
   "trXMaxF" : pc1xmax
   "trYMinF" : pc2ymin
   "trYMaxF" : pc2ymax
 end getvalues
 do i=0, ncol-1
   setvalues spec3_1(i)
     "xyLineColors"     : tr_col(i)
     "xyLineThicknessF" : 1.5
   end setvalues
 end do

 plot_0 = create "lplot" xyPlotClass wks
   "xyCoordData" : series_0

   "tmXBOn" : False
   "tmXTOn" : False
   "tmYLOn" : False
   "tmYROn" : False

   "vpXF" : 0.47
   "vpYF" : 0.9
   "vpWidthF"  : 0.30
   "vpHeightF" : 0.30

   "trXMinF" : pc1xmin
   "trXMaxF" : pc1xmax
   "trYMinF" : pc2ymin
   "trYMaxF" : pc2ymax
 end create
 getvalues plot_0
   "xyCoordDataSpec" : spec3_0
 end getvalues
 do i=0, dimsizes(series_0)-1
   setvalues spec3_0(i)
     "xyLineColor"      : ind_grey
     "xyLineThicknessF" : 1.8
   end setvalues
 end do


 draw(plot_0)

 draw(plot_1)


;------ text ------------------------------------------------

; tid = create "Text" textItemClass wks
;   "txString" : sprintf("%6.2f",evc(i))
;   "txPosXF"  : 0.81
;   "txPosYF"  : 0.875
;   "txFont"   : 22
;   "txFontHeightF" : 0.018
; end create
; draw(tid)

;-------- label bar ----------------------------

 lbstr = new(ncol-1,integer)
 do icol=0, ncol-2
   lbstr(icol) = icol + 1980
 end do

 label = create "lab" labelBarClass wks
   "lbAutoManage" : False
   "lbOrientation" : "vertical"
   "lbBoxCount"   : ncol
   "lbFillColors" : tr_col
   "lbMonoFillPattern" : True
;   "lbLabelsOn"     : True
   "lbLabelStride"  : 5
   "lbLabelStrings"     : lbstr
   "lbLabelFont"        : 21
   "lbLabelFontHeightF" : 0.011
   "lbLabelPosition"    : "Right"
   "lbLabelAlignment"   : "InteriorEdges"
   "lbPerimOn" : False
   "lbTitleOn" : False
;   "lbTitleString"      : ""
;   "lbTitlePosition"    : "Right"
;   "lbTitleDirection"   : "Across"
;   "lbTitleFont"        : 21
;   "lbTitleFontHeightF" : 0.013
   "vpXF" : 0.80
   "vpYF" : 0.90
   "vpWidthF"  : 0.04
   "vpHeightF" : 0.30
 end create
 draw(label)

;-------------------------------------------------

 frame(wks)

end



