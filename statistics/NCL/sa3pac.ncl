begin

 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

 setvalues NhlGetWorkspaceObjectId()
   "wsMaximumSize" : 55000000
 end setvalues

;--------------------------------------------

wks = create "1sttry" psWorkstationClass defaultapp
      "wkPSFileName" : "sa3pac.ps"
      "wkColorMap"   : "hotres"
end create

;---------------------------------------

 file1 = "res/sa3pacps.nc"
 f1    = addfile(file1,"r")
 k1    = f1->k1
 o1    = f1->o1
 p01   = f1->PS_pic

 p1 = new((/dimsizes(k1),dimsizes(o1)/),float)
 do i=0, dimsizes(k1)-1
   p1(i,:) = p01(:,i)
 end do

;---------------------------------- cn level

 cnlev = new(7,float)
 cnfil = new(8,integer)

 do i=0,6
  cnlev(i) = 10.^(0.5 * i - 1)
 end do
 cnfil = (/1,4,9,16,36,64,81,100/)*2

;------------------------------------------
 toplength1 = (/150,75,50,37.5,30,20/)
; minorl1    =
 topno1     = 1./toplength1
; minorn1    = 1./minorl1

 toplength2 = (/-0.5,-1,-2,-5,-10,-20,-50,50,20,10,5,2,1,0.5/)
 minorl2    = (/-3,-4,-5,-6,-7,-8,-9,-30,-40,40,30,9,8,7,6,5,4,3/)
 topno2     = 1./12./toplength2
 minorn2    = 1./12./minorl2
 toplength2 = fabs(toplength2)
;------------------------------------------

 field = new(2,graphic)
 field(0) = create "w" scalarFieldClass defaultapp
   "sfDataArray" : p1
   "sfYArray" : k1
   "sfXArray" : o1
 end create

 yf = (/0.95,0.50/)
 str = (/"PSD of SST_Pacific ocean",""/)

 xmin = (/-0.3,-0.045/)
 xmax = (/0.3,0.045/)
 ymin = (/0,0/)
 ymax = (/0.10,0.10/)
 xton = (/"False","True"/)
 xtick = (/0.1,0.015/)
  do i=0, 1
  cont1 = create "contour" contourPlotClass wks
   "cnScalarFieldData" : field(0)
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
   "tmXBTickSpacingF" : xtick(i)
   "tmXBMinorPerMajor" : 4
   "tmXBMajorLengthF" : -0.007
   "tmXBMinorLengthF" : -0.003
   "tmXBLabelFont" : 21
   "tmXBLabelFontHeightF" : 0.012
   "tmXBFormat" : "@*+^sg"

   "tmXUseBottom" : False
   "tmXTLabelsOn" : "True"
   "tmXTOn" : xton(i)
   "tmXTMode" : "EXPLICIT"
   "tmXTValues" : topno2
   "tmXTLabels" : toplength2
   "tmXTMinorValues" : minorn2
   "tmXTLabelFont" : 21
   "tmXTLabelFontHeightF" : 0.013
   "tmXTMajorLengthF" : -0.008
   "tmXTMinorLengthF" : -0.004

   "tmYLMode" : "Manual"
   "tmYLTickSpacingF" : 0.02
   "tmYLMinorPerMajor" : 3
   "tmYLMinorLengthF" : -0.003
   "tmYLMajorLengthF" : -0.007
   "tmYLLabelFont" : 21
   "tmYLLabelFontHeightF" : 0.012
;   "tmYLFormat" : "@*+^sg"

   "tmYUseLeft"   : False
   "tmYRLabelsOn" : True
   "tmYROn" : True
   "tmYRMode" : "EXPLICIT"
   "tmYRValues" : topno1
   "tmYRLabels" : toplength1
;   "tmYRMinorValues" : minorn1
   "tmYRLabelFont" : 21
   "tmYRLabelFontHeightF" : 0.013
   "tmYRMajorLengthF" : -0.008
   "tmYRMinorLengthF" : -0.004

   "vpXF" : 0.15
   "vpYF" : yf(i)
   "vpWidthF" : 0.7
   "vpHeightF" : 0.35

   "trXMinF" : xmin(i)
   "trXMaxF" : xmax(i)
   "trYMinF" : ymin(i)
   "trYMaxF" : ymax(i)

   "tiMainString" : str(i)
   "tiMainFont"  : 22
   "tiMainFontHeightF" : 0.017

   "tiXAxisString" : ":F33: w:F21: (cycle/mon)"
   "tiXAxisFont"  : 21
   "tiXAxisFontHeightF" : 0.015

   "tiYAxisString" : "k (cycle/deg)"
   "tiYAxisFont"  : 21
   "tiYAxisFontHeightF" : 0.015
  end create
  draw(cont1)
  delete(cont1)
  end do


   labels = log10(cnlev)
   label = create "lab" labelBarClass wks
     "lbOrientation" : "horizontal"
     "lbBoxCount" : 8
     "lbFillColors" : cnfil
     "lbMonoFillPattern" : True
     "lbLabelsOn" : True
     "lbLabelStrings" : labels
     "lbAutoManage" : False
     "lbLabelFont" : 21
     "lbLabelFontHeightF" : 0.012
     "lbLabelPosition" : "Bottom";,Right"
     "lbLabelAlignment" : "InteriorEdges"
     "lbLabelAngleF" : 0.
     "lbPerimOn" : False
     "lbTitleString" : "Log:B:10:N:(PSD)"
     "lbTitlePosition" : "Right"
     "lbTitleDirection" : "Across"
     "lbTitleFont" : 21
     "lbTitleFontHeightF" : 0.015
     "vpXF" : 0.50
     "vpYF" : 0.08
     "vpWidthF" : 0.35
     "vpHeightF" : 0.06
   end create
   draw(label)

;-------------------------------------------------

frame(wks)

end



