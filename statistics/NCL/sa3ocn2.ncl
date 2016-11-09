begin

 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

 setvalues NhlGetWorkspaceObjectId()
   "wsMaximumSize" : 55000000
 end setvalues

;--------------------------------------------

wks = create "1sttry" psWorkstationClass defaultapp
      "wkPSFileName" : "sa3ocn2.ps"
      "wkColorMap"   : "hotres"
end create

;---------------------------------------

 file1 = "res/sa3ocn2ps.nc"
 f1    = addfile(file1,"r")
 l1    = f1->l1
 o1    = f1->o1
 p01   = f1->Pac_pic
 p02   = f1->Atl_pic

 p1 = new((/dimsizes(l1),dimsizes(o1)/),float)
 p2 = new((/dimsizes(l1),dimsizes(o1)/),float)
 do i=0, dimsizes(l1)-1
   p1(i,:) = p01(:,i)
   p2(i,:) = p02(:,i)
 end do

;---------------------------------- cn level

 cnlev = new(7,float)
 cnfil = new(8,integer)

 do i=0,6
  cnlev(i) = 10.^(0.5 * i - 1)
 end do
 cnfil = (/1,4,9,16,36,64,81,100/)*2

;------------------------------------------
 toplength1 = (/80,40,(80/3.),20,16/)
; minorl1    = (/(80/3.),
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
   "sfYArray" : l1
   "sfXArray" : o1
 end create
 field(1) = create "w" scalarFieldClass defaultapp
   "sfDataArray" : p2
   "sfYArray" : l1
   "sfXArray" : o1
 end create

 yf = (/0.92,0.47/)
 str = (/"PSD of SST_Pacific ocean","PSD of SST_Atlantic ocean"/)

 xmin = (/-0.05,-0.05/)
 xmax = (/0.05,0.05/)
 ymin = (/0,0/)
 ymax = (/0.15,0.15/)
 xton = (/"True","True"/)
 xtick = (/0.01,0.01/)
  do i=0, 1
  cont1 = create "contour" contourPlotClass wks
   "cnScalarFieldData" : field(i)
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
   "tmYLTickSpacingF" : 0.05
   "tmYLMinorPerMajor" : 4
   "tmYLMinorLengthF" : -0.003
   "tmYLMajorLengthF" : -0.007
   "tmYLLabelFont" : 21
   "tmYLLabelFontHeightF" : 0.012
   "tmYLFormat" : "@*+^sg"

   "tmYUseLeft"   : False
   "tmYRLabelsOn" : True
   "tmYROn" : True
   "tmYRMode" : "EXPLICIT"
   "tmYRValues" : topno1
   "tmYRLabels" : toplength1
   "tmYRLabelFont" : 21
   "tmYRLabelFontHeightF" : 0.013
   "tmYRMajorLengthF" : -0.008
   "tmYRMinorLengthF" : -0.004

   "vpXF" : 0.15
   "vpYF" : yf(i)
   "vpWidthF" : 0.7
   "vpHeightF" : 0.30

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

   "tiYAxisString" : "l (cycle/deg)"
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



