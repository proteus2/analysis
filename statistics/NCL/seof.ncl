begin

 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

;--------------------------------------------

wks = create "1sttry" psWorkstationClass defaultapp
      "wkPSFileName" : "seof.ps"
      "wkColorMap"   : "BlWhRe"
end create

j = 0
;---------------------------------------

 file1 = "res/seof.nc"
 f1    = addfile(file1,"r")
 eof1  = f1->EOF_m1
 eof2  = f1->EOF_m4
 eof3  = f1->EOF_m7
 eof4  = f1->EOF_m10
 pc1   = f1->PC
 ev1   = f1->e_value_100
 err1  = f1->err_100
 x   = f1->x
 y   = f1->y
 t   = f1->t
 com = f1->com_ev

 eof = new((/4,dimsizes(com),dimsizes(y),dimsizes(x)/),float)
 eof(0,:,:,:) = eof1(:,:,:)
 eof(1,:,:,:) = eof2(:,:,:)
 eof(2,:,:,:) = eof3(:,:,:)
 eof(3,:,:,:) = eof4(:,:,:)

;---------------------------------- cn level

 cnlev = new(8,float)
 cnfil = new(9,integer)

 cnlev = (/-5,-3,-1,-0.2,0.2,1,3,5/)*0.01
 cnfil = (/-35,-16,-8,-3,0,3,8,16,35/)
 cnfil = cnfil + 52

;------------------------------------------

 field = new((/2,4/),graphic)
 do i=0, 3
 field(j,i) = create "w" scalarFieldClass defaultapp
   "sfDataArray" : eof(i,j,:,:)
   "sfYArray" : y
   "sfXArray" : x
   "sfMissingValueV" : 1.e32
 end create
 end do

 xf = (/0.10,0.39,0.10,0.39/)
 yf = (/0.95,0.95,0.62,0.62/)
 mstr = (/"Jan","Apr","Jul","Oct"/)
 xstr = (/"","","",""/)
 ystr = (/"","","",""/)
 do i=0, 3
  cont1 = create "contour" contourPlotClass wks
   "cnScalarFieldData" : field(j,i)
   "cnFillOn"          :  True           ; color fill
   "cnLevelSelectionMode" : "ExplicitLevels"
   "cnLevels" : cnlev
   "cnFillColors" : cnfil
   "cnLinesOn" : False
   "cnInfoLabelOn" : False
   "cnLineLabelPlacementMode" : "Constant"
   "cnLineLabelFont" : 21
   "cnLineLabelFontHeightF" : 0.011
   "cnLineLabelsOn" : False
   "cnMissingValFillPattern" : 11
   "cnMissingValFillScaleF"  : 0.8
   "cnMissingValFillColor"   : "black"

   "tmXBMode" : "Manual"
   "tmXBTickSpacingF" : 30
   "tmXBMinorPerMajor" : 2
   "tmXBMajorLengthF" : -0.007
   "tmXBMinorLengthF" : -0.004
   "tmXBLabelFont" : 21
   "tmXBLabelFontHeightF" : 0.011
   "tmXBFormat" : "@*+^sg"

   "tmYLMode" : "Manual"
   "tmYLTickSpacingF" : 20
   "tmYLMinorPerMajor" : 1
   "tmYLMinorLengthF" : -0.004
   "tmYLMajorLengthF" : -0.007
   "tmYLLabelFont" : 21
   "tmYLLabelFontHeightF" : 0.011
   "tmYLFormat" : "@*+^sg"

   "vpWidthF" : 0.24
   "vpHeightF" : 0.25

   "tiMainString" : mstr(i)
   "tiMainFont"  : 22
   "tiMainFontHeightF" : 0.017

   "tiXAxisString" : xstr(i)
   "tiXAxisFont"  : 21
   "tiXAxisFontHeightF" : 0.015

   "tiYAxisString" : ystr(i)
   "tiYAxisFont"  : 21
   "tiYAxisFontHeightF" : 0.015
  end create
  plotover = create "plot" logLinPlotClass wks
    "vpXF" : xf(i)
    "vpYF" : yf(i)
    "vpWidthF" : 0.24
    "vpHeightF" : 0.25
    
    "trXMinF" : -70
    "trXMaxF" : 20
    "trYMinF" : -40
    "trYMaxF" : 60
  end create
  overlay(plotover,cont1)
  draw(plotover)
  delete(plotover)
  delete(cont1)
  end do


  data = new(4,graphic)
  data(0) = create "00" coordArraysClass defaultapp
    "caXArray" : com
    "caYArray" : ev1
  end create
  data(1) = create "00" coordArraysClass defaultapp
    "caXArray" : com
    "caYArray" : ev1+err1
  end create
  data(2) = create "00" coordArraysClass defaultapp
    "caXArray" : com
    "caYArray" : ev1-err1
  end create

  d0 = t
  d0(:) = 0.
  data0 = create "00" coordArraysClass defaultapp
    "caXArray" : t
    "caYArray" : d0
  end create

  ; e-values
  new_index = NhlNewMarker(wks,"-",1,0.,0.,1.,6.,0.)
  mind = (/16,new_index,new_index/)
  do i=0, 2
  plotv = create "plot0" xyPlotClass wks
    "xyCoordData" : data(i)

    "trYMaxF" : 20
    "trYMinF" : 0
    "trXMaxF" : 6
    "trXMinF" : 0

    "tmXBMode" : "Manual"
    "tmXBTickSpacingF" : 1
    "tmXBMinorPerMajor" : 0
    "tmXBLabelFont" : 21
    "tmXBLabelFontHeightF" : 0.011
    "tmXBMajorLengthF" :  -0.008
    "tmXBMinorLengthF" :  -0.004
    "tmXBFormat" : "@*+^sg"

    "tmYLMinorLengthF" : -0.004
    "tmYLMajorLengthF" : -0.008
    "tmYLMode" : "MANUAL"
    "tmYLTickSpacingF"  : 5
    "tmYLMinorPerMajor" : 4
    "tmYLLabelFont" : 21
    "tmYLLabelFontHeightF" : 0.011
;    "tmYLFormat" : "@*+^sg"

    "tiMainFont" : 22
    "tiMainFontHeightF" : 0.017
    "tiMainString" : "eigen-values (%)"
    "tiXAxisFont" : 21
    "tiXAxisFontHeightF" : 0.017
    "tiXAxisString" : "rank"
    "tiYAxisFont" : 21
    "tiYAxisFontHeightF" : 0.017
    "tiYAxisString" : ""

    "vpYF" : 0.90
    "vpXF" : 0.70
    "vpHeightF" : 0.20
    "vpWidthF" : 0.18
  end create
  getvalues plotv
    "xyCoordDataSpec" : spec
  end getvalues
  setvalues spec
    "xyMarkLineMode" : "Markers"
    "xyMarkerSizeF"  : 0.005
    "xyMarker"       : mind(i)
  end setvalues
  if (j .eq. 0) then
    draw(plotv)
  end if
  delete(plotv)
  end do


   label = create "lab" labelBarClass wks
     "lbOrientation" : "horizontal"
     "lbBoxCount" : 9
     "lbFillColors" : cnfil
     "lbMonoFillPattern" : True
     "lbLabelsOn" : True
     "lbLabelStrings" : cnlev*100.
     "lbAutoManage" : False
     "lbLabelFont" : 21
     "lbLabelFontHeightF" : 0.013
     "lbLabelPosition" : "Bottom";,Right"
     "lbLabelAlignment" : "InteriorEdges"
     "lbLabelAngleF" : 0.
     "lbPerimOn" : False
     "lbTitleString" : "(10:S:-2:N:)"
     "lbTitlePosition" : "Right";,"Top"
     "lbTitleDirection" : "Across"
     "lbTitleFont" : 21
     "lbTitleFontHeightF" : 0.014
     "vpXF" : 0.55
     "vpYF" : 1.01
     "vpWidthF"  : 0.40
     "vpHeightF" : 0.05
   end create
   draw(label)


;   frame(wks)



  specpc1 = specx_anal(pc1(j,:),1,7,0.1)

  pcdata = create "00" coordArraysClass defaultapp
    "caXArray" : t
    "caYArray" : pc1(j,:)
  end create
  specdata = create "00" coordArraysClass defaultapp
    "caXArray" : 1./specpc1@frq
    "caYArray" : specpc1@spcx * specpc1@frq / 0.05
  end create

  xstrpc = (/"time (yr)"/)
  ystrpc = (/"PC1","PC2"/)
  yfpc = 0.25
  plot0 = create "plot0" xyPlotClass wks
    "xyCoordData" : pcdata

    "trYMaxF" : 20
    "trYMinF" : -20
    "trXMaxF" : 100
    "trXMinF" : 0

    "tmXBMode" : "Manual"
    "tmXBTickSpacingF" : 20 
    "tmXBMinorPerMajor" : 3
    "tmXBLabelFont" : 21
    "tmXBLabelFontHeightF" : 0.011
    "tmXBMajorLengthF" :  -0.008
    "tmXBMinorLengthF" :  -0.004
    "tmXBFormat" : "@*+^sg"

    "tmYLMinorLengthF" : -0.004
    "tmYLMajorLengthF" : -0.008
    "tmYLMode" : "MANUAL"
    "tmYLTickSpacingF"  : 10
    "tmYLMinorPerMajor" : 4
    "tmYLLabelFont" : 21
    "tmYLLabelFontHeightF" : 0.011
;    "tmYLFormat" : "@*+^sg"

    "tiMainFont" : 22
    "tiMainFontHeightF" : 0.020
    "tiMainString" : ""
    "tiXAxisFont" : 21
    "tiXAxisFontHeightF" : 0.017
    "tiXAxisString" : xstrpc
    "tiYAxisFont" : 21
    "tiYAxisFontHeightF" : 0.017
    "tiYAxisString" : ystrpc(j)

    "vpYF" : yfpc
    "vpXF" : 0.15
    "vpHeightF" : 0.12
    "vpWidthF" : 0.55
  end create
  draw(plot0)
  setvalues plot0
    "xyCoordData" : data0
  end setvalues
  draw(plot0)

  plot1 = create "plot0" xyPlotClass wks
    "xyCoordData" : specdata

    "trYMaxF" : 1000
    "trYMinF" : 0
    "trXMaxF" : 30
    "trXMinF" : 1
    "trXLog"  : True

    "tmXBTickSpacingF" : 1
    "tmXBMinorPerMajor" : 8
    "tmXBMode" : "MANUAL"
    "tmXBLabelFont" : 21
    "tmXBLabelFontHeightF" : 0.011
    "tmXBMajorLengthF" :  -0.008
    "tmXBMinorLengthF" :  -0.004
    "tmXBFormat" : "@*+^sg"

    "tmYLOn" : False
    "tmYROn" : False
    "tmYLTickSpacingF" : 200
    "tmYLMinorPerMajor" : 3
    "tmYLMode" : "MANUAL"
    "tmYLLabelFont" : 21
    "tmYLLabelFontHeightF" : 0.011
    "tmYLMinorLengthF" : -0.004
    "tmYLMajorLengthF" : -0.008

    "vpYF" : yfpc
    "vpXF" : 0.75
    "vpHeightF" : 0.12
    "vpWidthF" : 0.15
  end create
  draw(plot1)

  delete(plot0)
  delete(plot1)

  frame(wks)

end



