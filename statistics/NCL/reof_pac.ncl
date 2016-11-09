begin

 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

;--------------------------------------------

wks = create "1sttry" psWorkstationClass defaultapp
      "wkPSFileName" : "reof_pac.ps"
      "wkColorMap"   : "BlWhRe"
end create

;---------------------------------------

 file1 = "reof_pac.nc"
 f1    = addfile(file1,"r")
 eof1  = f1->REOF
 pc1   = f1->PC
 x   = f1->x
 y   = f1->y
 t   = f1->t
 com = f1->com_reof

;---------------------------------- cn level

 cnlev = new(10,float)
 cnfil = new(11,integer)

 cnlev = (/-1,-0.7,-0.5,-0.3,-0.1,0.1,0.3,0.5,0.7,1/)*0.1  ;*10
 cnfil = (/-35,-23,-12,-7,-3,0,3,7,12,23,35/)
 cnfil = cnfil + 52

;------------------------------------------

 field = new(3,graphic)
 do i=0, 2
 field(i) = create "w" scalarFieldClass defaultapp
   "sfDataArray" : eof1(i,:,:)
   "sfYArray" : y
   "sfXArray" : x
   "sfMissingValueV" : 1.e32
 end create
 end do

 xf = (/0.10,0.10,0.55/)
 yf = (/0.65,0.30,0.30/)
 mstr = (/"REOF1","REOF2","REOF3"/)
 ystr = (/"lat (deg)","lat (deg)",""/)
 do i=0, 2
  cont1 = create "contour" contourPlotClass wks
   "cnScalarFieldData" : field(i)
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
   "tmYLTickSpacingF" : 10
   "tmYLMinorPerMajor" : 1
   "tmYLMinorLengthF" : -0.004
   "tmYLMajorLengthF" : -0.007
   "tmYLLabelFont" : 21
   "tmYLLabelFontHeightF" : 0.011
   "tmYLFormat" : "@*+^sg"

   "vpXF" : xf(i)
   "vpYF" : yf(i)
   "vpWidthF" : 0.30
   "vpHeightF" : 0.20

   "tiMainString" : mstr(i)
   "tiMainFont"  : 22
   "tiMainFontHeightF" : 0.017

   "tiXAxisString" : "lon (deg)"
   "tiXAxisFont"  : 21
   "tiXAxisFontHeightF" : 0.015

   "tiYAxisString" : ystr(i)
   "tiYAxisFont"  : 21
   "tiYAxisFontHeightF" : 0.015
  end create
  plotover = create "plot" logLinPlotClass wks
    "vpXF" : xf(i)
    "vpYF" : yf(i)
    "vpWidthF" : 0.30
    "vpHeightF" : 0.20
    
    "trXMinF" : 120
    "trXMaxF" : 260
    "trYMinF" : -30
    "trYMaxF" : 30
  end create
  overlay(plotover,cont1)
  draw(plotover)
  delete(plotover)
  delete(cont1)
  end do


   label = create "lab" labelBarClass wks
     "lbOrientation" : "horizontal"
     "lbBoxCount" : 11
     "lbFillColors" : cnfil
     "lbMonoFillPattern" : True
     "lbLabelsOn" : True
     "lbLabelStrings" : cnlev*100.
     "lbAutoManage" : False
     "lbLabelFont" : 21
     "lbLabelFontHeightF" : 0.015
     "lbLabelPosition" : "Bottom";,Right"
     "lbLabelAlignment" : "InteriorEdges"
     "lbLabelAngleF" : 0.
     "lbPerimOn" : False
     "lbTitleString" : "(10:S:-2:N:)"
     "lbTitlePosition" : "Right";,"Top"
     "lbTitleDirection" : "Across"
     "lbTitleFont" : 21
     "lbTitleFontHeightF" : 0.015
     "vpXF" : 0.48
     "vpYF" : yf(0) - 0.15
     "vpWidthF"  : 0.50
     "vpHeightF" : 0.05
   end create
   draw(label)


   frame(wks)



  specpc1 = specx_anal(pc1(0,:),1,7,0.1)
  specpc2 = specx_anal(pc1(1,:),1,7,0.1)
  specpc3 = specx_anal(pc1(2,:),1,7,0.1)


  pcdata = new(3,graphic)
  do i=0, 2
  pcdata(i) = create "00" coordArraysClass defaultapp
    "caXArray" : t
    "caYArray" : pc1(i,:)
  end create
  end do
  specdata = new(3,graphic)
  specdata(0) = create "00" coordArraysClass defaultapp
    "caXArray" : 1./specpc1@frq / 12.
    "caYArray" : specpc1@spcx * specpc1@frq / 0.015
  end create
  specdata(1) = create "00" coordArraysClass defaultapp
    "caXArray" : 1./specpc2@frq / 12.
    "caYArray" : specpc2@spcx * specpc1@frq / 0.015
  end create
  specdata(2) = create "00" coordArraysClass defaultapp
    "caXArray" : 1./specpc3@frq / 12.
    "caYArray" : specpc3@spcx * specpc1@frq / 0.015
  end create


  d0 = t
  d0(:) = 0.
  data0 = create "00" coordArraysClass defaultapp
    "caXArray" : t
    "caYArray" : d0
  end create

  yfpc = (/0.90,0.65,0.40/)
  mstrpc = (/"PC1","PC2","PC3"/)
  do i=0, 2
  plot0 = create "plot0" xyPlotClass wks
    "xyCoordData" : pcdata(i)

    "trYMaxF" : 15
    "trYMinF" : -15
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
    "tmYLTickSpacingF"  : 5
    "tmYLMinorPerMajor" : 4
    "tmYLLabelFont" : 21
    "tmYLLabelFontHeightF" : 0.011
;    "tmYLFormat" : "@*+^sg"

    "tiMainFont" : 22
    "tiMainFontHeightF" : 0.020
    "tiMainString" : ""
    "tiXAxisFont" : 21
    "tiXAxisFontHeightF" : 0.017
    "tiXAxisString" : "time (yr)"
    "tiYAxisFont" : 21
    "tiYAxisFontHeightF" : 0.017
    "tiYAxisString" : mstrpc(i)

    "vpYF" : yfpc(i)
    "vpXF" : 0.10
    "vpHeightF" : 0.15
    "vpWidthF" : 0.60
  end create
  draw(plot0)
  setvalues plot0
    "xyCoordData" : data0
  end setvalues
  draw(plot0)

  plot1 = create "plot0" xyPlotClass wks
    "xyCoordData" : specdata(i)

    "trYMaxF" : 1000
    "trYMinF" : 0
    "trXMaxF" : 20
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

    "vpYF" : yfpc(i)
    "vpXF" : 0.75
    "vpHeightF" : 0.15
    "vpWidthF" : 0.20
  end create
  draw(plot1)

  delete(plot0)
  delete(plot1)
  end do


  frame(wks)

end



