begin

 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

;--------------------------------------------

wks = create "1sttry" psWorkstationClass defaultapp
      "wkPSFileName" : "condeof.ps"
      "wkColorMap"   : "BlWhRe"
end create

;---------------------------------------

 file1 = "res/cond_eof.nc"
 f1    = addfile(file1,"r")
 eof1  = f1->EOF
 pc1   = f1->PC
 ev1   = f1->e_value_100
 err1  = f1->err_100
 x   = f1->x
 y   = f1->y
 t   = f1->t
 com = f1->com_ev

 file2 = "res/cond_eof_ts.nc"
 f2    = addfile(file2,"r")
 ts    = f2->TS

;---------------------------------- cn level

 cnlev = new(10,float)
 cnfil = new(11,integer)

 cnlev = (/-1,-0.7,-0.5,-0.3,-0.1,0.1,0.3,0.5,0.7,1/)*0.1
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

 xf = (/0.10,0.40,0.70/)
 yf = (/0.60,0.60,0.60/)
 mstr = (/"EOF1","EOF2","EOF3"/)
 ystr = (/"lat (deg)","",""/)
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
   "vpWidthF" : 0.25
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
    "vpWidthF" : 0.25
    "vpHeightF" : 0.20
    
    "trXMinF" : -70
    "trXMaxF" : 20
    "trYMinF" : -30
    "trYMaxF" : 30
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

; e-values
  new_index = NhlNewMarker(wks,"-",1,0.,0.,1.,6.,0.)
  mind = (/16,new_index,new_index/)
  do i=0, 2
  plotv = create "plot0" xyPlotClass wks
    "xyCoordData" : data(i)

    "trYMaxF" : 30
    "trYMinF" : 0
    "trXMaxF" : 10
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
    "tiMainFontHeightF" : 0.020
    "tiMainString" : ""
    "tiXAxisFont" : 21
    "tiXAxisFontHeightF" : 0.017
    "tiXAxisString" : "rank"
    "tiYAxisFont" : 21
    "tiYAxisFontHeightF" : 0.017
    "tiYAxisString" : "eigen-values (%)"

    "vpYF" : 0.98
    "vpXF" : 0.1
    "vpHeightF" : 0.20
    "vpWidthF" : 0.8
  end create
  getvalues plotv
    "xyCoordDataSpec" : spec
  end getvalues
  setvalues spec
    "xyMarkLineMode" : "Markers"
    "xyMarkerSizeF"  : 0.005
    "xyMarker"       : mind(i)
  end setvalues
  draw(plotv)
  delete(plotv)
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
     "vpXF" : 0.25
     "vpYF" : yf(0) - 0.35
     "vpWidthF"  : 0.50
     "vpHeightF" : 0.05
   end create
   draw(label)


   frame(wks)



  spects  = specx_anal(ts(:),1,7,0.1)
  specpc1 = specx_anal(pc1(0,:),1,7,0.1)
  specpc2 = specx_anal(pc1(1,:),1,7,0.1)
  specpc3 = specx_anal(pc1(2,:),1,7,0.1)


  d0 = t
  d0(:) = 0.
  data0 = create "00" coordArraysClass defaultapp
    "caXArray" : t
    "caYArray" : d0
  end create

  pcdata = new(4,graphic)
  pcdata(0) = create "00" coordArraysClass defaultapp
    "caXArray" : t
    "caYArray" : ts
  end create
  do i=1, 3
  pcdata(i) = create "00" coordArraysClass defaultapp
    "caXArray" : t
    "caYArray" : pc1(i-1,:)
  end create
  end do

  specdata = new(4,graphic)
  specdata(0) = create "00" coordArraysClass defaultapp
    "caXArray" : 1./spects@frq / 12.
    "caYArray" : spects@spcx * spects@frq / 0.01 * 100
  end create
  specdata(1) = create "00" coordArraysClass defaultapp
    "caXArray" : 1./specpc1@frq / 12.
    "caYArray" : specpc1@spcx * specpc1@frq / 0.01 / 1.2
  end create
  specdata(2) = create "00" coordArraysClass defaultapp
    "caXArray" : 1./specpc2@frq / 12.
    "caYArray" : specpc2@spcx * specpc1@frq / 0.01 / 1.2
  end create
  specdata(3) = create "00" coordArraysClass defaultapp
    "caXArray" : 1./specpc3@frq / 12.
    "caYArray" : specpc3@spcx * specpc1@frq / 0.01 / 1.2
  end create


; PC and spec
  yfpc = (/0.95,0.65,0.45,0.25/)
  mstrpc = (/"Nino3.4 index","PC1","PC2","PC3"/)
  xstrpc = (/"time (yr)","","","time (yr)"/)
  ymax = (/2,20,20,20/)
  ymin = -ymax
  ytick = (/1,10,10,10/)
  do i=0, 3
  plot0 = create "plot0" xyPlotClass wks
    "xyCoordData" : pcdata(i)

    "trYMaxF" : ymax(i)
    "trYMinF" : ymin(i)
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
    "tmYLTickSpacingF"  : ytick(i)
    "tmYLMinorPerMajor" : 4
    "tmYLLabelFont" : 21
    "tmYLLabelFontHeightF" : 0.011
    "tmYLFormat" : "@*+^sg"

    "tiMainFont" : 22
    "tiMainFontHeightF" : 0.020
    "tiMainString" : ""
    "tiXAxisFont" : 21
    "tiXAxisFontHeightF" : 0.017
    "tiXAxisString" : xstrpc(i)
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



