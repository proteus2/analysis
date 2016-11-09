begin

 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

;--------------------------------------------

wks = create "1sttry" psWorkstationClass defaultapp
      "wkPSFileName" : "svd_pac_uas.ps"
      "wkColorMap"   : "BlWhRe"
end create

;---------------------------------------

 file1 = "res/svd_pac_uas_field1.nc"
 f1    = addfile(file1,"r")
 eof1  = f1->EOF
 pc1   = f1->PC
 var1  = f1->var100
 x1  = f1->x
 y1  = f1->y
 t   = f1->t
 com = f1->com_eof

 file2 = "res/svd_pac_uas_field2.nc"
 f2    = addfile(file2,"r")
 eof2  = f2->EOF
 pc2   = f2->PC
 var2  = f2->var100
 x2  = f2->x
 y2  = f2->y

 file3 = "res/svd_pac_uas_scf.nc"
 f3    = addfile(file3,"r")
 scf   = f3->SCF100
 corr  = f3->corr

;---------------------------------- cn level

 cnlev = new(10,float)
 cnfil = new(11,integer)

 cnlev = (/-1,-0.7,-0.5,-0.3,-0.1,0.1,0.3,0.5,0.7,1/)*0.1
 cnfil = (/-35,-23,-12,-7,-3,0,3,7,12,23,35/)
 cnfil = cnfil + 52

;------------------------------------------

 field = new((/2,2/),graphic)
 do i=0, 1
 field(i,0) = create "w" scalarFieldClass defaultapp
   "sfDataArray" : eof1(i,:,:)
   "sfYArray" : y1
   "sfXArray" : x1
   "sfMissingValueV" : 1.e32
 end create
 field(i,1) = create "w" scalarFieldClass defaultapp
   "sfDataArray" : eof2(i,:,:)
   "sfYArray" : y2
   "sfXArray" : x2
   "sfMissingValueV" : 1.e32
 end create
 end do

 xf = (/(/0.10,0.10/),(/0.55,0.55/)/)
 yf = (/(/0.65,0.30/),(/0.65,0.30/)/)
 mstr = (/(/"Pac. EOF1 ("+var1(0)+"%)","uas EOF1 ("+var2(0)+"%)"/), \
          (/"Pac. EOF2 ("+var1(1)+"%)","uas EOF2 ("+var2(1)+"%)"/) /)
 do j=0, 1
 do i=0, 1
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
   "tmYLTickSpacingF" : 10
   "tmYLMinorPerMajor" : 1
   "tmYLMinorLengthF" : -0.004
   "tmYLMajorLengthF" : -0.007
   "tmYLLabelFont" : 21
   "tmYLLabelFontHeightF" : 0.011
   "tmYLFormat" : "@*+^sg"

   "vpXF" : xf(j,i)
   "vpYF" : yf(j,i)
   "vpWidthF" : 0.35
   "vpHeightF" : 0.20

   "tiMainString" : mstr(j,i)
   "tiMainFont"  : 22
   "tiMainFontHeightF" : 0.017

   "tiXAxisString" : ""
   "tiXAxisFont"  : 21
   "tiXAxisFontHeightF" : 0.015

   "tiYAxisString" : ""
   "tiYAxisFont"  : 21
   "tiYAxisFontHeightF" : 0.015
  end create
  plotover = create "plot" logLinPlotClass wks
    "vpXF" : xf(j,i)
    "vpYF" : yf(j,i)
    "vpWidthF" : 0.35
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
  end do


  ; SCF
  datascf = create "00" coordArraysClass defaultapp
    "caXArray" : com
    "caYArray" : scf
  end create
  datacor = create "00" coordArraysClass defaultapp
    "caXArray" : com
    "caYArray" : corr * 100.
  end create
  plots = create "plot0" xyPlotClass wks
    "xyCoordData" : datascf

    "trYMaxF" : 60
    "trYMinF" : 0
    "trXMaxF" : 10
    "trXMinF" : 0

    "tmXBMode" : "Manual"
    "tmXBTickSpacingF"  : 1
    "tmXBMinorPerMajor" : 0
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
    "tiXAxisString" : "rank"
    "tiYAxisFont" : 21
    "tiYAxisFontHeightF" : 0.017
    "tiYAxisString" : "SCF, Corr. (%)"

    "vpYF" : 0.98
    "vpXF" : 0.1
    "vpHeightF" : 0.20
    "vpWidthF" : 0.8
  end create
  getvalues plots
    "xyCoordDataSpec" : spec
  end getvalues
  setvalues spec
    "xyMarkLineMode" : "MarkLines"
    "xyMarkerSizeF"  : 0.005
    "xyMarker"       : 16
  end setvalues
  draw(plots)

  setvalues plots
    "xyCoordData" : datacor
  end setvalues
  getvalues plots
    "xyCoordDataSpec" : spec
  end getvalues
  setvalues spec
    "xyMarkLineMode" : "Markers"
    "xyMarkerSizeF"  : 0.015
    "xyMarker"       : 5
  end setvalues
  draw(plots)



   frame(wks)



  specpc1 = specx_anal(pc1(0,:),1,7,0.1)
  specpc2 = specx_anal(pc2(0,:),1,7,0.1)
  specpc3 = specx_anal(pc1(1,:),1,7,0.1)
  specpc4 = specx_anal(pc2(1,:),1,7,0.1)

  pcdata = new((/2,2/),graphic)
  do i=0, 1
  pcdata(i,0) = create "00" coordArraysClass defaultapp
    "caXArray" : t
    "caYArray" : pc1(i,:)
  end create
  pcdata(i,1) = create "00" coordArraysClass defaultapp
    "caXArray" : t 
    "caYArray" : pc2(i,:)
  end create
  end do

  specdata = new((/2,2/),graphic)
  specdata(0,0) = create "00" coordArraysClass defaultapp
    "caXArray" : 1./specpc1@frq / 12.
    "caYArray" : specpc1@spcx * specpc1@frq / 0.015
  end create
  specdata(0,1) = create "00" coordArraysClass defaultapp
    "caXArray" : 1./specpc2@frq / 12.
    "caYArray" : specpc2@spcx * specpc1@frq / 0.015 / 3.
  end create
  specdata(1,0) = create "00" coordArraysClass defaultapp
    "caXArray" : 1./specpc3@frq / 12.
    "caYArray" : specpc3@spcx * specpc1@frq / 0.015
  end create 
  specdata(1,1) = create "00" coordArraysClass defaultapp
    "caXArray" : 1./specpc4@frq / 12.
    "caYArray" : specpc4@spcx * specpc1@frq / 0.015 / 3.
  end create


  d0 = t
  d0(:) = 0.
  data0 = create "00" coordArraysClass defaultapp
    "caXArray" : t
    "caYArray" : d0
  end create

  yfpc = (/(/0.95,0.73/),(/0.43,0.21/)/)
  mstrpc = (/(/"PC1 (Pac.)","PC1 (uas)"/),(/"PC2 (Pac.)","PC2 (uas)"/)/)
  xstrpc = (/(/"","time (yr)"/),(/"","time (yr)"/)/)
  do j=0, 1
  do i=0, 1
  plot0 = create "plot0" xyPlotClass wks
    "xyCoordData" : pcdata(j,i)

    "trYMaxF" : 30
    "trYMinF" : -30
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
    "tiXAxisString" : xstrpc(j,i)
    "tiYAxisFont" : 21
    "tiYAxisFontHeightF" : 0.017
    "tiYAxisString" : mstrpc(j,i)

    "vpYF" : yfpc(j,i)
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
    "xyCoordData" : specdata(j,i)

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

    "vpYF" : yfpc(j,i)
    "vpXF" : 0.75
    "vpHeightF" : 0.15
    "vpWidthF" : 0.20
  end create
  draw(plot1)

  delete(plot0)
  delete(plot1)
  end do
  end do


  frame(wks)

end



