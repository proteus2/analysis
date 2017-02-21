begin

;------------------------

 setvalues NhlGetWorkspaceObjectId()
   "wsMaximumSize" : 55000000
 end setvalues

 load "/usr/local/ncarg/lib/ncarg/nclscripts/csm/gsn_code.ncl"

 wks = create "ncgm" ncgmWorkstationClass defaultapp
     "wkMetaName" : "ncgm/movie60.ncgm"
     "wkColorMap" : "BlWhRe"
     "wkBackgroundColor" : (/1.,1.,1./)
 end create

 wi = 0.84         ; x * 7
 he = 0.36         ; z * 6

 yf = 0.93 - (0.72-he)

;---------------------------------- cn level
 npos = 6
 nneg = 6

 cnlev = new((npos+nneg+1),float)
 cnfil = new((npos+nneg+2),integer)

; cnlev = (/-40,-32,-25,-16,-10,-5,-3,-1,0,1,3,5,10,16,25,32,40/)
 cnlev = (/-40,-28,-18,-10,-4,-1,0,1,4,10,18,28,40/)
; cnfil = (/-49,-35,-30,-23,-17,-12,-7,-4,-2,0,4,7,11,18,23,30,35,49/)
 cnfil = (/-49,-35,-23,-12,-7,-4,-2,0,4,7,12,23,35,49/)
 cnfil = cnfil + 52

;-------------------------------- reading

 do n=1,240

 time = n * 120  ; 60

 time1= time 

 hour1= time1 / 3600
 min1 = (time1 - hour1*3600)/60

 fdir = "/export10/kyh/ing/x12src/nc/"

 fnam1 = "w" +sprinti("%0.5i",time1) + ".nc"
 f1 = addfile(fdir+fnam1,"r")
 w1 = f1->prt
 x = f1->x
 z = f1->z

 fnam2 = "pt"+sprinti("%0.5i",time1) + ".nc"
 f2 = addfile(fdir+fnam2,"r")
 pt1 = f2->prt
 pt1 = log(pt1)

 fnam5 = "qc"+sprinti("%0.5i",time1) + ".nc"
 f5 = addfile(fdir+fnam5,"r")
 qc1 = f5->prt

; fnam7 = "ri"+sprinti("%0.5i",time1) + ".nc"
; f7 = addfile(fdir+fnam7,"r")
; ri1 = f7->prt

; ri1(0:49,550:900) = 10.

;--------------------------------------------- field

  field1 = create "w" scalarFieldClass defaultapp
    "sfDataArray" : w1
    "sfYArray" : z
    "sfXArray" : x
  end create

  field2 = create "pt" scalarFieldClass defaultapp
    "sfDataArray" : pt1
    "sfYArray" : z
    "sfXArray" : x
  end create

  field5 = create "qc" scalarFieldClass defaultapp
    "sfDataArray" : qc1
    "sfYArray" : z
    "sfXArray" : x
  end create

;  field7 = create "ri" scalarFieldClass defaultapp
;    "sfDataArray" : ri1
;    "sfYArray" : z
;    "sfXArray" : x
;  end create


;--------------------------------------------- contour
; wprt1  
;----------------------------------------------------
  
  cont1 = create "contour" contourPlotClass wks
   "cnScalarFieldData" : field1
   "cnFillOn"        :  True           ; color fill
   "cnLevelSelectionMode" : "ExplicitLevels"
   "cnLevels"     : cnlev 
   "cnFillColors" : cnfil
   "cnLinesOn"    : False 
   "cnInfoLabelOn" : False
   "cnLineLabelPlacementMode" : "Constant"
   "cnLineLabelFont" : 21
   "cnLineLabelFontHeightF" : 0.011
   "cnLineLabelsOn" : True

   "tmXBMode" : "Manual"
   "tmXBTickSpacingF" : 100
   "tmXBMinorPerMajor" : 4
   "tmXTMinorPerMajor" : 4
   "tmXBMajorLengthF" : -0.012
   "tmXBMinorLengthF" : -0.006
   "tmXBLabelFont" : 21
   "tmXBLabelFontHeightF" : 0.018

   "tmYLMode" : "Manual"
   "tmYLTickSpacingF" : 10
   "tmYLMinorPerMajor" : 4
   "tmYLMinorLengthF" : -0.008
   "tmYLMajorLengthF" : -0.016
   "tmYLLabelFont" : 21
   "tmYLLabelFontHeightF" : 0.020

   "vpXF" : 0.11
   "vpYF" : yf
   "vpWidthF"  : wi
   "vpHeightF" : he

   "trXMinF" : 0 
   "trXMaxF" : 1200
   "trYMinF" : 0
   "trYMaxF" : 60

   "tiMainString" : hour1+" h  "+min1+" m "
   "tiMainFont"  : 22
   "tiMainFontHeightF" : 0.025

   "tiXAxisString" : "X (km)"
   "tiXAxisFont"  : 21
   "tiXAxisFontHeightF" : 0.022
;   "tiXAxisOn" : True 

   "tiYAxisString" : "Z (km)"
   "tiYAxisFont"  : 21
   "tiYAxisFontHeightF" : 0.022

   "pmLabelBarDisplayMode" : "Conditional"
   "pmLabelBarSide" : "Bottom"
   "lbOrientation"  : "Horizontal"
   "pmLabelBarWidthF"  : 0.65
   "pmLabelBarHeightF" : 0.11
   "pmLabelBarOrthogonalPosF" : 0.05
   "pmLabelBarParallelPosF"   : 0.55
   "lbPerimOn"    : False
   "lbAutoManage" : False
   "lbLabelFontHeightF" : 0.018
   "lbLabelFontAspectF" : 3.0
   "lbTitleString"   : "(m/s)"
   "lbTitlePosition" : "Right"
   "lbTitleFontHeightF" : 0.015
   "lbTitleDirection" : "Across"
  end create
  draw(cont1)

;---------------------------------------------------
; ptprt1
;---------------------------------------------------

  cont2 = create "contour" contourPlotClass wks
   "cnScalarFieldData" : field2
   "cnFillOn"          :  False
   "cnLevelSelectionMode" : "ManualLevels"
   "cnMaxLevelCount"      : 50
   "cnLevelSpacingF"      : 0.1
   "cnLinesOn"            : True
   "cnMonoLineThickness"  : True
   "cnLineThicknessF" : 0.5
   "cnLineLabelsOn" : False
   "cnInfoLabelOn"  : False

   "tmXBMode" : "Manual"
   "tmXBTickSpacingF" : 100
   "tmXBMinorPerMajor" : 4
   "tmXTMinorPerMajor" : 4
   "tmXBMajorLengthF" : -0.012
   "tmXBMinorLengthF" : -0.006
   "tmXBLabelFont" : 21
   "tmXBLabelFontHeightF" : 0.018

   "tmYLMode" : "Manual"
   "tmYLTickSpacingF" : 10
   "tmYLMinorPerMajor" : 4
   "tmYLMinorLengthF" : -0.008
   "tmYLMajorLengthF" : -0.016
   "tmYLLabelFont" : 21
   "tmYLLabelFontHeightF" : 0.020

   "vpXF" : 0.11
   "vpYF" : yf
   "vpWidthF"  : wi
   "vpHeightF" : he

   "trXMinF" : 0
   "trXMaxF" : 1200
   "trYMinF" : 0
   "trYMaxF" : 60

   "tiMainString" : hour1+" h  "+min1+" m "
   "tiMainFont"  : 22
   "tiMainFontHeightF" : 0.025

   "tiXAxisString" : "X (km)"
   "tiXAxisFont"  : 21
   "tiXAxisFontHeightF" : 0.022
;   "tiXAxisOn" : True 

   "tiYAxisString" : "Z (km)"
   "tiYAxisFont"  : 21
   "tiYAxisFontHeightF" : 0.022

  end create
  draw(cont2)

;----------------------------------------------------
; qc1
;----------------------------------------------------

  setvalues cont2
   "cnScalarFieldData" : field5
   "cnFillOn"        :  False           ; color fill
   "cnLevelSelectionMode" : "ExplicitLevels"
   "cnMonoLineDashPattern" : True
   "cnLineDashPattern" : 5
   "cnLevels" : (/0.0001/)
   "cnLinesOn" : True
   "cnMonoLineThickness" : True
   "cnLineThicknessF" : 1.5
   "cnInfoLabelOn" : False
   "cnLineLabelsOn" : False
   "cnConstFLabelOn" : False
  end setvalues
  draw(cont2)

;----------------------------------------------------
; ri1
;----------------------------------------------------

;  setvalues cont2
;   "cnScalarFieldData" : field7
;   "cnFillOn"        :  False           ; color fill
;   "cnLevelSelectionMode" : "ExplicitLevels"
;   "cnLevels" : (/0.25/)
;   "cnLinesOn" : True
;   "cnLineDashPattern" : 0
;   "cnMonoLineThickness" : True
;   "cnLineThicknessF" : 1.5
;   "cnInfoLabelOn" : False
;   "cnLineLabelsOn" : False
;  end setvalues
;  draw(cont2)

;----------------------------------------------------


   frame(wks)

 if (hour1 .eq. time1/3600.) then
   print(n)
 end if

 end do

 end



