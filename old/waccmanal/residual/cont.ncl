load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"

begin

 setvalues NhlGetWorkspaceObjectId()
  "wsMaximumSize" : 55000000
 end setvalues


 mon = 1

;--------------------------------------- read base data

   fdir = "res/"

   fnam1 = "control1.nc"
   fnam2 = "sgwdclm2.nc"

   missv = 1.e20

   f1 = addfile(fdir+fnam1,"r")
   lat= f1->lat
   z  = f1->z

   ny = dimsizes(lat)
   nz = dimsizes(z)

   var1v = f1->vsepd(mon-1,:,:)
   var1w = f1->wsepd(mon-1,:,:)*1.e3

   f2 = addfile(fdir+fnam2,"r")
   var2v = f2->vsepd(mon-1,:,:)
   var2w = f2->wsepd(mon-1,:,:)*1.e3

   varv_d = var2v - var1v
   varw_d = var2w - var1w


  vardraw = new((/2,nz,ny/),float,missv)
  vardraw(0,:,:) = (/var1v/)
  vardraw(1,:,:) = (/var1w/)

;---------------------------------------------------
 wks = create "ps" psWorkstationClass defaultapp
     "wkPSFileName"   : "cont.eps"
     "wkPSFormat"     : "eps"
     "wkPSResolution" : 600
     "wkColorMap"        : "rbw_l0_128";"BlWhRe"
     "wkBackgroundColor" : (/1.,1.,1./)
 end create
 colind = NhlNewColor(wks,0.3,0.3,0.3)

;--------------------------------------- cn level

; ncn = 25
; cnlev = (/-9,-8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10, \
;           11,12,13,14,15/)*0.1
;
; cnfil = new(ncn+1,integer)
; cnfil = (/ 48, 48, 48, 48, 48,  0,  0,  0,  0,  0, \
;             0,  0,  0,  0,  0, 82, 82, 82, 82, 82, \
;            88, 88, 88, 88, 88, 95/)

 ncn = 17
 cnlev = (/-10,-5,-2,-1,-0.5,-0.2,-0.1,-0.05,0, \
           0.05,0.1,0.2,0.5,1,2,5,10/)
 cnfil = new(ncn+1,integer)
 cnfil = 0

;--------------------------------------- field

  field = new(2,graphic)

  do n=0, 1
  field(n) = create "w" scalarFieldClass defaultapp
    "sfDataArray" : vardraw(n,:,:)
    "sfYArray" : z
    "sfXArray" : lat
    "sfMissingValueV" : missv
  end create
  end do
;print( max(uavg(0:43,:)) )
;print( min(uavg(0:43,:)) )

;--------------------------------------- contour
; wprt1  
;---------------------------------------

  do n=0, 1

  cont1 = create "contour" contourPlotClass wks
   "cnScalarFieldData" : field(n)
   "cnFillOn"          : True
   "cnLinesOn"         : True
   "cnLevelSelectionMode" : "ExplicitLevels"
   "cnLevels"             : cnlev
   "cnFillColors"         : cnfil
   "cnLineLabelsOn"           : True
   "cnLineLabelPlacementMode" : "Constant"
   "cnLineLabelFont"          : 21
   "cnLineLabelFontHeightF"   : 0.011
   "cnInfoLabelOn"            : False

   "tmXBMode" : "Manual"
   "tmXBTickSpacingF"  : 30
   "tmXBMinorPerMajor" : 2
   "tmXBMajorLengthF"  : -0.008
   "tmXBMinorLengthF"  : -0.005
   "tmXBLabelFont"        : 21
   "tmXBLabelFontHeightF" : 0.015

   "tmYLMode" : "Manual"
   "tmYLTickSpacingF"  : 10
   "tmYLMinorPerMajor" : 1
   "tmYLMajorLengthF"  : -0.008
   "tmYLMinorLengthF"  : -0.005
   "tmYLLabelFont"        : 21
   "tmYLLabelFontHeightF" : 0.015

   "tiMainFont"         : 22
   "tiMainFontHeightF"  : 0.022

   "tiXAxisString"      : "lat (deg)"
   "tiXAxisFont"        : 21
   "tiXAxisFontHeightF" : 0.020

   "tiYAxisString"      : "height (km)"
   "tiYAxisFont"        : 21
   "tiYAxisFontHeightF" : 0.020

    "vpXF" : 0.2
    "vpYF" : 0.95
    "vpWidthF"  : 0.5
    "vpHeightF" : 0.3
  end create
  getvalues cont1
    "cnLevels" : levels
  end getvalues
  nc = dimsizes(levels)
  ndash = new(nc,integer)
   do i=0,nc-1
     if (levels(i) .ge. 0.0) then
       ndash(i) = 0
     else
       ndash(i) = 1
     end if
   end do
  setvalues cont1
    "cnMonoLineDashPattern" : False
    "cnLineDashPatterns" : ndash
;    "cnLineThicknessF"   : 1.8
  end setvalues
  plotover = create "plot" logLinPlotClass wks
    "vpXF" : 0.2
    "vpYF" : 0.95 - n*0.4
    "vpWidthF"  : 0.4
    "vpHeightF" : 0.3

    "trXMinF" : -90
    "trXMaxF" :  90
    "trYMinF" : 0
    "trYMaxF" : 40;120
  end create

  overlay(plotover,cont1)

  draw(plotover)

  delete(plotover)

  end do


   label = create "lab" labelBarClass wks
     "lbAutoManage"      : False
     "lbOrientation"     : "horizontal"
     "lbBoxCount"        : dimsizes(cnfil)
     "lbFillColors"      : cnfil
     "lbMonoFillPattern" : True
     "lbPerimOn"         : False
     "lbLabelsOn"         : True
     "lbLabelStrings"     : cnlev
     "lbLabelFont"        : 21
     "lbLabelFontHeightF" : 0.011
     "lbLabelPosition"    : "Bottom";,Right"
     "lbLabelAlignment"   : "InteriorEdges"
     "lbLabelAngleF"      : 0.
     "lbLabelStride"      : 2
     "lbTitleString"      : "(m s:S:-1:N:)"
     "lbTitlePosition"    : "Right";,"Top"
     "lbTitleDirection"   : "Across"
     "lbTitleFont"        : 21
     "lbTitleFontHeightF" : 0.012
     "vpXF"      : 0.15
     "vpYF"      : 0.55
     "vpWidthF"  : 0.60
     "vpHeightF" : 0.03
   end create
;   draw(label)

   frame(wks)


 end

