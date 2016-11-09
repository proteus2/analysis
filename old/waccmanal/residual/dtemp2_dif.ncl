load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"

begin

 setvalues NhlGetWorkspaceObjectId()
  "wsMaximumSize" : 55000000
 end setvalues


do i22 = 0, 0
do i33 = 1, 12


 ;  0-EPD  ;  1-BGWD  ;  2-GWDO  ;  3-GWDC
 idvar = i22;0
 cvar = (/"EPD","BGWD","GWDO","GWDC"/)

 mon = i33;1
 cmon = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)

 cnlev = (/-5.,-2.,-1.,-0.5,-0.2,-0.1,-0.05, \
           0.,0.05,0.1,0.2,0.5,1.,2.,5./)

 cnfil = (/ 5,5,5,6,7,8,9, 0, 0 ,14,15,16,17,18,18,18/)

 mainstr_tail = " _ "+cvar(idvar)+"  ("+cmon(mon-1)+", LSGWDC-CTL)"
 mainstr = (/"dT[V*,W*] "/)+mainstr_tail

 ofname = "dtemp2_"+cvar(idvar)+"-"+sprinti("%2.2i",mon)+"-dif"

 xstr = (/"","lat (deg)"/)

 shading = (/True,False/)
 contour = (/False,True/)

;--------------------------------------- read base data

 fdir = "res/"

 fnam1 = "control1.nc"
 fnam2 = "lsgwdclm.nc"

 missv = 1.e20

 f1 = addfile(fdir+fnam1,"r")
 f2 = addfile(fdir+fnam2,"r")

 lat= f1->lat
 z  = f1->z

 ny = dimsizes(lat)
 nz = dimsizes(z)

 if (idvar .eq. 0) then
   var1v = f1->dtvepd(mon-1,:,:)
   var2v = f2->dtvepd(mon-1,:,:)
   var1w = f1->dtwepd(mon-1,:,:)
   var2w = f2->dtwepd(mon-1,:,:)
 end if
 if (idvar .eq. 1) then
   var1v = f1->dtvgwd(mon-1,:,:)
   var2v = f2->dtvgwd(mon-1,:,:)
   var1w = f1->dtwgwd(mon-1,:,:)
   var2w = f2->dtwgwd(mon-1,:,:)
 end if
 if (idvar .eq. 2) then
   var1v = f1->dtvgwo(mon-1,:,:)
   var2v = f2->dtvgwo(mon-1,:,:)
   var1w = f1->dtwgwo(mon-1,:,:)
   var2w = f2->dtwgwo(mon-1,:,:)
 end if
 if (idvar .eq. 3) then
   var2v = f2->dtvsgw(mon-1,:,:)
   var1v = var2v*0.
   var2w = f2->dtwsgw(mon-1,:,:)
   var1w = var2w*0.
 end if


 vardraw = new((/2,nz,ny/),float,missv)

 vardraw(0,:,:) = (/var1v+var1w/)
 vardraw(1,:,:) = (/var2v+var2w-var1v-var1w/)

;---------------------------------------------------

 wks = create "ps" psWorkstationClass defaultapp
     "wkPSFileName"   : ofname+".eps"
     "wkPSFormat"     : "eps"
     "wkPSResolution" : 600
     "wkColorMap"        : "posneg_2"
     "wkBackgroundColor" : (/1.,1.,1./)
 end create
 colind = NhlNewColor(wks,0.5,0.5,0.5)

;--------------------------------------- cn level

; print("contour levels :  "+cnlev+ " (K/day)")

; cnfil = new(dimsizes(cnlev)+1,integer,0)
; do i=0, dimsizes(cnlev)-1
;   cnfil(i) = 6
;   if (cnlev(i) .gt. -1.) then
;     cnfil(i) = 8
;   end if
;   if (cnlev(i) .gt. -0.1) then
;     cnfil(i) = 9
;   end if
;   if (cnlev(i) .gt. -0.05) then
;     cnfil(i) = 0
;   end if
;   if (cnlev(i) .gt. 0.05) then
;     cnfil(i) = 14
;   end if
;   if (cnlev(i) .gt. 0.1) then
;     cnfil(i) = 15
;   end if
;   if (cnlev(i) .gt. 1.) then
;     cnfil(i) = 17
;   end if
; end do
; cnfil(dimsizes(cnlev)) = 17

;--------------------------------------- field

  field = new(2,graphic)

  do j=0, 1
  field(j) = create "w" scalarFieldClass defaultapp
    "sfDataArray" : vardraw(j,:,:) * 86400.
    "sfYArray" : z
    "sfXArray" : lat
    "sfMissingValueV" : missv
  end create
  end do

;--------------------------------------- contour
; wprt1  
;---------------------------------------

  do j=0, 1
  do n=0, 0

  cont1 = create "contour" contourPlotClass wks
   "cnScalarFieldData" : field(j)
   "cnFillOn"          : shading(j)
   "cnLinesOn"         : contour(j)
   "cnLevelSelectionMode" : "ExplicitLevels"
   "cnLevels"             : cnlev
   "cnFillColors"         : cnfil
   "cnLineLabelsOn"           : True
   "cnLineLabelPlacementMode" : "Constant"
   "cnLineLabelFont"          : 21
   "cnLineLabelFontHeightF"   : 0.010
;   "cnLineLabelFontColor"     : colind
   "cnInfoLabelOn"            : False

   "tmXBMode" : "Manual"
   "tmXBTickSpacingF"  : 30
   "tmXBMinorPerMajor" : 2
   "tmXBMajorLengthF"  : -0.008
   "tmXBMinorLengthF"  : -0.005
   "tmXBLabelFont"        : 21
   "tmXBLabelFontHeightF" : 0.014

   "tmYLMode" : "Manual"
   "tmYLTickSpacingF"  : 10
   "tmYLMinorPerMajor" : 1
   "tmYLMajorLengthF"  : -0.008
   "tmYLMinorLengthF"  : -0.005
   "tmYLLabelFont"        : 21
   "tmYLLabelFontHeightF" : 0.014

   "tiMainString"       : mainstr(n)
   "tiMainFont"         : 22
   "tiMainFontHeightF"  : 0.018

   "tiXAxisString"      : xstr(n)
   "tiXAxisFont"        : 21
   "tiXAxisFontHeightF" : 0.017

   "tiYAxisString"      : "height (km)"
   "tiYAxisFont"        : 21
   "tiYAxisFontHeightF" : 0.017

    "vpXF" : 0.2
    "vpYF" : 0.90 - n*0.4
    "vpWidthF"  : 0.4
    "vpHeightF" : 0.3
  end create
  getvalues cont1
    "cnLevels" : levels
    "vpXF" : xfr
    "vpYF" : yfr
    "vpWidthF"  : wir
    "vpHeightF" : her
  end getvalues
  if (j .eq. 1) then
  ncr = dimsizes(levels)
  ndash = new(ncr,integer)
  ncol  = new(ncr,integer)
   do i=0, ncr-1
     if (levels(i) .ge. 0.0) then
       ndash(i) = 0
     else
       ndash(i) = 11
     end if
     if (levels(i) .eq. 0.0) then
       ncol(i) = colind
     else
       ncol(i) = 1
     end if
   end do
  setvalues cont1
    "cnMonoLineDashPattern" : False
    "cnLineDashPatterns"    : ndash
    "cnMonoLineColor"       : False
    "cnLineColors"          : ncol
    "cnLineThicknessF"      : 1.8
  end setvalues
  end if

  plotover = create "plot" logLinPlotClass wks
    "trXMinF" : -90
    "trXMaxF" :  90
    "trYMinF" : 0
    "trYMaxF" : 40

    "vpXF" : xfr
    "vpYF" : yfr
    "vpWidthF"  : wir
    "vpHeightF" : her
  end create

  overlay(plotover,cont1)

  draw(plotover)

  delete(plotover)

  end do
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
     "vpYF"      : 0.12
     "vpWidthF"  : 0.45
     "vpHeightF" : 0.03
   end create
;   draw(label)

   frame(wks)


end do
end do


 end

