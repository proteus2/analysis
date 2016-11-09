load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"

begin

 setvalues NhlGetWorkspaceObjectId()
  "wsMaximumSize" : 55000000
 end setvalues


do i11 = 0, 2
do i22 = 1, 1 ;3
do i33 = 1, 12


 ;  0-CTL  ;  1-SGWDC  ;  2-LSGWDC
 idexp = i11;0
 cexp = (/"CTL","SGWDC","LSGWDC"/)

 ;  0-EPD  ;  1-BGWD  ;  2-GWDO  ;  3-GWDC
 idvar = i22;0
 cvar = (/"EPD","BGWD","GWDO","GWDC"/)

 mon = i33;1
 cmon = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)

 if (idvar .eq. 0) then
   dcn = 200.
 end if
 if (idvar .eq. 1) then
   dcn = 20.
 end if
 if (idvar .eq. 2) then
   dcn = 20.
 end if
 if (idvar .eq. 3) then
   dcn = 20.
 end if

 mainstr = "Chi _ "+cvar(idvar)+"  ("+cmon(mon-1)+", "+cexp(idexp)+")"

 ofname = "chi_"+cvar(idvar)+"-"+sprinti("%2.2i",mon)+"-"+cexp(idexp)

;--------------------------------------- read base data

 fdir = "res/"

 if (idexp .eq. 0) then
   fnam1 = "control1.nc"
 end if
 if (idexp .eq. 1) then
   fnam1 = "sgwdclm2.nc"
 end if
 if (idexp .eq. 2) then
   fnam1 = "lsgwdclm.nc"
 end if

 missv = 1.e20

 f1 = addfile(fdir+fnam1,"r")

 lat= f1->lat
 z  = f1->z

 ny = dimsizes(lat)
 nz = dimsizes(z)

 if (idvar .eq. 0) then
   var1 = f1->chiepd(mon-1,:,:)
 end if
 if (idvar .eq. 1) then
   var1 = f1->chigwd(mon-1,:,:)
 end if
 if (idvar .eq. 2) then
   var1 = f1->chigwo(mon-1,:,:)
 end if
 if (idvar .eq. 3) then
   if (idexp .ne. 0) then
     var1 = f1->chigwc(mon-1,:,:)
   else
     var1 = f1->chigwd(mon-1,:,:)*0.
   end if
 end if


 vardraw = new((/1,nz,ny/),float,missv)
 vardraw(0,:,:) = (/var1/)

;---------------------------------------------------

 wks = create "ps" psWorkstationClass defaultapp
     "wkPSFileName"   : ofname+".eps"
     "wkPSFormat"     : "eps"
     "wkPSResolution" : 600
     "wkColorMap"        : "rbw_l0_128"
     "wkBackgroundColor" : (/1.,1.,1./)
 end create
 colind = NhlNewColor(wks,0.5,0.5,0.5)

;--------------------------------------- cn level

 ncn = 80
 cnlev = new(ncn,float)
 do i=0, ncn-1
   cnlev(i) = (i+0.5-ncn/2.)*dcn
 end do
 print("contour interval :  "+dcn+" m^2/s")

 cnfil = new(ncn+1,integer,0)
; cnfil = (/ 48, 48, 48, 48, 48,  0,  0,  0,  0,  0, \
;             0,  0,  0,  0,  0, 82, 82, 82, 82, 82, \
;            88, 88, 88, 88, 88, 95/)

;--------------------------------------- field

  field = new(1,graphic)

  do n=0, 0
  field(n) = create "w" scalarFieldClass defaultapp
    "sfDataArray" : vardraw(n,:,:)
    "sfYArray" : z
    "sfXArray" : lat
    "sfMissingValueV" : missv
  end create
  end do

;--------------------------------------- contour
; wprt1  
;---------------------------------------

  do n=0, 0

  cont1 = create "contour" contourPlotClass wks
   "cnScalarFieldData" : field(n)
;   "cnFillOn"          : True
   "cnLinesOn"         : True
   "cnLevelSelectionMode" : "ExplicitLevels"
   "cnLevels"             : cnlev
;   "cnFillColors"         : cnfil
   "cnLineLabelsOn"           : False ;True
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
   "tmXBLabelFontHeightF" : 0.014

   "tmYLMode" : "Manual"
   "tmYLTickSpacingF"  : 10
   "tmYLMinorPerMajor" : 1
   "tmYLMajorLengthF"  : -0.008
   "tmYLMinorLengthF"  : -0.005
   "tmYLLabelFont"        : 21
   "tmYLLabelFontHeightF" : 0.014

   "tiMainString"       : mainstr
   "tiMainFont"         : 22
   "tiMainFontHeightF"  : 0.018

   "tiXAxisString"      : "lat (deg)"
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
  ncr = dimsizes(levels)
  ndash = new(ncr,integer)
  ncol  = new(ncr,integer)
   do i=0, ncr-1
     if (levels(i) .ge. 0.0) then
       ndash(i) = 0
     else
       ndash(i) = 13
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
  end setvalues

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


end do
end do
end do


 end

