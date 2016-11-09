load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"

begin

 setvalues NhlGetWorkspaceObjectId()
  "wsMaximumSize" : 55000000
 end setvalues


do i11 = 2, 2
do i22 = 0, 3


 ;  0-CTL  ;  1-SGWDC  ;  2-LSGWDC
 idexp = i11;0
 cexp = (/"CTL","SGWDC","LSGWDC"/)

 ;  0-EPD  ;  1-BGWD  ;  2-GWDO  ;  3-GWDC
 idvar = i22;0
 cvar = (/"EPD","BGWD","GWDO","GWDC"/)

 ;  7- 13 km  ;  8- 15 km  ;  9- 17 km
 iheight = 8
 cheight = (/"13","15","17"/)

 mon = (/0,1,2,3,4,5,6,7,8,9,10,11,12,13/)
 cmon = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)

 cnlev = (/-5.,-2.,-1.,-0.5,-0.2,-0.1,-0.05, \
           0.,0.05,0.1,0.2,0.5,1.,2.,5./)

 if (idvar .ne. 0) then
   cnlev = cnlev * 0.1
 end if


 mainstr_tail = " _ "+cvar(idvar)+"  ("+cheight(iheight-7)+" km, "+cexp(idexp)+")"
 mainstr = (/"dT[V*,W*] "/)+mainstr_tail

 ofname = "dtemp"+cheight(iheight-7)+"_t_"+cvar(idvar)+"-"+cexp(idexp)

 xstr = (/"month"/)

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

 ny = dimsizes(lat)

 if (idvar .eq. 0) then
   var1v = f1->dtvepd(:,iheight,:)
   var1w = f1->dtwepd(:,iheight,:)
 end if
 if (idvar .eq. 1) then
   var1v = f1->dtvgwd(:,iheight,:)
   var1w = f1->dtwgwd(:,iheight,:)
 end if
 if (idvar .eq. 2) then
   var1v = f1->dtvgwo(:,iheight,:)
   var1w = f1->dtwgwo(:,iheight,:)
 end if
 if (idvar .eq. 3) then
   if (idexp .ne. 0) then
     var1v = f1->dtvgwc(:,iheight,:)
     var1w = f1->dtwgwc(:,iheight,:)
   else
     var1v = f1->dtvgwd(:,iheight,:)*0.
     var1w = f1->dtwgwd(:,iheight,:)*0.
   end if
 end if


 vardraw = new((/1,ny,14/),float,missv)
 do j=0, ny-1
 do n=1, 12
   vardraw(0,j,n) = (/var1v(n-1,j)+var1w(n-1,j)/)
 end do
 end do
 vardraw(0,:,0) = vardraw(0,:,12)
 vardraw(0,:,13) = vardraw(0,:,1)

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

 cnfil = new(dimsizes(cnlev)+1,integer,0)
 do i=0, dimsizes(cnlev)-1
   cnfil(i) = 6
   if (cnlev(i) .gt. -1.) then
     cnfil(i) = 8
   end if
   if (cnlev(i) .gt. -0.1) then
     cnfil(i) = 9
   end if
   if (cnlev(i) .gt. -0.05) then
     cnfil(i) = 0
   end if
   if (cnlev(i) .gt. 0.05) then
     cnfil(i) = 14
   end if
   if (cnlev(i) .gt. 0.1) then
     cnfil(i) = 15
   end if
   if (cnlev(i) .gt. 1.) then
     cnfil(i) = 17
   end if
 end do
 cnfil(dimsizes(cnlev)) = 17

;--------------------------------------- field

  field = new(1,graphic)

  do n=0, 0
  field(n) = create "w" scalarFieldClass defaultapp
    "sfDataArray" : vardraw(n,:,:) * 86400.
    "sfYArray" : lat
    "sfXArray" : mon
    "sfMissingValueV" : missv
  end create
  end do

;--------------------------------------- contour
; wprt1  
;---------------------------------------

  do n=0, 0

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
   "cnLineLabelFontHeightF"   : 0.010
;   "cnLineLabelFontColor"     : colind
   "cnInfoLabelOn"            : False

   "tmXBMode" : "Manual"
   "tmXBTickSpacingF"  : 1
   "tmXBMinorPerMajor" : 0
   "tmXBMajorLengthF"  : -0.008
   "tmXBMinorLengthF"  : -0.005
   "tmXBLabelFont"        : 21
   "tmXBLabelFontHeightF" : 0.014
   "tmXBTickStartF"       : 1
   "tmXBFormat"           : "@*+^sg"

   "tmYLMode" : "Manual"
   "tmYLTickSpacingF"  : 30
   "tmYLMinorPerMajor" : 2
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

   "tiYAxisString"      : "lat (deg)"
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
  end setvalues

  plotover = create "plot" logLinPlotClass wks
    "trXMinF" : 0.5
    "trXMaxF" : 12.5
    "trYMinF" : -90
    "trYMaxF" : 90

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


 end

