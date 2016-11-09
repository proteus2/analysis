load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin

 setvalues NhlGetWorkspaceObjectId()
  "wsMaximumSize" : 55000000
 end setvalues


 expname = "a"

 nmon = 4+24 ;84 ;24

;--------------------------------------- read base data
   fnam1 = "../dat/mfct_"+expname+".nc"

   f1  = addfile(fnam1,"r")
   mfect1 = f1->unspecified  (ispan(4 ,nmon-1,12),0,:,:)
   mfwct1 = f1->unspecified_1(ispan(4 ,nmon-1,12),0,:,:) *(-1.)
   mfect2 = f1->unspecified  (ispan(10,nmon-1,12),0,:,:)
   mfwct2 = f1->unspecified_1(ispan(10,nmon-1,12),0,:,:) *(-1.)
   lat = f1->latitude
   lon = f1->longitude

   fnam2 = "../dat/mfct0_"+expname+".nc"

   f2 = addfile(fnam2,"r")
   mfeun1 = f2->unspecified  (ispan(4 ,nmon-1,12),0,:,:)
   mfwun1 = f2->unspecified_1(ispan(4 ,nmon-1,12),0,:,:) *(-1.)
   mfeun2 = f2->unspecified  (ispan(10,nmon-1,12),0,:,:)
   mfwun2 = f2->unspecified_1(ispan(10,nmon-1,12),0,:,:) *(-1.)

   if (dimsizes(ispan(4 ,nmon-1,12)) .ne. 1) then
     mf_avg = dim_avg_n_Wrap( \
               (/ (/mfect1, mfwct1, mfect2, mfwct2/), \
                  (/mfeun1, mfwun1, mfeun2, mfwun2/) /), 2 )*1.e3
   else
     mf_avg = (/ (/mfect1, mfwct1, mfect2, mfwct2/), \
                 (/mfeun1, mfwun1, mfeun2, mfwun2/) /)*1.e3
   end if
mf_avg = mf_avg*3.

   mf_zavg = dim_avg(mf_avg)

;---------------------------------------------------
 wks = create "ps" psWorkstationClass defaultapp
     "wkPSFileName"   : "momfall.eps"
     "wkPSFormat"     : "eps"
     "wkPSResolution" : 600
     "wkColorMap"        : "WhBlGrYeRe" ;"rbw_l0_128"
     "wkBackgroundColor" : (/1.,1.,1./)
 end create
 colind = NhlNewColor(wks,0.95,0.95,0.95)
 colind2 = NhlNewColor(wks,0.50,0.50,0.50)

;--------------------------------------- cn level
 cnlev = (/0.2,0.5,1.0,2.0,5.0,10.0/)*0.6
 cnfil = (/colind, 7,11,15,20,77,95/)

;--------------------------------------- field
  field = new((/2,4/),graphic)
  do np=0, 1
  do n=0, 3
  field(np,n) = create "w" scalarFieldClass defaultapp
    "sfDataArray" : mf_avg(np,n,:,:)
    "sfYArray" : lat
    "sfXArray" : lon
  end create
  end do
  end do

  pair = new((/2,4/),graphic)
  do np=0, 1
  do n=0, 3
  pair(np,n) = create "w" coordArraysClass defaultapp
    "caXArray" : mf_zavg(np,n,:)
    "caYArray" : lat
  end create
  end do
  end do

;--------------------------------------- contour
; wprt1  
;---------------------------------------

  xf = (/0.1,0.45,0.80/)
  yf = (/0.95,0.72,0.48,0.25/)
  wi = 0.3
  he = 0.16
  mstr = (/ (/"E-MF (JAN)", "W-MF (JAN)", "E-MF (JUL)", "W-MF (JUL)"/), \
            (/"(unfiltered)", "", "", ""/) /)
  lcol = (/1,colind2/)

do np=0, 1

  do n=0, 3

  map = create "mp" mapPlotClass wks
   "vpXF"      : xf(np)
   "vpYF"      : yf(n)
   "vpWidthF"  : wi
   "vpHeightF" : he

;   "mpProjection"    : "LambertConformal"
   "mpFillOn"        : False
;   "mpLabelsOn"      : True;False
   "mpOutlineOn"     : True
;   "mpGeophysicalLineColor" : "Orange4"
   "mpGridAndLimbOn" : False
   "mpPerimOn"       : True

   "mpShapeMode"       : "FreeAspect"
   "mpDataBaseVersion" : "Ncarg4_0"

   "mpLimitMode"    : "LatLon"
   "mpCenterLonF"   : 180.
   "mpMinLonF"      : 0.
   "mpMaxLonF"      : 360.
   "mpMinLatF"      : -90.
   "mpMaxLatF"      : 90.
  end create

  cont1 = create "contour" contourPlotClass wks
   "cnScalarFieldData" : field(np,n)
   "cnFillOn"          : True
   "cnLinesOn"         : False
   "cnLevelSelectionMode" : "ExplicitLevels"
   "cnLevels"             : cnlev 
   "cnFillColors"         : cnfil
   "cnLineLabelsOn"           : True
   "cnLineLabelPlacementMode" : "Constant"
   "cnLineLabelFont"          : 21
   "cnLineLabelFontHeightF"   : 0.011
   "cnInfoLabelOn"            : False
  end create

  overlay(map,cont1)

  draw(map)


  rim = create "contour" contourPlotClass wks
   "cnNoDataLabelOn" : False

   "tmXBMode" : "Manual"
   "tmXBTickSpacingF"  : 60
   "tmXBMinorPerMajor" : 2
   "tmXBMajorLengthF"  : -0.005
   "tmXBMinorLengthF"  : -0.003
   "tmXBLabelFont"        : 21
   "tmXBLabelFontHeightF" : 0.009

   "tmYLMode" : "Manual"
   "tmYLTickSpacingF"  : 30
   "tmYLMinorPerMajor" : 2
   "tmYLMajorLengthF"  : -0.005
   "tmYLMinorLengthF"  : -0.003
   "tmYLLabelFont"        : 21
   "tmYLLabelFontHeightF" : 0.009

   "tiMainFont"         : 22
   "tiMainFontHeightF"  : 0.015
   "tiMainString"       : mstr(np,n)

   "tiXAxisString"      : ""
   "tiXAxisFont"        : 21
   "tiXAxisFontHeightF" : 0.013

   "tiYAxisString"      : ""
   "tiYAxisFont"        : 21
   "tiYAxisFontHeightF" : 0.013

   "vpXF"      : xf(np)
   "vpYF"      : yf(n)
   "vpWidthF"  : wi
   "vpHeightF" : he

   "trXMinF" : 0
   "trXMaxF" : 360
   "trYMinF" : -90
   "trYMaxF" : 90
  end create
  draw(rim)

  np_op = np-1
  if (np .eq. 0) then
    np_op = 1
  end if
  plot = create "lplot" xyPlotClass wks
   "xyCoordData" : pair(np_op,n)

   "trYMaxF" : 90
   "trYMinF" : -90
   "trXMaxF" : 4
   "trXMinF" : 0

   "tmXBMode" : "MANUAL"
   "tmXBTickSpacingF" : 1
   "tmXBMinorPerMajor" : 1
   "tmXBLabelFont" : 21
   "tmXBLabelFontHeightF" : 0.009
   "tmXBMajorLengthF" : -0.005
   "tmXBMinorLengthF" : -0.003
   "tmXBFormat" : "@*+^sg"

   "tmYLMode" : "Manual"
   "tmYLTickSpacingF"  : 30
   "tmYLMinorPerMajor" : 2
   "tmYLMajorLengthF"  : -0.005
   "tmYLMinorLengthF"  : -0.003
   "tmYLLabelFont"        : 21
   "tmYLLabelFontHeightF" : 0.009

   "tiMainOn"  : False
   "tiXAxisFont" : 21
   "tiXAxisFontHeightF" : 0.012
   "tiXAxisString" : "(mPa)"
   "tiYAxisOn" : False

   "vpYF" : yf(n)
   "vpXF" : xf(2)
   "vpHeightF" : he
   "vpWidthF"  : 0.12
  end create
  getvalues plot
   "xyCoordDataSpec" : spec
  end getvalues
  setvalues spec
   "xyLineColor" : lcol(np_op)
  end setvalues
  draw(plot)


  delete(map)
  delete(rim)
  delete(plot)

  end do

end do  ; np


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
     "lbLabelFontHeightF" : 0.010
     "lbLabelPosition"    : "Bottom"
     "lbLabelAlignment"   : "InteriorEdges"
     "lbLabelStride"      : 1
     "lbTitleString"      : "(mPa)"
     "lbTitlePosition"    : "Right"
     "lbTitleDirection"   : "Across"
     "lbTitleFont"        : 21
     "lbTitleFontHeightF" : 0.012
     "vpXF"      : (xf(0)+xf(1))*0.5
     "vpYF"      : yf(3)-he-0.03
     "vpWidthF"  : wi*1.5
     "vpHeightF" : 0.03
   end create
   draw(label)


   frame(wks)


 end



