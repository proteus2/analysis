load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

xspace = 170  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
axisfontratio = 1.2

begin

 expname = (/"fd","y"/)

 nmon = 4+144
 cmon = (/"January","July"/)
 imon = 4+(/0,6/)
;--------------------------------------- read base data

   fnam = "../dat/mf8_"+expname+".nc"
   f1 = addfile(fnam(0),"r")
   f2 = addfile(fnam(1),"r")

   mfe1_1 = dim_avg_n( f1->unspecified(ispan(imon(0),nmon,12),0,:,:), 0 )
   mfe2_1 = dim_avg_n( f1->field420   (ispan(imon(0),nmon,12),0,:,:), 0 )
   mfe3_1 = dim_avg_n( f2->field420   (ispan(imon(0),nmon,12),0,:,:), 0 )
   mfe2_1 = mfe2_1 + mfe1_1

   mfe1_7 = dim_avg_n( f1->unspecified(ispan(imon(1),nmon,12),0,:,:), 0 )
   mfe2_7 = dim_avg_n( f1->field420   (ispan(imon(1),nmon,12),0,:,:), 0 )
   mfe3_7 = dim_avg_n( f2->field420   (ispan(imon(1),nmon,12),0,:,:), 0 )
   mfe2_7 = mfe2_7 + mfe1_7

   mfw1_1 = dim_avg_n( f1->unspecified_1(ispan(imon(0),nmon,12),0,:,:), 0 )*(-1.)
   mfw2_1 = dim_avg_n( f1->field422     (ispan(imon(0),nmon,12),0,:,:), 0 )
   mfw3_1 = dim_avg_n( f2->field422     (ispan(imon(0),nmon,12),0,:,:), 0 )
   mfw2_1 = mfw2_1 + mfw1_1

   mfw1_7 = dim_avg_n( f1->unspecified_1(ispan(imon(1),nmon,12),0,:,:), 0 )*(-1.)
   mfw2_7 = dim_avg_n( f1->field422     (ispan(imon(1),nmon,12),0,:,:), 0 )
   mfw3_7 = dim_avg_n( f2->field422     (ispan(imon(1),nmon,12),0,:,:), 0 )
   mfw2_7 = mfw2_7 + mfw1_7

   mfe = (/ (/mfe1_1, mfe1_7/), \
            (/mfe2_1, mfe2_7/), \
            (/mfe3_1, mfe3_7/) /)*1.e3
   mfw = (/ (/mfw1_1, mfw1_7/), \
            (/mfw2_1, mfw2_7/), \
            (/mfw3_1, mfw3_7/) /)*1.e3

   lat = f1->latitude
   lon = f1->longitude

;--------------------------------------- contour

  field = new((/2,3,2/),graphic)
  pair  = new((/2,3,2/),graphic)

  do idr=0, 1
  do iv=0, 2
  do is=0, 1

  field(0,iv,is) = create "w" scalarFieldClass defaultapp
    "sfDataArray" : mfe(iv,is,:,:)
    "sfYArray" : lat
    "sfXArray" : lon
  end create
  field(1,iv,is) = create "w" scalarFieldClass defaultapp
    "sfDataArray" : mfw(iv,is,:,:)
    "sfYArray" : lat
    "sfXArray" : lon
  end create

  pair(0,iv,is) = create "w" coordArraysClass defaultapp
    "caXArray" : dim_avg(mfe(iv,is,:,:))
    "caYArray" : lat
  end create
  pair(1,iv,is) = create "w" coordArraysClass defaultapp
    "caXArray" : dim_avg(mfw(iv,is,:,:))
    "caYArray" : lat
  end create

  end do
  end do
  end do


  sdr = (/"Eastward","Westward"/)
  do idr=0, 1

wks = cr_wks("mf100_"+idr+"-"+expname(0),600,0)
setvalues wks
 "wkColorMap" : "WhBlGrYeRe"
end setvalues
colind = NhlNewColor(wks,0.95,0.95,0.95)
colgry = NhlNewColor(wks,0.50,0.50,0.50)

 cnlev = (/0.25,0.5,1.0,2.0,4.0,8.0/)
 cnfil = (/colind, 7,11,15,20,77,95/)

  xfitv = 83
  wi = 56
  he = 38
  xf = 11 + xfitv*(/0,1/)
  yf = 30 + (8+he)*(/0,1,2/)

  xstr = (/"","","lon (deg)"/)
  ystr = (/"lat (deg)",""/)
  txtstr = sdr(idr) + (/" CGW", " CGW + NGW", " NGW (CTL)"/)

  do iv=0, 2
  do is=0, 1

  xf2 = xf(is)+wi+6
  wi2 = 10
  xstr2 = (/"","","mean (mPa)"/)

  map = create "mp" mapPlotClass wks
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
  set_size_mp(map,xf(is),yf(iv),wi,he)

  cont1 = create "contour" contourPlotClass wks
   "cnScalarFieldData" : field(idr,iv,is)
   "cnFillOn"          : True
   "cnLinesOn"         : False
   "cnLevelSelectionMode" : "ExplicitLevels"
   "cnLevels"             : cnlev
   "cnFillColors"         : cnfil
   "cnLineLabelsOn"       : False
   "cnInfoLabelOn"        : False
  end create

  overlay(map,cont1)

  draw(map)


  rim = create "contour" contourPlotClass wks
   "cnNoDataLabelOn" : False

   "tmXBMode" : "Manual"
   "tmXBTickSpacingF"  : 60
   "tmXBMinorPerMajor" : 2
   "tmYLMode" : "Manual"
   "tmYLTickSpacingF"  : 30
   "tmYLMinorPerMajor" : 2

   "tiXAxisString" : xstr(iv)
   "tiYAxisString" : ystr(is)

   "trXMinF" : 0
   "trXMaxF" : 360
   "trYMinF" : -90
   "trYMaxF" : 90
  end create
  set_size(rim,xf(is),yf(iv),wi,he)

  draw(rim)

  dr_txt2(wks,txtstr(iv),xf(is)+1.2,yf(iv)+he-2,22,6,0.,"CenterLeft",1)

  plot =  create "lplot" xyPlotClass wks
   "xyCoordData" : pair(idr,iv,is)

   "trYMaxF" : 90
   "trYMinF" : -90
   "trXMaxF" : 4
   "trXMinF" : 0

   "tmXBMode" : "MANUAL"
   "tmXBTickSpacingF" : 2
   "tmXBMinorPerMajor" : 1
   "tmXBFormat" : "@*+^sg"
   "tmYLMode" : "Manual"
   "tmYLTickSpacingF"  : 30
   "tmYLMinorPerMajor" : 2

   "tiXAxisString" : xstr2(iv)
   "tiYAxisOn" : False
  end create
  set_size(plot,xf2,yf(iv),wi2,he)

  draw(plot)


  delete(map)
  delete(rim)
  delete(plot)


  end do  ; is
  end do  ; iv

  dr_txt(wks,cmon(0),xf(0)+wi/2,yf(0)-4,22,10,0.)
  dr_txt(wks,cmon(1),xf(1)+wi/2,yf(0)-4,22,10,0.)


   label = create "lab" labelBarClass wks
     "lbOrientation"     : "horizontal"
     "lbBoxCount"        : dimsizes(cnfil)
     "lbFillColors"      : cnfil
     "lbMonoFillPattern" : True
     "lbPerimOn"         : False
     "lbLabelsOn"         : True
     "lbLabelStrings"     : cnlev
     "lbLabelPosition"    : "Bottom"
     "lbLabelAlignment"   : "InteriorEdges"
     "lbLabelStride"      : 1
     "lbTitleString"      : "(mPa)"
     "lbTitlePosition"    : "Right"
     "lbTitleDirection"   : "Across"
   end create
   set_size_lb(label,60,yf(2)+he+10,70,5)

   draw(label)


  dr_txt(wks,"(a)",xf(0)-8,yf(0)+1,25,8,0.)
  dr_txt(wks,"(b)",xf(0)-8,yf(1)+1,25,8,0.)
  dr_txt(wks,"(c)",xf(0)-8,yf(2)+1,25,8,0.)


   frame(wks)

  end do  ; idr


 end



