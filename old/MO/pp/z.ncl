load "$UTIL/ncl/header.ncl"

xspace = 170  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 1.2

begin

 expname = (/"y","y"/)

 nmon = 4+144
 cmon = (/"Jan","Jul"/)
 imon = 4+(/0,6/)
;--------------------------------------- read base data

   fnam = "../dat/mf8_"+expname+".nc"
   f1 = addfile(fnam(0),"r")
   f2 = addfile(fnam(1),"r")

   mfe1_1 = dim_avg_n( f1->unspecified(ispan(imon(0),nmon-1,12),0,:,:), 2 )
   mfe2_1 = dim_avg_n( f1->field420   (ispan(imon(0),nmon-1,12),0,:,:), 2 )
   mfe3_1 = dim_avg_n( f2->field420   (ispan(imon(0),nmon-1,12),0,:,:), 2 )
;   mfe2_1 = mfe2_1 + mfe1_1

   mfe1_7 = dim_avg_n( f1->unspecified(ispan(imon(1),nmon-1,12),0,:,:), 2 )
   mfe2_7 = dim_avg_n( f1->field420   (ispan(imon(1),nmon-1,12),0,:,:), 2 )
   mfe3_7 = dim_avg_n( f2->field420   (ispan(imon(1),nmon-1,12),0,:,:), 2 )
;   mfe2_7 = mfe2_7 + mfe1_7

   mfw1_1 = dim_avg_n( f1->unspecified_1(ispan(imon(0),nmon-1,12),0,:,:), 2 )*(-1.)
   mfw2_1 = dim_avg_n( f1->field422     (ispan(imon(0),nmon-1,12),0,:,:), 2 )
   mfw3_1 = dim_avg_n( f2->field422     (ispan(imon(0),nmon-1,12),0,:,:), 2 )
;   mfw2_1 = mfw2_1 + mfw1_1

   mfw1_7 = dim_avg_n( f1->unspecified_1(ispan(imon(1),nmon-1,12),0,:,:), 2 )*(-1.)
   mfw2_7 = dim_avg_n( f1->field422     (ispan(imon(1),nmon-1,12),0,:,:), 2 )
   mfw3_7 = dim_avg_n( f2->field422     (ispan(imon(1),nmon-1,12),0,:,:), 2 )
;   mfw2_7 = mfw2_7 + mfw1_7

   mfe = (/ (/mfe1_1, mfe1_7/), \
            (/mfe2_1, mfe2_7/), \
            (/mfe3_1, mfe3_7/) /)*1.e3
   mfw = (/ (/mfw1_1, mfw1_7/), \
            (/mfw2_1, mfw2_7/), \
            (/mfw3_1, mfw3_7/) /)*1.e3

   lat = f1->latitude
   lon = f1->longitude
   time = ispan(0,11,1)

;--------------------------------------- contour

  field = new((/2,3,2/),graphic)
  pair  = new((/2,3,2/),graphic)
print(max(mfe(0,0,:,:))+"  "+max(mfe(0,1,:,:)))
print(max(mfw(0,0,:,:))+"  "+max(mfw(0,1,:,:)))

  do idr=0, 1
  do iv=0, 2
  do is=0, 1

;  pair(0,iv,is) = create "w" coordArraysClass defaultapp
;    "caXArray" : dim_avg(mfe(iv,is,:,:))
;    "caYArray" : time
;  end create
;  pair(1,iv,is) = create "w" coordArraysClass defaultapp
;    "caXArray" : dim_avg(mfw(iv,is,:,:))
;    "caYArray" : time
;  end create

  end do
  end do
  end do


  do is=0, 1

wks = cr_wks("fig/z_"+cmon(is)+"-"+expname(0),600,0)
gsn_define_colormap(wks,"StepSeq25")
colind = NhlNewColor(wks,0.4,0.4,0.4)

; cnlev = (/0.25,0.5,1.0,2.0,4.0,8.0/)
; cnfil = (/colind, 7,11,15,20,77,95/)
 cnlev = (/2 ,4 ,6,8,10,20,40,60,80,100,120,140,160,180,200,220,240,260/)*0.01
 cnfil = (/11,10,9,8,7 ,6 ,5 ,4 ,3 ,2  ,26 ,25 ,24 ,23 ,22 ,21 ,19 ,18 ,17/)

  xfitv = 82
  wi = 56
  he = 38
  xf = 12 + xfitv*(/0,1/)
  yf = 30 + (9+he)*(/0,1,2/)

  mstr = (/ (/"Eastward-MF","",""/), (/"Westward-MF","",""/) /)
  xstr = (/"","","lon (deg)"/)
  ystr = (/"lat (deg)",""/)

  do iv=0, 2
  do idr=0, 1

  xf2 = xf(idr)+wi+7
  wi2 = 10
  xstr2 = (/"","","avg. (mPa)"/)
  txtstr = (/"CGW (SC)", "CGW + BGW (SC)", "BGW (CTL)"/)

 if (idr .eq. 0) then
   contour1 = cnshadeplot(wks,lat,time,mfe(iv,is,:,:),cnlev,cnfil,"")
 else
   contour1 = cnshadeplot(wks,lat,time,mfw(iv,is,:,:),cnlev,cnfil,"")
 end if
 setvalues contour1
  "pmLabelBarOrthogonalPosF" : 0.06
  "pmLabelBarParallelPosF"   : 0.7
  "pmLabelBarKeepAspect"     : False
  "pmLabelBarWidthF"         : 2.2
  "pmLabelBarHeightF"        : 0.12
  "lbLabelAutoStride" : False
  "lbLabelFontHeightF" : 0.03
  "lbTitlePosition"    : "Right"
  "lbTitleDirection"   : "Across"
  "lbTitleOffsetF"     : -0.02
  "lbTitleFontHeightF" : 0.04
 end setvalues

 set_size(contour1,xf(idr),yf(iv),wi,he)
 axis_range(contour1,-30,30,-999,999,False,False)
 axis_str(contour1,mstr(idr,iv),"lat (deg)",ystr(idr))

 cnfinalize(contour1,wks)

 draw(contour1)

  dr_txt2(wks,txtstr(iv),xf(idr)+1.2,yf(iv)+he-2,22,6,0.,"CenterLeft",1)

;  plot =  create "lplot" xyPlotClass wks
;   "xyCoordData" : pair(idr,iv,is)
;
;   "trYMaxF" : 90
;   "trYMinF" : -90
;   "trXMaxF" : 4
;   "trXMinF" : 0
;
;   "tmXBMode" : "MANUAL"
;   "tmXBTickSpacingF" : 2
;   "tmXBMinorPerMajor" : 1
;   "tmXBFormat" : "@*+^sg"
;   "tmYLMode" : "Manual"
;   "tmYLTickSpacingF"  : 30
;   "tmYLMinorPerMajor" : 2
;
;   "tiXAxisString" : xstr2(iv)
;   "tiYAxisOn" : False
;  end create
;  set_size(plot,xf2,yf(iv),wi2,he)
;
;  draw(plot)


  end do  ; idr
  end do  ; iv

;  dr_txt(wks,cmon(0),xf(0)+wi/2,yf(0)-4,22,10,0.)
;  dr_txt(wks,cmon(1),xf(1)+wi/2,yf(0)-4,22,10,0.)


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
   set_size_lb(label,35,yf(2)+he+15,120,8)
   setvalues label
     "lbLabelFontHeightF" : 0.012
     "lbTitleFont"        : 21
     "lbTitleFontHeightF" : 0.015
   end setvalues

   draw(label)


  dr_txt(wks,"(a)",xf(0)-8,yf(0)-4,25,12,0.)
  dr_txt(wks,"(b)",xf(0)-8,yf(1)-4,25,12,0.)
  dr_txt(wks,"(c)",xf(0)-8,yf(2)-4,25,12,0.)


   frame(wks)

  end do  ; is


 end



