begin

;------------------------

 setvalues NhlGetWorkspaceObjectId()
   "wsMaximumSize" : 55000000
 end setvalues
;---------------------------------------------- wks

 load "/usr/local/ncarg/lib/ncarg/nclscripts/csm/gsn_code.ncl"

 wks = create "contour" psWorkstationClass defaultapp
     "wkPSFileName" : "fig6.ps"
     "wkColorMap" : "hotres"
     "wkDeviceLowerY" :  80
     "wkDeviceUpperY" : 770
 end create
 asp0 = 1.3125
 rf = (770-80) / 540.
 rx = asp0 / rf
 ry = asp0 * rf


 zk = (/1,2,4,5/)
 xf = (/0.1,0.51,0.1,0.51/)
 yf = (/0.95,0.95,0.58,0.58/)
 xst = (/"","","k (10:S:-4:N: cycle m:S:-1:N:)", \
               "k (10:S:-4:N: cycle m:S:-1:N:)"/)
 yst = (/":F33:w:F22: (10:S:-3:N: cycle s:S:-1:N:)","", \
         ":F33:w:F22: (10:S:-3:N: cycle s:S:-1:N:)",""/)

 nkkall = 601
 nooall = 421
 nkk = 241
 noo = 181

 f0 = addfile("/usr/users/kyh/result4/2ndgen/105jul/UNbar/unbar.nc","r")
 ub = f0->ubar(100:,0)    ; NETCDF bottom
ub = 0.
 fnam1 = "skoz.nc"
 f1 = addfile(fnam1,"r")
 d1 = f1->PSD(:,(nooall-noo)/2:(nooall-noo)/2+noo-1,:nkk-1)
 k1 = f1->k(:nkk-1)
 o1 = f1->o((nooall-noo)/2:(nooall-noo)/2+noo-1)
 z1 = f1->z

 d1 = log10(d1)  ;(d1+1.e-12)

 oh = new((/4,noo,nkk/),float)
 do n=0, noo-1
 do i=0, nkk-1
   oh(:,n,i) = o1(n) - ub(zk(:)+1) * k1(i)         ; intrinsic f
 end do
 end do

;---------------------------------- cn level

 cnlev = new(7,float)
 cnfil = new(8,integer)

 do i=0,6
;  cnlev(i) = 0.5 * i + 3.5
  cnlev(i) = 8.1 - (6-i)*0.4
 end do
 cnfil = (/1,4,8,20,40,70,95,120/)*2

;--------------------------------------------- field

  field = new(4,graphic)
  fieldu = new(4,graphic)
 do i=0, 3
  field(i) = create "w" scalarFieldClass defaultapp
    "sfDataArray" : d1(zk(i),:,:)
    "sfYArray" : o1*1000.
    "sfXArray" : k1*10000.
    "sfMissingValueV" : 1.e32
  end create
  fieldu(i) = create "w" scalarFieldClass defaultapp
    "sfDataArray" : oh(i,:,:)
    "sfYArray" : o1*1000.
    "sfXArray" : k1*10000.
    "sfMissingValueV" : 1.e32
  end create
 end do

;--------------------------------------------- contour
; wprt1  
;----------------------------------------------------
  
  cont1 = create "contour" contourPlotClass wks

   "cnScalarFieldData" : field(0)
   "cnFillOn"        :  True           ; color fill
   "cnLevelSelectionMode" : "ExplicitLevels"
   "cnLevels" : cnlev 
   "cnFillColors" : cnfil
   "cnLinesOn" : False 
   "cnInfoLabelOn" : False
   "cnLineLabelPlacementMode" : "Constant"
   "cnLineLabelFont" : 21
   "cnLineLabelFontHeightF" : 0.011
   "cnLineLabelsOn" : True

   "tmXBMode" : "Manual"
   "tmXBTickSpacingF" : 0.5
   "tmXBMinorPerMajor" : 4
   "tmXTMinorPerMajor" : 3
   "tmXBMajorLengthF" : -0.006
   "tmXBMinorLengthF" : -0.003
   "tmXBLabelFont" : 21
   "tmXBLabelFontHeightF" : 0.010
   "tmXBLabelFontAspectF" : rx

;  "tmXBMode" : "EXPLICIT"
;  "tmXBValues" : topno1
;  "tmXBLabels" : toplength1
;  "tmXBMinorValues" : minorn1
;  "tmXBLabelFont" : 21
;  "tmXBLabelFontHeightF" : 0.011
;  "tmXBMajorLengthF" : -0.007
;  "tmXBMinorLengthF" : -0.004


   "tmYLMode" : "Manual"
   "tmYLTickSpacingF" : 1
   "tmYLMinorPerMajor" : 3
   "tmYLMinorLengthF" : -0.003
   "tmYLMajorLengthF" : -0.007
   "tmYLLabelFont" : 21
   "tmYLLabelFontHeightF" : 0.010
   "tmYLLabelFontAspectF" : rx
   "tmYLFormat" : "@*+^sg"

   "vpXF" : xf(0)
   "vpYF" : yf(0)
   "vpWidthF" : 0.35
   "vpHeightF" : 0.28

   "trXMinF" : k1(0)*100000 
   "trXMaxF" : 1.5
   "trYMinF" : -3 
   "trYMaxF" : 3

   "tiMainString" : "z = "+floattoint(z1(zk(0)))+" km"
   "tiMainFont"  : 22
   "tiMainFontHeightF" : 0.016
   "tiMainFontAspectF" : rx
   "tiXAxisString" : xst(0)
   "tiXAxisFont"  : 22
   "tiXAxisFontHeightF" : 0.013
   "tiXAxisFontAspectF" : rx
   "tiYAxisString" : yst(0)
   "tiYAxisFont"  : 22
   "tiYAxisFontHeightF" : 0.013*rf
   "tiYAxisFontAspectF" : ry
  end create
  draw(cont1)

  setvalues cont1
    "cnScalarFieldData" : fieldu(0)
    "cnFillOn"          :  False          ; color fill
    "cnLevels"          : (/0.0/)
    "cnLinesOn"         : True
    "cnLineColor"       : 0
    "cnLineLabelsOn"    : False
    "cnMonoLineThickness" : True
    "cnLineThicknessF"  : 2.5
  end setvalues
  draw(cont1)

 do i=1, 3
  setvalues cont1
    "cnScalarFieldData" : field(i)
    "cnFillOn"          : True           ; color fill
    "cnLevels"          : cnlev
    "cnLinesOn"         : False
    "vpXF" : xf(i)
    "vpYF" : yf(i)
    "tiMainString" : "z = "+floattoint(z1(zk(i)))+" km"
    "tiXAxisString" : xst(i)
    "tiYAxisString" : yst(i)
  end setvalues
  draw(cont1)
  setvalues cont1
    "cnScalarFieldData" : fieldu(i)
    "cnFillOn"          :  False          ; color fill
    "cnLevels"          : (/0.0/)
    "cnLinesOn"         : True
    "cnLineColor"       : 0
    "cnLineLabelsOn"    : False
    "cnLineThicknessF"  : 2.5
  end setvalues
  draw(cont1)
 end do


   labels = cnlev
   label = create "lab" labelBarClass wks
     "lbOrientation" : "horizontal"
     "lbBoxCount" : 8
     "lbFillColors" : cnfil
     "lbMonoFillPattern" : True
     "lbLabelsOn" : True
     "lbLabelStrings" : labels
     "lbAutoManage" : False
     "lbLabelFont" : 21
     "lbLabelFontHeightF" : 0.011
     "lbLabelFontAspectF" : rx
     "lbLabelPosition" : "Bottom";,Right"
     "lbLabelAlignment" : "InteriorEdges"
     "lbLabelAngleF" : 0.
     "lbPerimOn" : False
     "lbTitleString" : "Log:B:10:N:(PSD)  "
     "lbTitlePosition" : "Left"
     "lbTitleDirection" : "Across"
     "lbTitleFont" : 21
     "lbTitleFontHeightF" : 0.012
     "lbTitleFontAspectF" : rx
     "vpXF" : 0.07
     "vpYF" : 0.21
     "vpWidthF" : 0.41
     "vpHeightF" : 0.035
   end create
   draw(label)

 txtid2 = create "Text" textItemClass wks
   "txPosXF"               :  0.71
   "txPosYF"               :  0.193
   "txFontHeightF"         :  0.012
   "txFontAspectF"         :  rx
   "txFont"                :  21
   "txString"              :  "PSD [(m:S:2:N: s:S:-4:N:) / (cycle m:S:-1:N: cycle s:S:-1:N:)]"
 end create
 draw(txtid2)

 txtid1 = create "Text" textItemClass wks
   "txPosXF"               :  xf(0) + 0.0
   "txPosYF"               :  yf(0) + 0.03
   "txFontHeightF"         :  0.020
   "txFontAspectF"         :  rx
   "txFont"                :  25
   "txString"              :  "(a)"
 end create
 draw(txtid1)
 setvalues txtid1
   "txPosXF"               :  xf(1) + 0.0
   "txString"              :  "(b)"
 end setvalues
 draw(txtid1)
 setvalues txtid1
   "txPosYF"               :  yf(2) + 0.03
   "txString"              :  "(d)"
 end setvalues
 draw(txtid1)
 setvalues txtid1
   "txPosXF"               :  xf(2) + 0.0
   "txString"              :  "(c)"
 end setvalues
 draw(txtid1)


  frame(wks)

 end



