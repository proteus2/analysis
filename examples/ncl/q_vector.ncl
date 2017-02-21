begin   ; synoptic chart

 setvalues NhlGetWorkspaceObjectId()
   "wsMaximumSize" : 55000000
 end setvalues

 load "/usr/local/ncarg/lib/ncarg/nclscripts/csm/gsn_code.ncl"

 wks = create "contour" psWorkstationClass defaultapp
     "wkPSFileName" : "q8501.ps"
     "wkColorMap" : "BlWhRe"
 end create

 inixf = 0.03
 iniyf = 1.00

 xf = (/0.,0.,0.,0.,0.5,0.5,0.5,0.5/) + inixf
 yf = iniyf - (/0.,0.25,0.5,0.75,0.,0.25,0.5,0.75/)
;-------------------------------------------------------------
 file1 = addfile("nc/qx.nc","r")
 x = file1->lon
 y = file1->lat
 t = file1->time
 qx = file1->Qx

 file2 = addfile("nc/qy.nc","r")
 qy = file2->Qy

 file3 = addfile("nc/div_q.nc","r")
 div = file3->div_Q

 q_vec = new((/2,dimsizes(t),dimsizes(y),dimsizes(x)/),float)
 q_div = new((/dimsizes(t),dimsizes(y),dimsizes(x)/),float)

 mag = (/3,3,2,1/)*1.E-11
 tn  = (/0,7,6/)
 ti  = (/0,0,8/)

 lev = 3     ; opt.
 tim = 1
 
 q_vec(0,:,:,:) = qx(:,lev,:,:)
 q_vec(1,:,:,:) = qy(:,lev,:,:)
 q_div(:,:,:) = div(:,lev,:,:)

;-------------------------------------------------------------
do i=0, tn(tim)
;-------------------------------------------------------------
 field1 = create "f1" vectorFieldClass defaultapp
  "vfDataArray" : q_vec(:,i+ti(tim),:,:)
  "vfXArray" : x
  "vfYArray" : y
 end create

 field3 = create "f1" scalarFieldClass defaultapp
  "sfDataArray" : q_div(i+ti(tim),:,:)
  "sfXArray" : x
  "sfYArray" : y
 end create

;-------------------------------------------------------------
 map = create "mp" mapPlotClass wks
  "vpXF"      : xf(i)
  "vpYF"      : yf(i)
  "vpWidthF"  : 0.45
  "vpHeightF" : 0.24

  "mpProjection"    : "LambertConformal"
  "mpLambertMeridianF" : 115
  "mpFillOn"        : False
  "mpLabelsOn"      : True;False
  "mpOutlineOn"     : True
  "mpGridAndLimbOn" : False
  "mpGridLineColor" : "gray"
  "mpGridLineDashPattern" : 1
  "mpGridLineDashSegLenF" : 0.2
  "mpPerimOn"       : True

  "mpShapeMode"     : "FixedAspectFitBB"
  "mpGeophysicalLineColor" : "Orange4"
  "mpDataBaseVersion" : "Ncarg4_0"

  "mpGridLonSpacingF" : 10
  "mpGridLatSpacingF" : 10
  "mpLimitMode"    : "LatLon"
  "mpMinLonF"      : 80.
  "mpMaxLonF"      : 150.
  "mpMinLatF"      : 20.
  "mpMaxLatF"      : 60.
 end create
 getvalues map
  "vpWidthF"  : wfg
  "vpHeightF" : hfg
  "vpXF"      : xfg
  "vpYF"      : yfg
 end getvalues

;------------------------------------------------------------

 vec1 = create "vector" vectorPlotClass wks
  "vcVectorFieldData" : field1
  "vcRefMagnitudeF" : mag(lev)
  "vcRefLengthF"    : 0.1
  "vcRefAnnoOn"     : False

  "tmXBOn" : False
  "tmYLOn" : False
  "tmXTOn" : False
  "tmYROn" : False

  "vpXF" : xfg
  "vpYF" : yfg
  "vpWidthF" : wfg
  "vpHeightF" : hfg

  "trXMinF" : 55
  "trXMaxF" : 175
  "trYMinF" : 10
  "trYMaxF" : 70

 end create

 cont3 = create "contour" contourPlotClass wks
  "cnScalarFieldData" : field3
  "cnLevelSelectionMode" : "ExplicitLevels"
  "cnLevels"       : (/-15,-10,-5,5,10,15/)*1.E-18
  "cnFillOn"       :  True
  "cnFillColors"   : (/-18,-12,-5,0,5,12,18/)+52
  "cnLinesOn"      : False
  "cnInfoLabelOn"  : False
  "cnLineLabelsOn" : False

  "cnSmoothingDistanceF"  : 0.001
;  "cnSmoothingTensionF"   : 0.1
  "cnSmoothingOn" : True

  "tmXBOn" : False
  "tmYLOn" : False
  "tmXTOn" : False
  "tmYROn" : False

  "vpXF" : xfg
  "vpYF" : yfg
  "vpWidthF" : wfg
  "vpHeightF" : hfg

  "trXMinF" : 55
  "trXMaxF" : 175
  "trYMinF" : 10
  "trYMaxF" : 70
 end create

 overlay(map,cont3) 
 draw(map)
 NhlRemoveOverlay(map,cont3,False)
 NhlAddOverlay(map,vec1,map)
 draw(map)


 delete(field1)
 delete(field3)
 delete(map)
 delete(vec1)
 delete(cont3)
;--------------------------------------------------------------
end do
;--------------------------------------------------------------
   label = create "lab" labelBarClass wks
     "lbOrientation" : "horizontal"
     "lbBoxCount" : 7
     "lbFillColors" : (/-18,-12,-5,0,5,12,18/)+52
     "lbMonoFillPattern" : True
     "lbLabelsOn" : True
     "lbLabelStrings" : (/-15,-10,-5,5,10,15/)
     "lbAutoManage" : False
     "lbLabelFont" : 21
     "lbLabelFontHeightF" : 0.010
     "lbLabelPosition" : "Bottom";,Right"
     "lbLabelAlignment" : "InteriorEdges"
     "lbLabelAngleF" : 0.
     "lbPerimOn" : False
     "lbTitleString" : "  (10:S:-18:N: kg m:S:-5:N: s:S:-1:N:)"
     "lbTitlePosition" : "Bottom";,"Top"
     "lbTitleDirection" : "Across"
     "lbTitleFont" : 21
     "lbTitleFontHeightF" : 0.012
     "vpXF" : 0.60
     "vpYF" : 0.15
     "vpWidthF" : 0.30
     "vpHeightF" : 0.05
   end create
if (tim .eq. 2) then
  draw(label)
end if


 frame(wks)

end
