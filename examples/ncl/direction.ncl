load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"

begin

  xf = (/0.1,0.5,0.1,0.5/)
  yf = (/0.95,0.95,0.50,0.50/)

;  nocal =
;  tropo = (/12,14,22,23,24,30,38,39,45,56,57,58,59,66,67,68,75,78,87,88,91/)


; read data
  datum = asciiread ("data97_",(/97,7/),"float")

  month  = new(97,integer)
  season = new(97,integer)
  angle  = new(97,float)
  pvel   = new(97,float)
  direction = new(97,integer)

  month(:)  = floattoint(datum(:,0)/10000)
  month(:)  = month(:) - month(:)/100*100
  angle(:)  = datum(:,5)
  pvel(:)   = datum(:,1)/datum(:,3) *1000./60.

  direction = 1
  do i=0, 96
    if ( (angle(i).ge. 45.) .and. (angle(i).lt. 135.) ) then
      direction(i) = 2
    end if
    if ( (angle(i).ge. 135.) .and. (angle(i).lt. 225.) ) then
      direction(i) = 3
    end if
    if ( (angle(i).ge. 225.) .and. (angle(i).lt. 315.) ) then
      direction(i) = 4
    end if
  end do

  pi = acos(-1.)
  angle(:) = (90. - angle(:)) /180.*pi

  arrow = new((/2,97/),float)
  arrow(0,:) = pvel(:) * cos(angle(:))
  arrow(1,:) = pvel(:) * sin(angle(:))

; season
  n1 = 0
  n2 = 0
  n3 = 0
  n4 = 0
  season = 4
  do i=0, 96
    if ( (month(i).ge. 3) .and. (month(i).le. 5) ) then
      season(i) = 1
      n1 = n1 + 1
    end if
    if ( (month(i).ge. 6) .and. (month(i).le. 8) ) then
      season(i) = 2
      n2 = n2 + 1
    end if
    if ( (month(i).ge. 9) .and. (month(i).le. 11) ) then
      season(i) = 3
      n3 = n3 + 1
    end if
  end do
  n4 = 97 - n1 - n2 - n3

  print(n1 +0)
  print(n2 +0)
  print(n3 +0)
  print(n4 +0)

  arrowall = new((/2,3,3,4,97/),float)
  noobs    = new((/4,97/),integer)
  colobs   = new((/4,97/),integer)
  direction1 = new(n1,integer)
  direction2 = new(n2,integer)
  direction3 = new(n3,integer)
  direction4 = new(n4,integer)

  n1 = 0
  n2 = 0
  n3 = 0
  n4 = 0
  arrowall = 1.e32
  do i=0, 96
    if ( season(i) .eq. 1 ) then         ; spring
      arrowall(:,1,1,0,n1) = arrow(:,i)
      noobs        (0,n1) = i
      direction1(n1) = direction(i)
      n1 = n1 + 1
    end if
    if ( season(i) .eq. 2 ) then         ; summer
      arrowall(:,1,1,1,n2) = arrow(:,i)
      noobs        (1,n2) = i
      direction2(n2) = direction(i)
      n2 = n2 + 1
    end if
    if ( season(i) .eq. 3 ) then         ; autumn
      arrowall(:,1,1,2,n3) = arrow(:,i)
      noobs        (2,n3) = i
      direction3(n3) = direction(i)
      n3 = n3 + 1
    end if
    if ( season(i) .eq. 4 ) then         ; winter
      arrowall(:,1,1,3,n4) = arrow(:,i)
      noobs        (3,n4) = i
      direction4(n4) = direction(i)
      n4 = n4 + 1
    end if
  end do
  na = (/n1,n2,n3,n4/)

; last position
  col = (/2,210,47/)  ;198

  do j=0, 3
  do i=0, na(j)-1
    colobs(j,i) = col(2)
;    do n=0, 22
;      if ( (noobs(j,i)+1) .eq. nocal(n) ) then
;        colobs(j,i) = col(0)
;      end if
;    end do
;    do n=0, 20
;      if ( (noobs(j,i)+1) .eq. tropo(n) ) then
;        colobs(j,i) = col(1)
;      end if
;    end do
  end do
  end do

; draw -----------------------------------------------------------------
  wks = create "hist" psWorkstationClass defaultapp
      "wkPSFileName" : "dir_season.ps"
      "wkColorMap"   : "rainbow+white+gray"
  end create

; draw circles
  x0 = new(1801,float)
  cir = new((/1801,1801/),float)

  do i=0, 1800
    x0(i) = i/10.-90.
  end do
  y0 = x0

  do j=0, 1800
  do i=0, 1800
    cir(j,i) = (x0(i)^2 + y0(j)^2)^0.5
  end do
  end do

  fieldcir = create "cir" scalarFieldClass defaultapp
    "sfDataArray" : cir
    "sfYArray" : y0
    "sfXArray" : x0
  end create

  do i=0, 3
  contcir = create "contour" contourPlotClass wks
   "cnScalarFieldData" : fieldcir
   "cnFillOn"          :  False
   "cnLevelSelectionMode" : "ExplicitLevels"
   "cnLevels" : (/20,40,60,80/)
   "cnLinesOn" : True
   "cnMonoLineThickness" : True
   "cnLineThicknessF" : 1.0
   "cnLineColor"   : 239
   "cnInfoLabelOn" : False
   "cnLineLabelPlacementMode" : "Constant"
   "cnLineLabelFont" : 21
   "cnLineLabelFontHeightF" : 0.011
   "cnLineLabelsOn" : False

   "tmXBOn" : False
   "tmXTOn" : False
   "tmYLOn" : False
   "tmYROn" : False
   "tmXBBorderOn" : False
   "tmXTBorderOn" : False
   "tmYLBorderOn" : False
   "tmYRBorderOn" : False

   "vpXF" : xf(i)
   "vpYF" : yf(i)
   "vpWidthF" : 0.35
   "vpHeightF" : 0.35

   "trXMinF" : -90
   "trXMaxF" : 90
   "trYMinF" : -90
   "trYMaxF" : 90
  end create

  draw(contcir)

  delete(contcir)
  end do


; draw axis
  do i=0, 3
  axi = create "lplot" xyPlotClass wks
    "trXMinF" : -90
    "trXMaxF" : 90
    "trYMinF" : -90
    "trYMaxF" : 90

    "tmXBMode"   : "Explicit"
    "tmXBValues" : (/-80,-60,-40,-20,20,40,60,80/)
    "tmXBLabels" : (/80,60,40,20,20,40,60,80/)
    "tmXBMajorLengthF" : -0.003
    "tmXBMinorLengthF" : -0.005
    "tmXBLabelFont" : 21
    "tmXBLabelFontHeightF" : 0.010
    "tmXBFormat" : "@*+^sg"

    "tmYLMode" : "Explicit"
    "tmYLValues" : (/-80,-60,-40,-20,20,40,60,80/)
    "tmYLLabels" : (/80,60,40,20,20,40,60,80/)
    "tmYLMinorLengthF" : -0.005
    "tmYLMajorLengthF" : -0.003
    "tmYLLabelFont" : 21
    "tmYLLabelFontHeightF" : 0.010

    "tmXBOn" : False
    "tmXTOn" : False
    "tmYLOn" : True
    "tmYROn" : False
    "tmXBBorderOn" : False
    "tmXTBorderOn" : False
    "tmYLBorderOn" : True
    "tmYRBorderOn" : False

    "vpXF" : xf(i)+0.175
    "vpYF" : yf(i)
    "vpWidthF" : 0.35
    "vpHeightF" : 0.35
  end create
  draw(axi)
  setvalues axi
    "tmXBOn" : True
    "tmYLOn" : False
    "tmXBBorderOn" : True
    "tmYLBorderOn" : False

    "vpXF" : xf(i)
    "vpYF" : yf(i)+0.175
    "vpWidthF" : 0.35
    "vpHeightF" : 0.35
  end setvalues
  draw(axi)

  delete(axi)
  end do


; draw vectors
  zero1 = new(n1,float)
  zero1 = 0.
  zero2 = new(n2,float)
  zero2 = 0.
  zero3 = new(n3,float)
  zero3 = 0.
  zero4 = new(n4,float)
  zero4 = 0.

  xarr = (/-100.,0.,100/)
  yarr = (/-100.,0.,100/)

  do j=0, 3
  do i=0, na(j)-1

  dir = create "vec" vectorFieldClass defaultapp
      "vfDataArray" : arrowall(:,:,:,j,i)
      "vfXArray"    : xarr
      "vfYArray"    : yarr
      "vfMissingUValueV" : 1.e32
      "vfMissingVValueV" : 1.e32
  end create

  vec = create "vector" vectorPlotClass wks
    "vcVectorFieldData" : dir
    "vcPositionMode" : "ArrowTail"
    "vcRefMagnitudeF" :  30
    "vcRefLengthF"    : 0.10
    "vcRefAnnoOn"     : False
    "vcLineArrowThicknessF" : 1.2
    "vcLineArrowColor"      : 1   ;colobs(j,i)
    "vcLineArrowHeadMaxSizeF" : 0.01
    "vcLineArrowHeadMinSizeF" : 0.01

    "trXMinF" : -90
    "trXMaxF" : 90
    "trYMinF" : -90
    "trYMaxF" : 90

    "tmXBOn" : False
    "tmXTOn" : False
    "tmYLOn" : False
    "tmYROn" : False
    "tmXBBorderOn" : False
    "tmXTBorderOn" : False
    "tmYLBorderOn" : False
    "tmYRBorderOn" : False

    "vpXF" : xf(j)
    "vpYF" : yf(j)
    "vpWidthF" : 0.35
    "vpHeightF" : 0.35
  end create
  draw(vec)

  delete(dir)
  delete(vec)

  end do
  end do


; text
  txtstr = (/"Spring","Summer","Autumn","Winter"/)
  do i=0, 3
  tid = create "Text" textItemClass wks
    "txString" : txtstr(i)
    "txPosXF"  : xf(i) + 0.02
    "txPosYF"  : yf(i)
    "txFont"   : 22
    "txFontHeightF" : 0.018
  end create
  draw(tid)
  delete(tid)
  end do

  frame(wks)

; histogram ------------------------------------------

;  wks2 = create "hist" psWorkstationClass defaultapp
;      "wkPSFileName" : "dir_season_histo.ps"
;      "wkColorMap"   : "gsdtol"
;  end create

;  direc1 = new(n1+n2+n3+n4,integer)
;  direc2 = new(n2+n3+n4,integer)
;  direc3 = new(n3+n4,integer)
;  direc4 = new(n4,integer)
;  do i=0, n4-1
;    direc4(i) = direction4(i)
;    direc3(i) = direction4(i)
;    direc2(i) = direction4(i)
;    direc1(i) = direction4(i)
;  end do
;  do i=0, n3-1
;    direc3(n4+i) = direction3(i)
;    direc2(n4+i) = direction3(i)
;    direc1(n4+i) = direction3(i)
;  end do
;  do i=0, n2-1
;    direc2(n4+n3+i) = direction2(i)
;    direc1(n4+n3+i) = direction2(i)
;  end do
;  do i=0, n1-1
;    direc1(n4+n3+n2+i) = direction1(i)
;  end do


;  res          = True
;  res@gsnDraw  = False
;  res@gsnFrame = False

;  res@gsnHistogramDiscreteBinValues = (/1,2,3,4/)
;  res@tiXAxisString = ":F21:propagation direction"
;  res@tiYAxisString = ""
;  res@tmXBMajorLengthF = 0.010
;  res@tmYLMajorLengthF = 0.010
;  res@tmYLMinorLengthF = 0.005
;  res@tmXBLabelFont      = 21
;  res@tmYLLabelFont      = 21
;  res@tmXBLabels         = (/"north","east","south","west"/)


;  res@gsFillColor = 30
;  plot1 = gsn_histogram(wks2,direc1,res)

;  res@gsFillColor = 20
;  plot2 = gsn_histogram(wks2,direc2,res)
;  overlay(plot1,plot2)

;  res@gsFillColor = 12
;  plot3 = gsn_histogram(wks2,direc3,res)
;  overlay(plot1,plot3)

;  res@gsFillColor = 5
;  plot4 = gsn_histogram(wks2,direc4,res)
;  overlay(plot1,plot4)

;  draw(plot1)


;  label = create "lab" labelBarClass wks2
;    "lbOrientation" : "vertical"
;    "lbBoxCount" : 1
;    "lbFillColors" : 30;(/5,12,20,30/)
;    "lbMonoFillPattern" : True
;    "lbLabelsOn" : True
;    "lbLabelStrings" : "Spring";(/"Winter","Autumn","Summer","Spring"/)
;    "lbAutoManage" : False
;    "lbLabelFont" : 21
;    "lbLabelFontHeightF" : 0.014
;    "lbLabelPosition" : "Bottom";,Right"
;    "lbLabelAlignment" : "BoxCenters";"InteriorEdges"
;    "lbLabelAngleF" : 0.
;    "lbPerimOn" : False
;    "lbTitleString" : ""
;    "lbTitlePosition" : "Right"
;    "lbTitleDirection" : "Across"
;    "lbTitleFont" : 21
;    "lbTitleFontHeightF" : 0.015
;    "vpXF" : 0.21
;    "vpYF" : 0.85
;    "vpWidthF" : 0.14
;    "vpHeightF" : 0.03
;  end create
;  draw(label)
;  setvalues label
;    "lbFillColors" : 20
;    "lbLabelStrings" : "Summer"
;    "vpXF" : 0.36
;  end setvalues
;  draw(label)
;  setvalues label
;    "lbFillColors" : 12
;    "lbLabelStrings" : "Autumn"
;    "vpXF" : 0.51
;  end setvalues
;  draw(label)
;  setvalues label
;    "lbFillColors" : 5
;    "lbLabelStrings" : "Winter"
;    "vpXF" : 0.66
;  end setvalues
;  draw(label)



;  frame(wks2)


end

