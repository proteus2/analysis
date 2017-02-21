load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"

begin

  datum = asciiread ("data97_",(/97,7/),"float")

  horil  = new(97,float)
  period = new(97,float)
  angle  = new(97,float)
  pvel   = new(97,float)
  direction = new(97,integer)

  horil(:)  = datum(:,1)
  period(:) = datum(:,3)
  angle(:)  = datum(:,5)
  pvel(:)   = horil(:)/period(:) *1000./60.

  avehoril  = 0.
  aveperiod = 0.
  avepvel   = 0.
  do i=0, 96
    avehoril  = avehoril  + horil(i)/97.
    aveperiod = aveperiod + period(i)/97.
    avepvel   = avepvel   + pvel(i)/97.
  end do

  do i=0, 96
    if ( (angle(i).ge. 315.) .or. (angle(i).lt. 45.) ) then
      direction(i) = 1
    end if
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


  wks = create "hist" psWorkstationClass defaultapp
      "wkPSFileName" : "histo.ps"
      "wkColorMap"   : "gsdtol"
  end create


;-----------------------------------------------------
  tid = create "Text" textItemClass wks
    "txString" : floattoint(avehoril*100)/100.+" km"
    "txPosXF" : 0.70
    "txPosYF" : 0.76
    "txFont" : 21
    "txFontHeightF" : 0.025
  end create
  draw(tid)

  res = True
  res@gsnHistogramClassIntervals = (/0,10,20,30,40,50,60/)
  res@gsFillColor = 15
  res@tiXAxisString = ":F21:horizontal wavelength (km)"
  res@tiYAxisString = ""
  res@tmXBMajorLengthF = 0.010
  res@tmYLMajorLengthF = 0.010
  res@tmYLMinorLengthF = 0.005
  res@tmXBLabelFont      = 21
  res@tmYLLabelFont      = 21
  res@trYMaxF = 40

  plot = gsn_histogram(wks,horil,res)

;-----------------------------------------------------
  tid = create "Text" textItemClass wks
    "txString" : floattoint(aveperiod*100)/100.+" min"
    "txPosXF" : 0.70 
    "txPosYF" : 0.76 
    "txFont" : 21 
    "txFontHeightF" : 0.025 
  end create 
  draw(tid) 

  res2 = True
  res2@gsnHistogramClassIntervals = (/0,5,10,15,20,25,30,35,40,45,50/)
  res2@gsFillColor = 15
  res2@tiXAxisString = ":F21:period (min)"
  res2@tiYAxisString = ""
  res2@tmXBMajorLengthF = 0.010
  res2@tmYLMajorLengthF = 0.010
  res2@tmYLMinorLengthF = 0.005
  res2@tmXBLabelFont      = 21
  res2@tmYLLabelFont      = 21
  res2@trYMaxF = 40
 
  plot2 = gsn_histogram(wks,period,res2)

;-----------------------------------------------------
  tid = create "Text" textItemClass wks
    "txString" : floattoint(avepvel*100)/100.+" m s:S:-1:N:"
    "txPosXF" : 0.70
    "txPosYF" : 0.76 
    "txFont" : 21 
    "txFontHeightF" : 0.025 
  end create 
  draw(tid) 

  res3 = True
  res3@gsnHistogramClassIntervals = (/0,10,20,30,40,50,60,70,80,90,100/)
  res3@gsFillColor = 15
  res3@tiXAxisString = ":F21:phase speed (ms:S:-1:N:)"
  res3@tiYAxisString = ""
  res3@tmXBMajorLengthF = 0.010
  res3@tmYLMajorLengthF = 0.010
  res3@tmYLMinorLengthF = 0.005
  res3@tmXBLabelFont      = 21
  res3@tmYLLabelFont      = 21
  res3@tmYLMinorPerMajor = 4
  res3@trYMaxF = 30

  plot3 = gsn_histogram(wks,pvel,res3)

;-----------------------------------------------------

  res4 = True
  res4@gsnHistogramDiscreteBinValues = (/1,2,3,4/)
  res4@gsFillColor = 15
  res4@tiXAxisString = ":F21:propagation direction"
  res4@tiYAxisString = ""
  res4@tmXBMajorLengthF = 0.010
  res4@tmYLMajorLengthF = 0.010
  res4@tmYLMinorLengthF = 0.005
  res4@tmXBLabelFont      = 21
  res4@tmYLLabelFont      = 21
  res4@tmXBLabels         = (/"north","east","south","west"/)
  res4@trYMaxF = 35

  plot4 = gsn_histogram(wks,direction,res4)




end

