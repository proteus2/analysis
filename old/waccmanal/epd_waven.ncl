
 begin
 
 f = addfile("/data3/kyh/analy/waccmanal/epf_yp-ctl.nc","r")
 epd0 = f->epd
 lat = f->lat
 lev = f->lev
 month = f->month
 wn = f->wn
 f = addfile("/data3/kyh/analy/waccmanal/epf_yp-lsgwdc.nc","r")
 epd1 = f->epd

 nfig = 10
 ny = 64
 for = new((/nfig,12/),float)

 z0 = 54  ; 85 mb 
 z1 = 54  ; 85 mb 
 print(lev(z0))
 print(lev(z1))
 y0 = 28  ; -9.767145
 y1 = 28  ; -9.767145
 y2 = 35  ;  9.767145 
 y3 = 35  ;  9.767145 
 print(lat(y0))
 print(lat(y1))
 print(lat(y2))
 print(lat(y3))

 do n = 0, 11
   for(0,n) = ( sum(epd0(n,10,z0:z1,y2:y3)) + sum(epd0(n,10,z0:z1,y0:y1)) )/2.
   for(1,n) = ( sum(epd0(n,0 ,z0:z1,y2:y3)) + sum(epd0(n,0 ,z0:z1,y0:y1)) )/2.
   for(2,n) = ( sum(epd0(n,1 ,z0:z1,y2:y3)) + sum(epd0(n,1 ,z0:z1,y0:y1)) )/2.
   for(3,n) = ( sum(epd0(n,2 ,z0:z1,y2:y3)) + sum(epd0(n,2 ,z0:z1,y0:y1)) )/2.
   for(4,n) = ( sum(epd0(n,3 ,z0:z1,y2:y3)) + sum(epd0(n,3, z0:z1,y0:y1)) )/2.
  
   for(5,n) = ( sum(epd1(n,10,z0:z1,y2:y3)) + sum(epd1(n,10,z0:z1,y0:y1)) )/2.
   for(6,n) = ( sum(epd1(n,0 ,z0:z1,y2:y3)) + sum(epd1(n,0 ,z0:z1,y0:y1)) )/2.
   for(7,n) = ( sum(epd1(n,1 ,z0:z1,y2:y3)) + sum(epd1(n,1 ,z0:z1,y0:y1)) )/2.
   for(8,n) = ( sum(epd1(n,2 ,z0:z1,y2:y3)) + sum(epd1(n,2 ,z0:z1,y0:y1)) )/2.
   for(9,n) = ( sum(epd1(n,3 ,z0:z1,y2:y3)) + sum(epd1(n,3 ,z0:z1,y0:y1)) )/2.
 end do
;----------------------------------------------------

 cmap =(/  (/1.000,1.000,1.000/), \
           (/0.000,0.000,0.000/), \
           (/0.600,0.000,0.000/), \
           (/0.000,0.000,0.600/), \
           (/0.600,0.600,0.600/) /)

 col  = (/1,1,1,1,1,3,3,3,3,3/)
 dash = (/0,1,2,3,4,0,1,2,3,4/)
 thick= (/3,1,1,1,1,3,1,1,1,1/)

 vpx = 0.1
 vpy = 0.9
 width = 0.5
 height = 0.50

 wks = create "plot_wmax" psWorkstationClass defaultapp
         "wkPSFileName"   : "epd_waven.eps"
         "wkPSFormat"     : "eps"
         "wkColorMap"     : cmap
 end create

 do i = 0, nfig-1
 dataid0 = create "xyData"coordArraysClass defaultapp
            "caXArray"      : month
            "caYArray"      : for(i,:)
 end create

 plotid = create "xyplot" xyPlotClass wks

  "xyCoordData" : dataid0

  "tmXBMode" : "Explicit"
  "tmXBValues" : (/1,2,3,4,5,6,7,8,9,10,11,12/)
  "tmXBLabels" : (/"Jan.","Feb.","Mar","Apr.","May","June","July","Agu.",\
                   "Sep.","Oct.","Nov.","Dec."/)
  "tmXBLabelFont" : 21
  "tmXBLabelFontHeightF" : 0.013
  "tmXBMajorLengthF" : -0.005
  "tmXBMinorLengthF" : -0.003

  "tmYLMode" : "MANUAL"
  "tmYLTickSpacingF" : 0.25
  "tmYLMinorPerMajor" : 4 
  "tmYLLabelFont" : 21
  "tmYLLabelFontHeightF" : 0.013
  "tmYLMajorLengthF" : -0.005
  "tmYLMinorLengthF" : -0.003

  "trYMinF" : -1.0 
  "trYMaxF" : 0.25
  "trXMinF" : 1 
  "trXMaxF" : 12 

  "tiMainFont" : 22
  "tiMainFontHeightF" : 0.016
  "tiMainString" : "EP flux divergence (10S plus 10N, 85 hPa)"
  "tiMainJust" : "CenterLeft"
  "tiMainPosition" : "Left"
  "tiXAxisFont" : 21
  "tiXAxisFontHeightF" : 0.016
  "tiXAxisString" : "Month"
  "tiYAxisFont" : 21
  "tiYAxisFontHeightF" : 0.016
  "tiYAxisString" : "EPD (m s:S:-1:N: day:S:-1:N:)" 
 
  "vpXF" : vpx
  "vpYF" : vpy
  "vpWidthF" : width 
  "vpHeightF" : height 
  end create

  getvalues plotid
    "xyCoordDataSpec" : dataspec
  end getvalues
  setvalues dataspec
   "xyLineColor" : col(i) 
   "xyDashPattern": dash(i)
   "xyLineThicknessF" : thick(i)
  end setvalues
  draw(plotid)

  getvalues wks
    "wkDefGraphicStyleId" : gsid
  end getvalues
  setvalues gsid
    "gsLineColor"      : 4
    "gsLineThicknessF" : 1
  end setvalues
  NhlDataPolyline(plotid,gsid,(/0,12/),(/0,0/))

  end do

  pid = create "Legend" legendClass wks
    "vpXF"          : vpx+width-0.05
    "vpYF"          : vpy
    "vpWidthF"      : 0.30
    "vpHeightF"     : 0.15
    "lgItemCount" : nfig
    "lgLabelStrings": (/"EPD (Total, CTL)","EPD (W1, CTL)","EPD (W2, CTL)",\
                        "EPD (W3, CTL)","EPD(W4 ,CTL)",\
                        "EPD (Total, RAY)","EPD (W1, RAY)","EPD (W2, RAY)",\
                        "EPD (W3, RAY)","EPD (W4, RAY)"/)
    "lgDashIndexes" : dash 
    "lgLineColors" : col 
    "lgLineLabelsOn" : False
    "lgMonoLineThickness" : False 
    "lgLineThicknesses" : thick 
    "lgBoxMinorExtentF" : 0.2
    "lgLabelFontHeightF" : 0.013
    "lgLabelFont" : 21
    "lgPerimOn" : False
  end create
  draw(pid)

  frame(wks)

  end
