begin

 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

 setvalues NhlGetWorkspaceObjectId()
  "wsMaximumSize" : 55000000
 end setvalues

;-------------------------------------------------------
 wks = create "WPwork"  psWorkstationClass defaultapp
  "wkPSFileName"                :  "map.ps"
  "wkColorMap"                  :  "gsdtol"
 end create

;-------------------------------------------------------
 gg = 0     ; 0 => larger scale,  1 => smaller scale
 err_sig = (/0,0,0,0,0,0,0,0,0,0/)
 wi0 = (/0.6,0.6,0.6/)
 he0 = (/0.6,0.6,0.6/)
 xf = (/0.1,0.1,0.1/)
 yf = (/0.95,0.95,0.95/)
 xmin = (/120,127,128/)
 xmax = (/135,131,130/)
 ymin = (/30,34,35/)
 ymax = (/45,38,37/)
 xtsp = (/5,1,0.5/)
 ytsp = (/5,1,0.5/)
 xtspm = (/4,0,4/)
 ytspm = (/4,0,4/)
 tit  = (/"longitude (deg)"/)

;-------------------------------------------------------
 mapid = create "map" mapPlotClass wks
  "vpXF"                        : xf(gg)
  "vpYF"                        : yf(gg)
  "vpWidthF"                    : wi0
  "vpHeightF"                   : he0
  "mpFillOn"                    : False
  "mpOutlineOn"                 : True
  "mpGridAndLimbOn"             : True ;False
  "mpPerimOn"                   : True
  "mpDataBaseVersion"           : "Ncarg4_1"  ; "RANGS_GSHHS"
  "mpLimitMode"                 : "LatLon"
  "mpMaxLatF"                   : ymax(gg)
  "mpMinLatF"                   : ymin(gg)
  "mpMaxLonF"                   : xmax(gg)
  "mpMinLonF"                   : xmin(gg)
;  "mpGeophysicalLineColor"      : "grey"
 end create
 draw(mapid)

 getvalues mapid
   "vpWidthF"  : wi_fix
   "vpHeightF" : he_fix
   "vpXF"      : xf_fix
   "vpYF"      : yf_fix
 end getvalues
 print(wi0(gg))
 print(he0(gg))
 print(wi_fix)
 print(he_fix)

;-------------------------------------------

 plot1 =  create "plot0" xyPlotClass wks
  "trYMaxF" : ymax(gg)
  "trYMinF" : ymin(gg)
  "trXMaxF" : xmax(gg)
  "trXMinF" : xmin(gg)

  "tmXBMode" : "MANUAL"
  "tmXBTickSpacingF"  : xtsp(gg)
  "tmXBMinorPerMajor" : xtspm(gg)
  "tmXBMajorLengthF" : -0.007
  "tmXBMinorLengthF" : -0.004
  "tmXBLabelFont"        : 21
  "tmXBLabelFontHeightF" : 0.012
  "tmXBAutoPrecision" : False
  "tmXBFormat" : "@*+^sg"

  "tmYLMode" : "MANUAL"
  "tmYLTickSpacingF"  : ytsp(gg)
  "tmYLMinorPerMajor" : ytspm(gg)
  "tmYLMajorLengthF" : -0.007
  "tmYLMinorLengthF" : -0.004
  "tmYLLabelFont"        : 21
  "tmYLLabelFontHeightF" : 0.012
  "tmYLFormat" : "@*+^sg"

  "tiMainFont"        : 21
  "tiMainFontHeightF" : 0.016
  "tiMainString"      : ""
  "tiXAxisFont"        : 21
  "tiXAxisFontHeightF" : 0.016
  "tiXAxisString"      : "longitude (deg)"
  "tiYAxisFont"        : 21
  "tiYAxisFontHeightF" : 0.016
  "tiYAxisString"      : "latitude (deg)"

  "vpYF" : yf_fix
  "vpXF" : xf_fix
  "vpHeightF" : he_fix
  "vpWidthF"  : wi_fix
 end create
 draw(plot1)

;----------------------------------------------------
 mkx = (/128.9/)
 mky = (/36.2/)
; mkx = (/125.1,128.34,124.40,127.03,129.23,126.81,125.27,126.26,126.34/)
; mky = (/32.07,38.15,37.58,37.09,36.02,35.1,34.41,33.17,34.33/)


 mres = True
 mres@gsMarkerIndex  = 5 ;14
 mres@gsMarkerSizeF  = 15.0
 mres@gsMarkerColor  = "black"
 mres@gsMarkerThicknessF = 3
 gsn_polymarker(wks,mapid,mkx,mky,mres)

;---------------------------------------------------
 frame(wks)

end
