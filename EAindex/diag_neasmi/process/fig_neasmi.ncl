load "header.ncl"

xspace = 170
yspace = 240
load "figsize.ncl"
load "attribute.ncl"

axisfontratio = 1.0

begin

;--------------------------------------- parameters
 year1 = stringtoint(getenv("year1"))
 year2 = stringtoint(getenv("year2"))
 yearf = stringtoint(getenv("yearf"))

;--------------------------------------- data reading
 if ( .not. isfilepresent("$fname") ) then
   print("The data file was not downloaded.")
   status_exit(20)
 end if
 fp = addfile("$fname.nc","r")
 ntf = dimsizes( fp->time )
 if ( mod(ntf,12) .ge. 6 .and. mod(ntf,12) .lt. 8 ) then
   ntf = ntf/12*12
 end if
 nt0 = dimsizes( ispan(8,ntf,12) )
 prcp = ( fp->precip(ispan(5,ntf-1,12),{18.75:51.25},{98.75:181.25}) + \
          fp->precip(ispan(6,ntf-1,12),{18.75:51.25},{98.75:181.25}) + \
          fp->precip(ispan(7,ntf-1,12),{18.75:51.25},{98.75:181.25}) )/3.
 time = fspan(1979,1979-1+nt0,nt0)
 lon = fp->lon({98.75:181.25})
 lat = fp->lat({18.75: 51.25})
 prcpm = dim_avg_n(prcp,0)
 prcpa = prcp - conform(prcp,prcpm,(/1,2/))
 prcpa!0 = "time"
 prcpa!1 = "lat"
 prcpa!2 = "lon"
 prcpa&time = time
 prcpa&lat  = lat
 prcpa&lon  = lon
 prcp0 = ( fp->precip(ispan(5,ntf-1,12),{-11.25:71.25},{58.75:181.25}) + \
           fp->precip(ispan(6,ntf-1,12),{-11.25:71.25},{58.75:181.25}) + \
           fp->precip(ispan(7,ntf-1,12),{-11.25:71.25},{58.75:181.25}) )/3.
 lon0 = fp->lon({ 58.75:181.25})
 lat0 = fp->lat({-11.25: 71.25})
 prcpm0 = dim_avg_n(prcp0,0)
 prcpa0 = prcp0 - conform(prcp0,prcpm0,(/1,2/))
 prcpa0!0 = "time"
 prcpa0!1 = "lat"
 prcpa0!2 = "lon"
 prcpa0&time = time
 prcpa0&lat  = lat0
 prcpa0&lon  = lon0

;--------------------------------------- process 0
 opteof      = True
 opteof@jopt = 0
 ev = eofunc_Wrap(prcpa(lat|:,lon|:,time|:),2,opteof)
 print(ev@pcvar)
 ev1 = ev(0,:,:)
 amv = avg( ev1({30.:50.},{110.:145.}) )
 ets = eofunc_ts_Wrap(prcpa(lat|:,lon|:,time|:),ev(0:1,:,:),False)
 print(ets@ts_mean)
 neasmi = ets(0,:)*amv

 evf = regCoef(ets(0,:), prcpa0(lat|:,lon|:,time|:)) * ets(0,yearf-1979)
 prcpa1 = prcpa0(yearf-1979,:,:)

;--------------------------------------- cn level
 cnlev = fspan(2,18,9)
 cnfil = (/0,4,5,6,8,9,10,11,12,13/)

 cnleva = fspan(-6,6,13)*0.75
 cnfila = (/14,15,16,17,18,19,20,20,21,22,23,24,25,26/)

;--------------------------------------- drawing
 wi = (/nt0*2,36/)
 he = (/20,30/)
 xf = (/10,10,50,90/)
 yf = (/10,42/)

 outfile = getenv("out1")
 wks = cr_wks(outfile,600,0)
 gsn_merge_colormaps(wks,"precip_11lev","precip_diff_12lev")
 colind = NhlNewColor(wks,0.6,0.6,0.6)
 plot = xylineplot(wks,time,neasmi,"D",0,"T",1.2,"C",1)
 xyaddmarker(plot,16,"S",0.5,"T",1.0,"C",1)
 xyaddline(plot,time,0.,"D",0,"T",0.6,"C",colind)
 set_size(plot,xf(0),yf(0),wi(0),he(0))
 axis_range(plot,year1-0.5,year2+0.5,-999,999,False,False)
 setvalues plot
  "tmXBMode" : "Manual"
  "tmXBTickStartF"    : 1975
  "tmXBTickSpacingF"  : 5
  "tmXBMinorPerMajor" : 4
 end setvalues
 axis_str(plot,"Northeast Asian Summer Monsoon Index","year","(mm day:S:-1:N:)")
 draw(plot)

 map = create "mp" mapPlotClass wks
  "mpProjection"     : "Mollweide" ;"Robinson"
  "mpCenterLonF"     : 125.
  "mpOutlineOn"      : True
  "mpGridAndLimbOn"  : False
  "mpPerimOn"        : True
  "mpPerimDrawOrder" : "PostDraw"

  "mpShapeMode"     : "FixedAspectFitBB"
  "mpGeophysicalLineColor" : "Orange4"
  "mpDataBaseVersion" : "Ncarg4_0"

  "mpLimitMode"    : "LatLon"
  "mpMinLonF"      : 85.
  "mpMaxLonF"      : 165.
  "mpMinLatF"      : 0.
  "mpMaxLatF"      : 50.
 end create
 set_size_mp(map,xf(1),yf(1),wi(1),he(1))

 contour1 = cnshadeplot(wks,lon0,lat0,prcpm0,cnlev,cnfil,"(mm day:S:-1:N:)")
 setvalues contour1
  "cnMissingValFillColor" : colind
  "lbTitlePosition"    : "Right"
  "lbTitleDirection"   : "Across"
  "lbTitleOffsetF"     : 0.10
 end setvalues
 set_size(contour1,xf(1),yf(1),wi(1),he(1))
 axis_str(contour1,"CMAP","","")

 overlay(map,contour1)

 draw(map)

 NhlRemoveOverlay(map,contour1,False)

 set_size_mp(map,xf(2),yf(1),wi(1),he(1))

; contour2 = cnshadeplot(wks,lon0,lat0,prcpa1,cnleva,cnfila,"(mm day:S:-1:N:)")
 contour2 = cnshadeplot(wks,lon0,lat0,evf   ,cnleva,cnfila,"(mm day:S:-1:N:)")
 setvalues contour2
  "cnMissingValFillColor" : colind
  "lbTitlePosition"    : "Right"
  "lbTitleDirection"   : "Across"
  "lbTitleOffsetF"     : 0.10
 end setvalues
 set_size(contour2,xf(2),yf(1),wi(1),he(1))
; axis_str(contour2,"Rainfall anomaly in JJA "+yearf,"","")
 axis_str(contour2,"CMA "+yearf,"","")

 overlay(map,contour2)

 draw(map)

 frame(wks)

 status_exit(10)

end

