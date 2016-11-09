load "$UTIL/ncl/header.ncl"

xspace = 170  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 0.85

begin

;--------------------------------------- parameters
 expname = (/"merra","jra55","cfsr","era-int","era-int_f"/)
 datname0= (/"MERRA","JRA55","CFSR","ERA-I"  ,"ERA-I-nr"/)
 datname = (/"MERRA","JRA55","CFSR","ERA-I"  ,"ERA-I_ml"/)
 uname = (/"u","U_GRD_GDS0_ISBL","UGRD","U_GDS0_ISBL","u"/)

 ie = 3
 y0 = 1979
 y1 = 1994;+7
 y2 = 1999;+7
 lat0 = 5

 wks = cr_wks("fig/u_tz-"+expname(ie),600,0)
 gsn_define_colormap(wks,"StepSeq25")
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- read data
 nt = (y2-y1+1)*12
 n0 = (y1-y0)*12

 fnam0 = "$DATD/clim/"+datname0(ie)+"/"+expname(ie)+".u_ypt.nc"

 latrng = (/-lat0,lat0/)

   f0 = addfile(fnam0,"r")

   latname = f0->$uname(ie)$!2
   lat = f0->$latname$({-lat0:lat0})
   dy05 = 0.5*(lat(1) - lat(0))
   wgt = where( lat .ge. latrng(0)+dy05 .and. lat .le. latrng(1)-dy05, 1., 0.)
   wgt = where( abs(lat-latrng(0)) .lt. dy05, 0.5*(1.+(lat-latrng(0))/dy05), wgt)
   wgt = where( abs(lat-latrng(1)) .lt. dy05, 0.5*(1.-(lat-latrng(1))/dy05), wgt)

   print(expname(ie)+"  "+abs(lat(0))+"  "+wgt(0))

   zname = f0->$uname(ie)$!1
   if (expname(ie) .eq. "era-int_f") then
     u0 = transpose( dim_avg_wgt_Wrap( tofloat( f0->$uname(ie)$(n0:n0+nt-1,{:25},{-lat0:lat0}) ), wgt,1) )
     pa = asciiread("/data11/data-arch/ERA-I-nr/plev/plev_ps1000_ERA-I",60,"float")
     zu = 18.6 + 6.6*log( 70./ pa(2:24) )
   end if
   if ( expname(ie) .eq. "era-int" .or. expname(ie) .eq. "jra55" .or. \
        expname(ie) .eq. "cfsr" ) then
     u0 = transpose( dim_avg_wgt_Wrap( tofloat( f0->$uname(ie)$(n0:n0+nt-1,{:100},{-lat0:lat0}) ), wgt,1) )
     zu = 18.6 + 6.6*log( 70./ f0->$zname$({:100}) )
   end if
   if (expname(ie) .eq. "merra") then
     u0 = transpose( dim_avg_wgt_Wrap( tofloat( f0->$uname(ie)$(n0:n0+nt-1,{100:},{-lat0:lat0}) ), wgt,1) )
     zu = 18.6 + 6.6*log( 70./ f0->$zname$({100:}) )
   end if

   t = y1 + (fspan(1,nt,nt)-0.5)/12.

;--------------------------------------- process 0

;--------------------------------------- cn level
 cnlev = ispan(-10,10,1)*8
; cnfil = (/2,3,4,5,6,7,8,9,10,11,0,0,26,25,24,23,22,21,20,19,18,17/)
 cnfil = (/17,18,19,20,21,22,23,24,25,26,0,0,11,10,9,8,7,6,5,4,3,2/)

;--------------------------------------- contour
 wi = 40
 he = 15
 xf = 22
 yf = 10

 mstr = datname(ie)
 xstr = "year"
 ystr = "height (km)"
 lbopt = "(m s:S:-1:N:)"
 txtstr = "U"

; contour1 = cnshadeplot(wks,t,zu,u0,cnlev,cnfil,lbopt)
 contour1 = cnlineplot(wks,t,zu,u0,cnlev,0)
 cnline_att(contour1,"D",(/0,0,1/),"T",0.8,"C",1)
 setvalues contour1
  "pmLabelBarOrthogonalPosF" : 0.12 ;+ 0.2
  "pmLabelBarParallelPosF"   : 0.50 ;+ 0.55
  "pmLabelBarKeepAspect"     : False
  "pmLabelBarWidthF"         : 0.60
  "pmLabelBarHeightF"        : 0.20
  "lbLabelStride"      : 2
  "lbLabelFontHeightF" : 0.016
  "lbTitlePosition"    : "Right"
  "lbTitleDirection"   : "Across"
  "lbTitleOffsetF"     : 0.05 + 0.05
  "lbTitleFontHeightF" : 0.018
 end setvalues

 set_size(contour1,xf,yf,wi,he)
 axis_range(contour1,y1,y2+1,19,33,False,False)
 axis_str(contour1,mstr,xstr,ystr)
 setvalues contour1
  "tmXBMode" : "Manual"
  "tmXBTickSpacingF"  : 1
  "tmXBMinorPerMajor" : 0
  "tiMainFontHeightF" : sz_mstr*axisfontratio*1.25
 end setvalues

 cnfinalize(contour1,wks)

 draw(contour1)

 dr_txt2(wks,txtstr,xf+1,yf+he-2,22,6,0.,"CenterLeft",1)


 frame(wks)

end

