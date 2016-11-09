load "$UTIL/ncl/header.ncl"

xspace = 83   ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 1.0

begin

;--------------------------------------- parameters
 expname = (/"uj"/)

 lat0 = 25.
 lat0avg = 10

 wks = cr_wks("fig/dt_my_m",600,0)
 gsn_define_colormap(wks,"precip3_16lev")
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- read data
   fnam1 = "$DATD/L60CGW-t/xm.uan"+expname(0)+"a.pe.nc"

   f1 = addfile(fnam1,"r")
   lat1 = f1->latitude({-lat0:lat0})

   pra = f1->dT_lsrain(:,{10.e3},{-lat0:lat0})*72.
;   prd = f1->precip_dc(:,{10.e3},{-lat0:lat0})*72.
;   prs = f1->precip_sc(:,{10.e3},{-lat0:lat0})*72.
;   prm = f1->precip_mc(:,{10.e3},{-lat0:lat0})*72.
;   prc = prd + prs + prm

   pr = pra ;- prc

   ny = dimsizes(lat1)
   nt = dimsizes(f1->t)

   t = fspan(0,13,14)

   uss = new((/ny,14/),float)
   do n=1, 12
     uss(:,n) = dim_avg_n(pr(ispan(n-1,nt-1,12),:),0)
   end do
   uss(:,0 ) = uss(:,12)
   uss(:,13) = uss(:,1 )
   uss!0 = "lat"
   uss&lat = lat1

;--------------------------------------- process 0

;--------------------------------------- cn level
 cnlev = (/1,2,3,4,5,6,7,8,9,10,11,12,13,14,15/)*0.5
 cnfil = (/3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18/)

;--------------------------------------- contour
 wi = 30
 he = 22
 xf = 12 + (/0,1/)*42
 yf = 10

 mstr = "Heating"
 xstr = "month"
 ystr = (/"lat (deg)","heating (K day:S:-1:N:)"/)


 contour1 = cnshadeplot(wks,t,lat1,abs(uss),999,cnfil,"(K day:S:-1:N:)")
 setvalues contour1
;  "cnFillMode" : "RasterFill"
  "pmLabelBarOrthogonalPosF" : 0.10
;  "pmLabelBarParallelPosF"   : 1.00
  "pmLabelBarKeepAspect"     : False
  "pmLabelBarWidthF"         : 0.6
  "pmLabelBarHeightF"        : 0.12
  "lbLabelAutoStride" : True
;  "lbLabelFontHeightF" : 0.03
  "lbTitlePosition"    : "Right"
  "lbTitleDirection"   : "Across"
  "lbTitleOffsetF"     : 0.02
;  "lbTitleFontHeightF" : 0.04
 end setvalues

 set_size(contour1,xf(0),yf,wi,he)
 axis_range(contour1,0.5,12.5,-lat0,lat0,False,False)
 axis_str(contour1,mstr,xstr,ystr(0))
 setvalues contour1
  "tmXBMode" : "Explicit" ;"Manual"
  "tmXBValues" : ispan(1,12,1)
  "tmXBLabels" : (/"J","F","M","A","M","J","J","A","S","O","N","D"/)
 end setvalues

 cnfinalize(contour1,wks)


 draw(contour1)

 plot = xylineplot(wks,t,dim_avg_n(uss({-lat0avg:lat0avg},:),0),"D",0,"T",1.5,"C",1)
 set_size(plot,xf(1),yf,wi,he)
 axis_range(plot,0.5,12.5,-999,999,False,False)
 axis_str(plot,"","",ystr(1))
 setvalues plot
  "tmXBMode" : "Explicit" ;"Manual"
  "tmXBValues" : ispan(1,12,1)
  "tmXBLabels" : (/"J","F","M","A","M","J","J","A","S","O","N","D"/)
 end setvalues
 draw(plot)

  dr_txt(wks,lat0avg+"N:S:o:NF29:-:F22:"+lat0avg+":S:o:N:S",xf(1)+wi-8,yf+2,22,6,0.)


 frame(wks)

end

