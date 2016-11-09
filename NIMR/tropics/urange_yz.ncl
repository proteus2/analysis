load "$UTIL/ncl/header.ncl"

xspace =  83  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 1.2

begin

;--------------------------------------- parameters
 expname = "uj"
 ppid    = "j" 
 varname = "u"

 y0 = 1953
 y1 = 1953
 y2 = 2006
 lat0 = 5.

 wks = cr_wks("fig/urange_yz",600,0)
 gsn_define_colormap(wks,"StepSeq25")
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- read data
 nt = (y2-y1+1)*12
 n1 = (y1-y0)*12

 fname = "$DATD/L60CGW-t/xm.uan"+expname+"a.p"+ppid+".nc"
 f0 = addfile(fname,"r")
 var = f0->$varname$(n1:n1+nt-1,:,:)
 zname = var!1
 z = var&$zname$/1.e3
 yname = var!2
 y = var&$yname$

 fname = "$DATD/L60CGW-t/n2_t.nc"
 f0 = addfile(fname,"r")
 nb = f0->N2(n1:n1+nt-1,:)
 nb = sqrt(nb)
 nzn = dimsizes(f0->hybrid_ht)
 c1 = nb*2.4e3/(2.*3.141592)
; var(:,15:,:) = var(:,15:,:) + conform_dims((/nt,nzn,dimsizes(y)/),c1,(/0,1/))

;--------------------------------------- process 0
 vmean = dim_avg_n( var, 0 )
 vsd = dim_stddev_n( var, 0 )

;--------------------------------------- cn level
 cnlev = ispan(-38,38,4)
 cnlev2 = ispan(1,19,2)
 cnlev = ispan(-19,19,2)
 cnlev2 = ispan(1,10,1)
; cnfil = (/2,3,4,5,6,7,8,9,10,11,0,0,26,25,24,23,22,21,20,19,18,17/)
 cnfil = (/17,18,19,20,21,22,23,24,25,26,0,11,10,9,8,7,6,5,4,3,2/)

;--------------------------------------- contour
 wi = 50
 he = 40
 xf = 10 + (/0,0,0,0,0,0,0,0/)*78
 yf = 10 + (/1,2,3,4,5,6,0,7/)*32

 xstr = "lat (deg)"
 ystr = "height (km)"
 lbopt = "(m s:S:-1:N:)"

 contour1 = cnshadeplot(wks,y,z,vmean,cnlev,cnfil,lbopt)
 cnaddline(contour1,wks,y,z,vsd,cnlev2,1,"D",0,"T",1,"C",1)
 setvalues contour1
  "pmLabelBarKeepAspect"     : False
  "pmLabelBarWidthF"         : 0.7
  "pmLabelBarHeightF"        : 0.15
  "lbLabelAutoStride"  : False
  "lbLabelFontHeightF" : 0.015
  "lbTitleDirection"   : "Across"
  "lbTitleOffsetF"     : 0.02
  "lbTitleFontHeightF" : 0.018
 end setvalues
 set_size(contour1,xf,yf,wi,he)
 axis_range(contour1,-999,999,-999,30,False,False)
 axis_str(contour1,"",xstr,ystr)

 cnfinalize(contour1,wks)

 draw(contour1)


 frame(wks)

end



