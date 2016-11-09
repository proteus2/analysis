load "$UTIL/ncl/header.ncl"

xspace =  83  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 1.2

begin

;--------------------------------------- parameters
 expname = (/"era-int","merra"/)
 dirname = (/"ERA-inter","MERRA"/)
 expname = (/"merra","merra"/)
 dirname = (/"MERRA","MERRA"/)

 tstt = 1979
 tend = 2006

 y1 = tstt ;+10
 y2 = y1 + 5

 wks = cr_wks("fig/t_yp-ra",600,0)
 gsn_define_colormap(wks,"StepSeq25")
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- read data
 ne = dimsizes(expname)
 n0 = (y1-tstt)*12
 n9 = (y2-tstt+1)*12 - 1

 fdir  = "../dat/clim/"+dirname
 fname = fdir+"/"+expname+".t_ypt."+tstt+"-"+tend+".nc"
 f1 = addfiles(fname,"r")

 lat1 = f1[0]->latitude
 lat2 = f1[1]->latitude

 z1 = 7.0*log(1.e3/f1[0]->levelist)
 z2 = 7.0*log(1.e3/f1[1]->levels)

 mm = (/ (/6,7,8/), (/12,1,2/) /)

do is=0, 1

 u1 = doubletofloat( \
         dim_avg_n( f1[0]->t(n0+ispan(mm(is,0)-1,n9,12),:,:), 0 ) + \
         dim_avg_n( f1[0]->t(n0+ispan(mm(is,1)-1,n9,12),:,:), 0 ) + \
         dim_avg_n( f1[0]->t(n0+ispan(mm(is,2)-1,n9,12),:,:), 0 ) )
 u1 = u1 / 3.
 u2 = dim_avg_n( f1[1]->t(n0+ispan(mm(is,0)-1,n9,12),:,:), 0 ) + \
      dim_avg_n( f1[1]->t(n0+ispan(mm(is,1)-1,n9,12),:,:), 0 ) + \
      dim_avg_n( f1[1]->t(n0+ispan(mm(is,2)-1,n9,12),:,:), 0 )
 u2 = u2 / 3.

;--------------------------------------- process 0

;--------------------------------------- cn level
 cnlev = fspan(-45,45,19)+245.
; cnfil = (/2,3,4,5,6,7,8,9,10,11,0,0,26,25,24,23,22,21,20,19,18,17/)
 cnfil = (/17,18,19,20,21,22,23,24,25,26,11,10,9,8,7,6,5,4,3,2/)
 cnlev2 = (/-20,-19,-18,-17,-16,-15,-14,-13,-12,-11,-10, \
            10,11,12,13,14,15,16,17,18,19,20/)*5.+245.

;--------------------------------------- contour
 wi = 32
 he = 40
 xf = 10 + (/0,1,2/)*40
 yf = 10 + (/0,1/)*52

 mstr  = new(ne,string)
 xstr  = new(ne,string)
 ystr  = new(2 ,string)
 lbopt = new((/ne,2/),string)

 mstr = (/"ERA-I","MERRA"/)
 xstr(:) = ""
 xstr(ne-1) = "lat (deg)"
 ystr(:) = ""
 ystr(0) = "log-p height (km)"
 lbopt(:,:) = ""
 lbopt(ne-1,0) = "(K)"

do ie=0, 1

 if (ie .eq. 0) then
   contour1 = cnshadeplot(wks,lat1,z1,u1,cnlev,cnfil,lbopt(ie,is))
   cnaddline(contour1,wks,lat1,z1,u1,cnlev2,0,"D",0,"T",0.5,"C",1)
 else
   contour1 = cnshadeplot(wks,lat2,z2,u2,cnlev,cnfil,lbopt(ie,is))
   cnaddline(contour1,wks,lat2,z2,u2,cnlev2,0,"D",0,"T",0.5,"C",1)
 end if
 setvalues contour1
  "pmLabelBarOrthogonalPosF" : 0.08
  "pmLabelBarParallelPosF"   : 1.15
  "pmLabelBarKeepAspect"     : False
  "pmLabelBarWidthF"         : 1.2
  "pmLabelBarHeightF"        : 0.12
  "lbLabelFontHeightF" : 0.026
  "lbTitlePosition"    : "Right"
  "lbTitleDirection"   : "Across"
  "lbTitleOffsetF"     : 0.02
  "lbTitleFontHeightF" : 0.032
 end setvalues

 set_size(contour1,xf(is),yf(ie),wi,he)
 axis_range(contour1,-90,90,0,50,False,False)
 axis_str(contour1,mstr(ie),xstr(ie),ystr(is))

 cnfinalize(contour1,wks)


 draw(contour1)

end do  ; ie

end do  ; is

 frame(wks)

end

