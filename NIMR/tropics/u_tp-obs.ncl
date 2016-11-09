load "$UTIL/ncl/header.ncl"

xspace = 170  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 1.2

begin

;--------------------------------------- parameters
 expname = (/"sonde","sonde_hr"/)

 y1 = 1979
 y2 = 2010

 wks = cr_wks("fig/u_tp-obs",600,0)
 gsn_define_colormap(wks,"StepSeq25")
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- read data
 ne = dimsizes(expname)

 fdir  = "/data18/kyh/qbo_obs"
 fname = fdir+"/qbo_"+expname+".nc"
 f1 = addfiles(fname,"r")

 u1 = f1[0]->u
 u2 = f1[1]->u

 t1 = f1[0]->t
 t2 = f1[1]->t

 z1 = 7.0*log(1.e3/f1[0]->p)
 z2 = 7.0*log(1.e3/f1[1]->p)

;   uss = new((/2,nz,12/),float)
;   do n=0, 11
;     uss(:,:,n) = dim_avg(u(:,:,ispan(n,nt-1,12)))
;   end do
;
;   udes = new(dimsizes(u),float)
;   do n=0, nt-1
;     udes(:,:,n) = u(:,:,n) - uss(:,:,mod(n,12))
;   end do

;--------------------------------------- process 0

;--------------------------------------- cn level
 cnlev = fspan(-40,40,21)
; cnfil = (/2,3,4,5,6,7,8,9,10,11,0,0,26,25,24,23,22,21,20,19,18,17/)
 cnfil = (/17,18,19,20,21,22,23,24,25,26,0,0,11,10,9,8,7,6,5,4,3,2/)

;--------------------------------------- contour
 wi = 72*2
 he = 20
 xf = 10 + (/0,0,0,0,0,0,0,0/)*78
 yf = 10 + (/1,2,3,4,5,6,7,8/)*32

 mstr  = new(ne,string)
 xstr  = new(ne,string)
 ystr  = new(ne,string)
 lbopt = new(ne,string)

 mstr = (/"Sonde data at three stations","Sonde data at Singapore (high res.)"/)
 xstr(:) = ""
 xstr(ne-1) = "year"
 ystr(:) = "height (km)"
 lbopt(:) = ""
 lbopt(ne-1) = "(m s:S:-1:N:)"

do ie=0, ne-1

 if ( ie .lt. 1 ) then
   contour0 = cnshadeplot(wks,t1,z1,u1,cnlev,cnfil,lbopt(ie))
   contour1 = cnshadeplot(wks,t1,z1,u1,cnlev,cnfil,lbopt(ie))
 else
   contour1 = cnshadeplot(wks,t2,z2,u2,cnlev,cnfil,lbopt(ie))
 end if
 setvalues contour1
  "pmLabelBarOrthogonalPosF" : 0.12
  "pmLabelBarKeepAspect"     : False
  "pmLabelBarWidthF"         : 0.6
  "pmLabelBarHeightF"        : 0.2
  "lbLabelFontHeightF" : 0.015
  "lbTitlePosition"    : "Right"
  "lbTitleDirection"   : "Across"
  "lbTitleOffsetF"     : 0.02
  "lbTitleFontHeightF" : 0.018
 end setvalues

 set_size(contour1,xf(ie),yf(ie),wi,he)
 axis_range(contour1,y1,y1+24,16,45,False,False)
 axis_str(contour1,mstr(ie),xstr(ie),ystr(ie))
 setvalues contour1
  "tmXBMode" : "Manual"
  "tmXBTickSpacingF"  : 1
;  "tmXBTickStartF"    : y1
  "tmXBMinorPerMajor" : 0
 end setvalues

 cnfinalize(contour1,wks)

 if ( ie .lt. 1 ) then
   setvalues contour0
    "pmLabelBarOrthogonalPosF" : 0.12
    "pmLabelBarKeepAspect"     : False
    "pmLabelBarWidthF"         : 0.6
    "pmLabelBarHeightF"        : 0.2
    "lbLabelFontHeightF" : 0.015
    "lbTitlePosition"    : "Right"
    "lbTitleDirection"   : "Across"
    "lbTitleOffsetF"     : 0.02
    "lbTitleFontHeightF" : 0.018
   end setvalues

   set_size(contour0,xf(ie),yf(0)-32,wi,he)
   axis_range(contour0,y1-24,y1,16,45,False,False)
   axis_str(contour0,mstr(ie),xstr(ie),ystr(ie))
   setvalues contour0
    "tmXBMode" : "Manual"
    "tmXBTickSpacingF"  : 1
;    "tmXBTickStartF"    : y1
    "tmXBMinorPerMajor" : 0
   end setvalues
   cnfinalize(contour0,wks)
   draw(contour0)
 end if


 draw(contour1)

end do


 frame(wks)

end

