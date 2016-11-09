load "$UTIL/ncl/header.ncl"

xspace =  83  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 1.2

begin

;--------------------------------------- parameters
 iexp = 0
 expname = (/"era-int","merra"/)
 cexp = (/"ERA-Interim","MERRA"/)

 y1 = 1979
 y2 = 2010
 lat0 = 5.

 wks = cr_wks("fig/u_tp-ra1",600,0)
 gsn_define_colormap(wks,"StepSeq25")
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- read data
 ne = dimsizes(expname)
 nt = (y2-y1+1)*12

 fdir  = "../dat"
 fname = fdir+"/"+expname+".u_ypt."+y1+"-"+y2+".nc"
 f1 = addfiles(fname,"r")

 z1 = 7.0*log(1.e3/f1[0]->levelist)
 z2 = 7.0*log(1.e3/f1[1]->levels)

 u1 = transpose( doubletofloat(dim_avg( f1[0]->u(:,:,{-lat0:lat0}) )) )
 u2 = transpose( dim_avg( f1[1]->u(:,:,{-lat0:lat0}) ) )

 t = y1 + (fspan(1,nt,nt)-0.5)/12.

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
 if (iexp .eq. 0) then
   u = u1
   z = z1
 end if
 if (iexp .eq. 1) then
   u = u2
   z = z2
 end if

;--------------------------------------- cn level
 cnlev = fspan(-40,40,21)
; cnfil = (/2,3,4,5,6,7,8,9,10,11,0,0,26,25,24,23,22,21,20,19,18,17/)
 cnfil = (/17,18,19,20,21,22,23,24,25,26,0,0,11,10,9,8,7,6,5,4,3,2/)

;--------------------------------------- contour
 wi = 72*2
 he = 20
 xf = 10 + (/0,0,0,0,0,0,0,0/)*78
 yf = 10 + (/0,1,2,3,4,5,6,7/)*28

 mstr = (/"Zonal wind :F21:(5:F34:0:F21:N:F30:-:F21:5:F34:0:F21:S, monthly and zonal mean, "+cexp(iexp)+")",""/)
 xstr = (/"","year"/)
 ystr = "height (km)"
 lbopt = (/"","(m s:S:-1:N:)"/)

 ys = y1

do il=0, 1

 ye = ys+(y2-y1+1)/2

 contour1 = cnshadeplot(wks,t,z,u,cnlev,cnfil,lbopt(il))
 setvalues contour1
  "pmLabelBarOrthogonalPosF" : 0.12
  "pmLabelBarKeepAspect"     : False
  "pmLabelBarWidthF"         : 0.5
  "pmLabelBarHeightF"        : 0.15
  "lbLabelFontHeightF" : 0.010
  "lbTitlePosition"    : "Right"
  "lbTitleDirection"   : "Across"
  "lbTitleOffsetF"     : 0.02
  "lbTitleFontHeightF" : 0.012
 end setvalues

 set_size(contour1,xf(il),yf(il),wi,he)
 axis_range(contour1,ys,ye,15,45,False,False)
 axis_str(contour1,mstr(il),xstr(il),ystr)
 setvalues contour1
  "tmXBMode" : "Manual"
  "tmXBTickSpacingF"  : 1
  "tmXBMinorPerMajor" : 0
 end setvalues

 cnfinalize(contour1,wks)


 draw(contour1)

 ys = ye

end do


 frame(wks)

end



