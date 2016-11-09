load "$UTIL/ncl/header.ncl"

xspace = 83   ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 1.2

begin

;--------------------------------------- parameters
 expname = (/"d","d"/)

 lat0 = 5.

 year1 = 1979
 year2 = 2010
 m1 = 0*12
 nmon0 = 12*21

 wks = cr_wks("fig/uss_tp-ra",600,0)
 gsn_define_colormap(wks,"StepSeq25")
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- read data
   fnam0 = "../dat/merra.u_ypt."+year1+"-"+year2+".nc"
   
   f3 = addfile(fnam0,"r")
   u3 = dim_avg( f3->u(m1:m1+nmon0-1,:,{-lat0:lat0}) )
   p  = f3->levels
   zp = 7.0*log(1.e3/p)

   nz3 = dimsizes(p)
   nt3 = nmon0

   ur = new((/nz3,nt3/),float)
   do n=0, nt3-1
     ur(:,n) = u3(n,:)
   end do

   t = fspan(0,13,14)

   uss0 = new((/2,nz3,14/),float)
   do n=1, 12
     uss0(0,:,n) = dim_avg(ur(:,ispan(n-1,nt3-1,12)))
   end do
   uss0(0,:,0 ) = uss0(0,:,12)
   uss0(0,:,13) = uss0(0,:,1 )

   fftu0 = ezfftf(uss0(0,:,1:12))
;   fftu0(:,:,1:) = 0.
   uss0(1,:,1:12) = ezfftb(fftu0,0.)
   uss0(1,:,0 ) = uss0(1,:,12)
   uss0(1,:,13) = uss0(1,:,1 )

;--------------------------------------- process 0

;--------------------------------------- cn level
 cnlev = fspan(-50,40,19)
; cnfil = (/2,3,4,5,6,7,8,9,10,11,0,0,26,25,24,23,22,21,20,19,18,17/)
 cnfil = (/17,18,19,20,21,22,23,24,25,26,0,0,11,10,9,8,7,6,5,4,3,2/)

;--------------------------------------- contour
 wi = 36
 he = 20
 xf = 12 + (/0,1,2/)*42
 yf = 10

 mstr = (/"MERRA",""/)
 xstr = (/"month"/)
 ystr = (/"log-p height (km)",""/)
 lbopt = (/"(m s:S:-1:N:)",""/)


do is=0, 1

 contour1 = cnshadeplot(wks,t,zp,uss0(is,:,:),cnlev,cnfil,lbopt(is))
 setvalues contour1
  "tiMainOffsetXF" : 0.078
  "pmLabelBarOrthogonalPosF" : 0.12
  "pmLabelBarParallelPosF"   : 1.1
  "pmLabelBarKeepAspect"     : False
  "pmLabelBarWidthF"         : 1.5
  "pmLabelBarHeightF"        : 0.2
;  "lbLabelAutoStride" : False
  "lbLabelFontHeightF" : 0.036
  "lbTitlePosition"    : "Right"
  "lbTitleDirection"   : "Across"
  "lbTitleOffsetF"     : 0.01
  "lbTitleFontHeightF" : 0.0432
 end setvalues

 set_size(contour1,xf(is),yf,wi,he)
 axis_range(contour1,0.5,12.5,35,65,False,False)
 axis_str(contour1,mstr(is),xstr,ystr(is))
 setvalues contour1
  "tmXBMode" : "Explicit" ;"Manual"
  "tmXBValues" : ispan(1,12,1)
  "tmXBLabels" : (/"J","F","M","A","M","J","J","A","S","O","N","D"/)
;  "tmXBTickSpacingF"  : 1
;  "tmXBTickStartF"    : 0
;  "tmXBMinorPerMajor" : 0
 end setvalues

 cnfinalize(contour1,wks)


 draw(contour1)

end do


 frame(wks)

end



