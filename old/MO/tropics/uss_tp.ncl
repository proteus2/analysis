load "$UTIL/ncl/header.ncl"

xspace = 83   ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 1.2

begin

;--------------------------------------- parameters
 expname = (/"fd","0"/)

 nmon = 4+144
 lat0 = 5.

 year1 = 1979
 year2 = 2010
 m1 = 0*12
 nmon0 = 12*12

 wks = cr_wks("uss_tp-"+expname(0),600,0)
 gsn_define_colormap(wks,"StepSeq25")
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- read data
   fnam1 = "../dat/uh-p_"+expname(0)+".nc"

   f1 = addfile(fnam1,"r")
   u1 = dim_avg( f1->u(4:nmon-1,:,{-lat0:lat0},0) )
   p1 = f1->p
   zp1 = 7.0*log(1.e3/p1)

   fnam2 = "../dat/uh-p_"+expname(1)+".nc"

   f2 = addfile(fnam2,"r")
   u2 = dim_avg( f2->u(4:nmon-1,:,{-lat0:lat0},0) )

   nz = dimsizes(zp1)
   nt = nmon-4

   fnam0 = "../dat/merra.u_ypt."+year1+"-"+year2+".nc"
   
   f3 = addfile(fnam0,"r")
   u3 = dim_avg( f3->u(m1:m1+nmon0-1,:,{-lat0:lat0}) )
   p  = f3->levels
   zp = 7.0*log(1.e3/p)

   nz3 = dimsizes(p)
   nt3 = nmon0

   u = new((/2,nz,nt/),float)
   do n=0, nt-1
     u(:,:,n) = (/u2(n,:),u1(n,:)/)
   end do

   ur = new((/nz3,nt3/),float)
   do n=0, nt3-1
     ur(:,n) = u3(n,:)
   end do

   t = fspan(0,13,14)

   uss = new((/2,2,nz,14/),float)
   do n=1, 12
     uss(0,:,:,n) = dim_avg(u(:,:,ispan(n-1,nt-1,12)))
   end do
   uss(0,:,:,0 ) = uss(0,:,:,12)
   uss(0,:,:,13) = uss(0,:,:,1 )

   uss0 = new((/2,nz3,14/),float)
   do n=1, 12
     uss0(0,:,n) = dim_avg(ur(:,ispan(n-1,nt3-1,12)))
   end do
   uss0(0,:,0 ) = uss0(0,:,12)
   uss0(0,:,13) = uss0(0,:,1 )

   fftu = ezfftf(uss(0,:,:,1:12))
   fftu(:,:,:,1:) = 0.
   uss(1,:,:,1:12) = ezfftb(fftu,0.)
   uss(1,:,:,0 ) = uss(1,:,:,12)
   uss(1,:,:,13) = uss(1,:,:,1 )

   fftu0 = ezfftf(uss0(0,:,1:12))
   fftu0(:,:,1:) = 0.
   uss0(1,:,1:12) = ezfftb(fftu0,0.)
   uss0(1,:,0 ) = uss0(1,:,12)
   uss0(1,:,13) = uss0(1,:,1 )

;--------------------------------------- process 0

;--------------------------------------- cn level
 cnlev = fspan(-50,50,21)
; cnfil = (/2,3,4,5,6,7,8,9,10,11,0,0,26,25,24,23,22,21,20,19,18,17/)
 cnfil = (/17,18,19,20,21,22,23,24,25,26,0,0,11,10,9,8,7,6,5,4,3,2/)

;--------------------------------------- contour
 wi = 30
 he = 24
 xf = 12 + (/0,1,2/)*38 + (/0,0,1/)*8
 yf = 10 + (/0,1,2,3,4/)*32

 mstr = (/ (/"CTL","SC","MERRA"/), (/"","",""/) /)
 xstr = (/"","month"/)
 ystr = (/"log-p height (km)","","log-p height (km)"/)
 lbopt = (/ (/"","",""/), (/"(m s:S:-1:N:)","",""/) /)


do is=0, 1
do iv=0, 2

 if (iv .le. 1) then
   contour1 = cnshadeplot(wks,t,zp1,uss(is,iv,:,:),cnlev,cnfil,lbopt(is,iv))
 else
   contour1 = cnshadeplot(wks,t,zp,uss0(is,:,:),cnlev,cnfil,lbopt(is,iv))
 end if
 setvalues contour1
  "pmLabelBarOrthogonalPosF" : 0.12
  "pmLabelBarParallelPosF"   : 1.8
  "pmLabelBarKeepAspect"     : False
  "pmLabelBarWidthF"         : 2.0
  "pmLabelBarHeightF"        : 0.12
  "lbLabelAutoStride" : False
  "lbLabelFontHeightF" : 0.03
  "lbTitlePosition"    : "Right"
  "lbTitleDirection"   : "Across"
  "lbTitleOffsetF"     : 0.02
  "lbTitleFontHeightF" : 0.04
 end setvalues

 set_size(contour1,xf(iv),yf(is),wi,he)
 axis_range(contour1,0.5,12.5,15,65,False,False)
 axis_str(contour1,mstr(is,iv),xstr(is),ystr(iv))
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
end do


 frame(wks)

end



