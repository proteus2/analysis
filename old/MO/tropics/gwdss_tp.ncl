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

 wks = cr_wks("fig/gwdss_tp-"+expname(0),600,0)
 gsn_define_colormap(wks,"StepSeq25")
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- read data
   fnam1 = "../dat/u_"+expname(0)+".nc"
   fnam1n = "../dat/ng_"+expname(0)+".nc"
   fnam1c = "../dat/cg_"+expname(0)+".nc"

   f1 = addfile(fnam1,"r")
   f1n = addfile(fnam1n,"r")
   f1c = addfile(fnam1c,"r")
   u1 = dim_avg( f1->u(4:nmon-1,:,{-lat0:lat0},0) )
   n1 = dim_avg( f1n->field424   (4:nmon-1,:,{-lat0:lat0},0) )*86400.
   c1 = dim_avg( f1c->unspecified(4:nmon-1,:,{-lat0:lat0},0) )*86400.
   zp1 = f1->hybrid_ht/1.e3

   fnam2 = "../dat/u_"+expname(1)+".nc"
   fnam2n = "../dat/ng_"+expname(1)+".nc"

   f2 = addfile(fnam2,"r")
   f2n = addfile(fnam2n,"r")
   u2 = dim_avg( f2->u(4:nmon-1,:,{-lat0:lat0},0) )
   n2 = dim_avg( f2n->field424(4:nmon-1,:,{-lat0:lat0},0) )*86400.

   nz = dimsizes(zp1)
   nt = nmon-4

   u = new((/2,nz,nt/),float)
   d = new((/2,nz,nt/),float)
   do n=0, nt-1
     u(:,:,n) = (/u2(n,:),u1(n,:)/)
     d(:,:,n) = (/n2(n,:),c1(n,:)+n1(n,:)/)
   end do

   t = fspan(0,13,14)

   uss = new((/2,2,nz,14/),float)
   dss = new((/2,2,nz,14/),float)
   do n=1, 12
     uss(0,:,:,n) = dim_avg(u(:,:,ispan(n-1,nt-1,12)))
     dss(0,:,:,n) = dim_avg(d(:,:,ispan(n-1,nt-1,12)))
   end do
   uss(0,:,:,0 ) = uss(0,:,:,12)
   uss(0,:,:,13) = uss(0,:,:,1 )
   dss(0,:,:,0 ) = dss(0,:,:,12)
   dss(0,:,:,13) = dss(0,:,:,1 )

   fftu = ezfftf(uss(0,:,:,1:12))
   fftu(:,:,:,1:) = 0.
   uss(1,:,:,1:12) = ezfftb(fftu,0.)
   uss(1,:,:,0 ) = uss(1,:,:,12)
   uss(1,:,:,13) = uss(1,:,:,1 )

   fftd = ezfftf(dss(0,:,:,1:12))
   fftd(:,:,:,1:) = 0.
   dss(1,:,:,1:12) = ezfftb(fftd,0.)
   dss(1,:,:,0 ) = dss(1,:,:,12)
   dss(1,:,:,13) = dss(1,:,:,1 )

   t2 = fspan(0.5,12.5,13)

   fss = new((/2,2,nz,13/),float)
   do n=0, 11
     fss(0,:,:,n) = (uss(0,:,:,n+1) - uss(0,:,:,n))/30.
   end do
   fss(0,:,:,12) = fss(0,:,:,0)

;--------------------------------------- process 0

;--------------------------------------- cn level
 cnlev = fspan(-50,40,19)
; cnfil = (/2,3,4,5,6,7,8,9,10,11,0,0,26,25,24,23,22,21,20,19,18,17/)
 cnfil = (/17,18,19,20,21,22,23,24,25,26,0,0,11,10,9,8,7,6,5,4,3,2/)

;--------------------------------------- contour
 wi = 30
 he = 20
 xf = 12 + (/0,1,2/)*38
 yf = 10 + (/0,1,2,3,4/)*32

 mstr = (/ (/"CTL","SC","SC - CTL"/), (/"","",""/) /)
 xstr = (/"","","month"/)
 ystr = (/"log-p height (km)",""/)
 lbopt = (/ (/"","","(m s:S:-1:N:)"/), (/"","",""/) /)


do is=0, 1
do iv=0, 2

 if (is .eq. 0) then
   if (iv .le. 1) then
     contour1 = cnshadeplot(wks,t,zp1,uss(is,iv,:,:),cnlev,cnfil,lbopt(is,iv))
     cnaddline(contour1,wks,t,zp1,dss(is,iv,:,:),cnlev*0.1,0, \
               "D",(/0,0,16/),"T",(/1.0,0.1,1.0/),"C",(/1,colind,1/))
   else
     contour1 = cnshadeplot(wks,t,zp1,uss(is,1,:,:)-uss(is,0,:,:),cnlev,cnfil,lbopt(is,iv))
     cnaddline(contour1,wks,t,zp1,dss(is,1,:,:)-dss(is,0,:,:),cnlev*0.1,0, \
               "D",(/0,0,16/),"T",(/1.0,0.1,1.0/),"C",(/1,colind,1/))
   end if
 else
   if (iv .le. 1) then
     contour1 = cnshadeplot(wks,t2,zp1,fss(0,iv,:,:),cnlev*0.02,cnfil,lbopt(is,iv))
;     cnaddline(contour1,wks,t2,zp1,dss(is,iv,:,:),cnlev*0.1,0, \
;               "D",(/0,0,16/),"T",(/1.0,0.1,1.0/),"C",(/1,colind,1/))
   else
     contour1 = cnshadeplot(wks,t2,zp1,fss(0,1,:,:)-fss(0,0,:,:),cnlev*0.02,cnfil,lbopt(is,iv))
;     cnaddline(contour1,wks,t2,zp1,dss(is,1,:,:)-dss(is,0,:,:),cnlev*0.1,0, \
;               "D",(/0,0,16/),"T",(/1.0,0.1,1.0/),"C",(/1,colind,1/))
   end if
 end if


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

 set_size(contour1,xf(is),yf(iv),wi,he)
 axis_range(contour1,0.5,12.5,15,65,False,False)
 axis_str(contour1,mstr(is,iv),xstr(iv),ystr(is))
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

