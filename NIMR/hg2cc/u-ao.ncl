load "$UTIL/ncl/header.ncl"

xspace = 170  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 1.2

begin

;--------------------------------------- parameters
 iyr = 126
 imon = 1
 nmon = 17*12
 lat0 = 5.

 year1 = 1979
 year2 = 2010
 y1 = 7
 nmon0 = 17*12

 wks = cr_wks("u-ao",600,0)
 gsn_define_colormap(wks,"StepSeq25")
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- read data
   fnam1 = "../dat/u-p_cc_pictl.nc"

   f1 = addfile(fnam1,"r")
   u1 = dim_avg( f1->ua(iyr*12+imon:iyr*12+imon+nmon-1,:,{-lat0:lat0}) )
   p1 = f1->plev/1.e2
   zp1 = 7.0*log(1.e3/p1)

   fnam2 = "../dat/u_b.nc"

   f2 = addfile(fnam2,"r")
   u2 = dim_avg( f2->u(imon:imon+nmon-1,:,{-lat0:lat0},0) )
;   p2 = f2->p
;   zp2 = 7.0*log(1.e3/p2)
   zp2 = f2->hybrid_ht/1.e3

   nz1 = dimsizes(zp1)
   nz2 = dimsizes(zp2)
   nt = nmon

   fnam0 = "../../MO/dat/merra.u_ypt."+year1+"-"+year2+".nc"
   
   f3 = addfile(fnam0,"r")
   u3 = dim_avg( f3->u(y1*12:y1*12+nmon0-1,:,{-lat0:lat0}) )
   p  = f3->levels
   zp = 7.0*log(1.e3/p)

   nz3 = dimsizes(p)
   nt3 = nmon0

   uu1 = new((/nz1,nt/),float)
   uu2 = new((/nz2,nt/),float)
   do n=0, nt-1
     uu1(:,n) = u1(n,:)
     uu2(:,n) = u2(n,:)
   end do

   ur = new((/nz3,nt3/),float)
   do n=0, nt3-1
     ur(:,n) = u3(n,:)
   end do

   t = fspan(0,nt-1,nt)/12. + iyr
   t2 = fspan(0,nt-1,nt)/12.
   t3 = fspan(0,nt3-1,nt3)/12. + y1

;--------------------------------------- process 0

;--------------------------------------- cn level
 cnlev = fspan(-50,50,21)
; cnfil = (/2,3,4,5,6,7,8,9,10,11,0,0,26,25,24,23,22,21,20,19,18,17/)
 cnfil = (/17,18,19,20,21,22,23,24,25,26,0,0,11,10,9,8,7,6,5,4,3,2/)

;--------------------------------------- contour
 xmax = 17+(/iyr,0,y1/)
 wi = 150 *17./25.
 he = 20
 xf = 12
 yf = 10 + (/0,1,2/)*36

 mstr = (/"PICTL_1860 (MO)","PICTL_1986 (NIMR)","MERRA from 1979"/)
 xstr = "year"
 ystr = (/"log-p height (km)","height (km)","log-p height (km)"/)
 lbopt = (/"","","(m s:S:-1:N:)"/)


do iv=0, 2

 if (iv .le. 1) then
   if (iv .eq. 0) then
     contour1 = cnshadeplot(wks,t,zp1,uu1,cnlev,cnfil,lbopt(iv))
   else
     contour1 = cnshadeplot(wks,t2,zp2,uu2,cnlev,cnfil,lbopt(iv))
   end if
 else
   contour1 = cnshadeplot(wks,t3,zp,ur,cnlev,cnfil,lbopt(iv))
 end if
 setvalues contour1
  "pmLabelBarOrthogonalPosF" : 0.12
  "pmLabelBarParallelPosF"   : 0.5
  "pmLabelBarKeepAspect"     : False
  "pmLabelBarWidthF"         : 0.5
  "pmLabelBarHeightF"        : 0.08
  "lbLabelAutoStride" : False
  "lbLabelFontHeightF" : 0.008
  "lbTitlePosition"    : "Right"
  "lbTitleDirection"   : "Across"
  "lbTitleOffsetF"     : 0.02
  "lbTitleFontHeightF" : 0.01
 end setvalues

 set_size(contour1,xf,yf(iv),wi,he)
 axis_range(contour1,-999,xmax(iv),15,55,False,False)
 axis_str(contour1,mstr(iv),xstr,ystr(iv))
;setvalues contour1
; "tmXBMode" : "Manual"
; "tmXBTickSpacingF"  : 12
; "tmXBTickStartF"    : 0
; "tmXBMinorPerMajor" : 11
;end setvalues

 cnfinalize(contour1,wks)


 draw(contour1)

end do


 frame(wks)

end



